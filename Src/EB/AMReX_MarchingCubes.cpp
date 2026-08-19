
#include <AMReX_Arena.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>

#include <AMReX_EB2_C.H>
#include <AMReX_MarchingCubes.H>
#include <AMReX_MFIter.H>
#include <AMReX_mc_jgt_table.H>

#include <algorithm>
#include <fstream>
#include <iomanip>

/*
 * http://thomas.lewiner.org/publication_page.php%EF%B9%96pubkey=marching_cubes_jgt.html
 *
 * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
 * @author  Math Dept, PUC-Rio
 * @version 0.2
 * @date    12/08/2002
 *
 * @brief   MarchingCubes Algorithm
 */

namespace amrex::MC {

namespace {

LookUpTable* h_table = nullptr;
LookUpTable* d_table = nullptr;

/**
 * Sign of a*c - b*d with a scale-free tie break: the products have the units
 * of the level set squared, so a tie is declared relative to their magnitude
 * rather than against an absolute epsilon.  Returns -1, 0 (tie) or 1.  Shared
 * by the MC33 face decider and the interior (tunnel) test so that both are
 * invariant under scaling of the implicit function.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int ambiguous_product_sign (Real a, Real c, Real b, Real d) noexcept
{
    Real const p = a*c;
    Real const q = b*d;
    Real const r = p - q;
    Real const tol = Real(8.0)*std::numeric_limits<Real>::epsilon()*(std::abs(p)+std::abs(q));
    if (std::abs(r) <= tol) { return 0; }
    return (r > 0.0_rt) ? 1 : -1;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool four_crossing_fluid_is_connected (Real const* levelset) noexcept
{
    // The tie break is scale free (see ambiguous_product_sign); |a|+|b| is
    // invariant under the index rotations/reversals neighboring cells apply
    // to a shared face.
    int const sign = ambiguous_product_sign(levelset[0], levelset[2], levelset[1], levelset[3]);
    if (sign == 0) {
        // MC33's test_face tie break always selects the positive material.
        // Expressing the decision in material terms makes it invariant under
        // the rotations/reversals used by neighboring cells on a shared face.
        return true;
    }
    return (levelset[0] > 0.0_rt) == (sign > 0);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
EB2::Type_t face_type (Real area, Real tolerance) noexcept
{
    if (area <= tolerance) {
        return EB2::Type::covered;
    } else if (area >= 1.0_rt-tolerance) {
        return EB2::Type::regular;
    } else {
        return EB2::Type::irregular;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real edge_intersection_fraction (Real lo, Real hi, Real exact) noexcept
{
    // Cleanup moves fluid nodes to exactly zero.  Those repaired crossings
    // belong at the node, rather than at the original STL intersection.
    if (lo == 0.0_rt) {
        return 0.0_rt;
    }
    if (hi == 0.0_rt) {
        return 1.0_rt;
    }
    return (exact < 0.0_rt || exact > 1.0_rt) ? lo/(lo-hi) : exact;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void consume_face_decision (Array4<int const> const& cell_data,
                            Box const& cell_box, int i, int j, int k,
                            int face, int& valid, int& connected) noexcept
{
    if (valid == 0 && cell_box.contains(i,j,k)) {
        int const bit = 1 << face;
        if ((cell_data(i,j,k,face_decision_valid_mask) & bit) != 0) {
            valid = 1;
            connected =
                (cell_data(i,j,k,face_fluid_connected_mask) & bit) != 0;
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void accumulate_face_decision (Array4<int const> const& cell_data,
                               Box const& cell_box, int i, int j, int k,
                               int face, int* error,
                               int& valid, int& connected) noexcept
{
    if (cell_box.contains(i,j,k)) {
        int const bit = 1 << face;
        int const valid_mask =
            cell_data(i,j,k,face_decision_valid_mask);
        if ((valid_mask & bit) != 0) {
            int const decision =
                (cell_data(i,j,k,face_fluid_connected_mask) & bit) != 0;
            if (valid != 0 && connected != decision) {
                Gpu::Atomic::AddNoRet(error,1);
            } else if (valid == 0) {
                connected = decision;
            }
            valid = 1;
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
GpuArray<int,2> resolved_face_decision (
    Array4<int const> const& cell_data, Box const& cell_box, int* error,
    int ilo, int jlo, int klo, int flo,
    int ihi, int jhi, int khi, int fhi) noexcept
{
    int valid = 0;
    int connected = 0;
    accumulate_face_decision(cell_data, cell_box, ilo, jlo, klo, flo,
                             error, valid, connected);
    accumulate_face_decision(cell_data, cell_box, ihi, jhi, khi, fhi,
                             error, valid, connected);
    return {valid,connected};
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool face_is_rejected (Real a, Real b, Real c, Real d,
                       Array4<int const> const& cell_data,
                       Box const& cell_box,
                       int ilo, int jlo, int klo, int flo,
                       int ihi, int jhi, int khi, int fhi) noexcept
{
    bool const f0 = a > 0.0_rt;
    bool const f1 = b > 0.0_rt;
    bool const f2 = c > 0.0_rt;
    bool const f3 = d > 0.0_rt;
    if (!(f0 == f2 && f1 == f3 && f0 != f1)) {
        return false;
    }

    int valid = 0;
    int connected = 0;
    consume_face_decision(cell_data, cell_box, ilo, jlo, klo, flo,
                          valid, connected);
    consume_face_decision(cell_data, cell_box, ihi, jhi, khi, fhi,
                          valid, connected);

    Real const values[4] = {a, b, c, d};
    bool const fluid_connected = valid != 0
        ? connected != 0
        : four_crossing_fluid_is_connected(values);
    return !fluid_connected;
}

//! Returns false when the face polygon is degenerate; the caller then routes
//! the face into the nodal repair set.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool cut_face_fraction (Real const* levelset, Real const* intersections,
                        Real& area,
                        Real& centroid_x, Real& centroid_y,
                        bool has_resolved_decision,
                        bool resolved_fluid_connected) noexcept
{
    // Safe sentinels for every early-error path.  The caller will reject the
    // geometry, but no downstream kernel may observe uninitialized storage.
    area = 0.0_rt;
    centroid_x = 0.0_rt;
    centroid_y = 0.0_rt;

    constexpr Real vertex_x[4] = {0.0_rt, 1.0_rt, 1.0_rt, 0.0_rt};
    constexpr Real vertex_y[4] = {0.0_rt, 0.0_rt, 1.0_rt, 1.0_rt};

    Real polygon_x[8];
    Real polygon_y[8];
    int polygon_size = 0;
    int crossing_count = 0;

    for (int n = 0; n < 4; ++n) {
        int const next = (n+1) % 4;
        bool const fluid = levelset[n] > 0.0_rt;
        bool const next_fluid = levelset[next] > 0.0_rt;

        if (fluid) {
            polygon_x[polygon_size] = vertex_x[n];
            polygon_y[polygon_size] = vertex_y[n];
            ++polygon_size;
        }
        if (fluid != next_fluid) {
            Real const alpha = edge_intersection_fraction(
                levelset[n], levelset[next], intersections[n]);
            polygon_x[polygon_size] =
                vertex_x[n] + alpha*(vertex_x[next]-vertex_x[n]);
            polygon_y[polygon_size] =
                vertex_y[n] + alpha*(vertex_y[next]-vertex_y[n]);
            ++polygon_size;
            ++crossing_count;
        }
    }

    if (crossing_count == 0) {
        area = (levelset[0] > 0.0_rt) ? 1.0_rt : 0.0_rt;
        centroid_x = 0.0_rt;
        centroid_y = 0.0_rt;
        return true;
    }

    if (crossing_count == 4) {
        Real edge_x[4];
        Real edge_y[4];
        for (int n = 0; n < 4; ++n) {
            int const next = (n+1) % 4;
            Real const alpha = edge_intersection_fraction(
                levelset[n], levelset[next], intersections[n]);
            edge_x[n] = vertex_x[n]
                + alpha*(vertex_x[next]-vertex_x[n]);
            edge_y[n] = vertex_y[n]
                + alpha*(vertex_y[next]-vertex_y[n]);
        }

        bool const fluid_connected = has_resolved_decision
            ? resolved_fluid_connected
            : four_crossing_fluid_is_connected(levelset);
        Real fluid_area = fluid_connected ? 1.0_rt : 0.0_rt;
        Real fluid_moment_x = fluid_connected ? 0.5_rt : 0.0_rt;
        Real fluid_moment_y = fluid_connected ? 0.5_rt : 0.0_rt;

        // The asymptotic decider selects which diagonal pair is connected.
        // With straight MC face segments, the other sign consists of the two
        // corner triangles. Add the fluid triangles or subtract the covered
        // triangles from the unit square.
        for (int n = 0; n < 4; ++n) {
            bool const fluid = levelset[n] > 0.0_rt;
            if (fluid != fluid_connected) {
                int const previous = (n+3) % 4;
                Real const ax = edge_x[n] - vertex_x[n];
                Real const ay = edge_y[n] - vertex_y[n];
                Real const bx = edge_x[previous] - vertex_x[n];
                Real const by = edge_y[previous] - vertex_y[n];
                Real const triangle_area =
                    0.5_rt*std::abs(ax*by-ay*bx);
                Real const triangle_centroid_x =
                    (vertex_x[n]+edge_x[n]+edge_x[previous])/3.0_rt;
                Real const triangle_centroid_y =
                    (vertex_y[n]+edge_y[n]+edge_y[previous])/3.0_rt;
                Real const sign = fluid_connected ? -1.0_rt : 1.0_rt;
                fluid_area += sign*triangle_area;
                fluid_moment_x += sign*triangle_area*triangle_centroid_x;
                fluid_moment_y += sign*triangle_area*triangle_centroid_y;
            }
        }

        if (fluid_area <= 0.0_rt || fluid_area >= 1.0_rt) {
            return false;
        }
        area = fluid_area;
        centroid_x = fluid_moment_x/fluid_area - 0.5_rt;
        centroid_y = fluid_moment_y/fluid_area - 0.5_rt;
        return true;
    }

    if (crossing_count != 2 || polygon_size < 3) {
        return false;
    }

    Real twice_area = 0.0_rt;
    Real centroid_x_numerator = 0.0_rt;
    Real centroid_y_numerator = 0.0_rt;
    for (int n = 0; n < polygon_size; ++n) {
        int const next = (n+1) % polygon_size;
        Real const cross = polygon_x[n]*polygon_y[next]
            - polygon_x[next]*polygon_y[n];
        twice_area += cross;
        centroid_x_numerator += (polygon_x[n]+polygon_x[next])*cross;
        centroid_y_numerator += (polygon_y[n]+polygon_y[next])*cross;
    }

    if (twice_area <= 0.0_rt) {
        return false;
    }

    area = 0.5_rt*twice_area;
    centroid_x = centroid_x_numerator/(3.0_rt*twice_area) - 0.5_rt;
    centroid_y = centroid_y_numerator/(3.0_rt*twice_area) - 0.5_rt;
    return true;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void process_cube (std::int8_t ipass, LookUpTable const* lut, int i, int j, int k,
                   Array4<Real const> const& sdf, Array4<int const> const& ex,
                   Array4<int const> const& ey, Array4<int const> const& ez,
                   GpuArray<Real*,3> const& pvrtx, Array4<int> const& ntri,
                   GpuArray<int*,3> const& ptri, int* error)
{
    std::uint8_t lut_entry = 0;
    Real cube[8] = { sdf(i  ,j  ,k  ), sdf(i+1,j  ,k  ), sdf(i+1,j+1,k  ), sdf(i  ,j+1,k  ),
                     sdf(i  ,j  ,k+1), sdf(i+1,j  ,k+1), sdf(i+1,j+1,k+1), sdf(i  ,j+1,k+1) };
    if (cube[0] > 0) { lut_entry +=   1; }
    if (cube[1] > 0) { lut_entry +=   2; }
    if (cube[2] > 0) { lut_entry +=   4; }
    if (cube[3] > 0) { lut_entry +=   8; }
    if (cube[4] > 0) { lut_entry +=  16; }
    if (cube[5] > 0) { lut_entry +=  32; }
    if (cube[6] > 0) { lut_entry +=  64; }
    if (cube[7] > 0) { lut_entry += 128; }

    int face_valid_mask = 0;
    int face_connected_mask = 0;

    auto add_triangle = [&] (const std::int8_t* trig, std::int8_t n, int v12 = -1) -> int
    {
        if (ipass == 0) { return n; }

        int r = 0;
        int tv[3] = {};

        for( int t = 0 ; t < 3*n ; t++ )
        {
            switch( trig[t] )
            {
            case  0 : tv[ t % 3 ] = ex(i  ,j  ,k  ,1) ; break ;
            case  1 : tv[ t % 3 ] = ey(i+1,j  ,k  ,1) ; break ;
            case  2 : tv[ t % 3 ] = ex(i  ,j+1,k  ,1) ; break ;
            case  3 : tv[ t % 3 ] = ey(i  ,j  ,k  ,1) ; break ;
            case  4 : tv[ t % 3 ] = ex(i  ,j  ,k+1,1) ; break ;
            case  5 : tv[ t % 3 ] = ey(i+1,j  ,k+1,1) ; break ;
            case  6 : tv[ t % 3 ] = ex(i  ,j+1,k+1,1) ; break ;
            case  7 : tv[ t % 3 ] = ey(i  ,j  ,k+1,1) ; break ;
            case  8 : tv[ t % 3 ] = ez(i  ,j  ,k  ,1) ; break ;
            case  9 : tv[ t % 3 ] = ez(i+1,j  ,k  ,1) ; break ;
            case 10 : tv[ t % 3 ] = ez(i+1,j+1,k  ,1) ; break ;
            case 11 : tv[ t % 3 ] = ez(i  ,j+1,k  ,1) ; break ;
            case 12 : tv[ t % 3 ] = v12 ; break ;
            default : break ;
            }

            if( tv[t%3] == -1 ) { *error = 1; return -1; }

            if( t%3 == 2 )
            {
                auto m = ntri(i,j,k,triangle_offset) + r;
                ptri[0][m] = tv[0];
                ptri[1][m] = tv[1];
                ptri[2][m] = tv[2];
                ++r;
            }
        }

        return r;
    };

    auto test_face = [&] (std::int8_t face) -> bool
    {
        Real A,B,C,D ;

        switch( face )
        {
        case -1 : case 1 :  A = cube[0] ;  B = cube[4] ;  C = cube[5] ;  D = cube[1] ;  break ;
        case -2 : case 2 :  A = cube[1] ;  B = cube[5] ;  C = cube[6] ;  D = cube[2] ;  break ;
        case -3 : case 3 :  A = cube[2] ;  B = cube[6] ;  C = cube[7] ;  D = cube[3] ;  break ;
        case -4 : case 4 :  A = cube[3] ;  B = cube[7] ;  C = cube[4] ;  D = cube[0] ;  break ;
        case -5 : case 5 :  A = cube[0] ;  B = cube[3] ;  C = cube[2] ;  D = cube[1] ;  break ;
        case -6 : case 6 :  A = cube[4] ;  B = cube[7] ;  C = cube[6] ;  D = cube[5] ;  break ;
        default : *error = 1 ;  A = B = C = D = 0 ;
        };

        Real const values[4] = {A,B,C,D};
        bool const fluid_connected = four_crossing_fluid_is_connected(values);
        int const bit = 1 << (std::abs(int(face))-1);
        face_valid_mask |= bit;
        if (fluid_connected) {
            face_connected_mask |= bit;
        }
        return fluid_connected == (face > 0);
    };

    auto test_interior_impl = [&] (std::int8_t _case, std::int8_t _config,
                                   std::int8_t _subconfig, std::int8_t s)
    {
        Real t, At=0, Bt=0, Ct=0, Dt=0, a, b ;
        std::int8_t  test =  0 ;
        std::int8_t  edge = -1 ; // reference edge of the triangulation

        switch( _case )
        {
        case  4 :
        case 10 :
            a = ( cube[4] - cube[0] ) * ( cube[6] - cube[2] ) - ( cube[7] - cube[3] ) * ( cube[5] - cube[1] ) ;
            b =  cube[2] * ( cube[4] - cube[0] ) + cube[0] * ( cube[6] - cube[2] )
                - cube[1] * ( cube[7] - cube[3] ) - cube[3] * ( cube[5] - cube[1] ) ;
            t = - b / (2*a) ;
            if( t<0 || t>1 ) { return s>0 ; }

            At = cube[0] + ( cube[4] - cube[0] ) * t ;
            Bt = cube[3] + ( cube[7] - cube[3] ) * t ;
            Ct = cube[2] + ( cube[6] - cube[2] ) * t ;
            Dt = cube[1] + ( cube[5] - cube[1] ) * t ;
            break ;

        case  6 :
        case  7 :
        case 12 :
        case 13 :
            // Lewiner's reference-edge shortcut: the interior test is taken on
            // the plane through the reference edge with At = 0.  With that
            // choice the tunnel alternatives 6.1.2, 7.4.2, 12.1.2 and 13.5.2
            // are never selected (Bt >= 0 is equivalent to the face test that
            // has already failed), so these cases always take the split
            // tiling.  That is conservative for the EB builder: split cells
            // hold two fluid corner groups and are rejected and repaired by
            // mark_cells_for_cleanup, exactly as tunnel cells would be.
            switch( _case ) // NOLINT(bugprone-switch-missing-default-case)
            {
            case  6 : edge = lut->test6 [_config][2] ; break ;
            case  7 : edge = lut->test7 [_config][4] ; break ;
            case 12 : edge = lut->test12[_config][3] ; break ;
            case 13 : edge = lut->tiling13_5_1[_config][_subconfig][0] ; break ;
            }
            switch( edge )
            {
            case  0 :
                t  = cube[0] / ( cube[0] - cube[1] ) ;
                At = 0 ;
                Bt = cube[3] + ( cube[2] - cube[3] ) * t ;
                Ct = cube[7] + ( cube[6] - cube[7] ) * t ;
                Dt = cube[4] + ( cube[5] - cube[4] ) * t ;
                break ;
            case  1 :
                t  = cube[1] / ( cube[1] - cube[2] ) ;
                At = 0 ;
                Bt = cube[0] + ( cube[3] - cube[0] ) * t ;
                Ct = cube[4] + ( cube[7] - cube[4] ) * t ;
                Dt = cube[5] + ( cube[6] - cube[5] ) * t ;
                break ;
            case  2 :
                t  = cube[2] / ( cube[2] - cube[3] ) ;
                At = 0 ;
                Bt = cube[1] + ( cube[0] - cube[1] ) * t ;
                Ct = cube[5] + ( cube[4] - cube[5] ) * t ;
                Dt = cube[6] + ( cube[7] - cube[6] ) * t ;
                break ;
            case  3 :
                t  = cube[3] / ( cube[3] - cube[0] ) ;
                At = 0 ;
                Bt = cube[2] + ( cube[1] - cube[2] ) * t ;
                Ct = cube[6] + ( cube[5] - cube[6] ) * t ;
                Dt = cube[7] + ( cube[4] - cube[7] ) * t ;
                break ;
            case  4 :
                t  = cube[4] / ( cube[4] - cube[5] ) ;
                At = 0 ;
                Bt = cube[7] + ( cube[6] - cube[7] ) * t ;
                Ct = cube[3] + ( cube[2] - cube[3] ) * t ;
                Dt = cube[0] + ( cube[1] - cube[0] ) * t ;
                break ;
            case  5 :
                t  = cube[5] / ( cube[5] - cube[6] ) ;
                At = 0 ;
                Bt = cube[4] + ( cube[7] - cube[4] ) * t ;
                Ct = cube[0] + ( cube[3] - cube[0] ) * t ;
                Dt = cube[1] + ( cube[2] - cube[1] ) * t ;
                break ;
            case  6 :
                t  = cube[6] / ( cube[6] - cube[7] ) ;
                At = 0 ;
                Bt = cube[5] + ( cube[4] - cube[5] ) * t ;
                Ct = cube[1] + ( cube[0] - cube[1] ) * t ;
                Dt = cube[2] + ( cube[3] - cube[2] ) * t ;
                break ;
            case  7 :
                t  = cube[7] / ( cube[7] - cube[4] ) ;
                At = 0 ;
                Bt = cube[6] + ( cube[5] - cube[6] ) * t ;
                Ct = cube[2] + ( cube[1] - cube[2] ) * t ;
                Dt = cube[3] + ( cube[0] - cube[3] ) * t ;
                break ;
            case  8 :
                t  = cube[0] / ( cube[0] - cube[4] ) ;
                At = 0 ;
                Bt = cube[3] + ( cube[7] - cube[3] ) * t ;
                Ct = cube[2] + ( cube[6] - cube[2] ) * t ;
                Dt = cube[1] + ( cube[5] - cube[1] ) * t ;
                break ;
            case  9 :
                t  = cube[1] / ( cube[1] - cube[5] ) ;
                At = 0 ;
                Bt = cube[0] + ( cube[4] - cube[0] ) * t ;
                Ct = cube[3] + ( cube[7] - cube[3] ) * t ;
                Dt = cube[2] + ( cube[6] - cube[2] ) * t ;
                break ;
            case 10 :
                t  = cube[2] / ( cube[2] - cube[6] ) ;
                At = 0 ;
                Bt = cube[1] + ( cube[5] - cube[1] ) * t ;
                Ct = cube[0] + ( cube[4] - cube[0] ) * t ;
                Dt = cube[3] + ( cube[7] - cube[3] ) * t ;
                break ;
            case 11 :
                t  = cube[3] / ( cube[3] - cube[7] ) ;
                At = 0 ;
                Bt = cube[2] + ( cube[6] - cube[2] ) * t ;
                Ct = cube[1] + ( cube[5] - cube[1] ) * t ;
                Dt = cube[0] + ( cube[4] - cube[0] ) * t ;
                break ;
            default : *error = 1;  break ;
            }
            break ;

        default : *error = 1;  break ;
        }

        if( At >= 0 ) { test ++ ; }
        if( Bt >= 0 ) { test += 2 ; }
        if( Ct >= 0 ) { test += 4 ; }
        if( Dt >= 0 ) { test += 8 ; }
        switch( test ) // NOLINT(bugprone-switch-missing-default-case)
        {
        case  0 : return s>0 ;
        case  1 : return s>0 ;
        case  2 : return s>0 ;
        case  3 : return s>0 ;
        case  4 : return s>0 ;
        case  5 : if( ambiguous_product_sign(At,Ct,Bt,Dt) <= 0 ) { return s>0 ; } break;
        case  6 : return s>0 ;
        case  7 : return s<0 ;
        case  8 : return s>0 ;
        case  9 : return s>0 ;
        case 10 : if( ambiguous_product_sign(At,Ct,Bt,Dt) >= 0 ) { return s>0 ; } break;
        case 11 : return s<0 ;
        case 12 : return s>0 ;
        case 13 : return s<0 ;
        case 14 : return s<0 ;
        case 15 : return s<0 ;
        }

        return s<0 ;
    };

    auto test_interior = [&] (std::int8_t c, std::int8_t config,
                              std::int8_t subconfig, std::int8_t s) -> bool
    {
        return test_interior_impl(c,config,subconfig,s);
    };

    auto add_c_vertex = [&] () -> int
    {
        if (ipass == 0) {
            ntri(i,j,k,interior_vertex_count) = 1;
            return -1;
        }

        Real u = 0 ;
        int   vid ;

        auto m = ntri(i,j,k,interior_vertex_offset);
        auto& vert_x  = pvrtx[0][m];
        auto& vert_y  = pvrtx[1][m];
        auto& vert_z  = pvrtx[2][m];

        vert_x = vert_y = vert_z = 0 ;

        auto update_vertex = [&] () {
            ++u ;
            vert_x  += pvrtx[0][vid];
            vert_y  += pvrtx[1][vid];
            vert_z  += pvrtx[2][vid];
        };


        // Computes the average of the intersection points of the cube
        vid = ex( i , j , k ,1) ;
        if( vid != -1 ) { update_vertex(); }
        vid = ey(i+1, j , k ,1) ;
        if( vid != -1 ) { update_vertex(); }
        vid = ex( i ,j+1, k ,1) ;
        if( vid != -1 ) { update_vertex(); }
        vid = ey( i , j , k ,1) ;
        if( vid != -1 ) { update_vertex(); }
        vid = ex( i , j ,k+1,1) ;
        if( vid != -1 ) { update_vertex(); }
        vid = ey(i+1, j ,k+1,1) ;
        if( vid != -1 ) { update_vertex(); }
        vid = ex( i ,j+1,k+1,1) ;
        if( vid != -1 ) { update_vertex(); }
        vid = ey( i , j ,k+1,1) ;
        if( vid != -1 ) { update_vertex(); }
        vid = ez( i , j , k ,1) ;
        if( vid != -1 ) { update_vertex(); }
        vid = ez(i+1, j , k ,1) ;
        if( vid != -1 ) { update_vertex(); }
        vid = ez(i+1,j+1, k ,1) ;
        if( vid != -1 ) { update_vertex(); }
        vid = ez( i ,j+1, k ,1) ;
        if( vid != -1 ) { update_vertex(); }

        vert_x  *= Real(1)/u ;
        vert_y  *= Real(1)/u ;
        vert_z  *= Real(1)/u ;

        return m;
    };

    int v12 = -1;
    std::int8_t _case   = lut->cases[lut_entry][0];
    std::int8_t _config = lut->cases[lut_entry][1];
    std::int8_t _subconfig = 0;
    int nt = 0;

    switch( _case )// NOLINT(bugprone-switch-missing-default-case)
    {
    case  0 :
        break ;

    case  1 :
        nt = add_triangle( lut->tiling1[_config], 1 ) ;
        break ;

    case  2 :
        nt = add_triangle( lut->tiling2[_config], 2 ) ;
        break ;

    case  3 :
        if( test_face( lut->test3[_config]) ) {
            nt = add_triangle( lut->tiling3_2[_config], 4 ) ; // 3.2
        } else {
            nt = add_triangle( lut->tiling3_1[_config], 2 ) ; // 3.1
        }
        break ;

    case  4 :
        if( test_interior( _case, _config, _subconfig, lut->test4[_config]) ) {
            nt = add_triangle( lut->tiling4_1[_config], 2 ) ; // 4.1.1
        } else {
            nt = add_triangle( lut->tiling4_2[_config], 6 ) ; // 4.1.2
        }
        break ;

    case  5 :
        nt = add_triangle( lut->tiling5[_config], 3 ) ;
        break ;

    case  6 :
        if( test_face( lut->test6[_config][0]) ) {
            nt = add_triangle( lut->tiling6_2[_config], 5 ) ; // 6.2
        } else {
            if( test_interior( _case, _config, _subconfig, lut->test6[_config][1]) ) {
                nt = add_triangle( lut->tiling6_1_1[_config], 3 ) ; // 6.1.1
            } else {
                v12 = add_c_vertex() ;
                nt = add_triangle( lut->tiling6_1_2[_config], 9 , v12) ; // 6.1.2
            }
        }
        break ;

    case  7 :
        if( test_face( lut->test7[_config][0] ) ) { _subconfig +=  1 ; }
        if( test_face( lut->test7[_config][1] ) ) { _subconfig +=  2 ; }
        if( test_face( lut->test7[_config][2] ) ) { _subconfig +=  4 ; }
        switch( _subconfig ) // NOLINT(bugprone-switch-missing-default-case)
        {
        case 0 :
            nt = add_triangle( lut->tiling7_1[_config], 3 ) ; break ;
        case 1 :
            nt = add_triangle( lut->tiling7_2[_config][0], 5 ) ; break ;
        case 2 :
            nt = add_triangle( lut->tiling7_2[_config][1], 5 ) ; break ;
        case 3 :
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling7_3[_config][0], 9, v12 ) ; break ;
        case 4 :
            nt = add_triangle( lut->tiling7_2[_config][2], 5 ) ; break ;
        case 5 :
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling7_3[_config][1], 9, v12 ) ; break ;
        case 6 :
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling7_3[_config][2], 9, v12 ) ; break ;
        case 7 :
            if( test_interior( _case, _config, _subconfig, lut->test7[_config][3]) ) {
                nt = add_triangle( lut->tiling7_4_2[_config], 9 ) ;
            } else {
                nt = add_triangle( lut->tiling7_4_1[_config], 5 ) ;
            }
            break ;
        };
        break ;

    case  8 :
        nt = add_triangle( lut->tiling8[_config], 2 ) ;
        break ;

    case  9 :
        nt = add_triangle( lut->tiling9[_config], 4 ) ;
        break ;

    case 10 :
        if( test_face( lut->test10[_config][0]) ) {
            if( test_face( lut->test10[_config][1]) ) {
                nt = add_triangle( lut->tiling10_1_1_[_config], 4 ) ; // 10.1.1
            } else {
                v12 = add_c_vertex() ;
                nt = add_triangle( lut->tiling10_2[_config], 8, v12 ) ; // 10.2
            }
        } else {
            if( test_face( lut->test10[_config][1]) ) {
                v12 = add_c_vertex() ;
                nt = add_triangle( lut->tiling10_2_[_config], 8, v12 ) ; // 10.2
            } else {
                if( test_interior( _case, _config, _subconfig, lut->test10[_config][2]) ) {
                    nt = add_triangle( lut->tiling10_1_1[_config], 4 ) ; // 10.1.1
                } else {
                    nt = add_triangle( lut->tiling10_1_2[_config], 8 ) ; // 10.1.2
                }
            }
        }
        break ;

    case 11 :
        nt = add_triangle( lut->tiling11[_config], 4 ) ;
        break ;

    case 12 :
        if( test_face( lut->test12[_config][0]) ) {
            if( test_face( lut->test12[_config][1]) ) {
                nt = add_triangle( lut->tiling12_1_1_[_config], 4 ) ; // 12.1.1
            } else {
                v12 = add_c_vertex() ;
                nt = add_triangle( lut->tiling12_2[_config], 8, v12 ) ; // 12.2
            }
        } else {
            if( test_face( lut->test12[_config][1]) ) {
                v12 = add_c_vertex() ;
                nt = add_triangle( lut->tiling12_2_[_config], 8, v12 ) ; // 12.2
            } else {
                if( test_interior( _case, _config, _subconfig, lut->test12[_config][2]) ) {
                    nt = add_triangle( lut->tiling12_1_1[_config], 4 ) ; // 12.1.1
                } else {
                    nt = add_triangle( lut->tiling12_1_2[_config], 8 ) ; // 12.1.2
                }
            }
        }
        break ;

    case 13 :
        if( test_face( lut->test13[_config][0] ) ) { _subconfig +=  1 ; }
        if( test_face( lut->test13[_config][1] ) ) { _subconfig +=  2 ; }
        if( test_face( lut->test13[_config][2] ) ) { _subconfig +=  4 ; }
        if( test_face( lut->test13[_config][3] ) ) { _subconfig +=  8 ; }
        if( test_face( lut->test13[_config][4] ) ) { _subconfig += 16 ; }
        if( test_face( lut->test13[_config][5] ) ) { _subconfig += 32 ; }
        switch( lut->subconfig13[_subconfig] ) // NOLINT(bugprone-switch-missing-default-case)
        {
        case 0 :/* 13.1 */
            nt = add_triangle( lut->tiling13_1[_config], 4 ) ; break ;

        case 1 :/* 13.2 */
            nt = add_triangle( lut->tiling13_2[_config][0], 6 ) ; break ;
        case 2 :/* 13.2 */
            nt = add_triangle( lut->tiling13_2[_config][1], 6 ) ; break ;
        case 3 :/* 13.2 */
            nt = add_triangle( lut->tiling13_2[_config][2], 6 ) ; break ;
        case 4 :/* 13.2 */
            nt = add_triangle( lut->tiling13_2[_config][3], 6 ) ; break ;
        case 5 :/* 13.2 */
            nt = add_triangle( lut->tiling13_2[_config][4], 6 ) ; break ;
        case 6 :/* 13.2 */
            nt = add_triangle( lut->tiling13_2[_config][5], 6 ) ; break ;

        case 7 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3[_config][0], 10, v12 ) ; break ;
        case 8 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3[_config][1], 10, v12 ) ; break ;
        case 9 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3[_config][2], 10, v12 ) ; break ;
        case 10 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3[_config][3], 10, v12 ) ; break ;
        case 11 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3[_config][4], 10, v12 ) ; break ;
        case 12 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3[_config][5], 10, v12 ) ; break ;
        case 13 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3[_config][6], 10, v12 ) ; break ;
        case 14 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3[_config][7], 10, v12 ) ; break ;
        case 15 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3[_config][8], 10, v12 ) ; break ;
        case 16 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3[_config][9], 10, v12 ) ; break ;
        case 17 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3[_config][10], 10, v12 ) ; break ;
        case 18 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3[_config][11], 10, v12 ) ; break ;

        case 19 :/* 13.4 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_4[_config][0], 12, v12 ) ; break ;
        case 20 :/* 13.4 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_4[_config][1], 12, v12 ) ; break ;
        case 21 :/* 13.4 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_4[_config][2], 12, v12 ) ; break ;
        case 22 :/* 13.4 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_4[_config][3], 12, v12 ) ; break ;

        case 23 :/* 13.5 */
            _subconfig = 0 ;
            if( test_interior( _case, _config, _subconfig, lut->test13[_config][6] ) ) {
                nt = add_triangle( lut->tiling13_5_1[_config][0], 6 ) ;
            } else {
                nt = add_triangle( lut->tiling13_5_2[_config][0], 10 ) ;
            }
            break ;
        case 24 :/* 13.5 */
            _subconfig = 1 ;
            if( test_interior( _case, _config, _subconfig, lut->test13[_config][6] ) ) {
                nt = add_triangle( lut->tiling13_5_1[_config][1], 6 ) ;
            } else {
                nt = add_triangle( lut->tiling13_5_2[_config][1], 10 ) ;
            }
            break ;
        case 25 :/* 13.5 */
            _subconfig = 2 ;
            if( test_interior( _case, _config, _subconfig, lut->test13[_config][6] ) ) {
                nt = add_triangle( lut->tiling13_5_1[_config][2], 6 ) ;
            } else {
                nt = add_triangle( lut->tiling13_5_2[_config][2], 10 ) ;
            }
            break ;
        case 26 :/* 13.5 */
            _subconfig = 3 ;
            if( test_interior( _case, _config, _subconfig, lut->test13[_config][6] ) ) {
                nt = add_triangle( lut->tiling13_5_1[_config][3], 6 ) ;
            } else {
                nt = add_triangle( lut->tiling13_5_2[_config][3], 10 ) ;
            }
            break ;

        case 27 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3_[_config][0], 10, v12 ) ; break ;
        case 28 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3_[_config][1], 10, v12 ) ; break ;
        case 29 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3_[_config][2], 10, v12 ) ; break ;
        case 30 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3_[_config][3], 10, v12 ) ; break ;
        case 31 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3_[_config][4], 10, v12 ) ; break ;
        case 32 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3_[_config][5], 10, v12 ) ; break ;
        case 33 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3_[_config][6], 10, v12 ) ; break ;
        case 34 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3_[_config][7], 10, v12 ) ; break ;
        case 35 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3_[_config][8], 10, v12 ) ; break ;
        case 36 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3_[_config][9], 10, v12 ) ; break ;
        case 37 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3_[_config][10], 10, v12 ) ; break ;
        case 38 :/* 13.3 */
            v12 = add_c_vertex() ;
            nt = add_triangle( lut->tiling13_3_[_config][11], 10, v12 ) ; break ;

        case 39 :/* 13.2 */
            nt = add_triangle( lut->tiling13_2_[_config][0], 6 ) ; break ;
        case 40 :/* 13.2 */
            nt = add_triangle( lut->tiling13_2_[_config][1], 6 ) ; break ;
        case 41 :/* 13.2 */
            nt = add_triangle( lut->tiling13_2_[_config][2], 6 ) ; break ;
        case 42 :/* 13.2 */
            nt = add_triangle( lut->tiling13_2_[_config][3], 6 ) ; break ;
        case 43 :/* 13.2 */
            nt = add_triangle( lut->tiling13_2_[_config][4], 6 ) ; break ;
        case 44 :/* 13.2 */
            nt = add_triangle( lut->tiling13_2_[_config][5], 6 ) ; break ;

        case 45 :/* 13.1 */
            nt = add_triangle( lut->tiling13_1_[_config], 4 ) ; break ;
        }
        break ;

    case 14 :
        nt = add_triangle( lut->tiling14[_config], 4 ) ;
        break ;
    };

    if (ipass == 0) {
        ntri(i,j,k,triangle_count) = nt;
        ntri(i,j,k,face_decision_valid_mask) = face_valid_mask;
        ntri(i,j,k,face_fluid_connected_mask) = face_connected_mask;
    }
}

}

void Initialize ()
{
    if (h_table == nullptr) {
        h_table = new LookUpTable{};
#ifdef AMREX_USE_GPU
        d_table = (LookUpTable*) The_Arena()->alloc(sizeof(LookUpTable));
        Gpu::htod_memcpy(d_table, h_table, sizeof(LookUpTable));
#else
        d_table = h_table;
#endif
    }
}

void Finalize ()
{
    if (h_table) {
        delete h_table;
        h_table = nullptr;
#ifdef AMREX_USE_GPU
        Gpu::streamSynchronizeAll();
        The_Arena()->free(d_table);
#endif
        d_table = nullptr;
    }
}

void Vertex::resize (int n)
{
    x.resize(n);
    y.resize(n);
    z.resize(n);
}

void Triangle::resize (int n)
{
    v1.resize(n);
    v2.resize(n);
    v3.resize(n);
}

void MCFab::defineEdgeIntersections (Box const& sdf_box)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Box const edge_box = amrex::enclosedCells(sdf_box, idim);
        m_edge_intersections[idim].resize(edge_box, 1);
        m_edge_intersections[idim].setVal<RunOn::Device>(invalid_edge_intersection);
    }
}

void marching_cubes (Geometry const& geom, FArrayBox& sdf_fab, MCFab& mc_fab, int* counters)
{
    BL_PROFILE("marching_cubes");

    // The prefix sums below index vertices (up to 3 per node) and triangles
    // (up to 12 per cell) with int, so bound the node count accordingly.
    AMREX_ALWAYS_ASSERT(sdf_fab.numPts() < Long(std::numeric_limits<int>::max())/12);

    // Exact zeros belong to the covered side. This lets the cleanup loop move
    // nodes to ON without introducing a new positive fluid sample.
    auto const& sdf = sdf_fab.array();

    Box const nbox = sdf_fab.box();
    Box cbox = amrex::enclosedCells(nbox);
    Box exbox = amrex::enclosedCells(nbox, 0);
    Box eybox = amrex::enclosedCells(nbox, 1);
    Box ezbox = amrex::enclosedCells(nbox, 2);

    BaseFab<int> ex_fab(exbox,2);
    BaseFab<int> ey_fab(eybox,2);
    BaseFab<int> ez_fab(ezbox,2);
    Array4<int> ex = ex_fab.array();
    Array4<int> ey = ey_fab.array();
    Array4<int> ez = ez_fab.array();
    if (mc_fab.m_edge_intersections[0].nComp() == 0
        || mc_fab.m_edge_intersections[1].nComp() == 0
        || mc_fab.m_edge_intersections[2].nComp() == 0
        || !mc_fab.m_edge_intersections[0].box().contains(exbox)
        || !mc_fab.m_edge_intersections[1].box().contains(eybox)
        || !mc_fab.m_edge_intersections[2].box().contains(ezbox))
    {
        // Standalone MC callers do not have an original STL to refine
        // against.  Keep their established sampled-SDF interpolation.
        mc_fab.defineEdgeIntersections(sdf_fab.box());
    }
    auto const exact_x = mc_fab.m_edge_intersections[0].const_array();
    auto const exact_y = mc_fab.m_edge_intersections[1].const_array();
    auto const exact_z = mc_fab.m_edge_intersections[2].const_array();
    BoxIndexer n_bi(nbox);

    auto nvx = Scan::PrefixSum<int>(int(nbox.numPts()),
                                    [=] AMREX_GPU_DEVICE (int m) {
                                        auto [i,j,k] = n_bi(m);
                                        int vx = 0, vy = 0, vz = 0;
                                        if (ex.contains(i,j,k)) {
                                            if ((sdf(i,j,k) > 0)
                                                != (sdf(i+1,j,k) > 0))
                                            {
                                                vx = 1;
                                            }
                                            ex(i,j,k,0) = vx;
                                        }
                                        if (ey.contains(i,j,k)) {
                                            if ((sdf(i,j,k) > 0)
                                                != (sdf(i,j+1,k) > 0))
                                            {
                                                vy = 1;
                                            }
                                            ey(i,j,k,0) = vy;
                                        }
                                        if (ez.contains(i,j,k)) {
                                            if ((sdf(i,j,k) > 0)
                                                != (sdf(i,j,k+1) > 0))
                                            {
                                                vz = 1;
                                            }
                                            ez(i,j,k,0) = vz;
                                        }
                                        return vx + vy + vz;
                                    },
                                    [=] AMREX_GPU_DEVICE (int m, int ps) {
                                        auto [i,j,k] = n_bi(m);
                                        // Component 1 holds the vertex index of the
                                        // crossing on this edge, or -1 when the edge has
                                        // no crossing.  The -1 sentinel is what
                                        // add_c_vertex and the triangle assembly test.
                                        if (ex.contains(i,j,k)) {
                                            if (ex(i,j,k,0)) { ex(i,j,k,1) = ps++; }
                                            else             { ex(i,j,k,1) = -1;   }
                                        }
                                        if (ey.contains(i,j,k)) {
                                            if (ey(i,j,k,0)) { ey(i,j,k,1) = ps++; }
                                            else             { ey(i,j,k,1) = -1;   }
                                        }
                                        if (ez.contains(i,j,k)) {
                                            if (ez(i,j,k,0)) { ez(i,j,k,1) = ps;   }
                                            else             { ez(i,j,k,1) = -1;   }
                                        }
                                    },
                                    Scan::Type::exclusive, Scan::retSum);

    Vertex vrtx;
    vrtx.resize(nvx);

    auto pvrtx = vrtx.dataPtrs();
    ParallelFor(nbox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (ex.contains(i,j,k) && ex(i,j,k,0)) {
            int m = ex(i,j,k,1);
            Real u = edge_intersection_fraction(
                sdf(i,j,k), sdf(i+1,j,k), exact_x(i,j,k));
            pvrtx[0][m] = Real(i) + u;
            pvrtx[1][m] = Real(j);
            pvrtx[2][m] = Real(k);
        }
        if (ey.contains(i,j,k) && ey(i,j,k,0)) {
            int m = ey(i,j,k,1);
            Real u = edge_intersection_fraction(
                sdf(i,j,k), sdf(i,j+1,k), exact_y(i,j,k));
            pvrtx[0][m] = Real(i);
            pvrtx[1][m] = Real(j) + u;
            pvrtx[2][m] = Real(k);
        }
        if (ez.contains(i,j,k) && ez(i,j,k,0)) {
            int m = ez(i,j,k,1);
            Real u = edge_intersection_fraction(
                sdf(i,j,k), sdf(i,j,k+1), exact_z(i,j,k));
            pvrtx[0][m] = Real(i);
            pvrtx[1][m] = Real(j);
            pvrtx[2][m] = Real(k) + u;
        }
    });

    LookUpTable const* lut = d_table;
    auto const& sdf_c = sdf_fab.const_array();
    auto const& ex_c = ex_fab.const_array();
    auto const& ey_c = ey_fab.const_array();
    auto const& ez_c = ez_fab.const_array();

    BaseFab<int> ntri_fab(cbox,num_cell_data_components);
    ntri_fab.setVal<RunOn::Device>(0);
    auto const& ntri = ntri_fab.array();
    GpuArray<int*,3> ptri{nullptr,nullptr,nullptr};

    int* const perror = counters + counter_invalid_triangles;

    BoxIndexer c_bi(cbox);

    int nttot = Scan::PrefixSum<int>(int(cbox.numPts()),
                                     [=] AMREX_GPU_DEVICE (int m) {
                                         auto [i,j,k] = c_bi(m);
                                         int err = 0;
                                         ntri(i,j,k,interior_vertex_count) = 0;
                                         process_cube(0,lut,i,j,k,sdf_c,ex_c,ey_c,ez_c,pvrtx,ntri,ptri,&err);
                                         if (err != 0) {
                                             Gpu::Atomic::AddNoRet(perror, 1);
                                         }
                                         return ntri(i,j,k,triangle_count);
                                     },
                                     [=] AMREX_GPU_DEVICE (int m, int ps) {
                                         auto [i,j,k] = c_bi(m);
                                         ntri(i,j,k,triangle_offset) = ps;
                                     },
                                     Scan::Type::exclusive, Scan::retSum);

    int const edge_vertex_count = nvx;
    int nvx_c = Scan::PrefixSum<int>(int(cbox.numPts()),
                                     [=] AMREX_GPU_DEVICE (int m) {
                                         auto [i,j,k] = c_bi(m);
                                         return ntri(i,j,k,interior_vertex_count);
                                     },
                                     [=] AMREX_GPU_DEVICE (int m, int ps) {
                                         auto [i,j,k] = c_bi(m);
                                         ntri(i,j,k,interior_vertex_offset) =
                                             edge_vertex_count + ps;
                                     },
                                     Scan::Type::exclusive, Scan::retSum);

    if (nvx_c > 0) {
        nvx += nvx_c;
        vrtx.resize(nvx);
        pvrtx = vrtx.dataPtrs();
    }

    Triangle tri;
    tri.resize(nttot);
    ptri = tri.dataPtrs();

    ParallelFor(cbox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        int err = 0;
        process_cube(1,lut,i,j,k,sdf_c,ex_c,ey_c,ez_c,pvrtx,ntri,ptri,&err);
        if (err != 0) {
            Gpu::Atomic::AddNoRet(perror, 1);
        }
    });

    // Shift vertices
    auto problo = geom.ProbLoArray();
    auto dx = geom.CellSizeArray();
    ParallelFor(nvx, [=] AMREX_GPU_DEVICE (int m)
    {
        pvrtx[0][m] = problo[0] + dx[0] * pvrtx[0][m];
        pvrtx[1][m] = problo[1] + dx[1] * pvrtx[1][m];
        pvrtx[2][m] = problo[2] + dx[2] * pvrtx[2][m];
    });

    // The scratch edge fabs are destroyed on return, so wait for the kernels
    // that read them.  Both passes accumulate into the same counter: an
    // unknown MC33 face id in pass 0 or a triangle referencing an edge without
    // a crossing in pass 1 are lookup-table invariants that the driver checks.
    Gpu::streamSynchronize();

    mc_fab.m_cell_data = std::move(ntri_fab);
    mc_fab.m_triangles = std::move(tri);
    mc_fab.m_vertices = std::move(vrtx);
}

void build_face_fractions (
    Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf_fab,
    FArrayBox& apx_fab, FArrayBox& apy_fab, FArrayBox& apz_fab,
    FArrayBox& fcx_fab, FArrayBox& fcy_fab, FArrayBox& fcz_fab,
    IArrayBox& rejected_x_fab, IArrayBox& rejected_y_fab,
    IArrayBox& rejected_z_fab, int* counters)
{
    BL_PROFILE("MC::build_face_fractions");

    AMREX_ALWAYS_ASSERT(mc_fab.m_cell_data.box().contains(bx));
    auto const rejected_x = rejected_x_fab.array();
    auto const rejected_y = rejected_y_fab.array();
    auto const rejected_z = rejected_z_fab.array();
    Box const rejected_xbox = rejected_x_fab.box();
    Box const rejected_ybox = rejected_y_fab.box();
    Box const rejected_zbox = rejected_z_fab.box();

    auto const sdf = sdf_fab.const_array();
    auto const cell_data = mc_fab.m_cell_data.const_array();
    auto const exact_x = mc_fab.m_edge_intersections[0].const_array();
    auto const exact_y = mc_fab.m_edge_intersections[1].const_array();
    auto const exact_z = mc_fab.m_edge_intersections[2].const_array();
    Box const cell_box = mc_fab.m_cell_data.box();
    auto const apx = apx_fab.array();
    auto const apy = apy_fab.array();
    auto const apz = apz_fab.array();
    auto const fcx = fcx_fab.array();
    auto const fcy = fcy_fab.array();
    auto const fcz = fcz_fab.array();

    // The two cells sharing a face resolved its MC33 ambiguity differently
    // (an invariant violation, fatal in the driver), and degenerate face
    // polygons that were marked in the rejected face arrays for nodal repair.
    int* const error = counters + counter_face_decision_errors;
    int* const degenerate = counters + counter_degenerate_faces;

    // Chombo's moment construction starts with boundary-face moments.  Build
    // those apertures from the same signed-distance edge intersections used
    // by marching cubes, so the six Cartesian patches and EB triangles close.
    // Edges traversed against their storage direction use 1 - exact; the
    // invalid_edge_intersection sentinel (-1) maps to 2, which
    // edge_intersection_fraction also treats as "no exact crossing".
    ParallelFor(amrex::surroundingNodes(bx,0),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real const face_levelset[4] = {
            sdf(i,j  ,k  ), sdf(i,j+1,k  ),
            sdf(i,j+1,k+1), sdf(i,j  ,k+1)
        };
        Real const intersections[4] = {
            exact_y(i,j,k), exact_z(i,j+1,k),
            1.0_rt-exact_y(i,j,k+1), 1.0_rt-exact_z(i,j,k)
        };
        auto const decision = resolved_face_decision(
            cell_data, cell_box, error, i-1,j,k,1, i,j,k,3);
        bool const ok = cut_face_fraction(face_levelset, intersections, apx(i,j,k),
                                          fcx(i,j,k,0), fcx(i,j,k,1),
                                          decision[0] != 0, decision[1] != 0);
        if (!ok && rejected_xbox.contains(i,j,k)) {
            rejected_x(i,j,k) = 1;
            Gpu::Atomic::AddNoRet(degenerate, 1);
        }
    });

    ParallelFor(amrex::surroundingNodes(bx,1),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real const face_levelset[4] = {
            sdf(i  ,j,k  ), sdf(i+1,j,k  ),
            sdf(i+1,j,k+1), sdf(i  ,j,k+1)
        };
        Real const intersections[4] = {
            exact_x(i,j,k), exact_z(i+1,j,k),
            1.0_rt-exact_x(i,j,k+1), 1.0_rt-exact_z(i,j,k)
        };
        auto const decision = resolved_face_decision(
            cell_data, cell_box, error, i,j-1,k,2, i,j,k,0);
        bool const ok = cut_face_fraction(face_levelset, intersections, apy(i,j,k),
                                          fcy(i,j,k,0), fcy(i,j,k,1),
                                          decision[0] != 0, decision[1] != 0);
        if (!ok && rejected_ybox.contains(i,j,k)) {
            rejected_y(i,j,k) = 1;
            Gpu::Atomic::AddNoRet(degenerate, 1);
        }
    });

    ParallelFor(amrex::surroundingNodes(bx,2),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real const face_levelset[4] = {
            sdf(i  ,j  ,k), sdf(i+1,j  ,k),
            sdf(i+1,j+1,k), sdf(i  ,j+1,k)
        };
        Real const intersections[4] = {
            exact_x(i,j,k), exact_y(i+1,j,k),
            1.0_rt-exact_x(i,j+1,k), 1.0_rt-exact_y(i,j,k)
        };
        auto const decision = resolved_face_decision(
            cell_data, cell_box, error, i,j,k-1,5, i,j,k,4);
        bool const ok = cut_face_fraction(face_levelset, intersections, apz(i,j,k),
                                          fcz(i,j,k,0), fcz(i,j,k,1),
                                          decision[0] != 0, decision[1] != 0);
        if (!ok && rejected_zbox.contains(i,j,k)) {
            rejected_z(i,j,k) = 1;
            Gpu::Atomic::AddNoRet(degenerate, 1);
        }
    });

}

void build_edge_centroids (
    Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf_fab,
    FArrayBox& ecx_fab, FArrayBox& ecy_fab, FArrayBox& ecz_fab)
{
    BL_PROFILE("MC::build_edge_centroids");

    Box const exbx = amrex::convert(bx,IntVect(0,1,1));
    Box const eybx = amrex::convert(bx,IntVect(1,0,1));
    Box const ezbx = amrex::convert(bx,IntVect(1,1,0));
    AMREX_ALWAYS_ASSERT(ecx_fab.box().contains(exbx));
    AMREX_ALWAYS_ASSERT(ecy_fab.box().contains(eybx));
    AMREX_ALWAYS_ASSERT(ecz_fab.box().contains(ezbx));

    auto const sdf = sdf_fab.const_array();
    auto const exact_x = mc_fab.m_edge_intersections[0].const_array();
    auto const exact_y = mc_fab.m_edge_intersections[1].const_array();
    auto const exact_z = mc_fab.m_edge_intersections[2].const_array();
    auto const ecx = ecx_fab.array();
    auto const ecy = ecy_fab.array();
    auto const ecz = ecz_fab.array();

    ParallelFor(exbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real const lo = sdf(i,j,k);
        Real const hi = sdf(i+1,j,k);
        bool const lo_fluid = lo > 0.0_rt;
        bool const hi_fluid = hi > 0.0_rt;
        if (lo_fluid && hi_fluid) {
            ecx(i,j,k) = 1.0_rt;
        } else if (!lo_fluid && !hi_fluid) {
            ecx(i,j,k) = -1.0_rt;
        } else {
            Real const cut = edge_intersection_fraction(
                lo, hi, exact_x(i,j,k));
            ecx(i,j,k) = lo_fluid
                ? 0.5_rt*cut-0.5_rt
                : 0.5_rt*cut;
        }
    });

    ParallelFor(eybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real const lo = sdf(i,j,k);
        Real const hi = sdf(i,j+1,k);
        bool const lo_fluid = lo > 0.0_rt;
        bool const hi_fluid = hi > 0.0_rt;
        if (lo_fluid && hi_fluid) {
            ecy(i,j,k) = 1.0_rt;
        } else if (!lo_fluid && !hi_fluid) {
            ecy(i,j,k) = -1.0_rt;
        } else {
            Real const cut = edge_intersection_fraction(
                lo, hi, exact_y(i,j,k));
            ecy(i,j,k) = lo_fluid
                ? 0.5_rt*cut-0.5_rt
                : 0.5_rt*cut;
        }
    });

    ParallelFor(ezbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real const lo = sdf(i,j,k);
        Real const hi = sdf(i,j,k+1);
        bool const lo_fluid = lo > 0.0_rt;
        bool const hi_fluid = hi > 0.0_rt;
        if (lo_fluid && hi_fluid) {
            ecz(i,j,k) = 1.0_rt;
        } else if (!lo_fluid && !hi_fluid) {
            ecz(i,j,k) = -1.0_rt;
        } else {
            Real const cut = edge_intersection_fraction(
                lo, hi, exact_z(i,j,k));
            ecz(i,j,k) = lo_fluid
                ? 0.5_rt*cut-0.5_rt
                : 0.5_rt*cut;
        }
    });
}

void build_cell_fractions (
    Box const& bx, Geometry const& geom, MCFab const& mc_fab,
    FArrayBox const& sdf_fab,
    FArrayBox& apx_fab, FArrayBox& apy_fab, FArrayBox& apz_fab,
    FArrayBox& vfrac_fab, FArrayBox& vcent_fab, FArrayBox& barea_fab,
    FArrayBox& bcent_fab, FArrayBox& bnorm_fab, int* counters)
{
    BL_PROFILE("MC::build_cell_fractions");

    AMREX_ALWAYS_ASSERT(mc_fab.m_cell_data.box().contains(bx));
    AMREX_ALWAYS_ASSERT(
        sdf_fab.box().contains(amrex::surroundingNodes(bx)));
    auto const dx = geom.CellSizeArray();
    Real const max_dx = amrex::max(dx[0],amrex::max(dx[1],dx[2]));
    Real const cubic_tolerance =
        16.0_rt*std::numeric_limits<Real>::epsilon()*max_dx;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::abs(dx[0]-dx[1]) <= cubic_tolerance
        && std::abs(dx[0]-dx[2]) <= cubic_tolerance,
        "Marching-cubes cut-cell fractions currently require dx == dy == dz");

    auto const cell_data = mc_fab.m_cell_data.const_array();
    auto const sdf = sdf_fab.const_array();
    auto const apx = apx_fab.array();
    auto const apy = apy_fab.array();
    auto const apz = apz_fab.array();
    auto const vfrac = vfrac_fab.array();
    auto const vcent = vcent_fab.array();
    auto const barea = barea_fab.array();
    auto const bcent = bcent_fab.array();
    auto const bnorm = bnorm_fab.array();

    auto const* tri_v1 = mc_fab.m_triangles.v1.data();
    auto const* tri_v2 = mc_fab.m_triangles.v2.data();
    auto const* tri_v3 = mc_fab.m_triangles.v3.data();
    auto const* vert_x = mc_fab.m_vertices.x.data();
    auto const* vert_y = mc_fab.m_vertices.y.data();
    auto const* vert_z = mc_fab.m_vertices.z.data();

    auto const problo = geom.ProbLoArray();
    // Area-vector closure, volume, centroid, and boundary-normal errors in
    // four consecutive counter slots.  Every quantity checked against
    // tolerance is dimensionless.
    static_assert(counter_volume_errors == counter_closure_errors + 1
                  && counter_centroid_errors == counter_closure_errors + 2
                  && counter_area_vector_errors == counter_closure_errors + 3);
    int* const errors = counters + counter_closure_errors;

#ifdef AMREX_USE_FLOAT
    constexpr Real tolerance = 2.e-5_rt;
#else
    constexpr Real tolerance = 2.e-12_rt;
#endif

    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        int const triangle_count =
            cell_data(i,j,k,CellDataComponent::triangle_count);
        int const triangle_offset =
            cell_data(i,j,k,CellDataComponent::triangle_offset);
        if (triangle_count <= 0) {
            return;
        }

        Real const cell_lo[3] = {
            problo[0] + Real(i)*dx[0],
            problo[1] + Real(j)*dx[1],
            problo[2] + Real(k)*dx[2]
        };

        // The coordinate-face patches and the EB triangles form a closed
        // surface.  Work in cell-local coordinates so volume is already a
        // volume fraction and first moments normalize directly to AMReX
        // centroid coordinates.
        Real eb_area_vector[3] = {0.0_rt, 0.0_rt, 0.0_rt};
        Real eb_area = 0.0_rt;
        Real eb_centroid_numerator[3] = {0.0_rt, 0.0_rt, 0.0_rt};

        for (int n = 0; n < triangle_count; ++n) {
            int const m = triangle_offset + n;
            int const iv1 = tri_v1[m];
            int const iv2 = tri_v2[m];
            int const iv3 = tri_v3[m];

            Real const p1[3] = {vert_x[iv1], vert_y[iv1], vert_z[iv1]};
            Real const p2[3] = {vert_x[iv2], vert_y[iv2], vert_z[iv2]};
            Real const p3[3] = {vert_x[iv3], vert_y[iv3], vert_z[iv3]};
            Real q1[3];
            Real q2[3];
            Real q3[3];
            for (int d = 0; d < 3; ++d) {
                q1[d] = (p1[d]-cell_lo[d])/dx[d] - 0.5_rt;
                q2[d] = (p2[d]-cell_lo[d])/dx[d] - 0.5_rt;
                q3[d] = (p3[d]-cell_lo[d])/dx[d] - 0.5_rt;
            }

            Real const q21[3] = {q2[0]-q1[0], q2[1]-q1[1], q2[2]-q1[2]};
            Real const q31[3] = {q3[0]-q1[0], q3[1]-q1[1], q3[2]-q1[2]};
            Real const local_av[3] = {
                0.5_rt*(q21[1]*q31[2]-q21[2]*q31[1]),
                0.5_rt*(q21[2]*q31[0]-q21[0]*q31[2]),
                0.5_rt*(q21[0]*q31[1]-q21[1]*q31[0])
            };
            for (int d = 0; d < 3; ++d) {
                eb_area_vector[d] += local_av[d];
            }
            Real const triangle_area = std::sqrt(
                local_av[0]*local_av[0]
                + local_av[1]*local_av[1]
                + local_av[2]*local_av[2]);
            eb_area += triangle_area;
            for (int d = 0; d < 3; ++d) {
                eb_centroid_numerator[d] +=
                    triangle_area*(q1[d]+q2[d]+q3[d])/3.0_rt;
            }
        }

        Real const expected_area_vector[3] = {
            apx(i,j,k)-apx(i+1,j,k),
            apy(i,j,k)-apy(i,j+1,k),
            apz(i,j,k)-apz(i,j,k+1)
        };
        Real const orientation_dot =
            eb_area_vector[0]*expected_area_vector[0]
            + eb_area_vector[1]*expected_area_vector[1]
            + eb_area_vector[2]*expected_area_vector[2];
        Real const orientation = (orientation_dot < 0.0_rt) ? -1.0_rt : 1.0_rt;
        if (std::abs(orientation*eb_area_vector[0]-expected_area_vector[0]) > tolerance
            || std::abs(orientation*eb_area_vector[1]-expected_area_vector[1]) > tolerance
            || std::abs(orientation*eb_area_vector[2]-expected_area_vector[2]) > tolerance)
        {
            Gpu::Atomic::AddNoRet(errors, 1);
            return;
        }

        Real volume = (
            0.5_rt*(apx(i,j,k)+apx(i+1,j,k)
                    +apy(i,j,k)+apy(i,j+1,k)
                    +apz(i,j,k)+apz(i,j,k+1)))/3.0_rt;
        Real first_moment[3] = {
            (apx(i+1,j,k)-apx(i,j,k))/8.0_rt,
            (apy(i,j+1,k)-apy(i,j,k))/8.0_rt,
            (apz(i,j,k+1)-apz(i,j,k))/8.0_rt
        };

        for (int n = 0; n < triangle_count; ++n) {
            int const m = triangle_offset + n;
            int const iv1 = tri_v1[m];
            int const iv2 = tri_v2[m];
            int const iv3 = tri_v3[m];

            Real const q1[3] = {
                (vert_x[iv1]-cell_lo[0])/dx[0] - 0.5_rt,
                (vert_y[iv1]-cell_lo[1])/dx[1] - 0.5_rt,
                (vert_z[iv1]-cell_lo[2])/dx[2] - 0.5_rt
            };
            Real const q2[3] = {
                (vert_x[iv2]-cell_lo[0])/dx[0] - 0.5_rt,
                (vert_y[iv2]-cell_lo[1])/dx[1] - 0.5_rt,
                (vert_z[iv2]-cell_lo[2])/dx[2] - 0.5_rt
            };
            Real const q3[3] = {
                (vert_x[iv3]-cell_lo[0])/dx[0] - 0.5_rt,
                (vert_y[iv3]-cell_lo[1])/dx[1] - 0.5_rt,
                (vert_z[iv3]-cell_lo[2])/dx[2] - 0.5_rt
            };
            Real const q21[3] = {q2[0]-q1[0], q2[1]-q1[1], q2[2]-q1[2]};
            Real const q31[3] = {q3[0]-q1[0], q3[1]-q1[1], q3[2]-q1[2]};
            Real const area_vector[3] = {
                orientation*0.5_rt*(q21[1]*q31[2]-q21[2]*q31[1]),
                orientation*0.5_rt*(q21[2]*q31[0]-q21[0]*q31[2]),
                orientation*0.5_rt*(q21[0]*q31[1]-q21[1]*q31[0])
            };
            Real const triangle_centroid[3] = {
                (q1[0]+q2[0]+q3[0])/3.0_rt,
                (q1[1]+q2[1]+q3[1])/3.0_rt,
                (q1[2]+q2[2]+q3[2])/3.0_rt
            };
            volume += (triangle_centroid[0]*area_vector[0]
                       + triangle_centroid[1]*area_vector[1]
                       + triangle_centroid[2]*area_vector[2])/3.0_rt;

            for (int d = 0; d < 3; ++d) {
                Real const quadratic =
                    q1[d]*q1[d] + q2[d]*q2[d] + q3[d]*q3[d]
                    + q1[d]*q2[d] + q2[d]*q3[d] + q3[d]*q1[d];
                first_moment[d] += area_vector[d]*quadratic/12.0_rt;
            }
        }

        bool nodal_plane_owner = true;
        bool nodal_plane_covered = true;
        bool has_on_node = false;
        for (int kk = 0; kk <= 1; ++kk) {
            for (int jj = 0; jj <= 1; ++jj) {
                for (int ii = 0; ii <= 1; ++ii) {
                    Real const value = sdf(i+ii,j+jj,k+kk);
                    nodal_plane_owner = nodal_plane_owner
                        && value >= 0.0_rt;
                    nodal_plane_covered = nodal_plane_covered
                        && value <= 0.0_rt;
                    has_on_node = has_on_node || value == 0.0_rt;
                }
            }
        }
        nodal_plane_owner = nodal_plane_owner && has_on_node;
        nodal_plane_covered = nodal_plane_covered && has_on_node;
        if (nodal_plane_covered) {
            // A coincident patch belongs to the fluid-side cell.  On the
            // solid side, the zero-area MC triangles are only a classification
            // artifact and the cell remains fully covered.
            vfrac(i,j,k) = 0.0_rt;
            return;
        }
        if (nodal_plane_owner) {
            if (eb_area <= tolerance) {
                // The surface only touches this cell at a node or edge.  It
                // has no measure inside the cell, so the fluid-side cell is
                // regular rather than a zero-area cut cell.  Write the full
                // set of "no boundary in this cell" values explicitly in case
                // a cut face still makes the cell single-valued downstream.
                vfrac(i,j,k) = 1.0_rt;
                barea(i,j,k) = 0.0_rt;
                for (int d = 0; d < 3; ++d) {
                    vcent(i,j,k,d) = 0.0_rt;
                    bcent(i,j,k,d) = 0.0_rt;
                    bnorm(i,j,k,d) = 0.0_rt;
                }
                return;
            }
            // The EB lies on one or more cell faces and this cell owns the
            // coincident patch. Its open volume is still the complete cell.
            volume = 1.0_rt;
            first_moment[0] = 0.0_rt;
            first_moment[1] = 0.0_rt;
            first_moment[2] = 0.0_rt;
        }

        if (volume <= tolerance || volume > 1.0_rt+tolerance
            || eb_area <= tolerance)
        {
            Gpu::Atomic::AddNoRet(errors+1, 1);
            return;
        }
        volume = amrex::min(amrex::max(volume,0.0_rt),1.0_rt);

        Real centroid[3] = {
            first_moment[0]/volume,
            first_moment[1]/volume,
            first_moment[2]/volume
        };
        if (centroid[0] < -0.5_rt-tolerance || centroid[0] > 0.5_rt+tolerance
            || centroid[1] < -0.5_rt-tolerance || centroid[1] > 0.5_rt+tolerance
            || centroid[2] < -0.5_rt-tolerance || centroid[2] > 0.5_rt+tolerance)
        {
            Gpu::Atomic::AddNoRet(errors+2, 1);
            return;
        }

        Real boundary_normal[3] = {
            orientation*eb_area_vector[0],
            orientation*eb_area_vector[1],
            orientation*eb_area_vector[2]
        };
        Real const area_vector_norm = std::sqrt(
            boundary_normal[0]*boundary_normal[0]
            + boundary_normal[1]*boundary_normal[1]
            + boundary_normal[2]*boundary_normal[2]);
        if (area_vector_norm <= tolerance) {
            Gpu::Atomic::AddNoRet(errors+3, 1);
            return;
        }
        for (Real& normal_component : boundary_normal) {
            normal_component /= area_vector_norm;
        }

        Real boundary_centroid[3] = {
            eb_centroid_numerator[0]/eb_area,
            eb_centroid_numerator[1]/eb_area,
            eb_centroid_numerator[2]/eb_area
        };
        if (boundary_centroid[0] < -0.5_rt-tolerance
            || boundary_centroid[0] > 0.5_rt+tolerance
            || boundary_centroid[1] < -0.5_rt-tolerance
            || boundary_centroid[1] > 0.5_rt+tolerance
            || boundary_centroid[2] < -0.5_rt-tolerance
            || boundary_centroid[2] > 0.5_rt+tolerance)
        {
            Gpu::Atomic::AddNoRet(errors+2, 1);
            return;
        }

        vfrac(i,j,k) = volume;
        for (int d = 0; d < 3; ++d) {
            vcent(i,j,k,d) = amrex::min(amrex::max(centroid[d],-0.5_rt),0.5_rt);
            bcent(i,j,k,d) = amrex::Clamp(
                boundary_centroid[d], -0.5_rt, 0.5_rt);
            bnorm(i,j,k,d) = boundary_normal[d];
        }
        barea(i,j,k) = eb_area;
    });

}

int build_cell_topology (Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf_fab,
                         EBCellFlagFab& cellflag, FArrayBox const& vfrac_fab,
                         FArrayBox const& apx_fab, FArrayBox const& apy_fab,
                         FArrayBox const& apz_fab)
{
    BL_PROFILE("MC::build_cell_topology");

    Box const bxg1 = amrex::grow(bx, 1);
    Box const nbxg1 = amrex::surroundingNodes(bxg1);
    IntVect const x_transverse_ghost =
        IntVect::TheUnitVector() - IntVect::TheDimensionVector(0);
    IntVect const y_transverse_ghost =
        IntVect::TheUnitVector() - IntVect::TheDimensionVector(1);
    IntVect const z_transverse_ghost =
        IntVect::TheUnitVector() - IntVect::TheDimensionVector(2);
    Box const valid_fxbx = amrex::surroundingNodes(bxg1, 0);
    Box const valid_fybx = amrex::surroundingNodes(bxg1, 1);
    Box const valid_fzbx = amrex::surroundingNodes(bxg1, 2);
    Box const fxbx = amrex::grow(valid_fxbx, x_transverse_ghost);
    Box const fybx = amrex::grow(valid_fybx, y_transverse_ghost);
    Box const fzbx = amrex::grow(valid_fzbx, z_transverse_ghost);

    AMREX_ALWAYS_ASSERT(sdf_fab.box().contains(nbxg1));
    AMREX_ALWAYS_ASSERT(cellflag.box().contains(bxg1));
    AMREX_ALWAYS_ASSERT(vfrac_fab.box().contains(bxg1));
    AMREX_ALWAYS_ASSERT(mc_fab.m_cell_data.box().contains(bx));
    AMREX_ALWAYS_ASSERT(apx_fab.box().contains(fxbx));
    AMREX_ALWAYS_ASSERT(apy_fab.box().contains(fybx));
    AMREX_ALWAYS_ASSERT(apz_fab.box().contains(fzbx));

    auto const sdf = sdf_fab.const_array();
    auto const cell = cellflag.array();
    auto const cell_data = mc_fab.m_cell_data.const_array();
    auto const vfrac = vfrac_fab.const_array();
    auto const apx = apx_fab.const_array();
    auto const apy = apy_fab.const_array();
    auto const apz = apz_fab.const_array();

    BaseFab<EB2::Type_t> fx_fab(fxbx);
    BaseFab<EB2::Type_t> fy_fab(fybx);
    BaseFab<EB2::Type_t> fz_fab(fzbx);
    fx_fab.setVal<RunOn::Device>(EB2::Type::regular);
    fy_fab.setVal<RunOn::Device>(EB2::Type::regular);
    fz_fab.setVal<RunOn::Device>(EB2::Type::regular);
    auto const fx = fx_fab.array();
    auto const fy = fy_fab.array();
    auto const fz = fz_fab.array();

#ifdef AMREX_USE_FLOAT
    constexpr Real tolerance = 2.e-5_rt;
#else
    constexpr Real tolerance = 2.e-12_rt;
#endif

    ParallelFor(bxg1, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        bool const all_faces_open =
            apx(i, j, k) >= 1.0_rt - tolerance && apx(i + 1, j, k) >= 1.0_rt - tolerance &&
            apy(i, j, k) >= 1.0_rt - tolerance && apy(i, j + 1, k) >= 1.0_rt - tolerance &&
            apz(i, j, k) >= 1.0_rt - tolerance && apz(i, j, k + 1) >= 1.0_rt - tolerance;

        if (vfrac(i, j, k) <= tolerance) {
            cell(i, j, k).setCovered();
        } else if (vfrac(i, j, k) >= 1.0_rt - tolerance && all_faces_open) {
            cell(i, j, k).setRegular();
        } else {
            cell(i, j, k).setSingleValued();
        }
    });

    ParallelFor(valid_fxbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        fx(i, j, k) = face_type(apx(i, j, k), tolerance);
    });
    ParallelFor(valid_fybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        fy(i, j, k) = face_type(apy(i, j, k), tolerance);
    });
    ParallelFor(valid_fzbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        fz(i, j, k) = face_type(apz(i, j, k), tolerance);
    });

    Gpu::DeviceScalar<int> error_count(0);
    int* const errors = error_count.dataPtr();
    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        int const triangle_count = cell_data(i, j, k, CellDataComponent::triangle_count);
        bool const has_triangles = triangle_count > 0;
        bool const source_has_fluid = sdf(i, j, k) > 0.0_rt || sdf(i + 1, j, k) > 0.0_rt ||
                                      sdf(i, j + 1, k) > 0.0_rt || sdf(i + 1, j + 1, k) > 0.0_rt ||
                                      sdf(i, j, k + 1) > 0.0_rt || sdf(i + 1, j, k + 1) > 0.0_rt ||
                                      sdf(i, j + 1, k + 1) > 0.0_rt ||
                                      sdf(i + 1, j + 1, k + 1) > 0.0_rt;
        bool const source_has_covered =
            sdf(i, j, k) <= 0.0_rt || sdf(i + 1, j, k) <= 0.0_rt || sdf(i, j + 1, k) <= 0.0_rt ||
            sdf(i + 1, j + 1, k) <= 0.0_rt || sdf(i, j, k + 1) <= 0.0_rt ||
            sdf(i + 1, j, k + 1) <= 0.0_rt || sdf(i, j + 1, k + 1) <= 0.0_rt ||
            sdf(i + 1, j + 1, k + 1) <= 0.0_rt;
        bool const source_is_cut = source_has_fluid && source_has_covered;
        if (has_triangles != source_is_cut) {
            Gpu::Atomic::AddNoRet(errors, 1);
        }
    });

    EBCellFlagFab cellflagtmp(cellflag.box(), 1);
    Elixir cellflagtmp_eli = cellflagtmp.elixir();
    EB2::set_connection_flags(bx, bxg1, cell, cellflagtmp.array(), fx, fy, fz);

    Gpu::streamSynchronize();
    return error_count.dataValue();
}

void mark_faces_for_cleanup (Box const& bx, MCFab const& mc_fab, FArrayBox const& sdf_fab,
                             IArrayBox& rejected_x_fab, IArrayBox& rejected_y_fab,
                             IArrayBox& rejected_z_fab, int* counters)
{
    BL_PROFILE("MC::mark_faces_for_cleanup");

    Box const xbx = amrex::surroundingNodes(bx, 0);
    Box const ybx = amrex::surroundingNodes(bx, 1);
    Box const zbx = amrex::surroundingNodes(bx, 2);
    AMREX_ALWAYS_ASSERT(mc_fab.m_cell_data.box().contains(bx));
    AMREX_ALWAYS_ASSERT(sdf_fab.box().contains(amrex::surroundingNodes(bx)));
    AMREX_ALWAYS_ASSERT(rejected_x_fab.box().contains(xbx));
    AMREX_ALWAYS_ASSERT(rejected_y_fab.box().contains(ybx));
    AMREX_ALWAYS_ASSERT(rejected_z_fab.box().contains(zbx));

    auto const cell_data = mc_fab.m_cell_data.const_array();
    Box const cell_box = mc_fab.m_cell_data.box();
    auto const sdf = sdf_fab.const_array();
    auto const rejected_x = rejected_x_fab.array();
    auto const rejected_y = rejected_y_fab.array();
    auto const rejected_z = rejected_z_fab.array();

    int* const count = counters + counter_face_rejections;

    ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        bool const rejected =
            face_is_rejected(sdf(i, j, k), sdf(i, j + 1, k), sdf(i, j + 1, k + 1), sdf(i, j, k + 1),
                             cell_data, cell_box, i - 1, j, k, 1, i, j, k, 3);
        if (rejected && rejected_x(i, j, k) == 0) {
            rejected_x(i, j, k) = 1;
            Gpu::Atomic::AddNoRet(count, 1);
        }
    });
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        bool const rejected =
            face_is_rejected(sdf(i, j, k), sdf(i + 1, j, k), sdf(i + 1, j, k + 1), sdf(i, j, k + 1),
                             cell_data, cell_box, i, j - 1, k, 2, i, j, k, 0);
        if (rejected && rejected_y(i, j, k) == 0) {
            rejected_y(i, j, k) = 1;
            Gpu::Atomic::AddNoRet(count, 1);
        }
    });
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        bool const rejected =
            face_is_rejected(sdf(i, j, k), sdf(i + 1, j, k), sdf(i + 1, j + 1, k), sdf(i, j + 1, k),
                             cell_data, cell_box, i, j, k - 1, 5, i, j, k, 4);
        if (rejected && rejected_z(i, j, k) == 0) {
            rejected_z(i, j, k) = 1;
            Gpu::Atomic::AddNoRet(count, 1);
        }
    });

}

void zero_nodes_for_cleanup (Box const& node_box, IArrayBox const& rejected_cells_fab,
                             IArrayBox const& rejected_x_fab, IArrayBox const& rejected_y_fab,
                             IArrayBox const& rejected_z_fab, FArrayBox& sdf_fab, int* counters)
{
    BL_PROFILE("MC::zero_nodes_for_cleanup");

    AMREX_ALWAYS_ASSERT(sdf_fab.box().contains(node_box));
    AMREX_ALWAYS_ASSERT(rejected_cells_fab.box().contains(
        amrex::enclosedCells(amrex::grow(node_box, 1))));
    AMREX_ALWAYS_ASSERT(rejected_x_fab.box().contains(amrex::convert(
        amrex::grow(node_box, IntVect(0, 1, 1)), IntVect::TheDimensionVector(0))));
    AMREX_ALWAYS_ASSERT(rejected_y_fab.box().contains(amrex::convert(
        amrex::grow(node_box, IntVect(1, 0, 1)), IntVect::TheDimensionVector(1))));
    AMREX_ALWAYS_ASSERT(rejected_z_fab.box().contains(amrex::convert(
        amrex::grow(node_box, IntVect(1, 1, 0)), IntVect::TheDimensionVector(2))));
    auto const rejected_cells = rejected_cells_fab.const_array();
    auto const rejected_x = rejected_x_fab.const_array();
    auto const rejected_y = rejected_y_fab.const_array();
    auto const rejected_z = rejected_z_fab.const_array();
    auto const sdf = sdf_fab.array();

    int* const changed = counters + counter_changed_nodes;
    ParallelFor(node_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if (sdf(i, j, k) <= 0.0_rt) {
            return;
        }

        bool rejected = false;
        for (int kk = 0; kk <= 1 && !rejected; ++kk) {
            for (int jj = 0; jj <= 1 && !rejected; ++jj) {
                for (int ii = 0; ii <= 1; ++ii) {
                    rejected = rejected || rejected_cells(i - ii, j - jj, k - kk) != 0;
                }
            }
        }
        for (int kk = 0; kk <= 1 && !rejected; ++kk) {
            for (int jj = 0; jj <= 1; ++jj) {
                rejected = rejected || rejected_x(i, j - jj, k - kk) != 0;
            }
        }
        for (int kk = 0; kk <= 1 && !rejected; ++kk) {
            for (int ii = 0; ii <= 1; ++ii) {
                rejected = rejected || rejected_y(i - ii, j, k - kk) != 0;
            }
        }
        for (int jj = 0; jj <= 1 && !rejected; ++jj) {
            for (int ii = 0; ii <= 1; ++ii) {
                rejected = rejected || rejected_z(i - ii, j - jj, k) != 0;
            }
        }
        if (rejected) {
            sdf(i, j, k) = 0.0_rt;
            Gpu::Atomic::AddNoRet(changed, 1);
        }
    });
}

void extend_domain_face_levelset (Box const& node_box, Box const& domain,
                                  GpuArray<int, 3> const& is_periodic, FArrayBox& sdf_fab,
                                  int* counters)
{
    BL_PROFILE("MC::extend_domain_face_levelset");

    AMREX_ALWAYS_ASSERT(sdf_fab.box().contains(node_box));
    Box const nodal_domain = amrex::surroundingNodes(domain);
    Box reference_box = node_box;
    for (int d = 0; d < 3; ++d) {
        if (!is_periodic[d]) {
            reference_box.setSmall(
                d, amrex::Clamp(node_box.smallEnd(d), nodal_domain.smallEnd(d),
                                nodal_domain.bigEnd(d)));
            reference_box.setBig(
                d, amrex::Clamp(node_box.bigEnd(d), nodal_domain.smallEnd(d),
                                nodal_domain.bigEnd(d)));
        }
    }
    AMREX_ALWAYS_ASSERT(sdf_fab.box().contains(reference_box));
    auto const sdf = sdf_fab.array();

    int const domlo_x = nodal_domain.smallEnd(0);
    int const domlo_y = nodal_domain.smallEnd(1);
    int const domlo_z = nodal_domain.smallEnd(2);
    int const domhi_x = nodal_domain.bigEnd(0);
    int const domhi_y = nodal_domain.bigEnd(1);
    int const domhi_z = nodal_domain.bigEnd(2);

    int* const changed = counters + counter_extended_nodes;
    ParallelFor(node_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        int const ii = is_periodic[0] ? i : amrex::Clamp(i, domlo_x, domhi_x);
        int const jj = is_periodic[1] ? j : amrex::Clamp(j, domlo_y, domhi_y);
        int const kk = is_periodic[2] ? k : amrex::Clamp(k, domlo_z, domhi_z);
        int const outside_directions = (ii != i) + (jj != j) + (kk != k);
        if (outside_directions != 1) {
            return;
        }
        Real const extended_value = sdf(ii, jj, kk);
        // A repaired (exact-zero) exterior node stays covered: re-extending a
        // fluid value onto it would undo the repair and could keep the repair
        // loop from converging.
        if (sdf(i, j, k) == 0.0_rt && extended_value > 0.0_rt) {
            return;
        }
        if (sdf(i, j, k) != extended_value) {
            sdf(i, j, k) = extended_value;
            Gpu::Atomic::AddNoRet(changed, 1);
        }
    });
}

void extend_domain_face_edge_intersections (Box const& domain,
                                            GpuArray<int, 3> const& is_periodic,
                                            MCFab& mc_fab)
{
    BL_PROFILE("MC::extend_domain_face_edge_intersections");

    GpuArray<int, 3> const domain_lo{
        domain.smallEnd(0), domain.smallEnd(1), domain.smallEnd(2)};
    GpuArray<int, 3> const domain_hi{
        domain.bigEnd(0), domain.bigEnd(1), domain.bigEnd(2)};
    for (int edge_direction = 0; edge_direction < 3; ++edge_direction) {
        auto const crossing = mc_fab.m_edge_intersections[edge_direction].array();
        Box const edge_box = mc_fab.m_edge_intersections[edge_direction].box();
        ParallelFor(edge_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int index[3] = {i, j, k};
            int reference[3] = {i, j, k};
            int outside_direction = -1;
            int outside_directions = 0;
            for (int d = 0; d < 3; ++d) {
                if (is_periodic[d]) {
                    continue;
                }
                int const lo = domain_lo[d];
                int const hi = domain_hi[d] + (d != edge_direction);
                reference[d] = amrex::Clamp(index[d], lo, hi);
                if (reference[d] != index[d]) {
                    outside_direction = d;
                    ++outside_directions;
                }
            }
            if (outside_directions == 1 && outside_direction != edge_direction) {
                crossing(i, j, k) = crossing(reference[0], reference[1], reference[2]);
            }
        });
    }
}

void mark_cells_for_cleanup (Box const& bx, MCFab const& mc_fab,
                             FArrayBox const& sdf_fab, FArrayBox const& vfrac_fab,
                             Real small_volfrac, IArrayBox& rejected_fab, int* counters)
{
    BL_PROFILE("MC::mark_cells_for_cleanup");

    AMREX_ALWAYS_ASSERT(mc_fab.m_cell_data.box().contains(bx));
    AMREX_ALWAYS_ASSERT(vfrac_fab.box().contains(bx));
    AMREX_ALWAYS_ASSERT(rejected_fab.box().contains(bx));
    AMREX_ALWAYS_ASSERT(sdf_fab.box().contains(amrex::surroundingNodes(bx)));

    auto const cell_data = mc_fab.m_cell_data.const_array();
    auto const sdf = sdf_fab.const_array();
    auto const vfrac = vfrac_fab.const_array();
    auto const rejected = rejected_fab.array();

    static_assert(counter_small_cell_rejections == counter_topology_rejections + 1);
    int* const counts = counters + counter_topology_rejections;

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        Real const cube[8] = {
            sdf(i, j, k),     sdf(i + 1, j, k),     sdf(i + 1, j + 1, k),     sdf(i, j + 1, k),
            sdf(i, j, k + 1), sdf(i + 1, j, k + 1), sdf(i + 1, j + 1, k + 1), sdf(i, j + 1, k + 1)};
        bool fluid[8];
        int nfluid = 0;
        for (int n = 0; n < 8; ++n) {
            fluid[n] = cube[n] > 0.0_rt;
            nfluid += fluid[n];
        }

        int const triangle_count = cell_data(i, j, k, CellDataComponent::triangle_count);
        bool const is_cut = nfluid > 0 && nfluid < 8;
        bool bad_topology = is_cut && triangle_count <= 0;

        // Count the fluid corner groups of the cell.  Corners are joined along
        // the 12 cube edges and across every ambiguous (four-crossing) face
        // whose MC33 decision says the two diagonal fluid corners are
        // connected: that decision is what the extracted surface and the face
        // apertures already use, so a group counted here is exactly a fluid
        // region that touches the cell faces.  Corner groups that MC33 joins
        // only through the cell interior (the tunnel tilings 4.1.2/10.1.2 and
        // the 7.3/10.2/12.2/13.x variants) are deliberately NOT merged: a
        // single-valued EB cell may hold one face-connected fluid region only,
        // so such cells are rejected and repaired.
        int corner_parent[8];
        for (int n = 0; n < 8; ++n) {
            corner_parent[n] = n;
        }
        auto find_root = [&] (int n) {
            while (corner_parent[n] != n) { n = corner_parent[n]; }
            return n;
        };
        auto join = [&] (int a, int b) {
            int const ra = find_root(a);
            int const rb = find_root(b);
            if (ra != rb) { corner_parent[rb] = ra; }
        };
        constexpr int edge_lo[12] = {0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3};
        constexpr int edge_hi[12] = {1, 2, 3, 0, 5, 6, 7, 4, 4, 5, 6, 7};
        for (int n = 0; n < 12; ++n) {
            if (fluid[edge_lo[n]] && fluid[edge_hi[n]]) {
                join(edge_lo[n], edge_hi[n]);
            }
        }
        // MC33 face ids 1..6 in cube-corner order, matching test_face()'s
        // A,B,C,D and the bits stored in m_cell_data.
        constexpr int face_corner[6][4] = {{0, 4, 5, 1}, {1, 5, 6, 2}, {2, 6, 7, 3},
                                           {3, 7, 4, 0}, {0, 3, 2, 1}, {4, 7, 6, 5}};
        int const valid_mask = cell_data(i, j, k, CellDataComponent::face_decision_valid_mask);
        int const connected_mask =
            cell_data(i, j, k, CellDataComponent::face_fluid_connected_mask);
        for (int f = 0; f < 6; ++f) {
            int const c0 = face_corner[f][0];
            int const c1 = face_corner[f][1];
            int const c2 = face_corner[f][2];
            int const c3 = face_corner[f][3];
            if (!(fluid[c0] == fluid[c2] && fluid[c1] == fluid[c3] && fluid[c0] != fluid[c1])) {
                continue; // not an ambiguous face
            }
            int const bit = 1 << f;
            Real const face_levelset[4] = {cube[c0], cube[c1], cube[c2], cube[c3]};
            bool const connected = ((valid_mask & bit) != 0)
                ? ((connected_mask & bit) != 0)
                : four_crossing_fluid_is_connected(face_levelset);
            if (connected) {
                if (fluid[c0]) { join(c0, c2); } else { join(c1, c3); }
            }
        }
        int fluid_components = 0;
        for (int n = 0; n < 8; ++n) {
            if (fluid[n]) {
                fluid_components += find_root(n) == n;
            }
        }
        bad_topology = bad_topology || fluid_components > 1;

        // A negative sentinel means geometry construction rejected the closed
        // boundary. Cover it and let the next MC pass rebuild its neighbors.
        bad_topology = bad_topology || (is_cut && vfrac(i, j, k) < 0.0_rt);
        bool const small_cell = is_cut && !bad_topology && vfrac(i, j, k) < small_volfrac;

        rejected(i, j, k) = bad_topology ? RejectionReason::invalid_topology
                                         : (small_cell ? RejectionReason::small_volume : 0);
        if (bad_topology) {
            Gpu::Atomic::AddNoRet(counts, 1);
        } else if (small_cell) {
            Gpu::Atomic::AddNoRet(counts + 1, 1);
        }
    });

}

void write_stl (std::string const& filename, LayoutData<MCFab> const& mc_fabs)
{
    BoxArray const& grids = mc_fabs.boxArray();
    // The EB may be built inside a ParallelContext sub-frame, so the token
    // chain runs over the sub-communicator like the rest of the builder.
    int const myproc = ParallelContext::MyProcSub();
    int const nprocs = ParallelContext::NProcsSub();

    std::ofstream ofs;

    if (myproc == 0) {
        ofs.open(filename);
        ofs << "solid Created by AMReX\n";
    }

#ifdef AMREX_USE_MPI
    if (myproc > 0) {
        int foo = 0;
        ParallelDescriptor::Recv(&foo, 1, myproc - 1, 100, ParallelContext::CommunicatorSub());
    }
#endif

    if (!ofs.is_open()) {
        ofs.open(filename, std::ios_base::app);
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ofs.good(),
                                     "Could not open marching-cubes STL output " + filename);
    ofs << std::setprecision(std::numeric_limits<Real>::max_digits10);

    for (MFIter mfi(mc_fabs); mfi.isValid(); ++mfi) {
        int const k = mfi.index();
        MCFab const* p = &mc_fabs[mfi];
        AMREX_ALWAYS_ASSERT(p->m_cell_data.box().contains(grids[k]));
        auto ntri = int(p->m_triangles.v1.size());
        amrex::ignore_unused(ntri);

#ifdef AMREX_USE_GPU
        Gpu::PinnedVector<int> tri_v1(ntri);
        Gpu::PinnedVector<int> tri_v2(ntri);
        Gpu::PinnedVector<int> tri_v3(ntri);
        auto nvert = p->m_vertices.x.size();
        Gpu::PinnedVector<Real> vert_x(nvert);
        Gpu::PinnedVector<Real> vert_y(nvert);
        Gpu::PinnedVector<Real> vert_z(nvert);
        BaseFab<int> cell_data(p->m_cell_data.box(), p->m_cell_data.nComp(), The_Pinned_Arena());
        Gpu::copyAsync(Gpu::deviceToHost, p->m_triangles.v1.begin(), p->m_triangles.v1.end(),
                       tri_v1.begin());
        Gpu::copyAsync(Gpu::deviceToHost, p->m_triangles.v2.begin(), p->m_triangles.v2.end(),
                       tri_v2.begin());
        Gpu::copyAsync(Gpu::deviceToHost, p->m_triangles.v3.begin(), p->m_triangles.v3.end(),
                       tri_v3.begin());
        Gpu::copyAsync(Gpu::deviceToHost, p->m_vertices.x.begin(), p->m_vertices.x.end(),
                       vert_x.begin());
        Gpu::copyAsync(Gpu::deviceToHost, p->m_vertices.y.begin(), p->m_vertices.y.end(),
                       vert_y.begin());
        Gpu::copyAsync(Gpu::deviceToHost, p->m_vertices.z.begin(), p->m_vertices.z.end(),
                       vert_z.begin());
        Gpu::dtoh_memcpy_async(cell_data.dataPtr(), p->m_cell_data.dataPtr(),
                               p->m_cell_data.nBytes(p->m_cell_data.box(), p->m_cell_data.nComp()));
        Gpu::streamSynchronize();
#else
        auto const& tri_v1 = p->m_triangles.v1;
        auto const& tri_v2 = p->m_triangles.v2;
        auto const& tri_v3 = p->m_triangles.v3;
        auto const& vert_x = p->m_vertices.x;
        auto const& vert_y = p->m_vertices.y;
        auto const& vert_z = p->m_vertices.z;
        auto const& cell_data = p->m_cell_data;
#endif
        auto const cell = cell_data.const_array();
        amrex::LoopOnCpu(grids[k], [&] (int i, int j, int kk) noexcept {
            int const count = cell(i, j, kk, 0);
            int const offset = cell(i, j, kk, 1);
            for (int n = 0; n < count; ++n) {
                int const itri = offset + n;
                AMREX_ASSERT(itri >= 0 && itri < ntri);
                auto iv1 = tri_v1[itri];
                auto iv2 = tri_v2[itri];
                auto iv3 = tri_v3[itri];
                XDim3 v1{.x = vert_x[iv1], .y = vert_y[iv1], .z = vert_z[iv1]};
                XDim3 v2{.x = vert_x[iv2], .y = vert_y[iv2], .z = vert_z[iv2]};
                XDim3 v3{.x = vert_x[iv3], .y = vert_y[iv3], .z = vert_z[iv3]};
                XDim3 vec1{.x = v2.x - v1.x, .y = v2.y - v1.y, .z = v2.z - v1.z};
                XDim3 vec2{.x = v3.x - v2.x, .y = v3.y - v2.y, .z = v3.z - v2.z};
                XDim3 norm{.x = vec1.y * vec2.z - vec1.z * vec2.y,
                           .y = vec1.z * vec2.x - vec1.x * vec2.z,
                           .z = vec1.x * vec2.y - vec1.y * vec2.x};
                auto tmp = std::sqrt(norm.x * norm.x + norm.y * norm.y + norm.z * norm.z);
                Real const edge_scale_sq =
                    amrex::max(vec1.x * vec1.x + vec1.y * vec1.y + vec1.z * vec1.z,
                               vec2.x * vec2.x + vec2.y * vec2.y + vec2.z * vec2.z);
                Real const degenerate_tolerance =
                    64.0_rt * std::numeric_limits<Real>::epsilon() * edge_scale_sq;
                if (tmp <= degenerate_tolerance) {
                    continue;
                }
                tmp = Real(1) / tmp;
                ofs << "facet normal " << norm.x * tmp << " " << norm.y * tmp << " " << norm.z * tmp
                    << "\n"
                    << "  outer loop\n"
                    << "    vertex " << v1.x << " " << v1.y << " " << v1.z << "\n"
                    << "    vertex " << v2.x << " " << v2.y << " " << v2.z << "\n"
                    << "    vertex " << v3.x << " " << v3.y << " " << v3.z << "\n"
                    << "  endloop\n"
                    << "endfacet\n";
            }
        });
    }

    if (myproc == nprocs - 1) {
        ofs << "endsolid Created by AMReX\n";
    }
    ofs.close();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!ofs.fail(),
                                     "Could not complete marching-cubes STL output " + filename);

#ifdef AMREX_USE_MPI
    if (myproc < nprocs - 1) {
        int foo = 0;
        ParallelDescriptor::Send(&foo, 1, myproc + 1, 100, ParallelContext::CommunicatorSub());
    }
    ParallelDescriptor::Barrier(ParallelContext::CommunicatorSub());
#endif
}
} // namespace amrex::MC
