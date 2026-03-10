
#include <AMReX_Arena.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParallelDescriptor.H>

#include <AMReX_MarchingCubes.H>
#include <AMReX_mc_jgt_table.H>

#include <fstream>

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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void process_cube (std::int8_t ipass, LookUpTable const* lut, int i, int j, int k,
                   Array4<Real const> const& sdf, Array4<int const> const& ex,
                   Array4<int const> const& ey, Array4<int const> const& ez,
                   GpuArray<Real*,6> const& pvrtx, Array4<int> const& ntri,
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
                auto m = ntri(i,j,k,1) + r;
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

        if( std::abs( A*C - B*D ) < std::numeric_limits<Real>::epsilon() ) {
            return face >= 0 ;
        } else {
            return Real(face) * A * ( A*C - B*D ) >= 0  ;  // face and A invert signs
        }
    };

    auto test_interior = [&] (std::int8_t _case, std::int8_t _config, std::int8_t _subconfig, std::int8_t s)
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
        case  5 : if( At * Ct - Bt * Dt <  FLT_EPSILON ) { return s>0 ; } break;
        case  6 : return s>0 ;
        case  7 : return s<0 ;
        case  8 : return s>0 ;
        case  9 : return s>0 ;
        case 10 : if( At * Ct - Bt * Dt >= FLT_EPSILON ) { return s>0 ; } break;
        case 11 : return s<0 ;
        case 12 : return s>0 ;
        case 13 : return s<0 ;
        case 14 : return s<0 ;
        case 15 : return s<0 ;
        }

        return s<0 ;
    };

    auto add_c_vertex = [&] () -> int
    {
        if (ipass == 0) { ntri(i,j,k,2) = 1; return -1; }

        Real u = 0 ;
        int   vid ;

        auto m = ntri(i,j,k,3);
        auto& vert_x  = pvrtx[0][m];
        auto& vert_y  = pvrtx[1][m];
        auto& vert_z  = pvrtx[2][m];
        auto& vert_nx = pvrtx[3][m];
        auto& vert_ny = pvrtx[4][m];
        auto& vert_nz = pvrtx[5][m];

        vert_x = vert_y = vert_z = vert_nx = vert_ny = vert_nz = 0 ;

        auto update_vertex = [&] () {
            ++u ;
            vert_x  += pvrtx[0][vid];
            vert_y  += pvrtx[1][vid];
            vert_z  += pvrtx[2][vid];
            vert_nx += pvrtx[3][vid];
            vert_ny += pvrtx[4][vid];
            vert_nz += pvrtx[5][vid];
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

        u = std::sqrt( vert_nx * vert_nx + vert_ny * vert_ny +vert_nz * vert_nz ) ;
        if( u > 0 )
        {
            vert_nx *= Real(1)/u ;
            vert_ny *= Real(1)/u ;
            vert_nz *= Real(1)/u ;
        }

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
        ntri(i,j,k,0) = nt;
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
    nx.resize(n);
    ny.resize(n);
    nz.resize(n);
}

void Triangle::resize (int n)
{
    v1.resize(n);
    v2.resize(n);
    v3.resize(n);
}

void marching_cubes (Geometry const& geom, FArrayBox& sdf_fab, MCFab& mc_fab)
{
    BL_PROFILE("marching_cubes");

    AMREX_ALWAYS_ASSERT(sdf_fab.numPts() < Long(std::numeric_limits<int>::max()));

    // Remove small numbers.
    auto const& sdf = sdf_fab.array();
    ParallelFor(sdf_fab.box(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        if (std::abs(sdf(i,j,k)) < std::numeric_limits<Real>::epsilon()) {
            sdf(i,j,k) = std::numeric_limits<Real>::epsilon();
        }
    });

    Box nbox = sdf_fab.box();
    nbox.grow(-1); // Shrink the box by 1 so that we can compute gradient of sdf
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
    BoxIndexer n_bi(nbox);

    auto nvx = Scan::PrefixSum<int>(int(nbox.numPts()),
                                    [=] AMREX_GPU_DEVICE (int m) {
                                        auto [i,j,k] = n_bi(m);
                                        int vx = 0, vy = 0, vz = 0;
                                        if (ex.contains(i,j,k)) {
                                            if (((sdf(i,j,k) < 0) && (sdf(i+1,j,k) > 0)) ||
                                                ((sdf(i,j,k) > 0) && (sdf(i+1,j,k) < 0)))
                                            {
                                                vx = 1;
                                            }
                                            ex(i,j,k,0) = vx;
                                        }
                                        if (ey.contains(i,j,k)) {
                                            if (((sdf(i,j,k) < 0) && (sdf(i,j+1,k) > 0)) ||
                                                ((sdf(i,j,k) > 0) && (sdf(i,j+1,k) < 0)))
                                            {
                                                vy = 1;
                                            }
                                            ey(i,j,k,0) = vy;
                                        }
                                        if (ez.contains(i,j,k)) {
                                            if (((sdf(i,j,k) < 0) && (sdf(i,j,k+1) > 0)) ||
                                                ((sdf(i,j,k) > 0) && (sdf(i,j,k+1) < 0)))
                                            {
                                                vz = 1;
                                            }
                                            ez(i,j,k,0) = vz;
                                        }
                                        return vx + vy + vz;
                                    },
                                    [=] AMREX_GPU_DEVICE (int m, int ps) {
                                        auto [i,j,k] = n_bi(m);
                                        if (ex.contains(i,j,k)) {
                                            ex(i,j,k,1) = ps;
                                            if (ex(i,j,k,0)) { ++ps; }
                                        }
                                        if (ey.contains(i,j,k)) {
                                            ey(i,j,k,1) = ps;
                                            if (ey(i,j,k,0)) { ++ps; }
                                        }
                                        if (ez.contains(i,j,k)) {
                                            ez(i,j,k,1) = ps;
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
            Real u = sdf(i,j,k) / (sdf(i,j,k) - sdf(i+1,j,k));
            pvrtx[0][m] = Real(i) + u;
            pvrtx[1][m] = Real(j);
            pvrtx[2][m] = Real(k);
            Real nx = (Real(1)-u) * (sdf(i+1,j  ,k  ) - sdf(i-1,j  ,k  ))
                +              u  * (sdf(i+2,j,  k  ) - sdf(i  ,j  ,k  ));
            Real ny = (Real(1)-u) * (sdf(i  ,j+1,k  ) - sdf(i  ,j-1,k  ))
                +              u  * (sdf(i+1,j+1,k  ) - sdf(i+1,j-1,k  ));
            Real nz = (Real(1)-u) * (sdf(i  ,j  ,k+1) - sdf(i  ,j  ,k-1))
                +              u  * (sdf(i+1,j  ,k+1) - sdf(i+1,j  ,k-1));
            Real norm = std::sqrt(nx*nx + ny*ny + nz*nz);
            if (norm > 0) {
                nx *= Real(1)/norm;
                ny *= Real(1)/norm;
                nz *= Real(1)/norm;
            }
            pvrtx[3][m] = nx;
            pvrtx[4][m] = ny;
            pvrtx[5][m] = nz;
        }
        if (ey.contains(i,j,k) && ey(i,j,k,0)) {
            int m = ey(i,j,k,1);
            Real u = sdf(i,j,k) / (sdf(i,j,k) - sdf(i,j+1,k));
            pvrtx[0][m] = Real(i);
            pvrtx[1][m] = Real(j) + u;
            pvrtx[2][m] = Real(k);
            Real nx = (Real(1)-u) * (sdf(i+1,j  ,k  ) - sdf(i-1,j  ,k  ))
                +              u  * (sdf(i+1,j+1,k  ) - sdf(i-1,j+1,k  ));
            Real ny = (Real(1)-u) * (sdf(i  ,j+1,k  ) - sdf(i  ,j-1,k  ))
                +              u  * (sdf(i  ,j+2,k  ) - sdf(i  ,j  ,k  ));
            Real nz = (Real(1)-u) * (sdf(i  ,j  ,k+1) - sdf(i  ,j  ,k-1))
                +              u  * (sdf(i  ,j+1,k+1) - sdf(i  ,j+1,k-1));
            Real norm = std::sqrt(nx*nx + ny*ny + nz*nz);
            if (norm > 0) {
                nx *= Real(1)/norm;
                ny *= Real(1)/norm;
                nz *= Real(1)/norm;
            }
            pvrtx[3][m] = nx;
            pvrtx[4][m] = ny;
            pvrtx[5][m] = nz;
        }
        if (ez.contains(i,j,k) && ez(i,j,k,0)) {
            int m = ez(i,j,k,1);
            Real u = sdf(i,j,k) / (sdf(i,j,k) - sdf(i,j,k+1));
            pvrtx[0][m] = Real(i);
            pvrtx[1][m] = Real(j);
            pvrtx[2][m] = Real(k) + u;
            Real nx = (Real(1)-u) * (sdf(i+1,j  ,k  ) - sdf(i-1,j  ,k  ))
                +              u  * (sdf(i+1,j,  k+1) - sdf(i-1,j  ,k+1));
            Real ny = (Real(1)-u) * (sdf(i  ,j+1,k  ) - sdf(i  ,j-1,k  ))
                +              u  * (sdf(i  ,j+1,k+1) - sdf(i  ,j-1,k+1));
            Real nz = (Real(1)-u) * (sdf(i  ,j  ,k+1) - sdf(i  ,j  ,k-1))
                +              u  * (sdf(i  ,j  ,k+2) - sdf(i  ,j  ,k  ));
            Real norm = std::sqrt(nx*nx + ny*ny + nz*nz);
            if (norm > 0) {
                nx *= Real(1)/norm;
                ny *= Real(1)/norm;
                nz *= Real(1)/norm;
            }
            pvrtx[3][m] = nx;
            pvrtx[4][m] = ny;
            pvrtx[5][m] = nz;
        }
    });

    LookUpTable const* lut = d_table;
    auto const& sdf_c = sdf_fab.const_array();
    auto const& ex_c = ex_fab.const_array();
    auto const& ey_c = ey_fab.const_array();
    auto const& ez_c = ez_fab.const_array();

    // 0: # tri, 1: partial sum of # tri, 2: # c_vertex, 3: partial sum of # c_vertex
    BaseFab<int> ntri_fab(cbox,4);
    auto const& ntri = ntri_fab.array();
    GpuArray<int*,3> ptri{nullptr,nullptr,nullptr};

    Gpu::Buffer<int> error({0});
    auto* perror = error.data();

    BoxIndexer c_bi(cbox);

    int nttot = Scan::PrefixSum<int>(int(cbox.numPts()),
                                     [=] AMREX_GPU_DEVICE (int m) {
                                         auto [i,j,k] = c_bi(m);
                                         int err = 0;
                                         ntri(i,j,k,2) = 0;
                                         process_cube(0,lut,i,j,k,sdf_c,ex_c,ey_c,ez_c,pvrtx,ntri,ptri,&err);
                                         if (err != 0) {
                                             Gpu::Atomic::AddNoRet(perror, 1);
                                         }
                                         return ntri(i,j,k,0);
                                     },
                                     [=] AMREX_GPU_DEVICE (int m, int ps) {
                                         auto [i,j,k] = c_bi(m);
                                         ntri(i,j,k,1) = ps;
                                     },
                                     Scan::Type::exclusive, Scan::retSum);

    auto* nerror = error.copyToHost();
    if (*nerror > 0) {
        amrex::Abort("Marching Cubes: invalid triangle");
    }

    int nvx_c = Scan::PrefixSum<int>(int(cbox.numPts()),
                                     [=] AMREX_GPU_DEVICE (int m) {
                                         auto [i,j,k] = c_bi(m);
                                         return ntri(i,j,k,2);
                                     },
                                     [=] AMREX_GPU_DEVICE (int m, int ps) {
                                         auto [i,j,k] = c_bi(m);
                                         ntri(i,j,k,3) = ps;
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

    Gpu::streamSynchronize();

    mc_fab.m_triangles = std::move(tri);
    mc_fab.m_vertices = std::move(vrtx); // We can probably release the memory used by nx, ny, nz.
}

void write_stl (std::string const& filename, std::map<int,std::unique_ptr<MCFab>> const& mc_fabs)
{
    int myproc = ParallelDescriptor::MyProc();
    int nprocs = ParallelDescriptor::NProcs();

    std::ofstream ofs;

    if (myproc == 0) {
        ofs.open(filename);
        ofs << "solid Created by AMReX\n";
    }

#ifdef AMREX_USE_MPI
    if (myproc > 0) {
        int foo = 0;
        ParallelDescriptor::Recv(&foo, 1, myproc-1, 100);
    }
#endif

    if (! ofs.is_open()) {
        ofs.open(filename, std::ios_base::app);
    }

    for (auto const& [k,p] : mc_fabs) {
        auto ntri = int(p->m_triangles.v1.size());

#ifdef AMREX_USE_GPU
        Gpu::PinnedVector<int> tri_v1(ntri);
        Gpu::PinnedVector<int> tri_v2(ntri);
        Gpu::PinnedVector<int> tri_v3(ntri);
        auto nvert = p->m_vertices.x.size();
        Gpu::PinnedVector<Real> vert_x(nvert);
        Gpu::PinnedVector<Real> vert_y(nvert);
        Gpu::PinnedVector<Real> vert_z(nvert);
        Gpu::copyAsync(Gpu::deviceToHost, p->m_triangles.v1.begin(), p->m_triangles.v1.end(), tri_v1.begin());
        Gpu::copyAsync(Gpu::deviceToHost, p->m_triangles.v2.begin(), p->m_triangles.v2.end(), tri_v2.begin());
        Gpu::copyAsync(Gpu::deviceToHost, p->m_triangles.v3.begin(), p->m_triangles.v3.end(), tri_v3.begin());
        Gpu::copyAsync(Gpu::deviceToHost, p->m_vertices.x.begin(), p->m_vertices.x.end(), vert_x.begin());
        Gpu::copyAsync(Gpu::deviceToHost, p->m_vertices.y.begin(), p->m_vertices.y.end(), vert_y.begin());
        Gpu::copyAsync(Gpu::deviceToHost, p->m_vertices.z.begin(), p->m_vertices.z.end(), vert_z.begin());
        Gpu::streamSynchronize();
#else
        auto const& tri_v1 = p->m_triangles.v1;
        auto const& tri_v2 = p->m_triangles.v2;
        auto const& tri_v3 = p->m_triangles.v3;
        auto const& vert_x = p->m_vertices.x;
        auto const& vert_y = p->m_vertices.y;
        auto const& vert_z = p->m_vertices.z;
#endif
        for (int itri = 0; itri < ntri; ++itri) {
            auto iv1 = tri_v1[itri];
            auto iv2 = tri_v2[itri];
            auto iv3 = tri_v3[itri];
            XDim3 v1 = { vert_x[iv1], vert_y[iv1], vert_z[iv1] };
            XDim3 v2 = { vert_x[iv2], vert_y[iv2], vert_z[iv2] };
            XDim3 v3 = { vert_x[iv3], vert_y[iv3], vert_z[iv3] };
            XDim3 vec1{v2.x-v1.x, v2.y-v1.y, v2.z-v1.z};
            XDim3 vec2{v3.x-v2.x, v3.y-v2.y, v3.z-v2.z};
            XDim3 norm{vec1.y*vec2.z-vec1.z*vec2.y,
                       vec1.z*vec2.x-vec1.x*vec2.z,
                       vec1.x*vec2.y-vec1.y*vec2.x};
            auto tmp = std::sqrt(norm.x*norm.x + norm.y*norm.y + norm.z*norm.z);
            if (tmp != 0) { tmp = Real(1) / tmp; }
            ofs << "facet normal " << norm.x*tmp << " " << norm.y*tmp << " " << norm.z*tmp << "\n"
                << "  outer loop\n"
                << "    vertex " << v1.x << " " << v1.y << " " << v1.z << "\n"
                << "    vertex " << v2.x << " " << v2.y << " " << v2.z << "\n"
                << "    vertex " << v3.x << " " << v3.y << " " << v3.z << "\n"
                << "  endloop\n"
                << "endfacet\n";
        }
    }

#ifdef AMREX_USE_MPI
    if (myproc < nprocs-1) {
        int foo = 0;
        ParallelDescriptor::Send(&foo, 1, myproc+1, 100);
    }
#endif

    if (myproc == nprocs-1) {
        ofs << "endsolid Created by AMReX\n";
    }
}

}
