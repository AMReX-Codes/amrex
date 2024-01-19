#ifndef AMREX_TRACERPARTICLE_MOD_K_H
#define AMREX_TRACERPARTICLE_MOD_K_H

#include <AMReX_Config.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_TracerParticles.H>
#include <cmath>

namespace amrex {

//
// cic_interpolate computes the velocity at the particle by interpolating from velocity
//   components all co-located at cell centers
//
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cic_interpolate (const P& p,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                      const amrex::Array4<amrex::Real const>&  uarr,
                      amrex::ParticleReal * val, int M = AMREX_SPACEDIM)
{
    int icomp = 0;
    int ncomp_per_array = M;
    int num_arrays = 1;
    IntVect is_nodal = amrex::IntVect::TheZeroVector();
    linear_interpolate_to_particle (p, plo, dxi, &uarr, val, &is_nodal, icomp, ncomp_per_array, num_arrays);
}

//
// cic_interpolate computes the velocity at the particle by interpolating from velocity
//   components all co-located at nodes (cell corners)
//
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cic_interpolate_nd (const P& p,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                         const amrex::Array4<amrex::Real const>&  uarr,
                         amrex::ParticleReal * val, int M = AMREX_SPACEDIM)
{
    int icomp = 0;
    int ncomp_per_array = M;
    int num_arrays = 1;
    IntVect is_nodal = amrex::IntVect::TheUnitVector();
    linear_interpolate_to_particle (p, plo, dxi, &uarr, val, &is_nodal, icomp, ncomp_per_array, num_arrays);
}

//
// mac_interpolate computes the velocity at the particle by interpolating from normal velocity components
//   located at the faces of the cell the particle lives in
//
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mac_interpolate (const P& p,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                      amrex::GpuArray<amrex::Array4<amrex::Real const>,AMREX_SPACEDIM> const& umac_arr,
                      amrex::ParticleReal * val)
{
    int icomp = 0;
    int ncomp_per_array = 1;
    int num_arrays = AMREX_SPACEDIM;
    IntVect is_nodal[AMREX_SPACEDIM];
    for (int d=0; d < AMREX_SPACEDIM; ++d) {
        is_nodal[d] = amrex::IntVect::TheZeroVector();
        is_nodal[d][d] = 1;
    }
    linear_interpolate_to_particle (p, plo, dxi, &umac_arr[0], val, &is_nodal[0], icomp, ncomp_per_array, num_arrays);
}


//
// linear_interpolate computes the data at the particle by interpolating from mesh data
//
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void linear_interpolate_to_particle (const P& p,
                                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                     const Array4<amrex::Real const>* data_arr,
                                     amrex::ParticleReal * val,
                                     const IntVect* is_nodal,
                                     int start_comp, int ncomp, int num_arrays)
{
    AMREX_ASSERT(val != nullptr);

    int ctr = 0;

    for (int d = 0; d < num_arrays; d++)
    {
        AMREX_D_TERM(amrex::Real lx = (Real(p.pos(0))-plo[0])*dxi[0] - static_cast<Real>(!is_nodal[d][0])*Real(0.5);,
                     amrex::Real ly = (Real(p.pos(1))-plo[1])*dxi[1] - static_cast<Real>(!is_nodal[d][1])*Real(0.5);,
                     amrex::Real lz = (Real(p.pos(2))-plo[2])*dxi[2] - static_cast<Real>(!is_nodal[d][2])*Real(0.5));

        // (i0,j0,k0) is the lower corner of the box needed for interpolation
        // i0 = (i-1) if particle is lower  than center of cell i
        // i0 = (i  ) if particle is higher than center of cell i
        AMREX_D_TERM(int const i0 = static_cast<int>(amrex::Math::floor(lx));,
                     int const j0 = static_cast<int>(amrex::Math::floor(ly));,
                     int const k0 = static_cast<int>(amrex::Math::floor(lz)));

        AMREX_D_TERM(amrex::Real const xint = lx - static_cast<Real>(i0);,
                     amrex::Real const yint = ly - static_cast<Real>(j0);,
                     amrex::Real const zint = lz - static_cast<Real>(k0));

        amrex::Real sx[] = {amrex::Real(1.0) - xint, xint};
#if (AMREX_SPACEDIM > 1)
        amrex::Real sy[] = {amrex::Real(1.0) - yint, yint};
#endif
#if (AMREX_SPACEDIM > 2)
        amrex::Real sz[] = {amrex::Real(1.0) - zint, zint};
#endif

        for (int comp = start_comp; comp < ncomp; ++comp) {
            val[ctr] = ParticleReal(0.0);
#if (AMREX_SPACEDIM > 2)
            for (int kk = 0; kk <=1; ++kk) {
#endif

#if (AMREX_SPACEDIM > 1)
                for (int jj = 0; jj <= 1; ++jj) {
#endif
                    for (int ii = 0; ii <= 1; ++ii) {
                        val[ctr] += static_cast<ParticleReal>((data_arr[d])(IntVect(AMREX_D_DECL(i0+ii, j0+jj, k0+kk)), comp) *
                                                                                    AMREX_D_TERM(sx[ii],*sy[jj],*sz[kk]));
            AMREX_D_TERM(},},});
            ctr++;
        } // ncomp
    } // d
}

//
// linear_interpolate computes the data at the particle by interpolating from mesh data
//
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void linear_interpolate_to_particle_z (const P& p,
                                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                       const Array4<amrex::Real const>* data_arr,
                                       const amrex::Array4<amrex::Real const>&  height_arr,
                                       amrex::ParticleReal * val,
                                       const IntVect* is_nodal,
                                       int start_comp, int ncomp, int num_arrays)
{
    AMREX_ASSERT(val != nullptr);

    int ctr = 0;

    for (int d = 0; d < num_arrays; d++)
    {
        AMREX_D_TERM(amrex::Real lx = (Real(p.pos(0))-plo[0])*dxi[0] - static_cast<Real>(!is_nodal[d][0])*Real(0.5);,
                     amrex::Real ly = (Real(p.pos(1))-plo[1])*dxi[1] - static_cast<Real>(!is_nodal[d][1])*Real(0.5);,
                     amrex::Real lz = (Real(p.pos(2))-plo[2])*dxi[2] - static_cast<Real>(!is_nodal[d][2])*Real(0.5));

        int const i0 = static_cast<int>(amrex::Math::floor(lx));
        int const j  = static_cast<int>(amrex::Math::floor(ly));
        int k = 0; // temporary for 2D
        int k0 = 0; // temporary for 2D

//        AMREX_D_TERM(amrex::Real xoff =  static_cast<Real>(!is_nodal[d][0])*Real(0.5);,
//                     amrex::Real yoff =  static_cast<Real>(!is_nodal[d][1])*Real(0.5);,
//                     amrex::Real zoff =  static_cast<Real>(!is_nodal[d][2])*Real(0.5));
//

#if (AMREX_SPACEDIM == 2)

        amrex::Real hlo_xlo = 0.25 * ( height_arr(i0                       , j                     , k)
                                   +   height_arr(i0 + (!(is_nodal[d][0])) , j                     , k)
                                   +   height_arr(i0                       , j + (!is_nodal[d][1]) , k)
                                   +   height_arr(i0 + (!(is_nodal[d][0])) , j + (!is_nodal[d][1]) , k) );

        amrex::Real hlo_xhi = 0.25 * ( height_arr(i0 + 1                       , j                    , k )
                                   +   height_arr(i0 + 1 + (!(is_nodal[d][0])) , j                    , k )
                                   +   height_arr(i0 + 1                       , j + (!is_nodal[d][1]), k )
                                   +   height_arr(i0 + 1 + (!(is_nodal[d][0])) , j + (!is_nodal[d][1]), k ) );


        amrex::Real const xint = lx - static_cast<Real>(i0);
        amrex::Real sx[] = { amrex::Real(1.) - xint, xint};
        amrex::Real height_at_px = sx[0] * hlo_xlo + sx[1] * hlo_xhi;

        int const j0 = (amrex::Real(p.pos(1)) >= height_at_px) ? j : j-1;

        amrex::Real y1 = 0.25 * ( height_arr(i0                       , j0                    , k )
                              +   height_arr(i0 + (!(is_nodal[d][0])) , j0                    , k )
                              +   height_arr(i0                       , j0 + (!is_nodal[d][1]), k )
                              +   height_arr(i0 + (!(is_nodal[d][0])) , j0 + (!is_nodal[d][1]), k ) );

        amrex::Real y2 = 0.25 * ( height_arr(i0                       , j0 + 1                    , k )
                              +   height_arr(i0 + (!(is_nodal[d][0])) , j0 + 1                    , k )
                              +   height_arr(i0                       , j0 + 1 + (!is_nodal[d][1]), k )
                              +   height_arr(i0 + (!(is_nodal[d][0])) , j0 + 1 + (!is_nodal[d][1]), k ) );

        amrex::Real y3 = 0.25 * ( height_arr(i0 + 1                       , j0                    , k )
                              +   height_arr(i0 + 1 + (!(is_nodal[d][0])) , j0                    , k )
                              +   height_arr(i0 + 1                       , j0 + (!is_nodal[d][1]), k )
                              +   height_arr(i0 + 1 + (!(is_nodal[d][0])) , j0 + (!is_nodal[d][1]), k ) );

        amrex::Real y4 = 0.25 * ( height_arr(i0 + 1                       , j0 + 1                    , k )
                              +   height_arr(i0 + 1 + (!(is_nodal[d][0])) , j0 + 1                    , k )
                              +   height_arr(i0 + 1                       , j0 + 1 + (!is_nodal[d][1]), k )
                              +   height_arr(i0 + 1 + (!(is_nodal[d][0])) , j0 + 1 + (!is_nodal[d][1]), k ) );

        amrex::Real hint_ilo = (p.pos(1) - y1) / (y2 - y1);
        amrex::Real hint_ihi = (p.pos(1) - y3) / (y4 - y3);
        amrex::Real sy[] = { amrex::Real(1.) - hint_ilo, hint_ilo,
                             amrex::Real(1.) - hint_ihi, hint_ihi};

#elif (AMREX_SPACEDIM == 3)

        int const j0 = static_cast<int>(amrex::Math::floor(ly));,
        amrex::Real const xint = lx - static_cast<Real>(i0);
        amrex::Real const yint = ly - static_cast<Real>(j0);
        amrex::Real sx[] = { amrex::Real(1.) - xint, xint};
        amrex::Real sy[] = { amrex::Real(1.) - yint, yint};

//        amrex::Real hlo_xloylo = 0.125
//            * ( height_arr(i0                    , j0                    , k                    )
//            +   height_arr(i0 + (!is_nodal[d][0]), j0                    , k                    )
//            +   height_arr(i0                    , j0 + (!is_nodal[d][1]), k                    )
//            +   height_arr(i0 + (!is_nodal[d][0]), j0 + (!is_nodal[d][1]), k                    )
//            +   height_arr(i0                    , j0                    , k + (!is_nodal[d][2]))
//            +   height_arr(i0 + (!is_nodal[d][0]), j0                    , k + (!is_nodal[d][2]))
//            +   height_arr(i0                    , j0 + (!is_nodal[d][1]), k + (!is_nodal[d][2]))
//            +   height_arr(i0 + (!is_nodal[d][0]), j0 + (!is_nodal[d][1]), k + (!is_nodal[d][2]))
//            ); 
//        amrex::Real hlo_xloyhi = 0.125 
//            * ( height_arr(i0                    , j0 + 1                    , k                    )
//            +   height_arr(i0 + (!is_nodal[d][0]), j0 + 1                    , k                    )
//            +   height_arr(i0                    , j0 + 1 + (!is_nodal[d][1]), k                    )
//            +   height_arr(i0 + (!is_nodal[d][0]), j0 + 1 + (!is_nodal[d][1]), k                    )
//            +   height_arr(i0                    , j0 + 1                    , k + (!is_nodal[d][2]))
//            +   height_arr(i0 + (!is_nodal[d][0]), j0 + 1                    , k + (!is_nodal[d][2]))
//            +   height_arr(i0                    , j0 + 1 + (!is_nodal[d][1]), k + (!is_nodal[d][2]))
//            +   height_arr(i0 + (!is_nodal[d][0]), j0 + 1 + (!is_nodal[d][1]), k + (!is_nodal[d][2]))
//            ); 
//        amrex::Real hlo_xhiylo = 0.125
//            * ( height_arr(i0 + 1                    , j0                    , k                    )
//            +   height_arr(i0 + 1 + (!is_nodal[d][0]), j0                    , k                    )
//            +   height_arr(i0 + 1                    , j0 + (!is_nodal[d][1]), k                    )
//            +   height_arr(i0 + 1 + (!is_nodal[d][0]), j0 + (!is_nodal[d][1]), k                    )
//            +   height_arr(i0 + 1                    , j0                    , k + (!is_nodal[d][2]))
//            +   height_arr(i0 + 1 + (!is_nodal[d][0]), j0                    , k + (!is_nodal[d][2]))
//            +   height_arr(i0 + 1                    , j0 + (!is_nodal[d][1]), k + (!is_nodal[d][2]))
//            +   height_arr(i0 + 1 + (!is_nodal[d][0]), j0 + (!is_nodal[d][1]), k + (!is_nodal[d][2]))
//            ); 
//        amrex::Real hlo_xhiyhi = 
//            * ( height_arr(i0 + 1                    , j0 + 1                    , k                    )
//            +   height_arr(i0 + 1 + (!is_nodal[d][0]), j0 + 1                    , k                    )
//            +   height_arr(i0 + 1                    , j0 + 1 + (!is_nodal[d][1]), k                    )
//            +   height_arr(i0 + 1 + (!is_nodal[d][0]), j0 + 1 + (!is_nodal[d][1]), k                    )
//            +   height_arr(i0 + 1                    , j0 + 1                    , k + (!is_nodal[d][2]))
//            +   height_arr(i0 + 1 + (!is_nodal[d][0]), j0 + 1                    , k + (!is_nodal[d][2]))
//            +   height_arr(i0 + 1                    , j0 + 1 + (!is_nodal[d][1]), k + (!is_nodal[d][2]))
//            +   height_arr(i0 + 1 + (!is_nodal[d][0]), j0 + 1 + (!is_nodal[d][1]), k + (!is_nodal[d][2]))
//            ); 
//        amrex::Real height_at_pxy = hlo_xloylo * sx[0] * sy[0] + hlo_xloyhi * sx[0] * sy[1]
//                                  + hlo_xhiylo * sx[1] * sy[0] + hlo_xhiyhi * sx[1] * sy[1]; 

        amrex::Real hlo[4];
        ilo = 0;
        amrex::Real height_at_pxy = 0.;
        for (ii = 0; ii < 2; ++ii) {
            for (jj = 0; jj < 2; ++jj) {
                hlo[ilo] = 0.125
                * ( height_arr(i0 + ii                    , j0 + jj                    , k                    )
                +   height_arr(i0 + ii + (!is_nodal[d][0]), j0 + jj                    , k                    )
                +   height_arr(i0 + ii                    , j0 + jj + (!is_nodal[d][1]), k                    )
                +   height_arr(i0 + ii + (!is_nodal[d][0]), j0 + jj + (!is_nodal[d][1]), k                    )
                +   height_arr(i0 + ii                    , j0 + jj                    , k + (!is_nodal[d][2]))
                +   height_arr(i0 + ii + (!is_nodal[d][0]), j0 + jj                    , k + (!is_nodal[d][2]))
                +   height_arr(i0 + ii                    , j0 + jj + (!is_nodal[d][1]), k + (!is_nodal[d][2]))
                +   height_arr(i0 + ii + (!is_nodal[d][0]), j0 + jj + (!is_nodal[d][1]), k + (!is_nodal[d][2]))
                );
                height_at_pxy += hlo[ilo] * sx[ii] * sy[jj];
                ++ilo;
            }
        }

        int const k0 = (amrex::Real(p.pos(2)) >= height_at_pxy ) ? k : k-1;
//
//        amrex::Real x1 = (static_cast<Real>(i  ) + xoff) * dx[0];
//        amrex::Real x2 = (static_cast<Real>(i  ) + xoff) * dx[0] - x1;
//        amrex::Real x3 = (static_cast<Real>(i+1) + xoff) * dx[0] - x1;
//        amrex::Real x4 = (static_cast<Real>(i+1) + xoff) * dx[0] - x1;
//
//        amrex::Real y1 = (static_cast<Real>(j  ) + yoff) * dx[1];
//        amrex::Real y2 = (static_cast<Real>(j+1) + yoff) * dx[1] - y1;
//        amrex::Real y3 = (static_cast<Real>(j  ) + yoff) * dx[1] - y1;
//        amrex::Real y4 = (static_cast<Real>(j+1) + yoff) * dx[1] - y1;
//
//        amrex::Real z1 = height_arr(i0  ,j0  ); // ONLY IF NODAL
//        amrex::Real z2 = height_arr(i0  ,j0+1); // ONLY IF NODAL
//        amrex::Real z3 = height_arr(i0+1,j0  ); // ONLY IF NODAL
//        amrex::Real z4 = height_arr(i0+1,j0+1); // ONLY IF NODAL
#endif
//
//        AMREX_D_TERM(amrex::Real x = p.pos(0) - x1;,
//                     amrex::Real y = p.pos(1) - y1;,
//                     amrex::Real z = p.pos(2) - z1;);
//
        for (int comp = start_comp; comp < ncomp; ++comp) {
            val[ctr] = amrex::ParticleReal(0.);
#if (AMREX_SPACEDIM == 2)
            // interpolate 4 points in y
            for (int jj = 0; jj <= 3; ++jj) {
                for (int ii = 0; ii <=1; ++ii) {
                    val[ctr] += static_cast<ParticleReal>( (data_arr[d])(i0+ii, j0+jj, k0 ,comp)*sx[ii]*sy[jj]);
                }
            }
#endif
//#if (AMREX_SPACEDIM == 2)
//            f1 = data_arr(i0  ,j0  ,k0);
//            f2 = data_arr(i0+1,j0  ,k0) - f1;
//            f3 = data_arr(i0  ,j0+1,k0) - f1;
//            f4 = data_arr(i0+1,j0+1,k0) - f1;
//
//            amrex::Real det = -x3*y2*(y3 - y4);
//
//            amrex::Real b = ( -y2*y3*f4 + y2*y4*f3         ) / det;
//            amrex::Real c = ( -x3*y3*f2 + x3*y4*f2         ) / det;
//            amrex::Real d = ( y3*f2 + y2*f4 - y4*f2 - y2*f3) / det;
//
//            val[ctr] = f1 + b*x + c*y + d*x*y;
//#elif (AMREX_SPACEDIM == 3)
//#endif
            ctr++;
        } // ncomp
    } // d
}

#if 0
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
bilinear_interp_for_3d (amrex::Real x, amrex::Real y,
                        amrex::Real const x2,amrex::Real const x3,amrex::Real const x4,
                        amrex::Real const y2,amrex::Real const y3,amrex::Real const y4,
                        amrex::Real const f2,amrex::Real const f3,amrex::Real const f4)
{
    // Here we interpolate between four arbitrary locations (with positions calculated as offsets from (x1,y1)),
    //      but we define the values relative to point 1, so we can set f1 = 0 in the formulae below
    // amrex::Real a = (f1*(x2 - x3)*x4*y2*y3 +
    //                 (f1*x3*(-x2 + x4)*y2 + f1*x2*(x3 - x4)*y3)*y4)/det;
    // amrex::Real b = ( f1*(-x2 + x3)*y2*y3 +
    //                   f4*(x2*y2*y3 - x3*y2*y3) +
    //                  (f1*x2*y2 - f1*x4*y2 - (f1 - f2)*(x3 - x4)*y3)*y4 +
    //                   f3*(-(x2*y2*y4) + x4*y2*y4) ) / det;
    // amrex::Real c = ( f1*x2*x3*y2 - f1*x2*x4*y2 - f1*x2*x3*y3 + f1*x3*x4*y3 - f2*x3*x4*y3 +
    //                   f4*x2*x3*(-y2 + y3) + f3*x2*x4*(y2 - y4) + f1*(x2 - x3)*x4*y4 +
    //                   f2*x3*x4*y4 )/det;
    // amrex::Real d = ( -(f1*x3*y2) + f1*x4*y2 + f1*x2*y3 - f1*x4*y3 +
    //                     f2*x4*y3 + f4*(x3*y2 - x2*y3) - f1*x2*y4 + f1*x3*y4 - f2*x3*y4 +
    //                     f3*(-(x4*y2) + x2*y4) ) / det;
    // amrex::Real f = a + b*x + c*y + d*x*y;

    amrex::Real det = x2*x4*y3*(y2 - y4) - x3*x4*y2*(y3 - y4) - x2*x3*(y2 - y3)*y4;

    // amrex::Real a = Real(0.0);

    amrex::Real b = ( f4*(x2*y2*y3 - x3*y2*y3) + (f2*(x3 - x4)*y3)*y4 + f3*(-(x2*y2*y4) + x4*y2*y4)          ) / det;
    amrex::Real c = ( -f2*x3*x4*y3 + f4*x2*x3*(-y2 + y3) + f3*x2*x4*(y2 - y4) + f2*x3*x4*y4                  ) /det;
    amrex::Real d = ( f2*x4*y3 + f4*(x3*y2 - x2*y3) - f2*x3*y4 + f3*(-(x4*y2) + x2*y4) ) / det;

    amrex::Real f = b*x + c*y + d*x*y;
    return f;
}
#endif

//
// cic_interpolate_mapped computes the velocity at the particle by interpolating from velocity
//   components all co-located at cell centers
// This version assumes that in 2D: node (i,j)   lives at physical location (i*dx, j*dx)
//                                  and cell (i,j) that the particle lives in is carried in (p.idata(0), p.idata(1))
//                           in 3D: node (i,j,k) lives at physical location (i*dx, j*dy, k*dz)
//                                  and cell (i,j,k) that the particle lives in is carried in (p.idata(0), p.idata(1), p.idata(2))
//
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cic_interpolate_mapped (const P& p,
                             const amrex::Array4<amrex::Real const>&  ucc_arr,
                             const amrex::Array4<amrex::Real const>&  loc_arr,
                             amrex::ParticleReal* val)
{
    AMREX_ASSERT(val != nullptr);

    // First identify the cell (i,j,k) the particle is in
    AMREX_D_TERM(int i = p.idata(0);,
                 int j = p.idata(1);,
                 int k = p.idata(2));

    AMREX_D_TERM(amrex::Real lx = amrex::Real(p.pos(0)) - 0.5*(loc_arr(IntVect(AMREX_D_DECL(i-1,j,k)),0) + loc_arr(IntVect(AMREX_D_DECL(i,j,k)),0));,
                 amrex::Real ly = amrex::Real(p.pos(1)) - 0.5*(loc_arr(IntVect(AMREX_D_DECL(i,j-1,k)),1) + loc_arr(IntVect(AMREX_D_DECL(i,j,k)),1));,
                 amrex::Real lz = amrex::Real(p.pos(2)) - 0.5*(loc_arr(IntVect(AMREX_D_DECL(i,j,k-1)),2) + loc_arr(IntVect(AMREX_D_DECL(i,j,k)),2)));

    // (i0,j0,k0) is the lower corner of the box needed for interpolation
    // i0 = (i-1) if particle is lower  than center of cell i
    // i0 = (i  ) if particle is higher than center of cell i
    AMREX_D_TERM(int i0 = (lx >= 0) ? i : i+1;,
                 int j0 = (ly >= 0) ? j : j+1;,
                 int k0 = (lz >= 0) ? k : k+1;);

    AMREX_D_TERM(amrex::Real const xint = 2.0*lx/(loc_arr(IntVect(AMREX_D_DECL(i0+2,j,k)),0) - loc_arr(IntVect(AMREX_D_DECL(i0,j,k)),0));,
                 amrex::Real const yint = 2.0*ly/(loc_arr(IntVect(AMREX_D_DECL(i,j0+2,k)),1) - loc_arr(IntVect(AMREX_D_DECL(i,j0,k)),1));,
                 amrex::Real const zint = 2.0*lz/(loc_arr(IntVect(AMREX_D_DECL(i,j,k0+2)),2) - loc_arr(IntVect(AMREX_D_DECL(i,j,k0)),2)));

    amrex::Real sx[] = {amrex::Real(1.0) - xint, xint};
    amrex::Real sy[] = {amrex::Real(1.0) - yint, yint};

#if (AMREX_SPACEDIM == 2)
//    amrex::Real u = bilinear_interp(x1,x2-x1,x3-x1,x4-x1,y1,y2-y1,y3-y1,y4-y1,f1,f2,f3,f4);
//    amrex::Real v = bilinear_interp(x1,x2-x1,x3-x1,x4-x1,y1,y2-y1,y3-y1,y4-y1,f1,f2,f3,f4);
#elif (AMREX_SPACEDIM == 3)
    amrex::Real sz[] = {amrex::Real(1.0) - zint, zint};
    amrex::Abort("cic_interpolate_mapped has not yet been implemented in 3D");
#endif // SPACEDIM
}

//
// cic_interpolate_mapped_z computes the velocity at the particle by interpolating from velocity
//     components co-located at cell centers
// This version assumes that in 2D: node (i,j)   lives at physical location (i*dx, height_arr()),
//                                  and we can easily calculate the "i" of the particle, while
//                                  "j" is carried in p.idata(0)
//                           in 3D: node (i,j,k) lives at physical location (i*dx, j*dy, height_arr())
//                                  and we can easily calculate the "i" and "j" of the particle, while
//                                  "k" is carried in p.idata(0)
//
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cic_interpolate_mapped_z (const P& p,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                               const amrex::Array4<amrex::Real const>&  ucc_arr,
                               const amrex::Array4<amrex::Real const>&  loc_arr,
                               amrex::ParticleReal* val)
{
    AMREX_ASSERT(val != nullptr);

#if 0
    amrex::Real lx = (amrex::Real(p.pos(0)) - plo[0]) * dxi[0] - amrex::Real(0.5);
    amrex::Real ly = (amrex::Real(p.pos(1)) - plo[1]) * dxi[1] - amrex::Real(0.5);

    // (i0,j0,k0) is the lower corner of the box needed for interpolation
    int const i0 = static_cast<int>(amrex::Math::floor(lx));
    int const j0 = static_cast<int>(amrex::Math::floor(ly));

#if (AMREX_SPACEDIM == 3)
    int const k  = p.idata(2);
    amrex::Real lz = amrex::Real(p.pos(2)) - 0.5*(loc_arr(IntVect(AMREX_D_DECL(i,j,k-1)),2) + loc_arr(IntVect(AMREX_D_DECL(i,j,k)),2));
    int k0 = (lz >= 0) ? k : k+1;
#endif // SPACEDIM

    // z_cc is the center of the cell that the particle is in
    amrex::Real z_cc = 0.125*(loc_arr(IntVect(AMREX_D_DECL(i0  ,j0  ,k-1)),2) + loc_arr(IntVect(AMREX_D_DECL(i0  ,j0  ,k)),2) +
                              loc_arr(IntVect(AMREX_D_DECL(i0  ,j0+1,k-1)),2) + loc_arr(IntVect(AMREX_D_DECL(i0  ,j0+1,k)),2) +
                              loc_arr(IntVect(AMREX_D_DECL(i0+1,j0  ,k-1)),2) + loc_arr(IntVect(AMREX_D_DECL(i0+1,j0  ,k)),2) +
                              loc_arr(IntVect(AMREX_D_DECL(i0+1,j0+1,k-1)),2) + loc_arr(IntVect(AMREX_D_DECL(i0+1,j0+1,k)),2) );


    // amrex::Real const xint = lx - static_cast<Real>(i);
    // amrex::Real const yint = ly - static_cast<Real>(j);
#if (AMREX_SPACEDIM == 3)
    // amrex::Real const zint = 2.0*lz/(loc_arr(IntVect(AMREX_D_DECL(i,j,k0+2)),2) - loc_arr(IntVect(AMREX_D_DECL(i,j,k0)),2));
#endif // SPACEDIM

    amrex::Real sx[] = {amrex::Real(1.0) - xint, xint};,
    amrex::Real sy[] = {amrex::Real(1.0) - yint, yint};,

    for (int d=0; d < M; ++d) {
        val[d] = ParticleReal(0.0);
#if (AMREX_SPACEDIM == 2)
        for (int jj = 0; jj <= 1; ++jj) {
            for (int ii = 0; ii <= 1; ++ii) {
                val[d] += static_cast<ParticleReal>(sx[ii] * sy[jj] * uccarr(IntVect(i0+ii,j0+jj)), d);
            } // ii
        } // jj
#elif (AMREX_SPACEDIM == 3)
        amrex::Real sz[] = {amrex::Real(1.0) - zint, zint};
        for (int kk = 0; kk <= 1; ++kk) {
            for (int jj = 0; jj <= 1; ++jj) {
                for (int ii = 0; ii <= 1; ++ii) {
                    val[d] += static_cast<ParticleReal>(sx[ii] * sy[jj] * sz[kk] * uccarr(IntVect(i0+ii,j0+jj,k0+kk),d));
                } // ii
            } // jj
        } // kk
#endif // SPACEDIM
    } // d
#endif // if 0
}

//
// mac_interpolate computes the velocity at the particle by interpolating from normal velocity components
//   located at the faces of the cell the particle lives in
// This version assumes that in 2D: node (i,j)   lives at physical location (i*dx, height_arr()),
//                                  and we can easily calculate the "i" of the particle, while
//                                  "j" is carried in p.idata(0)
//                           in 3D: node (i,j,k) lives at physical location (i*dx, j*dy, height_arr())
//                                  and we can easily calculate the "i" and "j" of the particle, while
//                                  "k" is carried in p.idata(0)
// We use the normal velocity components at cell faces for interpolation
//
#if (AMREX_SPACEDIM > 1)
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mac_interpolate_mapped_z (const P& p,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                               amrex::GpuArray<amrex::Array4<amrex::Real const>,AMREX_SPACEDIM> const& umac_arr,
                               const amrex::Array4<amrex::Real const>&  height_arr,
                               amrex::ParticleReal * val)
{
    AMREX_ASSERT(val != nullptr);

#if (AMREX_SPACEDIM == 2)
    // Particle is in cell (i,j)
    int const i = static_cast<int>((Real(p.pos(0))-plo[0])*dxi[0]);
    int const j  = p.idata(0);
    int k = 0;

    // ************************************************************************************
    // For x-vel
    // ************************************************************************************
    int d = 0;

    // If doing x-vel, i0 = i
    amrex::Real lx = (Real(p.pos(0))-plo[0])*dxi[0] - static_cast<Real>(d != 0)*amrex::Real(0.5);
    int         i0 = static_cast<int>(amrex::Math::floor(lx));

    amrex::Real height_of_u_lo_iside = 0.5*(height_arr(IntVect(i  ,j)) + height_arr(IntVect(i  ,j+1)));
    amrex::Real height_of_u_hi_iside = 0.5*(height_arr(IntVect(i+1,j)) + height_arr(IntVect(i+1,j+1)));

    amrex::Real const xint = lx - static_cast<Real>(i0);
    amrex::Real height_of_u_at_p = xint * height_of_u_lo_iside + (1.0-xint) * height_of_u_hi_iside;
    amrex::Real val_of_u_at_p_j = xint * umac_arr[0](i,j,k) + (1.0-xint) * umac_arr[0](i+1,j,k);

    amrex::Real ly = amrex::Real(p.pos(1)) - height_of_u_at_p;

    int jother = (ly >= 0) ? j+1 : j-1;

    // Now that we know whether to use {(i,j) and (i,j+1)} OR {(i,j) and (i,j-1)} ...
    amrex::Real height_of_u_lo_iside_other_j = 0.5*(height_arr(IntVect(i  ,jother)) + height_arr(IntVect(i  ,jother+1)));
    amrex::Real height_of_u_hi_iside_other_j = 0.5*(height_arr(IntVect(i+1,jother)) + height_arr(IntVect(i+1,jother+1)));

    amrex::Real height_of_u_at_p_other_j = xint * height_of_u_lo_iside_other_j + (1.0-xint) * height_of_u_hi_iside_other_j;

    amrex::Real val_of_u_at_p_other_j = xint * umac_arr[0](i,jother,k) + (1.0-xint) * umac_arr[0](i+1,jother,k);

    amrex::Real dy   = std::abs(height_of_u_at_p_other_j - height_of_u_at_p);
    amrex::Real yint = std::abs(height_of_u_at_p_other_j - p.pos(1)) / dy;

    val[0] = ParticleReal(yint * val_of_u_at_p_j + (1.0-yint) * val_of_u_at_p_other_j);

    // ************************************************************************************
    // For y-vel
    // ************************************************************************************

    d = 1;
    // If doing y-vel, i0 = i or i-1
    lx = (Real(p.pos(0))-plo[0])*dxi[0] - static_cast<Real>(d != 0)*amrex::Real(0.5);
    i0 = static_cast<int>(amrex::Math::floor(lx));

    amrex::Real dx = 1.0 / dxi[0];

    amrex::Real x1 = (i0+0.5) * dx;
    amrex::Real x2 = (i0+0.5) * dx - x1;
    amrex::Real x3 = (i0+1.5) * dx - x1;
    amrex::Real x4 = (i0+1.5) * dx - x1;

    amrex::Real y1 = 0.5*(height_arr(IntVect(i0  ,j  )) + height_arr(IntVect(i0+1,j  )));
    amrex::Real y2 = 0.5*(height_arr(IntVect(i0  ,j+1)) + height_arr(IntVect(i0+1,j+1))) - y1;
    amrex::Real y3 = 0.5*(height_arr(IntVect(i0+1,j  )) + height_arr(IntVect(i0+1,j  ))) - y1;
    amrex::Real y4 = 0.5*(height_arr(IntVect(i0+1,j+1)) + height_arr(IntVect(i0+1,j+1))) - y1;

    amrex::Real f1 = umac_arr[1](i0  ,j  ,k);
    amrex::Real f2 = umac_arr[1](i0  ,j+1,k) - f1;
    amrex::Real f3 = umac_arr[1](i0+1,j  ,k) - f1;
    amrex::Real f4 = umac_arr[1](i0+1,j+1,k) - f1;

    // amrex::Print() << "F " << f1 << " " << f2 << " " << f3 << " " << f4 << std::endl;

    // val[1] = bilinear_interp_2d (p.pos(0)-x1, p.pos(1)-y1,
    //                              x2, x3, x4,
    //                              y2, y3, y4,
    //                              f2, f3, f4);
    val[1] += f1;

#elif (AMREX_SPACEDIM == 3)

    // Particle is in cell (i,j)
    int const i = static_cast<int>((Real(p.pos(0))-plo[0])*dxi[0]);
    int const j = static_cast<int>((Real(p.pos(1))-plo[1])*dxi[1]);
    int const k  = p.idata(0);

    // ************************************************************************************
    // For x-vel
    // ************************************************************************************
    int d = 0;

    // If doing x-vel, i0 = i
    amrex::Real lx = (Real(p.pos(0))-plo[0])*dxi[0] - static_cast<Real>(d != 0)*amrex::Real(0.5);
    int         i0 = static_cast<int>(amrex::Math::floor(lx));

    amrex::Real ly = (Real(p.pos(1))-plo[1])*dxi[1] - static_cast<Real>(d != 0)*amrex::Real(0.5);
    int         j0 = static_cast<int>(amrex::Math::floor(ly));

    amrex::Real height_of_u_lo_iside = Real(0.25)*(height_arr(i  ,j  ,k) + height_arr(i  ,j  ,k+1) +
                                                   height_arr(i  ,j+1,k) + height_arr(i  ,j+1,k+1) );
    amrex::Real height_of_u_hi_iside = Real(0.25)*(height_arr(i+1,j  ,k) + height_arr(i+1,j  ,k+1) +
                                                   height_arr(i+1,j+1,k) + height_arr(i+1,j+1,k+1) );

    amrex::Real const xint = lx - static_cast<Real>(i0);
    amrex::Real height_of_u_at_p = xint * height_of_u_lo_iside + (1.0-xint) * height_of_u_hi_iside;
    amrex::Real val_of_u_at_p_k = xint * umac_arr[0](i,j,k) + (1.0-xint) * umac_arr[0](i+1,j,k);

    amrex::Real lz = amrex::Real(p.pos(2)) - height_of_u_at_p;

    int kother = (lz >= 0) ? k+1 : k-1;

    // Now that we know whether to use {(i,j,k) and (i,j,k+1)} OR {(i,j,k) and (i,j,k-1)} ...
    amrex::Real height_of_u_lo_iside_other_k = 0.5*(height_arr(IntVect(i  ,j,kother)) + height_arr(IntVect(i  ,j,kother+1)));
    amrex::Real height_of_u_hi_iside_other_k = 0.5*(height_arr(IntVect(i+1,j,kother)) + height_arr(IntVect(i+1,j,kother+1)));

    amrex::Real height_of_u_at_p_other_k = xint * height_of_u_lo_iside_other_k + (1.0-xint) * height_of_u_hi_iside_other_k;

    amrex::Real val_of_u_at_p_other_k = xint * umac_arr[0](i,j,kother) + (1.0-xint) * umac_arr[0](i+1,j,kother);

    amrex::Real dz   = std::abs(height_of_u_at_p_other_k - height_of_u_at_p);
    amrex::Real zint = std::abs(height_of_u_at_p_other_k - p.pos(2)) / dz;

    val[0] = ParticleReal(zint * val_of_u_at_p_k + (1.0-zint) * val_of_u_at_p_other_k);

    // ************************************************************************************
    // For y-vel
    // ************************************************************************************
    // THIS IS JUST A HACK FOR NOW
    val[1] = 0.0;

    // ************************************************************************************
    // For z-vel
    // ************************************************************************************

    d = 2;
    // If doing z-vel, i0 = i or i-1
    lx = (Real(p.pos(0))-plo[0])*dxi[0] - static_cast<Real>(d != 0)*amrex::Real(0.5);
    i0 = static_cast<int>(amrex::Math::floor(lx));

    amrex::Real dx = 1.0 / dxi[0];

    amrex::Real x1 = (i0+0.5) * dx;
    amrex::Real x2 = (i0+0.5) * dx - x1;
    amrex::Real x3 = (i0+1.5) * dx - x1;
    amrex::Real x4 = (i0+1.5) * dx - x1;

    amrex::Real z1 = 0.5*(height_arr(IntVect(i0  ,j,k  )) + height_arr(IntVect(i0+1,j,k  )));
    amrex::Real z2 = 0.5*(height_arr(IntVect(i0  ,j,k+1)) + height_arr(IntVect(i0+1,j,k+1))) - z1;
    amrex::Real z3 = 0.5*(height_arr(IntVect(i0+1,j,k  )) + height_arr(IntVect(i0+1,j,k  ))) - z1;
    amrex::Real z4 = 0.5*(height_arr(IntVect(i0+1,j,k+1)) + height_arr(IntVect(i0+1,j,k+1))) - z1;

    amrex::Real f1 = umac_arr[1](i0  ,j,k  ,2);
    amrex::Real f2 = umac_arr[1](i0  ,j,k+1,2) - f1;
    amrex::Real f3 = umac_arr[1](i0+1,j,k  ,2) - f1;
    amrex::Real f4 = umac_arr[1](i0+1,j,k+1,2) - f1;

    amrex::Print() << "F " << f1 << " " << f2 << " " << f3 << " " << f4 << std::endl;

    // val[2] = bilinear_interp_2d (p.pos(0)-x1, p.pos(2)-z1,
    //                              x2, x3, x4,
    //                              z2, z3, z4,
    //                              f2, f3, f4);
    val[2] += f1;

#endif // SPACEDIM
}
#endif // if (SPACEDIM > 1)

}  // namespace amrex
#endif  // include guard
