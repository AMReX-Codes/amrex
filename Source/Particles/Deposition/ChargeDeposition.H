#ifndef CHARGEDEPOSITION_H_
#define CHARGEDEPOSITION_H_

#include "ShapeFactors.H"

/* \brief Charge Deposition for thread thread_num
 * /param xp, yp, zp   : Pointer to arrays of particle positions.
 * \param wp           : Pointer to array of particle weights.
 * \param rho_arr      : Array4 of charge density, either full array or tile.
 * \param np_to_depose : Number of particles for which current is deposited.
 * \param dx           : 3D cell size
 * \param xyzmin       : Physical lower bounds of domain.
 * \param lo           : Index lower bounds of domain.
 * /param q            : species charge.
 */
template <int depos_order>
void doChargeDepositionShapeN(const amrex::Real * const xp, 
                              const amrex::Real * const yp, 
                              const amrex::Real * const zp,
                              const amrex::Real * const wp,
                              const amrex::Array4<amrex::Real>& rho_arr,
                              const long np_to_depose,
                              const std::array<amrex::Real,3>& dx,
                              const std::array<amrex::Real, 3> xyzmin,
                              const amrex::Dim3 lo,
                              const amrex::Real q)
{
    const amrex::Real dxi = 1.0/dx[0];
    const amrex::Real dzi = 1.0/dx[2];
#if (AMREX_SPACEDIM == 2)
    const amrex::Real invvol = dxi*dzi;
#elif (defined WARPX_DIM_3D)
    const amrex::Real dyi = 1.0/dx[1];
    const amrex::Real invvol = dxi*dyi*dzi;
#endif

    const amrex::Real xmin = xyzmin[0];
    const amrex::Real ymin = xyzmin[1];
    const amrex::Real zmin = xyzmin[2];

    // Loop over particles and deposit into rho_arr
    amrex::ParallelFor(
        np_to_depose,
        [=] AMREX_GPU_DEVICE (long ip) {
            // --- Get particle quantities
            const amrex::Real wq = q*wp[ip]*invvol;

            // --- Compute shape factors
            // x direction
            // Get particle position in grid coordinates
#if (defined WARPX_DIM_RZ)
            const amrex::Real r = std::sqrt(xp[ip]*xp[ip] + yp[ip]*yp[ip]);
            const amrex::Real x = (r - xmin)*dxi;
#else
            const amrex::Real x = (xp[ip] - xmin)*dxi;
#endif
            // Compute shape factors for node-centered quantities
            amrex::Real AMREX_RESTRICT sx[depos_order + 1];
            // i: leftmost grid point (node-centered) that the particle touches
            const int i = compute_shape_factor<depos_order>(sx,  x);
                     
#if (defined WARPX_DIM_3D)
            // y direction
            const amrex::Real y = (yp[ip] - ymin)*dyi;
            amrex::Real AMREX_RESTRICT sy[depos_order + 1];
            const int j = compute_shape_factor<depos_order>(sy,  y);
#endif
            // z direction
            const amrex::Real z = (zp[ip] - zmin)*dzi;
            amrex::Real AMREX_RESTRICT sz[depos_order + 1];
            const int k = compute_shape_factor<depos_order>(sz,  z);

            // Deposit charge into rho_arr
#if (defined WARPX_DIM_2D) || (defined WARPX_DIM_RZ)
            for (int iz=0; iz<=depos_order; iz++){
                for (int ix=0; ix<=depos_order; ix++){
                    amrex::Gpu::Atomic::Add(
                        &rho_arr(lo.x+i+ix, lo.y+k+iz, 0), 
                        sx[ix]*sz[iz]*wq);
                }
            }
#elif (defined WARPX_DIM_3D)
            for (int iz=0; iz<=depos_order; iz++){
                for (int iy=0; iy<=depos_order; iy++){
                    for (int ix=0; ix<=depos_order; ix++){
                        amrex::Gpu::Atomic::Add(
                            &rho_arr(lo.x+i+ix, lo.y+j+iy, lo.z+k+iz),
                            sx[ix]*sy[iy]*sz[iz]*wq);
                    }
                }
            }
#endif
        }
        );
}

#endif // CHARGEDEPOSITION_H_
