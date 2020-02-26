/* Copyright 2019 Axel Huebl, David Grote, Maxence Thevenet
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef CHARGEDEPOSITION_H_
#define CHARGEDEPOSITION_H_

#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/ShapeFactors.H"


/* \brief Charge Deposition for thread thread_num
 * /param GetPosition : A functor for returning the particle position.
 * \param wp           : Pointer to array of particle weights.
 * \param ion_lev      : Pointer to array of particle ionization level. This is
                         required to have the charge of each macroparticle
                         since q is a scalar. For non-ionizable species,
                         ion_lev is a null pointer.
 * \param rho_arr      : Array4 of charge density, either full array or tile.
 * \param np_to_depose : Number of particles for which current is deposited.
 * \param dx           : 3D cell size
 * \param xyzmin       : Physical lower bounds of domain.
 * \param lo           : Index lower bounds of domain.
 * /param q            : species charge.
 */
template <int depos_order>
void doChargeDepositionShapeN(const GetParticlePosition& GetPosition,
                              const amrex::ParticleReal * const wp,
                              const int * const ion_lev,
                              const amrex::Array4<amrex::Real>& rho_arr,
                              const long np_to_depose,
                              const std::array<amrex::Real,3>& dx,
                              const std::array<amrex::Real, 3> xyzmin,
                              const amrex::Dim3 lo,
                              const amrex::Real q)
{
    // Whether ion_lev is a null pointer (do_ionization=0) or a real pointer
    // (do_ionization=1)
    const bool do_ionization = ion_lev;
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
            amrex::Real wq = q*wp[ip]*invvol;
            if (do_ionization){
                wq *= ion_lev[ip];
            }

            amrex::ParticleReal xp, yp, zp;
            GetPosition(ip, xp, yp, zp);

            // --- Compute shape factors
            // x direction
            // Get particle position in grid coordinates
#if (defined WARPX_DIM_RZ)
            const amrex::Real r = std::sqrt(xp*xp + yp*yp);
            const amrex::Real x = (r - xmin)*dxi;
#else
            const amrex::Real x = (xp - xmin)*dxi;
#endif
            // Compute shape factors for node-centered quantities
            amrex::Real sx[depos_order + 1];
            // i: leftmost grid point (node-centered) that the particle touches
            const int i = compute_shape_factor<depos_order>(sx,  x);

#if (defined WARPX_DIM_3D)
            // y direction
            const amrex::Real y = (yp - ymin)*dyi;
            amrex::Real sy[depos_order + 1];
            const int j = compute_shape_factor<depos_order>(sy,  y);
#endif
            // z direction
            const amrex::Real z = (zp - zmin)*dzi;
            amrex::Real sz[depos_order + 1];
            const int k = compute_shape_factor<depos_order>(sz,  z);

            // Deposit charge into rho_arr
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
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
