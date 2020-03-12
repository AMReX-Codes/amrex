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
 * \param rho_fab      : FArrayBox of charge density, either full array or tile.
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
                              amrex::FArrayBox& rho_fab,
                              const long np_to_depose,
                              const std::array<amrex::Real,3>& dx,
                              const std::array<amrex::Real, 3> xyzmin,
                              const amrex::Dim3 lo,
                              const amrex::Real q,
                              const long n_rz_azimuthal_modes)
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

    amrex::Array4<amrex::Real> const& rho_arr = rho_fab.array();
    amrex::IntVect const rho_type = rho_fab.box().type();

    constexpr int zdir = (AMREX_SPACEDIM - 1);
    constexpr int NODE = amrex::IndexType::NODE;
    constexpr int CELL = amrex::IndexType::CELL;

    // Loop over particles and deposit into rho_fab
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
            const amrex::Real rp = std::sqrt(xp*xp + yp*yp);
            amrex::Real costheta;
            amrex::Real sintheta;
            if (rp > 0.) {
                costheta = xp/rp;
                sintheta = yp/rp;
            } else {
                costheta = 1.;
                sintheta = 0.;
            }
            const Complex xy0 = Complex{costheta, sintheta};
            const amrex::Real x = (rp - xmin)*dxi;
#else
            const amrex::Real x = (xp - xmin)*dxi;
#endif

            // Compute shape factor along x
            // i: leftmost grid point that the particle touches
            amrex::Real sx[depos_order + 1];
            int i;
            if (rho_type[0] == NODE) {
                i = compute_shape_factor<depos_order>(sx, x);
            } else if (rho_type[0] == CELL) {
                i = compute_shape_factor<depos_order>(sx, x - 0.5);
            }

#if (defined WARPX_DIM_3D)
            // y direction
            const amrex::Real y = (yp - ymin)*dyi;
            amrex::Real sy[depos_order + 1];
            int j;
            if (rho_type[1] == NODE) {
                j = compute_shape_factor<depos_order>(sy, y);
            } else if (rho_type[1] == CELL) {
                j = compute_shape_factor<depos_order>(sy, y - 0.5);
            }
#endif
            // z direction
            const amrex::Real z = (zp - zmin)*dzi;
            amrex::Real sz[depos_order + 1];
            int k;
            if (rho_type[zdir] == NODE) {
                k = compute_shape_factor<depos_order>(sz, z);
            } else if (rho_type[zdir] == CELL) {
                k = compute_shape_factor<depos_order>(sz, z - 0.5);
            }

            // Deposit charge into rho_arr
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
            for (int iz=0; iz<=depos_order; iz++){
                for (int ix=0; ix<=depos_order; ix++){
                    amrex::Gpu::Atomic::Add(
                        &rho_arr(lo.x+i+ix, lo.y+k+iz, 0, 0),
                        sx[ix]*sz[iz]*wq);
#if (defined WARPX_DIM_RZ)
                    Complex xy = xy0; // Throughout the following loop, xy takes the value e^{i m theta}
                    for (int imode=1 ; imode < n_rz_azimuthal_modes ; imode++) {
                        amrex::Gpu::Atomic::Add( &rho_arr(lo.x+i+ix, lo.y+k+iz, 0, 2*imode-1), sx[ix]*sz[iz]*wq*xy.real());
                        amrex::Gpu::Atomic::Add( &rho_arr(lo.x+i+ix, lo.y+k+iz, 0, 2*imode  ), sx[ix]*sz[iz]*wq*xy.imag());
                        xy = xy*xy0;
                    }
#endif
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
