/* Copyright 2019 Axel Huebl, David Grote, Maxence Thevenet
 * Remi Lehe, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef CURRENTDEPOSITION_H_
#define CURRENTDEPOSITION_H_

#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/ShapeFactors.H"
#include "Utils/WarpX_Complex.H"

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

/**
 * \brief Current Deposition for thread thread_num
 * /param GetPosition : A functor for returning the particle position.
 * \param wp           : Pointer to array of particle weights.
 * \param uxp uyp uzp  : Pointer to arrays of particle momentum.
 * \param ion_lev      : Pointer to array of particle ionization level. This is
                         required to have the charge of each macroparticle
                         since q is a scalar. For non-ionizable species,
                         ion_lev is a null pointer.
 * \param jx_fab       : FArrayBox of current density, either full array or tile.
 * \param jy_fab       : FArrayBox of current density, either full array or tile.
 * \param jz_fab       : FArrayBox of current density, either full array or tile.
 * \param np_to_depose : Number of particles for which current is deposited.
 * \param dt           : Time step for particle level
 * \param dx           : 3D cell size
 * \param xyzmin       : Physical lower bounds of domain.
 * \param lo           : Index lower bounds of domain.
 * /param q            : species charge.
 */
template <int depos_order>
void doDepositionShapeN(const GetParticlePosition& GetPosition,
                        const amrex::ParticleReal * const wp,
                        const amrex::ParticleReal * const uxp,
                        const amrex::ParticleReal * const uyp,
                        const amrex::ParticleReal * const uzp,
                        const int * const ion_lev,
                        amrex::FArrayBox& jx_fab,
                        amrex::FArrayBox& jy_fab,
                        amrex::FArrayBox& jz_fab,
                        const long np_to_depose, const amrex::Real dt,
                        const std::array<amrex::Real,3>& dx,
                        const std::array<amrex::Real,3>& xyzmin,
                        const amrex::Dim3 lo,
                        const amrex::Real q,
                        const long n_rz_azimuthal_modes)
{
    // Whether ion_lev is a null pointer (do_ionization=0) or a real pointer
    // (do_ionization=1)
    const bool do_ionization = ion_lev;
    const amrex::Real dxi = 1.0/dx[0];
    const amrex::Real dzi = 1.0/dx[2];
    const amrex::Real dts2dx = 0.5*dt*dxi;
    const amrex::Real dts2dz = 0.5*dt*dzi;
#if (AMREX_SPACEDIM == 2)
    const amrex::Real invvol = dxi*dzi;
#elif (defined WARPX_DIM_3D)
    const amrex::Real dyi = 1.0/dx[1];
    const amrex::Real dts2dy = 0.5*dt*dyi;
    const amrex::Real invvol = dxi*dyi*dzi;
#endif

    const amrex::Real xmin = xyzmin[0];
    const amrex::Real ymin = xyzmin[1];
    const amrex::Real zmin = xyzmin[2];

    const amrex::Real clightsq = 1.0/PhysConst::c/PhysConst::c;

    amrex::Array4<amrex::Real> const& jx_arr = jx_fab.array();
    amrex::Array4<amrex::Real> const& jy_arr = jy_fab.array();
    amrex::Array4<amrex::Real> const& jz_arr = jz_fab.array();
    amrex::IntVect const jx_type = jx_fab.box().type();
    amrex::IntVect const jy_type = jy_fab.box().type();
    amrex::IntVect const jz_type = jz_fab.box().type();

    constexpr int zdir = (AMREX_SPACEDIM - 1);
    constexpr int NODE = amrex::IndexType::NODE;
    constexpr int CELL = amrex::IndexType::CELL;

    // Loop over particles and deposit into jx_fab, jy_fab and jz_fab
    amrex::ParallelFor(
        np_to_depose,
        [=] AMREX_GPU_DEVICE (long ip) {
            // --- Get particle quantities
            const amrex::Real gaminv = 1.0/std::sqrt(1.0 + uxp[ip]*uxp[ip]*clightsq
                                                         + uyp[ip]*uyp[ip]*clightsq
                                                         + uzp[ip]*uzp[ip]*clightsq);
            amrex::Real wq  = q*wp[ip];
            if (do_ionization){
                wq *= ion_lev[ip];
            }

            amrex::ParticleReal xp, yp, zp;
            GetPosition(ip, xp, yp, zp);

            const amrex::Real vx  = uxp[ip]*gaminv;
            const amrex::Real vy  = uyp[ip]*gaminv;
            const amrex::Real vz  = uzp[ip]*gaminv;
            // wqx, wqy wqz are particle current in each direction
#if (defined WARPX_DIM_RZ)
            // In RZ, wqx is actually wqr, and wqy is wqtheta
            // Convert to cylinderical at the mid point
            const amrex::Real xpmid = xp - 0.5*dt*vx;
            const amrex::Real ypmid = yp - 0.5*dt*vy;
            const amrex::Real rpmid = std::sqrt(xpmid*xpmid + ypmid*ypmid);
            amrex::Real costheta;
            amrex::Real sintheta;
            if (rpmid > 0.) {
                costheta = xpmid/rpmid;
                sintheta = ypmid/rpmid;
            } else {
                costheta = 1.;
                sintheta = 0.;
            }
            const Complex xy0 = Complex{costheta, sintheta};
            const amrex::Real wqx = wq*invvol*(+vx*costheta + vy*sintheta);
            const amrex::Real wqy = wq*invvol*(-vx*sintheta + vy*costheta);
#else
            const amrex::Real wqx = wq*invvol*vx;
            const amrex::Real wqy = wq*invvol*vy;
#endif
            const amrex::Real wqz = wq*invvol*vz;

            // --- Compute shape factors
            // x direction
            // Get particle position after 1/2 push back in position
#if (defined WARPX_DIM_RZ)
            const amrex::Real xmid = (rpmid - xmin)*dxi;
#else
            const amrex::Real xmid = (xp - xmin)*dxi - dts2dx*vx;
#endif
            // j_j[xyz] leftmost grid point in x that the particle touches for the centering of each current
            // sx_j[xyz] shape factor along x for the centering of each current
            // There are only two possible centerings, node or cell centered, so at most only two shape factor
            // arrays will be needed.
            amrex::Real sx_node[depos_order + 1];
            amrex::Real sx_cell[depos_order + 1];
            int j_node;
            int j_cell;
            if (jx_type[0] == NODE || jy_type[0] == NODE || jz_type[0] == NODE) {
                j_node = compute_shape_factor<depos_order>(sx_node, xmid);
            }
            if (jx_type[0] == CELL || jy_type[0] == CELL || jz_type[0] == CELL) {
                j_cell = compute_shape_factor<depos_order>(sx_cell, xmid - 0.5);
            }
            const amrex::Real (&sx_jx)[depos_order + 1] = ((jx_type[0] == NODE) ? sx_node : sx_cell);
            const amrex::Real (&sx_jy)[depos_order + 1] = ((jy_type[0] == NODE) ? sx_node : sx_cell);
            const amrex::Real (&sx_jz)[depos_order + 1] = ((jz_type[0] == NODE) ? sx_node : sx_cell);
            int const j_jx = ((jx_type[0] == NODE) ? j_node : j_cell);
            int const j_jy = ((jy_type[0] == NODE) ? j_node : j_cell);
            int const j_jz = ((jz_type[0] == NODE) ? j_node : j_cell);

#if (defined WARPX_DIM_3D)
            // y direction
            const amrex::Real ymid = (yp - ymin)*dyi - dts2dy*vy;
            amrex::Real sy_node[depos_order + 1];
            amrex::Real sy_cell[depos_order + 1];
            int k_node;
            int k_cell;
            if (jx_type[1] == NODE || jy_type[1] == NODE || jz_type[1] == NODE) {
                k_node = compute_shape_factor<depos_order>(sy_node, ymid);
            }
            if (jx_type[1] == CELL || jy_type[1] == CELL || jz_type[1] == CELL) {
                k_cell = compute_shape_factor<depos_order>(sy_cell, ymid - 0.5);
            }
            const amrex::Real (&sy_jx)[depos_order + 1] = ((jx_type[1] == NODE) ? sy_node : sy_cell);
            const amrex::Real (&sy_jy)[depos_order + 1] = ((jy_type[1] == NODE) ? sy_node : sy_cell);
            const amrex::Real (&sy_jz)[depos_order + 1] = ((jz_type[1] == NODE) ? sy_node : sy_cell);
            int const k_jx = ((jx_type[1] == NODE) ? k_node : k_cell);
            int const k_jy = ((jy_type[1] == NODE) ? k_node : k_cell);
            int const k_jz = ((jz_type[1] == NODE) ? k_node : k_cell);
#endif

            // z direction
            const amrex::Real zmid = (zp - zmin)*dzi - dts2dz*vz;
            amrex::Real sz_node[depos_order + 1];
            amrex::Real sz_cell[depos_order + 1];
            int l_node;
            int l_cell;
            if (jx_type[zdir] == NODE || jy_type[zdir] == NODE || jz_type[zdir] == NODE) {
                l_node = compute_shape_factor<depos_order>(sz_node, zmid);
            }
            if (jx_type[zdir] == CELL || jy_type[zdir] == CELL || jz_type[zdir] == CELL) {
                l_cell = compute_shape_factor<depos_order>(sz_cell, zmid - 0.5);
            }
            const amrex::Real (&sz_jx)[depos_order + 1] = ((jx_type[zdir] == NODE) ? sz_node : sz_cell);
            const amrex::Real (&sz_jy)[depos_order + 1] = ((jy_type[zdir] == NODE) ? sz_node : sz_cell);
            const amrex::Real (&sz_jz)[depos_order + 1] = ((jz_type[zdir] == NODE) ? sz_node : sz_cell);
            int const l_jx = ((jx_type[zdir] == NODE) ? l_node : l_cell);
            int const l_jy = ((jy_type[zdir] == NODE) ? l_node : l_cell);
            int const l_jz = ((jz_type[zdir] == NODE) ? l_node : l_cell);

            // Deposit current into jx_arr, jy_arr and jz_arr
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
            for (int iz=0; iz<=depos_order; iz++){
                for (int ix=0; ix<=depos_order; ix++){
                    amrex::Gpu::Atomic::Add(
                        &jx_arr(lo.x+j_jx+ix, lo.y+l_jx+iz, 0, 0),
                        sx_jx[ix]*sz_jx[iz]*wqx);
                    amrex::Gpu::Atomic::Add(
                        &jy_arr(lo.x+j_jy+ix, lo.y+l_jy+iz, 0, 0),
                        sx_jy[ix]*sz_jy[iz]*wqy);
                    amrex::Gpu::Atomic::Add(
                        &jz_arr(lo.x+j_jz+ix, lo.y+l_jz+iz, 0, 0),
                        sx_jz[ix]*sz_jz[iz]*wqz);
#if (defined WARPX_DIM_RZ)
                    Complex xy = xy0; // Note that xy is equal to e^{i m theta}
                    for (int imode=1 ; imode < n_rz_azimuthal_modes ; imode++) {
                        // The factor 2 on the weighting comes from the normalization of the modes
                        amrex::Gpu::Atomic::Add( &jx_arr(lo.x+j_jx+ix, lo.y+l_jx+iz, 0, 2*imode-1), 2.*sx_jx[ix]*sz_jx[iz]*wqx*xy.real());
                        amrex::Gpu::Atomic::Add( &jx_arr(lo.x+j_jx+ix, lo.y+l_jx+iz, 0, 2*imode  ), 2.*sx_jx[ix]*sz_jx[iz]*wqx*xy.imag());
                        amrex::Gpu::Atomic::Add( &jy_arr(lo.x+j_jy+ix, lo.y+l_jy+iz, 0, 2*imode-1), 2.*sx_jy[ix]*sz_jy[iz]*wqy*xy.real());
                        amrex::Gpu::Atomic::Add( &jy_arr(lo.x+j_jy+ix, lo.y+l_jy+iz, 0, 2*imode  ), 2.*sx_jy[ix]*sz_jy[iz]*wqy*xy.imag());
                        amrex::Gpu::Atomic::Add( &jz_arr(lo.x+j_jz+ix, lo.y+l_jz+iz, 0, 2*imode-1), 2.*sx_jz[ix]*sz_jz[iz]*wqz*xy.real());
                        amrex::Gpu::Atomic::Add( &jz_arr(lo.x+j_jz+ix, lo.y+l_jz+iz, 0, 2*imode  ), 2.*sx_jz[ix]*sz_jz[iz]*wqz*xy.imag());
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
                            &jx_arr(lo.x+j_jx+ix, lo.y+k_jx+iy, lo.z+l_jx+iz),
                            sx_jx[ix]*sy_jx[iy]*sz_jx[iz]*wqx);
                        amrex::Gpu::Atomic::Add(
                            &jy_arr(lo.x+j_jy+ix, lo.y+k_jy+iy, lo.z+l_jy+iz),
                            sx_jy[ix]*sy_jy[iy]*sz_jy[iz]*wqy);
                        amrex::Gpu::Atomic::Add(
                            &jz_arr(lo.x+j_jz+ix, lo.y+k_jz+iy, lo.z+l_jz+iz),
                            sx_jz[ix]*sy_jz[iy]*sz_jz[iz]*wqz);
                    }
                }
            }
#endif
        }
        );
}

/**
 * \brief Esirkepov Current Deposition for thread thread_num
 *
 * /param GetPosition : A functor for returning the particle position.
 * \param wp           : Pointer to array of particle weights.
 * \param uxp uyp uzp  : Pointer to arrays of particle momentum.
 * \param ion_lev      : Pointer to array of particle ionization level. This is
                         required to have the charge of each macroparticle
                         since q is a scalar. For non-ionizable species,
                         ion_lev is a null pointer.
 * \param Jx_arr       : Array4 of current density, either full array or tile.
 * \param Jy_arr       : Array4 of current density, either full array or tile.
 * \param Jz_arr       : Array4 of current density, either full array or tile.
 * \param np_to_depose : Number of particles for which current is deposited.
 * \param dt           : Time step for particle level
 * \param dx           : 3D cell size
 * \param xyzmin       : Physical lower bounds of domain.
 * \param lo           : Index lower bounds of domain.
 * \param q            : species charge.
 * \param n_rz_azimuthal_modes: Number of azimuthal modes when using RZ geometry
 */
template <int depos_order>
void doEsirkepovDepositionShapeN (const GetParticlePosition& GetPosition,
                                  const amrex::ParticleReal * const wp,
                                  const amrex::ParticleReal * const uxp,
                                  const amrex::ParticleReal * const uyp,
                                  const amrex::ParticleReal * const uzp,
                                  const int * ion_lev,
                                  const amrex::Array4<amrex::Real>& Jx_arr,
                                  const amrex::Array4<amrex::Real>& Jy_arr,
                                  const amrex::Array4<amrex::Real>& Jz_arr,
                                  const long np_to_depose,
                                  const amrex::Real dt,
                                  const std::array<amrex::Real,3>& dx,
                                  const std::array<amrex::Real, 3> xyzmin,
                                  const amrex::Dim3 lo,
                                  const amrex::Real q,
                                  const long n_rz_azimuthal_modes)
{
    using namespace amrex;

    // Whether ion_lev is a null pointer (do_ionization=0) or a real pointer
    // (do_ionization=1)
    bool const do_ionization = ion_lev;
    Real const dxi = 1.0_rt / dx[0];
    Real const dtsdx0 = dt*dxi;
    Real const xmin = xyzmin[0];
#if (defined WARPX_DIM_3D)
    Real const dyi = 1.0_rt / dx[1];
    Real const dtsdy0 = dt*dyi;
    Real const ymin = xyzmin[1];
#endif
    Real const dzi = 1.0_rt / dx[2];
    Real const dtsdz0 = dt*dzi;
    Real const zmin = xyzmin[2];

#if (defined WARPX_DIM_3D)
    Real const invdtdx = 1.0_rt / (dt*dx[1]*dx[2]);
    Real const invdtdy = 1.0_rt / (dt*dx[0]*dx[2]);
    Real const invdtdz = 1.0_rt / (dt*dx[0]*dx[1]);
#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
    Real const invdtdx = 1.0_rt / (dt*dx[2]);
    Real const invdtdz = 1.0_rt / (dt*dx[0]);
    Real const invvol = 1.0_rt / (dx[0]*dx[2]);
#endif

#if (defined WARPX_DIM_RZ)
    Complex const I = Complex{0., 1.};
#endif

    Real const clightsq = 1.0_rt / ( PhysConst::c * PhysConst::c );

    // Loop over particles and deposit into Jx_arr, Jy_arr and Jz_arr
    amrex::ParallelFor(
        np_to_depose,
        [=] AMREX_GPU_DEVICE (long const ip) {

            // --- Get particle quantities
        Real const gaminv = 1.0/std::sqrt(1.0 + uxp[ip]*uxp[ip]*clightsq
                                              + uyp[ip]*uyp[ip]*clightsq
                                              + uzp[ip]*uzp[ip]*clightsq);

            // wqx, wqy wqz are particle current in each direction
            Real wq = q*wp[ip];
            if (do_ionization){
                wq *= ion_lev[ip];
            }

            ParticleReal xp, yp, zp;
            GetPosition(ip, xp, yp, zp);

            Real const wqx = wq*invdtdx;
#if (defined WARPX_DIM_3D)
            Real const wqy = wq*invdtdy;
#endif
            Real const wqz = wq*invdtdz;

            // computes current and old position in grid units
#if (defined WARPX_DIM_RZ)
            Real const xp_mid = xp - 0.5_rt * dt*uxp[ip]*gaminv;
            Real const yp_mid = yp - 0.5_rt * dt*uyp[ip]*gaminv;
            Real const xp_old = xp - dt*uxp[ip]*gaminv;
            Real const yp_old = yp - dt*uyp[ip]*gaminv;
            Real const rp_new = std::sqrt(xp*xp
                                        + yp*yp);
            Real const rp_mid = std::sqrt(xp_mid*xp_mid + yp_mid*yp_mid);
            Real const rp_old = std::sqrt(xp_old*xp_old + yp_old*yp_old);
            Real costheta_new, sintheta_new;
            if (rp_new > 0._rt) {
                costheta_new = xp/rp_new;
                sintheta_new = yp/rp_new;
            } else {
                costheta_new = 1.;
                sintheta_new = 0.;
            }
            amrex::Real costheta_mid, sintheta_mid;
            if (rp_mid > 0._rt) {
                costheta_mid = xp_mid/rp_mid;
                sintheta_mid = yp_mid/rp_mid;
            } else {
                costheta_mid = 1.;
                sintheta_mid = 0.;
            }
            amrex::Real costheta_old, sintheta_old;
            if (rp_old > 0._rt) {
                costheta_old = xp_old/rp_old;
                sintheta_old = yp_old/rp_old;
            } else {
                costheta_old = 1.;
                sintheta_old = 0.;
            }
            const Complex xy_new0 = Complex{costheta_new, sintheta_new};
            const Complex xy_mid0 = Complex{costheta_mid, sintheta_mid};
            const Complex xy_old0 = Complex{costheta_old, sintheta_old};
            Real const x_new = (rp_new - xmin)*dxi;
            Real const x_old = (rp_old - xmin)*dxi;
#else
            Real const x_new = (xp - xmin)*dxi;
            Real const x_old = x_new - dtsdx0*uxp[ip]*gaminv;
#endif
#if (defined WARPX_DIM_3D)
            Real const y_new = (yp - ymin)*dyi;
            Real const y_old = y_new - dtsdy0*uyp[ip]*gaminv;
#endif
            Real const z_new = (zp - zmin)*dzi;
            Real const z_old = z_new - dtsdz0*uzp[ip]*gaminv;

#if (defined WARPX_DIM_RZ)
            Real const vy = (-uxp[ip]*sintheta_mid + uyp[ip]*costheta_mid)*gaminv;
#elif (defined WARPX_DIM_XZ)
            Real const vy = uyp[ip]*gaminv;
#endif

            // Shape factor arrays
            // Note that there are extra values above and below
            // to possibly hold the factor for the old particle
            // which can be at a different grid location.
            Real sx_new[depos_order + 3] = {0.};
            Real sx_old[depos_order + 3] = {0.};
#if (defined WARPX_DIM_3D)
            Real sy_new[depos_order + 3] = {0.};
            Real sy_old[depos_order + 3] = {0.};
#endif
            Real sz_new[depos_order + 3] = {0.};
            Real sz_old[depos_order + 3] = {0.};

            // --- Compute shape factors
            // Compute shape factors for position as they are now and at old positions
            // [ijk]_new: leftmost grid point that the particle touches
            const int i_new = compute_shape_factor<depos_order>(sx_new+1, x_new);
            const int i_old = compute_shifted_shape_factor<depos_order>(sx_old, x_old, i_new);
#if (defined WARPX_DIM_3D)
            const int j_new = compute_shape_factor<depos_order>(sy_new+1, y_new);
            const int j_old = compute_shifted_shape_factor<depos_order>(sy_old, y_old, j_new);
#endif
            const int k_new = compute_shape_factor<depos_order>(sz_new+1, z_new);
            const int k_old = compute_shifted_shape_factor<depos_order>(sz_old, z_old, k_new);

            // computes min/max positions of current contributions
            int dil = 1, diu = 1;
            if (i_old < i_new) dil = 0;
            if (i_old > i_new) diu = 0;
#if (defined WARPX_DIM_3D)
            int djl = 1, dju = 1;
            if (j_old < j_new) djl = 0;
            if (j_old > j_new) dju = 0;
#endif
            int dkl = 1, dku = 1;
            if (k_old < k_new) dkl = 0;
            if (k_old > k_new) dku = 0;

#if (defined WARPX_DIM_3D)

            for (int k=dkl; k<=depos_order+2-dku; k++) {
                for (int j=djl; j<=depos_order+2-dju; j++) {
                    amrex::Real sdxi = 0.;
                    for (int i=dil; i<=depos_order+1-diu; i++) {
                        sdxi += wqx*(sx_old[i] - sx_new[i])*((sy_new[j] + 0.5*(sy_old[j] - sy_new[j]))*sz_new[k] +
                                                             (0.5*sy_new[j] + 1./3.*(sy_old[j] - sy_new[j]))*(sz_old[k] - sz_new[k]));
                        amrex::Gpu::Atomic::Add( &Jx_arr(lo.x+i_new-1+i, lo.y+j_new-1+j, lo.z+k_new-1+k), sdxi);
                    }
                }
            }
            for (int k=dkl; k<=depos_order+2-dku; k++) {
                for (int i=dil; i<=depos_order+2-diu; i++) {
                    amrex::Real sdyj = 0.;
                    for (int j=djl; j<=depos_order+1-dju; j++) {
                        sdyj += wqy*(sy_old[j] - sy_new[j])*((sz_new[k] + 0.5*(sz_old[k] - sz_new[k]))*sx_new[i] +
                                                             (0.5*sz_new[k] + 1./3.*(sz_old[k] - sz_new[k]))*(sx_old[i] - sx_new[i]));
                        amrex::Gpu::Atomic::Add( &Jy_arr(lo.x+i_new-1+i, lo.y+j_new-1+j, lo.z+k_new-1+k), sdyj);
                    }
                }
            }
            for (int j=djl; j<=depos_order+2-dju; j++) {
                for (int i=dil; i<=depos_order+2-diu; i++) {
                    amrex::Real sdzk = 0.;
                    for (int k=dkl; k<=depos_order+1-dku; k++) {
                        sdzk += wqz*(sz_old[k] - sz_new[k])*((sx_new[i] + 0.5*(sx_old[i] - sx_new[i]))*sy_new[j] +
                                                             (0.5*sx_new[i] + 1./3.*(sx_old[i] - sx_new[i]))*(sy_old[j] - sy_new[j]));
                        amrex::Gpu::Atomic::Add( &Jz_arr(lo.x+i_new-1+i, lo.y+j_new-1+j, lo.z+k_new-1+k), sdzk);
                    }
                }
            }

#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)

            for (int k=dkl; k<=depos_order+2-dku; k++) {
                amrex::Real sdxi = 0.;
                for (int i=dil; i<=depos_order+1-diu; i++) {
                    sdxi += wqx*(sx_old[i] - sx_new[i])*(sz_new[k] + 0.5*(sz_old[k] - sz_new[k]));
                    amrex::Gpu::Atomic::Add( &Jx_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 0), sdxi);
#if (defined WARPX_DIM_RZ)
                    Complex xy_mid = xy_mid0; // Throughout the following loop, xy_mid takes the value e^{i m theta}
                    for (int imode=1 ; imode < n_rz_azimuthal_modes ; imode++) {
                        // The factor 2 comes from the normalization of the modes
                        const Complex djr_cmplx = 2._rt *sdxi*xy_mid;
                        amrex::Gpu::Atomic::Add( &Jx_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 2*imode-1), djr_cmplx.real());
                        amrex::Gpu::Atomic::Add( &Jx_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 2*imode), djr_cmplx.imag());
                        xy_mid = xy_mid*xy_mid0;
                    }
#endif
                }
            }
            for (int k=dkl; k<=depos_order+2-dku; k++) {
                for (int i=dil; i<=depos_order+2-diu; i++) {
                    Real const sdyj = wq*vy*invvol*((sz_new[k] + 0.5_rt * (sz_old[k] - sz_new[k]))*sx_new[i] +
                                                           (0.5_rt * sz_new[k] + 1._rt / 3._rt *(sz_old[k] - sz_new[k]))*(sx_old[i] - sx_new[i]));
                    amrex::Gpu::Atomic::Add( &Jy_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 0), sdyj);
#if (defined WARPX_DIM_RZ)
                    Complex xy_new = xy_new0;
                    Complex xy_mid = xy_mid0;
                    Complex xy_old = xy_old0;
                    // Throughout the following loop, xy_ takes the value e^{i m theta_}
                    for (int imode=1 ; imode < n_rz_azimuthal_modes ; imode++) {
                        // The factor 2 comes from the normalization of the modes
                        // The minus sign comes from the different convention with respect to Davidson et al.
                        const Complex djt_cmplx = -2._rt * I*(i_new-1 + i + xmin*dxi)*wq*invdtdx/(amrex::Real)imode*
                                                  (sx_new[i]*sz_new[k]*(xy_new - xy_mid) + sx_old[i]*sz_old[k]*(xy_mid - xy_old));
                        amrex::Gpu::Atomic::Add( &Jy_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 2*imode-1), djt_cmplx.real());
                        amrex::Gpu::Atomic::Add( &Jy_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 2*imode), djt_cmplx.imag());
                        xy_new = xy_new*xy_new0;
                        xy_mid = xy_mid*xy_mid0;
                        xy_old = xy_old*xy_old0;
                    }
#endif
                }
            }
            for (int i=dil; i<=depos_order+2-diu; i++) {
                Real sdzk = 0.;
                for (int k=dkl; k<=depos_order+1-dku; k++) {
                    sdzk += wqz*(sz_old[k] - sz_new[k])*(sx_new[i] + 0.5_rt * (sx_old[i] - sx_new[i]));
                    amrex::Gpu::Atomic::Add( &Jz_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 0), sdzk);
#if (defined WARPX_DIM_RZ)
                    Complex xy_mid = xy_mid0; // Throughout the following loop, xy_mid takes the value e^{i m theta}
                    for (int imode=1 ; imode < n_rz_azimuthal_modes ; imode++) {
                        // The factor 2 comes from the normalization of the modes
                        const Complex djz_cmplx = 2._rt * sdzk * xy_mid;
                        amrex::Gpu::Atomic::Add( &Jz_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 2*imode-1), djz_cmplx.real());
                        amrex::Gpu::Atomic::Add( &Jz_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0, 2*imode), djz_cmplx.imag());
                        xy_mid = xy_mid*xy_mid0;
                    }
#endif
                }
            }

#endif
        }
        );
}

#endif // CURRENTDEPOSITION_H_
