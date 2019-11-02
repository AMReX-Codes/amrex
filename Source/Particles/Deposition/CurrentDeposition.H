#ifndef CURRENTDEPOSITION_H_
#define CURRENTDEPOSITION_H_

#include "ShapeFactors.H"
#include <WarpX_Complex.H>

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

/* \brief Current Deposition for thread thread_num
 * /param xp, yp, zp   : Pointer to arrays of particle positions.
 * \param wp           : Pointer to array of particle weights.
 * \param uxp uyp uzp  : Pointer to arrays of particle momentum.
 * \param ion_lev      : Pointer to array of particle ionization level. This is
                         required to have the charge of each macroparticle
                         since q is a scalar. For non-ionizable species,
                         ion_lev is a null pointer.
 * \param jx_arr       : Array4 of current density, either full array or tile.
 * \param jy_arr       : Array4 of current density, either full array or tile.
 * \param jz_arr       : Array4 of current density, either full array or tile.
 * \param np_to_depose : Number of particles for which current is deposited.
 * \param dt           : Time step for particle level
 * \param dx           : 3D cell size
 * \param tbx          : Box of domain for jx.
 * \param tby          : Box of domain for jy.
 * \param tbz          : Box of domain for jz.
 * \param depos_lev    : Refinement level being used for the deposition.
 * \param lo           : Index lower bounds of domain.
 * /param q            : species charge.
 */
template <int depos_order>
void doDepositionShapeN(const amrex::ParticleReal * const xp,
                        const amrex::ParticleReal * const yp,
                        const amrex::ParticleReal * const zp,
                        const amrex::ParticleReal * const wp,
                        const amrex::ParticleReal * const uxp,
                        const amrex::ParticleReal * const uyp,
                        const amrex::ParticleReal * const uzp,
                        const int * const ion_lev,
                        const amrex::Array4<amrex::Real>& jx_arr,
                        const amrex::Array4<amrex::Real>& jy_arr,
                        const amrex::Array4<amrex::Real>& jz_arr,
                        const long np_to_depose, const amrex::Real dt,
                        const std::array<amrex::Real,3>& dx,
                        const amrex::Box& tbx,
                        const amrex::Box& tby,
                        const amrex::Box& tbz,
                        const int depos_lev,
                        const amrex::Dim3 lo,
                        const amrex::Real q)
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

    const std::array<amrex::Real, 3> xyzmin_jx = WarpX::LowerCornerWithCentering(tbx, depos_lev);
    const std::array<amrex::Real, 3> xyzmin_jy = WarpX::LowerCornerWithCentering(tby, depos_lev);
    const std::array<amrex::Real, 3> xyzmin_jz = WarpX::LowerCornerWithCentering(tbz, depos_lev);
    const amrex::Real xmin_jx = xyzmin_jx[0];
    const amrex::Real xmin_jy = xyzmin_jy[0];
    const amrex::Real xmin_jz = xyzmin_jz[0];
#if (defined WARPX_DIM_3D)
    const amrex::Real ymin_jx = xyzmin_jx[1];
    const amrex::Real ymin_jy = xyzmin_jy[1];
    const amrex::Real ymin_jz = xyzmin_jz[1];
#endif
    const amrex::Real zmin_jx = xyzmin_jx[2];
    const amrex::Real zmin_jy = xyzmin_jy[2];
    const amrex::Real zmin_jz = xyzmin_jz[2];

    const amrex::Real clightsq = 1.0/PhysConst::c/PhysConst::c;

    // Loop over particles and deposit into jx_arr, jy_arr and jz_arr
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
            const amrex::Real vx  = uxp[ip]*gaminv;
            const amrex::Real vy  = uyp[ip]*gaminv;
            const amrex::Real vz  = uzp[ip]*gaminv;
            // wqx, wqy wqz are particle current in each direction
#if (defined WARPX_DIM_RZ)
            // In RZ, wqx is actually wqr, and wqy is wqtheta
            // Convert to cylinderical at the mid point
            const amrex::Real xpmid = xp[ip] - 0.5*dt*vx;
            const amrex::Real ypmid = yp[ip] - 0.5*dt*vy;
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
            const amrex::Real xmid_jx = (rpmid - xmin_jx)*dxi;
            const amrex::Real xmid_jy = (rpmid - xmin_jy)*dxi;
            const amrex::Real xmid_jz = (rpmid - xmin_jz)*dxi;
#else
            const amrex::Real xmid_jx = (xp[ip] - xmin_jx)*dxi - dts2dx*vx;
            const amrex::Real xmid_jy = (xp[ip] - xmin_jy)*dxi - dts2dx*vx;
            const amrex::Real xmid_jz = (xp[ip] - xmin_jz)*dxi - dts2dx*vx;
#endif
            // j_[xyz]: leftmost grid point in x that the particle touches for the centering of each current
            // sx_[xyz] shape factor along x for the centering of each current
            // There are only two possible centerings, node or cell centered, so at most only two shape factor
            // arrays will be needed.
            // sx_jx is always calculated. sx_2 only if needed (when Jy or Jz have different centerings than Jx).
            amrex::Real sx_jx[depos_order + 1];
            amrex::Real sx_2[depos_order + 1];
            const int j_jx = compute_shape_factor<depos_order>(sx_jx, xmid_jx);
            // If Jy has the same centering as Jx, use the same index and shape factor, otherwise
            // calculate and use the different shape factor.
            const int j_jy = (tby.type(0) == tbx.type(0) ? j_jx : compute_shape_factor<depos_order>(sx_2, xmid_jy));
            const amrex::Real (&sx_jy)[depos_order + 1] = (tby.type(0) == tbx.type(0) ? sx_jx : sx_2);
            // If Jz has the same centering as Jx, use the same index and shape factor, otherwise
            // calculate and use the different shape factor, unless it has been already calculated for Jy.
            const int j_jz = (tbz.type(0) == tbx.type(0) ? j_jx :
                             (tbz.type(0) == tby.type(0) ? j_jy :
                                 compute_shape_factor<depos_order>(sx_2, xmid_jz)));
            const amrex::Real (&sx_jz)[depos_order + 1] = (tbz.type(0) == tbx.type(0) ? sx_jx : sx_2);

#if (defined WARPX_DIM_3D)
            // y direction
            const amrex::Real ymid_jx = (yp[ip] - ymin_jx)*dxi - dts2dx*vx;
            const amrex::Real ymid_jy = (yp[ip] - ymin_jy)*dxi - dts2dx*vx;
            const amrex::Real ymid_jz = (yp[ip] - ymin_jz)*dxi - dts2dx*vx;
            amrex::Real sy_jx[depos_order + 1];
            amrex::Real sy_2[depos_order + 1];
            const int k_jx = compute_shape_factor<depos_order>(sy_jx, ymid_jx);
            const int k_jy = (tby.type(1) == tbx.type(1) ? k_jx : compute_shape_factor<depos_order>(sy_2, ymid_jy));
            const amrex::Real (&sy_jy)[depos_order + 1] = (tby.type(1) == tbx.type(1) ? sy_jx : sy_2);
            const int k_jz = (tbz.type(1) == tbx.type(1) ? k_jx :
                             (tbz.type(1) == tby.type(1) ? k_jy :
                                 compute_shape_factor<depos_order>(sy_2, ymid_jz)));
            const amrex::Real (&sy_jz)[depos_order + 1] = (tbz.type(1) == tbx.type(1) ? sy_jx : sy_2);
#endif
            // z direction
            const amrex::Real zmid_jx = (zp[ip] - zmin_jx)*dxi - dts2dx*vx;
            const amrex::Real zmid_jy = (zp[ip] - zmin_jy)*dxi - dts2dx*vx;
            const amrex::Real zmid_jz = (zp[ip] - zmin_jz)*dxi - dts2dx*vx;
            amrex::Real sz_jx[depos_order + 1];
            amrex::Real sz_2[depos_order + 1];
            const int l_jx = compute_shape_factor<depos_order>(sz_jx, zmid_jx);
            const int l_jy = (tby.type(2) == tbx.type(2) ? l_jx : compute_shape_factor<depos_order>(sz_2, zmid_jy));
            const amrex::Real (&sz_jy)[depos_order + 1] = (tby.type(2) == tbx.type(2) ? sz_jx : sz_2);
            const int l_jz = (tbz.type(2) == tbx.type(2) ? l_jx :
                             (tbz.type(2) == tby.type(2) ? l_jy :
                                 compute_shape_factor<depos_order>(sz_2, zmid_jz)));
            const amrex::Real (&sz_jz)[depos_order + 1] = (tbz.type(2) == tbx.type(2) ? sz_jx : sz_2);

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

/* \brief Esirkepov Current Deposition for thread thread_num
 * \param xp, yp, zp   : Pointer to arrays of particle positions.
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
void doEsirkepovDepositionShapeN (const amrex::ParticleReal * const xp,
                                  const amrex::ParticleReal * const yp,
                                  const amrex::ParticleReal * const zp,
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
            Real const wqx = wq*invdtdx;
#if (defined WARPX_DIM_3D)
            Real const wqy = wq*invdtdy;
#endif
            Real const wqz = wq*invdtdz;

            // computes current and old position in grid units
#if (defined WARPX_DIM_RZ)
            Real const xp_mid = xp[ip] - 0.5_rt * dt*uxp[ip]*gaminv;
            Real const yp_mid = yp[ip] - 0.5_rt * dt*uyp[ip]*gaminv;
            Real const xp_old = xp[ip] - dt*uxp[ip]*gaminv;
            Real const yp_old = yp[ip] - dt*uyp[ip]*gaminv;
            Real const rp_new = std::sqrt(xp[ip]*xp[ip] + yp[ip]*yp[ip]);
            Real const rp_mid = std::sqrt(xp_mid*xp_mid + yp_mid*yp_mid);
            Real const rp_old = std::sqrt(xp_old*xp_old + yp_old*yp_old);
            Real costheta_new, sintheta_new;
            if (rp_new > 0._rt) {
                costheta_new = xp[ip]/rp_new;
                sintheta_new = yp[ip]/rp_new;
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
            Real const x_new = (xp[ip] - xmin)*dxi;
            Real const x_old = x_new - dtsdx0*uxp[ip]*gaminv;
#endif
#if (defined WARPX_DIM_3D)
            Real const y_new = (yp[ip] - ymin)*dyi;
            Real const y_old = y_new - dtsdy0*uyp[ip]*gaminv;
#endif
            Real const z_new = (zp[ip] - zmin)*dzi;
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
