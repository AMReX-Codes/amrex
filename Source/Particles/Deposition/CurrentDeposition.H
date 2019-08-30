#ifndef CURRENTDEPOSITION_H_
#define CURRENTDEPOSITION_H_

#include "ShapeFactors.H"

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
 * \param xyzmin       : Physical lower bounds of domain.
 * \param lo           : Index lower bounds of domain.
 * \param stagger_shift: 0 if nodal, 0.5 if staggered.
 * /param q            : species charge.
 */
template <int depos_order>
void doDepositionShapeN(const amrex::Real * const xp, 
                        const amrex::Real * const yp, 
                        const amrex::Real * const zp,
                        const amrex::Real * const wp,
                        const amrex::Real * const uxp,
                        const amrex::Real * const uyp,
                        const amrex::Real * const uzp,
                        const int * const ion_lev,
                        const amrex::Array4<amrex::Real>& jx_arr,
                        const amrex::Array4<amrex::Real>& jy_arr,
                        const amrex::Array4<amrex::Real>& jz_arr,
                        const long np_to_depose, const amrex::Real dt,
                        const std::array<amrex::Real,3>& dx,
                        const std::array<amrex::Real, 3> xyzmin,
                        const amrex::Dim3 lo,
                        const amrex::Real stagger_shift, 
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

    const amrex::Real xmin = xyzmin[0];
    const amrex::Real ymin = xyzmin[1];
    const amrex::Real zmin = xyzmin[2];
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
            const amrex::Real xmid = (rpmid-xmin)*dxi;
#else
            const amrex::Real xmid = (xp[ip]-xmin)*dxi-dts2dx*vx;
#endif
            // Compute shape factors for node-centered quantities
            amrex::Real sx [depos_order + 1];
            // j: leftmost grid point (node-centered) that the particle touches
            const int j  = compute_shape_factor<depos_order>(sx,  xmid);
            // Compute shape factors for cell-centered quantities
            amrex::Real sx0[depos_order + 1];
            // j0: leftmost grid point (cell-centered) that the particle touches
            const int j0 = compute_shape_factor<depos_order>(sx0, xmid-stagger_shift);
                     
#if (defined WARPX_DIM_3D)
            // y direction
            const amrex::Real ymid= (yp[ip]-ymin)*dyi-dts2dy*vy;
            amrex::Real sy [depos_order + 1];
            const int k  = compute_shape_factor<depos_order>(sy,  ymid);
            amrex::Real sy0[depos_order + 1];
            const int k0 = compute_shape_factor<depos_order>(sy0, ymid-stagger_shift);
#endif
            // z direction
            const amrex::Real zmid= (zp[ip]-zmin)*dzi-dts2dz*vz;
            amrex::Real sz [depos_order + 1];
            const int l  = compute_shape_factor<depos_order>(sz,  zmid);
            amrex::Real sz0[depos_order + 1];
            const int l0 = compute_shape_factor<depos_order>(sz0, zmid-stagger_shift);

            // Deposit current into jx_arr, jy_arr and jz_arr
#if (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
            for (int iz=0; iz<=depos_order; iz++){
                for (int ix=0; ix<=depos_order; ix++){
                    amrex::Gpu::Atomic::Add(
                        &jx_arr(lo.x+j0+ix, lo.y+l +iz, 0), 
                        sx0[ix]*sz [iz]*wqx);
                    amrex::Gpu::Atomic::Add(
                        &jy_arr(lo.x+j +ix, lo.y+l +iz, 0), 
                        sx [ix]*sz [iz]*wqy);
                    amrex::Gpu::Atomic::Add(
                        &jz_arr(lo.x+j +ix, lo.y+l0+iz, 0), 
                        sx [ix]*sz0[iz]*wqz);
                }
            }
#elif (defined WARPX_DIM_3D)
            for (int iz=0; iz<=depos_order; iz++){
                for (int iy=0; iy<=depos_order; iy++){
                    for (int ix=0; ix<=depos_order; ix++){
                        amrex::Gpu::Atomic::Add(
                            &jx_arr(lo.x+j0+ix, lo.y+k +iy, lo.z+l +iz),
                            sx0[ix]*sy [iy]*sz [iz]*wqx);
                        amrex::Gpu::Atomic::Add(
                            &jy_arr(lo.x+j +ix, lo.y+k0+iy, lo.z+l +iz), 
                            sx [ix]*sy0[iy]*sz [iz]*wqy);
                        amrex::Gpu::Atomic::Add(
                            &jz_arr(lo.x+j +ix, lo.y+k +iy, lo.z+l0+iz),
                            sx [ix]*sy [iy]*sz0[iz]*wqz);
                    }
                }
            }
#endif
        }
        );
}

/* \brief Esirkepov Current Deposition for thread thread_num
 * /param xp, yp, zp   : Pointer to arrays of particle positions.
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
 * /param q            : species charge.
 */
template <int depos_order>
void doEsirkepovDepositionShapeN (const amrex::Real * const xp,
                                  const amrex::Real * const yp,
                                  const amrex::Real * const zp,
                                  const amrex::Real * const wp,
                                  const amrex::Real * const uxp,
                                  const amrex::Real * const uyp,
                                  const amrex::Real * const uzp,
                                  const int * ion_lev,
                                  const amrex::Array4<amrex::Real>& Jx_arr,
                                  const amrex::Array4<amrex::Real>& Jy_arr,
                                  const amrex::Array4<amrex::Real>& Jz_arr,
                                  const long np_to_depose,
                                  const amrex::Real dt,
                                  const std::array<amrex::Real,3>& dx,
                                  const std::array<amrex::Real, 3> xyzmin,
                                  const amrex::Dim3 lo,
                                  const amrex::Real q)
{
    // Whether ion_lev is a null pointer (do_ionization=0) or a real pointer
    // (do_ionization=1)
    const bool do_ionization = ion_lev;
    const amrex::Real dxi = 1.0/dx[0];
    const amrex::Real dtsdx0 = dt*dxi;
    const amrex::Real xmin = xyzmin[0];
#if (defined WARPX_DIM_3D)
    const amrex::Real dyi = 1.0/dx[1];
    const amrex::Real dtsdy0 = dt*dyi;
    const amrex::Real ymin = xyzmin[1];
#endif
    const amrex::Real dzi = 1.0/dx[2];
    const amrex::Real dtsdz0 = dt*dzi;
    const amrex::Real zmin = xyzmin[2];

#if (defined WARPX_DIM_3D)
    const amrex::Real invdtdx = 1.0/(dt*dx[1]*dx[2]);
    const amrex::Real invdtdy = 1.0/(dt*dx[0]*dx[2]);
    const amrex::Real invdtdz = 1.0/(dt*dx[0]*dx[1]);
#elif (defined WARPX_DIM_XZ) || (defined WARPX_DIM_RZ)
    const amrex::Real invdtdx = 1.0/(dt*dx[2]);
    const amrex::Real invdtdz = 1.0/(dt*dx[0]);
    const amrex::Real invvol = 1.0/(dx[0]*dx[2]);
#endif

    const amrex::Real clightsq = 1.0/PhysConst::c/PhysConst::c;

    // Loop over particles and deposit into Jx_arr, Jy_arr and Jz_arr
    amrex::ParallelFor( 
        np_to_depose,
        [=] AMREX_GPU_DEVICE (long ip) {

            // --- Get particle quantities
            const amrex::Real gaminv = 1.0/std::sqrt(1.0 + uxp[ip]*uxp[ip]*clightsq
                                                         + uyp[ip]*uyp[ip]*clightsq
                                                         + uzp[ip]*uzp[ip]*clightsq);

            // wqx, wqy wqz are particle current in each direction
            amrex::Real wq = q*wp[ip];
            if (do_ionization){
                wq *= ion_lev[ip];
            }
            const amrex::Real wqx = wq*invdtdx;
#if (defined WARPX_DIM_3D)
            const amrex::Real wqy = wq*invdtdy;
#endif
            const amrex::Real wqz = wq*invdtdz;

            // computes current and old position in grid units
#if (defined WARPX_DIM_RZ)
            const amrex::Real r_new = std::sqrt(xp[ip]*xp[ip] + yp[ip]*yp[ip]);
            const amrex::Real r_old = std::sqrt((xp[ip] - dt*uxp[ip]*gaminv)*(xp[ip] - dt*uxp[ip]*gaminv) +
                                                (yp[ip] - dt*uyp[ip]*gaminv)*(yp[ip] - dt*uyp[ip]*gaminv));
            const amrex::Real x_new = (r_new - xmin)*dxi;
            const amrex::Real x_old = (r_old - xmin)*dxi;
#else
            const amrex::Real x_new = (xp[ip] - xmin)*dxi;
            const amrex::Real x_old = x_new - dtsdx0*uxp[ip]*gaminv;
#endif
#if (defined WARPX_DIM_3D)        
            const amrex::Real y_new = (yp[ip] - ymin)*dyi;
            const amrex::Real y_old = y_new - dtsdy0*uyp[ip]*gaminv;
#endif
            const amrex::Real z_new = (zp[ip] - zmin)*dzi;
            const amrex::Real z_old = z_new - dtsdz0*uzp[ip]*gaminv;

#if (defined WARPX_DIM_RZ)
            amrex::Real costheta;
            amrex::Real sintheta;
            if (r_new > 0.) {
                costheta = xp[ip]/r_new;
                sintheta = yp[ip]/r_new;
            } else {
                costheta = 1.;
                sintheta = 0.;
            }
            const amrex::Real vy = (-uxp[ip]*sintheta + uyp[ip]*costheta)*gaminv;
#elif (defined WARPX_DIM_XZ)
            const amrex::Real vy = uyp[ip]*gaminv;
#endif

            // Shape factor arrays
            // Note that there are extra values above and below
            // to possibly hold the factor for the old particle
            // which can be at a different grid location.
            amrex::Real sx_new[depos_order + 3] = {0.};
            amrex::Real sx_old[depos_order + 3] = {0.};
#if (defined WARPX_DIM_3D)
            amrex::Real sy_new[depos_order + 3] = {0.};
            amrex::Real sy_old[depos_order + 3] = {0.};
#endif
            amrex::Real sz_new[depos_order + 3] = {0.};
            amrex::Real sz_old[depos_order + 3] = {0.};

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
                    amrex::Gpu::Atomic::Add( &Jx_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0), sdxi);
                }
            }
            for (int k=dkl; k<=depos_order+2-dku; k++) {
                for (int i=dil; i<=depos_order+2-diu; i++) {
                    const amrex::Real sdyj = wq*vy*invvol*((sz_new[k] + 0.5*(sz_old[k] - sz_new[k]))*sx_new[i] +
                                                           (0.5*sz_new[k] + 1./3.*(sz_old[k] - sz_new[k]))*(sx_old[i] - sx_new[i]));
                    amrex::Gpu::Atomic::Add( &Jy_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0), sdyj);
                }
            }
            for (int i=dil; i<=depos_order+2-diu; i++) {
                amrex::Real sdzk = 0.;
                for (int k=dkl; k<=depos_order+1-dku; k++) {
                    sdzk += wqz*(sz_old[k] - sz_new[k])*(sx_new[i] + 0.5*(sx_old[i] - sx_new[i]));
                    amrex::Gpu::Atomic::Add( &Jz_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0), sdzk);
                }
            }


#endif
        }
        );
}

#endif // CURRENTDEPOSITION_H_
