#ifndef CURRENTDEPOSITION_H_
#define CURRENTDEPOSITION_H_

using namespace amrex;

// Compute shape factor and return index of leftmost cell where
// particle writes.
// Specialized templates are defined below for orders 1, 2 and 3.
template <int depos_order>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor(Real* const sx, Real xint);

// Compute shape factor for order 1.
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <1> (Real* const sx, Real xmid){
    int j = (int) xmid;
    Real xint = xmid-j;
    sx[0] = 1.0 - xint;
    sx[1] = xint;
    return j;
}

// Compute shape factor for order 2.
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <2> (Real* const sx, Real xmid){
    int j = (int) (xmid+0.5);
    Real xint = xmid-j;
    sx[0] = 0.5*(0.5-xint)*(0.5-xint);
    sx[1] = 0.75-xint*xint;
    sx[2] = 0.5*(0.5+xint)*(0.5+xint);
    // index of the leftmost cell where particle deposits
    return j-1;
}

// Compute shape factor for order 3.
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <3> (Real* const sx, Real xmid){
    int j = (int) xmid;
    Real xint = xmid-j;
    sx[0] = 1.0/6.0*(1.0-xint)*(1.0-xint)*(1.0-xint);
    sx[1] = 2.0/3.0-xint*xint*(1-xint/2.0);
    sx[2] = 2.0/3.0-(1-xint)*(1-xint)*(1.0-0.5*(1-xint));
    sx[3] = 1.0/6.0*xint*xint*xint;
    // index of the leftmost cell where particle deposits
    return j-1;
}

/* \brief Current Deposition for thread thread_num
 * /param xp, yp, zp   : Pointer to arrays of particle positions.
 * \param wp           : Pointer to array of particle weights.
 * \param uxp uyp uzp  : Pointer to arrays of particle momentum.
 * \param jx_arr       : Array4 of current density, either full array or tile.
 * \param jy_arr       : Array4 of current density, either full array or tile.
 * \param jz_arr       : Array4 of current density, either full array or tile.
 * \param offset       : Index of first particle for which current is deposited
 * \param np_to_depose : Number of particles for which current is deposited.
                        Particles [offset,offset+np_tp_depose] deposit current.
 * \param dt           : Time step for particle level
 * \param dx           : 3D cell size
 * \param xyzmin       : Physical lower bounds of domain.
 * \param lo           : Index lower bounds of domain.
 * \param stagger_shift: 0 if nodal, 0.5 if staggered.
 * /param q            : species charge.
 */
template <int depos_order>
void doDepositionShapeN(const Real * const xp, const Real * const yp, const Real * const zp,
                        const Real * const wp, const Real * const uxp,
                        const Real * const uyp, const Real * const uzp,
                        const amrex::Array4<amrex::Real>& jx_arr, 
                        const amrex::Array4<amrex::Real>& jy_arr, 
                        const amrex::Array4<amrex::Real>& jz_arr,
                        const long offset, const long np_to_depose, 
                        const amrex::Real dt, const std::array<amrex::Real,3>& dx,
                        const std::array<Real, 3> xyzmin,
                        const Dim3 lo,
                        const amrex::Real stagger_shift, 
                        const amrex::Real q)
{
    const Real dxi = 1.0/dx[0];
    const Real dzi = 1.0/dx[2];
    const Real dts2dx = 0.5*dt*dxi;
    const Real dts2dz = 0.5*dt*dzi;
#if (AMREX_SPACEDIM == 2)
    const Real invvol = dxi*dzi;
#else // (AMREX_SPACEDIM == 3)
    const Real dyi = 1.0/dx[1];
    const Real dts2dy = 0.5*dt*dyi;
    const Real invvol = dxi*dyi*dzi;
#endif

    const Real xmin = xyzmin[0];
    const Real ymin = xyzmin[1];
    const Real zmin = xyzmin[2];
    const Real clightsq = 1.0/PhysConst::c/PhysConst::c;

    // Loop over particles and deposit into jx_arr, jy_arr and jz_arr
    ParallelFor( np_to_depose,
                 [=] AMREX_GPU_DEVICE (long ip) {
                     // --- Get particle quantities
                     const Real gaminv = 1.0/std::sqrt(1.0 + uxp[ip]*uxp[ip]*clightsq
						       + uyp[ip]*uyp[ip]*clightsq
						       + uzp[ip]*uzp[ip]*clightsq);
                     const Real wq  = q*wp[ip];
                     const Real vx  = uxp[ip]*gaminv;
                     const Real vy  = uyp[ip]*gaminv;
                     const Real vz  = uzp[ip]*gaminv;
                     // wqx, wqy wqz are particle current in each direction 
                     const Real wqx = wq*invvol*vx;
                     const Real wqy = wq*invvol*vy;
                     const Real wqz = wq*invvol*vz;

                     // --- Compute shape factors
                     // x direction
                     // Get particle position after 1/2 push back in position
                     const Real xmid = (xp[ip]-xmin)*dxi-dts2dx*vx;
                     // Compute shape factors for node-centered quantities
                     Real AMREX_RESTRICT sx [depos_order + 1];
                     // j: leftmost grid point (node-centered) that the particle touches
                     const int j  = compute_shape_factor<depos_order>(sx,  xmid);
                     // Compute shape factors for cell-centered quantities
                     Real AMREX_RESTRICT sx0[depos_order + 1];
                     // j0: leftmost grid point (cell-centered) that the particle touches
                     const int j0 = compute_shape_factor<depos_order>(sx0, xmid-stagger_shift);
                     
#if (AMREX_SPACEDIM == 3)
                     // y direction
                     const Real ymid= (yp[ip]-ymin)*dyi-dts2dy*vy;
                     Real AMREX_RESTRICT sy [depos_order + 1];
		     const int k  = compute_shape_factor<depos_order>(sy,  ymid);
                     Real AMREX_RESTRICT sy0[depos_order + 1];
		     const int k0 = compute_shape_factor<depos_order>(sy0, ymid-stagger_shift);
#endif
                     // z direction
                     const Real zmid= (zp[ip]-zmin)*dzi-dts2dz*vz;
                     Real AMREX_RESTRICT sz [depos_order + 1];
		     const int l  = compute_shape_factor<depos_order>(sz,  zmid);
                     Real AMREX_RESTRICT sz0[depos_order + 1];
		     const int l0 = compute_shape_factor<depos_order>(sz0, zmid-stagger_shift);

                     // Deposit current into jx_arr, jy_arr and jz_arr
#if (AMREX_SPACEDIM == 2)
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
#else // (AMREX_SPACEDIM == 3)
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

// Compute shape factor and return index of leftmost cell where
// particle writes.
// Specialized templates are defined below for orders 1, 2 and 3.
template <int depos_order>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor (Real* const sx, const Real x_old, const int i_new);

// Compute shape factor for order 1.
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor <1> (Real* const sx, const Real x_old, const int i_new){
    const int i = (int) x_old;
    const int i_shift = i - i_new;
    const Real xint = x_old - i;
    sx[1+i_shift] = 1.0 - xint;
    sx[2+i_shift] = xint;
    return i;
}

// Compute shape factor for order 2.
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor <2> (Real* const sx, const Real x_old, const int i_new){
    const int i = (int) (x_old+0.5);
    const int i_shift = i - (i_new + 1);
    const Real xint = x_old - i;
    sx[1+i_shift] = 0.5*(0.5-xint)*(0.5-xint);
    sx[2+i_shift] = 0.75-xint*xint;
    sx[3+i_shift] = 0.5*(0.5+xint)*(0.5+xint);
    // index of the leftmost cell where particle deposits
    return i-1;
}

// Compute shape factor for order 3.
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shifted_shape_factor <3> (Real* const sx, const Real x_old, const int i_new){
    const int i = (int) x_old;
    const int i_shift = i - (i_new + 1);
    const Real xint = x_old - i;
    sx[1+i_shift] = 1.0/6.0*(1.0-xint)*(1.0-xint)*(1.0-xint);
    sx[2+i_shift] = 2.0/3.0-xint*xint*(1-xint/2.0);
    sx[3+i_shift] = 2.0/3.0-(1-xint)*(1-xint)*(1.0-0.5*(1-xint));
    sx[4+i_shift] = 1.0/6.0*xint*xint*xint;
    // index of the leftmost cell where particle deposits
    return i-1;
}

/* \brief Esirkepov Current Deposition for thread thread_num
 * /param xp, yp, zp   : Pointer to arrays of particle positions.
 * \param wp           : Pointer to array of particle weights.
 * \param uxp uyp uzp  : Pointer to arrays of particle momentum.
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
void doEsirkepovDepositionShapeN (const Real * const xp, const Real * const yp, const Real * const zp,
                                  const Real * const wp, const Real * const uxp,
                                  const Real * const uyp, const Real * const uzp,
                                  const amrex::Array4<amrex::Real>& Jx_arr,
                                  const amrex::Array4<amrex::Real>& Jy_arr,
                                  const amrex::Array4<amrex::Real>& Jz_arr,
                                  const long np_to_depose,
                                  const amrex::Real dt, const std::array<amrex::Real,3>& dx,
                                  const std::array<Real, 3> xyzmin,
                                  const Dim3 lo,
                                  const amrex::Real q)
{
    const Real dxi = 1.0/dx[0];
    const Real dtsdx0 = dt*dxi;
    const Real xmin = xyzmin[0];
#if (AMREX_SPACEDIM == 3)
    const Real dyi = 1.0/dx[1];
    const Real dtsdy0 = dt*dyi;
    const Real ymin = xyzmin[1];
#endif
    const Real dzi = 1.0/dx[2];
    const Real dtsdz0 = dt*dzi;
    const Real zmin = xyzmin[2];

#if (AMREX_SPACEDIM == 3)
    const Real invdtdx = 1.0/(dt*dx[1]*dx[2]);
    const Real invdtdy = 1.0/(dt*dx[0]*dx[2]);
    const Real invdtdz = 1.0/(dt*dx[0]*dx[1]);
#elif (AMREX_SPACEDIM == 2)
    const Real invdtdx = 1.0/(dt*dx[2]);
    const Real invdtdz = 1.0/(dt*dx[0]);
    const Real invvol = 1.0/(dx[0]*dx[2]);
#endif

    const Real clightsq = 1.0/PhysConst::c/PhysConst::c;

    // Loop over particles and deposit into Jx_arr, Jy_arr and Jz_arr
    ParallelFor( np_to_depose,
                 [=] AMREX_GPU_DEVICE (long ip) {

                     // --- Get particle quantities
                     const Real gaminv = 1.0/std::sqrt(1.0 + uxp[ip]*uxp[ip]*clightsq
                                                           + uyp[ip]*uyp[ip]*clightsq
                                                           + uzp[ip]*uzp[ip]*clightsq);

                     // wqx, wqy wqz are particle current in each direction
                     const Real wq = q*wp[ip];
                     const Real wqx = wq*invdtdx;
#if (AMREX_SPACEDIM == 3)
                     const Real wqy = wq*invdtdy;
#endif
                     const Real wqz = wq*invdtdz;

                     // computes current and old position in grid units
                     const Real x_new = (xp[ip] - xmin)*dxi;
                     const Real x_old = x_new - dtsdx0*uxp[ip]*gaminv;
#if (AMREX_SPACEDIM == 3)        
                     const Real y_new = (yp[ip] - ymin)*dyi;
                     const Real y_old = y_new - dtsdy0*uyp[ip]*gaminv;
#endif
                     const Real z_new = (zp[ip] - zmin)*dzi;
                     const Real z_old = z_new - dtsdz0*uzp[ip]*gaminv;

                     // Shape factor arrays
                     // Note that there are extra values above and below
                     // to possibly hold the factor for the old particle
                     // which can be at a different grid location.
                     Real AMREX_RESTRICT sx_new[depos_order + 3] = {0.};
                     Real AMREX_RESTRICT sx_old[depos_order + 3] = {0.};
#if (AMREX_SPACEDIM == 3)
                     Real AMREX_RESTRICT sy_new[depos_order + 3] = {0.};
                     Real AMREX_RESTRICT sy_old[depos_order + 3] = {0.};
#endif
                     Real AMREX_RESTRICT sz_new[depos_order + 3] = {0.};
                     Real AMREX_RESTRICT sz_old[depos_order + 3] = {0.};

                     // --- Compute shape factors
                     // Compute shape factors for position as they are now and at old positions
                     // [ijk]_new: leftmost grid point that the particle touches
                     const int i_new = compute_shape_factor<depos_order>(sx_new+1, x_new);
                     const int i_old = compute_shifted_shape_factor<depos_order>(sx_old, x_old, i_new);
#if (AMREX_SPACEDIM == 3)
                     const int j_new = compute_shape_factor<depos_order>(sy_new+1, y_new);
                     const int j_old = compute_shifted_shape_factor<depos_order>(sy_old, y_old, j_new);
#endif 
                     const int k_new = compute_shape_factor<depos_order>(sz_new+1, z_new);
                     const int k_old = compute_shifted_shape_factor<depos_order>(sz_old, z_old, k_new);

                     // computes min/max positions of current contributions
                     int dil = 1, diu = 1;
                     if (i_old < i_new) dil = 0;
                     if (i_old > i_new) diu = 0;
#if (AMREX_SPACEDIM == 3)
                     int djl = 1, dju = 1;
                     if (j_old < j_new) djl = 0;
                     if (j_old > j_new) dju = 0;
#endif
                     int dkl = 1, dku = 1;
                     if (k_old < k_new) dkl = 0;
                     if (k_old > k_new) dku = 0;

#if (AMREX_SPACEDIM == 3)

                     for (int k=dkl; k<=depos_order+2-dku; k++) {
                         for (int j=djl; j<=depos_order+2-dju; j++) {
                             Real sdxi = 0.;
                             for (int i=dil; i<=depos_order+1-diu; i++) {
                                 sdxi += wqx*(sx_old[i] - sx_new[i])*((sy_new[j] + 0.5*(sy_old[j] - sy_new[j]))*sz_new[k] +
                                         (0.5*sy_new[j] + 1./3.*(sy_old[j] - sy_new[j]))*(sz_old[k] - sz_new[k]));
                                 amrex::Gpu::Atomic::Add( &Jx_arr(lo.x+i_new-1+i, lo.y+j_new-1+j, lo.z+k_new-1+k), sdxi);
                             }
                         }
                     }
                     for (int k=dkl; k<=depos_order+2-dku; k++) {
                         for (int i=dil; i<=depos_order+2-diu; i++) {
                             Real sdyj = 0.;
                             for (int j=djl; j<=depos_order+1-dju; j++) {
                                 sdyj += wqy*(sy_old[j] - sy_new[j])*((sz_new[k] + 0.5*(sz_old[k] - sz_new[k]))*sx_new[i] +
                                         (0.5*sz_new[k] + 1./3.*(sz_old[k] - sz_new[k]))*(sx_old[i] - sx_new[i]));
                                 amrex::Gpu::Atomic::Add( &Jy_arr(lo.x+i_new-1+i, lo.y+j_new-1+j, lo.z+k_new-1+k), sdyj);
                             }
                         }
                     }
                     for (int j=djl; j<=depos_order+2-dju; j++) {
                         for (int i=dil; i<=depos_order+2-diu; i++) {
                             Real sdzk = 0.;
                             for (int k=dkl; k<=depos_order+1-dku; k++) {
                                 sdzk += wqz*(sz_old[k] - sz_new[k])*((sx_new[i] + 0.5*(sx_old[i] - sx_new[i]))*sy_new[j] +
                                         (0.5*sx_new[i] + 1./3.*(sx_old[i] - sx_new[i]))*(sy_old[j] - sy_new[j]));
                                 amrex::Gpu::Atomic::Add( &Jz_arr(lo.x+i_new-1+i, lo.y+j_new-1+j, lo.z+k_new-1+k), sdzk);
                             }
                         }
                     }

#elif (AMREX_SPACEDIM == 2)

                    for (int k=dkl; k<=depos_order+2-dku; k++) {
                        Real sdxi = 0.;
                        for (int i=dil; i<=depos_order+1-diu; i++) {
                            sdxi += wqx*(sx_old[i] - sx_new[i])*(sz_new[k] + 0.5*(sz_old[k] - sz_new[k]));
                            amrex::Gpu::Atomic::Add( &Jx_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0), sdxi);
                        }
                    }
                    for (int k=dkl; k<=depos_order+2-dku; k++) {
                        for (int i=dil; i<=depos_order+2-diu; i++) {
                            const Real sdyj = wq*uyp[ip]*gaminv*invvol*((sz_new[k] + 0.5*(sz_old[k] - sz_new[k]))*sx_new[i] +
                                        (0.5*sz_new[k] + 1./3.*(sz_old[k] - sz_new[k]))*(sx_old[i] - sx_new[i]));
                            amrex::Gpu::Atomic::Add( &Jy_arr(lo.x+i_new-1+i, lo.y+k_new-1+k, 0), sdyj);
                        }
                    }
                    for (int i=dil; i<=depos_order+2-diu; i++) {
                        Real sdzk = 0.;
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
