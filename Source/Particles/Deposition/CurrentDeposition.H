#ifndef CURRENTDEPOSITION_H_
#define CURRENTDEPOSITION_H_

using namespace amrex;

// Compute lower bound for current deposition shape factor:
// order 1: 0
// order 2: -1
// order 3: -1
// non-type template parameter depos_order is the order for current deposition
// shape factor. Specialized templates are defined below for orders 1, 2 and 3.
// Currently, the shape factor is the same in all directions.
template <int depos_order>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_shape_factor_lo(amrex::Dim3& shape_lo) {};

// Compute lower bound for order 1.
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_shape_factor_lo <1> (amrex::Dim3& shape_lo) {
    shape_lo = {0, 0, 0};
}

// Compute lower bound for order 2.
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_shape_factor_lo <2> (amrex::Dim3& shape_lo) {
    shape_lo = {-1, -1, -1};
}

// Compute lower bound for order 3.
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_shape_factor_lo <3> (amrex::Dim3& shape_lo) {
    shape_lo = {-1, -1, -1};
}

// Compute shape factor and return index of cell containing particle.
// Specialized templates are defined below for orders 1, 2 and 3.
template <int depos_order>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor(Real* const sx, Real xint) {return 0;};

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
    return j;
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
    return j;
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
    BL_PROFILE("PICSAR::CurrentDeposition");
    
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

    Dim3 shape_lo;
    compute_shape_factor_lo<depos_order>(shape_lo);
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
                     const Real xmid = (Real) (xp[ip]-xmin)*dxi-dts2dx*vx;
                     // Compute shape factors for node-centered quantities
                     Real AMREX_RESTRICT sx [depos_order + 1];
		     const int j  = compute_shape_factor<depos_order>(sx,  xmid);
                     // Compute shape factors for cell-centered quantities
                     Real AMREX_RESTRICT sx0[depos_order + 1];
		     const int j0 = compute_shape_factor<depos_order>(sx0, xmid-stagger_shift);
                     
#if (AMREX_SPACEDIM == 3)
                     // y direction
                     const Real ymid= (Real) (yp[ip]-ymin)*dyi-dts2dy*vy;
                     Real AMREX_RESTRICT sy [depos_order + 1];
		     const int k  = compute_shape_factor<depos_order>(sy,  ymid);
                     Real AMREX_RESTRICT sy0[depos_order + 1];
		     const int k0 = compute_shape_factor<depos_order>(sy0, ymid-stagger_shift);
#endif
                     // z direction
                     const Real zmid=(Real) (zp[ip]-zmin)*dzi-dts2dz*vz;
                     Real AMREX_RESTRICT sz [depos_order + 1];
		     const int l  = compute_shape_factor<depos_order>(sz,  zmid);
                     Real AMREX_RESTRICT sz0[depos_order + 1];
		     const int l0 = compute_shape_factor<depos_order>(sz0, zmid-stagger_shift);

                     // Deposit current into jx_arr, jy_arr and jz_arr
#if (AMREX_SPACEDIM == 2)
                     for (int iz=0; iz<=depos_order; iz++){
                         for (int ix=0; ix<=depos_order; ix++){
                             amrex::Gpu::Atomic::Add(
                                 &jx_arr(lo.x+j0+ix+shape_lo.x, 
                                         lo.y+l +iz+shape_lo.y, 
                                         0), 
                                 sx0[ix]*sz [iz]*wqx);
                             amrex::Gpu::Atomic::Add(
                                 &jy_arr(lo.x+j +ix+shape_lo.x, 
                                         lo.y+l +iz+shape_lo.y, 
                                         0), 
                                 sx [ix]*sz [iz]*wqy);
                             amrex::Gpu::Atomic::Add(
                                 &jz_arr(lo.x+j +ix+shape_lo.x, 
                                         lo.y+l0+iz+shape_lo.y, 
                                         0), 
                                 sx [ix]*sz0[iz]*wqz);
                         }
                     }
#else // (AMREX_SPACEDIM == 3)
                     for (int iz=0; iz<=depos_order; iz++){
                         for (int iy=0; iy<=depos_order; iy++){
                             for (int ix=0; ix<=depos_order; ix++){
                                 amrex::Gpu::Atomic::Add(
                                     &jx_arr(lo.x+j0+ix+shape_lo.x, 
                                             lo.y+k +iy+shape_lo.y, 
                                             lo.z+l +iz+shape_lo.z), 
                                     sx0[ix]*sy [iy]*sz [iz]*wqx);
                                 amrex::Gpu::Atomic::Add(
                                     &jy_arr(lo.x+j +ix+shape_lo.x, 
                                             lo.y+k0+iy+shape_lo.y, 
                                             lo.z+l +iz+shape_lo.z), 
                                     sx [ix]*sy0[iy]*sz [iz]*wqy);
                                 amrex::Gpu::Atomic::Add(
                                     &jz_arr(lo.x+j +ix+shape_lo.x, 
                                             lo.y+k +iy+shape_lo.y, 
                                             lo.z+l0+iz+shape_lo.z), 
                                     sx [ix]*sy [iy]*sz0[iz]*wqz);
                             }
                         }
                     }
#endif
                 }
        );
}

#endif // CURRENTDEPOSITION_H_
