#ifndef CURRENTDEPOSITION_H_
#define CURRENTDEPOSITION_H_

using namespace amrex;


template <int depos_order>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_shape_factor_lo(amrex::Dim3& shape_lo) {};

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_shape_factor_lo <1> (amrex::Dim3& shape_lo) {
    shape_lo = {0, 0, 0};
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_shape_factor_lo <2> (amrex::Dim3& shape_lo) {
    shape_lo = {-1, -1, -1};
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_shape_factor_lo <3> (amrex::Dim3& shape_lo) {
    shape_lo = {-1, -1, -1};
}

template <int depos_order>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor(Real* sx, Real xint) {return 0;};

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <1> (Real* sx, Real xmid){
    int j = (int) xmid;
    Real xint = xmid-j;
    sx[0] = 1.0 - xint;
    sx[1] = xint;
    return j;
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <2> (Real* sx, Real xmid){
    int j = (int) (xmid+0.5);
    Real xint = xmid-j;
    sx[0] = 0.5*(0.5-xint)*(0.5-xint);
    sx[1] = 0.75-xint*xint;
    sx[2] = 0.5*(0.5+xint)*(0.5+xint);
    return j;
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_shape_factor <3> (Real* sx, Real xmid){
    int j = (int) xmid;
    Real xint = xmid-j;
    sx[0] = 1.0/6.0*(1.0-xint)*(1.0-xint)*(1.0-xint);
    sx[1] = 2.0/3.0-xint*xint*(1-xint/2.0);
    sx[2] = 2.0/3.0-(1-xint)*(1-xint)*(1.0-0.5*(1-xint));
    sx[3] = 1.0/6.0*xint*xint*xint;
    return j;
}

template <int depos_order>
void doDepositionShapeN(Real* xp, Real* yp, Real* zp,
                        Real* wp, Real* uxp,
                        Real* uyp, Real* uzp,
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

    Real dxi = 1.0/dx[0];
    Real dzi = 1.0/dx[2];
    Real dts2dx = 0.5*dt*dxi;
    Real dts2dz = 0.5*dt*dzi;
#if (AMREX_SPACEDIM == 2)
    Real invvol = dxi*dzi;
#else // (AMREX_SPACEDIM == 3)
    Real dyi = 1.0/dx[1];
    Real dts2dy = 0.5*dt*dyi;
    Real invvol = dxi*dyi*dzi;
#endif

    const Real xmin = xyzmin[0];
    const Real ymin = xyzmin[1];
    const Real zmin = xyzmin[2];

    Dim3 shape_lo;
    compute_shape_factor_lo<depos_order>(shape_lo);
    Real clightsq = 1.0/PhysConst::c/PhysConst::c;

    ParallelFor( np_to_depose,
                 [=] AMREX_GPU_DEVICE (long ip) {

                     // macro-particle current in all 3D (even for 2D runs)
                     // wqx, wqy and wqz are the particle current in each 
                     // direction
                     Real gaminv = 1.0/std::sqrt(1.0 + uxp[ip]*uxp[ip]*clightsq
                                                 + uyp[ip]*uyp[ip]*clightsq
                                                 + uzp[ip]*uzp[ip]*clightsq);
                     Real wq = q*wp[ip];
                     Real vx = uxp[ip]*gaminv;
                     Real wqx = wq*invvol*vx;
                     Real vy = uyp[ip]*gaminv;
                     Real wqy = wq*invvol*vy;
                     Real vz = uzp[ip]*gaminv;
                     Real wqz = wq*invvol*vz;

                     // --- x direction
                     // Get particle position after 1/2 push back in position
                     Real xmid= (Real) (xp[ip]-xmin)*dxi-dts2dx*vx;
                     // Initialize arrays for node-centered and cell-centered
                     // shape factors
                     Real sx[depos_order + 1], sx0[depos_order + 1];
                     // Compute shape factors for node-centered quantities
                     int j  = compute_shape_factor<depos_order>(sx, xmid);
                     // Compute shape factors for cell-centered quantities
                     int j0 = compute_shape_factor<depos_order>(sx0, xmid-stagger_shift);
                     
#if (AMREX_SPACEDIM == 3)
                     // --- y direction
                     Real ymid= (Real) (yp[ip]-ymin)*dyi-dts2dy*vy;
                     Real sy[depos_order + 1], sy0[depos_order + 1];
                     int k  = compute_shape_factor<depos_order>(sy, ymid);
                     int k0 = compute_shape_factor<depos_order>(sy0, ymid-stagger_shift);
#endif
                     // --- z direction
                     Real zmid=(Real) (zp[ip]-zmin)*dzi-dts2dz*vz;
                     Real sz[depos_order + 1], sz0[depos_order + 1];
                     int l  = compute_shape_factor<depos_order>(sz, zmid);
                     int l0 = compute_shape_factor<depos_order>(sz0, zmid-stagger_shift);

#if (AMREX_SPACEDIM == 2)
                     for (int iz=0; iz<=depos_order; iz++){
                         for (int ix=0; ix<=depos_order; ix++){
                             jx_arr(lo.x+j0+ix+shape_lo.x, lo.y+l +iz+shape_lo.y, 0) += sx0[ix]*sz [iz]*wqx;
                             jy_arr(lo.x+j +ix+shape_lo.x, lo.y+l +iz+shape_lo.y, 0) += sx [ix]*sz [iz]*wqy;
                             jz_arr(lo.x+j +ix+shape_lo.x, lo.y+l0+iz+shape_lo.y, 0) += sx [ix]*sz0[iz]*wqz;
                         }
                     }
#else // (AMREX_SPACEDIM == 3)
                     for (int iz=0; iz<=depos_order; iz++){
                         for (int iy=0; iy<=depos_order; iy++){
                             for (int ix=0; ix<=depos_order; ix++){
			         amrex::Gpu::Atomic::Add(&jx_arr(lo.x+j0+ix+shape_lo.x, lo.y+k +iy+shape_lo.y, lo.z+l +iz+shape_lo.z), sx0[ix]*sy [iy]*sz [iz]*wqx);
			         amrex::Gpu::Atomic::Add(&jy_arr(lo.x+j +ix+shape_lo.x, lo.y+k0+iy+shape_lo.y, lo.z+l +iz+shape_lo.z), sx [ix]*sy0[iy]*sz [iz]*wqy);
			         amrex::Gpu::Atomic::Add(&jz_arr(lo.x+j +ix+shape_lo.x, lo.y+k +iy+shape_lo.y, lo.z+l0+iz+shape_lo.z), sx [ix]*sy [iy]*sz0[iz]*wqz);
			         // amrex::Gpu::Atomic::Add(&jx(j0+joff,k +koff,l +loff), sx0[joff]*sy [koff]*sz [loff]*wqx);
                                 // jx_arr(lo.x+j0+ix+shape_lo.x, lo.y+k +iy+shape_lo.y, lo.z+l +iz+shape_lo.z) += sx0[ix]*sy [iy]*sz [iz]*wqx;
                                 // jy_arr(lo.x+j +ix+shape_lo.x, lo.y+k0+iy+shape_lo.y, lo.z+l +iz+shape_lo.z) += sx [ix]*sy0[iy]*sz [iz]*wqy;
                                 // jz_arr(lo.x+j +ix+shape_lo.x, lo.y+k +iy+shape_lo.y, lo.z+l0+iz+shape_lo.z) += sx [ix]*sy [iy]*sz0[iz]*wqz;
                             }
                         }
                     }
#endif
                 }
        );
}

#endif // CURRENTDEPOSITION_H_
