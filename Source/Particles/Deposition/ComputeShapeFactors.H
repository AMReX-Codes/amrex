#ifndef COMPUTESHAPEFACTORS_H_
#define COMPUTESHAPEFACTORS_H_

//#include <memory>
//#include <AMReX_Particles.H>
//#include <AMReX_AmrCore.H>
//#include <MultiParticleContainer.H>
//#include <WarpXParticleContainer.H>
//#include <AMReX_AmrParGDB.H>
//#include <WarpX_f.H>
//#include <WarpX.H>

//#include <AMReX_FArrayBox.H>
//#include <WarpXParticleContainer.H>
//#include <AMReX_Particles.H>
//#include <WarpX.H> 
//#include <AMReX_AmrCore.H>

//#include <AMReX_Particles.H>
//#include <AMReX_AmrCore.H>
//#include <MultiParticleContainer.H>
//#include <WarpXParticleContainer.H>
//#include <AMReX_AmrParGDB.H>
// #include <WarpX_f.H>
// #include <WarpX.H>
//#include "WarpXParticleContainer.H"
//#include <AMReX_FArrayBox.H>
//#include <AMReX_Dim3.H>

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
void compute_shape_factor(Real* sx, Real xint) {};

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_shape_factor <1> (Real* sx, Real xint){
    sx[0] = 1.0 - xint;
    sx[1] = xint;
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_shape_factor <2> (Real* sx, Real xint){
    Print()<<"compute_shape_factor <2>\n";
    sx[0] = 0.5*(0.5-xint)*(0.5-xint);
    sx[1] = 0.75-xint*xint;
    sx[2] = 0.5*(0.5+xint)*(0.5+xint);
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_shape_factor <3> (Real* sx, Real xint){
    Print()<<"compute_shape_factor <2>\n";
    sx[0] = 1.0/6.0*(1.0-xint)*(1.0-xint)*(1.0-xint);
    sx[1] = 2.0/3.0-xint*xint*(1-xint/2.0);
    sx[2] = 2.0/3.0-(1-xint)*(1-xint)*(1.0-0.5*(1-xint));
    sx[3] = 1.0/6.0*xint*xint*xint;
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

    Print()<<"toto\n";

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
                     // particle position in cell unit (1/2 timestep back)
                     // Real x = (Real) (xp[ip]-xmin)*dxi;
                     // particle position in cell unit
                     // Real xmid=x-dts2dx*vx;
                     Real xmid= (Real) (xp[ip]-xmin)*dxi-dts2dx*vx;
                     // x index of cell containing particle (cell-centered)
                     int j = (int) xmid;
                     // x index of cell containing particle (node-centered)
                     int j0= (int) (xmid-stagger_shift);
                     // for the cell-centered quantities, compute current on 
                     // each point
                     Real xint = xmid-j;
                     Real sx[depos_order + 1], sx0[depos_order + 1];
                     // Real sx[2] = {1.0 - xint, xint};
                     compute_shape_factor<depos_order>(sx, xint);
                     // for the node-centered quantities, compute current on 
                     // each point
                     xint = xmid-stagger_shift-j0;
                     compute_shape_factor<depos_order>(sx0, xint);
                     // Real sx0[2] = {1.0 - xint, xint};
                     // compute_shape_factor<WarpX::nox>(sx0, xint);
                     
#if (AMREX_SPACEDIM == 3)
                     // --- y direction
                     Real y = (Real) (yp[ip]-ymin)*dyi;
                     Real ymid=y-dts2dy*vy;
                     int k = (int) ymid;
                     int k0= (int) (ymid-stagger_shift);
                     Real yint = ymid-k;
                     Real sy[depos_order + 1], sy0[depos_order + 1];
                     compute_shape_factor<depos_order>(sy, yint);
                     // Real sy[2] = {1.0 - yint, yint};
                     yint = ymid-stagger_shift-k0;
                     compute_shape_factor<depos_order>(sy0, yint);
                     //Real sy0[2] = {1.0 - yint, yint};
                     // compute_shape_factor<WarpX::noy>(shape_lo, shape_hi, sy, yint);
#endif
                     // --- z direction
                     Real z = (Real) (zp[ip]-zmin)*dzi;
                     Real zmid=z-dts2dz*vz;
                     int l = (int) zmid;
                     int l0= (int) (zmid-stagger_shift);
                     Real zint = zmid-l;
                     Real sz[depos_order + 1], sz0[depos_order + 1];
                     compute_shape_factor<depos_order>(sz, zint);
                     // Real sz[2] = {1.0 - zint, zint};
                     zint = zmid-stagger_shift-l0;
                     compute_shape_factor<depos_order>(sz0, zint);
                     // Real sz0[2] = {1.0 - zint, zint};
                     // compute_shape_factor<WarpX::noz>(shape_lo, shape_hi, sz, zint);
/*
  #if (AMREX_SPACEDIM == 2)
  for (int iz=shape_lo.y; iz<=shape_hi.y; iz++){
  for (int ix=shape_lo.x; ix<=shape_hi.x; ix++){
  jx_arr(lo.x+j0+ix  ,lo.y+l +iz  ,0) += sx0[ix]*sz [iz]*wqx;
  jy_arr(lo.x+j +ix  ,lo.y+l +iz  ,0) += sx [ix]*sz [iz]*wqy;
  jz_arr(lo.x+j +ix  ,lo.y+l0+iz  ,0) += sx [ix]*sz0[iz]*wqz;
  }
  }
  #else // (AMREX_SPACEDIM == 3)
  for (int iz=shape_lo.z; iz<=shape_hi.z; iz++){
  for (int iy=shape_lo.y; iy<=shape_hi.y; iy++){
  for (int ix=shape_lo.x; ix<=shape_hi.x; ix++){
  jx_arr(lo.x+j0+ix  ,lo.y+k +iy  ,lo.z+l +iz  ) += sx0[ix]*sy [iy]*sz [iz]*wqx;
  jy_arr(lo.x+j +ix  ,lo.y+k0+iy  ,lo.z+l +iz  ) += sx [ix]*sy0[iy]*sz [iz]*wqy;
  jz_arr(lo.x+j +ix  ,lo.y+k +iy  ,lo.z+l0+iz  ) += sx [ix]*sy [iy]*sz0[iz]*wqz;
  }
  }
  }
  #endif
*/
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
                                 jx_arr(lo.x+j0+ix+shape_lo.x, lo.y+k +iy+shape_lo.y, lo.z+l +iz+shape_lo.z) += sx0[ix]*sy [iy]*sz [iz]*wqx;
                                 jy_arr(lo.x+j +ix+shape_lo.x, lo.y+k0+iy+shape_lo.y, lo.z+l +iz+shape_lo.z) += sx [ix]*sy0[iy]*sz [iz]*wqy;
                                 jz_arr(lo.x+j +ix+shape_lo.x, lo.y+k +iy+shape_lo.y, lo.z+l0+iz+shape_lo.z) += sx [ix]*sy [iy]*sz0[iz]*wqz;
                             }
                         }
                     }
#endif
                 }
        );
}

#endif // COMPUTESHAPEFACTORS_H_
