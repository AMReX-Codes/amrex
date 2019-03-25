#ifndef _compute_flux_3d_H_
#define _compute_flux_3d_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AmrCoreAdv_F.H>

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_INLINE
void flux_x(int i, int j, int k,
            Array4<Real> const& phi,
            Array4<Real> const& vx,
            Array4<Real> const& px,
            Array4<Real> const& slope,
            const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    px(i,j,k) = ( (vx(i,j,k) < 0) ? 
                phi(i  ,j,k) - slope(i  ,j,k)*(0.5 + 0.5*dtdx[0]*vx(i,j,k)) : 
                phi(i-1,j,k) + slope(i-1,j,k)*(0.5 - 0.5*dtdx[0]*vx(i,j,k)) );
}

AMREX_GPU_DEVICE
AMREX_INLINE
void flux_y(int i, int j, int k,
            Array4<Real> const& phi,
            Array4<Real> const& vy,
            Array4<Real> const& py,
            Array4<Real> const& slope,
            const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    py(i,j,k) = ( (vy(i,j,k) < 0) ? 
                phi(i,j  ,k) - slope(i,j  ,k)*(0.5 + 0.5*dtdx[0]*vy(i,j,k)) : 
                phi(i,j-1,k) + slope(i,j-1,k)*(0.5 - 0.5*dtdx[0]*vy(i,j,k)) );
}

AMREX_GPU_DEVICE
AMREX_INLINE
void flux_z(int i, int j, int k,
            Array4<Real> const& phi,
            Array4<Real> const& vz,
            Array4<Real> const& pz,
            Array4<Real> const& slope,
            const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    pz(i,j,k) = ( (vz(i,j,k) < 0) ? 
                phi(i,j,k  ) - slope(i,j,k  )*(0.5 + 0.5*dtdx[0]*vz(i,j,k)) : 
                phi(i,j,k-1) + slope(i,j,k-1)*(0.5 - 0.5*dtdx[0]*vz(i,j,k)) );
}

AMREX_GPU_DEVICE
AMREX_INLINE
void flux_xy(int i, int j, int k, 
             AMREX_D_DECL(Array4<Real> const& vx, 
                          Array4<Real> const& vy,
                          Array4<Real> const& vz),
             AMREX_D_DECL(Array4<Real> const& px,
                          Array4<Real> const& py,
                          Array4<Real> const& pz),
             Array4<Real> const& pxy,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    pxy(i,j,k) = ( (vx(i,j,k) < 0) ? 
                 px(i,j,k) - dtdx[1]/3.0 * ( 0.5*(vy(i,  j+1,k) + vy(i  ,j,k)) * (py(i  ,j+1,k) - py(i  ,j,k))) : 
                 px(i,j,k) - dtdx[1]/3.0 * ( 0.5*(vy(i-1,j+1,k) + vy(i-1,j,k)) * (py(i-1,j+1,k) - py(i-1,j,k))) );
}

AMREX_GPU_DEVICE
AMREX_INLINE
void flux_xz(int i, int j, int k,
             AMREX_D_DECL(Array4<Real> const& vx, 
                          Array4<Real> const& vy,
                          Array4<Real> const& vz),
             AMREX_D_DECL(Array4<Real> const& px,
                          Array4<Real> const& py,
                          Array4<Real> const& pz),
             Array4<Real> const& pxz,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    pxz(i,j,k) = ( (vx(i,j,k) < 0) ? 
                 px(i,j,k) - dtdx[2]/3.0 * ( 0.5*(vz(i,  j,k+1) + vz(i  ,j,k)) * (pz(i  ,j,k+1) - pz(i  ,j,k))) : 
                 px(i,j,k) - dtdx[2]/3.0 * ( 0.5*(vz(i-1,j,k+1) + vz(i-1,j,k)) * (pz(i-1,j,k+1) - pz(i-1,j,k))) );
}

AMREX_GPU_DEVICE
AMREX_INLINE
void flux_yx(int i, int j, int k,
             AMREX_D_DECL(Array4<Real> const& vx, 
                          Array4<Real> const& vy,
                          Array4<Real> const& vz),
             AMREX_D_DECL(Array4<Real> const& px,
                          Array4<Real> const& py,
                          Array4<Real> const& pz),
             Array4<Real> const& pyx,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    pyx(i,j,k) = ( (vy(i,j,k) < 0) ? 
                 py(i,j,k) - dtdx[0]/3.0 * ( 0.5*(vx(i+1,j  ,k) + vx(i,j  ,k)) * (px(i+1,j  ,k) - px(i,j  ,k))) : 
                 py(i,j,k) - dtdx[0]/3.0 * ( 0.5*(vx(i+1,j-1,k) + vx(i,j-1,k)) * (px(i+1,j-1,k) - px(i,j-1,k))) );
}

AMREX_GPU_DEVICE
AMREX_INLINE
void flux_yz(int i, int j, int k,
             AMREX_D_DECL(Array4<Real> const& vx, 
                          Array4<Real> const& vy,
                          Array4<Real> const& vz),
             AMREX_D_DECL(Array4<Real> const& px,
                          Array4<Real> const& py,
                          Array4<Real> const& pz),
             Array4<Real> const& pyz,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    pyz(i,j,k) = ( (vy(i,j,k) < 0) ? 
                 py(i,j,k) - dtdx[2]/3.0 * ( 0.5*(vz(i,  j,k+1) + vz(i,j  ,k)) * (pz(i,j  ,k+1) - pz(i,j  ,k))) : 
                 py(i,j,k) - dtdx[2]/3.0 * ( 0.5*(vz(i,j-1,k+1) + vz(i,j-1,k)) * (pz(i,j-1,k+1) - pz(i,j-1,k))) );
}

AMREX_GPU_DEVICE
AMREX_INLINE
void flux_zx(int i, int j, int k,
             AMREX_D_DECL(Array4<Real> const& vx, 
                          Array4<Real> const& vy,
                          Array4<Real> const& vz),
             AMREX_D_DECL(Array4<Real> const& px,
                          Array4<Real> const& py,
                          Array4<Real> const& pz),
             Array4<Real> const& pzx,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    pzx(i,j,k) = ( (vz(i,j,k) < 0) ? 
                 pz(i,j,k) - dtdx[0]/3.0 * ( 0.5*(vx(i+1,j,k  ) + vx(i,j,k  )) * (px(i+1,j,k  ) - px(i,j,k  ))) : 
                 pz(i,j,k) - dtdx[0]/3.0 * ( 0.5*(vx(i+1,j,k-1) + vx(i,j,k-1)) * (px(i+1,j,k-1) - px(i,j,k-1))) );
}

AMREX_GPU_DEVICE
AMREX_INLINE
void flux_zy(int i, int j, int k,
             AMREX_D_DECL(Array4<Real> const& vx, 
                          Array4<Real> const& vy,
                          Array4<Real> const& vz),
             AMREX_D_DECL(Array4<Real> const& px,
                          Array4<Real> const& py,
                          Array4<Real> const& pz),
             Array4<Real> const& pzy,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    pzy(i,j,k) = ( (vz(i,j,k) < 0) ? 
                 pz(i,j,k) - dtdx[1]/3.0 * ( 0.5*(vy(i,j+1,k  ) + vy(i,j,k  )) * (py(i,j+1,k  ) - py(i,j,k  ))) : 
                 pz(i,j,k) - dtdx[1]/3.0 * ( 0.5*(vy(i,j+1,k-1) + vy(i,j,k-1)) * (py(i,j+1,k-1) - py(i,j,k-1))) );
}

AMREX_GPU_DEVICE
AMREX_INLINE
void combine_flux_x(int i, int j, int k,
                    AMREX_D_DECL(Array4<Real> const& vx, 
                                 Array4<Real> const& vy,
                                 Array4<Real> const& vz),
                    Array4<Real> const& px,
                    Array4<Real> const& pyz,
                    Array4<Real> const& pzy,
                    Array4<Real> const& fx,
                    const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    px(i,j,k) = ( (vx(i,j,k) < 0) ? 
                px(i,j,k) - 0.5*dtdx[1] * ( 0.5*(vy(i  ,j+1,k  ) + vy(i  ,j,k)) * (pyz(i  ,j+1,k  )-pyz(i  ,j,k)))
                          - 0.5*dtdx[2] * ( 0.5*(vz(i  ,j  ,k+1) + vz(i  ,j,k)) * (pzy(i  ,j  ,k+1)-pzy(i  ,j,k))) :
                px(i,j,k) - 0.5*dtdx[1] * ( 0.5*(vy(i-1,j+1,k  ) + vy(i-1,j,k)) * (pyz(i-1,j+1,k  )-pyz(i-1,j,k)))
                          - 0.5*dtdx[2] * ( 0.5*(vz(i-1,j  ,k+1) + vz(i-1,j,k)) * (pzy(i-1,j  ,k+1)-pzy(i-1,j,k))) );

    fx(i,j,k) = vx(i,j,k)*px(i,j,k);
}

AMREX_GPU_DEVICE
AMREX_INLINE
void combine_flux_y(int i, int j, int k,
                    AMREX_D_DECL(Array4<Real> const& vx, 
                                 Array4<Real> const& vy,
                                 Array4<Real> const& vz),
                    Array4<Real> const& py,
                    Array4<Real> const& pxz,
                    Array4<Real> const& pzx,
                    Array4<Real> const& fy,
                    const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    py(i,j,k) = ( (vy(i,j,k) < 0) ? 
                py(i,j,k) - 0.5*dtdx[0] * ( 0.5*(vx(i+1,j  ,k  ) + vx(i,j  ,k)) * (pxz(i+1,j  ,k  )-pxz(i,j  ,k)))
                          - 0.5*dtdx[2] * ( 0.5*(vz(i,  j  ,k+1) + vz(i,j  ,k)) * (pzx(i,  j  ,k+1)-pzx(i,j  ,k))) :
                py(i,j,k) - 0.5*dtdx[0] * ( 0.5*(vx(i+1,j-1,k  ) + vx(i,j-1,k)) * (pxz(i+1,j-1,k  )-pxz(i,j-1,k)))
                          - 0.5*dtdx[2] * ( 0.5*(vz(i  ,j-1,k+1) + vz(i,j-1,k)) * (pzx(i  ,j-1,k+1)-pzx(i,j-1,k))) );

    fy(i,j,k) = vy(i,j,k)*py(i,j,k);
}

AMREX_GPU_DEVICE
AMREX_INLINE
void combine_flux_z(int i, int j, int k,
                    AMREX_D_DECL(Array4<Real> const& vx, 
                                 Array4<Real> const& vy,
                                 Array4<Real> const& vz),
                    Array4<Real> const& pz,
                    Array4<Real> const& pxy,
                    Array4<Real> const& pyx,
                    Array4<Real> const& fz,
                    const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    pz(i,j,k) = ( (vz(i,j,k) < 0) ? 
                pz(i,j,k) - 0.5*dtdx[0] * ( 0.5*(vx(i+1,j  ,k  ) + vx(i,j,k  )) * (pxy(i+1,j  ,k  )-pxy(i,j,k  )))
                          - 0.5*dtdx[1] * ( 0.5*(vy(i,  j+1,k  ) + vy(i,j,k  )) * (pyx(i,  j+1,k  )-pyx(i,j,k  ))) :
                pz(i,j,k) - 0.5*dtdx[0] * ( 0.5*(vx(i+1,j  ,k-1) + vx(i,j,k-1)) * (pxy(i+1,j  ,k-1)-pxy(i,j,k-1)))
                          - 0.5*dtdx[1] * ( 0.5*(vy(i  ,j+1,k-1) + vy(i,j,k-1)) * (pyx(i  ,j+1,k-1)-pyx(i,j,k-1))) );

    fz(i,j,k) = vz(i,j,k)*pz(i,j,k);
}

#endif
