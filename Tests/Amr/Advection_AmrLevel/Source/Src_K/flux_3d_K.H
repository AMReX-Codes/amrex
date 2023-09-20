#ifndef flux_3d_K_H_
#define flux_3d_K_H_

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_x(int i, int j, int k,
            amrex::Array4<amrex::Real> const& px,
            amrex::Array4<amrex::Real const> const& phi,
            amrex::Array4<amrex::Real const> const& vx,
            amrex::Array4<amrex::Real const> const& slope,
            amrex::Real hdtdx)
{
    px(i,j,k) = ( (vx(i,j,k) < 0) ?
                phi(i  ,j,k) - slope(i  ,j,k)*(0.5 + hdtdx*vx(i,j,k)) :
                phi(i-1,j,k) + slope(i-1,j,k)*(0.5 - hdtdx*vx(i,j,k)) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_y(int i, int j, int k,
            amrex::Array4<amrex::Real> const& py,
            amrex::Array4<amrex::Real const> const& phi,
            amrex::Array4<amrex::Real const> const& vy,
            amrex::Array4<amrex::Real const> const& slope,
            amrex::Real hdtdy)
{
    py(i,j,k) = ( (vy(i,j,k) < 0) ?
                phi(i,j  ,k) - slope(i,j  ,k)*(0.5 + hdtdy*vy(i,j,k)) :
                phi(i,j-1,k) + slope(i,j-1,k)*(0.5 - hdtdy*vy(i,j,k)) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_z(int i, int j, int k,
            amrex::Array4<amrex::Real> const& pz,
            amrex::Array4<amrex::Real const> const& phi,
            amrex::Array4<amrex::Real const> const& vz,
            amrex::Array4<amrex::Real const> const& slope,
            amrex::Real hdtdz)
{
    pz(i,j,k) = ( (vz(i,j,k) < 0) ?
                phi(i,j,k  ) - slope(i,j,k  )*(0.5 + hdtdz*vz(i,j,k)) :
                phi(i,j,k-1) + slope(i,j,k-1)*(0.5 - hdtdz*vz(i,j,k)) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_xy(int i, int j, int k,
             amrex::Array4<amrex::Real> const& pxy,
             amrex::Array4<amrex::Real const> const& vx,
             amrex::Array4<amrex::Real const> const& vy,
             amrex::Array4<amrex::Real const> const& px,
             amrex::Array4<amrex::Real const> const& py,
             amrex::Real tdtdy)
{
    pxy(i,j,k) = ( (vx(i,j,k) < 0) ?
                 px(i,j,k) - tdtdy * ( 0.5*(vy(i,  j+1,k) + vy(i  ,j,k)) * (py(i  ,j+1,k) - py(i  ,j,k))) :
                 px(i,j,k) - tdtdy * ( 0.5*(vy(i-1,j+1,k) + vy(i-1,j,k)) * (py(i-1,j+1,k) - py(i-1,j,k))) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_xz(int i, int j, int k,
             amrex::Array4<amrex::Real> const& pxz,
             amrex::Array4<amrex::Real const> const& vx,
             amrex::Array4<amrex::Real const> const& vz,
             amrex::Array4<amrex::Real const> const& px,
             amrex::Array4<amrex::Real const> const& pz,
             amrex::Real tdtdz)
{
    pxz(i,j,k) = ( (vx(i,j,k) < 0) ?
                 px(i,j,k) - tdtdz * ( 0.5*(vz(i,  j,k+1) + vz(i  ,j,k)) * (pz(i  ,j,k+1) - pz(i  ,j,k))) :
                 px(i,j,k) - tdtdz * ( 0.5*(vz(i-1,j,k+1) + vz(i-1,j,k)) * (pz(i-1,j,k+1) - pz(i-1,j,k))) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_yx(int i, int j, int k,
             amrex::Array4<amrex::Real> const& pyx,
             amrex::Array4<amrex::Real const> const& vx,
             amrex::Array4<amrex::Real const> const& vy,
             amrex::Array4<amrex::Real const> const& px,
             amrex::Array4<amrex::Real const> const& py,
             amrex::Real tdtdx)
{
    pyx(i,j,k) = ( (vy(i,j,k) < 0) ?
                 py(i,j,k) - tdtdx * ( 0.5*(vx(i+1,j  ,k) + vx(i,j  ,k)) * (px(i+1,j  ,k) - px(i,j  ,k))) :
                 py(i,j,k) - tdtdx * ( 0.5*(vx(i+1,j-1,k) + vx(i,j-1,k)) * (px(i+1,j-1,k) - px(i,j-1,k))) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_yz(int i, int j, int k,
             amrex::Array4<amrex::Real> const& pyz,
             amrex::Array4<amrex::Real const> const& vy,
             amrex::Array4<amrex::Real const> const& vz,
             amrex::Array4<amrex::Real const> const& py,
             amrex::Array4<amrex::Real const> const& pz,
             amrex::Real tdtdz)
{
    pyz(i,j,k) = ( (vy(i,j,k) < 0) ?
                 py(i,j,k) - tdtdz * ( 0.5*(vz(i,  j,k+1) + vz(i,j  ,k)) * (pz(i,j  ,k+1) - pz(i,j  ,k))) :
                 py(i,j,k) - tdtdz * ( 0.5*(vz(i,j-1,k+1) + vz(i,j-1,k)) * (pz(i,j-1,k+1) - pz(i,j-1,k))) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_zx(int i, int j, int k,
             amrex::Array4<amrex::Real> const& pzx,
             amrex::Array4<amrex::Real const> const& vx,
             amrex::Array4<amrex::Real const> const& vz,
             amrex::Array4<amrex::Real const> const& px,
             amrex::Array4<amrex::Real const> const& pz,
             amrex::Real tdtdx)
{
    pzx(i,j,k) = ( (vz(i,j,k) < 0) ?
                 pz(i,j,k) - tdtdx * ( 0.5*(vx(i+1,j,k  ) + vx(i,j,k  )) * (px(i+1,j,k  ) - px(i,j,k  ))) :
                 pz(i,j,k) - tdtdx * ( 0.5*(vx(i+1,j,k-1) + vx(i,j,k-1)) * (px(i+1,j,k-1) - px(i,j,k-1))) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void flux_zy(int i, int j, int k,
             amrex::Array4<amrex::Real> const& pzy,
             amrex::Array4<amrex::Real const> const& vy,
             amrex::Array4<amrex::Real const> const& vz,
             amrex::Array4<amrex::Real const> const& py,
             amrex::Array4<amrex::Real const> const& pz,
             amrex::Real tdtdy)
{
    pzy(i,j,k) = ( (vz(i,j,k) < 0) ?
                 pz(i,j,k) - tdtdy * ( 0.5*(vy(i,j+1,k  ) + vy(i,j,k  )) * (py(i,j+1,k  ) - py(i,j,k  ))) :
                 pz(i,j,k) - tdtdy * ( 0.5*(vy(i,j+1,k-1) + vy(i,j,k-1)) * (py(i,j+1,k-1) - py(i,j,k-1))) );
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void create_flux_x(int i, int j, int k,
                   amrex::Array4<amrex::Real> const& fx,
                   amrex::Array4<amrex::Real const> const& vx,
                   amrex::Array4<amrex::Real const> const& vy,
                   amrex::Array4<amrex::Real const> const& vz,
                   amrex::Array4<amrex::Real const> const& px,
                   amrex::Array4<amrex::Real const> const& pyz,
                   amrex::Array4<amrex::Real const> const& pzy,
                   amrex::Real hdtdy, amrex::Real hdtdz)
{
    amrex::Real f = ( (vx(i,j,k) < 0) ?
                px(i,j,k) - hdtdy * ( 0.5*(vy(i  ,j+1,k  ) + vy(i  ,j,k)) * (pyz(i  ,j+1,k  )-pyz(i  ,j,k)))
                          - hdtdz * ( 0.5*(vz(i  ,j  ,k+1) + vz(i  ,j,k)) * (pzy(i  ,j  ,k+1)-pzy(i  ,j,k))) :
                px(i,j,k) - hdtdy * ( 0.5*(vy(i-1,j+1,k  ) + vy(i-1,j,k)) * (pyz(i-1,j+1,k  )-pyz(i-1,j,k)))
                          - hdtdz * ( 0.5*(vz(i-1,j  ,k+1) + vz(i-1,j,k)) * (pzy(i-1,j  ,k+1)-pzy(i-1,j,k))) );

    fx(i,j,k) = vx(i,j,k)*f;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void create_flux_y(int i, int j, int k,
                   amrex::Array4<amrex::Real> const& fy,
                   amrex::Array4<amrex::Real const> const& vx,
                   amrex::Array4<amrex::Real const> const& vy,
                   amrex::Array4<amrex::Real const> const& vz,
                   amrex::Array4<amrex::Real const> const& py,
                   amrex::Array4<amrex::Real const> const& pxz,
                   amrex::Array4<amrex::Real const> const& pzx,
                   amrex::Real hdtdx, amrex::Real hdtdz)
{
    amrex::Real f = ( (vy(i,j,k) < 0) ?
                py(i,j,k) - hdtdx * ( 0.5*(vx(i+1,j  ,k  ) + vx(i,j  ,k)) * (pxz(i+1,j  ,k  )-pxz(i,j  ,k)))
                          - hdtdz * ( 0.5*(vz(i,  j  ,k+1) + vz(i,j  ,k)) * (pzx(i,  j  ,k+1)-pzx(i,j  ,k))) :
                py(i,j,k) - hdtdx * ( 0.5*(vx(i+1,j-1,k  ) + vx(i,j-1,k)) * (pxz(i+1,j-1,k  )-pxz(i,j-1,k)))
                          - hdtdz * ( 0.5*(vz(i  ,j-1,k+1) + vz(i,j-1,k)) * (pzx(i  ,j-1,k+1)-pzx(i,j-1,k))) );

    fy(i,j,k) = vy(i,j,k)*f;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void create_flux_z(int i, int j, int k,
                   amrex::Array4<amrex::Real> const& fz,
                   amrex::Array4<amrex::Real const> const& vx,
                   amrex::Array4<amrex::Real const> const& vy,
                   amrex::Array4<amrex::Real const> const& vz,
                   amrex::Array4<amrex::Real const> const& pz,
                   amrex::Array4<amrex::Real const> const& pxy,
                   amrex::Array4<amrex::Real const> const& pyx,
                   amrex::Real hdtdx, amrex::Real hdtdy)
{
    amrex::Real f = ( (vz(i,j,k) < 0) ?
                pz(i,j,k) - hdtdx * ( 0.5*(vx(i+1,j  ,k  ) + vx(i,j,k  )) * (pxy(i+1,j  ,k  )-pxy(i,j,k  )))
                          - hdtdy * ( 0.5*(vy(i,  j+1,k  ) + vy(i,j,k  )) * (pyx(i,  j+1,k  )-pyx(i,j,k  ))) :
                pz(i,j,k) - hdtdx * ( 0.5*(vx(i+1,j  ,k-1) + vx(i,j,k-1)) * (pxy(i+1,j  ,k-1)-pxy(i,j,k-1)))
                          - hdtdy * ( 0.5*(vy(i  ,j+1,k-1) + vy(i,j,k-1)) * (pyx(i  ,j+1,k-1)-pyx(i,j,k-1))) );

    fz(i,j,k) = vz(i,j,k)*f;
}

#endif /* flux_3d_K_H_ */
