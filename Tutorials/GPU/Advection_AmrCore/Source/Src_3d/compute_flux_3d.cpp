
#include <AmrCoreAdv_F.H>

using namespace amrex;

AMREX_GPU_DEVICE
void flux_x(Box const& bx, 
            const FArrayBox& state,
            const FArrayBox& velx,
            FArrayBox& phix,
            const FArrayBox& slope,
            const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);

    const auto vx  = velx.view(lo);
    const auto phi = state.view(lo); 
    const auto px  = phix.view(lo);
    const auto del = slope.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                px(i,j,k) = ( (vx(i,j,k) < 0) ? 
                               phi(i  ,j,k) - del(i  ,j,k)*(0.5 + 0.5*dtdx[0]*vx(i,j,k)) : 
                               phi(i-1,j,k) + del(i-1,j,k)*(0.5 - 0.5*dtdx[0]*vx(i,j,k)) );
            }
        }
    }
}

AMREX_GPU_DEVICE
void flux_y(Box const& bx, 
            const FArrayBox& state,
            const FArrayBox& vely,
            FArrayBox& phiy,
            const FArrayBox& slope,
            const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);

    const auto vy  = vely.view(lo);
    const auto phi = state.view(lo); 
    const auto py  = phiy.view(lo);
    const auto del = slope.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                py(i,j,k) = ( (vy(i,j,k) < 0) ? 
                               phi(i,j  ,k) - del(i,j  ,k)*(0.5 + 0.5*dtdx[0]*vy(i,j,k)) : 
                               phi(i,j-1,k) + del(i,j-1,k)*(0.5 - 0.5*dtdx[0]*vy(i,j,k)) );
            }
        }
    }
}

AMREX_GPU_DEVICE
void flux_z(Box const& bx, 
            const FArrayBox& state,
            const FArrayBox& velz,
            FArrayBox& phiz,
            const FArrayBox& slope,
            const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);

    const auto vz  = velz.view(lo);
    const auto phi = state.view(lo); 
    const auto pz  = phiz.view(lo);
    const auto del = slope.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                pz(i,j,k) = ( (vz(i,j,k) < 0) ? 
                               phi(i,j,k  ) - del(i,j,k  )*(0.5 + 0.5*dtdx[0]*vz(i,j,k)) : 
                               phi(i,j,k-1) + del(i,j,k-1)*(0.5 - 0.5*dtdx[0]*vz(i,j,k)) );
            }
        }
    }
}

AMREX_GPU_DEVICE
void flux_xy(Box const& bx,
             AMREX_D_DECL(const FArrayBox& velx, 
                          const FArrayBox& vely, 
                          const FArrayBox& velz),
             AMREX_D_DECL(const FArrayBox& phix,
                          const FArrayBox& phiy,
                          const FArrayBox& phiz),
             FArrayBox& phix_y,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);

    const auto vx  = velx.view(lo);
    const auto vy  = vely.view(lo);
    const auto px  = phix.view(lo);
    const auto py  = phiy.view(lo);
    const auto pxy = phix_y.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                pxy(i,j,k) = ( (vx(i,j,k) < 0) ? 
                               px(i,j,k) - dtdx[1]/3.0 * ( 0.5*(vy(i,  j+1,k) + vy(i  ,j,k)) * (py(i  ,j+1,k) - py(i  ,j,k))) : 
                               px(i,j,k) - dtdx[1]/3.0 * ( 0.5*(vy(i-1,j+1,k) + vy(i-1,j,k)) * (py(i-1,j+1,k) - py(i-1,j,k))) );
            }
        }
    }
}

AMREX_GPU_DEVICE
void flux_xz(Box const& bx, 
             AMREX_D_DECL(const FArrayBox& velx, 
                          const FArrayBox& vely, 
                          const FArrayBox& velz),
             AMREX_D_DECL(const FArrayBox& phix,
                          const FArrayBox& phiy,
                          const FArrayBox& phiz),
             FArrayBox& phix_z,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);

    const auto vx  = velx.view(lo);
    const auto vz  = velz.view(lo);
    const auto px  = phix.view(lo);
    const auto pz  = phiz.view(lo);
    const auto pxz = phix_z.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                pxz(i,j,k) = ( (vx(i,j,k) < 0) ? 
                               px(i,j,k) - dtdx[2]/3.0 * ( 0.5*(vz(i,  j,k+1) + vz(i  ,j,k)) * (pz(i  ,j,k+1) - pz(i  ,j,k))) : 
                               px(i,j,k) - dtdx[2]/3.0 * ( 0.5*(vz(i-1,j,k+1) + vz(i-1,j,k)) * (pz(i-1,j,k+1) - pz(i-1,j,k))) );
            }
        }
    }
}

AMREX_GPU_DEVICE
void flux_yx(Box const& bx, 
             AMREX_D_DECL(const FArrayBox& velx, 
                          const FArrayBox& vely, 
                          const FArrayBox& velz),
             AMREX_D_DECL(const FArrayBox& phix,
                          const FArrayBox& phiy,
                          const FArrayBox& phiz),
             FArrayBox& phiy_x,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);

    const auto vx  = velx.view(lo);
    const auto vy  = vely.view(lo);
    const auto py  = phiy.view(lo);
    const auto px  = phix.view(lo);
    const auto pyx = phiy_x.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                pyx(i,j,k) = ( (vy(i,j,k) < 0) ? 
                               py(i,j,k) - dtdx[0]/3.0 * ( 0.5*(vx(i+1,j  ,k) + vx(i,j  ,k)) * (px(i+1,j  ,k) - px(i,j  ,k))) : 
                               py(i,j,k) - dtdx[0]/3.0 * ( 0.5*(vx(i+1,j-1,k) + vx(i,j-1,k)) * (px(i+1,j-1,k) - px(i,j-1,k))) );
            }
        }
    }
};

AMREX_GPU_DEVICE
void flux_yz(Box const& bx,
             AMREX_D_DECL(const FArrayBox& velx,
                          const FArrayBox& vely,
                          const FArrayBox& velz),
             AMREX_D_DECL(const FArrayBox& phix,
                          const FArrayBox& phiy,
                          const FArrayBox& phiz),
             FArrayBox& phiy_z,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);

    const auto vz  = velz.view(lo);
    const auto vy  = vely.view(lo);
    const auto py  = phiy.view(lo);
    const auto pz  = phiz.view(lo);
    const auto pyz = phiy_z.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                pyz(i,j,k) = ( (vy(i,j,k) < 0) ? 
                               py(i,j,k) - dtdx[2]/3.0 * ( 0.5*(vz(i,  j,k+1) + vz(i,j  ,k)) * (pz(i,j  ,k+1) - pz(i,j  ,k))) : 
                               py(i,j,k) - dtdx[2]/3.0 * ( 0.5*(vz(i,j-1,k+1) + vz(i,j-1,k)) * (pz(i,j-1,k+1) - pz(i,j-1,k))) );
            }
        }
    }
}

AMREX_GPU_DEVICE
void flux_zx(Box const& bx,
             AMREX_D_DECL(const FArrayBox& velx,
                          const FArrayBox& vely,
                          const FArrayBox& velz),
             AMREX_D_DECL(const FArrayBox& phix,
                          const FArrayBox& phiy,
                          const FArrayBox& phiz),
             FArrayBox& phiz_x,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);

    const auto vx  = velx.view(lo);
    const auto vz  = velz.view(lo);
    const auto pz  = phiz.view(lo);
    const auto px  = phix.view(lo);
    const auto pzx = phiz_x.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                pzx(i,j,k) = ( (vz(i,j,k) < 0) ? 
                               pz(i,j,k) - dtdx[0]/3.0 * ( 0.5*(vx(i+1,j,k  ) + vx(i,j,k  )) * (px(i+1,j,k  ) - px(i,j,k  ))) : 
                               pz(i,j,k) - dtdx[0]/3.0 * ( 0.5*(vx(i+1,j,k-1) + vx(i,j,k-1)) * (px(i+1,j,k-1) - px(i,j,k-1))) );
            }
        }
    }
}

AMREX_GPU_DEVICE
void flux_zy(Box const& bx,
             AMREX_D_DECL(const FArrayBox& velx,
                          const FArrayBox& vely,
                          const FArrayBox& velz),
             AMREX_D_DECL(const FArrayBox& phix,
                          const FArrayBox& phiy,
                          const FArrayBox& phiz),
             FArrayBox& phiz_y,
             const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{

    const auto len = length(bx);
    const auto lo  = lbound(bx);

    const auto vy  = vely.view(lo);
    const auto vz  = velz.view(lo);
    const auto pz  = phiz.view(lo);
    const auto py  = phiy.view(lo);
    const auto pzy = phiz_y.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {

                pzy(i,j,k) = ( (vz(i,j,k) < 0) ? 
                               pz(i,j,k) - dtdx[1]/3.0 * ( 0.5*(vy(i,j+1,k  ) + vy(i,j,k  )) * (py(i,j+1,k  ) - py(i,j,k  ))) : 
                               pz(i,j,k) - dtdx[1]/3.0 * ( 0.5*(vy(i,j+1,k-1) + vy(i,j,k-1)) * (py(i,j+1,k-1) - py(i,j,k-1))) );

            }
        }
    }
}

AMREX_GPU_DEVICE
void combine_flux_x(Box const& bx,
                    const FArrayBox& velx,
                    const FArrayBox& vely,
                    const FArrayBox& velz,
                    FArrayBox& phix,
                    const FArrayBox& phiy_z,
                    const FArrayBox& phiz_y,
                    FArrayBox& flxx,
                    const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);

    const auto vx  = velx.view(lo);
    const auto vy  = vely.view(lo);
    const auto vz  = velz.view(lo);
    const auto px  = phix.view(lo);
    const auto pyz = phiy_z.view(lo);
    const auto pzy = phiz_y.view(lo);
    const auto fx  = flxx.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {

                px(i,j,k) = ( (vx(i,j,k) < 0) ? 
                              px(i,j,k) - 0.5*dtdx[1] * ( 0.5*(vy(i  ,j+1,k  ) + vy(i  ,j,k)) * (pyz(i  ,j+1,k  )-pyz(i  ,j,k)))
                                        - 0.5*dtdx[2] * ( 0.5*(vz(i  ,j  ,k+1) + vz(i  ,j,k)) * (pzy(i  ,j  ,k+1)-pzy(i  ,j,k))) :
                              px(i,j,k) - 0.5*dtdx[1] * ( 0.5*(vy(i-1,j+1,k  ) + vy(i-1,j,k)) * (pyz(i-1,j+1,k  )-pyz(i-1,j,k)))
                                        - 0.5*dtdx[2] * ( 0.5*(vz(i-1,j  ,k+1) + vz(i-1,j,k)) * (pzy(i-1,j  ,k+1)-pzy(i-1,j,k))) );

                fx(i,j,k) = vx(i,j,k)*px(i,j,k);

            }
        }
    }
}

AMREX_GPU_DEVICE
void combine_flux_y(Box const& bx,
                    const FArrayBox& velx,
                    const FArrayBox& vely,
                    const FArrayBox& velz,
                    FArrayBox& phiy,
                    const FArrayBox& phix_z,
                    const FArrayBox& phiz_x,
                    FArrayBox& flxy,
                    const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);

    const auto vx  = velx.view(lo);
    const auto vy  = vely.view(lo);
    const auto vz  = velz.view(lo);
    const auto py  = phiy.view(lo);
    const auto pxz = phix_z.view(lo);
    const auto pzx = phiz_x.view(lo);
    const auto fy  = flxy.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                py(i,j,k) = ( (vy(i,j,k) < 0) ? 
                              py(i,j,k) - 0.5*dtdx[0] * ( 0.5*(vx(i+1,j  ,k  ) + vx(i,j  ,k)) * (pxz(i+1,j  ,k  )-pxz(i,j  ,k)))
                                        - 0.5*dtdx[2] * ( 0.5*(vz(i,  j  ,k+1) + vz(i,j  ,k)) * (pzx(i,  j  ,k+1)-pzx(i,j  ,k))) :
                              py(i,j,k) - 0.5*dtdx[0] * ( 0.5*(vx(i+1,j-1,k  ) + vx(i,j-1,k)) * (pxz(i+1,j-1,k  )-pxz(i,j-1,k)))
                                        - 0.5*dtdx[2] * ( 0.5*(vz(i  ,j-1,k+1) + vz(i,j-1,k)) * (pzx(i  ,j-1,k+1)-pzx(i,j-1,k))) );

                fy(i,j,k) = vy(i,j,k)*py(i,j,k);
            }
        }
    }
}

AMREX_GPU_DEVICE
void combine_flux_z(Box const& bx,
                    const FArrayBox& velx,
                    const FArrayBox& vely,
                    const FArrayBox& velz,
                    FArrayBox& phiz,
                    const FArrayBox& phix_y,
                    const FArrayBox& phiy_x,
                    FArrayBox& flxz,
                    const GpuArray<Real, AMREX_SPACEDIM>& dtdx)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);

    const auto vx  = velx.view(lo);
    const auto vy  = vely.view(lo);
    const auto vz  = velz.view(lo);
    const auto pz  = phiz.view(lo);
    const auto pxy = phix_y.view(lo);
    const auto pyx = phiy_x.view(lo);
    const auto fz  = flxz.view(lo);

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            for (int i = 0; i < len.x; ++i) {
                pz(i,j,k) = ( (vz(i,j,k) < 0) ? 
                              pz(i,j,k) - 0.5*dtdx[0] * ( 0.5*(vx(i+1,j  ,k  ) + vx(i,j,k  )) * (pxy(i+1,j  ,k  )-pxy(i,j,k  )))
                                        - 0.5*dtdx[1] * ( 0.5*(vy(i,  j+1,k  ) + vy(i,j,k  )) * (pyx(i,  j+1,k  )-pyx(i,j,k  ))) :
                              pz(i,j,k) - 0.5*dtdx[0] * ( 0.5*(vx(i+1,j  ,k-1) + vx(i,j,k-1)) * (pxy(i+1,j  ,k-1)-pxy(i,j,k-1)))
                                        - 0.5*dtdx[1] * ( 0.5*(vy(i  ,j+1,k-1) + vy(i,j,k-1)) * (pyx(i  ,j+1,k-1)-pyx(i,j,k-1))) );

                fz(i,j,k) = vz(i,j,k)*pz(i,j,k);
            }
        }
    }
}

