#ifndef _compute_flux_3d_cpp_
#define _compute_flux_3d_cpp_

#include <AmrCoreAdv_F.H>

using namespace amrex;

void compute_flux_3d(Box const& bx, 
                     GpuArray<Real, AMREX_SPACEDIM>& dtdx, 
                     const FArrayBox& phifab, 
                     AMREX_D_DECL(const FArrayBox& velx, 
                                  const FArrayBox& vely,
                                  const FArrayBox& velz),
                     AMREX_D_DECL(FArrayBox& flxx,
                                  FArrayBox& flxy,
                                  FArrayBox& flxz))
{
    GpuArray<Real, AMREX_SPACEDIM> hdtdx;
    GpuArray<Real, AMREX_SPACEDIM> tdtdx;
    for (int i=0; i<AMREX_SPACEDIM; ++i)
    {
        hdtdx[i] = 0.5*dtdx[i];
        tdtdx[i] = dtdx[i] / 3.0;
    }

    const Box nbx = amrex::grow(bx, 1);

    Gpu::AsyncFab fabslope (nbx, 1);
    FArrayBox* slope = fabslope.fabPtr();

/*
    call slopex(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)
*/

    Gpu::AsyncFab fabphix (nbx, 1);
    Gpu::AsyncFab fabphiy (nbx, 1);
    Gpu::AsyncFab fabphiz (nbx, 1);
    FArrayBox* phix = fabphix.fabPtr();
    FArrayBox* phiy = fabphiy.fabPtr();
    FArrayBox* phiz = fabphiz.fabPtr();

    const Box xfbx = amrex::growLo(nbx, 0, -1);
    const Box yfbx = amrex::growLo(nbx, 1, -1);
    const Box zfbx = amrex::growLo(nbx, 2, -1);

    AMREX_LAUNCH_DEVICE_LAMBDA(xfbx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);

        const auto vx  = velx.view(lo);
        const auto phi = phifab.view(lo); 
        const auto px  = phix->view(lo);
        const auto del = slope->view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    px(i,j,k) = ( (vx(i,j,k) < 0) ? 
                                   phi(i  ,j,k) - del(i  ,j,k)*(0.5 + hdtdx[0]*vx(i,j,k)) : 
                                   phi(i-1,j,k) - del(i-1,j,k)*(0.5 + hdtdx[0]*vx(i,j,k)) );
                }
            }
        }
    });

/*
    call slopey(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)
*/                

    AMREX_LAUNCH_DEVICE_LAMBDA(yfbx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);

        const auto vy  = vely.view(lo);
        const auto phi = phifab.view(lo); 
        const auto py  = phiy->view(lo);
        const auto del = slope->view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    py(i,j,k) = ( (vy(i,j,k) < 0) ? 
                                   phi(i,j  ,k) - del(i,j  ,k)*(0.5 + hdtdx[0]*vy(i,j,k)) : 
                                   phi(i,j-1,k) - del(i,j-1,k)*(0.5 + hdtdx[0]*vy(i,j,k)) );
                }
            }
        }
    });

/*
    call slopez(glo, ghi, &
                phi, ph_lo, ph_hi, &
                slope, glo, ghi)
*/                

    AMREX_LAUNCH_DEVICE_LAMBDA(zfbx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);

        const auto vz  = velz.view(lo);
        const auto phi = phifab.view(lo); 
        const auto pz  = phiz->view(lo);
        const auto del = slope->view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    pz(i,j,k) = ( (vz(i,j,k) < 0) ? 
                                   phi(i,j,k  ) - del(i,j,k  )*(0.5 + hdtdx[0]*vz(i,j,k)) : 
                                   phi(i,j,k-1) - del(i,j,k-1)*(0.5 + hdtdx[0]*vz(i,j,k)) );
                }
            }
        }
    });

//    !!!!!!!!!!!!!!!!!!!!
//    ! transverse terms
//    !!!!!!!!!!!!!!!!!!!!

    Gpu::AsyncFab fabphix_y (nbx, 1);
    Gpu::AsyncFab fabphix_z (nbx, 1);
    FArrayBox* phix_y = fabphix_y.fabPtr();
    FArrayBox* phix_z = fabphix_z.fabPtr();

    Box xyfbx = amrex::grow(bx, 2, 1);
    xyfbx.growHi(0, 1); 
    Box xzfbx = amrex::grow(bx, 1, 1);
    xzfbx.growHi(0, 1); 

    // update phi on x faces by adding in y-transverse terms
    AMREX_LAUNCH_DEVICE_LAMBDA(xyfbx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);

        const auto vx  = velx.view(lo);
        const auto vy  = vely.view(lo);
        const auto px  = phix->view(lo);
        const auto py  = phiy->view(lo);
        const auto pxy = phix_y->view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    pxy(i,j,k) = ( (vx(i,j,k) < 0) ? 
                                   px(i,j,k) - tdtdx[1] * ( 0.5*(vy(i,  j+1,k) + vy(i  ,j,k)) * (py(i  ,j+1,k)-py(i  ,j,k))) : 
                                   px(i,j,k) - tdtdx[1] * ( 0.5*(vy(i-1,j+1,k) + vy(i-1,j,k)) * (py(i-1,j+1,k)-py(i-1,j,k))) );
                }
            }
        }
    });

    // update phi on x faces by adding in z-transverse terms
    AMREX_LAUNCH_DEVICE_LAMBDA(xzfbx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);

        const auto vx  = velx.view(lo);
        const auto vz  = velz.view(lo);
        const auto px  = phix->view(lo);
        const auto pz  = phiz->view(lo);
        const auto pxz = phix_z->view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    pxz(i,j,k) = ( (vx(i,j,k) < 0) ? 
                                   px(i,j,k) - tdtdx[2] * ( 0.5*(vz(i,  j,k+1) + vz(i  ,j,k)) * (pz(i  ,j,k+1)-pz(i  ,j,k))) : 
                                   px(i,j,k) - tdtdx[2] * ( 0.5*(vz(i-1,j,k+1) + vz(i-1,j,k)) * (pz(i-1,j,k+1)-pz(i-1,j,k))) );
                }
            }
        }
    });

    Gpu::AsyncFab fabphiy_x (nbx, 1);
    Gpu::AsyncFab fabphiy_z (nbx, 1);
    FArrayBox* phiy_x = fabphiy_x.fabPtr();
    FArrayBox* phiy_z = fabphiy_z.fabPtr();

    Box yxfbx = amrex::grow(bx, 2, 1);
    yxfbx.growHi(1, 1); 
    Box yzfbx = amrex::grow(bx, 0, 1);
    yzfbx.growHi(1, 1); 

    // update phi on y faces by adding in x-transverse terms
    AMREX_LAUNCH_DEVICE_LAMBDA(yxfbx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);

        const auto vx  = velx.view(lo);
        const auto vy  = vely.view(lo);
        const auto py  = phiy->view(lo);
        const auto px  = phix->view(lo);
        const auto pyx = phiy_x->view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    pyx(i,j,k) = ( (vy(i,j,k) < 0) ? 
                                   py(i,j,k) - tdtdx[0] * ( 0.5*(vx(i+1,j  ,k) + vx(i,j  ,k)) * (px(i+1,j  ,k)-px(i,j  ,k))) : 
                                   py(i,j,k) - tdtdx[0] * ( 0.5*(vx(i+1,j-1,k) + vx(i,j-1,k)) * (px(i+1,j-1,k)-px(i,j-1,k))) );
                }
            }
        }
    });

    // update phi on y faces by adding in z-transverse terms
    AMREX_LAUNCH_DEVICE_LAMBDA(yzfbx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);

        const auto vz  = velz.view(lo);
        const auto vy  = vely.view(lo);
        const auto py  = phiy->view(lo);
        const auto pz  = phiz->view(lo);
        const auto pyz = phiy_z->view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    pyz(i,j,k) = ( (vy(i,j,k) < 0) ? 
                                   py(i,j,k) - tdtdx[2] * ( 0.5*(vz(i,  j,k+1) + vz(i,j  ,k)) * (pz(i,j  ,k+1)-pz(i,j  ,k))) : 
                                   py(i,j,k) - tdtdx[2] * ( 0.5*(vz(i,j-1,k+1) + vz(i,j-1,k)) * (pz(i,j-1,k+1)-pz(i,j-1,k))) );
                }
            }
        }
    });

    Gpu::AsyncFab fabphiz_x (nbx, 1);
    Gpu::AsyncFab fabphiz_y (nbx, 1);
    FArrayBox* phiz_x = fabphiz_x.fabPtr();
    FArrayBox* phiz_y = fabphiz_y.fabPtr();

    Box zxfbx = amrex::grow(bx, 1, 1);
    zxfbx.growHi(2, 1); 
    Box zyfbx = amrex::grow(bx, 0, 1);
    zyfbx.growHi(2, 1); 

    // update phi on z faces by adding in x-transverse terms
    AMREX_LAUNCH_DEVICE_LAMBDA(zxfbx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);

        const auto vx  = velx.view(lo);
        const auto vz  = velz.view(lo);
        const auto pz  = phiz->view(lo);
        const auto px  = phix->view(lo);
        const auto pzx = phiz_x->view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    pzx(i,j,k) = ( (vz(i,j,k) < 0) ? 
                                   pz(i,j,k) - tdtdx[0] * ( 0.5*(vx(i+1,j,k  ) + vx(i,j,k  )) * (px(i+1,j,k  )-px(i,j,k  ))) : 
                                   pz(i,j,k) - tdtdx[0] * ( 0.5*(vx(i+1,j,k-1) + vx(i,j,k-1)) * (px(i+1,j,k-1)-px(i,j,k-1))) );
                }
            }
        }
    });

    // update phi on z faces by adding in y-transverse terms
    AMREX_LAUNCH_DEVICE_LAMBDA(zyfbx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);

        const auto vx  = velx.view(lo);
        const auto vz  = velz.view(lo);
        const auto pz  = phiz->view(lo);
        const auto py  = phiy->view(lo);
        const auto pzy = phiz_y->view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    pzy(i,j,k) = ( (vz(i,j,k) < 0) ? 
                                   pz(i,j,k) - tdtdx[1] * ( 0.5*(vx(i,j+1,k  ) + vx(i,j,k-1)) * (py(i,j+1,k  )-py(i,j,k  ))) : 
                                   pz(i,j,k) - tdtdx[1] * ( 0.5*(vx(i,j+1,k-1) + vx(i,j,k  )) * (py(i,j+1,k-1)-py(i,j,k-1))) );
                }
            }
        }
    });

//    !!!!!!!!!!!!!!!!!!!!
//    ! final edge states
//    !!!!!!!!!!!!!!!!!!!!

    const Box xbx = amrex::growHi(bx, 0, 1);
    const Box ybx = amrex::growHi(bx, 1, 1);
    const Box zbx = amrex::growHi(bx, 2, 1);

    // update phi on x faces by adding in yz and zy transverse terms
    AMREX_LAUNCH_DEVICE_LAMBDA(xbx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);

        const auto vx  = velx.view(lo);
        const auto vy  = vely.view(lo);
        const auto vz  = velz.view(lo);
        const auto px  = phix->view(lo);
        const auto pyz = phiy_z->view(lo);
        const auto pzy = phiz_y->view(lo);
        const auto fx  = flxx.view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    px(i,j,k) = ( (vx(i,j,k) < 0) ? 
                                  px(i,j,k) - hdtdx[1] * ( 0.5*(vy(i,j+1,k  ) + vy(i,j,k)) * (pyz(i,j+1,k  )-pyz(i,j,k)))
                                - px(i,j,k) - hdtdx[2] * ( 0.5*(vz(i,j  ,k+1) + vz(i,j,k)) * (pzy(i,j  ,k+1)-pzy(i,j,k)))  :
                                  px(i,j,k) - hdtdx[1] * ( 0.5*(vy(i-1,j+1,k  ) + vy(i-1,j,k)) * (pyz(i-1,j+1,k  )-pyz(i-1,j,k)))
                                - px(i,j,k) - hdtdx[2] * ( 0.5*(vz(i-1,j  ,k+1) + vz(i-1,j,k)) * (pzy(i-1,j  ,k+1)-pzy(i-1,j,k))) );

                    fx(i,j,k) = vx(i,j,k)*px(i,j,k);
                }
            }
        }
    });


    // update phi on y faces by adding in xz and zx transverse terms
    AMREX_LAUNCH_DEVICE_LAMBDA(ybx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);

        const auto vx  = velx.view(lo);
        const auto vy  = vely.view(lo);
        const auto vz  = velz.view(lo);
        const auto py  = phiy->view(lo);
        const auto pxz = phix_z->view(lo);
        const auto pzx = phiz_x->view(lo);
        const auto fy  = flxy.view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    py(i,j,k) = ( (vy(i,j,k) < 0) ? 
                                  py(i,j,k) - hdtdx[0] * ( 0.5*(vx(i+1,j,k  ) + vx(i,j,k)) * (pxz(i+1,j,k  )-pxz(i,j,k)))
                                - py(i,j,k) - hdtdx[2] * ( 0.5*(vz(i,  j,k+1) + vz(i,j,k)) * (pzx(i,  j,k+1)-pzx(i,j,k)))  :
                                  py(i,j,k) - hdtdx[0] * ( 0.5*(vy(i+1,j-1,k  ) + vy(i,j-1,k)) * (pxz(i+1,j-1,k  )-pxz(i,j-1,k)))
                                - py(i,j,k) - hdtdx[2] * ( 0.5*(vz(i  ,j-1,k+1) + vz(i,j-1,k)) * (pzx(i  ,j-1,k+1)-pzx(i,j-1,k))) );

                    fy(i,j,k) = vy(i,j,k)*py(i,j,k);
                }
            }
        }
    });

    // update phi on z faces by adding in xy and yx transverse terms
    AMREX_LAUNCH_DEVICE_LAMBDA(zbx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);

        const auto vx  = velx.view(lo);
        const auto vy  = vely.view(lo);
        const auto vz  = velz.view(lo);
        const auto pz  = phiz->view(lo);
        const auto pxy = phix_y->view(lo);
        const auto pyx = phiy_x->view(lo);
        const auto fz  = flxz.view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    pz(i,j,k) = ( (vz(i,j,k) < 0) ? 
                                  pz(i,j,k) - hdtdx[0] * ( 0.5*(vx(i+1,j  ,k) + vx(i,j,k)) * (pxy(i+1,j  ,k)-pxy(i,j,k)))
                                - pz(i,j,k) - hdtdx[1] * ( 0.5*(vy(i,  j+1,k) + vy(i,j,k)) * (pyx(i,  j+1,k)-pyx(i,j,k)))  :
                                  pz(i,j,k) - hdtdx[0] * ( 0.5*(vx(i+1,j  ,k-1) + vx(i,j,k-1)) * (pxy(i+1,j  ,k-1)-pxy(i,j,k-1)))
                                - pz(i,j,k) - hdtdx[1] * ( 0.5*(vy(i  ,j+1,k-1) + vy(i,j,k-1)) * (pyx(i  ,j+1,k-1)-pyx(i,j,k-1))) );

                    fz(i,j,k) = vz(i,j,k)*pz(i,j,k);
                }
            }
        }
    });

    fabphix.clear();
    fabphix_y.clear();
    fabphix_z.clear();
    fabphiy.clear();
    fabphiy_x.clear();
    fabphiy_z.clear();
    fabphiz.clear();
    fabphiz_x.clear();
    fabphiz_y.clear();
    fabslope.clear();

    // ADD SYNCHS AS NEEDED

}

#endif
