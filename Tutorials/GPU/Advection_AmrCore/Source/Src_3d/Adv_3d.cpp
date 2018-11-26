#ifndef _Adv_3d_cpp_
#define _Adv_3d_cpp_

#include <AmrCoreAdv_F.H>
#include <AMReX_FArrayBox.H>

using namespace amrex;

void advect(const Real time, Box const& bx,
            const FArrayBox& statein,
            FArrayBox& stateout,
            AMREX_D_DECL(const FArrayBox& xvel,
                         const FArrayBox& yvel,
                         const FArrayBox& zvel),
            AMREX_D_DECL(FArrayBox& fx,
                         FArrayBox& fy,
                         FArrayBox& fz),
            const GpuArray<Real, AMREX_SPACEDIM>& dx, const Real dt)
{

    GpuArray<Real, AMREX_SPACEDIM> dtdx; 
    for (int i=0; i<AMREX_SPACEDIM; ++i)
    {
        dtdx[i] = dt/dx[i]; 
    }

    const Box nbx = amrex::grow(bx, 1);

/*  NOT USED HERE. MOVE INTO COMPUTE_FLUX_3D?
    Gpu::DeviceFab fphix   (nbx, 1);
    Gpu::DeviceFab fphix_y (nbx, 1);
    Gpu::DeviceFab fphix_z (nbx, 1);
    Gpu::DeviceFab fphiy   (nbx, 1);
    Gpu::DeviceFab fphiy_x (nbx, 1);
    Gpu::DeviceFab fphiy_z (nbx, 1);
    Gpu::DeviceFab fphiz   (nbx, 1);
    Gpu::DeviceFab fphiz_x (nbx, 1);
    Gpu::DeviceFab fphiz_y (nbx, 1);
    Gpu::DeviceFab fslope  (nbx, 1);

    FArrayBox* phix   = fphix.fabPtr();
    FArrayBox* phix_y = fphix_y.fabPtr();
    FArrayBox* phix_z = fphix_z.fabPtr();
    FArrayBox* phiy   = fphiy.fabPtr();
    FArrayBox* phiy_x = fphiy_x.fabPtr();
    FArrayBox* phiy_z = fphiy_z.fabPtr();
    FArrayBox* phiz   = fphiz.fabPtr();
    FArrayBox* phiz_x = fphiz_x.fabPtr();
    FArrayBox* phiz_y = fphiz_y.fabPtr();
    FArrayBox* slope  = fslope.fabPtr();
*/

/* CHECK CFL. ADD MODIFIED MAX LOOP HERE. 
    Real umax = ;
    Real vmax = ;
    Real wmax = ;

    umax = maxval(abs(xvel))
    vmax = maxval(abs(yvel))
    wmax = maxval(abs(zvel))

    if ( umax*dt .ge. dx(1) .or. &
         vmax*dt .ge. dx(2) .or. &
         wmax*dt .ge. dx(3) ) then
       print *, "umax = ", umax, ", vmax = ", vmax, ", wmax = ", wmax, ", dt = ", dt, ", dx = ", dx
       call bl_error("CFL violation. Use smaller adv.cfl.")
    end if
*/

//  call a function to compute flux

    compute_flux_3d(bx, dtdx,
                    statein,
                    AMREX_D_DECL(xvel, yvel, zvel),
                    AMREX_D_DECL(fx, fy, fz));
  
/*
    compute_flux_3d(lo, hi, dt, dx, 
                    uin, ui_lo, ui_hi, 
                    vx, vx_lo, vx_hi, 
                    vy, vy_lo, vy_hi,
                    vz, vz_lo, vz_h,i
                    fx, fx_lo, fx_hi,
                    fy, fy_lo, fy_hi,
                    fz, fz_lo, fz_hi,
                    phix, phix_y, phix_z, 
                    phiy, phiy_x, phiy_z, 
                    phiz, phiz_x, phiz_y, 
                    slope, glo, ghi);
*/

/*  NOT USED HERE. MOVE TO COMPUTE_FLUX_3D
    fphix.clear();
    fphix_y.clear();
    fphix_z.clear();
    fphiy.clear();
    fphiy_x.clear();
    fphiy_z.clear();
    fphiz.clear();
    fphiz_x.clear();
    fphiz_y.clear();
    fslope.clear();
*/

    // Do a conservative update
    AMREX_LAUNCH_DEVICE_LAMBDA(bx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);

        const auto uin  = statein.view(lo);
        const auto uout = stateout.view(lo);
        const auto flxx = fx.view(lo);
        const auto flxy = fy.view(lo);
        const auto flxz = fz.view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    uout(i,j,k) = uin(i,j,k) + 
                        ( (flxx(i,j,k) - flxx(i+1,j,k)) * dtdx[0] 
                        + (flxy(i,j,k) - flxy(i,j+1,k)) * dtdx[1] 
                        + (flxz(i,j,k) - flxz(i,j,k+1)) * dtdx[2] );
                }
            }
        }
    });

    // Scale by face area in order to correctly reflux
    const Box xbx = amrex::growHi(bx, 0, 1); 
    const Box ybx = amrex::growHi(bx, 1, 1); 
    const Box zbx = amrex::growHi(bx, 2, 1); 

    AMREX_LAUNCH_DEVICE_LAMBDA(xbx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);
        const auto flxx = fx.view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    flxx(i,j,k) = flxx(i,j,k) * (dt * dx[1]*dx[2]);
                }
            }
        }
    });

    AMREX_LAUNCH_DEVICE_LAMBDA(ybx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);
        const auto flxy = fy.view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    flxy(i,j,k) = flxy(i,j,k) * (dt * dx[0]*dx[2]);
                }
            }
        }
    });

    AMREX_LAUNCH_DEVICE_LAMBDA(zbx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);
        const auto flxz = fz.view(lo);

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                for (int i = 0; i < len.x; ++i) {
                    flxz(i,j,k) = flxz(i,j,k) * (dt * dx[0]*dx[1]);
                }
            }
        }
    });

/*
    POSSIBLE ALTERNATIVE TO THREE LAUNCHES AND THREE BOX BUILDS?
    SIMPLE EXAMPLE. MAY BE MORE USEFUL FOR COMPUTE_FLUX_3D AND 
    ITS 12 SEPARATE LOOPS WITH UNIQUE LOOP RANGES. 

    WORTH IT?: 
    1) READABLE? (MIRRORED REFLECTION OF PREVIOUS METHOD).
    2) LAUNCH LATENCY VS. NON-COALESCED MEMORY ACCESSES (ACROSS LOOPS).

    CAN (SOMEHOW) THESE LOOPS BE COMBINED TO BE CLEANER?
    (E.G. A "SHARED" LOOP AND A "REMAINDER" LOOP?
     IFS PLACED ON OUTER LOOP LEVELS TO PREVENT INTERFERENCE WITH SIMD?)

    // One box with hi(i)+1 for all i.
    const Box abx = amrex::growHi(bx, 0); 
                       abx.growHi(1);
                       abx.growHi(2);

    AMREX_LAUNCH_DEVICE_LAMBDA(abx, tbx
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);

        // Alter lo to start at the appropriate place.

        const auto flxx = fx.view(lo);
        const auto flxy = fy.view(lo);
        const auto flxz = fz.view(lo);

        for         (int k = 0; k < (len.z-1); ++k) {
            for     (int j = 0; j < (len.y-1); ++j) {
                for (int i = 0; i < len.x;     ++i) {
                    flxx(i,j,k) = flxx(i,j,k) * (dt * dx[1]*dx[2]);
                }
            }
        }

        for         (int k = 0; k < (len.z-1); ++k) {
            for     (int j = 0; j < len.y;     ++j) {
                for (int i = 0; i < (len.x-1); ++i) {
                    flxy(i,j,k) = flxy(i,j,k) * (dt * dx[0]*dx[2]);
                }
            }
        }

        for         (int k = 0; k < len.z;     ++k) {
            for     (int j = 0; j < (len.y-1); ++j) {
                for (int i = 0; i < (len.x-1); ++i) {
                    flxz(i,j,k) = flxz(i,j,k) * (dt * dx[0]*dx[1]);
                }
            }
        }

    }
*/

}

#endif
