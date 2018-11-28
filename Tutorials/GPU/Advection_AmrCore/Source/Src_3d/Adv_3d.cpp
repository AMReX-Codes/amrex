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

/*
// CHECK CFL. ADD MODIFIED MAX LOOP HERE. 
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
}

#endif
