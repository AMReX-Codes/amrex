#ifndef face_velocity_cpp
#define face_velocity_cpp

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

using namespace amrex;

void get_face_velocity(const int level, const Real time,
                       AMREX_D_DECL(FArrayBox& xvel,
                                    FArrayBox& yvel,
                                    FArrayBox& zvel),
                       GeometryData const& geom)
{

// Create FArrayBox for psi. Base it on AmrCoreAdv.cpp call of get_face_velocity.

    Box xbx = xvel.box();
    Box ybx = yvel.box();
    Box psi_box = Box(IntVect(AMREX_D_DECL(std::min(xbx.smallEnd(0)-1, ybx.smallEnd(0)-1), 
                                           std::min(xbx.smallEnd(1)-1, ybx.smallEnd(0)-1), 0)),
                      IntVect(AMREX_D_DECL(std::max(xbx.bigEnd(0)  ,   ybx.bigEnd(0)+1),
                                           std::max(xbx.bigEnd(1)+1,   ybx.bigEnd(1)  ), 0)));

/*  Setup for size of psi FArrayBox in Fortran.
    plo(1) = min(vx_l1-1, vy_l1-1)
    plo(2) = min(vx_l2-1, vy_l2-1)
    phi(1) = max(vx_h1  , vy_h1+1)
    phi(2) = max(vx_h2+1, vy_h2  )
*/

    Gpu::DeviceFab psi_dfab(psi_box, 1);
    FArrayBox* psifab = psi_dfab.fabPtr();

    // Calculate psi
    AMREX_LAUNCH_DEVICE_LAMBDA(psi_box, tbx,
    {
//        Use C++ M_PI instead.
//        const Real M_PI = 3.141592653589793238462643383279502884197;

        const auto len = length(tbx);
        const auto lo  = lbound(tbx);
        const auto psi = psifab->view(lo); 
        const Real* AMREX_RESTRICT prob_lo = geom.ProbLo();
        const Real* AMREX_RESTRICT dx      = geom.CellSize(); 

        for     (int j = 0; j < len.y; ++j) {
            Real y = dx[1]*(0.5 + j) + prob_lo[1]; 
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                Real x = dx[0]*(0.5 + i) + prob_lo[0];
                psi(i,j,0) = pow(sin(M_PI*x), 2) * pow(sin(M_PI*y), 2)
                           * cos(M_PI*time/2.0) * 1.0/M_PI; 
                }
            }
    });

    // Change box and calculate x velocity.
    AMREX_LAUNCH_DEVICE_LAMBDA(xbx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);
        const auto vx  = xvel.view(lo);

        const auto psi_lo = lbound(tbx);
        const auto psi = psifab->view(psi_lo);

        const Real* AMREX_RESTRICT prob_lo = geom.ProbLo();
        const Real* AMREX_RESTRICT dx      = geom.CellSize(); 

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                Real y = dx[1] * (0.5 + j) + prob_lo[1];
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    Real x = dx[0] * i + prob_lo[0];
                    vx(i,j,k) = -( (psi(i,j+1,0)+psi(i-1,j+1,0)) - (psi(i,j-1,0)+psi(i-1,j-1,0)) ) * (0.25/dx[1]); 
                }
            }
        }
    });

    // Change box and calculate y velocity.
    AMREX_LAUNCH_DEVICE_LAMBDA(ybx, tbx,
    {
        const auto len = length(tbx);
        const auto lo  = lbound(tbx);
        const auto vy  = yvel.view(lo); 

        const auto psi_lo = lbound(tbx);
        const auto psi = psifab->view(psi_lo);

        const Real* AMREX_RESTRICT prob_lo = geom.ProbLo();
        const Real* AMREX_RESTRICT dx      = geom.CellSize(); 

        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                Real y = dx[1] * j + prob_lo[1];
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    Real x = dx[0] * (0.5 + i) + prob_lo[0];
                    vy(i,j,k) = -( (psi(i+1,j,0)+psi(i+1,j-1,0)) - (psi(i-1,j,0)+psi(i-1,j-1,0)) ) * (0.25/dx[1]); 
                }
            }
        }
    });

    phi_dfab.clear();

    // If 3d, set z velocity to 1.
#if (AMREX_SPACEDIM == 3)
    zvel.setVal(1.0);
#endif

}

#endif
