#ifndef GetFaceVelocity_cpp
#define GetFaceVelocity_cpp

#include <AMReX_FArrayBox.H>

void get_face_velocity(const int level, const Real time,
                       AMREX_D_DECL(FArrayBox& xvel,
                                    FArrayBox& yvel,
                                    FArrayBox& zvel),
                       GeometryData const& geom)
{

// Create FArrayBox for psi. Base it on AmrCoreAdv.cpp call of get_face_velocity.

    FArrayBox psi;
    Box xbx = xvel.box();
    Box ybx = yvel.box();
    Box psi_box = ({min(xbx.lo(0)-1, ybx.lo(0)-1), min(xbx.lo(1)-1, ybx.lo(0)-1), 0}
                   {max(xbx.hi(0)  , ybx.hi(0)+1), max(xbx.hi(1)+1, ybx.hi(1)   , 0});

    psi.resize(psi_box);

/*  Setup for size of psi FArrayBox in Fortran.
    plo(1) = min(vx_l1-1, vy_l1-1)
    plo(2) = min(vx_l2-1, vy_l2-1)
    phi(1) = max(vx_h1  , vy_h1+1)
    phi(2) = max(vx_h2+1, vy_h2  )
*/

    // Calculate psi
    const Box& box = psi_box;
    AMREX_LAUNCH_DEVICE_LAMBDA(box, tbox,
    {
        get_psi(time, tbox, psi, geom);
    });

    // Change box and calculate x velocity.
    box = xbx; 
    AMREX_LAUNCH_DEVICE_LAMBDA(box, tbox,
    {
        get_xvel(time, tbox, xvel, psi, geom);
    });

    // Change box and calculate y velocity.
    box = ybx;
    AMREX_LAUNCH_DEVICE_LAMBDA(box, tbox,
    {
        get_yvel(time, tbox, yvel, psi, geom);
    });

#ifdef (AMREX_SPACEDIM == 3)
    zvel.setVal(1.0);
#endif

}

void get_psi(const Real time, Box const& bx, FArrayBox& psi, GeometryData const& geom)
{
    const Real M_PI = 3.141592653589793238462643383279502884197;

    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto psi = psifab.view(lo); 
    const Real* AMREX_RESTRICT prob_lo = geom.ProbLo();
    const Real* AMREX_RESTRICT dx      = geom.CellSize(); 

    for     (int j = 0; j < len.y; ++j) {
        Real y = dx[1]*(0.5 + j) + prob_lo[1]; 
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            Real x = dx[0]*(0.5 + i) + prob_lo[0];
            psi(i,j,0) = pow(sin(M_PI*x), 2) * pow(sin(M_PI*y), 2)
                       * cos(M_PI*time/2.0) * 1/M_PI;
            }
        }
}

void get_xvel(const Real time, Box const& bx, FArrayBox& xvel,
              FArrayBox& psi, GeometryData const& geom)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto vx  = xvel.view(lo); 
    const Real* AMREX_RESTRICT prob_lo = geom.ProbLo();
    const Real* AMREX_RESTRICT dx      = geom.CellSize(); 

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            Real y = dx[1] * (0.5 + j) + prob_lo[1];
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                Real x = dx[0] * i + prob_lo[0];
                vx(i,j,k) = -( (psi(i,j+1)+psi(i-1,j+1)) - (psi(i,j-1)+psi(i-1,j-1)) ) * (0.25/dx[1]); 
            }
        }
    }
}

void get_yvel(const Real time, Box const& bx, FArrayBox& yvel,
              FArrayBox& psi, GeometryData const& geom)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto vy  = yvel.view(lo); 
    const Real* AMREX_RESTRICT prob_lo = geom.ProbLo();
    const Real* AMREX_RESTRICT dx      = geom.CellSize(); 

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            Real y = dx[1] * j + prob_lo[1];
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                Real x = dx[0] * (0.5 + i) + prob_lo[0];
                vy(i,j,k) = -( (psi(i+1,j)+psi(i+1,j-1)) - (psi(i-1,j)+psi(i-1,j-1)) ) * (0.25/dx[1]); 
            }
        }
    }
}

#endif
