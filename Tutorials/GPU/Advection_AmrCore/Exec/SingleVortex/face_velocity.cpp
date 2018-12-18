
#include <face_velocity.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AmrCoreAdv_F.H>

using namespace amrex;

AMREX_GPU_DEVICE
void get_face_velocity_psi(Box const& bx,
                           const Real time,
                           FArrayBox& psifab,
                           GeometryData const& geom)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto psi = psifab.view(lo); 
    const Real* AMREX_RESTRICT prob_lo = geom.ProbLo();
    const Real* AMREX_RESTRICT dx      = geom.CellSize(); 

    for     (int j = 0; j < len.y; ++j) {
        Real y = dx[1]*(0.5 + j+lo.y) + prob_lo[1]; 
        AMREX_PRAGMA_SIMD
        for (int i = 0; i < len.x; ++i) {
            Real x = dx[0]*(0.5 + i+lo.x) + prob_lo[0];
            psi(i,j,0) = pow(sin(M_PI*x), 2) * pow(sin(M_PI*y), 2)
                       * cos(M_PI*time/2.0) * 1.0/M_PI; 
        }
    }
}

AMREX_GPU_DEVICE
void get_face_velocity_x(Box const& bx,
                         FArrayBox& xvelfab,
                         const FArrayBox& psifab,
                         GeometryData const& geom)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto vx  = xvelfab.view(lo);
    const auto psi = psifab.view(IntVect{lo.x, lo.y, 0});

    const Real* AMREX_RESTRICT prob_lo = geom.ProbLo();
    const Real* AMREX_RESTRICT dx      = geom.CellSize(); 

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                vx(i,j,k) = -( (psi(i,j+1,0)+psi(i-1,j+1,0)) - (psi(i,j-1,0)+psi(i-1,j-1,0)) ) * (0.25/dx[1]); 
            }
        }
    }
}

AMREX_GPU_DEVICE
void get_face_velocity_y(Box const& bx,
                         FArrayBox& yvelfab,
                         const FArrayBox& psifab,
                         GeometryData const& geom)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto vy  = yvelfab.view(lo); 
    const auto psi = psifab.view(IntVect{lo.x, lo.y, 0});

    const Real* AMREX_RESTRICT prob_lo = geom.ProbLo();
    const Real* AMREX_RESTRICT dx      = geom.CellSize(); 

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = 0; i < len.x; ++i) {
                vy(i,j,k) =  ( (psi(i+1,j,0)+psi(i+1,j-1,0)) - (psi(i-1,j,0)+psi(i-1,j-1,0)) ) * (0.25/dx[0]); 
            }
        }
    }
}


AMREX_GPU_DEVICE
void get_face_velocity_z(Box const& bx,
                         FArrayBox& zvelfab,
                         const FArrayBox& psifab,
                         GeometryData const& geom)
{
    // Set z velocity to 1.
    zvelfab.setVal(1.0,bx,0);
}

