#ifndef FACE_VELOCITY_3D_H_
#define FACE_VELOCITY_3D_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Array.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_psi(amrex::Box const& bx,
                           const amrex::Real time,
                           amrex::Array4<amrex::Real> const& psi,
                           amrex::Real const dx,
                           amrex::Real const dy,
                           amrex::Real const /*dz*/,
                           amrex::Real const prob_lo_x,
                           amrex::Real const prob_lo_y,
                           amrex::Real const /*prob_lo_z*/)
{
    using namespace amrex;
    constexpr Real PI = 3.141592653589793238462643383279502884197;

    // Extract box bounds
    const auto plo = lbound(bx);
    const auto phi = ubound(bx);

    // Compute streamfunction psi
    for     (int j = plo.y; j <= phi.y; ++j)
    {
        Real y = ( (Real)j + 0.5 )*dy + prob_lo_y;
        AMREX_PRAGMA_SIMD
        for (int i = plo.x; i <= phi.x; ++i)
        {
            Real x = ( (Real)i + 0.5 )*dx + prob_lo_x;
            psi(i,j,0) = std::pow(std::sin(PI*x), 2) * std::pow(std::sin(PI*y), 2)
                       * std::cos(PI*time/2.0) * 1.0/PI;
        }
    }

    return;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_x(int i, int j, int k,
                         amrex::Array4<amrex::Real>       const& vx,
                         amrex::Array4<const amrex::Real> const& psi,
                         const amrex::Real dy)
{
    vx(i,j,k) = -( (psi(i,j+1,0)+psi(i-1,j+1,0)) - (psi(i,j-1,0)+psi(i-1,j-1,0)) ) * (0.25/dy);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_y(int i, int j, int k,
                         amrex::Array4<amrex::Real>       const& vy,
                         amrex::Array4<const amrex::Real> const& psi,
                         const amrex::Real dx)
{
    vy(i,j,k) =  ( (psi(i+1,j,0)+psi(i+1,j-1,0)) - (psi(i-1,j,0)+psi(i-1,j-1,0)) ) * (0.25/dx);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void get_face_velocity_z(int i, int j, int k,
                         amrex::Array4<amrex::Real>       const& vz)
{
    vz(i,j,k) = 1.0;
}

AMREX_GPU_HOST
AMREX_FORCE_INLINE
void get_face_velocity(const amrex::Real time,
                       amrex::FArrayBox& vx,
                       amrex::FArrayBox& vy,
                       amrex::FArrayBox& vz,
                       amrex::GpuArray<amrex::Real,3> dx,
                       amrex::GpuArray<amrex::Real,3> prob_lo)
{
    using namespace amrex;

    // Extract boxes out of vx and vy fabs
    GpuArray<Box, 2> nbx;
    nbx[0] = vx.box();
    nbx[1] = vy.box();
    const Box& ngbxx = nbx[0];
    const Box& ngbxy = nbx[1];

    // Construct GpuArray for velocities
    GpuArray<Array4<Real>, 3> vel{ vx.array(), vy.array(), vz.array() };

    // Construct box for streamfunction psi
    const Box& psibox = Box(IntVect(AMREX_D_DECL(std::min( ngbxx.smallEnd(0)-1, ngbxy.smallEnd(0)-1 ),
                                                 std::min( ngbxx.smallEnd(1)-1, ngbxy.smallEnd(1)-1 ),
                                                 0)),
                            IntVect(AMREX_D_DECL(std::max( ngbxx.bigEnd(0),   ngbxy.bigEnd(0)+1 ),
                                                 std::max( ngbxx.bigEnd(1)+1, ngbxy.bigEnd(1)   ),
                                                 0)));

    // Construct fab for streamfunction psi and extract pointer to array
    FArrayBox psifab( psibox, 1 );
    Elixir psieli = psifab.elixir();
    Array4<Real> psi = psifab.array();

    // Compute streamfunction psi (can run on GPU)
    amrex::launch(psibox,
    [=] AMREX_GPU_DEVICE (const Box& tbx)
    {
        get_face_velocity_psi(tbx, time, psi,
                              dx[0], dx[1], dx[2],
                              prob_lo[0], prob_lo[1], prob_lo[2]);
    });

    // Compute face velocities
    amrex::ParallelFor
        (vx.box(), vy.box(), vz.box(),
         [=] AMREX_GPU_DEVICE (int i, int j, int k)
         {
             get_face_velocity_x(i, j, k, vel[0], psi, dx[1]);
         },
         [=] AMREX_GPU_DEVICE (int i, int j, int k)
         {
             get_face_velocity_y(i, j, k, vel[1], psi, dx[0]);
         },
         [=] AMREX_GPU_DEVICE (int i, int j, int k)
         {
             get_face_velocity_z(i, j, k, vel[2]);
         });

    return;
}

#endif
