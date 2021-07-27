#ifndef FACE_VELOCITY_3D_H_
#define FACE_VELOCITY_3D_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

void get_face_velocity(const amrex::Real time,
                       amrex::FArrayBox const& vx,
                       amrex::FArrayBox const& vy,
                       amrex::FArrayBox const& vz,
                       const amrex::Real* dx, const amrex::Real* prob_lo)
{
    using namespace amrex;

    // Extract bounds from vx and vy fabs
    Dim3 xlo = lbound(vx.box());
    Dim3 xhi = ubound(vx.box());
    Dim3 ylo = lbound(vy.box());
    Dim3 yhi = ubound(vy.box());

    // Define box bounds for streamfunction psi
    Dim3 plo;
    plo.x = min( xlo.x-1, ylo.x-1 );
    plo.y = min( xlo.y-1, ylo.y-1 );
    Dim3 phi;
    phi.x = max( xhi.x,   yhi.x+1 );
    phi.y = max( xhi.y+1, yhi.y   );

    // Construct box for streamfunction psi
    const Box& psibox = Box( plo, phi );

    // Construct fab for streamfunction psi and extract pointer to array
    // TODO: Convert to GPUArray
    FArrayBox psifab( psibox, 1 );

    Array4<Real>& psi = psifab.array();

    // TODO: use amrex::launch to perform on GPU
    get_face_velocity_psi(psibox, time,
                          psi,
                          dx, prob_lo);

    // TODO: convert into amrex::ParallelFor
    //
    // x velocity
    for         (int k = xlo.z; j <= xhi.z; ++k)
        for     (int j = xlo.y; j <= xhi.y; ++j)
            for (int i = xlo.x; i <= xhi.x; ++i)
                get_face_velocity_x(i, j, k,
                                    vx, psi,
                                    dx[1]);
    // y velocity
    for         (int k = ylo.z; j <= yhi.z; ++k)
        for     (int j = ylo.y; j <= yhi.y; ++j)
            for (int i = ylo.x; i <= yhi.x; ++i)
                get_face_velocity_y(i, j, k,
                                    vy, psi,
                                    dx[0]);

    // z velocity
    vz.setVal(1.0);

    // call bl_deallocate(psi)

    return;
}

// TODO: Convert to GPU kernel
void get_face_velocity_psi(amrex::Box const& bx,
                           const amrex::Real time,
                           amrex::Array4<amrex::Real> const& psi,
                           amrex::Real const& dx, amrex::Real const& prob_lo)
{
    using namespace amrex;
    constexpr Real PI = 3.141592653589793238462643383279502884197;

    // Extract box bounds
    const auto plo = lbound(bx);
    const auto phi = ubound(bx);

    // Compute streamfunction psi
    for     (int j = plo.y; j <= phi.y; ++j)
    {
        Real y = ( (Real)j + 0.5 )*dx[1] + prob_lo[1];
        for (int i = plo.x; i <= phi.x; ++i)
        {
            Real x = ( (Real)i + 0.5 )*dx[0] + prob_lo[0];
            psi(i,j,0) = std::pow(std::sin(PI*x), 2) * std::pow(std::sin(PI*y), 2)
                       * std::cos(PI*time/2.0) * 1.0/PI;
        }
    }

    return;
}

void get_face_velocity_x(int i, int j, int k,
                         amrex::Array4<amrex::Real>       const& vx,
                         amrex::Array4<const amrex::Real> const& psi,
                         const amrex::Real dy)
{
    vx(i,j,k) = -( (psi(i,j+1,0)+psi(i-1,j+1,0)) - (psi(i,j-1,0)+psi(i-1,j-1,0)) ) * (0.25/dy);
}

void get_face_velocity_y(int i, int j, int k,
                         amrex::Array4<amrex::Real>       const& vx,
                         amrex::Array4<const amrex::Real> const& psi,
                         amrex::FArrayBox& vy,
                         const amrex::Real dx)
{
    vy(i,j,k) =  ( (psi(i+1,j,0)+psi(i+1,j-1,0)) - (psi(i-1,j,0)+psi(i-1,j-1,0)) ) * (0.25/dx);
}

#endif
