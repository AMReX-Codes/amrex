#ifndef INIT_PROB_K_H_
#define INIT_PROB_K_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void actual_init_prob (int i, int j, int k,
                       amrex::GpuArray<amrex::Array4<amrex::Real>,3> const& rhs,
                       amrex::GpuArray<amrex::Array4<amrex::Real>,3> const& sol,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& problo,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
                       amrex::Real alpha, amrex::Real beta)
{
    using namespace amrex;

    constexpr Real pi = amrex::Math::pi<Real>();

    Real xnd = problo[0] + Real(i)*dx[0];
    Real ynd = problo[1] + Real(j)*dx[1];
    Real xcc = xnd + Real(0.5)*dx[0];
    Real ycc = ynd + Real(0.5)*dx[1];
#if (AMREX_SPACEDIM == 3)
    Real znd = problo[2] + Real(k)*dx[2];
    Real zcc = znd + Real(0.5)*dx[2];
#endif

    if (sol[0].contains(i,j,k)) {
        Real x = xcc;
        Real y = ynd;
        Real Ex = std::sin(pi*x) * std::sin(Real(2.5)*pi*y);
#if (AMREX_SPACEDIM == 3)
        Real z = znd;
        Ex *= std::sin(Real(2.0)*pi*z + Real(1./3.)*pi);
#endif
        sol[0](i,j,k) = Ex;
    }

    if (sol[1].contains(i,j,k)) {
        Real x = xnd;
        Real y = ycc;
        Real Ey = std::cos(Real(2.5)*pi*x) * std::sin(Real(3.)*pi*y);
#if (AMREX_SPACEDIM == 3)
        Real z = znd;
        Ey *= std::sin(Real(4.)*pi*z + Real(0.25)*pi);
#endif
        sol[1](i,j,k) = Ey;
    }

    if (sol[2].contains(i,j,k)) {
        Real x = xnd;
        Real y = ynd;
        Real Ez = std::cos(Real(3.5)*pi*x) * std::sin(Real(3.5)*pi*y);
#if (AMREX_SPACEDIM == 3)
        Real z = zcc;
        Ez *= std::sin(Real(4.)*pi*z + Real(1./6.)*pi);
#endif
        sol[2](i,j,k) = Ez;
    }

    if (rhs[0].contains(i,j,k)) {
        Real x = xcc;
        Real y = ynd;
#if (AMREX_SPACEDIM == 2)
        Real cce = Real(-7.5)*pi*pi*std::sin(Real(2.5)*pi*x)*std::cos(Real(3.)*pi*y)
            + Real(6.25)*pi*pi*std::sin(pi*x)*std::sin(Real(2.5)*pi*y);
#else
        Real z = znd;
        Real cce = Real(-7.5)*pi*pi*std::sin(Real(2.5)*pi*x)*std::cos(Real(3.)*pi*y)*std::sin(Real(4.)*pi*z+Real(0.25)*pi)
            + Real(6.25)*pi*pi*std::sin(pi*x)*std::sin(Real(2.5)*pi*y)*std::sin(Real(2.)*pi*z+Real(1./3.)*pi)
            - Real(14.)*pi*pi*std::sin(Real(3.5)*pi*x)*std::sin(Real(3.5)*pi*y)*std::cos(Real(4.)*pi*z+Real(1./6.)*pi)
            + Real(4.)*pi*pi*std::sin(pi*x)*std::sin(Real(2.5)*pi*y)*std::sin(Real(2.)*pi*z+Real(1./3.)*pi);
#endif
        rhs[0](i,j,k) = alpha*cce + beta*sol[0](i,j,k);
    }

    if (rhs[1].contains(i,j,k)) {
        Real x = xnd;
        Real y = ycc;
#if (AMREX_SPACEDIM == 2)
        Real cce = Real(6.25)*pi*pi*std::cos(Real(2.5)*pi*x)*std::sin(Real(3.)*pi*y)
            + Real(2.5)*pi*pi*std::cos(pi*x)*std::cos(Real(2.5)*pi*y);
#else
        Real z = znd;
        Real cce = Real(6.25)*pi*pi*std::cos(Real(2.5)*pi*x)*std::sin(Real(3.)*pi*y)*std::sin(Real(4.)*pi*z+Real(0.25)*pi)
            + Real(2.5)*pi*pi*std::cos(pi*x)*std::cos(Real(2.5)*pi*y)*std::sin(Real(2.)*pi*z+Real(1./3.)*pi)
            + Real(14.)*pi*pi*std::cos(Real(3.5)*pi*x)*std::cos(Real(3.5)*pi*y)*std::cos(Real(4.)*pi*z+Real(1./6.)*pi)
            + Real(16.)*pi*pi*std::cos(Real(2.5)*pi*x)*std::sin(Real(3.)*pi*y)*std::sin(Real(4.)*pi*z+Real(0.25)*pi);
#endif
        rhs[1](i,j,k) = alpha*cce + beta*sol[1](i,j,k);
    }

    if (rhs[2].contains(i,j,k)) {
        Real x = xnd;
        Real y = ynd;
#if (AMREX_SPACEDIM == 2)
        Real cce = Real(24.5)*pi*pi*std::cos(Real(3.5)*pi*x)*std::sin(Real(3.5)*pi*y);
#else
        Real z = zcc;
        Real cce = Real(24.5)*pi*pi*std::cos(Real(3.5)*pi*x)*std::sin(Real(3.5)*pi*y)*std::sin(Real(4.)*pi*z+Real(1./6.)*pi)
            + Real(2.)*pi*pi*std::cos(pi*x)*std::sin(Real(2.5)*pi*y)*std::cos(Real(2.)*pi*z+Real(1./3.)*pi)
            + Real(12.)*pi*pi*std::cos(Real(2.5)*pi*x)*std::cos(Real(3.)*pi*y)*std::cos(Real(4.)*pi*z+Real(0.25)*pi);
#endif
        rhs[2](i,j,k) = alpha*cce + beta*sol[2](i,j,k);
    }
}

#endif
