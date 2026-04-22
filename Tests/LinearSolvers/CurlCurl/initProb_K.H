#ifndef INIT_PROB_K_H_
#define INIT_PROB_K_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void actual_init_prob (int i, int j, int k,
                       amrex::GpuArray<amrex::Array4<amrex::Real>,3> const& rhs,
                       amrex::GpuArray<amrex::Array4<amrex::Real>,3> const& sol,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& problo,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
                       amrex::Real alpha, amrex::Real beta,
                       amrex::Array4<amrex::Real> const& alpha_node)
{
    using namespace amrex;

    constexpr Real pi = amrex::Math::pi<Real>();

    Real xnd = problo[0] + Real(i)*dx[0];
    Real xcc = xnd + Real(0.5)*dx[0];
#if (AMREX_SPACEDIM > 1)
    Real ynd = problo[1] + Real(j)*dx[1];
    Real ycc = ynd + Real(0.5)*dx[1];
#endif
#if (AMREX_SPACEDIM == 3)
    Real znd = problo[2] + Real(k)*dx[2];
    Real zcc = znd + Real(0.5)*dx[2];
#endif

    auto get_alpha = [&] (AMREX_D_DECL(Real x, Real y, Real z))
    {
        Real a = AMREX_D_TERM(std::cos(Real(2.0)*pi*x),
                            * std::cos(Real(3.0)*pi*y),
                            * std::cos(Real(4.0)*pi*z));
        return Real(1.0) + Real(0.3)*a;
    };

    if (alpha_node && alpha_node.contains(i,j,k)) {
        alpha_node(i,j,k) = get_alpha(AMREX_D_DECL(xnd,ynd,znd));
    }

    if (sol[0].contains(i,j,k)) {
        Real x = xcc;
        Real Ex = std::sin(pi*x);
#if (AMREX_SPACEDIM > 1)
        Real y = ynd;
        Ex *= std::sin(Real(2.5)*pi*y);
#endif
#if (AMREX_SPACEDIM == 3)
        Real z = znd;
        Ex *= std::sin(Real(2.0)*pi*z + Real(1./3.)*pi);
#endif
        sol[0](i,j,k) = Ex;
    }

    if (sol[1].contains(i,j,k)) {
        Real x = xnd;
        Real Ey = std::cos(Real(2.5)*pi*x);
#if (AMREX_SPACEDIM > 1)
        Real y = ycc;
        Ey *= std::sin(Real(3.)*pi*y);
#endif
#if (AMREX_SPACEDIM == 3)
        Real z = znd;
        Ey *= std::sin(Real(4.)*pi*z + Real(0.25)*pi);
#endif
        sol[1](i,j,k) = Ey;
    }

    if (sol[2].contains(i,j,k)) {
        Real x = xnd;
        Real Ez = std::cos(Real(3.5)*pi*x);
#if (AMREX_SPACEDIM > 1)
        Real y = ynd;
        Ez *= std::sin(Real(3.5)*pi*y);
#endif
#if (AMREX_SPACEDIM == 3)
        Real z = zcc;
        Ez *= std::sin(Real(4.)*pi*z + Real(1./6.)*pi);
#endif
        sol[2](i,j,k) = Ez;
    }

    if (rhs[0].contains(i,j,k)) {
        Real x = xcc;
#if (AMREX_SPACEDIM > 1)
        Real y = ynd;
#endif
#if (AMREX_SPACEDIM == 1)
        Real cce = 0;
#elif (AMREX_SPACEDIM == 2)
        Real cce = Real(-7.5)*pi*pi*std::sin(Real(2.5)*pi*x)*std::cos(Real(3.)*pi*y)
            + Real(6.25)*pi*pi*std::sin(pi*x)*std::sin(Real(2.5)*pi*y);
#else
        Real z = znd;
        Real cce = Real(-7.5)*pi*pi*std::sin(Real(2.5)*pi*x)*std::cos(Real(3.)*pi*y)*std::sin(Real(4.)*pi*z+Real(0.25)*pi)
            + Real(6.25)*pi*pi*std::sin(pi*x)*std::sin(Real(2.5)*pi*y)*std::sin(Real(2.)*pi*z+Real(1./3.)*pi)
            - Real(14.)*pi*pi*std::sin(Real(3.5)*pi*x)*std::sin(Real(3.5)*pi*y)*std::cos(Real(4.)*pi*z+Real(1./6.)*pi)
            + Real(4.)*pi*pi*std::sin(pi*x)*std::sin(Real(2.5)*pi*y)*std::sin(Real(2.)*pi*z+Real(1./3.)*pi);
#endif

        Real a = (alpha_node) ? get_alpha(AMREX_D_DECL(x,y,z)) : alpha;
        rhs[0](i,j,k) = a*cce + beta*sol[0](i,j,k);

#if (AMREX_SPACEDIM >= 2)
        if (alpha_node) {
#if (AMREX_SPACEDIM == 2)
            Real rhs2 = Real(-3.)*pi*std::cos(Real(2.)*pi*x)*std::sin(Real(3.)*pi*y) *
                (Real(-2.5)*pi*std::sin(Real(2.5)*pi*x)*std::sin(Real(3.)*pi*y)
                 + Real(-2.5)*pi*std::sin(pi*x)*std::cos(Real(2.5)*pi*y));
#elif (AMREX_SPACEDIM == 3)
            Real rhs2 = Real(-3.)*pi*std::cos(Real(2)*pi*x)*std::sin(Real(3)*pi*y)*std::cos(Real(4)*pi*z) *
                (Real(-2.5)*pi*std::sin(Real(2.5)*pi*x)*std::sin(Real(3)*pi*y)*std::sin(Real(4)*pi*z+Real(0.25)*pi) +
                 Real(-2.5)*pi*std::sin(pi*x)*std::cos(Real(2.5)*pi*y)*std::sin(Real(2)*pi*z+Real(1./3.)*pi))
                + Real(4.)*pi*std::cos(Real(2)*pi*x)*std::cos(Real(3)*pi*y)*std::sin(Real(4)*pi*z) *
                (Real(2)*pi*std::sin(pi*x)*std::sin(Real(2.5)*pi*y)*std::cos(Real(2)*pi*z+Real(1./3.)*pi) +
                 Real(3.5)*pi*std::sin(Real(3.5)*pi*x)*std::sin(Real(3.5)*pi*y)*std::sin(Real(4)*pi*z+Real(1./6.)*pi));
#endif
            rhs[0](i,j,k) += Real(0.3)*rhs2;
        }
#endif
    }

    if (rhs[1].contains(i,j,k)) {
        Real x = xnd;
#if (AMREX_SPACEDIM > 1)
        Real y = ycc;
#endif
#if (AMREX_SPACEDIM == 1)
        Real cce = Real(6.25)*pi*pi*std::cos(Real(2.5)*pi*x);
#elif (AMREX_SPACEDIM == 2)
        Real cce = Real(6.25)*pi*pi*std::cos(Real(2.5)*pi*x)*std::sin(Real(3.)*pi*y)
            + Real(2.5)*pi*pi*std::cos(pi*x)*std::cos(Real(2.5)*pi*y);
#else
        Real z = znd;
        Real cce = Real(6.25)*pi*pi*std::cos(Real(2.5)*pi*x)*std::sin(Real(3.)*pi*y)*std::sin(Real(4.)*pi*z+Real(0.25)*pi)
            + Real(2.5)*pi*pi*std::cos(pi*x)*std::cos(Real(2.5)*pi*y)*std::sin(Real(2.)*pi*z+Real(1./3.)*pi)
            + Real(14.)*pi*pi*std::cos(Real(3.5)*pi*x)*std::cos(Real(3.5)*pi*y)*std::cos(Real(4.)*pi*z+Real(1./6.)*pi)
            + Real(16.)*pi*pi*std::cos(Real(2.5)*pi*x)*std::sin(Real(3.)*pi*y)*std::sin(Real(4.)*pi*z+Real(0.25)*pi);
#endif

        Real a = (alpha_node) ? get_alpha(AMREX_D_DECL(x,y,z)) : alpha;
        rhs[1](i,j,k) = a*cce + beta*sol[1](i,j,k);

        if (alpha_node) {
#if (AMREX_SPACEDIM == 1)
            Real rhs2 = Real(-5.0)*pi*pi* std::sin(Real(2.0)*pi*x)* std::sin(Real(2.5)*pi*x);
#elif (AMREX_SPACEDIM == 2)
            Real rhs2 = Real(2)*pi*std::sin(Real(2.)*pi*x)*std::cos(Real(3.)*pi*y) *
                (Real(-2.5)*pi*std::sin(Real(2.5)*pi*x)*std::sin(Real(3.)*pi*y)
                 + Real(-2.5)*pi*std::sin(pi*x)*std::cos(Real(2.5)*pi*y));
#elif (AMREX_SPACEDIM == 3)
            Real rhs2 = Real(-4.)*pi*std::cos(Real(2)*pi*x)*std::cos(Real(3.)*pi*y)*std::sin(Real(4.)*pi*z) *
                (Real(3.5)*pi*std::cos(Real(3.5)*pi*x)*std::cos(Real(3.5)*pi*y)*std::sin(Real(4.)*pi*z+Real(1./6.)*pi) +
                 Real(-4.)*pi*std::cos(Real(2.5)*pi*x)*std::sin(Real(3.)*pi*y)*std::cos(Real(4.)*pi*z+Real(0.25)*pi))
                + Real(2)*pi*std::sin(Real(2)*pi*x)*std::cos(Real(3.)*pi*y)*std::cos(Real(4.)*pi*z) *
                (Real(-2.5)*pi*std::sin(Real(2.5)*pi*x)*std::sin(Real(3.)*pi*y)*std::sin(Real(4.)*pi*z+Real(0.25)*pi) +
                 Real(-2.5)*pi*std::sin(pi*x)*std::cos(Real(2.5)*pi*y)*std::sin(Real(2)*pi*z+Real(1./3.)*pi));
#endif
            rhs[1](i,j,k) += Real(0.3)*rhs2;
        }
    }

    if (rhs[2].contains(i,j,k)) {
        Real x = xnd;
#if (AMREX_SPACEDIM > 1)
        Real y = ynd;
#endif
#if (AMREX_SPACEDIM == 1)
        Real cce = Real(12.25)*pi*pi*std::cos(Real(3.5)*pi*x);
#elif (AMREX_SPACEDIM == 2)
        Real cce = Real(24.5)*pi*pi*std::cos(Real(3.5)*pi*x)*std::sin(Real(3.5)*pi*y);
#else
        Real z = zcc;
        Real cce = Real(24.5)*pi*pi*std::cos(Real(3.5)*pi*x)*std::sin(Real(3.5)*pi*y)*std::sin(Real(4.)*pi*z+Real(1./6.)*pi)
            + Real(2.)*pi*pi*std::cos(pi*x)*std::sin(Real(2.5)*pi*y)*std::cos(Real(2.)*pi*z+Real(1./3.)*pi)
            + Real(12.)*pi*pi*std::cos(Real(2.5)*pi*x)*std::cos(Real(3.)*pi*y)*std::cos(Real(4.)*pi*z+Real(0.25)*pi);
#endif

        Real a = (alpha_node) ? get_alpha(AMREX_D_DECL(x,y,z)) : alpha;
        rhs[2](i,j,k) = a*cce + beta*sol[2](i,j,k);

        if (alpha_node) {
#if (AMREX_SPACEDIM == 1)
            Real rhs2 = Real(-7.0)*pi*pi*std::sin(Real(2.0)*pi*x)*std::sin(Real(3.5)*pi*x);
#elif (AMREX_SPACEDIM == 2)
            Real rhs2 = Real(-2)*pi*Real(3.5)*pi *
                std::sin(Real(2.)*pi*x)*std::cos(Real(3.)*pi*y)*
                std::sin(Real(3.5)*pi*x)*std::sin(Real(3.5)*pi*y)
                + Real(3)*pi*Real(3.5)*pi *
                std::cos(Real(2.)*pi*x)*std::sin(Real(3.)*pi*y)*
                std::cos(Real(3.5)*pi*x)*std::cos(Real(3.5)*pi*y);
#elif (AMREX_SPACEDIM == 3)
            Real rhs2 = Real(-2.)*pi*std::sin(Real(2.)*pi*x)*std::cos(Real(3.)*pi*y)*std::cos(Real(4.)*pi*z) *
                (Real(2.)*pi*std::sin(pi*x)*std::sin(Real(2.5)*pi*y)*std::cos(Real(2.)*pi*z+Real(1./3.)*pi) +
                 Real(3.5)*pi*std::sin(Real(3.5)*pi*x)*std::sin(Real(3.5)*pi*y)*std::sin(Real(4.)*pi*z+Real(1./6.)*pi))
                + Real(3.)*pi*std::cos(Real(2.)*pi*x)*std::sin(Real(3.)*pi*y)*std::cos(Real(4.)*pi*z) *
                (Real(3.5)*pi*std::cos(Real(3.5)*pi*x)*std::cos(Real(3.5)*pi*y)*std::sin(Real(4.)*pi*z+Real(1./6.)*pi) -
                 Real(4.)*pi*std::cos(Real(2.5)*pi*x)*std::sin(Real(3.)*pi*y)*std::cos(Real(4.)*pi*z+Real(0.25)*pi));
#endif
            rhs[2](i,j,k) += Real(0.3)*rhs2;
        }
    }
}

#if (AMREX_SPACEDIM == 1)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void actual_init_prob_sph1d (int i, int j, int k,
                             amrex::GpuArray<amrex::Array4<amrex::Real>,3> const& rhs,
                             amrex::GpuArray<amrex::Array4<amrex::Real>,3> const& sol,
                             amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& problo,
                             amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
                             amrex::Real alpha, amrex::Real beta,
                             amrex::IntVect const& ndhi,
                             amrex::Array4<amrex::Real> const& alpha_node)
{
    using namespace amrex;

    constexpr Real pi = amrex::Math::pi<Real>();

    auto alpha_val = [=] AMREX_GPU_DEVICE (Real r) noexcept
    {
        return Real(1.0) + Real(0.3)*std::cos(Real(2.0)*pi*r);
    };

    auto dalpha_val = [=] AMREX_GPU_DEVICE (Real r) noexcept
    {
        return Real(-0.6)*pi*std::sin(Real(2.0)*pi*r);
    };

    Real xnd = problo[0] + Real(i)*dx[0];
    Real xcc = xnd + Real(0.5)*dx[0];

    if (alpha_node && alpha_node.contains(i,j,k)) {
        alpha_node(i,j,k) = alpha_val(xnd);
    }

    if (sol[0].contains(i,j,k)) {
        Real x = xcc;
        Real Ex = std::sin(pi*x);
        sol[0](i,j,k) = Ex;
    }

    if (sol[1].contains(i,j,k)) {
        Real x = xnd;
        Real Ey = std::sin(Real(2.0)*pi*x);
        sol[1](i,j,k) = Ey;
    }

    if (sol[2].contains(i,j,k)) {
        Real x = xnd;
        Real Ez = std::sin(Real(3.0)*pi*x);
        sol[2](i,j,k) = Ez;
    }

    if (rhs[0].contains(i,j,k)) {
        rhs[0](i,j,k) = beta*sol[0](i,j,k);
    }

    if (rhs[1].contains(i,j,k)) {
        Real x = xnd;
        if ((x == 0) || (i == ndhi[0])) {
            rhs[1](i,j,k) = 0;
        } else {
            auto s = Real(2.0)*pi;
            Real cce = -Real(2.0)*s*std::cos(s*x)/x + s*s*std::sin(s*x);
            Real alpha_loc = alpha;
            Real dalpha = Real(0.0);
            if (alpha_node) {
                alpha_loc = alpha_val(x);
                dalpha = dalpha_val(x);
            }
            Real rhs_val = alpha_loc*cce + beta*sol[1](i,j,k);
            if (alpha_node) {
                // variable-alpha correction = +alpha'(r)*(curl E)_phi = -alpha'(r)*(dEy/dr + Ey/r)
                rhs_val -= dalpha * (s*std::cos(s*x) + std::sin(s*x)/x);
            }
            rhs[1](i,j,k) = rhs_val;
        }
    }

    if (rhs[2].contains(i,j,k)) {
        Real x = xnd;
        if ((x == 0) || (i == ndhi[0])) {
            rhs[2](i,j,k) = 0;
        } else {
            auto s = Real(3.0)*pi;
            Real cce = -Real(2.0)*s*std::cos(s*x)/x + s*s*std::sin(s*x);
            Real alpha_loc = alpha;
            Real dalpha = Real(0.0);
            if (alpha_node) {
                alpha_loc = alpha_val(x);
                dalpha = dalpha_val(x);
            }
            Real rhs_val = alpha_loc*cce + beta*sol[2](i,j,k);
            if (alpha_node) {
                // variable-alpha correction = +alpha'(r)*(curl E)_theta = -alpha'(r)*(dEz/dr + Ez/r)
                rhs_val -= dalpha * (s*std::cos(s*x) + std::sin(s*x)/x);
            }
            rhs[2](i,j,k) = rhs_val;
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void actual_init_prob_cyl1d (int i, int j, int k,
                             amrex::GpuArray<amrex::Array4<amrex::Real>,3> const& rhs,
                             amrex::GpuArray<amrex::Array4<amrex::Real>,3> const& sol,
                             amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& problo,
                             amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
                             amrex::Real alpha, amrex::Real beta,
                             amrex::IntVect const& ndhi,
                             amrex::Array4<amrex::Real> const& alpha_node)
{
    using namespace amrex;

    constexpr Real pi = amrex::Math::pi<Real>();

    auto alpha_val = [=] AMREX_GPU_DEVICE (Real r) noexcept
    {
        return Real(1.0) + Real(0.3)*std::cos(Real(2.0)*pi*r);
    };

    auto dalpha_val = [=] AMREX_GPU_DEVICE (Real r) noexcept
    {
        return Real(-0.6)*pi*std::sin(Real(2.0)*pi*r);
    };

    Real xnd = problo[0] + Real(i)*dx[0];
    Real xcc = xnd + Real(0.5)*dx[0];

    if (alpha_node && alpha_node.contains(i,j,k)) {
        alpha_node(i,j,k) = alpha_val(xnd);
    }

    if (sol[0].contains(i,j,k)) {
        Real x = xcc;
        Real Ex = std::sin(pi*x);
        sol[0](i,j,k) = Ex;
    }

    if (sol[1].contains(i,j,k)) {
        Real x = xnd;
        Real Ey = std::sin(Real(2.0)*pi*x);
        sol[1](i,j,k) = Ey;
    }

    if (sol[2].contains(i,j,k)) {
        Real x = xnd;
        Real Ez = std::cos(Real(2.5)*pi*x);
        sol[2](i,j,k) = Ez;
    }

    if (rhs[0].contains(i,j,k)) {
        rhs[0](i,j,k) = beta*sol[0](i,j,k);
    }

    if (rhs[1].contains(i,j,k)) {
        Real x = xnd;
        if ((x == 0) || (i == ndhi[0])) {
            rhs[1](i,j,k) = 0;
        } else {
            auto s = Real(2.0)*pi;
            Real cce = std::sin(s*x)*(s*s+Real(1.0)/(x*x)) - std::cos(s*x)*s/x;
            Real alpha_loc = alpha;
            Real dalpha = Real(0.0);
            if (alpha_node) {
                alpha_loc = alpha_val(x);
                dalpha = dalpha_val(x);
            }
            Real rhs_val = alpha_loc*cce + beta*sol[1](i,j,k);
            if (alpha_node) {
                // variable-alpha correction = +alpha'(r)*(curl E)_z = -alpha'(r)*(dEy/dr + Ey/r)
                rhs_val -= dalpha * (s*std::cos(s*x) + std::sin(s*x)/x);
            }
            rhs[1](i,j,k) = rhs_val;
        }
    }

    if (rhs[2].contains(i,j,k)) {
        if (i == ndhi[0]) {
            rhs[2](i,j,k) = 0;
        } else {
            Real x = xnd;
            Real s = Real(2.5)*pi;
            Real cce;
            if (x == 0) {
                cce = 2*s*s;
            } else {
                cce = std::sin(s*x)*s/x + std::cos(s*x)*s*s;
            }
            Real alpha_loc = alpha;
            Real dalpha = Real(0.0);
            if (alpha_node) {
                alpha_loc = alpha_val(x);
                dalpha = dalpha_val(x);
            }
            Real rhs_val = alpha_loc*cce + beta*sol[2](i,j,k);
            if (alpha_node) {
                // variable-alpha correction = +alpha'(r)*(curl E)_theta = -alpha'(r)*dE_z/dr
                rhs_val += dalpha * s * std::sin(s*x);
            }
            rhs[2](i,j,k) = rhs_val;
        }
    }
}

#endif

#endif
