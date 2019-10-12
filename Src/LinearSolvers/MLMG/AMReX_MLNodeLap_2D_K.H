#ifndef AMREX_MLNODELAP_2D_K_H_
#define AMREX_MLNODELAP_2D_K_H_

namespace amrex {

//
// masks
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_nodal_mask (int i, int j, int k, Array4<int> const& nmsk,
                             Array4<int const> const& cmsk) noexcept
{
    int s = cmsk(i-1,j-1,k) + cmsk(i  ,j-1,k)
        +   cmsk(i-1,j  ,k) + cmsk(i  ,j  ,k);
    if (s == 4*crse_cell) {
        nmsk(i,j,k) = crse_node;
    }
    else if (s == 4*fine_cell) {
        nmsk(i,j,k) = fine_node;
    } else {
        nmsk(i,j,k) = crse_fine_node;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_dirichlet_mask (Box const& bx, Array4<int> const& dmsk,
                                 Array4<int const> const& omsk, Box const& dom,
                                 Array<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                 Array<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        dmsk(i,j,0) = (omsk(i-1,j-1,0) == 1 or omsk(i,j-1,0) == 1 or
                       omsk(i-1,j  ,0) == 1 or omsk(i,j  ,0) == 1);
    }}

    const auto domlo = amrex::lbound(dom);
    const auto domhi = amrex::ubound(dom);

    if (bclo[0] == LinOpBCType::Dirichlet and lo.x == domlo.x) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(lo.x,j,0) = 1;
        }
    }

    if (bchi[0] == LinOpBCType::Dirichlet and hi.x == domhi.x) {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(hi.x,j,0) = 1;
        }
    }

    if (bclo[1] == LinOpBCType::Dirichlet and lo.y == domlo.y) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,lo.y,0) = 1;
        }
    }

    if (bchi[1] == LinOpBCType::Dirichlet and hi.y == domhi.y) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,hi.y,0) = 1;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_dot_mask (Box const& bx, Array4<Real> const& dmsk,
                           Array4<int const> const& omsk, Box const& dom,
                           Array<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                           Array<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        dmsk(i,j,0) = static_cast<Real>(omsk(i,j,0));
    }}

    const auto domlo = amrex::lbound(dom);
    const auto domhi = amrex::ubound(dom);

    if ((bclo[0] == LinOpBCType::Neumann or bclo[0] == LinOpBCType::inflow)
        and lo.x == domlo.x)
    {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(lo.x,j,0) *= 0.5;
        }
    }

    if ((bchi[0] == LinOpBCType::Neumann or bchi[0] == LinOpBCType::inflow)
        and hi.x == domhi.x)
    {
        for (int j = lo.y; j <= hi.y; ++j) {
            dmsk(hi.x,j,0) *= 0.5;
        }
    }

    if ((bclo[1] == LinOpBCType::Neumann or bclo[1] == LinOpBCType::inflow)
        and lo.y == domlo.y)
    {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,lo.y,0) *= 0.5;
        }
    }

    if ((bchi[1] == LinOpBCType::Neumann or bchi[1] == LinOpBCType::inflow)
        and hi.y == domhi.y)
    {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            dmsk(i,hi.y,0) *= 0.5;
        }
    }
}

//
// coeffs
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_avgdown_coeff_x (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine) noexcept
{
    Real a = fine(2*i  ,2*j,k) + fine(2*i  ,2*j+1,k);
    Real b = fine(2*i+1,2*j,k) + fine(2*i+1,2*j+1,k);
    crse(i,j,k) = a*b/(a+b);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_avgdown_coeff_y (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine) noexcept
{
    Real a = fine(2*i,2*j  ,k) + fine(2*i+1,2*j  ,k);
    Real b = fine(2*i,2*j+1,k) + fine(2*i+1,2*j+1,k);
    crse(i,j,k) = a*b/(a+b);
}

//
// bc
//

//
// operator
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_ha (Box const& bx, Array4<Real> const& y, Array4<Real const> const& x,
                       Array4<Real const> const& sx, Array4<Real const> const& sy,
                       Array4<int const> const& msk, bool is_rz,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1./6.)*dxinv[0]*dxinv[0];
    Real facy = (1./6.)*dxinv[1]*dxinv[1];
    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            y(i,j,k) = 0.0;
        } else {
            y(i,j,k) = x(i-1,j-1,k)*(facx*sx(i-1,j-1,k)+facy*sy(i-1,j-1,k))
                   +   x(i+1,j-1,k)*(facx*sx(i  ,j-1,k)+facy*sy(i  ,j-1,k))
                   +   x(i-1,j+1,k)*(facx*sx(i-1,j  ,k)+facy*sy(i-1,j  ,k))
                   +   x(i+1,j+1,k)*(facx*sx(i  ,j  ,k)+facy*sy(i  ,j  ,k))
                   +   x(i-1,j,k)*(2.0*facx*(sx(i-1,j-1,k)+sx(i-1,j,k))
                                 -     facy*(sy(i-1,j-1,k)+sx(i-1,j,k)))
                   +   x(i+1,j,k)*(2.0*facx*(sx(i  ,j-1,k)+sx(i  ,j,k))
                                 -     facy*(sy(i  ,j-1,k)+sx(i  ,j,k)))
                   +   x(i,j-1,k)*(   -facx*(sx(i-1,j-1,k)+sx(i,j-1,k))
                                  +2.0*facy*(sy(i-1,j-1,k)+sy(i,j-1,k)))
                   +   x(i,j+1,k)*(   -facx*(sx(i-1,j  ,k)+sx(i,j  ,k))
                                  +2.0*facy*(sy(i-1,j  ,k)+sy(i,j  ,k)))
                   +   x(i,j,k)*(-2.0)*(facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
                                       +facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));
            if (is_rz) {
                Real fp = facy / static_cast<Real>(2*i+1);
                Real fm = facy / static_cast<Real>(2*i-1);
                y(i,j,k) += (fm*sy(i-1,j  ,k)-fp*sy(i,j  ,k))*(x(i,j+1,k)-x(i,j,k))
                          + (fm*sy(i-1,j-1,k)-fp*sy(i,j-1,k))*(x(i,j-1,k)-x(i,j,k));
            }
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_aa (Box const& bx, Array4<Real> const& y, Array4<Real const> const& x,
                       Array4<Real const> const& sig, Array4<int const> const& msk,
                       bool is_rz, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/6.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/6.0)*dxinv[1]*dxinv[1];
    Real fxy = facx + facy;
    Real f2xmy = 2.0*facx - facy;
    Real fmx2y = 2.0*facy - facx;

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            y(i,j,k) = 0.0;
        } else {
            y(i,j,k) = x(i-1,j-1,k)*fxy*sig(i-1,j-1,k)
                   +   x(i+1,j-1,k)*fxy*sig(i  ,j-1,k)
                   +   x(i-1,j+1,k)*fxy*sig(i-1,j  ,k)
                   +   x(i+1,j+1,k)*fxy*sig(i  ,j  ,k)
                   +   x(i-1,j,k)*f2xmy*(sig(i-1,j-1,k)+sig(i-1,j,k))
                   +   x(i+1,j,k)*f2xmy*(sig(i  ,j-1,k)+sig(i  ,j,k))
                   +   x(i,j-1,k)*fmx2y*(sig(i-1,j-1,k)+sig(i,j-1,k))
                   +   x(i,j+1,k)*fmx2y*(sig(i-1,j  ,k)+sig(i,j  ,k))
                   +   x(i,j,k)*(-2.0)*fxy*(sig(i-1,j-1,k)+sig(i,j-1,k)
                                            +sig(i-1,j,k)+sig(i,j,k));
            if (is_rz) {
                Real fp = facy / static_cast<Real>(2*i+1);
                Real fm = facy / static_cast<Real>(2*i-1);
                y(i,j,k) += (fm*sig(i-1,j  ,k)-fp*sig(i,j  ,k))*(x(i,j+1,k)-x(i,j,k))
                          + (fm*sig(i-1,j-1,k)-fp*sig(i,j-1,k))*(x(i,j-1,k)-x(i,j,k));
            }
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_ha (Box const& bx, Array4<Real> const& x, Array4<Real const> const& sx,
                           Array4<Real const> const& sy, Array4<int const> const& msk,
                           GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/6.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/6.0)*dxinv[1]*dxinv[1];

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (!msk(i,j,k)) {
            x(i,j,k) /= (-2.0)*(facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
                               +facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_normalize_aa (Box const& bx, Array4<Real> const& x, Array4<Real const> const& sig,
                           Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/6.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/6.0)*dxinv[1]*dxinv[1];
    Real fxy = facx + facy;

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (!msk(i,j,k)) {
            x(i,j,k) /= (-2.0)*fxy*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_ha (Box const& bx, Array4<Real> const& sol, Array4<Real const> const& Ax,
                        Array4<Real const> const& rhs, Array4<Real const> const& sx,
                        Array4<Real const> const& sy, Array4<int const> const& msk,
                        GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = -2.0 * (1.0/6.0)*dxinv[0]*dxinv[0];
    Real facy = -2.0 * (1.0/6.0)*dxinv[1]*dxinv[1];

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            sol(i,j,k) += (2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
                / (facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
                +  facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_jacobi_aa (Box const& bx, Array4<Real> const& sol, Array4<Real const> const& Ax,
                        Array4<Real const> const& rhs, Array4<Real const> const& sig,
                        Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real fac = -2.0 * (1.0/6.0)*(dxinv[0]*dxinv[0] + dxinv[1]*dxinv[1]);

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            sol(i,j,k) += (2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
                / (fac*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k)));
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_ha (Box const& bx, Array4<Real> const& sol,
                              Array4<Real const> const& rhs, Array4<Real const> const& sx,
                              Array4<Real const> const& sy, Array4<int const> const& msk,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              bool is_rz) noexcept
{
    Real facx = (1.0/6.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/6.0)*dxinv[1]*dxinv[1];

    amrex::Loop(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            Real s0 = (-2.0)*(facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
                             +facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));

            Real Ax = sol(i-1,j-1,k)*(facx*sx(i-1,j-1,k)+facy*sy(i-1,j-1,k))
                    + sol(i+1,j-1,k)*(facx*sx(i  ,j-1,k)+facy*sy(i  ,j-1,k))
                    + sol(i-1,j+1,k)*(facx*sx(i-1,j  ,k)+facy*sy(i-1,j  ,k))
                    + sol(i+1,j+1,k)*(facx*sx(i  ,j  ,k)+facy*sy(i  ,j  ,k))
                    + sol(i-1,j,k)*(2.0*facx*(sx(i-1,j-1,k)+sx(i-1,j,k))
                                  -     facy*(sy(i-1,j-1,k)+sx(i-1,j,k)))
                    + sol(i+1,j,k)*(2.0*facx*(sx(i  ,j-1,k)+sx(i  ,j,k))
                                  -     facy*(sy(i  ,j-1,k)+sx(i  ,j,k)))
                    + sol(i,j-1,k)*(   -facx*(sx(i-1,j-1,k)+sx(i,j-1,k))
                                   +2.0*facy*(sy(i-1,j-1,k)+sy(i,j-1,k)))
                    + sol(i,j+1,k)*(   -facx*(sx(i-1,j  ,k)+sx(i,j  ,k))
                                   +2.0*facy*(sy(i-1,j  ,k)+sy(i,j  ,k)))
                    + sol(i,j,k)*s0;

            if (is_rz) {
                Real fp = facy / static_cast<Real>(2*i+1);
                Real fm = facy / static_cast<Real>(2*i-1);
                Real frzlo = fm*sy(i-1,j-1,k)-fp*sy(i,j-1,k);
                Real frzhi = fm*sy(i-1,j  ,k)-fp*sy(i,j  ,k);
                s0 += - frzhi - frzlo;
                Ax += frzhi*(sol(i,j+1,k)-sol(i,j,k))
                    + frzlo*(sol(i,j-1,k)-sol(i,j,k));
            }

            sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_aa (Box const& bx, Array4<Real> const& sol,
                              Array4<Real const> const& rhs, Array4<Real const> const& sig,
                              Array4<int const> const& msk,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              bool is_rz) noexcept
{
    Real facx = (1.0/6.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/6.0)*dxinv[1]*dxinv[1];
    Real fxy = facx + facy;
    Real f2xmy = 2.0*facx - facy;
    Real fmx2y = 2.0*facy - facx;

    amrex::Loop(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else {
            Real s0 = (-2.0)*fxy*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));
            Real Ax =   sol(i-1,j-1,k)*fxy*sig(i-1,j-1,k)
                      + sol(i+1,j-1,k)*fxy*sig(i  ,j-1,k)
                      + sol(i-1,j+1,k)*fxy*sig(i-1,j  ,k)
                      + sol(i+1,j+1,k)*fxy*sig(i  ,j  ,k)
                      + sol(i-1,j,k)*f2xmy*(sig(i-1,j-1,k)+sig(i-1,j,k))
                      + sol(i+1,j,k)*f2xmy*(sig(i  ,j-1,k)+sig(i  ,j,k))
                      + sol(i,j-1,k)*fmx2y*(sig(i-1,j-1,k)+sig(i,j-1,k))
                      + sol(i,j+1,k)*fmx2y*(sig(i-1,j  ,k)+sig(i,j  ,k))
                      + sol(i,j,k)*s0;

            if (is_rz) {
                Real fp = facy / static_cast<Real>(2*i+1);
                Real fm = facy / static_cast<Real>(2*i-1);
                Real frzlo = fm*sig(i-1,j-1,k)-fp*sig(i,j-1,k);
                Real frzhi = fm*sig(i-1,j  ,k)-fp*sig(i,j  ,k);
                s0 += - frzhi - frzlo;
                Ax += frzhi*(sol(i,j+1,k)-sol(i,j,k))
                    + frzlo*(sol(i,j-1,k)-sol(i,j,k));
            }

            sol(i,j,k) += (rhs(i,j,k) - Ax) / s0;
        }
    });
}

//
// restriction
//

//
// interpolation
//

//
// rhs & u
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu (int i, int j, int k, Array4<Real> const& rhs, Array4<Real const> const& vel,
                   Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                   bool is_rz) noexcept
{
    Real facx = 0.5*dxinv[0];
    Real facy = 0.5*dxinv[1];

    if (msk(i,j,k)) {
        rhs(i,j,k) = 0.0;
    } else {
        rhs(i,j,k) = facx*(-vel(i-1,j-1,k,0) + vel(i,j-1,k,0)
                           -vel(i-1,j  ,k,0) + vel(i,j  ,k,0))
                   + facy*(-vel(i-1,j-1,k,1) - vel(i,j-1,k,1)
                           +vel(i-1,j  ,k,1) + vel(i,j  ,k,1));
        if (is_rz) {
            Real fm = facy / (6*i-3);
            Real fp = facy / (6*i+3);
            rhs(i,j,k) += fm*(vel(i-1,j,k,1)-vel(i-1,j-1,k,1))
                        - fp*(vel(i  ,j,k,1)-vel(i  ,j-1,k,1));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_rhcc (int i, int j, int k, Array4<Real> const& rhs, Array4<Real const> const& rhcc,
                   Array4<int const> const& msk) noexcept
{
    if (msk(i,j,k)) {
        rhs(i,j,k) = 0.0;
    } else {
        rhs(i,j,k) = 0.25 * (rhcc(i-1,j-1,k)+rhcc(i,j-1,k)+rhcc(i-1,j,k)+rhcc(i,j,k));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu (int i, int j, int k, Array4<Real> const& u, Array4<Real const> const& p,
                     Array4<Real const> const& sig, Real dxi, Real dyi, bool is_rz) noexcept
{
    Real facx = 0.5*dxi;
    Real facy = 0.5*dyi;
    u(i,j,k,0) -= sig(i,j,k)*facx*(-p(i,j,k)+p(i+1,j,k)-p(i,j+1,k)+p(i+1,j+1,k));
    u(i,j,k,1) -= sig(i,j,k)*facy*(-p(i,j,k)-p(i+1,j,k)+p(i,j+1,k)+p(i+1,j+1,k));
    if (is_rz) {
        Real frz = sig(i,j,k)*facy / (6*i+3);
        u(i,j,k,1) += frz*(p(i,j,k)-p(i+1,j,k)-p(i,j+1,k)+p(i+1,j+1,k));
    }
}

//
// residual
//

//
// sync residual
//

}
#endif
