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

template <typename T>
void mlndlap_bc_doit (Box const& vbx, Array4<T> const& a, Box const& domain,
                      GpuArray<bool,AMREX_SPACEDIM> const& bflo,
                      GpuArray<bool,AMREX_SPACEDIM> const& bfhi) noexcept
{
    if (domain.strictly_contains(vbx)) return;

    int offset = domain.cellCentered() ? 0 : 1;

    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);
    Box const& sbox = amrex::grow(vbx,1);

    Box xlo_face = domain;
    xlo_face.setSmall(0,dlo.x-1);
    xlo_face.setBig  (0,dlo.x-1);
    xlo_face &= sbox;
    Box ylo_face = domain;
    ylo_face.setSmall(1,dlo.y-1);
    ylo_face.setBig  (1,dlo.y-1);
    ylo_face &= sbox;
    int xoffset = vbx.length(0)+1;
    int yoffset = vbx.length(1)+1;

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
    xlo_face, txbxlo,
    {
        auto lo = amrex::lbound(txbxlo);
        auto hi = amrex::ubound(txbxlo);
        if (lo.x == dlo.x-1 and bflo[0]) {
            for (int j = lo.y; j <= hi.y; ++j) {
                a(dlo.x-1,j,0) = a(dlo.x+offset,j,0);
            }
        }
        if (lo.x+xoffset == dhi.x+1 and bfhi[0]) {
            for (int j = lo.y; j <= hi.y; ++j) {
                a(dhi.x+1,j,0) = a(dhi.x-offset,j,0);
            }
        }
    },
    ylo_face, tybxlo,
    {
        auto lo = amrex::lbound(tybxlo);
        auto hi = amrex::ubound(tybxlo);
        if (lo.y == dlo.y-1 and bflo[1]) {
            for (int i = lo.x; i <= hi.x; ++i) {
                a(i,dlo.y-1,0) = a(i,dlo.y+offset,0);
            }
        }
        if (lo.y+yoffset == dhi.y+1 and bfhi[1]) {
            for (int i = lo.x; i <= hi.x; ++i) {
                a(i,dhi.y+1,0) = a(i,dhi.y-offset,0);
            }
        }
    });

    const auto lo = amrex::lbound(sbox);
    const auto hi = amrex::ubound(sbox);

    AMREX_HOST_DEVICE_FOR_1D ( 4, icorner,
    {
        switch (icorner) {
        case 0: {
            // xlo & ylo
            if (lo.x == dlo.x-1 and lo.y == dlo.y-1) {
                if (bflo[0]) {
                    a(dlo.x-1,dlo.y-1,0) = a(dlo.x+offset,dlo.y-1,0);
                } else if (bflo[1]) {
                    a(dlo.x-1,dlo.y-1,0) = a(dlo.x-1,dlo.y+offset,0);
                }
            }
            break;
        }
        case 1: {
            // xhi & ylo
            if (hi.x == dhi.x+1 and lo.y == dlo.y-1) {
                if (bfhi[0]) {
                    a(dhi.x+1,dlo.y-1,0) = a(dhi.x-offset,dlo.y-1,0);
                } else if (bflo[1]) {
                    a(dhi.x+1,dlo.y-1,0) = a(dhi.x+1,dlo.y+offset,0);
                }
            }
            break;
        }
        case 2: {
            // xlo & yhi
            if (lo.x == dlo.x-1 and hi.y == dhi.y+1) {
                if (bflo[0]) {
                    a(dlo.x-1,dhi.y+1,0) = a(dlo.x+offset,dhi.y+1,0);
                } else if (bfhi[1]) {
                    a(dlo.x-1,dhi.y+1,0) = a(dlo.x-1,dhi.y-offset,0);
                }
            }
            break;
        }
        case 3: {
            // xhi & yhi
            if (hi.x == dhi.x+1 and hi.y == dhi.y+1) {
                if (bfhi[0]) {
                    a(dhi.x+1,dhi.y+1,0) = a(dhi.x-offset,dhi.y+1,0);
                } else if (bfhi[1]) {
                    a(dhi.x+1,dhi.y+1,0) = a(dhi.x+1,dhi.y-offset,0);
                }
            }
            break;
        }
        default: {}
        }
    });
}

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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction (int i, int j, int k, Array4<Real> const& crse,
                          Array4<Real const> const& fine, Array4<int const> const& msk) noexcept
{
    int ii = i*2;
    int jj = j*2;
    int kk = 0;
    if (msk(ii,jj,kk)) {
        crse(i,j,k) = 0.0;
    } else {
        crse(i,j,k) = (1./16.)*(fine(ii-1,jj-1,kk) + 2.*fine(ii  ,jj-1,kk) +    fine(ii+1,jj-1,kk)
                           + 2.*fine(ii-1,jj  ,kk) + 4.*fine(ii  ,jj  ,kk) + 2.*fine(ii+1,jj  ,kk)
                              + fine(ii-1,jj+1,kk) + 2.*fine(ii  ,jj+1,kk) +    fine(ii+1,jj+1,kk));
    }
}

//
// interpolation
//

namespace {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_x (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int ic, int jc) noexcept
    {
        Real w1 = sig(i-1,j-1,0) + sig(i-1,j,0);
        Real w2 = sig(i  ,j-1,0) + sig(i  ,j,0);
        return (w1*crse(ic,jc,0)+w2*crse(ic+1,jc,0))/(w1+w2);
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_y (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int ic, int jc) noexcept
    {
        Real w1 = sig(i-1,j-1,0) + sig(i,j-1,0);
        Real w2 = sig(i-1,j  ,0) + sig(i,j  ,0);
        return (w1*crse(ic,jc,0)+w2*crse(ic,jc+1,0))/(w1+w2);
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_face_xy (Array4<Real const> const& crse, Array4<Real const> const& sig,
                            int i, int j, int ic, int jc) noexcept
    {
        Real w1 = sig(i-1,j-1,0) + sig(i-1,j,0);
        Real w2 = sig(i  ,j-1,0) + sig(i  ,j,0);
        Real w3 = sig(i-1,j-1,0) + sig(i,j-1,0);
        Real w4 = sig(i-1,j  ,0) + sig(i,j  ,0);
        return (w1 * aa_interp_line_y(crse,sig,i-1,j  ,ic  ,jc  ) +
                w2 * aa_interp_line_y(crse,sig,i+1,j  ,ic+1,jc  ) +
                w3 * aa_interp_line_x(crse,sig,i  ,j-1,ic  ,jc  ) +
                w4 * aa_interp_line_x(crse,sig,i  ,j+1,ic  ,jc+1)) / (w1+w2+w3+w4);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_aa (int i, int j, int, Array4<Real> const& fine,
                           Array4<Real const> const& crse, Array4<Real const> const& sig,
                           Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,0)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        if (i_is_odd and j_is_odd) {
            // Node on a X-Y face
            fine(i,j,0) += aa_interp_face_xy(crse,sig,i,j,ic,jc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,0) += aa_interp_line_x(crse,sig,i,j,ic,jc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,0) += aa_interp_line_y(crse,sig,i,j,ic,jc);
        } else {
            // Node coincident with coarse node
            fine(i,j,0) += crse(ic,jc,0);
        }
    }
}

namespace {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real ha_interp_face_xy (Array4<Real const> const& crse,
                            Array4<Real const> const& sigx, Array4<Real const> const& sigy,
                            int i, int j, int ic, int jc) noexcept
    {
        Real w1 = sigx(i-1,j-1,0) + sigx(i-1,j,0);
        Real w2 = sigx(i  ,j-1,0) + sigx(i  ,j,0);
        Real w3 = sigy(i-1,j-1,0) + sigy(i,j-1,0);
        Real w4 = sigy(i-1,j  ,0) + sigy(i,j  ,0);
        return (w1 * aa_interp_line_y(crse,sigy,i-1,j  ,ic  ,jc  ) +
                w2 * aa_interp_line_y(crse,sigy,i+1,j  ,ic+1,jc  ) +
                w3 * aa_interp_line_x(crse,sigx,i  ,j-1,ic  ,jc  ) +
                w4 * aa_interp_line_x(crse,sigx,i  ,j+1,ic  ,jc+1)) / (w1+w2+w3+w4);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_ha (int i, int j, int,
                           Array4<Real> const& fine, Array4<Real const> const& crse,
                           Array4<Real const> const& sigx, Array4<Real const> const& sigy,
                           Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,0)) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        bool i_is_odd = (ic*2 != i);
        bool j_is_odd = (jc*2 != j);
        if (i_is_odd and j_is_odd) {
            // Node on a X-Y face
            fine(i,j,0) += ha_interp_face_xy(crse,sigx,sigy,i,j,ic,jc);
        } else if (i_is_odd) {
            // Node on X line
            fine(i,j,0) += aa_interp_line_x(crse,sigx,i,j,ic,jc);
        } else if (j_is_odd) {
            // Node on Y line
            fine(i,j,0) += aa_interp_line_y(crse,sigy,i,j,ic,jc);
        } else {
            // Node coincident with coarse node
            fine(i,j,0) += crse(ic,jc,0);
        }
    }
}

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

//
// RAP
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil (Box const& bx, Array4<Real> const& sten,
                          Array4<Real const> const& sigma,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = (1.0/6.0)*dxinv[0]*dxinv[0];
    Real facy = (1.0/6.0)*dxinv[1]*dxinv[1];
    Real fxy = facx + facy;
    Real f2xmy = 2.0*facx - facy;
    Real fmx2y = 2.0*facy - facx;

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        sten(i,j,k,1) = f2xmy*(sigma(i,j-1,k)+sigma(i,j,k));
        sten(i,j,k,2) = fmx2y*(sigma(i-1,j,k)+sigma(i,j,k));
        sten(i,j,k,3) = fxy*sigma(i,j,k);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil_s0 (int i, int j, int k, Array4<Real> const& sten) noexcept
{
    sten(i,j,k,0) = -(sten(i-1,j  ,k,1) + sten(i  ,j  ,k,1)
                    + sten(i  ,j-1,k,2) + sten(i  ,j  ,k,2)
                    + sten(i-1,j-1,k,3) + sten(i  ,j-1,k,3)
                    + sten(i-1,j  ,k,3) + sten(i  ,j  ,k,3));
    sten(i,j,k,4) = 1.0 / (std::abs(sten(i-1,j  ,k,1)) + std::abs(sten(i,j  ,k,1))
                         + std::abs(sten(i  ,j-1,k,2)) + std::abs(sten(i,j  ,k,2))
                         + std::abs(sten(i-1,j-1,k,3)) + std::abs(sten(i,j-1,k,3))
                         + std::abs(sten(i-1,j  ,k,3)) + std::abs(sten(i,j  ,k,3)) + eps);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_sten (Box const& bx, Array4<Real> const& y, Array4<Real const> const& x,
                         Array4<Real const> const& sten, Array4<int const> const& msk) noexcept
{
    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            y(i,j,k) = 0.0;
        } else {
            y(i,j,k) = x(i-1,j-1,k)*sten(i-1,j-1,k,3)
                +      x(i  ,j-1,k)*sten(i  ,j-1,k,2)
                +      x(i+1,j-1,k)*sten(i  ,j-1,k,3)
                +      x(i-1,j  ,k)*sten(i-1,j  ,k,1)
                +      x(i  ,j  ,k)*sten(i  ,j  ,k,0)
                +      x(i+1,j  ,k)*sten(i  ,j  ,k,1)
                +      x(i-1,j+1,k)*sten(i-1,j  ,k,3)
                +      x(i  ,j+1,k)*sten(i  ,j  ,k,2)
                +      x(i+1,j+1,k)*sten(i  ,j  ,k,3);
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_sten (Box const& bx, Array4<Real> const& sol,
                                Array4<Real const> const& rhs,
                                Array4<Real const> const& sten,
                                Array4<int const> const& msk) noexcept
{
    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = 0.0;
        } else if (sten(i,j,k,0) != 0.0) {
            Real Ax = sol(i-1,j-1,k)*sten(i-1,j-1,k,3)
                +     sol(i  ,j-1,k)*sten(i  ,j-1,k,2)
                +     sol(i+1,j-1,k)*sten(i  ,j-1,k,3)
                +     sol(i-1,j  ,k)*sten(i-1,j  ,k,1)
                +     sol(i  ,j  ,k)*sten(i  ,j  ,k,0)
                +     sol(i+1,j  ,k)*sten(i  ,j  ,k,1)
                +     sol(i-1,j+1,k)*sten(i-1,j  ,k,3)
                +     sol(i  ,j+1,k)*sten(i  ,j  ,k,2)
                +     sol(i+1,j+1,k)*sten(i  ,j  ,k,3);
            sol(i,j,k) += (rhs(i,j,k) - Ax) / sten(i,j,k,0);
        }
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_interpadd_rap (int i, int j, int, Array4<Real> const& fine,
                            Array4<Real const> const& crse, Array4<Real const> const& sten,
                            Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,0) and sten(i,j,0,0) != 0.0) {
        int ic = amrex::coarsen(i,2);
        int jc = amrex::coarsen(j,2);
        bool ieven = ic*2 == i;
        bool jeven = jc*2 == j;
        Real fv;
        if (ieven and jeven) {
            fv = crse(ic,jc,0);
        } else if (ieven) {
            Real wym = std::abs(sten(i,j-1,0,2));
            Real wyp = std::abs(sten(i,j  ,0,2));
            fv = (wym*crse(ic,jc,0) + wyp*crse(ic,jc+1,0)) / (wym+wyp+eps);
        } else if (jeven) {
            Real wxm = std::abs(sten(i-1,j,0,1));
            Real wxp = std::abs(sten(i  ,j,0,1));
            fv = (wxm*crse(ic,jc,0) + wxp*crse(ic+1,jc,0)) / (wxm+wxp+eps);
        } else {
            Real wxm = std::abs(sten(i-1,j  ,0,1)) /
                (std::abs(sten(i-1,j-1,0,3))+std::abs(sten(i-1,j  ,0,3))+eps);
            Real wxp = std::abs(sten(i  ,j  ,0,1)) /
                (std::abs(sten(i  ,j-1,0,3))+std::abs(sten(i  ,j  ,0,3))+eps);
            Real wym = std::abs(sten(i  ,j-1,0,2)) /
                (std::abs(sten(i-1,j-1,0,3))+std::abs(sten(i  ,j-1,0,3))+eps);
            Real wyp = std::abs(sten(i  ,j  ,0,2)) /
                (std::abs(sten(i-1,j  ,0,3))+std::abs(sten(i  ,j  ,0,3))+eps);
            Real wmm = std::abs(sten(i-1,j-1,0,3)) * (1.0 + wxm + wym);
            Real wpm = std::abs(sten(i,j-1,0,3)) * (1.0 + wxp + wym);
            Real wmp = std::abs(sten(i-1,j,0,3)) *(1.0 + wxm + wyp);
            Real wpp = std::abs(sten(i,j,0,3)) * (1.0 + wxp + wyp);
            fv = (wmm*crse(ic,jc,0) + wpm*crse(ic+1,jc,0)
                  + wmp*crse(ic,jc+1,0) + wpp*crse(ic+1,jc+1,0))
                / (wmm+wpm+wmp+wpp+eps);
        }

        fine(i,j,0) += fv;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction_rap (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine, Array4<Real const> const& sten,
                              Array4<int const> const& msk) noexcept
{
    int ii = i*2;
    int jj = j*2;
    if (msk(ii,jj,0)) {
        crse(i,j,0) = 0.0;
    } else {

        Real cv = fine(ii,jj,0)
            + fine(ii-1,jj  ,0)*std::abs(sten(ii-1,jj  ,0,1))
            /                  (std::abs(sten(ii-2,jj  ,0,1))
                               +std::abs(sten(ii-1,jj  ,0,1))+eps)
            + fine(ii+1,jj  ,0)*std::abs(sten(ii  ,jj  ,0,1))
            /                  (std::abs(sten(ii  ,jj  ,0,1))
                               +std::abs(sten(ii+1,jj  ,0,1))+eps)
            + fine(ii  ,jj-1,0)*std::abs(sten(ii  ,jj-1,0,2))
            /                  (std::abs(sten(ii  ,jj-2,0,2))
                               +std::abs(sten(ii  ,jj-1,0,2))+eps)
            + fine(ii  ,jj+1,0)*std::abs(sten(ii  ,jj  ,0,2))
            /                  (std::abs(sten(ii  ,jj  ,0,2))
                               +std::abs(sten(ii  ,jj+1,0,2))+eps);

        Real wxp = std::abs(sten(ii-1,jj-1,0,1))
            /     (std::abs(sten(ii-1,jj-2,0,3))
                  +std::abs(sten(ii-1,jj-1,0,3))+eps);
        Real wyp = std::abs(sten(ii-1,jj-1,0,2))
            /     (std::abs(sten(ii-2,jj-1,0,3))
                  +std::abs(sten(ii-1,jj-1,0,3))+eps);
        Real wpp = std::abs(sten(ii-1,jj-1,0,3))*(1.0+wxp+wyp);
        cv +=           wpp*sten(ii-1,jj-1,0,4)*fine(ii-1,jj-1,0);

        Real wxm = std::abs(sten(ii  ,jj-1,0,1))
            /     (std::abs(sten(ii  ,jj-2,0,3))
                  +std::abs(sten(ii  ,jj-1,0,3))+eps);
        wyp      = std::abs(sten(ii+1,jj-1,0,2))
            /     (std::abs(sten(ii  ,jj-1,0,3))
                  +std::abs(sten(ii+1,jj-1,0,3))+eps);
        Real wmp = std::abs(sten(ii  ,jj-1,0,3))*(1.0 + wxm + wyp);
        cv +=           wmp*sten(ii+1,jj-1,0,4)*fine(ii+1,jj-1,0);

        wxp      = std::abs(sten(ii-1,jj+1,0,1))
            /     (std::abs(sten(ii-1,jj  ,0,3))
                  +std::abs(sten(ii-1,jj+1,0,3))+eps);
        Real wym = std::abs(sten(ii-1,jj  ,0,2))
            /     (std::abs(sten(ii-2,jj  ,0,3))
                  +std::abs(sten(ii-1,jj  ,0,3))+eps);
        Real wpm = std::abs(sten(ii-1,jj  ,0,3)) * (1.0 + wxp + wym);
        cv +=           wpm*sten(ii-1,jj+1,0,4)*fine(ii-1,jj+1,0);

        wxm      = std::abs(sten(ii  ,jj+1,0,1))
            /     (std::abs(sten(ii  ,jj  ,0,3))
                  +std::abs(sten(ii  ,jj+1,0,3))+eps);
        wym      = std::abs(sten(ii+1,jj  ,0,2))
            /     (std::abs(sten(ii  ,jj  ,0,3))
                  +std::abs(sten(ii+1,jj  ,0,3))+eps);
        Real wmm = std::abs(sten(ii  ,jj  ,0,3)) * (1.0 + wxm + wym);
        cv +=           wmm*sten(ii+1,jj+1,0,4)*fine(ii+1,jj+1,0);

        crse(i,j,0) = cv * 0.25;
    }
}

}
#endif
