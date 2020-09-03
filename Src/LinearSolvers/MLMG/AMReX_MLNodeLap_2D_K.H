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
                                 GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                 GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int j = lo.y; j <= hi.y; ++j) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        if (!dmsk(i,j,0)) {
            dmsk(i,j,0) = (omsk(i-1,j-1,0) == 1 or omsk(i,j-1,0) == 1 or
                           omsk(i-1,j  ,0) == 1 or omsk(i,j  ,0) == 1);
        }
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
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                           GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_zero_fine (int i, int j, int, Array4<Real> const& phi,
                        Array4<int const> const& msk, int fine_flag) noexcept
{
    // Testing if the node is covered by a fine level in computing
    // coarse sync residual
    if (msk(i-1,j-1,0) == fine_flag and
        msk(i  ,j-1,0) == fine_flag and
        msk(i-1,j  ,0) == fine_flag and
        msk(i  ,j  ,0) == fine_flag)
    {
        phi(i,j,0) = Real(0.0);
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_avgdown_coeff (int i, int j, int k, Array4<Real> const& crse,
                              Array4<Real const> const& fine, int idir) noexcept
{   
    if (idir == 1) {
        Real a = fine(2*i  ,j,k);
        Real b = fine(2*i+1,j,k);
        crse(i,j,k) = 2.0*a*b/(a+b);
    } else {
        Real a = fine(i,2*j  ,k);
        Real b = fine(i,2*j+1,k);
        crse(i,j,k) = 2.0*a*b/(a+b);
    }
}

//
// bc
//

template <typename T>
void mlndlap_bc_doit (Box const& vbx, Array4<T> const& a, Box const& domain,
                      GpuArray<bool,AMREX_SPACEDIM> const& bflo,
                      GpuArray<bool,AMREX_SPACEDIM> const& bfhi) noexcept
{
    Box gdomain = domain;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (not bflo[idim]) gdomain.growLo(idim,1);
        if (not bfhi[idim]) gdomain.growHi(idim,1);
    }

    if (gdomain.strictly_contains(vbx)) return;

    int offset = domain.cellCentered() ? 0 : 1;

    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    Box const& sbox = amrex::grow(vbx,1);
    AMREX_HOST_DEVICE_FOR_3D(sbox, i, j, k,
    {
        if (not gdomain.contains(IntVect(i,j))) {
            // xlo & ylo
            if (i == dlo.x-1 and j == dlo.y-1 and (bflo[0] or bflo[1]))
            {
                if (bflo[0] and bflo[1])
                {
                    a(i,j,k) = a(i+1+offset, j+1+offset, k);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
            }
            // xhi & ylo
            else if (i == dhi.x+1 and j == dlo.y-1 and (bfhi[0] or bflo[1]))
            {
                if (bfhi[0] and bflo[1])
                {
                    a(i,j,k) = a(i-1-offset, j+1+offset, k);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bflo[1])
                {
                    a(i,j,k) = a(i, j+1+offset, k);
                }
            }
            // xlo & yhi
            else if (i == dlo.x-1 and j == dhi.y+1 and (bflo[0] or bfhi[1]))
            {
                if (bflo[0] and bfhi[1])
                {
                    a(i,j,k) = a(i+1+offset, j-1-offset, k);
                }
                else if (bflo[0])
                {
                    a(i,j,k) = a(i+1+offset, j, k);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
            }
            // xhi & yhi
            else if (i == dhi.x+1 and j == dhi.y+1 and (bfhi[0] or bfhi[1]))
            {
                if (bfhi[0] and bfhi[1])
                {
                    a(i,j,k) = a(i-1-offset, j-1-offset, k);
                }
                else if (bfhi[0])
                {
                    a(i,j,k) = a(i-1-offset, j, k);
                }
                else if (bfhi[1])
                {
                    a(i,j,k) = a(i, j-1-offset, k);
                }
            }
            else if (i == dlo.x-1 and bflo[0])
            {
                a(i,j,k) = a(i+1+offset, j, k);
            }
            else if (i == dhi.x+1 and bfhi[0])
            {
                a(i,j,k) = a(i-1-offset, j, k);
            }
            else if (j == dlo.y-1 and bflo[1])
            {
                a(i,j,k) = a(i, j+1+offset, k);
            }
            else if (j == dhi.y+1 and bfhi[1])
            {
                a(i,j,k) = a(i, j-1-offset, k);
            }
        }
    });
}

//
// operator
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_ha (int i, int j, int k, Array4<Real> const& y, Array4<Real const> const& x,
                       Array4<Real const> const& sx, Array4<Real const> const& sy,
                       Array4<int const> const& msk, bool is_rz,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        Real facx = (1./6.)*dxinv[0]*dxinv[0];
        Real facy = (1./6.)*dxinv[1]*dxinv[1];
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
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_aa (int i, int j, int k, Array4<Real> const& y, Array4<Real const> const& x,
                       Array4<Real const> const& sig, Array4<int const> const& msk,
                       bool is_rz, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        y(i,j,k) = Real(0.0);
    } else {
        Real facx = (1.0/6.0)*dxinv[0]*dxinv[0];
        Real facy = (1.0/6.0)*dxinv[1]*dxinv[1];
        Real fxy = facx + facy;
        Real f2xmy = 2.0*facx - facy;
        Real fmx2y = 2.0*facy - facx;
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
    Real facx = -Real(2.0 * (1.0/6.0))*dxinv[0]*dxinv[0];
    Real facy = -Real(2.0 * (1.0/6.0))*dxinv[1]*dxinv[1];

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = Real(0.0);
        } else {
            sol(i,j,k) += Real(2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
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
    Real fac = -Real(2.0 * (1.0/6.0))*(dxinv[0]*dxinv[0] + dxinv[1]*dxinv[1]);

    amrex::LoopConcurrent(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = Real(0.0);
        } else {
            sol(i,j,k) += Real(2.0/3.0) * (rhs(i,j,k) - Ax(i,j,k))
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
    Real facx = Real(1.0/6.0)*dxinv[0]*dxinv[0];
    Real facy = Real(1.0/6.0)*dxinv[1]*dxinv[1];

    amrex::Loop(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = Real(0.0);
        } else {
            Real s0 = Real(-2.0)*(facx*(sx(i-1,j-1,k)+sx(i,j-1,k)+sx(i-1,j,k)+sx(i,j,k))
                                 +facy*(sy(i-1,j-1,k)+sy(i,j-1,k)+sy(i-1,j,k)+sy(i,j,k)));

            Real Ax = sol(i-1,j-1,k)*(facx*sx(i-1,j-1,k)+facy*sy(i-1,j-1,k))
                    + sol(i+1,j-1,k)*(facx*sx(i  ,j-1,k)+facy*sy(i  ,j-1,k))
                    + sol(i-1,j+1,k)*(facx*sx(i-1,j  ,k)+facy*sy(i-1,j  ,k))
                    + sol(i+1,j+1,k)*(facx*sx(i  ,j  ,k)+facy*sy(i  ,j  ,k))
                    + sol(i-1,j,k)*(Real(2.0)*facx*(sx(i-1,j-1,k)+sx(i-1,j,k))
                                        -     facy*(sy(i-1,j-1,k)+sx(i-1,j,k)))
                    + sol(i+1,j,k)*(Real(2.0)*facx*(sx(i  ,j-1,k)+sx(i  ,j,k))
                                        -     facy*(sy(i  ,j-1,k)+sx(i  ,j,k)))
                    + sol(i,j-1,k)*(   -facx*(sx(i-1,j-1,k)+sx(i,j-1,k))
                             +Real(2.0)*facy*(sy(i-1,j-1,k)+sy(i,j-1,k)))
                    + sol(i,j+1,k)*(   -facx*(sx(i-1,j  ,k)+sx(i,j  ,k))
                             +Real(2.0)*facy*(sy(i-1,j  ,k)+sy(i,j  ,k)))
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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void tridiagonal_solve (Array1D<Real,0,31>& a_ls, Array1D<Real,0,31>& b_ls, Array1D<Real,0,31>& c_ls,
                        Array1D<Real,0,31>& r_ls, Array1D<Real,0,31>& u_ls, Array1D<Real,0,31>& gam,
                        int ilen ) noexcept
{
    Real bet = b_ls(0);
    u_ls(0) = r_ls(0) / bet;

    for (int i = 1; i <= ilen - 1; i++) {
        gam(i) = c_ls(i-1) / bet;
        bet = b_ls(i) - a_ls(i)*gam(i);
        if (bet == 0) amrex::Abort(">>>TRIDIAG FAILED");
        u_ls(i) = (r_ls(i)-a_ls(i)*u_ls(i-1)) / bet;
    }
    for (int i = ilen-2; i >= 0; i--) {
        u_ls(i) = u_ls(i) - gam(i+1)*u_ls(i+1);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_gauss_seidel_with_line_solve_aa (Box const& bx, Array4<Real> const& sol,
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

    if (is_rz) {
	amrex::Abort("mlndlap_gauss_seidel_with_line_solve_aa is not implemented in r-z 2D ");
    } 

    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    int idir = -1;
    int ilen = 33;
    int k = 0;
    if (dxinv[0] <= dxinv[1])  {
        idir = 1;
        ilen = hi.y - lo.y + 1;
    }
    if (dxinv[1] <= dxinv[0]) {
        idir = 0;
        ilen = hi.x - lo.x + 1;
    }

    if (ilen > 32) amrex::Abort("mlndlap_gauss_seidel_with_line_solve_aa is hard-wired to be no longer than 32");

    Array1D<Real,0,31> a_ls,b_ls,c_ls,u_ls,r_ls,gam;

    if (idir == 1) {
        for (int i = lo.x; i <= hi.x; ++i)
        {
            for (int j = lo.y; j <= hi.y; ++j)
            {
                if (msk(i,j,k)) {
                    a_ls(j-lo.y) = 0.;
                    b_ls(j-lo.y) = 1.;
                    c_ls(j-lo.y) = 0.;
                    u_ls(j-lo.y) = 0.;
                    r_ls(j-lo.y) = 0.;
                }
                else
                {
                    Real s0 = (-2.0)*fxy*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));

                    Real Ax = sol(i-1,j-1,k)*fxy*sig(i-1,j-1,k)
                            + sol(i+1,j-1,k)*fxy*sig(i  ,j-1,k)
                            + sol(i-1,j+1,k)*fxy*sig(i-1,j  ,k)
                            + sol(i+1,j+1,k)*fxy*sig(i  ,j  ,k)
                            + sol(i-1,j,k)*f2xmy*(sig(i-1,j-1,k)+sig(i-1,j,k))
                            + sol(i+1,j,k)*f2xmy*(sig(i  ,j-1,k)+sig(i  ,j,k));

                    a_ls(j-lo.y) = fmx2y*(sig(i-1,j-1,k)+sig(i,j-1,k));
                    b_ls(j-lo.y) = s0;
                    c_ls(j-lo.y) = fmx2y*(sig(i-1,j  ,k)+sig(i,j  ,k));          
                    u_ls(j-lo.y) = 0.;
                    r_ls(j-lo.y) = rhs(i,j,k) - Ax;
                }
            }
            tridiagonal_solve(a_ls, b_ls, c_ls, r_ls, u_ls, gam, ilen);

            for (int j = lo.y; j <= hi.y; ++j)
            {
                sol(i,j,k) = u_ls(j-lo.y);
            }
        }
    } else if (idir == 0) {
        for (int j = lo.y ;j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                if (msk(i,j,k)) {
                    a_ls(i-lo.x) = 0.;
                    b_ls(i-lo.x) = 1.;
                    c_ls(i-lo.x) = 0.;
                    u_ls(i-lo.x) = 0.;
                    r_ls(i-lo.x) = 0.;
                }
                else
                {
                    Real s0 = (-2.0)*fxy*(sig(i-1,j-1,k)+sig(i,j-1,k)+sig(i-1,j,k)+sig(i,j,k));

                    Real Ax = sol(i-1,j-1,k)*fxy*sig(i-1,j-1,k)
                            + sol(i+1,j-1,k)*fxy*sig(i  ,j-1,k)
                            + sol(i-1,j+1,k)*fxy*sig(i-1,j  ,k)
                            + sol(i+1,j+1,k)*fxy*sig(i  ,j  ,k)
                            + sol(i,j-1,k)*fmx2y*(sig(i-1,j-1,k)+sig(i,j-1,k))
                            + sol(i,j+1,k)*fmx2y*(sig(i-1,j  ,k)+sig(i,j  ,k));

                    a_ls(i-lo.x) = f2xmy*(sig(i-1,j-1,k)+sig(i-1,j,k));
                    b_ls(i-lo.x) = s0;
                    c_ls(i-lo.x) = f2xmy*(sig(i  ,j-1,k)+sig(i  ,j,k));
                    u_ls(i-lo.x) = 0.;
                    r_ls(i-lo.x) = rhs(i,j,k) - Ax;

                }
            }
            tridiagonal_solve(a_ls, b_ls, c_ls, r_ls, u_ls, gam, ilen);

            for (int i = lo.x; i <= hi.x; ++i)
            {
                sol(i,j,k) = u_ls(i-lo.x);
            }
        }
    } else {
        amrex::Abort("mlndlap_gauss_seidel_with_line_solve_aa is wrong direction.");
    }

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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_restriction (int i, int j, int k, Array4<Real> const& crse,
                          Array4<Real const> const& fine, Array4<int const> const& msk, int idir) noexcept
{
    int kk = 0;
    if (idir == 1) {
        int ii = i*2;
        int jj = j;
        if (msk(ii,jj,kk)) {
            crse(i,j,k) = 0.0;
        } else {
            crse(i,j,k) = (1./4.)*(fine(ii-1,jj,kk) + 2.*fine(ii,jj,kk) + fine(ii+1,jj,kk));
        }
    } else if (idir == 0) {
        int ii = i;
        int jj = j*2;
        if (msk(ii,jj,kk)) {
            crse(i,j,k) = 0.0;
        } else {
            crse(i,j,k) = (1./4.)*(fine(ii,jj-1,kk) + 2.*fine(ii,jj,kk) + fine(ii,jj+1,kk));
        }
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
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real aa_interp_line_y (Array4<Real const> const& crse, Array4<Real const> const& sig,
                           int i, int j, int ic, int jc) noexcept
    {
        Real w1 = sig(i-1,j-1,0) + sig(i,j-1,0);
        Real w2 = sig(i-1,j  ,0) + sig(i,j  ,0);
        return (w1*crse(ic,jc,0)+w2*crse(ic,jc+1,0))/(w1+w2);
    }

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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_semi_interpadd_aa (int i, int j, int, Array4<Real> const& fine,
                           Array4<Real const> const& crse, Array4<Real const> const& sig,
                           Array4<int const> const& msk, int idir) noexcept
{
    if (idir == 1) {
        if (!msk(i,j,0)) {
            int ic = amrex::coarsen(i,2);
            int jc = j;
            bool i_is_odd = (ic*2 != i);
            if (i_is_odd) {
                // Node on X line
                fine(i,j,0) += aa_interp_line_x(crse,sig,i,j,ic,jc);
            } else {
                //Node coincident with coarse node
                fine(i,j,0) += crse(ic,jc,0);
            }
        }
    } else if (idir == 0 ) { 
        if (!msk(i,j,0)) { 
            int ic = i;
            int jc = amrex::coarsen(j,2);
            bool j_is_odd = (ic*2 != i); 
            if (j_is_odd) {
                // Node on Y line
                fine(i,j,0) += aa_interp_line_y(crse,sig,i,j,ic,jc);
            } else {
                //Node coincident with coarse node
                fine(i,j,0) += crse(ic,jc,0);
            }
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
                   Box const& nodal_domain,
                   GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                   GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi,
                   bool is_rz) noexcept
{
    Real facx = 0.5*dxinv[0];
    Real facy = 0.5*dxinv[1];

    const auto domlo = amrex::lbound(nodal_domain);
    const auto domhi = amrex::ubound(nodal_domain);

    if (msk(i,j,k)) {
        rhs(i,j,k) = 0.0;
    } else {

        Real zero_ilo = 1.0;
        Real zero_ihi = 1.0;
        Real zero_jlo = 1.0;
        Real zero_jhi = 1.0;

        // The nodal divergence operator should not see the tangential velocity
        //     at an inflow face
        if ((bclo[0] == LinOpBCType::Neumann or bclo[0] == LinOpBCType::inflow)
            and i == domlo.x) zero_ilo = 0.0;
        if ((bchi[0] == LinOpBCType::Neumann or bchi[0] == LinOpBCType::inflow)
            and i == domhi.x) zero_ihi = 0.0;
        if ((bclo[1] == LinOpBCType::Neumann or bclo[1] == LinOpBCType::inflow)
            and j == domlo.y) zero_jlo = 0.0;
        if ((bchi[1] == LinOpBCType::Neumann or bchi[1] == LinOpBCType::inflow)
            and j == domhi.y) zero_jhi = 0.0; 

        rhs(i,j,k) = facx*(-vel(i-1,j-1,k,0)*zero_jlo + vel(i,j-1,k,0)*zero_jlo
                           -vel(i-1,j  ,k,0)*zero_jhi + vel(i,j  ,k,0)*zero_jhi)
                   + facy*(-vel(i-1,j-1,k,1)*zero_ilo - vel(i,j-1,k,1)*zero_ihi
                           +vel(i-1,j  ,k,1)*zero_ilo + vel(i,j  ,k,1)*zero_ihi);
        if (is_rz) {
            // Here we assume we can't have inflow in the radial direction
            Real fm = facy / static_cast<Real>(6*i-3);
            Real fp = facy / static_cast<Real>(6*i+3);
            rhs(i,j,k) += fm*(vel(i-1,j,k,1)-vel(i-1,j-1,k,1))
                        - fp*(vel(i  ,j,k,1)-vel(i  ,j-1,k,1));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndlap_rhcc (int i, int j, int k, Array4<Real const> const& rhcc,
                   Array4<int const> const& msk) noexcept
{
    Real r;
    if (msk(i,j,k)) {
        r = 0.0;
    } else {
        r = 0.25 * (rhcc(i-1,j-1,k)+rhcc(i,j-1,k)+rhcc(i-1,j,k)+rhcc(i,j,k));
    }
    return r;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu (int i, int j, int k, Array4<Real> const& u, Array4<Real const> const& p,
                     Array4<Real const> const& sig,  GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                     bool is_rz) noexcept
{
    Real facx = 0.5*dxinv[0];
    Real facy = 0.5*dxinv[1];
    u(i,j,k,0) -= sig(i,j,k)*facx*(-p(i,j,k)+p(i+1,j,k)-p(i,j+1,k)+p(i+1,j+1,k));
    u(i,j,k,1) -= sig(i,j,k)*facy*(-p(i,j,k)-p(i+1,j,k)+p(i,j+1,k)+p(i+1,j+1,k));
    if (is_rz) {
        Real frz = sig(i,j,k)*facy / static_cast<Real>(6*i+3);
        u(i,j,k,1) += frz*(p(i,j,k)-p(i+1,j,k)-p(i,j+1,k)+p(i+1,j+1,k));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_compute_fine_contrib (int i, int j, int, Box const& fvbx,
                                        Array4<Real> const& frh, Array4<Real const> const& vel,
                                        GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                        Box const& nodal_domain,
                                        GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                                        GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi,
                                        bool is_rz) noexcept
{
    const auto domlo = amrex::lbound(nodal_domain);
    const auto domhi = amrex::ubound(nodal_domain);

    IntVect iv(i,j);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv))
    {
        Real zero_ilo = 1.0;
        Real zero_ihi = 1.0;
        Real zero_jlo = 1.0;
        Real zero_jhi = 1.0;

        // The nodal divergence operator should not see the tangential velocity
        //     at an inflow face
        if ((bclo[0] == LinOpBCType::Neumann or bclo[0] == LinOpBCType::inflow)
            and i == domlo.x) zero_ilo = 0.0;
        if ((bchi[0] == LinOpBCType::Neumann or bchi[0] == LinOpBCType::inflow)
            and i == domhi.x) zero_ihi = 0.0;
        if ((bclo[1] == LinOpBCType::Neumann or bclo[1] == LinOpBCType::inflow)
            and j == domlo.y) zero_jlo = 0.0;
        if ((bchi[1] == LinOpBCType::Neumann or bchi[1] == LinOpBCType::inflow)
            and j == domhi.y) zero_jhi = 0.0; 

        Real facx = Real(0.5)*dxinv[0];
        Real facy = Real(0.5)*dxinv[1];

        frh(i,j,0) = facx*(-vel(i-1,j-1,0,0)*zero_jlo+vel(i,j-1,0,0)*zero_jlo-vel(i-1,j,0,0)*zero_jhi+vel(i,j,0,0)*zero_jhi)
            +        facy*(-vel(i-1,j-1,0,1)*zero_ilo-vel(i,j-1,0,1)*zero_ihi+vel(i-1,j,0,1)*zero_ilo+vel(i,j,0,1)*zero_ihi);

        if (is_rz) {
            // Here we assume we can't have inflow in the radial direction
            Real fm = facy / static_cast<Real>(6*i-3);
            Real fp = facy / static_cast<Real>(6*i+3);
            frh(i,j,0) += fm*(vel(i-1,j,0,1)-vel(i-1,j-1,0,1))
                        - fp*(vel(i  ,j,0,1)-vel(i  ,j-1,0,1));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_add_fine_contrib (int i, int j, int /*k*/, Box const& fvbx,
                                    Array4<Real> const& rhs, Array4<Real const> const& frh,
                                    Array4<int const> const& msk) noexcept
{
    constexpr Real rfd = Real(0.25);
    constexpr Real chip = Real(0.5);
    constexpr Real chip2 = Real(0.25);

    int ii = 2*i;
    int jj = 2*j;
    IntVect iv(ii,jj);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv) and msk(ii,jj,0))
    {
        rhs(i,j,0) +=
            rfd*(frh(ii,jj,0)
                 + chip*(frh(ii-1,jj,0)+frh(ii+1,jj,0)+frh(ii,jj-1,0)+frh(ii,jj+1,0))
                 + chip2*(frh(ii-1,jj-1,0)+frh(ii+1,jj-1,0)+frh(ii-1,jj+1,0)+frh(ii+1,jj+1,0)));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_rhcc_fine_contrib (int i, int j, int, Box const& fvbx,
                                Array4<Real> const& rhs, Array4<Real const> const& cc,
                                Array4<int const> const& msk) noexcept
{
    constexpr Real w1 = Real(9./64.);
    constexpr Real w2 = Real(3./64.);
    constexpr Real w3 = Real(1./64.);

    int ii = 2*i;
    int jj = 2*j;
    IntVect iv(ii,jj);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv) and msk(ii,jj,0))
    {
        rhs(i,j,0) += w1*(cc(ii-1,jj-1,0)+cc(ii  ,jj-1,0)+cc(ii-1,jj  ,0)+cc(ii  ,jj  ,0))
                    + w2*(cc(ii-2,jj-1,0)+cc(ii+1,jj-1,0)+cc(ii-2,jj  ,0)+cc(ii+1,jj  ,0)
                         +cc(ii-1,jj-2,0)+cc(ii  ,jj-2,0)+cc(ii-1,jj+1,0)+cc(ii  ,jj+1,0))
                    + w3*(cc(ii-2,jj-2,0)+cc(ii+1,jj-2,0)+cc(ii-2,jj+1,0)+cc(ii+1,jj+1,0));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_cf_contrib (int i, int j, int, Array4<Real> const& rhs,
                              Array4<Real const> const& vel, Array4<Real const> const& fc,
                              Array4<Real const> const& rhcc, Array4<int const> const& dmsk,
                              Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                              bool is_rz, GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Box const& nddom, GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bclo,
                              GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bchi,
                              bool neumann_doubling) noexcept
{
    if (!dmsk(i,j,0) and ndmsk(i,j,0) == crse_fine_node) {
        Real facx = Real(0.5) * dxinv[0];
        Real facy = Real(0.5) * dxinv[1];
        Real r = fc(i,j,0);
        if (rhcc) {
            r += Real(0.25)*( (Real(1.)-ccmsk(i-1,j-1,0)) * rhcc(i-1,j-1,0)
                            + (Real(1.)-ccmsk(i  ,j-1,0)) * rhcc(i  ,j-1,0)
                            + (Real(1.)-ccmsk(i-1,j  ,0)) * rhcc(i-1,j  ,0)
                            + (Real(1.)-ccmsk(i  ,j  ,0)) * rhcc(i  ,j  ,0));
        }
        r += (Real(1.)-ccmsk(i-1,j-1,0)) * (-facx*vel(i-1,j-1,0,0) - facy*vel(i-1,j-1,0,1))
           + (Real(1.)-ccmsk(i  ,j-1,0)) * ( facx*vel(i  ,j-1,0,0) - facy*vel(i  ,j-1,0,1))
           + (Real(1.)-ccmsk(i-1,j  ,0)) * (-facx*vel(i-1,j  ,0,0) + facy*vel(i-1,j  ,0,1))
           + (Real(1.)-ccmsk(i  ,j  ,0)) * ( facx*vel(i  ,j  ,0,0) + facy*vel(i  ,j  ,0,1));
        if (is_rz) {
            Real fm = facy / static_cast<Real>(6*i-3);
            Real fp = facy / static_cast<Real>(6*i+3);
            r += fm*((Real(1.)-ccmsk(i-1,j  ,0))*vel(i-1,j  ,0,1)
                    -(Real(1.)-ccmsk(i-1,j-1,0))*vel(i-1,j-1,0,1))
               - fp*((Real(1.)-ccmsk(i  ,j  ,0))*vel(i  ,j  ,0,1)
                    -(Real(1.)-ccmsk(i  ,j-1,0))*vel(i  ,j-1,0,1));
        }

        if (neumann_doubling) {
            const auto ndlo = amrex::lbound(nddom);
            const auto ndhi = amrex::ubound(nddom);
            if (i == ndlo.x and ( bclo[0] == LinOpBCType::Neumann or
                                  bclo[0] == LinOpBCType::inflow)) {
                r *= Real(2.);
            } else if (i== ndhi.x and ( bchi[0] == LinOpBCType::Neumann or
                                        bchi[0] == LinOpBCType::inflow)) {
                r *= Real(2.);
            }

            if (j == ndlo.y and ( bclo[1] == LinOpBCType::Neumann or
                                  bclo[1] == LinOpBCType::inflow)) {
                r *= Real(2.);
            } else if (j == ndhi.y and ( bchi[1] == LinOpBCType::Neumann or
                                         bchi[1] == LinOpBCType::inflow)) {
                r *= Real(2.);
            }
        }

        rhs(i,j,0) = r;
    }
}

//
// residual
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_crse_resid (int i, int j, int k, Array4<Real> const& resid,
                         Array4<Real const> const& rhs, Array4<int const> const& msk,
                         Box const& nddom, GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bclo,
                         GpuArray<LinOpBCType,AMREX_SPACEDIM> const& bchi,
                         bool neumann_doubling) noexcept
{
    if ((msk(i-1,j-1,k  ) == 0 or
         msk(i  ,j-1,k  ) == 0 or
         msk(i-1,j  ,k  ) == 0 or
         msk(i  ,j  ,k  ) == 0) and
        (msk(i-1,j-1,k  ) == 0 or
         msk(i  ,j-1,k  ) == 0 or
         msk(i-1,j  ,k  ) == 0 or
         msk(i  ,j  ,k  ) == 0))
    {
        Real fac = Real(1.0);
        if (neumann_doubling) {
            const auto ndlo = amrex::lbound(nddom);
            const auto ndhi = amrex::ubound(nddom);
            if (i == ndlo.x and ( bclo[0] == LinOpBCType::Neumann or
                                  bclo[0] == LinOpBCType::inflow)) {
                fac *= Real(2.);
            } else if (i== ndhi.x and ( bchi[0] == LinOpBCType::Neumann or
                                        bchi[0] == LinOpBCType::inflow)) {
                fac *= Real(2.);
            }
            if (j == ndlo.y and ( bclo[1] == LinOpBCType::Neumann or
                                  bclo[1] == LinOpBCType::inflow)) {
                fac *= Real(2.);
            } else if (j == ndhi.y and ( bchi[1] == LinOpBCType::Neumann or
                                         bchi[1] == LinOpBCType::inflow)) {
                fac *= Real(2.);
            }
        }

        resid(i,j,k) = (rhs(i,j,k) - resid(i,j,k)) * fac;
    } else {
        resid(i,j,k) = Real(0.);
    }
}

//
// sync residual
//

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_fine_Ax (int i, int j, int, Box const& fvbx, Array4<Real> const& Ax,
                          Array4<Real const> const& x, Array4<Real const> const& sig,
                          bool is_rz,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    IntVect iv(i,j);
    if (fvbx.contains(iv) and !fvbx.strictly_contains(iv)) {
        Real facx = Real(1./6.)*dxinv[0]*dxinv[0];
        Real facy = Real(1./6.)*dxinv[1]*dxinv[1];
        Real fxy = facx + facy;
        Real f2xmy = Real(2.)*facx - facy;
        Real fmx2y = Real(2.)*facy - facx;
        Ax(i,j,0) = x(i-1,j-1,0)*fxy*sig(i-1,j-1,0)
            +       x(i+1,j-1,0)*fxy*sig(i  ,j-1,0)
            +       x(i-1,j+1,0)*fxy*sig(i-1,j  ,0)
            +       x(i+1,j+1,0)*fxy*sig(i  ,j  ,0)
            +       x(i-1,j,0)*f2xmy*(sig(i-1,j-1,0)+sig(i-1,j  ,0))
            +       x(i+1,j,0)*f2xmy*(sig(i  ,j-1,0)+sig(i  ,j  ,0))
            +       x(i,j-1,0)*fmx2y*(sig(i-1,j-1,0)+sig(i  ,j-1,0))
            +       x(i,j+1,0)*fmx2y*(sig(i-1,j  ,0)+sig(i  ,j  ,0))
            +       x(i,j,0)*Real(-2.)*fxy*(sig(i-1,j-1,0)+sig(i,j-1,0)
                                           +sig(i-1,j,0)+sig(i,j,0));
        if (is_rz) {
            Real fp = facy / static_cast<Real>(2*i+1);
            Real fm = facy / static_cast<Real>(2*i-1);
            Ax(i,j,0) += (fm*sig(i-1,j  ,0)-fp*sig(i,j  ,0))*(x(i,j+1,0)-x(i,j,0))
                       + (fm*sig(i-1,j-1,0)-fp*sig(i,j-1,0))*(x(i,j-1,0)-x(i,j,0));
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_fine_contrib (int i, int j, int k, Array4<Real> const& f,
                               Array4<Real const> const& Ax,
                               Array4<int const> const& msk) noexcept
{
    constexpr Real rfd = Real(0.25);
    constexpr Real chip = Real(0.5);
    constexpr Real chip2 = Real(0.25);

    int ii = 2*i;
    int jj = 2*j;
    int kk = 2*k;
    if (msk(ii,jj,kk)) {
        f(i,j,0) += rfd*(Ax(ii,jj,0)
                 + chip*(Ax(ii-1,jj,0)+Ax(ii+1,jj,0)+Ax(ii,jj-1,0)+Ax(ii,jj+1,0))
                 + chip2*(Ax(ii-1,jj-1,0)+Ax(ii+1,jj-1,0)+Ax(ii-1,jj+1,0)+Ax(ii+1,jj+1,0)));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_res_cf_contrib (int i, int j, int, Array4<Real> const& res,
                             Array4<Real const> const& phi, Array4<Real const> const& rhs,
                             Array4<Real const> const& sig, Array4<int const> const& dmsk,
                             Array4<int const> const& ndmsk, Array4<int const> const& ccmsk,
                             Array4<Real const> const& fc,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv, Box const& nddom,
                             bool is_rz,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                             GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi,
                             bool neumann_doubling) noexcept
{
    if (!dmsk(i,j,0) and ndmsk(i,j,0) == crse_fine_node) {
        Real facx = Real(1./6.)*dxinv[0]*dxinv[0];
        Real facy = Real(1./6.)*dxinv[1]*dxinv[1];

        Real fp, fm;
        if (is_rz) {
            fp = facy / static_cast<Real>(2*i+1);
            fm = facy / static_cast<Real>(2*i-1);
        }

        Real Ax = Real(0.);
        if (ccmsk(i-1,j-1,0) == crse_cell) {
            Ax += sig(i-1,j-1,0)*(facx*(Real(2.)*(phi(i-1,j  ,0)-phi(i  ,j  ,0))
                                           +     (phi(i-1,j-1,0)-phi(i  ,j-1,0)))
                                + facy*(Real(2.)*(phi(i  ,j-1,0)-phi(i  ,j  ,0))
                                           +     (phi(i-1,j-1,0)-phi(i-1,j  ,0))));
            if (is_rz) {
                Ax += fm*sig(i-1,j-1,0)*(phi(i,j-1,0)-phi(i,j,0));
            }
        }
        if (ccmsk(i,j-1,0) == crse_cell) {
            Ax += sig(i,j-1,0)*(facx*(Real(2.)*(phi(i+1,j  ,0)-phi(i  ,j  ,0))
                                         +     (phi(i+1,j-1,0)-phi(i  ,j-1,0)))
                              + facy*(Real(2.)*(phi(i  ,j-1,0)-phi(i  ,j  ,0))
                                         +     (phi(i+1,j-1,0)-phi(i+1,j  ,0))));
            if (is_rz) {
                Ax -= fp*sig(i,j-1,0)*(phi(i,j-1,0)-phi(i,j,0));
            }
        }
        if (ccmsk(i-1,j,0) == crse_cell) {
            Ax += sig(i-1,j,0)*(facx*(Real(2.)*(phi(i-1,j  ,0)-phi(i  ,j  ,0))
                                         +     (phi(i-1,j+1,0)-phi(i  ,j+1,0)))
                              + facy*(Real(2.)*(phi(i  ,j+1,0)-phi(i  ,j  ,0))
                                         +     (phi(i-1,j+1,0)-phi(i-1,j  ,0))));
            if (is_rz) {
                Ax += fm*sig(i-1,j,0)*(phi(i,j+1,0)-phi(i,j,0));
            }
        }
        if (ccmsk(i,j,0) == crse_cell) {
            Ax += sig(i,j,0)*(facx*(Real(2.)*(phi(i+1,j  ,0)-phi(i  ,j  ,0))
                                      +      (phi(i+1,j+1,0)-phi(i  ,j+1,0)))
                            + facy*(Real(2.)*(phi(i  ,j+1,0)-phi(i  ,j  ,0))
                                      +      (phi(i+1,j+1,0)-phi(i+1,j  ,0))));
            if (is_rz) {
                Ax -= fp*sig(i,j,0)*(phi(i,j+1,0)-phi(i,j,0));
            }
        }

        Real Axf = fc(i,j,0);
        if (neumann_doubling) {
            const auto ndlo = amrex::lbound(nddom);
            const auto ndhi = amrex::ubound(nddom);

            if (i == ndlo.x and (bclo[0] == LinOpBCType::Neumann or
                                 bclo[0] == LinOpBCType::inflow)) {
                Axf *= Real(2.);
            } else if (i== ndhi.x and (bchi[0] == LinOpBCType::Neumann or
                                       bchi[0] == LinOpBCType::inflow)) {
                Axf *= Real(2.);
            }

            if (j == ndlo.y and (bclo[1] == LinOpBCType::Neumann or
                                 bclo[1] == LinOpBCType::inflow)) {
                Axf *= Real(2.);
            } else if (j == ndhi.y and (bchi[1] == LinOpBCType::Neumann or
                                        bchi[1] == LinOpBCType::inflow)) {
                Axf *= Real(2.);
            }
        }

        res(i,j,0) = rhs(i,j,0) - (Ax + Axf);
    }
}

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
    sten(i,j,k,4) = 1.0 / (amrex::Math::abs(sten(i-1,j  ,k,1)) + amrex::Math::abs(sten(i,j  ,k,1))
                         + amrex::Math::abs(sten(i  ,j-1,k,2)) + amrex::Math::abs(sten(i,j  ,k,2))
                         + amrex::Math::abs(sten(i-1,j-1,k,3)) + amrex::Math::abs(sten(i,j-1,k,3))
                         + amrex::Math::abs(sten(i-1,j  ,k,3)) + amrex::Math::abs(sten(i,j  ,k,3)) + eps);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_stencil_rap (int i, int j, int, Array4<Real> const& csten,
                          Array4<Real const> const& fsten) noexcept
{
    constexpr int k = 0;

#if 0
    auto interp_from_mm_to = [&fsten] (int i_, int j_) -> Real {
        Real wxm = amrex::Math::abs(fsten(i_-1,j_  ,0,1))/(amrex::Math::abs(fsten(i_-1,j_-1,0,3))+amrex::Math::abs(fsten(i_-1,j_  ,0,3))+eps);
        Real wym = amrex::Math::abs(fsten(i_  ,j_-1,0,2))/(amrex::Math::abs(fsten(i_-1,j_-1,0,3))+amrex::Math::abs(fsten(i_  ,j_-1,0,3))+eps);
        Real wmm = amrex::Math::abs(fsten(i_-1,j_-1,0,3)) * (Real(1.) + wxm + wym);
        return wmm * fsten(i_,j_,0,4);
    };
#endif

    auto interp_from_mp_to = [&fsten] (int i_, int j_) -> Real {
        Real wxm = amrex::Math::abs(fsten(i_-1,j_  ,0,1))/(amrex::Math::abs(fsten(i_-1,j_-1,0,3))+amrex::Math::abs(fsten(i_-1,j_  ,0,3))+eps);
        Real wyp = amrex::Math::abs(fsten(i_  ,j_  ,0,2))/(amrex::Math::abs(fsten(i_-1,j_  ,0,3))+amrex::Math::abs(fsten(i_  ,j_  ,0,3))+eps);
        Real wmp = amrex::Math::abs(fsten(i_-1,j_  ,0,3)) *(Real(1.) + wxm + wyp);
        return wmp * fsten(i_,j_,0,4);
    };

    auto interp_from_pm_to = [&fsten] (int i_, int j_) -> Real {
        Real wxp = amrex::Math::abs(fsten(i_  ,j_  ,0,1))/(amrex::Math::abs(fsten(i_  ,j_-1,0,3))+amrex::Math::abs(fsten(i_  ,j_  ,0,3))+eps);
        Real wym = amrex::Math::abs(fsten(i_  ,j_-1,0,2))/(amrex::Math::abs(fsten(i_-1,j_-1,0,3))+amrex::Math::abs(fsten(i_  ,j_-1,0,3))+eps);
        Real wpm = amrex::Math::abs(fsten(i_  ,j_-1,0,3)) * (Real(1.) + wxp + wym);
        return wpm * fsten(i_,j_,0,4);
    };

    auto interp_from_pp_to = [&fsten] (int i_, int j_) -> Real {
        Real wxp = amrex::Math::abs(fsten(i_  ,j_  ,0,1))/(amrex::Math::abs(fsten(i_  ,j_-1,0,3))+amrex::Math::abs(fsten(i_  ,j_  ,0,3))+eps);
        Real wyp = amrex::Math::abs(fsten(i_  ,j_  ,0,2))/(amrex::Math::abs(fsten(i_-1,j_  ,0,3))+amrex::Math::abs(fsten(i_  ,j_  ,0,3))+eps);
        Real wpp = amrex::Math::abs(fsten(i_  ,j_  ,0,3)) * (Real(1.) + wxp + wyp);
        return wpp * fsten(i_,j_,0,4);
    };

    auto interp_from_m0_to = [&fsten] (int i_, int j_) -> Real {
        return amrex::Math::abs(fsten(i_-1,j_,0,1))/(amrex::Math::abs(fsten(i_-1,j_,0,1))+amrex::Math::abs(fsten(i_,j_,0,1))+eps);
    };

    auto interp_from_p0_to = [&fsten] (int i_, int j_) -> Real {
        return amrex::Math::abs(fsten(i_,j_,0,1))/(amrex::Math::abs(fsten(i_-1,j_,0,1))+amrex::Math::abs(fsten(i_,j_,0,1))+eps);
    };

    auto interp_from_0m_to = [&fsten] (int i_, int j_) -> Real {
        return amrex::Math::abs(fsten(i_,j_-1,0,2))/(amrex::Math::abs(fsten(i_,j_-1,0,2))+amrex::Math::abs(fsten(i_,j_,0,2))+eps);
    };

    auto interp_from_0p_to = [&fsten] (int i_, int j_) -> Real {
        return amrex::Math::abs(fsten(i_,j_,0,2))/(amrex::Math::abs(fsten(i_,j_-1,0,2))+amrex::Math::abs(fsten(i_,j_,0,2))+eps);
    };

#if 0
    auto Amm = [&fsten] (int i_, int j_) -> Real {
        return fsten(i_-1,j_-1,0,3);
    };
#endif

    auto A0m = [&fsten] (int i_, int j_) -> Real {
        return fsten(i_,j_-1,0,2);
    };

    auto Apm = [&fsten] (int i_, int j_) -> Real {
        return fsten(i_,j_-1,0,3);
    };

    auto Am0 = [&fsten] (int i_, int j_) -> Real {
        return fsten(i_-1,j_,0,1);
    };

    auto A00 = [&fsten] (int i_, int j_) -> Real {
        return fsten(i_,j_,0,0);
    };

    auto Ap0 = [&fsten] (int i_, int j_) -> Real {
        return fsten(i_,j_,0,1);
    };

    auto Amp = [&fsten] (int i_, int j_) -> Real {
        return fsten(i_-1,j_,0,3);
    };

    auto A0p = [&fsten] (int i_, int j_) -> Real {
        return fsten(i_,j_,0,2);
    };

    auto App = [&fsten] (int i_, int j_) -> Real {
        return fsten(i_,j_,0,3);
    };

#if 0
    auto restrict_from_mm_to = [&fsten] (int ii_, int jj_) -> Real {
        Real wxp = amrex::Math::abs(fsten(ii_-1,jj_-1,0,1))/(amrex::Math::abs(fsten(ii_-1,jj_-2,0,3))+amrex::Math::abs(fsten(ii_-1,jj_-1,0,3))+eps);
        Real wyp = amrex::Math::abs(fsten(ii_-1,jj_-1,0,2))/(amrex::Math::abs(fsten(ii_-2,jj_-1,0,3))+amrex::Math::abs(fsten(ii_-1,jj_-1,0,3))+eps);
        Real wpp = amrex::Math::abs(fsten(ii_-1,jj_-1,0,3))*(Real(1.)+wxp+wyp);
        return wpp * fsten(ii_-1,jj_-1,0,4);
    };
#endif

    auto restrict_from_0m_to = [&fsten] (int ii_, int jj_) -> Real {
        return amrex::Math::abs(fsten(ii_,jj_-1,0,2))/(amrex::Math::abs(fsten(ii_,jj_-2,0,2))+amrex::Math::abs(fsten(ii_,jj_-1,0,2))+eps);
    };

    auto restrict_from_pm_to = [&fsten] (int ii_, int jj_) -> Real {
        Real wxm = amrex::Math::abs(fsten(ii_  ,jj_-1,0,1))/(amrex::Math::abs(fsten(ii_,jj_-2,0,3))+amrex::Math::abs(fsten(ii_  ,jj_-1,0,3))+eps);
        Real wyp = amrex::Math::abs(fsten(ii_+1,jj_-1,0,2))/(amrex::Math::abs(fsten(ii_,jj_-1,0,3))+amrex::Math::abs(fsten(ii_+1,jj_-1,0,3))+eps);
        Real wmp = amrex::Math::abs(fsten(ii_  ,jj_-1,0,3)) *(Real(1.) + wxm + wyp);
        return wmp * fsten(ii_+1,jj_-1,0,4);
    };

    auto restrict_from_m0_to = [&fsten] (int ii_, int jj_) -> Real {
        return amrex::Math::abs(fsten(ii_-1,jj_,0,1))/(amrex::Math::abs(fsten(ii_-2,jj_,0,1))+amrex::Math::abs(fsten(ii_-1,jj_,0,1))+eps);
    };

    auto restrict_from_p0_to = [&fsten] (int ii_, int jj_) -> Real {
        return amrex::Math::abs(fsten(ii_,jj_,0,1))/(amrex::Math::abs(fsten(ii_,jj_,0,1))+amrex::Math::abs(fsten(ii_+1,jj_,0,1))+eps);
    };

    auto restrict_from_mp_to = [&fsten] (int ii_, int jj_) -> Real {
        Real wxp = amrex::Math::abs(fsten(ii_-1,jj_+1,0,1))/(amrex::Math::abs(fsten(ii_-1,jj_,0,3))+amrex::Math::abs(fsten(ii_-1,jj_+1,0,3))+eps);
        Real wym = amrex::Math::abs(fsten(ii_-1,jj_  ,0,2))/(amrex::Math::abs(fsten(ii_-2,jj_,0,3))+amrex::Math::abs(fsten(ii_-1,jj_  ,0,3))+eps);
        Real wpm = amrex::Math::abs(fsten(ii_-1,jj_  ,0,3)) * (Real(1.) + wxp + wym);
        return wpm * fsten(ii_-1,jj_+1,0,4);
    };

    auto restrict_from_0p_to = [&fsten] (int ii_, int jj_) -> Real {
        return amrex::Math::abs(fsten(ii_,jj_,0,2))/(amrex::Math::abs(fsten(ii_,jj_,0,2))+amrex::Math::abs(fsten(ii_,jj_+1,0,2))+eps);
    };

    auto restrict_from_pp_to = [&fsten] (int ii_, int jj_) -> Real {
        Real wxm = amrex::Math::abs(fsten(ii_  ,jj_+1,0,1))/(amrex::Math::abs(fsten(ii_  ,jj_  ,0,3))+amrex::Math::abs(fsten(ii_  ,jj_+1,0,3))+eps);
        Real wym = amrex::Math::abs(fsten(ii_+1,jj_  ,0,2))/(amrex::Math::abs(fsten(ii_  ,jj_  ,0,3))+amrex::Math::abs(fsten(ii_+1,jj_  ,0,3))+eps);
        Real wmm = amrex::Math::abs(fsten(ii_  ,jj_  ,0,3)) * (Real(1.) + wxm + wym);
        return wmm * fsten(ii_+1,jj_+1,0,4);
    };

    int ii = 2*i;
    int jj = 2*j;
    Array2D<Real,-1,1,-1,1> ap, p;

    // csten(i,j,k,1)
    p(-1,-1) = interp_from_pp_to(ii+1,jj-1);
    p( 0,-1) = interp_from_0p_to(ii+2,jj-1);
    p(-1, 0) = interp_from_p0_to(ii+1,jj  );
    p( 0, 0) = Real(1.);
    p(-1, 1) = interp_from_pm_to(ii+1,jj+1);
    p( 0, 1) = interp_from_0m_to(ii+2,jj+1);

    ap(0,-1) = Ap0(ii,jj-1)*p(-1,-1) + App(ii,jj-1)*p(-1,0);
    ap(1,-1) = A00(ii+1,jj-1)*p(-1,-1) + Ap0(ii+1,jj-1)*p(0,-1)
        +      A0p(ii+1,jj-1)*p(-1,0) + App(ii+1,jj-1)*p(0,0);
    ap(0,0) = Apm(ii,jj)*p(-1,-1) + Ap0(ii,jj)*p(-1,0) + App(ii,jj)*p(-1,1);
    ap(1,0) = A0m(ii+1,jj)*p(-1,-1) + Apm(ii+1,jj)*p(0,-1)
        +     A00(ii+1,jj)*p(-1,0) + Ap0(ii+1,jj)*p(0,0)
        +     A0p(ii+1,jj)*p(-1,1) + App(ii+1,jj)*p(0,1);
    ap(0,1) = Apm(ii,jj+1)*p(-1,0) + Ap0(ii,jj+1)*p(-1,1);
    ap(1,1) = A0m(ii+1,jj+1)*p(-1,0) + Apm(ii+1,jj+1)*p(0,0)
        +     A00(ii+1,jj+1)*p(-1,1) + Ap0(ii+1,jj+1)*p(0,1);

    csten(i,j,k,1) = Real(0.25)*(restrict_from_0m_to(ii,jj)*ap(0,-1)
                            + restrict_from_pm_to(ii,jj)*ap(1,-1)
                            + ap(0,0)
                            + restrict_from_p0_to(ii,jj)*ap(1,0)
                            + restrict_from_0p_to(ii,jj)*ap(0,1)
                            + restrict_from_pp_to(ii,jj)*ap(1,1));

    // csten(i,j,k,2)
    p(-1,-1) = interp_from_pp_to(ii-1,jj+1);
    p( 0,-1) = interp_from_0p_to(ii  ,jj+1);
    p( 1,-1) = interp_from_mp_to(ii+1,jj+1);
    p(-1, 0) = interp_from_p0_to(ii-1,jj+2);
    p( 0, 0) = Real(1.);
    p( 1, 0) = interp_from_m0_to(ii+1,jj+2);

    ap(-1,0) = A0p(ii-1,jj)*p(-1,-1) + App(ii-1,jj)*p(0,-1);
    ap(0,0) = Amp(ii,jj)*p(-1,-1) + A0p(ii,jj)*p(0,-1) + App(ii,jj)*p(1,-1);
    ap(1,0) = Amp(ii+1,jj)*p(0,-1) + A0p(ii+1,jj)*p(1,-1);
    ap(-1,1) = A00(ii-1,jj+1)*p(-1,-1) + Ap0(ii-1,jj+1)*p(0,-1)
        +      A0p(ii-1,jj+1)*p(-1,0) + App(ii-1,jj+1)*p(0,0);
    ap(0,1) = Am0(ii,jj+1)*p(-1,-1) + A00(ii,jj+1)*p(0,-1) + Ap0(ii,jj+1)*p(1,-1)
        +     Amp(ii,jj+1)*p(-1,0) + A0p(ii,jj+1)*p(0,0) + App(ii,jj+1)*p(1,0);
    ap(1,1) = Am0(ii+1,jj+1)*p(0,-1) + A00(ii+1,jj+1)*p(1,-1)
        +     Amp(ii+1,jj+1)*p(0,0) + A0p(ii+1,jj+1)*p(1,0);

    csten(i,j,k,2) = Real(0.25)*(restrict_from_m0_to(ii,jj)*ap(-1,0)
                            + ap(0,0)
                            + restrict_from_p0_to(ii,jj)*ap(1,0)
                            + restrict_from_mp_to(ii,jj)*ap(-1,1)
                            + restrict_from_0p_to(ii,jj)*ap(0,1)
                            + restrict_from_pp_to(ii,jj)*ap(1,1));

    // csten(i,j,k,3)
    p(-1,-1) = interp_from_pp_to(ii+1,jj+1);
    p( 0,-1) = interp_from_0p_to(ii+2,jj+1);
    p(-1, 0) = interp_from_p0_to(ii+1,jj+2);
    p( 0, 0) = Real(1.);

    ap(0,0) = App(ii,jj)*p(-1,-1);
    ap(1,0) = A0p(ii+1,jj)*p(-1,-1) + App(ii+1,jj)*p(0,-1);
    ap(0,1) = Ap0(ii,jj+1)*p(-1,-1) + App(ii,jj+1)*p(-1,0);
    ap(1,1) = A00(ii+1,jj+1)*p(-1,-1) + Ap0(ii+1,jj+1)*p(0,-1)
        +     A0p(ii+1,jj+1)*p(-1,0) + App(ii+1,jj+1)*p(0,0);

    Real cross1 = Real(0.25)*(ap(0,0)
                           + restrict_from_p0_to(ii,jj)*ap(1,0)
                           + restrict_from_0p_to(ii,jj)*ap(0,1)
                           + restrict_from_pp_to(ii,jj)*ap(1,1));

    p(0,-1) = interp_from_0p_to(ii,jj+1);
    p(1,-1) = interp_from_mp_to(ii+1,jj+1);
    p(0, 0) = Real(1.);
    p(1, 0) = interp_from_m0_to(ii+1,jj+2);

    ap(-1,0) = Amp(ii+1,jj)*p(0,-1) + A0p(ii+1,jj)*p(1,-1);
    ap( 0,0) = Amp(ii+2,jj)*p(1,-1);
    ap(-1,1) = Am0(ii+1,jj+1)*p(0,-1) + A00(ii+1,jj+1)*p(1,-1) + Amp(ii+1,jj+1)*p(0,0)
        + A0p(ii+1,jj+1)*p(1,0);
    ap( 0,1) = Am0(ii+2,jj+1)*p(1,-1) + Amp(ii+2,jj+1)*p(1,0);

    Real cross2 = Real(0.25)*(ap(0,0)
                           + restrict_from_m0_to(ii+2,jj)*ap(-1,0)
                           + restrict_from_mp_to(ii+2,jj)*ap(-1,1)
                           + restrict_from_0p_to(ii+2,jj)*ap( 0,1));

    csten(i,j,k,3) = Real(0.5)*(cross1+cross2);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_adotx_sten (int i, int j, int k, Array4<Real> const& y, Array4<Real const> const& x,
                         Array4<Real const> const& sten, Array4<int const> const& msk) noexcept
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
            Real wym = amrex::Math::abs(sten(i,j-1,0,2));
            Real wyp = amrex::Math::abs(sten(i,j  ,0,2));
            fv = (wym*crse(ic,jc,0) + wyp*crse(ic,jc+1,0)) / (wym+wyp+eps);
        } else if (jeven) {
            Real wxm = amrex::Math::abs(sten(i-1,j,0,1));
            Real wxp = amrex::Math::abs(sten(i  ,j,0,1));
            fv = (wxm*crse(ic,jc,0) + wxp*crse(ic+1,jc,0)) / (wxm+wxp+eps);
        } else {
            Real wxm = amrex::Math::abs(sten(i-1,j  ,0,1)) /
                (amrex::Math::abs(sten(i-1,j-1,0,3))+amrex::Math::abs(sten(i-1,j  ,0,3))+eps);
            Real wxp = amrex::Math::abs(sten(i  ,j  ,0,1)) /
                (amrex::Math::abs(sten(i  ,j-1,0,3))+amrex::Math::abs(sten(i  ,j  ,0,3))+eps);
            Real wym = amrex::Math::abs(sten(i  ,j-1,0,2)) /
                (amrex::Math::abs(sten(i-1,j-1,0,3))+amrex::Math::abs(sten(i  ,j-1,0,3))+eps);
            Real wyp = amrex::Math::abs(sten(i  ,j  ,0,2)) /
                (amrex::Math::abs(sten(i-1,j  ,0,3))+amrex::Math::abs(sten(i  ,j  ,0,3))+eps);
            Real wmm = amrex::Math::abs(sten(i-1,j-1,0,3)) * (1.0 + wxm + wym);
            Real wpm = amrex::Math::abs(sten(i,j-1,0,3)) * (1.0 + wxp + wym);
            Real wmp = amrex::Math::abs(sten(i-1,j,0,3)) *(1.0 + wxm + wyp);
            Real wpp = amrex::Math::abs(sten(i,j,0,3)) * (1.0 + wxp + wyp);
            fv = (wmm*crse(ic,jc,0) + wpm*crse(ic+1,jc,0)
                  + wmp*crse(ic,jc+1,0) + wpp*crse(ic+1,jc+1,0))
                / (wmm+wpm+wmp+wpp+eps);
        }

        fine(i,j,0) += fv;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_restriction_rap (int i, int j, int /*k*/, Array4<Real> const& crse,
                              Array4<Real const> const& fine, Array4<Real const> const& sten,
                              Array4<int const> const& msk) noexcept
{
    int ii = i*2;
    int jj = j*2;
    if (msk(ii,jj,0)) {
        crse(i,j,0) = 0.0;
    } else {

        Real cv = fine(ii,jj,0)
            + fine(ii-1,jj  ,0)*amrex::Math::abs(sten(ii-1,jj  ,0,1))
            /                  (amrex::Math::abs(sten(ii-2,jj  ,0,1))
                               +amrex::Math::abs(sten(ii-1,jj  ,0,1))+eps)
            + fine(ii+1,jj  ,0)*amrex::Math::abs(sten(ii  ,jj  ,0,1))
            /                  (amrex::Math::abs(sten(ii  ,jj  ,0,1))
                               +amrex::Math::abs(sten(ii+1,jj  ,0,1))+eps)
            + fine(ii  ,jj-1,0)*amrex::Math::abs(sten(ii  ,jj-1,0,2))
            /                  (amrex::Math::abs(sten(ii  ,jj-2,0,2))
                               +amrex::Math::abs(sten(ii  ,jj-1,0,2))+eps)
            + fine(ii  ,jj+1,0)*amrex::Math::abs(sten(ii  ,jj  ,0,2))
            /                  (amrex::Math::abs(sten(ii  ,jj  ,0,2))
                               +amrex::Math::abs(sten(ii  ,jj+1,0,2))+eps);

        Real wxp = amrex::Math::abs(sten(ii-1,jj-1,0,1))
            /     (amrex::Math::abs(sten(ii-1,jj-2,0,3))
                  +amrex::Math::abs(sten(ii-1,jj-1,0,3))+eps);
        Real wyp = amrex::Math::abs(sten(ii-1,jj-1,0,2))
            /     (amrex::Math::abs(sten(ii-2,jj-1,0,3))
                  +amrex::Math::abs(sten(ii-1,jj-1,0,3))+eps);
        Real wpp = amrex::Math::abs(sten(ii-1,jj-1,0,3))*(1.0+wxp+wyp);
        cv +=           wpp*sten(ii-1,jj-1,0,4)*fine(ii-1,jj-1,0);

        Real wxm = amrex::Math::abs(sten(ii  ,jj-1,0,1))
            /     (amrex::Math::abs(sten(ii  ,jj-2,0,3))
                  +amrex::Math::abs(sten(ii  ,jj-1,0,3))+eps);
        wyp      = amrex::Math::abs(sten(ii+1,jj-1,0,2))
            /     (amrex::Math::abs(sten(ii  ,jj-1,0,3))
                  +amrex::Math::abs(sten(ii+1,jj-1,0,3))+eps);
        Real wmp = amrex::Math::abs(sten(ii  ,jj-1,0,3))*(1.0 + wxm + wyp);
        cv +=           wmp*sten(ii+1,jj-1,0,4)*fine(ii+1,jj-1,0);

        wxp      = amrex::Math::abs(sten(ii-1,jj+1,0,1))
            /     (amrex::Math::abs(sten(ii-1,jj  ,0,3))
                  +amrex::Math::abs(sten(ii-1,jj+1,0,3))+eps);
        Real wym = amrex::Math::abs(sten(ii-1,jj  ,0,2))
            /     (amrex::Math::abs(sten(ii-2,jj  ,0,3))
                  +amrex::Math::abs(sten(ii-1,jj  ,0,3))+eps);
        Real wpm = amrex::Math::abs(sten(ii-1,jj  ,0,3)) * (1.0 + wxp + wym);
        cv +=           wpm*sten(ii-1,jj+1,0,4)*fine(ii-1,jj+1,0);

        wxm      = amrex::Math::abs(sten(ii  ,jj+1,0,1))
            /     (amrex::Math::abs(sten(ii  ,jj  ,0,3))
                  +amrex::Math::abs(sten(ii  ,jj+1,0,3))+eps);
        wym      = amrex::Math::abs(sten(ii+1,jj  ,0,2))
            /     (amrex::Math::abs(sten(ii  ,jj  ,0,3))
                  +amrex::Math::abs(sten(ii+1,jj  ,0,3))+eps);
        Real wmm = amrex::Math::abs(sten(ii  ,jj  ,0,3)) * (1.0 + wxm + wym);
        cv +=           wmm*sten(ii+1,jj+1,0,4)*fine(ii+1,jj+1,0);

        crse(i,j,0) = cv * 0.25;
    }
}

#ifdef AMREX_USE_EB

namespace {
    constexpr int i_S_x     = 0;
    constexpr int i_S_y     = 1;
    constexpr int i_S_x2    = 2;
    constexpr int i_S_y2    = 3;
    constexpr int i_S_xy    = 4;
    constexpr int n_Sintg   = 5;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_connection (int i, int j, int, Array4<Real> const& conn,
                             Array4<Real const> const& intg, Array4<Real const> const& vol,
                             Array4<EBCellFlag const> const& flag) noexcept
{
    if (flag(i,j,0).isCovered()) {
        for (int n = 0; n < 6; ++n) conn(i,j,0,n) = Real(0.);
    } else if (flag(i,j,0).isRegular() or vol(i,j,0) >= almostone) {
        for (int n = 0; n < 6; ++n) conn(i,j,0,n) = Real(1.);
    } else {
        // Note that these are normalized so that they equal 1 in the case of a regular cell

        conn(i,j,0,0) = Real(3.)*(Real(0.25)*vol(i,j,0) + intg(i,j,0,i_S_y2) - intg(i,j,0,i_S_y));
        conn(i,j,0,1) = Real(6.)*(Real(0.25)*vol(i,j,0) - intg(i,j,0,i_S_y2));
        conn(i,j,0,2) = Real(3.)*(Real(0.25)*vol(i,j,0) + intg(i,j,0,i_S_y2) + intg(i,j,0,i_S_y));

        conn(i,j,0,3) = Real(3.)*(Real(0.25)*vol(i,j,0) + intg(i,j,0,i_S_x2) - intg(i,j,0,i_S_x));
        conn(i,j,0,4) = Real(6.)*(Real(0.25)*vol(i,j,0) - intg(i,j,0,i_S_x2));
        conn(i,j,0,5) = Real(3.)*(Real(0.25)*vol(i,j,0) + intg(i,j,0,i_S_x2) + intg(i,j,0,i_S_x));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_stencil_eb (int i, int j, int, Array4<Real> const& sten,
                             Array4<Real const> const& sig, Array4<Real const> const& conn,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = Real(1./6.)*dxinv[0]*dxinv[0];
    Real facy = Real(1./6.)*dxinv[1]*dxinv[1];

    sten(i,j,0,1) = Real(2.)*facx*(sig(i,j-1,0)*conn(i,j-1,0,2)+sig(i,j,0)*conn(i,j,0,0))
                            -facy*(sig(i,j-1,0)*conn(i,j-1,0,4)+sig(i,j,0)*conn(i,j,0,4));
    sten(i,j,0,2) = Real(2.)*facy*(sig(i-1,j,0)*conn(i-1,j,0,5)+sig(i,j,0)*conn(i,j,0,3))
                            -facx*(sig(i-1,j,0)*conn(i-1,j,0,1)+sig(i,j,0)*conn(i,j,0,1));
    sten(i,j,0,3) = (facx*conn(i,j,0,1)+facy*conn(i,j,0,4))*sig(i,j,0);
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_divu_eb (int i, int j, int, Array4<Real> const& rhs, Array4<Real const> const& vel,
                      Array4<Real const> const& vfrac, Array4<Real const> const& intg,
                      Array4<int const> const& msk, GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                      Box const& nodal_domain,
                      GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bclo,
                      GpuArray<LinOpBCType, AMREX_SPACEDIM> const& bchi) noexcept
{
    Real facx = Real(0.5)*dxinv[0];
    Real facy = Real(0.5)*dxinv[1];

    const auto domlo = amrex::lbound(nodal_domain);
    const auto domhi = amrex::ubound(nodal_domain);

    if (!msk(i,j,0)) {

        Real zero_ilo = 1.0;
        Real zero_ihi = 1.0;
        Real zero_jlo = 1.0;
        Real zero_jhi = 1.0;

        // The nodal divergence operator should not see the tangential velocity
        //     at an inflow face
        if ((bclo[0] == LinOpBCType::Neumann or bclo[0] == LinOpBCType::inflow)
            and i == domlo.x) zero_ilo = 0.0;
        if ((bchi[0] == LinOpBCType::Neumann or bchi[0] == LinOpBCType::inflow)
            and i == domhi.x) zero_ihi = 0.0;
        if ((bclo[1] == LinOpBCType::Neumann or bclo[1] == LinOpBCType::inflow)
            and j == domlo.y) zero_jlo = 0.0;
        if ((bchi[1] == LinOpBCType::Neumann or bchi[1] == LinOpBCType::inflow)
            and j == domhi.y) zero_jhi = 0.0; 

        rhs(i,j,0) = facx*(-vel(i-1,j-1,0,0)*(vfrac(i-1,j-1,0)+Real(2.)*intg(i-1,j-1,0,1))*zero_jlo
                           +vel(i  ,j-1,0,0)*(vfrac(i  ,j-1,0)+Real(2.)*intg(i  ,j-1,0,1))*zero_jlo
                           -vel(i-1,j  ,0,0)*(vfrac(i-1,j  ,0)-Real(2.)*intg(i-1,j  ,0,1))*zero_jhi
                           +vel(i  ,j  ,0,0)*(vfrac(i  ,j  ,0)-Real(2.)*intg(i  ,j  ,0,1))*zero_jhi)
                   + facy*(-vel(i-1,j-1,0,1)*(vfrac(i-1,j-1,0)+Real(2.)*intg(i-1,j-1,0,0))*zero_ilo
                           -vel(i  ,j-1,0,1)*(vfrac(i  ,j-1,0)-Real(2.)*intg(i  ,j-1,0,0))*zero_ihi
                           +vel(i-1,j  ,0,1)*(vfrac(i-1,j  ,0)+Real(2.)*intg(i-1,j  ,0,0))*zero_ilo
                           +vel(i  ,j  ,0,1)*(vfrac(i  ,j  ,0)-Real(2.)*intg(i  ,j  ,0,0))*zero_ihi);
    } else {
        rhs(i,j,0) = Real(0.);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_mknewu_eb (int i, int j, int, Array4<Real> const& u, Array4<Real const> const& p,
                        Array4<Real const> const& sig, Array4<Real const> const& vfrac,
                        Array4<Real const> const& intg, GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = Real(0.5)*dxinv[0];
    Real facy = Real(0.5)*dxinv[1];
    if (vfrac(i,j,0) == Real(0.)) {
        u(i,j,0,0) = u(i,j,0,1) = Real(0.);
    } else {
        Real dpdx = facx*(-p(i,j,0)+p(i+1,j,0)-p(i,j+1,0)+p(i+1,j+1,0));
        Real dpdy = facy*(-p(i,j,0)-p(i+1,j,0)+p(i,j+1,0)+p(i+1,j+1,0));
        Real dpp = (p(i,j,0)+p(i+1,j+1,0)-p(i+1,j,0)-p(i,j+1,0))/vfrac(i,j,0);
        u(i,j,0,0) -= sig(i,j,0)*(dpdx + dxinv[0]*intg(i,j,0,1)*dpp);
        u(i,j,0,1) -= sig(i,j,0)*(dpdy + dxinv[1]*intg(i,j,0,0)*dpp);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mlndlap_rhcc_eb (int i, int j, int, Array4<Real const> const& rhcc,
                      Array4<Real const> const& vfrac, Array4<Real const> const& intg,
                      Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,0)) {
        return
            rhcc(i  ,j  ,0)*(Real(0.25)*vfrac(i  ,j  ,0)-intg(i  ,j  ,0,i_S_x)-intg(i  ,j  ,0,i_S_y)+intg(i  ,j  ,0,i_S_xy)) +
            rhcc(i-1,j  ,0)*(Real(0.25)*vfrac(i-1,j  ,0)+intg(i-1,j  ,0,i_S_x)-intg(i-1,j  ,0,i_S_y)-intg(i-1,j  ,0,i_S_xy)) +
            rhcc(i-1,j-1,0)*(Real(0.25)*vfrac(i-1,j-1,0)+intg(i-1,j-1,0,i_S_x)+intg(i-1,j-1,0,i_S_y)+intg(i-1,j-1,0,i_S_xy)) +
            rhcc(i  ,j-1,0)*(Real(0.25)*vfrac(i  ,j-1,0)-intg(i  ,j-1,0,i_S_x)+intg(i  ,j-1,0,i_S_y)-intg(i  ,j-1,0,i_S_xy));
    } else {
        return Real(0.);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_integral (int i, int j, int, Array4<Real> const& intg) noexcept
{
    intg(i,j,0,i_S_x ) = Real(0.);
    intg(i,j,0,i_S_y ) = Real(0.);
    intg(i,j,0,i_S_x2) = Real(1./12.);
    intg(i,j,0,i_S_y2) = Real(1./12.);
    intg(i,j,0,i_S_xy) = Real(0.);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndlap_set_integral_eb (int i, int j, int, Array4<Real> const& intg,
                              Array4<EBCellFlag const> const& flag, Array4<Real const> const& vol,
                              Array4<Real const> const& ax, Array4<Real const> const& ay,
                              Array4<Real const> const& bcen) noexcept
{
    if (flag(i,j,0).isCovered()) {
        intg(i,j,0,i_S_x ) = Real(0.);
        intg(i,j,0,i_S_y ) = Real(0.);
        intg(i,j,0,i_S_x2) = Real(0.);
        intg(i,j,0,i_S_y2) = Real(0.);
        intg(i,j,0,i_S_xy) = Real(0.);
    } else if (flag(i,j,0).isRegular() or vol(i,j,0) >= almostone) {
        intg(i,j,0,i_S_x ) = Real(0.);
        intg(i,j,0,i_S_y ) = Real(0.);
        intg(i,j,0,i_S_x2) = Real(1./12.);
        intg(i,j,0,i_S_y2) = Real(1./12.);
        intg(i,j,0,i_S_xy) = Real(0.);
    } else {
        Real axm = ax(i,j,0);
        Real axp = ax(i+1,j,0);
        Real aym = ay(i,j,0);
        Real ayp = ay(i,j+1,0);

        Real apnorm = std::sqrt((axm-axp)*(axm-axp) + (aym-ayp)*(aym-ayp));
        if (apnorm == Real(0.)) {
            amrex::Abort("amrex_mlndlap_set_integral: we are in trouble");
        }

        Real apnorminv = Real(1.)/apnorm;
        Real anrmx = (axm-axp) * apnorminv;  // pointing to the wall
        Real anrmy = (aym-ayp) * apnorminv;

        Real bcx = bcen(i,j,0,0);
        Real bcy = bcen(i,j,0,1);

        Real Sx, Sy, Sx2, Sy2, Sxy;
        if (anrmx == Real(0.)) {
            Sx = Real(0.);
            Sx2 = Real(1./24.)*(axm+axp);
            Sxy = Real(0.);
        } else if (anrmy == Real(0.)) {
            Sx  = Real(1./8.) *(axp-axm) + anrmx*Real(0.5)*(bcx*bcx);
            Sx2 = Real(1./24.)*(axp+axm) + anrmx*Real(1./3.)*(bcx*bcx*bcx);
            Sxy = Real(0.);
        } else {
            Real xmin, xmax;
            if (anrmx > Real(0.)) {
                xmin = Real(-0.5) + amrex::min(aym,ayp);
                xmax = Real(-0.5) + amrex::max(aym,ayp);
            } else {
                xmin = Real(0.5) - amrex::max(aym,ayp);
                xmax = Real(0.5) - amrex::min(aym,ayp);
            }
            Real xmin3 = xmin*xmin*xmin;
            Real xmin4 = xmin3*xmin;
            Real xmax3 = xmax*xmax*xmax;
            Real xmax4 = xmax3*xmax;
            Sx  = Real(1./8.) *(axp-axm) + (anrmx/amrex::Math::abs(anrmy))*Real(1./6.) *(xmax3-xmin3);
            Sx2 = Real(1./24.)*(axp+axm) + (anrmx/amrex::Math::abs(anrmy))*Real(1./12.)*(xmax4-xmin4);

            Real kk = -anrmx/anrmy;
            Real bb = bcy-kk*bcx;
            Sxy = Real(1./8.)*kk*kk*(xmax4-xmin4) + Real(1./3.)*kk*bb*(xmax3-xmin3)
                + (Real(0.25)*bb*bb-Real(1./16.))*(xmax*xmax-xmin*xmin);
            Sxy = amrex::Math::copysign(Sxy, anrmy);
        }

        if (anrmy == Real(0.)) {
            Sy = Real(0.);
            Sy2 = Real(1./24.)*(aym+ayp);
        } else if (anrmx == Real(0.)) {
            Sy  = Real(1./8.) *(ayp-aym) + anrmy*Real(0.5)*(bcy*bcy);
            Sy2 = Real(1./24.)*(ayp+aym) + anrmy*Real(1./3.)*(bcy*bcy*bcy);
        } else {
            Real ymin, ymax;
            if (anrmy > Real(0.)) {
                ymin = Real(-0.5) + amrex::min(axm,axp);
                ymax = Real(-0.5) + amrex::max(axm,axp);
            } else {
                ymin = Real(0.5) - amrex::max(axm,axp);
                ymax = Real(0.5) - amrex::min(axm,axp);
            }
            Real ymin3 = ymin*ymin*ymin;
            Real ymin4 = ymin3*ymin;
            Real ymax3 = ymax*ymax*ymax;
            Real ymax4 = ymax3*ymax;
            Sy  = Real(1./8.) *(ayp-aym) + (anrmy/amrex::Math::abs(anrmx))*Real(1./6.) *(ymax3-ymin3);
            Sy2 = Real(1./24.)*(ayp+aym) + (anrmy/amrex::Math::abs(anrmx))*Real(1./12.)*(ymax4-ymin4);
        }

        intg(i,j,0,i_S_x ) = Sx;
        intg(i,j,0,i_S_y ) = Sy;
        intg(i,j,0,i_S_x2) = Sx2;
        intg(i,j,0,i_S_y2) = Sy2;
        intg(i,j,0,i_S_xy) = Sxy;
    }
}

#endif

}
#endif
