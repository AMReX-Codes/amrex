#ifndef AMREX_MLNODEABECLAP_2D_K_H_
#define AMREX_MLNODEABECLAP_2D_K_H_

namespace amrex {

inline void
mlndabeclap_gauss_seidel_aa (Box const& bx, Array4<Real> const& sol,
                             Array4<Real const> const& rhs,
                             Real scalarA, Real scalarB,
                             Array4<Real const> const& alphacf,
                             Array4<Real const> const& betacf,
                             Array4<int const> const& msk,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = Real(1.0/6.0)*dxinv[0]*dxinv[0];
    Real facy = Real(1.0/6.0)*dxinv[1]*dxinv[1];
    Real fxy = facx + facy;
    Real f2xmy = Real(2.0)*facx - facy;
    Real fmx2y = Real(2.0)*facy - facx;

    amrex::Loop(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = Real(0.0);
        } else {
            Real s0 = (-Real(2.0))*fxy*(betacf(i-1,j-1,k)+betacf(i,j-1,k)+betacf(i-1,j,k)+betacf(i,j,k));
            Real lap =  sol(i-1,j-1,k)*fxy*betacf(i-1,j-1,k)
                      + sol(i+1,j-1,k)*fxy*betacf(i  ,j-1,k)
                      + sol(i-1,j+1,k)*fxy*betacf(i-1,j  ,k)
                      + sol(i+1,j+1,k)*fxy*betacf(i  ,j  ,k)
                      + sol(i-1,j,k)*f2xmy*(betacf(i-1,j-1,k)+betacf(i-1,j,k))
                      + sol(i+1,j,k)*f2xmy*(betacf(i  ,j-1,k)+betacf(i  ,j,k))
                      + sol(i,j-1,k)*fmx2y*(betacf(i-1,j-1,k)+betacf(i,j-1,k))
                      + sol(i,j+1,k)*fmx2y*(betacf(i-1,j  ,k)+betacf(i,j  ,k))
                      + sol(i,j,k)*s0;
            Real Ax = scalarA*alphacf(i,j,k)*sol(i,j,k) - scalarB*lap;

            sol(i,j,k) += (rhs(i,j,k) - Ax) / (scalarA*alphacf(i,j,k)-scalarB*s0);
        }
    });
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
mlndabeclap_jacobi_aa (int i, int j, int k, Array4<Real> const& sol,
                       Real lap, Array4<Real const> const& rhs,
                       Real scalarA, Real scalarB,
                       Array4<Real const> const& alphacf,
                       Array4<Real const> const& betacf,
                       Array4<int const> const& msk,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        sol(i,j,k) = Real(0.0);
    } else {
        Real fac = -Real(2.0/6.0)*(dxinv[0]*dxinv[0] + dxinv[1]*dxinv[1]);
        Real s0 = fac*(betacf(i-1,j-1,k)+betacf(i,j-1,k)+betacf(i-1,j,k)+betacf(i,j,k));
        Real Ax = scalarA*alphacf(i,j,k)*sol(i,j,k) - scalarB*lap;

        sol(i,j,k) += Real(2.0/3.0) * (rhs(i,j,k) - Ax)
            / (scalarA*alphacf(i,j,k)-scalarB*s0);
    }

}

}

#endif
