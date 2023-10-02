#ifndef AMREX_MLNODEABECLAP_3D_K_H_
#define AMREX_MLNODEABECLAP_3D_K_H_

namespace amrex {

inline void
mlndabeclap_gauss_seidel_aa (Box const& bx, Array4<Real> const& sol,
                             Array4<Real const> const& rhs,
                             Real alpha, Real beta,
                             Array4<Real const> const& acf,
                             Array4<Real const> const& bcf,
                             Array4<int const> const& msk,
                             GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    Real facx = Real(1.0/36.0)*dxinv[0]*dxinv[0];
    Real facy = Real(1.0/36.0)*dxinv[1]*dxinv[1];
    Real facz = Real(1.0/36.0)*dxinv[2]*dxinv[2];
    Real fxyz = facx + facy + facz;
    Real fmx2y2z = -facx + Real(2.0)*facy + Real(2.0)*facz;
    Real f2xmy2z = Real(2.0)*facx - facy + Real(2.0)*facz;
    Real f2x2ymz = Real(2.0)*facx + Real(2.0)*facy - facz;
    Real f4xm2ym2z = Real(4.0)*facx - Real(2.0)*facy - Real(2.0)*facz;
    Real fm2x4ym2z = -Real(2.0)*facx + Real(4.0)*facy - Real(2.0)*facz;
    Real fm2xm2y4z = -Real(2.0)*facx - Real(2.0)*facy + Real(4.0)*facz;

    amrex::LoopOnCpu(bx, [=] (int i, int j, int k) noexcept
    {
        if (msk(i,j,k)) {
            sol(i,j,k) = Real(0.0);
        } else {
            Real s0 = Real(-4.0)*fxyz*(bcf(i-1,j-1,k-1)+bcf(i,j-1,k-1)+bcf(i-1,j,k-1)+bcf(i,j,k-1)
                                      +bcf(i-1,j-1,k  )+bcf(i,j-1,k  )+bcf(i-1,j,k  )+bcf(i,j,k  ));
            Real lap = sol(i,j,k)*s0
                + fxyz*(sol(i-1,j-1,k-1)*bcf(i-1,j-1,k-1)
                      + sol(i+1,j-1,k-1)*bcf(i  ,j-1,k-1)
                      + sol(i-1,j+1,k-1)*bcf(i-1,j  ,k-1)
                      + sol(i+1,j+1,k-1)*bcf(i  ,j  ,k-1)
                      + sol(i-1,j-1,k+1)*bcf(i-1,j-1,k  )
                      + sol(i+1,j-1,k+1)*bcf(i  ,j-1,k  )
                      + sol(i-1,j+1,k+1)*bcf(i-1,j  ,k  )
                      + sol(i+1,j+1,k+1)*bcf(i  ,j  ,k  ))
                + fmx2y2z*(sol(i  ,j-1,k-1)*(bcf(i-1,j-1,k-1)+bcf(i,j-1,k-1))
                         + sol(i  ,j+1,k-1)*(bcf(i-1,j  ,k-1)+bcf(i,j  ,k-1))
                         + sol(i  ,j-1,k+1)*(bcf(i-1,j-1,k  )+bcf(i,j-1,k  ))
                         + sol(i  ,j+1,k+1)*(bcf(i-1,j  ,k  )+bcf(i,j  ,k  )))
                + f2xmy2z*(sol(i-1,j  ,k-1)*(bcf(i-1,j-1,k-1)+bcf(i-1,j,k-1))
                         + sol(i+1,j  ,k-1)*(bcf(i  ,j-1,k-1)+bcf(i  ,j,k-1))
                         + sol(i-1,j  ,k+1)*(bcf(i-1,j-1,k  )+bcf(i-1,j,k  ))
                         + sol(i+1,j  ,k+1)*(bcf(i  ,j-1,k  )+bcf(i  ,j,k  )))
                + f2x2ymz*(sol(i-1,j-1,k  )*(bcf(i-1,j-1,k-1)+bcf(i-1,j-1,k))
                         + sol(i+1,j-1,k  )*(bcf(i  ,j-1,k-1)+bcf(i  ,j-1,k))
                         + sol(i-1,j+1,k  )*(bcf(i-1,j  ,k-1)+bcf(i-1,j  ,k))
                         + sol(i+1,j+1,k  )*(bcf(i  ,j  ,k-1)+bcf(i  ,j  ,k)))
                + f4xm2ym2z*(sol(i-1,j,k)*(bcf(i-1,j-1,k-1)+bcf(i-1,j,k-1)+bcf(i-1,j-1,k)+bcf(i-1,j,k))
                           + sol(i+1,j,k)*(bcf(i  ,j-1,k-1)+bcf(i  ,j,k-1)+bcf(i  ,j-1,k)+bcf(i  ,j,k)))
                + fm2x4ym2z*(sol(i,j-1,k)*(bcf(i-1,j-1,k-1)+bcf(i,j-1,k-1)+bcf(i-1,j-1,k)+bcf(i,j-1,k))
                           + sol(i,j+1,k)*(bcf(i-1,j  ,k-1)+bcf(i,j  ,k-1)+bcf(i-1,j  ,k)+bcf(i,j  ,k)))
                + fm2xm2y4z*(sol(i,j,k-1)*(bcf(i-1,j-1,k-1)+bcf(i,j-1,k-1)+bcf(i-1,j,k-1)+bcf(i,j,k-1))
                           + sol(i,j,k+1)*(bcf(i-1,j-1,k  )+bcf(i,j-1,k  )+bcf(i-1,j,k  )+bcf(i,j,k  )));
            Real Ax = alpha*acf(i,j,k)*sol(i,j,k) - beta*lap;

            sol(i,j,k) += (rhs(i,j,k) - Ax) / (alpha*acf(i,j,k)-beta*s0);
        }
    });
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
mlndabeclap_jacobi_aa (int i, int j, int k, Array4<Real> const& sol,
                       Real lap, Array4<Real const> const& rhs,
                       Real alpha, Real beta,
                       Array4<Real const> const& acf,
                       Array4<Real const> const& bcf,
                       Array4<int const> const& msk,
                       GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    if (msk(i,j,k)) {
        sol(i,j,k) = Real(0.0);
    } else {
        Real fxyz = Real(-4.0 / 36.0)*(dxinv[0]*dxinv[0] +
                                       dxinv[1]*dxinv[1] +
                                       dxinv[2]*dxinv[2]);
        Real s0 = fxyz*(bcf(i-1,j-1,k-1)+bcf(i,j-1,k-1)+bcf(i-1,j,k-1)+bcf(i,j,k-1)
                       +bcf(i-1,j-1,k  )+bcf(i,j-1,k  )+bcf(i-1,j,k  )+bcf(i,j,k));
        Real Ax = alpha*acf(i,j,k)*sol(i,j,k) - beta*lap;

        sol(i,j,k) += Real(2.0/3.0) * (rhs(i,j,k) - Ax)
            / (alpha*acf(i,j,k)-beta*s0);
    }
}

}

#endif
