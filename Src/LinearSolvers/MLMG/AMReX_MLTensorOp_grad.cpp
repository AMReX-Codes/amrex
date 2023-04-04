#include <AMReX_MLTensorOp.H>
#include <AMReX_MLTensor_K.H>

namespace amrex {

void
MLTensorOp::compFlux (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes,
                       MultiFab& sol, Location loc) const
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(amrlev, fluxes, sol, loc);
#else
    BL_PROFILE("MLTensorOp::compFlux()");

    const int mglev = 0;
    const int ncomp = getNComp();
    MLABecLaplacian::compFlux(amrlev, fluxes, sol, loc);

    MLMGBndry const* bndry = m_bndry_sol[amrlev].get();
    applyBCTensor(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution, bndry);

    const auto& bcondloc = *m_bcondloc[amrlev][mglev];

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const Box& domain = m_geom[amrlev][mglev].growPeriodicDomain(1);
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    Array<MultiFab,AMREX_SPACEDIM> const& etamf = m_b_coeffs[amrlev][mglev];
    Array<MultiFab,AMREX_SPACEDIM> const& kapmf = m_kappa[amrlev][mglev];
    Real bscalar = m_b_scalar;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        Array<FArrayBox,AMREX_SPACEDIM> fluxfab_tmp;

        for (MFIter mfi(sol, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<Real const> const vfab = sol.const_array(mfi);
            AMREX_D_TERM(Array4<Real const> const etaxfab = etamf[0].const_array(mfi);,
                         Array4<Real const> const etayfab = etamf[1].const_array(mfi);,
                         Array4<Real const> const etazfab = etamf[2].const_array(mfi););
            AMREX_D_TERM(Array4<Real const> const kapxfab = kapmf[0].const_array(mfi);,
                         Array4<Real const> const kapyfab = kapmf[1].const_array(mfi);,
                         Array4<Real const> const kapzfab = kapmf[2].const_array(mfi););
            AMREX_D_TERM(Box const xbx = mfi.nodaltilebox(0);,
                         Box const ybx = mfi.nodaltilebox(1);,
                         Box const zbx = mfi.nodaltilebox(2););
            AMREX_D_TERM(fluxfab_tmp[0].resize(xbx,AMREX_SPACEDIM);,
                         fluxfab_tmp[1].resize(ybx,AMREX_SPACEDIM);,
                         fluxfab_tmp[2].resize(zbx,AMREX_SPACEDIM););
            AMREX_D_TERM(Elixir fxeli = fluxfab_tmp[0].elixir();,
                         Elixir fyeli = fluxfab_tmp[1].elixir();,
                         Elixir fzeli = fluxfab_tmp[2].elixir(););
            AMREX_D_TERM(Array4<Real> const fxfab = fluxfab_tmp[0].array();,
                         Array4<Real> const fyfab = fluxfab_tmp[1].array();,
                         Array4<Real> const fzfab = fluxfab_tmp[2].array(););

            if (domain.strictly_contains(mfi.tilebox())) {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
                ( xbx, txbx,
                  {
                      mltensor_cross_terms_fx(txbx,fxfab,vfab,etaxfab,kapxfab,dxinv);
                  }
                , ybx, tybx,
                  {
                      mltensor_cross_terms_fy(tybx,fyfab,vfab,etayfab,kapyfab,dxinv);
                  }
                , zbx, tzbx,
                  {
                      mltensor_cross_terms_fz(tzbx,fzfab,vfab,etazfab,kapzfab,dxinv);
                  }
                );
            } else {
                const auto & bdcv = bcondloc.bndryConds(mfi);

                Array2D<BoundCond,0,2*AMREX_SPACEDIM,0,AMREX_SPACEDIM> bct;
                for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                    for (OrientationIter face; face; ++face) {
                        Orientation ori = face();
                        bct(ori,icomp) = bdcv[icomp][ori];
                    }
                }

                const auto& bvxlo = (*bndry)[Orientation(0,Orientation::low )].array(mfi);
                const auto& bvylo = (*bndry)[Orientation(1,Orientation::low )].array(mfi);
                const auto& bvxhi = (*bndry)[Orientation(0,Orientation::high)].array(mfi);
                const auto& bvyhi = (*bndry)[Orientation(1,Orientation::high)].array(mfi);
#if (AMREX_SPACEDIM == 3)
                const auto& bvzlo = (*bndry)[Orientation(2,Orientation::low )].array(mfi);
                const auto& bvzhi = (*bndry)[Orientation(2,Orientation::high)].array(mfi);
#endif

                AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
                ( xbx, txbx,
                  {
                      mltensor_cross_terms_fx(txbx,fxfab,vfab,etaxfab,kapxfab,dxinv,
                                              bvxlo, bvxhi, bct, dlo, dhi);
                  }
                , ybx, tybx,
                  {
                      mltensor_cross_terms_fy(tybx,fyfab,vfab,etayfab,kapyfab,dxinv,
                                              bvylo, bvyhi, bct, dlo, dhi);
                  }
                , zbx, tzbx,
                  {
                      mltensor_cross_terms_fz(tzbx,fzfab,vfab,etazfab,kapzfab,dxinv,
                                              bvzlo, bvzhi, bct, dlo, dhi);
                  }
                );
            }

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const Box& nbx = mfi.nodaltilebox(idim);
                Array4<Real      > dst = fluxes[idim]->array(mfi);
                Array4<Real const> src = fluxfab_tmp[idim].const_array();
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (nbx, ncomp, i, j, k, n,
                {
                    dst(i,j,k,n) += bscalar*src(i,j,k,n);
                });
            }

        }
    }
#endif
}



void
MLTensorOp::compVelGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes, // NOLINT(readability-convert-member-functions-to-static)
                       MultiFab& sol, Location /*loc*/) const
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(amrlev,fluxes,sol);
#else
    BL_PROFILE("MLTensorOp::compVelGrad()");

    const int mglev = 0;

    MLMGBndry const* bndry = m_bndry_sol[amrlev].get();
    applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution, bndry);
    applyBCTensor(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution, bndry);

    const auto& bcondloc = *m_bcondloc[amrlev][mglev];

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const Box& domain = m_geom[amrlev][mglev].growPeriodicDomain(1);
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sol, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real const> const vfab = sol.const_array(mfi);
        AMREX_D_TERM(Box const xbx = mfi.nodaltilebox(0);,
                     Box const ybx = mfi.nodaltilebox(1);,
                     Box const zbx = mfi.nodaltilebox(2);)
        AMREX_D_TERM(Array4<Real> const fxfab = fluxes[0]->array(mfi);,
                     Array4<Real> const fyfab = fluxes[1]->array(mfi);,
                     Array4<Real> const fzfab = fluxes[2]->array(mfi);)

// The derivatives are put in the array with the following order:
// component: 0    ,  1    ,  2    ,  3    ,  4    , 5    ,  6    ,  7    ,  8
// in 2D:     dU/dx,  dV/dx,  dU/dy,  dV/dy
// in 3D:     dU/dx,  dV/dx,  dW/dx,  dU/dy,  dV/dy, dW/dy,  dU/dz,  dV/dz,  dW/dz

        if (domain.strictly_contains(mfi.tilebox())) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
            ( xbx, txbx,
              {
                  mltensor_vel_grads_fx(txbx,fxfab,vfab,dxinv);
              }
            , ybx, tybx,
              {
                  mltensor_vel_grads_fy(tybx,fyfab,vfab,dxinv);
              }
            , zbx, tzbx,
              {
                  mltensor_vel_grads_fz(tzbx,fzfab,vfab,dxinv);
              }
            );
        } else {
            const auto & bdcv = bcondloc.bndryConds(mfi);

            Array2D<BoundCond,0,2*AMREX_SPACEDIM,0,AMREX_SPACEDIM> bct;
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                for (OrientationIter face; face; ++face) {
                    Orientation ori = face();
                    bct(ori,icomp) = bdcv[icomp][ori];
                }
            }

            const auto& bvxlo = (*bndry)[Orientation(0,Orientation::low )].array(mfi);
            const auto& bvylo = (*bndry)[Orientation(1,Orientation::low )].array(mfi);
            const auto& bvxhi = (*bndry)[Orientation(0,Orientation::high)].array(mfi);
            const auto& bvyhi = (*bndry)[Orientation(1,Orientation::high)].array(mfi);
#if (AMREX_SPACEDIM == 3)
            const auto& bvzlo = (*bndry)[Orientation(2,Orientation::low )].array(mfi);
            const auto& bvzhi = (*bndry)[Orientation(2,Orientation::high)].array(mfi);
#endif
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM
            ( xbx, txbx,
              {
                  mltensor_vel_grads_fx(txbx,fxfab,vfab,dxinv,bvxlo,bvxhi,bct,dlo,dhi);
              }
            , ybx, tybx,
              {
                  mltensor_vel_grads_fy(tybx,fyfab,vfab,dxinv,bvylo,bvyhi,bct,dlo,dhi);
              }
            , zbx, tzbx,
              {
                  mltensor_vel_grads_fz(tzbx,fzfab,vfab,dxinv,bvzlo,bvzhi,bct,dlo,dhi);
              }
            );
        }
    }
#endif
}

}
