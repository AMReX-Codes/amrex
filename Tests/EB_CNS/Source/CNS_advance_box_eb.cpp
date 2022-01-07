
#include <CNS.H>
#include <CNS_hydro_K.H>
#include <CNS_hydro_eb_K.H>
#include <CNS_divop_K.H>
#include <CNS_diffusion_eb_K.H>

#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>

#if (AMREX_SPACEDIM == 2)
#include <AMReX_EBMultiFabUtil_2D_C.H>
#elif (AMREX_SPACEDIM == 3)
#include <AMReX_EBMultiFabUtil_3D_C.H>
#endif

using namespace amrex;

void
CNS::compute_dSdt_box_eb (const Box& bx,
                          Array4<Real const> const& sfab,
                          Array4<Real      > const& dsdtfab,
                          std::array<FArrayBox*, AMREX_SPACEDIM> const& flux,
                          Array4<EBCellFlag const> const& flag,
                          Array4<Real       const> const& vfrac,
                          AMREX_D_DECL(
                          Array4<Real       const> const& apx,
                          Array4<Real       const> const& apy,
                          Array4<Real       const> const& apz),
                          AMREX_D_DECL(
                          Array4<Real       const> const& fcx,
                          Array4<Real       const> const& fcy,
                          Array4<Real       const> const& fcz),
                          Array4<Real       const> const& bcent,
                          int as_crse,
                          Array4<Real            > const& drho_as_crse,
                          Array4<int        const> const& rrflag_as_crse,
                          int as_fine,
                          Array4<Real            > const& dm_as_fine,
                          Array4<int        const> const& lev_mask,
                          Real dt)
{
    BL_PROFILE("CNS::compute_dSdt_box_eb()");

    const Box& bxg1 = amrex::grow(bx,1);
    const Box& bxg2 = amrex::grow(bx,2);
    const Box& bxg3 = amrex::grow(bx,3);
    const Box& bxg4 = amrex::grow(bx,4);
    const Box& bxg5 = amrex::grow(bx,5);

    const auto dxinv = geom.InvCellSizeArray();

    // Quantities for redistribution
    FArrayBox divc,optmp,redistwgt,delta_m;
    divc.resize(bxg2,NEQNS);
    optmp.resize(bxg2,NEQNS);
    delta_m.resize(bxg1,NEQNS);
    redistwgt.resize(bxg2,1);

    // Set to zero just in case
    divc.setVal<RunOn::Device>(0.0);
    optmp.setVal<RunOn::Device>(0.0);
    delta_m.setVal<RunOn::Device>(0.0);
    redistwgt.setVal<RunOn::Device>(0.0);

    // Primitive variables
    FArrayBox qtmp;
    qtmp.resize(bxg5, NPRIM);

    // Slopes
    FArrayBox slopetmp;
    slopetmp.resize(bxg4,NEQNS);

    FArrayBox diff_coeff;

    if (do_visc == 1)
    {
       diff_coeff.resize(bxg4, NCOEF);
    }

    FArrayBox flux_tmp[AMREX_SPACEDIM];
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        flux_tmp[idim].resize(amrex::surroundingNodes(bxg3,idim),NCONS);
        flux_tmp[idim].setVal<RunOn::Device>(0.);
    }

    Parm const* lparm = d_parm;

    AMREX_D_TERM(auto const& fxfab = flux_tmp[0].array();,
                 auto const& fyfab = flux_tmp[1].array();,
                 auto const& fzfab = flux_tmp[2].array(););

    auto const& q = qtmp.array();

    GpuArray<Real,3> weights;
    weights[0] = 0.;
    weights[1] = 1.;
    weights[2] = 0.5;

    // Initialize dm_as_fine to 0
    if (as_fine)
    {
        amrex::ParallelFor(bxg1, NCONS,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
           dm_as_fine(i,j,k,n) = 0.;
        });
    }

    amrex::ParallelFor(bxg5,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cns_ctoprim(i, j, k, sfab, q, *lparm);
    });

    if (do_visc == 1)
    {
       auto const& coefs = diff_coeff.array();
       if(use_const_visc == 1 ) {
          amrex::ParallelFor(bxg4,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              cns_constcoef_eb(i, j, k, flag, coefs, *lparm);
          });
       } else {
          amrex::ParallelFor(bxg4,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              cns_diffcoef_eb(i, j, k, q, flag, coefs, *lparm);
          });
       }
    }

    auto const& slope = slopetmp.array();

    auto l_plm_iorder = plm_iorder;
    auto l_plm_theta = plm_theta;

    // x-direction
    int cdir = 0;
    const Box& xslpbx = amrex::grow(bxg3, cdir, 1);
    amrex::ParallelFor(xslpbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cns_slope_eb_x(i, j, k, slope, q, flag, l_plm_iorder, l_plm_theta);
    });

    const Box& xflxbx = amrex::surroundingNodes(bxg3,cdir);
    amrex::ParallelFor(xflxbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cns_riemann_x(i, j, k, fxfab, slope, q, *lparm);
        for (int n = NEQNS; n < NCONS; ++n) fxfab(i,j,k,n) = Real(0.0);
    });


    if (do_visc == 1)
    {
        auto const& coefs = diff_coeff.array();
        amrex::ParallelFor(xflxbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cns_diff_eb_x(i, j, k, q, coefs, flag, dxinv, weights, fxfab);
        });
    }


    // y-direction
    cdir = 1;
    const Box& yslpbx = amrex::grow(bxg3, cdir, 1);
    amrex::ParallelFor(yslpbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cns_slope_eb_y(i, j, k, slope, q, flag, l_plm_iorder, l_plm_theta);
    });

    const Box& yflxbx = amrex::surroundingNodes(bxg3,cdir);
    amrex::ParallelFor(yflxbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cns_riemann_y(i, j, k, fyfab, slope, q, *lparm);
        for (int n = NEQNS; n < NCONS; ++n) fyfab(i,j,k,n) = Real(0.0);
    });

    if(do_visc == 1)
    {
       auto const& coefs = diff_coeff.array();
       amrex::ParallelFor(yflxbx,
       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
           cns_diff_eb_y(i, j, k, q, coefs, flag, dxinv, weights, fyfab);
       });
    }

#if (AMREX_SPACEDIM == 3)
    // z-direction
    cdir = 2;
    const Box& zslpbx = amrex::grow(bxg3, cdir, 1);
    amrex::ParallelFor(zslpbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cns_slope_eb_z(i, j, k, slope, q, flag, l_plm_iorder, l_plm_theta);
    });

    const Box& zflxbx = amrex::surroundingNodes(bxg3,cdir);
    amrex::ParallelFor(zflxbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        cns_riemann_z(i, j, k, fzfab, slope, q, *lparm);
        for (int n = NEQNS; n < NCONS; ++n) fzfab(i,j,k,n) = Real(0.0);
    });

    if(do_visc == 1)
    {
       auto const& coefs = diff_coeff.array();
       amrex::ParallelFor(zflxbx,
       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
           cns_diff_eb_z(i, j, k, q, coefs, flag, dxinv, weights, fzfab);
       });
    }
#endif

    // These are the fluxes we computed above -- they live at face centers
    AMREX_D_TERM(auto const& fx_in_arr = flux_tmp[0].array();,
                 auto const& fy_in_arr = flux_tmp[1].array();,
                 auto const& fz_in_arr = flux_tmp[2].array(););

    // These are the fluxes on face centroids -- they are defined in eb_compute_div
    //    and are the fluxes that go into the flux registers
    AMREX_D_TERM(auto const& fx_out_arr = flux[0]->array();,
                 auto const& fy_out_arr = flux[1]->array();,
                 auto const& fz_out_arr = flux[2]->array(););

    auto const& blo = bx.smallEnd();
    auto const& bhi = bx.bigEnd();

    // Because we are going to redistribute, we put the divergence into divc
    //    rather than directly into dsdtfab
    auto const& divc_arr = divc.array();

    bool l_do_visc = do_visc;
    auto l_eb_weights_type = eb_weights_type;

    auto const& coefs = diff_coeff.array();
    auto const& redistwgt_arr = redistwgt.array();

    amrex::ParallelFor(bxg2, NEQNS,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
       // This does the divergence but not the redistribution -- we will do that later
       // We do compute the weights here though
       eb_compute_div(i,j,k,n,blo,bhi,q,divc_arr,
                      AMREX_D_DECL(fx_in_arr ,fy_in_arr ,fz_in_arr),
                      AMREX_D_DECL(fx_out_arr,fy_out_arr,fz_out_arr),
                      flag, vfrac, bcent, coefs, redistwgt_arr,
                      AMREX_D_DECL(apx, apy, apz),
                      AMREX_D_DECL(fcx, fcy, fcz), dxinv, *lparm, l_eb_weights_type, l_do_visc);
    });

    auto const& optmp_arr = optmp.array();
    auto const& del_m_arr = delta_m.array();

    // Now do redistribution
    cns_flux_redistribute(bx,dsdtfab,divc_arr,optmp_arr,del_m_arr,redistwgt_arr,vfrac,flag,
                          as_crse, drho_as_crse, rrflag_as_crse, as_fine, dm_as_fine, lev_mask, dt);

    if (gravity != Real(0.0))
    {
        const Real g = gravity;
        const int irho = URHO;
#if (AMREX_SPACEDIM == 2)
        const int imz = UMY;
#elif (AMREX_SPACEDIM == 3)
        const int imz = UMZ;
#endif
        const int irhoE = UEDEN;
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            dsdtfab(i,j,k,imz  ) += g * sfab(i,j,k,irho);
            dsdtfab(i,j,k,irhoE) += g * sfab(i,j,k,imz);
        });
    }

    Gpu::streamSynchronize();
}
