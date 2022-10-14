
#include "CNS.H"
#include "CNS_hydro_K.H"
#include "CNS_diffusion_K.H"
#include "CNS_K.H"

using namespace amrex;

Real
CNS::advance (Real time, Real dt, int iteration, int ncycle)
{
    BL_PROFILE("CNS::advance()");

    for (int i = 0; i < num_state_data_types; ++i) {
        state[i].allocOldData();
        state[i].swapTimeLevels(dt);
    }

    FluxRegister* fr_as_crse = nullptr;
    if (do_reflux && level < parent->finestLevel()) {
        CNS& fine_level = getLevel(level+1);
        fr_as_crse = fine_level.flux_reg.get();
    }

    FluxRegister* fr_as_fine = nullptr;
    if (do_reflux && level > 0) {
        fr_as_fine = flux_reg.get();
    }

    if (fr_as_crse) {
        fr_as_crse->setVal(Real(0.0));
    }

    RK(rk_order, State_Type, time, dt, iteration, ncycle,
       // Given state S, compute dSdt. dtsub is needed for flux register operations
       [&] (int /*stage*/, MultiFab& dSdt, MultiFab const& S,
            Real /*t*/, Real dtsub) {
           compute_dSdt(S, dSdt, dtsub, fr_as_crse, fr_as_fine);
       },
       // Optional. In case if there is anything needed after each RK substep.
       [&] (int /*stage*/, MultiFab& S) { computeTemp(S,0); });

    return dt;
}

void
CNS::compute_dSdt (const MultiFab& S, MultiFab& dSdt, Real dt,
                   FluxRegister* fr_as_crse, FluxRegister* fr_as_fine)
{
    BL_PROFILE("CNS::compute_dSdt()");

    const auto dx = geom.CellSizeArray();
    const auto dxinv = geom.InvCellSizeArray();
    const int ncomp = NUM_STATE;
    const int neqns = 5;
    const int ncons = 7;
    const int nprim = 8;
    const int ncoef = 3;

    Array<MultiFab,AMREX_SPACEDIM> fluxes;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        fluxes[idim].define(amrex::convert(S.boxArray(),IntVect::TheDimensionVector(idim)),
                            S.DistributionMap(), ncomp, 0);
    }

    Parm const* lparm = d_parm;

    FArrayBox qtmp, slopetmp, diff_coeff;
    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& sfab = S.array(mfi);
        auto const& dsdtfab = dSdt.array(mfi);
        AMREX_D_TERM(auto const& fxfab = fluxes[0].array(mfi);,
                     auto const& fyfab = fluxes[1].array(mfi);,
                     auto const& fzfab = fluxes[2].array(mfi););

        const Box& bxg2 = amrex::grow(bx,2);
        qtmp.resize(bxg2, nprim);
        Elixir qeli = qtmp.elixir();
        Elixir dcoeff_eli;
        auto const& q = qtmp.array();

        if (do_visc)
        {
           diff_coeff.resize(bxg2, ncoef);
           dcoeff_eli = diff_coeff.elixir();
        }

        amrex::ParallelFor(bxg2,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cns_ctoprim(i, j, k, sfab, q, *lparm);
        });

        if (do_visc)
        {
           auto const& coefs = diff_coeff.array();
           if(use_const_visc == 1 ) {
              amrex::ParallelFor(bxg2,
              [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                  cns_constcoef(i, j, k, q, coefs, *lparm);
              });
           } else {
              amrex::ParallelFor(bxg2,
              [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                  cns_diffcoef(i, j, k, q, coefs, *lparm);
              });
           }
        }


        const Box& bxg1 = amrex::grow(bx,1);
        slopetmp.resize(bxg1,neqns);
        Elixir slopeeli = slopetmp.elixir();
        auto const& slope = slopetmp.array();

        // x-direction
        int cdir = 0;
        const Box& xslpbx = amrex::grow(bx, cdir, 1);
        amrex::ParallelFor(xslpbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cns_slope_x(i, j, k, slope, q);
        });
        const Box& xflxbx = amrex::surroundingNodes(bx,cdir);
        amrex::ParallelFor(xflxbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cns_riemann_x(i, j, k, fxfab, slope, q, *lparm);
            for (int n = neqns; n < ncons; ++n) fxfab(i,j,k,n) = Real(0.0);
        });

        if (do_visc)
        {
           auto const& coefs = diff_coeff.array();
           amrex::ParallelFor(xflxbx,
           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
               cns_diff_x(i, j, k, q, coefs, dxinv, fxfab, *lparm);
           });
        }

        // y-direction
        cdir = 1;
        const Box& yslpbx = amrex::grow(bx, cdir, 1);
        amrex::ParallelFor(yslpbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cns_slope_y(i, j, k, slope, q);
        });
        const Box& yflxbx = amrex::surroundingNodes(bx,cdir);
        amrex::ParallelFor(yflxbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cns_riemann_y(i, j, k, fyfab, slope, q, *lparm);
            for (int n = neqns; n < ncons; ++n) fyfab(i,j,k,n) = Real(0.0);
        });

        if (do_visc)
        {
           auto const& coefs = diff_coeff.array();
           amrex::ParallelFor(yflxbx,
           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
               cns_diff_y(i, j, k, q, coefs, dxinv, fyfab, *lparm);
           });
        }

        // z-direction
        cdir = 2;
        const Box& zslpbx = amrex::grow(bx, cdir, 1);
        amrex::ParallelFor(zslpbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cns_slope_z(i, j, k, slope, q);
        });
        const Box& zflxbx = amrex::surroundingNodes(bx,cdir);
        amrex::ParallelFor(zflxbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cns_riemann_z(i, j, k, fzfab, slope, q, *lparm);
            for (int n = neqns; n < ncons; ++n) fzfab(i,j,k,n) = Real(0.0);
        });

        if (do_visc)
        {
           auto const& coefs = diff_coeff.array();
           amrex::ParallelFor(zflxbx,
           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
               cns_diff_z(i, j, k, q, coefs, dxinv, fzfab, *lparm);
           });
        }

        // don't have to do this, but we could
        qeli.clear(); // don't need them anymore
        slopeeli.clear();

        if (do_visc) {
           dcoeff_eli.clear();
        }

        amrex::ParallelFor(bx, ncons,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            cns_flux_to_dudt(i, j, k, n, dsdtfab, AMREX_D_DECL(fxfab,fyfab,fzfab), dxinv);
        });

        if (gravity != Real(0.0)) {
            const Real g = gravity;
            const int irho = Density;
            const int imz = Zmom;
            const int irhoE = Eden;
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dsdtfab(i,j,k,imz) += g * sfab(i,j,k,irho);
                dsdtfab(i,j,k,irhoE) += g * sfab(i,j,k,imz);
            });
        }
    }

    if (fr_as_crse) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            const Real dA = (idim == 0) ? dx[1]*dx[2] : ((idim == 1) ? dx[0]*dx[2] : dx[0]*dx[1]);
            const Real scale = -dt*dA;
            fr_as_crse->CrseInit(fluxes[idim], idim, 0, 0, NUM_STATE, scale, FluxRegister::ADD);
        }
    }

    if (fr_as_fine) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            const Real dA = (idim == 0) ? dx[1]*dx[2] : ((idim == 1) ? dx[0]*dx[2] : dx[0]*dx[1]);
            const Real scale = dt*dA;
            fr_as_fine->FineAdd(fluxes[idim], idim, 0, 0, NUM_STATE, scale);
        }
    }
}
