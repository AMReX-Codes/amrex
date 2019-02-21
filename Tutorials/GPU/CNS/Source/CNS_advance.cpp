
#include "CNS.H"
#include "CNS_hydro_K.H"
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

    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab dSdt(grids,dmap,NUM_STATE,0,MFInfo(),Factory());
    MultiFab Sborder(grids,dmap,NUM_STATE,NUM_GROW,MFInfo(),Factory());

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
        fr_as_crse->setVal(0.0);
    }

    // RK2 stage 1
    FillPatch(*this, Sborder, NUM_GROW, time, State_Type, 0, NUM_STATE);
    compute_dSdt(Sborder, dSdt, 0.5*dt, fr_as_crse, fr_as_fine);
    // U^* = U^n + dt*dUdt^n
    MultiFab::LinComb(S_new, 1.0, Sborder, 0, dt, dSdt, 0, 0, NUM_STATE, 0);
    computeTemp(S_new,0);

    // RK2 stage 2
    // After fillpatch Sborder = U^n+dt*dUdt^n
    FillPatch(*this, Sborder, NUM_GROW, time+dt, State_Type, 0, NUM_STATE);
    compute_dSdt(Sborder, dSdt, 0.5*dt, fr_as_crse, fr_as_fine);
    // S_new = 0.5*(Sborder+S_old) = U^n + 0.5*dt*dUdt^n
    MultiFab::LinComb(S_new, 0.5, Sborder, 0, 0.5, S_old, 0, 0, NUM_STATE, 0);
    // S_new += 0.5*dt*dSdt
    MultiFab::Saxpy(S_new, 0.5*dt, dSdt, 0, 0, NUM_STATE, 0);
    // We now have S_new = U^{n+1} = (U^n+0.5*dt*dUdt^n) + 0.5*dt*dUdt^*
    computeTemp(S_new,0);
    
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

    Array<MultiFab,AMREX_SPACEDIM> fluxes;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        fluxes[idim].define(amrex::convert(S.boxArray(),IntVect::TheDimensionVector(idim)),
                            S.DistributionMap(), ncomp, 0);
    }

    FArrayBox qtmp, slopetmp;
    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& sfab = S.array(mfi);
        auto const& dsdtfab = dSdt.array(mfi);
        AMREX_D_TERM(auto const& fxfab = fluxes[0].array(mfi);,
                     auto const& fyfab = fluxes[1].array(mfi);,
                     auto const& fzfab = fluxes[2].array(mfi););

        const Box& bxg2 = amrex::grow(bx,2);
        AsyncFab q(qtmp,bxg2, nprim);
        auto const& qfab = q.array();

        amrex::ParallelFor(bxg2,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            cns_ctoprim(i, j, k, sfab, qfab);
        });

        const Box& bxg1 = amrex::grow(bx,1);
        AsyncFab slope(slopetmp,bxg1,neqns);
        auto const& slopefab = slope.array();

        // x-direction
        int cdir = 0;
        const Box& xslpbx = amrex::grow(bx, cdir, 1);
        amrex::ParallelFor(xslpbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            cns_slope_x(i, j, k, slopefab, qfab);
        });
        const Box& xflxbx = amrex::surroundingNodes(bx,cdir);
        amrex::ParallelFor(xflxbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            cns_riemann_x(i, j, k, fxfab, slopefab, qfab);
            for (int n = neqns; n < ncons; ++n) fxfab(i,j,k,n) = 0.0;
        });

        // y-direction
        cdir = 1;
        const Box& yslpbx = amrex::grow(bx, cdir, 1);
        amrex::ParallelFor(yslpbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            cns_slope_y(i, j, k, slopefab, qfab);
        });
        const Box& yflxbx = amrex::surroundingNodes(bx,cdir);
        amrex::ParallelFor(yflxbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            cns_riemann_y(i, j, k, fyfab, slopefab, qfab);
            for (int n = neqns; n < ncons; ++n) fyfab(i,j,k,n) = 0.0;
        });

        // z-direction
        cdir = 2;
        const Box& zslpbx = amrex::grow(bx, cdir, 1);
        amrex::ParallelFor(zslpbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            cns_slope_z(i, j, k, slopefab, qfab);
        });
        const Box& zflxbx = amrex::surroundingNodes(bx,cdir);
        amrex::ParallelFor(zflxbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            cns_riemann_z(i, j, k, fzfab, slopefab, qfab);
            for (int n = neqns; n < ncons; ++n) fzfab(i,j,k,n) = 0.0;
        });

        // don't have to do this, but we could
        q.clear(); // don't need them anymore
        slope.clear();

        amrex::ParallelFor(bx, ncons,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
        {
            cns_flux_to_dudt(i, j, k, n, dsdtfab, AMREX_D_DECL(fxfab,fyfab,fzfab), dxinv);
        });

        if (gravity != 0.0) {
            const Real g = gravity;
            const int irho = Density;
            const int imz = Zmom;
            const int irhoE = Eden;
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
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


