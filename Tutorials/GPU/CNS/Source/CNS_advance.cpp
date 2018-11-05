
#include <CNS.H>
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
//    if (do_reflux && level < parent->finestLevel()) {
//        CNS& fine_level = getLevel(level+1);
//        fr_as_crse = fine_level.flux_reg.get();
//    }

    FluxRegister* fr_as_fine = nullptr;
//    if (do_reflux && level > 0) {
//        fr_as_fine = flux_reg.get();
//    }

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
                   FluxRegister* /*fr_as_crse*/, FluxRegister* /*fr_as_fine*/)
{
    BL_PROFILE("CNS::compute_dSdt()");

    const auto dx = geom.CellSizeArray();
    const auto dxinv = geom.InvCellSizeArray();
    const int ncomp = dSdt.nComp();
    const int neqns = 5;
    const int ncons = 7;
    const int nprim = 8;

    Array<MultiFab,AMREX_SPACEDIM> fluxes;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        fluxes[idim].define(amrex::convert(S.boxArray(),IntVect::TheDimensionVector(idim)),
                            S.DistributionMap(), ncomp, 0);
    }

    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox const* sfab = &(S[mfi]);
        FArrayBox * dsdtfab = &(dSdt[mfi]);
        AMREX_D_TERM(FArrayBox* fxfab = &(fluxes[0][mfi]);,
                     FArrayBox* fyfab = &(fluxes[1][mfi]);,
                     FArrayBox* fzfab = &(fluxes[2][mfi]););

        const Box& bxg2 = amrex::grow(bx,2);
        Gpu::DeviceFab q(bxg2, nprim);
        FArrayBox* qfab = q.fabPtr();

        AMREX_LAUNCH_DEVICE_LAMBDA ( bxg2, tbx,
        {
            cns_ctoprim(tbx, *sfab, *qfab);
        });

        const Box& nbxg2 = amrex::surroundingNodes(bxg2);
        Gpu::DeviceFab slope(nbxg2,neqns);
        FArrayBox* slopefab = slope.fabPtr();

        // x-direction
        int cdir = 0;
        const Box& xbx   = amrex::surroundingNodes(bx,cdir);
        const Box& xbxg1 = amrex::grow(xbx, cdir, 1);
        AMREX_LAUNCH_DEVICE_LAMBDA ( xbxg1, tbx,
        {
            cns_slope_x(tbx, *slopefab, *qfab);
        });
        AMREX_LAUNCH_DEVICE_LAMBDA (xbx, tbx,
        {
            cns_riemann_x(tbx, *fxfab, *slopefab, *qfab);
            fxfab->setVal(0.0, tbx, neqns, (fxfab->nComp()-neqns));
        });

        // y-direction
        cdir = 1;
        const Box& ybx   = amrex::surroundingNodes(bx,cdir);
        const Box& ybxg1 = amrex::grow(ybx, cdir, 1);
        AMREX_LAUNCH_DEVICE_LAMBDA ( ybxg1, tbx,
        {
            cns_slope_y(tbx, *slopefab, *qfab);
        });
        AMREX_LAUNCH_DEVICE_LAMBDA (ybx, tbx,
        {
            cns_riemann_y(tbx, *fyfab, *slopefab, *qfab);
            fyfab->setVal(0.0, tbx, neqns, (fyfab->nComp()-neqns));
        });

        // z-direction
        cdir = 2;
        const Box& zbx   = amrex::surroundingNodes(bx,cdir);
        const Box& zbxg1 = amrex::grow(zbx, cdir, 1);
        AMREX_LAUNCH_DEVICE_LAMBDA ( zbxg1, tbx,
        {
            cns_slope_z(tbx, *slopefab, *qfab);
        });
        AMREX_LAUNCH_DEVICE_LAMBDA (zbx, tbx,
        {
            cns_riemann_z(tbx, *fzfab, *slopefab, *qfab);
            fzfab->setVal(0.0, tbx, neqns, (fzfab->nComp()-neqns));
        });

        q.clear(); // don't need them anymore
        slope.clear();

        AMREX_LAUNCH_DEVICE_LAMBDA ( bx, tbx,
        {
            cns_flux_to_dudt(tbx, *dsdtfab, AMREX_D_DECL(*fxfab,*fyfab,*fzfab), dxinv);
        });
    }
}


