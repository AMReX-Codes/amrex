#include <AMReX_MGT_Solver.H>
#include <AMReX_stencil_types.H>

#include <WarpX.H>
#include <WarpX_f.H>

namespace
{
    const std::string level_prefix {"Level_"};
}

using namespace amrex;

class NoOpPhysBC
    : public amrex::PhysBCFunctBase
{
public:
    NoOpPhysBC () {}
    virtual ~NoOpPhysBC () {}
    virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real time) override { }
    using amrex::PhysBCFunctBase::FillBoundary;
};

void
WarpX::EvolveES (int numsteps) {

    amrex::Print() << "Running in electrostatic mode \n";

    BL_PROFILE("WarpX::EvolveES()");
    Real cur_time = t_new[0];
    static int last_plot_file_step = 0;
    static int last_check_file_step = 0;

    int numsteps_max;
    if (numsteps < 0) {  // Note that the default argument is numsteps = -1
        numsteps_max = max_step;
    } else {
        numsteps_max = std::min(istep[0]+numsteps, max_step);
    }

    bool max_time_reached = false;

    // nodal storage for thee electrostatic case
    const int num_levels = max_level + 1;
    Vector<std::unique_ptr<MultiFab> > rhoNodal(num_levels);
    Vector<std::unique_ptr<MultiFab> > phiNodal(num_levels);
    Vector<std::array<std::unique_ptr<MultiFab>, 3> > eFieldNodal(num_levels);
    const int ng = 1;
    for (int lev = 0; lev <= max_level; lev++) {
        BoxArray nba = boxArray(lev);
        nba.surroundingNodes();
        rhoNodal[lev].reset(new MultiFab(nba, dmap[lev], 1, ng));
        phiNodal[lev].reset(new MultiFab(nba, dmap[lev], 1, 2));

        eFieldNodal[lev][0].reset(new MultiFab(nba, dmap[lev], 1, ng));
        eFieldNodal[lev][1].reset(new MultiFab(nba, dmap[lev], 1, ng));
        eFieldNodal[lev][2].reset(new MultiFab(nba, dmap[lev], 1, ng));
    }

    const int lev = 0;
    for (int step = istep[0]; step < numsteps_max && cur_time < stop_time; ++step)
    {

        // Start loop on time steps
        amrex::Print() << "\nSTEP " << step+1 << " starts ...\n";
#ifdef WARPX_USE_PY
        if (warpx_py_beforestep) warpx_py_beforestep();
#endif

        // At initialization, particles have p^{n-1/2} and x^{n-1/2}.
        // Beyond one step, particles have p^{n-1/2} and x^{n}.
        if (is_synchronized) {
            // on first step, push X by 0.5*dt
            mypc->PushXES(0.5*dt[lev]);
            UpdatePlasmaInjectionPosition(0.5*dt[lev]);
            mypc->Redistribute();
            mypc->DepositCharge(rhoNodal);
            computePhi(rhoNodal, phiNodal);
            computeE(eFieldNodal, phiNodal);
            is_synchronized = false;
        }

        mypc->FieldGatherES(eFieldNodal, gather_masks);

        const std::string& ppltfile = amrex::Concatenate("particles", istep[0], 5);
        auto& pc = mypc->GetParticleContainer(0);
        pc.WriteAsciiFile(ppltfile);

        // Evolve particles to p^{n+1/2} and x^{n+1}
        mypc->EvolveES(eFieldNodal, rhoNodal, cur_time, dt[lev]);

#ifdef WARPX_USE_PY
        if (warpx_py_particleinjection) warpx_py_particleinjection();
        if (warpx_py_particlescraper) warpx_py_particlescraper();
        if (warpx_py_beforedeposition) warpx_py_beforedeposition();
#endif
        mypc->DepositCharge(rhoNodal);
#ifdef WARPX_USE_PY
        if (warpx_py_afterdeposition) warpx_py_afterdeposition();

        if (warpx_py_beforeEsolve) warpx_py_beforeEsolve();
#endif
        computePhi(rhoNodal, phiNodal);
        computeE(eFieldNodal, phiNodal);
#ifdef WARPX_USE_PY
        if (warpx_py_afterEsolve) warpx_py_afterEsolve();
#endif

        if (cur_time + dt[0] >= stop_time - 1.e-3*dt[0] || step == numsteps_max-1) {
            // on last step, push by only 0.5*dt to synchronize all at n+1/2
            mypc->PushXES(-0.5*dt[lev]);
            UpdatePlasmaInjectionPosition(-0.5*dt[lev]);
            is_synchronized = true;
        }

        mypc->Redistribute();

        ++istep[0];

        cur_time += dt[0];

        bool to_make_plot = (plot_int > 0) && ((step+1) % plot_int == 0);

        amrex::Print()<< "STEP " << step+1 << " ends." << " TIME = "
                      << cur_time << " DT = " << dt[0] << "\n";

        // sync up time
        for (int i = 0; i <= finest_level; ++i) {
            t_new[i] = cur_time;
        }

        if (to_make_plot) {
            // replace with ES field Gather
            mypc->DepositCharge(rhoNodal);
            computePhi(rhoNodal, phiNodal);
            phiNodal[0]->FillBoundary(Geom(0).periodicity());
            computeE(eFieldNodal, phiNodal);
            mypc->FieldGatherES(eFieldNodal, gather_masks);
            last_plot_file_step = step+1;
            WritePlotFileES(rhoNodal, phiNodal, eFieldNodal);
        }

        if (check_int > 0 && (step+1) % check_int == 0) {
            last_check_file_step = step+1;
            WriteCheckPointFile();
        }

        if (cur_time >= stop_time - 1.e-3*dt[0]) {
            max_time_reached = true;
            break;
        }

#ifdef WARPX_USE_PY
        if (warpx_py_afterstep) warpx_py_afterstep();
#endif
        // End loop on time steps
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step && (max_time_reached || istep[0] >= max_step)) {
        WritePlotFileES(rhoNodal, phiNodal, eFieldNodal);
    }

    if (check_int > 0 && istep[0] > last_check_file_step && (max_time_reached || istep[0] >= max_step)) {
        WriteCheckPointFile();
    }
}

void WarpX::zeroOutBoundary(amrex::MultiFab& input_data,
                            amrex::MultiFab& bndry_data,
                            const FabArray<BaseFab<int> >& mask) const {
    bndry_data.setVal(0.0, 1);
    for (MFIter mfi(input_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        WRPX_ZERO_OUT_BNDRY(bx.loVect(), bx.hiVect(),
                            input_data[mfi].dataPtr(),
                            bndry_data[mfi].dataPtr(),
                            mask[mfi].dataPtr());
    }
    bndry_data.FillBoundary();
}

void
WarpX::fixRHSForSolve(Vector<std::unique_ptr<MultiFab> >& rhs,
                      const Vector<std::unique_ptr<FabArray<BaseFab<int> > > >& masks) const {
    int num_levels = rhs.size();
    for (int lev = 0; lev < num_levels; ++lev) {
        MultiFab& fine_rhs = *rhs[lev];
        const FabArray<BaseFab<int> >& mask = *masks[lev];
        const BoxArray& fine_ba = fine_rhs.boxArray();
        const DistributionMapping& fine_dm = fine_rhs.DistributionMap();
        MultiFab fine_bndry_data(fine_ba, fine_dm, 1, 1);
        zeroOutBoundary(fine_rhs, fine_bndry_data, mask);
    }
}

void WarpX::getLevelMasks(Vector<std::unique_ptr<FabArray<BaseFab<int> > > >& masks,
                          const int nnodes) {
    int num_levels = grids.size();
    BL_ASSERT(num_levels == dmap.size());

    int covered = 0;
    int notcovered = 1;
    int physbnd = 1;
    int interior = 0;

    for (int lev = 0; lev < num_levels; ++lev) {
        BoxArray nba = grids[lev];
        nba.surroundingNodes();

        FabArray<BaseFab<int> > tmp_mask(nba, dmap[lev], 1, nnodes);
        tmp_mask.BuildMask(geom[lev].Domain(), geom[lev].periodicity(),
                           covered, notcovered, physbnd, interior);
        masks[lev].reset(new FabArray<BaseFab<int> >(nba, dmap[lev], 1, 0));
        for (MFIter mfi(tmp_mask); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            WRPX_BUILD_MASK(bx.loVect(), bx.hiVect(),
                            tmp_mask[mfi].dataPtr(), (*masks[lev])[mfi].dataPtr(), &nnodes);
        }
    }
}


void WarpX::computePhi(const Vector<std::unique_ptr<MultiFab> >& rho,
                             Vector<std::unique_ptr<MultiFab> >& phi) const {


    int num_levels = rho.size();
    Vector<std::unique_ptr<MultiFab> > rhs(num_levels);
    for (int lev = 0; lev < num_levels; ++lev) {
        phi[lev]->setVal(0.0, 2);
        rhs[lev].reset(new MultiFab(rho[lev]->boxArray(), dmap[lev], 1, 0));
        MultiFab::Copy(*rhs[lev], *rho[lev], 0, 0, 1, 0);
        rhs[lev]->mult(-1.0/PhysConst::ep0, 0);
    }

    fixRHSForSolve(rhs, masks);

    bool nodal = true;
    bool have_rhcc = false;
    int  nc = 0;
    int Ncomp = 1;
    int stencil = ND_CROSS_STENCIL;
    int verbose = 0;
    Vector<int> mg_bc(2*AMREX_SPACEDIM, 1); // this means Dirichlet
    Real rel_tol = 1.0e-14;
    Real abs_tol = 1.0e-14;

    Vector<Geometry>            level_geom(1);
    Vector<BoxArray>            level_grids(1);
    Vector<DistributionMapping> level_dm(1);
    Vector<MultiFab*>           level_phi(1);
    Vector<MultiFab*>           level_rhs(1);

    for (int lev = 0; lev < num_levels; ++lev) {
        level_phi[0]   = phi[lev].get();
        level_rhs[0]   = rhs[lev].get();
        level_geom[0]  = geom[lev];
        level_grids[0] = grids[lev];
        level_dm[0]    = dmap[lev];

        MGT_Solver solver(level_geom, mg_bc.dataPtr(), level_grids,
                          level_dm, nodal,
                          stencil, have_rhcc, nc, Ncomp, verbose);

        solver.set_nodal_const_coefficients(1.0);

        solver.solve_nodal(level_phi, level_rhs, rel_tol, abs_tol);

        if (lev < num_levels-1) {

            NoOpPhysBC cphysbc, fphysbc;
#if AMREX_SPACEDIM == 3
            int lo_bc[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
            int hi_bc[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
#else
            int lo_bc[] = {BCType::int_dir, BCType::int_dir};
            int hi_bc[] = {BCType::int_dir, BCType::int_dir};
#endif
            Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));
            NodeBilinear mapper;

            amrex::InterpFromCoarseLevel(*phi[lev+1], 0.0, *phi[lev],
                                         0, 0, 1, geom[lev], geom[lev+1],
                                         cphysbc, fphysbc,
                                         IntVect(AMREX_D_DECL(2, 2, 2)), &mapper, bcs);
        }
    }

    for (int lev = 0; lev < num_levels; ++lev) {
        const Geometry& gm = geom[lev];
        phi[lev]->FillBoundary(gm.periodicity());
    }
}

void WarpX::computeE(Vector<std::array<std::unique_ptr<MultiFab>, 3> >& E,
                     const Vector<std::unique_ptr<MultiFab> >& phi) const {

    const int num_levels = E.size();
    for (int lev = 0; lev < num_levels; ++lev) {
        const auto& gm = GetInstance().Geom(lev);
        const Real* dx = gm.CellSize();
        for (MFIter mfi(*phi[lev]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();

            WRPX_COMPUTE_E_NODAL(bx.loVect(), bx.hiVect(),
                                 (*phi[lev] )[mfi].dataPtr(),
                                 (*E[lev][0])[mfi].dataPtr(),
                                 (*E[lev][1])[mfi].dataPtr(),
#if AMREX_SPACEDIM == 3
                                 (*E[lev][2])[mfi].dataPtr(),
#endif
                                 dx);
        }

        E[lev][0]->FillBoundary(gm.periodicity());
        E[lev][1]->FillBoundary(gm.periodicity());
#if AMREX_SPACEDIM == 3
        E[lev][2]->FillBoundary(gm.periodicity());
#endif
    }
}
