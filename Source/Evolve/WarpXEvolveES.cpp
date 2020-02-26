/* Copyright 2019 Andrew Myers, Axel Huebl, David Bizzozero
 * David Grote, Maxence Thevenet, Remi Lehe
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "FortranInterface/WarpX_f.H"


namespace
{
    const std::string level_prefix {"Level_"};
}

#ifdef WARPX_DO_ELECTROSTATIC
using namespace amrex;

void
WarpX::EvolveES (int numsteps) {

    amrex::Print() << "Running in electrostatic mode \n";

    WARPX_PROFILE("WarpX::EvolveES()");
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
            mypc->PushX(0.5*dt[lev]);
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
            mypc->PushX(-0.5*dt[lev]);
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
#endif // WARPX_DO_ELECTROSTATIC
