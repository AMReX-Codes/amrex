
#include <cmath>
#include <limits>

#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpX_f.H>
#ifdef WARPX_USE_PY
#include <WarpX_py.H>
#endif

#ifdef BL_USE_SENSEI_INSITU
#include <AMReX_AmrMeshInSituBridge.H>
#endif

using namespace amrex;

void
WarpX::Evolve (int numsteps) {
    BL_PROFILE_REGION("WarpX::Evolve()");

#ifdef WARPX_DO_ELECTROSTATIC
    if (do_electrostatic) {
        EvolveES(numsteps);
    } else {
      EvolveEM(numsteps);
    }
#else
    EvolveEM(numsteps);
#endif // WARPX_DO_ELECTROSTATIC

}

void
WarpX::EvolveEM (int numsteps)
{
    BL_PROFILE("WarpX::EvolveEM()");

    Real cur_time = t_new[0];
    static int last_plot_file_step = 0;
    static int last_check_file_step = 0;
    static int last_insitu_step = 0;

    int numsteps_max;
    if (numsteps < 0) {  // Note that the default argument is numsteps = -1
        numsteps_max = max_step;
    } else {
        numsteps_max = std::min(istep[0]+numsteps, max_step);
    }

    bool max_time_reached = false;
    Real walltime, walltime_start = amrex::second();
    for (int step = istep[0]; step < numsteps_max && cur_time < stop_time; ++step)
    {
        Real walltime_beg_step = amrex::second();

	// Start loop on time steps
        amrex::Print() << "\nSTEP " << step+1 << " starts ...\n";
#ifdef WARPX_USE_PY
        if (warpx_py_beforestep) warpx_py_beforestep();
#endif

        if (costs[0] != nullptr)
        {
#ifdef WARPX_USE_PSATD
            amrex::Abort("LoadBalance for PSATD: TODO");
#endif

            if (step > 0 && (step+1) % load_balance_int == 0)
            {
                LoadBalance();
		// Reset the costs to 0
		for (int lev = 0; lev <= finest_level; ++lev) {
		  costs[lev]->setVal(0.0);
		}
            }

            for (int lev = 0; lev <= finest_level; ++lev) {
	      // Perform running average of the costs
	      // (Giving more importance to most recent costs)
	      (*costs[lev].get()).mult( (1. - 2./load_balance_int) );
            }
        }

        // At the beginning, we have B^{n} and E^{n}.
        // Particles have p^{n} and x^{n}.
        // is_synchronized is true.
        if (is_synchronized) {
            FillBoundaryE();
            FillBoundaryB();
            UpdateAuxilaryData();
            // on first step, push p by -0.5*dt
            for (int lev = 0; lev <= finest_level; ++lev) {
                mypc->PushP(lev, -0.5*dt[lev],
                            *Efield_aux[lev][0],*Efield_aux[lev][1],*Efield_aux[lev][2],
                            *Bfield_aux[lev][0],*Bfield_aux[lev][1],*Bfield_aux[lev][2]);
            }
            is_synchronized = false;
        } else {
           // Beyond one step, we have E^{n} and B^{n}.
           // Particles have p^{n-1/2} and x^{n}.
            FillBoundaryE();
            FillBoundaryB();
            UpdateAuxilaryData();
        }

        if (do_subcycling == 0 || finest_level == 0) {
            OneStep_nosub(cur_time);
        } else if (do_subcycling == 1 && finest_level == 1) {
            OneStep_sub1(cur_time);
        } else {
            amrex::Print() << "Error: do_subcycling = " << do_subcycling << std::endl;
            amrex::Abort("Unsupported do_subcycling type");
        }

#ifdef WARPX_USE_PY
        if (warpx_py_beforeEsolve) warpx_py_beforeEsolve();
#endif
        if (cur_time + dt[0] >= stop_time - 1.e-3*dt[0] || step == numsteps_max-1) {
            // At the end of last step, push p by 0.5*dt to synchronize
            UpdateAuxilaryData();
            for (int lev = 0; lev <= finest_level; ++lev) {
                mypc->PushP(lev, 0.5*dt[lev],
                    *Efield_aux[lev][0],*Efield_aux[lev][1],*Efield_aux[lev][2],
                    *Bfield_aux[lev][0],*Bfield_aux[lev][1],*Bfield_aux[lev][2]);
            }
            is_synchronized = true;
        }
#ifdef WARPX_USE_PY
        if (warpx_py_afterEsolve) warpx_py_afterEsolve();
#endif

        for (int lev = 0; lev <= max_level; ++lev) {
            ++istep[lev];
        }

	cur_time += dt[0];

        bool to_make_plot = (plot_int > 0) && ((step+1) % plot_int == 0);

        bool do_insitu = ((step+1) >= insitu_start) &&
             (insitu_int > 0) && ((step+1) % insitu_int == 0);

        bool move_j = is_synchronized || to_make_plot || do_insitu;
        // If is_synchronized we need to shift j too so that next step we can evolve E by dt/2.
        // We might need to move j because we are going to make a plotfile.

	int num_moved = MoveWindow(move_j);
        
        if (max_level == 0) {
            int num_redistribute_ghost = num_moved + 1;
            mypc->RedistributeLocal(num_redistribute_ghost);
        }
        else {
            mypc->Redistribute();
        }

	bool to_sort = (sort_int > 0) && ((step+1) % sort_int == 0);
	if (to_sort) {
	    amrex::Print() << "re-sorting particles \n";
	    mypc->SortParticlesByCell();
	}

        amrex::Print()<< "STEP " << step+1 << " ends." << " TIME = " << cur_time
                      << " DT = " << dt[0] << "\n";
        Real walltime_end_step = amrex::second();
        walltime = walltime_end_step - walltime_start;
        amrex::Print()<< "Walltime = " << walltime
                      << " s; This step = " << walltime_end_step-walltime_beg_step
                      << " s; Avg. per step = " << walltime/(step+1) << " s\n";

	// sync up time
	for (int i = 0; i <= max_level; ++i) {
	    t_new[i] = cur_time;
	}

        if (do_boosted_frame_diagnostic) {
            std::unique_ptr<MultiFab> cell_centered_data = nullptr;
            if (WarpX::do_boosted_frame_fields) {
                cell_centered_data = GetCellCenteredData();
            }
            myBFD->writeLabFrameData(cell_centered_data.get(), *mypc, geom[0], cur_time, dt[0]);
        }

	if (to_make_plot || do_insitu)
        {
            FillBoundaryE();
            FillBoundaryB();
            UpdateAuxilaryData();

            for (int lev = 0; lev <= finest_level; ++lev) {
                mypc->FieldGather(lev,
                                  *Efield_aux[lev][0],*Efield_aux[lev][1],*Efield_aux[lev][2],
                                  *Bfield_aux[lev][0],*Bfield_aux[lev][1],*Bfield_aux[lev][2]);
            }

	    last_plot_file_step = step+1;
	    last_insitu_step = step+1;

        if (to_make_plot)
    	    WritePlotFile();

        if (do_insitu)
            UpdateInSitu();
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

    bool write_plot_file = plot_int > 0 && istep[0] > last_plot_file_step && (max_time_reached || istep[0] >= max_step);

    bool do_insitu = (insitu_start >= istep[0]) && (insitu_int > 0) &&
        (istep[0] > last_insitu_step) && (max_time_reached || istep[0] >= max_step);

    if (write_plot_file || do_insitu)
    {
        FillBoundaryE();
        FillBoundaryB();
        UpdateAuxilaryData();

        for (int lev = 0; lev <= finest_level; ++lev) {
            mypc->FieldGather(lev,
                              *Efield_aux[lev][0],*Efield_aux[lev][1],*Efield_aux[lev][2],
                              *Bfield_aux[lev][0],*Bfield_aux[lev][1],*Bfield_aux[lev][2]);
        }

        if (write_plot_file)
            WritePlotFile();

        if (do_insitu)
            UpdateInSitu();
    }

    if (check_int > 0 && istep[0] > last_check_file_step && (max_time_reached || istep[0] >= max_step)) {
	WriteCheckPointFile();
    }

    if (do_boosted_frame_diagnostic) {
        myBFD->Flush(geom[0]);
    }

#ifdef BL_USE_SENSEI_INSITU
    insitu_bridge->finalize();
#endif
}

/* /brief Perform one PIC iteration, without subcycling
*  i.e. all levels/patches use the same timestep (that of the finest level)
*  for the field advance and particle pusher.
*/
void
WarpX::OneStep_nosub (Real cur_time)
{
    // Push particle from x^{n} to x^{n+1}
    //               from p^{n-1/2} to p^{n+1/2}
    // Deposit current j^{n+1/2}
    // Deposit charge density rho^{n}
#ifdef WARPX_USE_PY
    if (warpx_py_particleinjection) warpx_py_particleinjection();
    if (warpx_py_particlescraper) warpx_py_particlescraper();
    if (warpx_py_beforedeposition) warpx_py_beforedeposition();
#endif
    PushParticlesandDepose(cur_time);
#ifdef WARPX_USE_PY
    if (warpx_py_afterdeposition) warpx_py_afterdeposition();
#endif

    SyncCurrent();

    SyncRho(rho_fp, rho_cp);

    // Push E and B from {n} to {n+1}
    // (And update guard cells immediately afterwards)
#ifdef WARPX_USE_PSATD
    PushPSATD(dt[0]);
    FillBoundaryE();
    FillBoundaryB();
#else
    EvolveF(0.5*dt[0], DtType::FirstHalf);
    FillBoundaryF();
    EvolveB(0.5*dt[0]); // We now have B^{n+1/2}
    FillBoundaryB();
    EvolveE(dt[0]); // We now have E^{n+1}
    FillBoundaryE();
    EvolveF(0.5*dt[0], DtType::SecondHalf);
    EvolveB(0.5*dt[0]); // We now have B^{n+1}
    if (do_pml) {
        DampPML();
        FillBoundaryE();
    }
    FillBoundaryB();
#endif
}

/* /brief Perform one PIC iteration, with subcycling
*  i.e. The fine patch uses a smaller timestep (and steps more often)
*  than the coarse patch, for the field advance and particle pusher.
*
* This version of subcycling only works for 2 levels and with a refinement
* ratio of 2.
* The particles and fields of the fine patch are pushed twice
* (with dt[coarse]/2) in this routine.
* The particles of the coarse patch and mother grid are pushed only once
* (with dt[coarse]). The fields on the coarse patch and mother grid
* are pushed in a way which is equivalent to pushing once only, with
* a current which is the average of the coarse + fine current at the 2
* steps of the fine grid.
*
*/
void
WarpX::OneStep_sub1 (Real curtime)
{
    // TODO: we could save some charge depositions

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(finest_level == 1, "Must have exactly two levels");
    const int fine_lev = 1;
    const int coarse_lev = 0;

    // i) Push particles and fields on the fine patch (first fine step)
    PushParticlesandDepose(fine_lev, curtime);
    RestrictCurrentFromFineToCoarsePatch(fine_lev);
    RestrictRhoFromFineToCoarsePatch(fine_lev);
    ApplyFilterandSumBoundaryJ(fine_lev, PatchType::fine);
    NodalSyncJ(fine_lev, PatchType::fine);
    ApplyFilterandSumBoundaryRho(fine_lev, PatchType::fine, 0, 2);
    NodalSyncRho(fine_lev, PatchType::fine, 0, 2);

    EvolveB(fine_lev, PatchType::fine, 0.5*dt[fine_lev]);
    EvolveF(fine_lev, PatchType::fine, 0.5*dt[fine_lev], DtType::FirstHalf);
    FillBoundaryB(fine_lev, PatchType::fine);
    FillBoundaryF(fine_lev, PatchType::fine);

    EvolveE(fine_lev, PatchType::fine, dt[fine_lev]);
    FillBoundaryE(fine_lev, PatchType::fine);

    EvolveB(fine_lev, PatchType::fine, 0.5*dt[fine_lev]);
    EvolveF(fine_lev, PatchType::fine, 0.5*dt[fine_lev], DtType::SecondHalf);

    if (do_pml) {
        DampPML(fine_lev, PatchType::fine);
        FillBoundaryE(fine_lev, PatchType::fine);
    }

    FillBoundaryB(fine_lev, PatchType::fine);

    // ii) Push particles on the coarse patch and mother grid.
    // Push the fields on the coarse patch and mother grid
    // by only half a coarse step (first half)
    PushParticlesandDepose(coarse_lev, curtime);
    StoreCurrent(coarse_lev);
    AddCurrentFromFineLevelandSumBoundary(coarse_lev);
    AddRhoFromFineLevelandSumBoundary(coarse_lev, 0, 1);

    EvolveB(fine_lev, PatchType::coarse, dt[fine_lev]);
    EvolveF(fine_lev, PatchType::coarse, dt[fine_lev], DtType::FirstHalf);
    FillBoundaryB(fine_lev, PatchType::coarse);
    FillBoundaryF(fine_lev, PatchType::coarse);

    EvolveE(fine_lev, PatchType::coarse, dt[fine_lev]);
    FillBoundaryE(fine_lev, PatchType::coarse);

    EvolveB(coarse_lev, PatchType::fine, 0.5*dt[coarse_lev]);
    EvolveF(coarse_lev, PatchType::fine, 0.5*dt[coarse_lev], DtType::FirstHalf);
    FillBoundaryB(coarse_lev, PatchType::fine);
    FillBoundaryF(coarse_lev, PatchType::fine);

    EvolveE(coarse_lev, PatchType::fine, 0.5*dt[coarse_lev]);
    FillBoundaryE(coarse_lev, PatchType::fine);

    // iii) Get auxiliary fields on the fine grid, at dt[fine_lev]
    UpdateAuxilaryData();

    // iv) Push particles and fields on the fine patch (second fine step)
    PushParticlesandDepose(fine_lev, curtime+dt[fine_lev]);
    RestrictCurrentFromFineToCoarsePatch(fine_lev);
    RestrictRhoFromFineToCoarsePatch(fine_lev);
    ApplyFilterandSumBoundaryJ(fine_lev, PatchType::fine);
    NodalSyncJ(fine_lev, PatchType::fine);
    ApplyFilterandSumBoundaryRho(fine_lev, PatchType::fine, 0, 2);
    NodalSyncRho(fine_lev, PatchType::fine, 0, 2);

    EvolveB(fine_lev, PatchType::fine, 0.5*dt[fine_lev]);
    EvolveF(fine_lev, PatchType::fine, 0.5*dt[fine_lev], DtType::FirstHalf);
    FillBoundaryB(fine_lev, PatchType::fine);
    FillBoundaryF(fine_lev, PatchType::fine);

    EvolveE(fine_lev, PatchType::fine, dt[fine_lev]);
    FillBoundaryE(fine_lev, PatchType::fine);

    EvolveB(fine_lev, PatchType::fine, 0.5*dt[fine_lev]);
    EvolveF(fine_lev, PatchType::fine, 0.5*dt[fine_lev], DtType::SecondHalf);

    if (do_pml) {
        DampPML(fine_lev, PatchType::fine);
        FillBoundaryE(fine_lev, PatchType::fine);
    }

    FillBoundaryB(fine_lev, PatchType::fine);
    FillBoundaryF(fine_lev, PatchType::fine);

    // v) Push the fields on the coarse patch and mother grid
    // by only half a coarse step (second half)
    RestoreCurrent(coarse_lev);
    AddCurrentFromFineLevelandSumBoundary(coarse_lev);
    AddRhoFromFineLevelandSumBoundary(coarse_lev, 1, 1);

    EvolveE(fine_lev, PatchType::coarse, dt[fine_lev]);
    FillBoundaryE(fine_lev, PatchType::coarse);

    EvolveB(fine_lev, PatchType::coarse, dt[fine_lev]);
    EvolveF(fine_lev, PatchType::coarse, dt[fine_lev], DtType::SecondHalf);

    if (do_pml) {
        DampPML(fine_lev, PatchType::coarse); // do it twice
        DampPML(fine_lev, PatchType::coarse);
        FillBoundaryE(fine_lev, PatchType::coarse);
    }

    FillBoundaryB(fine_lev, PatchType::coarse);
    FillBoundaryF(fine_lev, PatchType::coarse);

    EvolveE(coarse_lev, PatchType::fine, 0.5*dt[coarse_lev]);
    FillBoundaryE(coarse_lev, PatchType::fine);

    EvolveB(coarse_lev, PatchType::fine, 0.5*dt[coarse_lev]);
    EvolveF(coarse_lev, PatchType::fine, 0.5*dt[coarse_lev], DtType::SecondHalf);

    if (do_pml) {
        DampPML(coarse_lev, PatchType::fine);
        FillBoundaryE(coarse_lev, PatchType::fine);
    }

    FillBoundaryB(coarse_lev, PatchType::fine);
}

void
WarpX::EvolveB (Real dt)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        EvolveB(lev, dt);
    }
}

void
WarpX::EvolveB (int lev, Real dt)
{
    BL_PROFILE("WarpX::EvolveB()");
    EvolveB(lev, PatchType::fine, dt);
    if (lev > 0)
    {
        EvolveB(lev, PatchType::coarse, dt);
    }
}

void
WarpX::EvolveB (int lev, PatchType patch_type, amrex::Real dt)
{
    const int patch_level = (patch_type == PatchType::fine) ? lev : lev-1;
    const std::array<Real,3>& dx = WarpX::CellSize(patch_level);
    const std::array<Real,3> dtsdx {dt/dx[0], dt/dx[1], dt/dx[2]};

    MultiFab *Ex, *Ey, *Ez, *Bx, *By, *Bz;
    if (patch_type == PatchType::fine)
    {
        Ex = Efield_fp[lev][0].get();
        Ey = Efield_fp[lev][1].get();
        Ez = Efield_fp[lev][2].get();
        Bx = Bfield_fp[lev][0].get();
        By = Bfield_fp[lev][1].get();
        Bz = Bfield_fp[lev][2].get();
    }
    else
    {
        Ex = Efield_cp[lev][0].get();
        Ey = Efield_cp[lev][1].get();
        Ez = Efield_cp[lev][2].get();
        Bx = Bfield_cp[lev][0].get();
        By = Bfield_cp[lev][1].get();
        Bz = Bfield_cp[lev][2].get();
    }

    MultiFab* cost = costs[lev].get();
    const IntVect& rr = (lev > 0) ? refRatio(lev-1) : IntVect::TheUnitVector();

    // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Bx, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Real wt = amrex::second();

        const Box& tbx  = mfi.tilebox(Bx_nodal_flag);
        const Box& tby  = mfi.tilebox(By_nodal_flag);
        const Box& tbz  = mfi.tilebox(Bz_nodal_flag);

        // Call picsar routine for each tile
        warpx_push_bvec(
		      tbx.loVect(), tbx.hiVect(),
		      tby.loVect(), tby.hiVect(),
		      tbz.loVect(), tbz.hiVect(),
		      BL_TO_FORTRAN_3D((*Ex)[mfi]),
		      BL_TO_FORTRAN_3D((*Ey)[mfi]),
		      BL_TO_FORTRAN_3D((*Ez)[mfi]),
		      BL_TO_FORTRAN_3D((*Bx)[mfi]),
		      BL_TO_FORTRAN_3D((*By)[mfi]),
		      BL_TO_FORTRAN_3D((*Bz)[mfi]),
		      &dtsdx[0], &dtsdx[1], &dtsdx[2],
		      &WarpX::maxwell_fdtd_solver_id);

        if (cost) {
            Box cbx = mfi.tilebox(IntVect{AMREX_D_DECL(0,0,0)});
            if (patch_type == PatchType::coarse) cbx.refine(rr);
            wt = (amrex::second() - wt) / cbx.d_numPts();\
            FArrayBox* costfab = cost->fabPtr(mfi);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( cbx, work_box,
            {
                costfab->plus(wt, work_box);
            });
        }
    }

    if (do_pml && pml[lev]->ok())
    {
        const auto& pml_B = (patch_type == PatchType::fine) ? pml[lev]->GetB_fp() : pml[lev]->GetB_cp();
        const auto& pml_E = (patch_type == PatchType::fine) ? pml[lev]->GetE_fp() : pml[lev]->GetE_cp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*pml_B[0], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const Box& tbx  = mfi.tilebox(Bx_nodal_flag);
            const Box& tby  = mfi.tilebox(By_nodal_flag);
            const Box& tbz  = mfi.tilebox(Bz_nodal_flag);

            WRPX_PUSH_PML_BVEC(
			     tbx.loVect(), tbx.hiVect(),
			     tby.loVect(), tby.hiVect(),
			     tbz.loVect(), tbz.hiVect(),
			     BL_TO_FORTRAN_3D((*pml_E[0])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_E[1])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_E[2])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_B[0])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_B[1])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_B[2])[mfi]),
			     &dtsdx[0], &dtsdx[1], &dtsdx[2],
			     &WarpX::maxwell_fdtd_solver_id);
        }
    }
}

void
WarpX::EvolveE (Real dt)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        EvolveE(lev, dt);
    }
}

void
WarpX::EvolveE (int lev, Real dt)
{
    BL_PROFILE("WarpX::EvolveE()");
    EvolveE(lev, PatchType::fine, dt);
    if (lev > 0)
    {
        EvolveE(lev, PatchType::coarse, dt);
    }
}

void
WarpX::EvolveE (int lev, PatchType patch_type, amrex::Real dt)
{
    const Real mu_c2_dt = (PhysConst::mu0*PhysConst::c*PhysConst::c) * dt;
    const Real c2dt = (PhysConst::c*PhysConst::c) * dt;

    int patch_level = (patch_type == PatchType::fine) ? lev : lev-1;
    const std::array<Real,3>& dx = WarpX::CellSize(patch_level);
    const std::array<Real,3> dtsdx_c2 {c2dt/dx[0], c2dt/dx[1], c2dt/dx[2]};

    MultiFab *Ex, *Ey, *Ez, *Bx, *By, *Bz, *jx, *jy, *jz, *F;
    if (patch_type == PatchType::fine)
    {
        Ex = Efield_fp[lev][0].get();
        Ey = Efield_fp[lev][1].get();
        Ez = Efield_fp[lev][2].get();
        Bx = Bfield_fp[lev][0].get();
        By = Bfield_fp[lev][1].get();
        Bz = Bfield_fp[lev][2].get();
        jx = current_fp[lev][0].get();
        jy = current_fp[lev][1].get();
        jz = current_fp[lev][2].get();
        F  = F_fp[lev].get();
    }
    else if (patch_type == PatchType::coarse)
    {
        Ex = Efield_cp[lev][0].get();
        Ey = Efield_cp[lev][1].get();
        Ez = Efield_cp[lev][2].get();
        Bx = Bfield_cp[lev][0].get();
        By = Bfield_cp[lev][1].get();
        Bz = Bfield_cp[lev][2].get();
        jx = current_cp[lev][0].get();
        jy = current_cp[lev][1].get();
        jz = current_cp[lev][2].get();
        F  = F_cp[lev].get();
    }

    MultiFab* cost = costs[lev].get();
    const IntVect& rr = (lev > 0) ? refRatio(lev-1) : IntVect::TheUnitVector();

    // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(*Ex, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        Real wt = amrex::second();

        const Box& tex  = mfi.tilebox(Ex_nodal_flag);
        const Box& tey  = mfi.tilebox(Ey_nodal_flag);
        const Box& tez  = mfi.tilebox(Ez_nodal_flag);

        // Call picsar routine for each tile
        warpx_push_evec(
		      tex.loVect(), tex.hiVect(),
		      tey.loVect(), tey.hiVect(),
		      tez.loVect(), tez.hiVect(),
		      BL_TO_FORTRAN_3D((*Ex)[mfi]),
		      BL_TO_FORTRAN_3D((*Ey)[mfi]),
		      BL_TO_FORTRAN_3D((*Ez)[mfi]),
		      BL_TO_FORTRAN_3D((*Bx)[mfi]),
		      BL_TO_FORTRAN_3D((*By)[mfi]),
		      BL_TO_FORTRAN_3D((*Bz)[mfi]),
		      BL_TO_FORTRAN_3D((*jx)[mfi]),
		      BL_TO_FORTRAN_3D((*jy)[mfi]),
		      BL_TO_FORTRAN_3D((*jz)[mfi]),
		      &mu_c2_dt,
		      &dtsdx_c2[0], &dtsdx_c2[1], &dtsdx_c2[2]);

        if (F)
        {
            warpx_push_evec_f(
			  tex.loVect(), tex.hiVect(),
			  tey.loVect(), tey.hiVect(),
			  tez.loVect(), tez.hiVect(),
			  BL_TO_FORTRAN_3D((*Ex)[mfi]),
			  BL_TO_FORTRAN_3D((*Ey)[mfi]),
			  BL_TO_FORTRAN_3D((*Ez)[mfi]),
			  BL_TO_FORTRAN_3D((*F)[mfi]),
			  &dtsdx_c2[0], &dtsdx_c2[1], &dtsdx_c2[2],
			  &WarpX::maxwell_fdtd_solver_id);
        }

        if (cost) {
            Box cbx = mfi.tilebox(IntVect{AMREX_D_DECL(0,0,0)});
            if (patch_type == PatchType::coarse) cbx.refine(rr);
            wt = (amrex::second() - wt) / cbx.d_numPts();
            FArrayBox* costfab = cost->fabPtr(mfi);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( cbx, work_box,
            {
                costfab->plus(wt, work_box);
            });
        }
    }

    if (do_pml && pml[lev]->ok())
    {
        if (F) pml[lev]->ExchangeF(patch_type, F);

        const auto& pml_B = (patch_type == PatchType::fine) ? pml[lev]->GetB_fp() : pml[lev]->GetB_cp();
        const auto& pml_E = (patch_type == PatchType::fine) ? pml[lev]->GetE_fp() : pml[lev]->GetE_cp();
        const auto& pml_F = (patch_type == PatchType::fine) ? pml[lev]->GetF_fp() : pml[lev]->GetF_cp();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*pml_E[0], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const Box& tex  = mfi.tilebox(Ex_nodal_flag);
            const Box& tey  = mfi.tilebox(Ey_nodal_flag);
            const Box& tez  = mfi.tilebox(Ez_nodal_flag);

            WRPX_PUSH_PML_EVEC(
			     tex.loVect(), tex.hiVect(),
			     tey.loVect(), tey.hiVect(),
			     tez.loVect(), tez.hiVect(),
			     BL_TO_FORTRAN_3D((*pml_E[0])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_E[1])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_E[2])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_B[0])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_B[1])[mfi]),
			     BL_TO_FORTRAN_3D((*pml_B[2])[mfi]),
			     &dtsdx_c2[0], &dtsdx_c2[1], &dtsdx_c2[2]);

            if (pml_F)
            {
                WRPX_PUSH_PML_EVEC_F(
				   tex.loVect(), tex.hiVect(),
				   tey.loVect(), tey.hiVect(),
				   tez.loVect(), tez.hiVect(),
				   BL_TO_FORTRAN_3D((*pml_E[0])[mfi]),
				   BL_TO_FORTRAN_3D((*pml_E[1])[mfi]),
				   BL_TO_FORTRAN_3D((*pml_E[2])[mfi]),
				   BL_TO_FORTRAN_3D((*pml_F   )[mfi]),
				   &dtsdx_c2[0], &dtsdx_c2[1], &dtsdx_c2[2],
				   &WarpX::maxwell_fdtd_solver_id);
            }
        }
    }
}

void
WarpX::EvolveF (Real dt, DtType dt_type)
{
    if (!do_dive_cleaning) return;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        EvolveF(lev, dt, dt_type);
    }
}

void
WarpX::EvolveF (int lev, Real dt, DtType dt_type)
{
    if (!do_dive_cleaning) return;

    EvolveF(lev, PatchType::fine, dt, dt_type);
    if (lev > 0) EvolveF(lev, PatchType::coarse, dt, dt_type);
}

void
WarpX::EvolveF (int lev, PatchType patch_type, Real dt, DtType dt_type)
{
    if (!do_dive_cleaning) return;

    BL_PROFILE("WarpX::EvolveF()");

    static constexpr Real mu_c2 = PhysConst::mu0*PhysConst::c*PhysConst::c;

    int patch_level = (patch_type == PatchType::fine) ? lev : lev-1;
    const auto& dx = WarpX::CellSize(patch_level);
    const std::array<Real,3> dtsdx {dt/dx[0], dt/dx[1], dt/dx[2]};

    MultiFab *Ex, *Ey, *Ez, *rho, *F;
    if (patch_type == PatchType::fine)
    {
        Ex = Efield_fp[lev][0].get();
        Ey = Efield_fp[lev][1].get();
        Ez = Efield_fp[lev][2].get();
        rho = rho_fp[lev].get();
        F = F_fp[lev].get();
    }
    else
    {
        Ex = Efield_cp[lev][0].get();
        Ey = Efield_cp[lev][1].get();
        Ez = Efield_cp[lev][2].get();
        rho = rho_cp[lev].get();
        F = F_cp[lev].get();
    }

    const int rhocomp = (dt_type == DtType::FirstHalf) ? 0 : 1;

    MultiFab src(rho->boxArray(), rho->DistributionMap(), 1, 0);
    ComputeDivE(src, 0, {Ex,Ey,Ez}, dx);
    MultiFab::Saxpy(src, -mu_c2, *rho, rhocomp, 0, 1, 0);
    MultiFab::Saxpy(*F, dt, src, 0, 0, 1, 0);

    if (do_pml && pml[lev]->ok())
    {
        const auto& pml_F = (patch_type == PatchType::fine) ? pml[lev]->GetF_fp() : pml[lev]->GetF_cp();
        const auto& pml_E = (patch_type == PatchType::fine) ? pml[lev]->GetE_fp() : pml[lev]->GetE_cp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*pml_F, TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const Box& bx = mfi.tilebox();
            WRPX_PUSH_PML_F(bx.loVect(), bx.hiVect(),
			  BL_TO_FORTRAN_ANYD((*pml_F   )[mfi]),
			  BL_TO_FORTRAN_ANYD((*pml_E[0])[mfi]),
			  BL_TO_FORTRAN_ANYD((*pml_E[1])[mfi]),
			  BL_TO_FORTRAN_ANYD((*pml_E[2])[mfi]),
			  &dtsdx[0], &dtsdx[1], &dtsdx[2]);
        }
    }
}

void
WarpX::DampPML ()
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        DampPML(lev);
    }
}

void
WarpX::DampPML (int lev)
{
    DampPML(lev, PatchType::fine);
    if (lev > 0) DampPML(lev, PatchType::coarse);
}

void
WarpX::DampPML (int lev, PatchType patch_type)
{
    if (!do_pml) return;

    BL_PROFILE("WarpX::DampPML()");

    if (pml[lev]->ok())
    {
        const auto& pml_E = (patch_type == PatchType::fine) ? pml[lev]->GetE_fp() : pml[lev]->GetE_cp();
        const auto& pml_B = (patch_type == PatchType::fine) ? pml[lev]->GetB_fp() : pml[lev]->GetB_cp();
        const auto& pml_F = (patch_type == PatchType::fine) ? pml[lev]->GetF_fp() : pml[lev]->GetF_cp();
        const auto& sigba = (patch_type == PatchType::fine) ? pml[lev]->GetMultiSigmaBox_fp()
                                                              : pml[lev]->GetMultiSigmaBox_cp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(*pml_E[0], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            const Box& tex  = mfi.tilebox(Ex_nodal_flag);
            const Box& tey  = mfi.tilebox(Ey_nodal_flag);
            const Box& tez  = mfi.tilebox(Ez_nodal_flag);
            const Box& tbx  = mfi.tilebox(Bx_nodal_flag);
            const Box& tby  = mfi.tilebox(By_nodal_flag);
            const Box& tbz  = mfi.tilebox(Bz_nodal_flag);

            WRPX_DAMP_PML(tex.loVect(), tex.hiVect(),
			    tey.loVect(), tey.hiVect(),
			    tez.loVect(), tez.hiVect(),
			    tbx.loVect(), tbx.hiVect(),
    			tby.loVect(), tby.hiVect(),
	    		tbz.loVect(), tbz.hiVect(),
		    	BL_TO_FORTRAN_3D((*pml_E[0])[mfi]),
			    BL_TO_FORTRAN_3D((*pml_E[1])[mfi]),
			    BL_TO_FORTRAN_3D((*pml_E[2])[mfi]),
			    BL_TO_FORTRAN_3D((*pml_B[0])[mfi]),
			    BL_TO_FORTRAN_3D((*pml_B[1])[mfi]),
			    BL_TO_FORTRAN_3D((*pml_B[2])[mfi]),
			    WRPX_PML_TO_FORTRAN(sigba[mfi]));

            if (pml_F) {
                const Box& tnd  = mfi.nodaltilebox();
                WRPX_DAMP_PML_F(tnd.loVect(), tnd.hiVect(),
			        BL_TO_FORTRAN_3D((*pml_F)[mfi]),
			        WRPX_PML_TO_FORTRAN(sigba[mfi]));
            }
        }
    }
}

void
WarpX::PushParticlesandDepose (Real cur_time)
{
    // Evolve particles to p^{n+1/2} and x^{n+1}
    // Depose current, j^{n+1/2}
    for (int lev = 0; lev <= finest_level; ++lev) {
        PushParticlesandDepose(lev, cur_time);
    }
}

void
WarpX::PushParticlesandDepose (int lev, Real cur_time)
{
    mypc->Evolve(lev,
                 *Efield_aux[lev][0],*Efield_aux[lev][1],*Efield_aux[lev][2],
                 *Bfield_aux[lev][0],*Bfield_aux[lev][1],*Bfield_aux[lev][2],
                 *current_fp[lev][0],*current_fp[lev][1],*current_fp[lev][2],
                 current_buf[lev][0].get(), current_buf[lev][1].get(), current_buf[lev][2].get(),
                 rho_fp[lev].get(), charge_buf[lev].get(),
                 Efield_cax[lev][0].get(), Efield_cax[lev][1].get(), Efield_cax[lev][2].get(),
                 Bfield_cax[lev][0].get(), Bfield_cax[lev][1].get(), Bfield_cax[lev][2].get(),
                 cur_time, dt[lev]);
}

void
WarpX::ComputeDt ()
{
    const Real* dx = geom[max_level].CellSize();
    Real deltat = 0.;

    if (maxwell_fdtd_solver_id == 0) {
      // CFL time step Yee solver
      deltat  = cfl * 1./( std::sqrt(AMREX_D_TERM(  1./(dx[0]*dx[0]),
                                                  + 1./(dx[1]*dx[1]),
                                                  + 1./(dx[2]*dx[2]))) * PhysConst::c );
    } else {
      // CFL time step CKC solver
#if (BL_SPACEDIM == 3)
      const Real delta = std::min(dx[0],std::min(dx[1],dx[2]));
#elif (BL_SPACEDIM == 2)
      const Real delta = std::min(dx[0],dx[1]);
#endif
      deltat = cfl*delta/PhysConst::c;
    }
    dt.resize(0);
    dt.resize(max_level+1,deltat);

    if (do_subcycling) {
        for (int lev = max_level-1; lev >= 0; --lev) {
            dt[lev] = dt[lev+1] * refRatio(lev)[0];
        }
    }

    if (do_electrostatic) {
        dt[0] = const_dt;
    }
}
