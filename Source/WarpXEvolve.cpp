
#include <cmath>
#include <limits>

#include <WarpX.H>
#include <WarpXConst.H>
#include <WarpX_f.H>
#include <WarpX_py.H>

using namespace amrex;

void
WarpX::Evolve (int numsteps)
{
    BL_PROFILE("WarpX::Evolve()");

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

    for (int step = istep[0]; step < numsteps_max && cur_time < stop_time; ++step)
    {
        if (warpx_py_print_step) {
            warpx_py_print_step(step);
        }

	// Start loop on time steps
        amrex::Print() << "\nSTEP " << step+1 << " starts ...\n";

        // At the beginning, we have B^{n-1/2} and E^{n-1/2}.
        // Particles have p^{n-1/2} and x^{n-1/2}.

        // Beyond one step, we have B^{n-1/2} and E^{n}.
        // Particles have p^{n-1/2} and x^{n}.
        // F for div E cleaning is at n-1/2.

        if (is_synchronized) {
            // on first step, push E and X by 0.5*dt
            FillBoundaryB();
            EvolveE(0.5*dt[0]);
            mypc->PushX(0.5*dt[0]);
            mypc->Redistribute();  // Redistribute particles
            is_synchronized = false;
        }

        FillBoundaryE();

        EvolveB(0.5*dt[0]); // We now B^{n}

        FillBoundaryB();

        UpdateAuxilaryData();

        // Push particle from x^{n} to x{n+1}
        // Deposit current j^{n+1/2}
        // Deposit charge density rho^{n}
        PushParticlesandDepose(cur_time);

        EvolveB(0.5*dt[0]); // We now B^{n+1/2}

        SyncCurrent();

        // xxxxx SyncRho and then evolve F if ...
        // We now have F^{n+1/2}

        // Fill B's ghost cells because of the next step of evolving E.
        FillBoundaryB();

        if (cur_time + dt[0] >= stop_time - 1.e-3*dt[0] || step == numsteps_max-1) {
            // on last step, push by only 0.5*dt to synchronize all at n+1/2
            EvolveE(0.5*dt[0]); // We now have E^{n+1/2}
            mypc->PushX(-0.5*dt[0]);
            is_synchronized = true;
        } else {
            EvolveE(dt[0]); // We now have E^{n+1}
        }
        
        mypc->Redistribute();  // Redistribute particles

         for (int lev = 0; lev <= max_level; ++lev) {
            ++istep[lev];
        }

	cur_time += dt[0];

        bool to_make_plot = (plot_int > 0) && ((step+1) % plot_int == 0);

        bool move_j = is_synchronized || to_make_plot;
        // If is_synchronized we need to shift j too so that next step we can evolve E by dt/2.
        // We might need to move j because we are going to make a plotfile.
	MoveWindow(move_j);

        amrex::Print()<< "STEP " << step+1 << " ends." << " TIME = " << cur_time
                      << " DT = " << dt[0] << "\n";

	// sync up time
	for (int i = 0; i <= max_level; ++i) {
	    t_new[i] = cur_time;
	}

	if (to_make_plot)
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
	    WritePlotFile();
	}

	if (check_int > 0 && (step+1) % check_int == 0) {
	    last_check_file_step = step+1;
	    WriteCheckPointFile();
	}

	if (cur_time >= stop_time - 1.e-3*dt[0]) {
	    max_time_reached = true;
	    break;
	}

	// End loop on time steps
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step && (max_time_reached || istep[0] >= max_step)) {
	WritePlotFile();
    }

    if (check_int > 0 && istep[0] > last_check_file_step && (max_time_reached || istep[0] >= max_step)) {
	WriteCheckPointFile();
    }
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

    // Parameters of the solver: order and mesh spacing
    const int norder = 2;

    int npatches = (lev == 0) ? 1 : 2;

    for (int ipatch = 0; ipatch < npatches; ++ipatch)
    {
        int patch_level = (ipatch == 0) ? lev : lev-1;
        const std::array<Real,3>& dx = WarpX::CellSize(patch_level);
        const std::array<Real,3> dtsdx {dt/dx[0], dt/dx[1], dt/dx[2]};
        
        MultiFab *Ex, *Ey, *Ez, *Bx, *By, *Bz;
        if (ipatch == 0)
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

        // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(*Bx,true); mfi.isValid(); ++mfi )
        {
            const Box& tbx  = mfi.tilebox(Bx_nodal_flag);
            const Box& tby  = mfi.tilebox(By_nodal_flag);
            const Box& tbz  = mfi.tilebox(Bz_nodal_flag);
            
            // Call picsar routine for each tile
            WRPX_PXR_PUSH_BVEC(
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
                &norder);
        }
    }

    if (do_pml && pml[lev]->ok())
    {
        pml[lev]->ComputePMLFactorsB(dt);

        for (int ipatch = 0; ipatch < npatches; ++ipatch)
        {
            const auto& pml_B = (ipatch==0) ? pml[lev]->GetB_fp() : pml[lev]->GetB_cp();
            const auto& pml_E = (ipatch==0) ? pml[lev]->GetE_fp() : pml[lev]->GetE_cp();
            const auto& sigba = (ipatch==0) ? pml[lev]->GetMultiSigmaBox_fp()
                                            : pml[lev]->GetMultiSigmaBox_cp();
            
#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(*pml_B[0],true); mfi.isValid(); ++mfi )
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
                    WRPX_PML_SIGMA_STAR_TO_FORTRAN(sigba[mfi]));
            }
        }
    }
}

void
WarpX::EvolveE (Real dt)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        EvolveE(lev, dt);
    }
}

void
WarpX::EvolveE (int lev, Real dt)
{
    BL_PROFILE("WarpX::EvolveE()");

    // Parameters of the solver: order and mesh spacing
    const int norder = 2;
    const Real mu_c2_dt = (PhysConst::mu0*PhysConst::c*PhysConst::c) * dt;
    const Real foo = (PhysConst::c*PhysConst::c) * dt;

    int npatches = (lev == 0) ? 1 : 2;

    for (int ipatch = 0; ipatch < npatches; ++ipatch)
    {
        int patch_level = (ipatch == 0) ? lev : lev-1;
        const std::array<Real,3>& dx = WarpX::CellSize(patch_level);
        const std::array<Real,3> dtsdx_c2 {foo/dx[0], foo/dx[1], foo/dx[2]};

        MultiFab *Ex, *Ey, *Ez, *Bx, *By, *Bz, *jx, *jy, *jz;
        if (ipatch == 0)
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
        }
        else
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
        }
        
        // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel
#endif
        for ( MFIter mfi(*Ex,true); mfi.isValid(); ++mfi )
        {
            const Box& tex  = mfi.tilebox(Ex_nodal_flag);
            const Box& tey  = mfi.tilebox(Ey_nodal_flag);
            const Box& tez  = mfi.tilebox(Ez_nodal_flag);

            // Call picsar routine for each tile
            WRPX_PXR_PUSH_EVEC(
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
                &dtsdx_c2[0], &dtsdx_c2[1], &dtsdx_c2[2],
                &norder);
        }
    }

    if (do_pml && pml[lev]->ok())
    {
        pml[lev]->ComputePMLFactorsE(dt);

        for (int ipatch = 0; ipatch < npatches; ++ipatch)
        {
            const auto& pml_B = (ipatch==0) ? pml[lev]->GetB_fp() : pml[lev]->GetB_cp();
            const auto& pml_E = (ipatch==0) ? pml[lev]->GetE_fp() : pml[lev]->GetE_cp();
            const auto& sigba = (ipatch==0) ? pml[lev]->GetMultiSigmaBox_fp()
                                            : pml[lev]->GetMultiSigmaBox_cp();

#ifdef _OPENMP
#pragma omp parallel
#endif
            for ( MFIter mfi(*pml_E[0],true); mfi.isValid(); ++mfi )
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
                    WRPX_PML_SIGMA_TO_FORTRAN(sigba[mfi]));
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
                 rho_fp[lev].get(), cur_time, dt[lev]);
}

void
WarpX::ComputeDt ()
{
    const Real* dx = geom[max_level].CellSize();
    const Real deltat  = cfl * 1./( std::sqrt(D_TERM(  1./(dx[0]*dx[0]),
                                                     + 1./(dx[1]*dx[1]),
                                                     + 1./(dx[2]*dx[2]))) * PhysConst::c );
    dt.resize(0);
    dt.resize(max_level+1,deltat);
}

void
WarpX::InjectPlasma (int num_shift, int dir)
{
    if(do_plasma_injection)
    {
        const int lev = 0;

        // particleBox encloses the cells where we generate particles
        Box particleBox = geom[lev].Domain();
        int domainLength = particleBox.length(dir);
        int sign = (num_shift < 0) ? -1 : 1;
        particleBox.shift(dir, sign*(domainLength - std::abs(num_shift)));
        particleBox &= geom[lev].Domain();

        for (int i = 0; i < num_injected_species; ++i) {
            int ispecies = injected_plasma_species[i];
            WarpXParticleContainer& pc = mypc->GetParticleContainer(ispecies);
            auto& ppc = dynamic_cast<PhysicalParticleContainer&>(pc);
            ppc.AddParticles(lev, particleBox);
        }
    }
}

