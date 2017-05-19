
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

        if (is_synchronized) {
            // on first step, push E and X by 0.5*dt
            FillBoundaryB();
            EvolveE(0.5*dt[0]);
            mypc->PushX(0.5*dt[0]);
            mypc->Redistribute();  // Redistribute particles
            is_synchronized = false;
        }

        EvolveB(0.5*dt[0]); // We now B^{n}

        if (WarpX::nox > 1 || WarpX::noz > 1 || WarpX::noz > 1) {
            FillBoundaryB();
            FillBoundaryE();
        }
        
        for (int lev = 0; lev <= finest_level; ++lev) {
	    // Evolve particles to p^{n+1/2} and x^{n+1}
	    // Depose current, j^{n+1/2}

            MultiFab* pjx_bnd = (lev == finest_level) ? nullptr : &(cfbndry[lev]->GetJx());
            MultiFab* pjy_bnd = (lev == finest_level) ? nullptr : &(cfbndry[lev]->GetJy());
            MultiFab* pjz_bnd = (lev == finest_level) ? nullptr : &(cfbndry[lev]->GetJz());

	    mypc->Evolve(lev,
			 *Efield[lev][0],*Efield[lev][1],*Efield[lev][2],
			 *Bfield[lev][0],*Bfield[lev][1],*Bfield[lev][2],
			 *current[lev][0],*current[lev][1],*current[lev][2],
                         pjx_bnd, pjy_bnd, pjz_bnd,
                         cur_time, dt[lev]);
        }

        EvolveB(0.5*dt[0]); // We now B^{n+1/2}

        SyncCurrent();

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

	if (to_make_plot) {
            if (WarpX::nox > 1 || WarpX::noz > 1 || WarpX::noz > 1) {
                FillBoundaryB();
                FillBoundaryE();
            }
            for (int lev = 0; lev <= finest_level; ++lev) {
                mypc->FieldGather(lev,
                                  *Efield[lev][0],*Efield[lev][1],*Efield[lev][2],
                                  *Bfield[lev][0],*Bfield[lev][1],*Bfield[lev][2]);
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
    AverageDownB();
}

void
WarpX::EvolveB (int lev, Real dt)
{
    BL_PROFILE("WarpX::EvolveB()");

    // Parameters of the solver: order and mesh spacing
    const int norder = 2;

    const std::array<Real,3>& dx = WarpX::CellSize(lev);
    const std::array<Real,3> dtsdx {dt/dx[0], dt/dx[1], dt/dx[2]};

    // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(*Bfield[lev][0],true); mfi.isValid(); ++mfi )
    {
        const Box& tbx  = mfi.tilebox(Bx_nodal_flag);
        const Box& tby  = mfi.tilebox(By_nodal_flag);
        const Box& tbz  = mfi.tilebox(Bz_nodal_flag);

	// Call picsar routine for each tile
	WRPX_PXR_PUSH_BVEC(
	    tbx.loVect(), tbx.hiVect(),
	    tby.loVect(), tby.hiVect(),
	    tbz.loVect(), tbz.hiVect(),
	    BL_TO_FORTRAN_3D((*Efield[lev][0])[mfi]),
	    BL_TO_FORTRAN_3D((*Efield[lev][1])[mfi]),
	    BL_TO_FORTRAN_3D((*Efield[lev][2])[mfi]),
	    BL_TO_FORTRAN_3D((*Bfield[lev][0])[mfi]),
	    BL_TO_FORTRAN_3D((*Bfield[lev][1])[mfi]),
	    BL_TO_FORTRAN_3D((*Bfield[lev][2])[mfi]),
	    &dtsdx[0], &dtsdx[1], &dtsdx[2],
	    &norder);
    }

    if (do_pml && lev == 0)
    {

        ComputePMLFactorsB(lev, dt);

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
                pml_sigma_star_fac1[0].data(),pml_sigma_star_fac1[0].lo(),pml_sigma_star_fac1[0].hi(),
                pml_sigma_star_fac2[0].data(),pml_sigma_star_fac2[0].lo(),pml_sigma_star_fac2[0].hi(),
                pml_sigma_star_fac1[1].data(),pml_sigma_star_fac1[1].lo(),pml_sigma_star_fac1[1].hi(),
                pml_sigma_star_fac2[1].data(),pml_sigma_star_fac2[1].lo(),pml_sigma_star_fac2[1].hi()
#if (BL_SPACEDIM == 3)
               ,pml_sigma_star_fac1[2].data(),pml_sigma_star_fac1[2].lo(),pml_sigma_star_fac1[2].hi(),
                pml_sigma_star_fac2[2].data(),pml_sigma_star_fac2[2].lo(),pml_sigma_star_fac2[2].hi()
#endif
                );
        }
    }
}

void
WarpX::EvolveE (Real dt)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        EvolveE(lev, dt);
    }
    AverageDownE();
}

void
WarpX::EvolveE (int lev, Real dt)
{
    BL_PROFILE("WarpX::EvolveE()");

    // Parameters of the solver: order and mesh spacing
    const int norder = 2;

    const std::array<Real,3>& dx = WarpX::CellSize(lev);
    const Real mu_c2_dt = (PhysConst::mu0*PhysConst::c*PhysConst::c) * dt;
    const Real foo = (PhysConst::c*PhysConst::c) * dt;
    const std::array<Real,3> dtsdx_c2 {foo/dx[0], foo/dx[1], foo/dx[2]};

  // Loop through the grids, and over the tiles within each grid
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(*Efield[lev][0],true); mfi.isValid(); ++mfi )
    {
	const Box& tex  = mfi.tilebox(Ex_nodal_flag);
	const Box& tey  = mfi.tilebox(Ey_nodal_flag);
	const Box& tez  = mfi.tilebox(Ez_nodal_flag);

        // Call picsar routine for each tile
	WRPX_PXR_PUSH_EVEC(
	    tex.loVect(), tex.hiVect(),
	    tey.loVect(), tey.hiVect(),
	    tez.loVect(), tez.hiVect(),
	    BL_TO_FORTRAN_3D((*Efield[lev][0])[mfi]),
	    BL_TO_FORTRAN_3D((*Efield[lev][1])[mfi]),
	    BL_TO_FORTRAN_3D((*Efield[lev][2])[mfi]),
	    BL_TO_FORTRAN_3D((*Bfield[lev][0])[mfi]),
	    BL_TO_FORTRAN_3D((*Bfield[lev][1])[mfi]),
	    BL_TO_FORTRAN_3D((*Bfield[lev][2])[mfi]),
	    BL_TO_FORTRAN_3D((*current[lev][0])[mfi]),
	    BL_TO_FORTRAN_3D((*current[lev][1])[mfi]),
	    BL_TO_FORTRAN_3D((*current[lev][2])[mfi]),
	    &mu_c2_dt,
	    &dtsdx_c2[0], &dtsdx_c2[1], &dtsdx_c2[2],
	    &norder);
    }

    if (do_pml && lev == 0)
    {

        ComputePMLFactorsE(lev, dt);

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
                pml_sigma_fac1[0].data(),pml_sigma_fac1[0].lo(),pml_sigma_fac1[0].hi(),
                pml_sigma_fac2[0].data(),pml_sigma_fac2[0].lo(),pml_sigma_fac2[0].hi(),
                pml_sigma_fac1[1].data(),pml_sigma_fac1[1].lo(),pml_sigma_fac1[1].hi(),
                pml_sigma_fac2[1].data(),pml_sigma_fac2[1].lo(),pml_sigma_fac2[1].hi()
#if (BL_SPACEDIM == 3)
               ,pml_sigma_fac1[2].data(),pml_sigma_fac1[2].lo(),pml_sigma_fac1[2].hi(),
                pml_sigma_fac2[2].data(),pml_sigma_fac2[2].lo(),pml_sigma_fac2[2].hi()
#endif
                );
        }
    }
}

void
WarpX::PushParticlesandDepose(int lev, Real cur_time)
{
    // Evolve particles to p^{n+1/2} and x^{n+1}
    // Depose current, j^{n+1/2}
    mypc->Evolve(lev,
		 *Efield[lev][0],*Efield[lev][1],*Efield[lev][2],
		 *Bfield[lev][0],*Bfield[lev][1],*Bfield[lev][2],
		 *current[lev][0],*current[lev][1],*current[lev][2],
                 nullptr, nullptr, nullptr, // xxxxx is this the right thing to do?
                 cur_time, dt[lev]);
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

