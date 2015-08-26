#include <SMC.H>
#include <SMC_F.H>

void
SMC::advance (int istep)
{
    compute_dUdt(U, istep);
    Utmp.setVal(0.0);
    set_dt(istep);
    update_rk3(0.0, Utmp, 1.0, U, dt, Uprime);
    reset_density(Utmp);

    compute_dUdt(Utmp);
    update_rk3(0.25, Utmp, 0.75, U, 0.25*dt, Uprime);
    reset_density(Utmp);

    compute_dUdt(Utmp);
    update_rk3(1./3., U, 2./3., Utmp, (2./3.)*dt, Uprime);
    reset_density(U);
}

void 
SMC::compute_dUdt (MultiFab& UU, int istep)
{
    bool update_courno = false;
    if (istep >=0 && fixed_dt <= 0.0) {
	if (istep%cfl_int == 1 || cfl_int <= 1) update_courno = true;
    }

    UU.FillBoundary_nowait();
    geom.FillPeriodicBoundary_nowait(UU);

    if (!overlap_comm_comp) {
	UU.FillBoundary_finish();
	geom.FillPeriodicBoundary_finish(UU);	
    }

    Uprime.setVal(0.0);

    int ng_ctoprim = (overlap_comm_comp) ? 0 : ngrow;

#ifdef _OPENMP
#pragma omp parallel reduction(max:courno)
#endif
    for (MFIter mfi(UU); mfi.isValid(); ++mfi) {
	const Box& tbx0 = mfi.tilebox();
	const Box& tbxg = mfi.growntilebox(ng_ctoprim);
	const Box&  bx  = mfi.validbox();

	ctoprim_3d(tbxg.loVect(), tbxg.hiVect(), bx.loVect(), bx.hiVect(),
		   UU[mfi].dataPtr(), Q[mfi].dataPtr(), ngrow, ngrow);

	if (update_courno) {
	    comp_courno_3d(tbx0.loVect(), tbx0.hiVect(), dx,
			   BL_TO_FORTRAN_3D(Q[mfi]), courno);
	}

	chemterm_3d(tbx0.loVect(), tbx0.hiVect(),
		    BL_TO_FORTRAN_3D(     Q[mfi]),
		    BL_TO_FORTRAN_3D(Uprime[mfi]));

	get_trans_prop_3d(tbxg.loVect(), tbxg.hiVect(), bx.loVect(), bx.hiVect(),
			  Q[mfi].dataPtr(), mu[mfi].dataPtr(), xi[mfi].dataPtr(), 
			  lam[mfi].dataPtr(), Ddiag[mfi].dataPtr(), ngrow);
    }

    if (overlap_comm_comp) {
	UU.FillBoundary_finish();
	geom.FillPeriodicBoundary_finish(UU);

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFGhostIter mfi(UU); mfi.isValid(); ++mfi) {
	    const Box& tbx = mfi.tilebox();
	    const Box&  bx = mfi.validbox();

	    ctoprim_3d(tbx.loVect(), tbx.hiVect(), bx.loVect(), bx.hiVect(),
		       UU[mfi].dataPtr(), Q[mfi].dataPtr(), ngrow, ngrow);
	    
	    get_trans_prop_3d(tbx.loVect(), tbx.hiVect(), bx.loVect(), bx.hiVect(),
			      Q[mfi].dataPtr(), mu[mfi].dataPtr(), xi[mfi].dataPtr(), 
			      lam[mfi].dataPtr(), Ddiag[mfi].dataPtr(), ngrow);
	}
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(UU); mfi.isValid(); ++mfi) {
	const Box& tbx = mfi.tilebox();
	const Box&  bx = mfi.validbox();

	hypterm_3d(tbx.loVect(), tbx.hiVect(), dx,
		   BL_TO_FORTRAN_3D(    UU[mfi]),
		   BL_TO_FORTRAN_3D(     Q[mfi]),
		   BL_TO_FORTRAN_3D(Uprime[mfi]));

	narrow_diffterm_3d(tbx.loVect(), tbx.hiVect(), dx,
			   BL_TO_FORTRAN_3D(     Q[mfi]),
			   BL_TO_FORTRAN_3D(Uprime[mfi]),
			   mu[mfi].dataPtr(), xi[mfi].dataPtr(), 
			   lam[mfi].dataPtr(), Ddiag[mfi].dataPtr());
    }


    if (update_courno) {
	ParallelDescriptor::ReduceRealMax(courno);
    }
}

void
SMC::set_dt(int istep)
{
    const Real max_dt_growth = 1.1;

    if (fixed_dt > 0.0) {
	dt = fixed_dt;
	if (ParallelDescriptor::IOProcessor()) {
	    std::cout << "\n" << "Setting fixed dt = " << dt << "\n" << std::endl;
	}
    } else {
	Real dtold = dt;
	dt = cfl/courno;

	if (ParallelDescriptor::IOProcessor()) {
	    std::cout << "CFL: dt = " << dt << "\n";
	}

	if (istep == 1) {
	    dt *= init_shrink;
	    if (ParallelDescriptor::IOProcessor()) {
		std::cout << "Limited by init_shrink: dt = " << dt << "\n";
	    }	    
	} else {
	    if (dt > dtold*max_dt_growth) {
		dt = dtold*max_dt_growth;
		if (ParallelDescriptor::IOProcessor()) {
		    std::cout << "Limited by dt_growth: dt = " << dt << "\n";
		}	    
	    }
	}

	if (stop_time > 0.0) {
	    if (t+dt > stop_time) {
		dt = stop_time - t;
		if (ParallelDescriptor::IOProcessor()) {
		    std::cout << "Limited by stop_time: dt = " << dt << "\n";
		}
	    }
	}

	if (ParallelDescriptor::IOProcessor()) std::cout << std::endl;
    }
}

// Ua = a*Ua + b*Ub + c*Up
void
SMC::update_rk3 (Real a, MultiFab& Ua, Real b, const MultiFab& Ub, Real c, const MultiFab& Uc)
{
    int ncomp = Ua.nComp();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Ua,true); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.tilebox();
	Ua[mfi].mult(a, bx, 0, ncomp);
	Ua[mfi].saxpy(b, Ub[mfi], bx, bx, 0, 0, ncomp);
	Ua[mfi].saxpy(c, Uc[mfi], bx, bx, 0, 0, ncomp);
    }
}

void
SMC::reset_density (MultiFab& UU)
{
    int ngrow = UU.nGrow();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(UU,true); mfi.isValid(); ++mfi) {
	const Box& tbx = mfi.tilebox();
	const Box&  bx = mfi.validbox();
	reset_rho_3d(tbx.loVect(), tbx.hiVect(), bx.loVect(), bx.hiVect(), 
		     UU[mfi].dataPtr(), ngrow);
    }
}

