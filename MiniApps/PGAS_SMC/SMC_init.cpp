
#include <SMC.H>
#include <SMC_F.H>

void
SMC::build_multifabs ()
{
    IntVect dlo(0, 0, 0);
    IntVect dhi(ncell[0]-1, ncell[1]-1, ncell[2]-1);
    Box bx(dlo, dhi);

    BoxArray ba(bx);
    IntVect mgs{max_grid_size};
    ba.maxSize(mgs);

    int dir = 2;
    while (ba.size() < ParallelDescriptor::NTeams()) {
	if (dir == 2) {
	    ba.maxSize(IntVect{mgs[0],mgs[1],mgs[2]/2});
	    dir = 1;
	} else if (dir == 1) {
	    ba.maxSize(IntVect{mgs[0],mgs[1]/2,mgs[2]/2});
	    dir = 0;
	} else {
	    ba.maxSize(IntVect{mgs[0]/2,mgs[1]/2,mgs[2]/2});
	    dir = 2;
	    mgs = IntVect{mgs[0]/2,mgs[1]/2,mgs[2]/2};
	}
    }

    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Number of boxes: " << ba.size() << std::endl;
    }

    RealBox real_box(prob_lo.dataPtr(), prob_hi.dataPtr());  // physical size
    int coord = 0; // Cartesian coordinates
    Array<int> is_per(3, 1); // triply periodic
    geom.define(bx, &real_box, 0, &is_per[0]);

         U.define(ba, ncons, ngrow, Fab_allocate);
      Utmp.define(ba, ncons, ngrow, Fab_allocate);
    Uprime.define(ba, ncons, 0    , Fab_allocate);
         Q.define(ba, nprim, ngrow, Fab_allocate);
        mu.define(ba, 1    , ngrow, Fab_allocate);
        xi.define(ba, 1    , ngrow, Fab_allocate);
       lam.define(ba, 1    , ngrow, Fab_allocate);
     Ddiag.define(ba, nspec, ngrow, Fab_allocate);

    Q.setVal(0.0);
}

void
SMC::init_from_scratch ()
{
    t = 0.0;
    dt = 1.e10;
    courno = -1.e10;
    for (int i=0; i<3; ++i) {
	dx[i] = (prob_hi[i]-prob_lo[i]) / Real(ncell[i]);
    }

#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (MFIter mfi(U,true); mfi.isValid(); ++mfi)
    {
	const Box& tbx  = mfi.tilebox();
	const Box&  bx = mfi.validbox();

	init_data_3d(tbx.loVect(), tbx.hiVect(), bx.loVect(), bx.hiVect(), 
		     U[mfi].dataPtr(), ngrow, dx, prob_lo.dataPtr(), prob_hi.dataPtr());
    }
}
