#include <Utility.H>
#include <ParallelDescriptor.H>
#include <MultiFab.H>
#include <PArray.H>
#include <ParmParse.H>

#include <algorithm>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

int
main (int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    BoxArray ba;

    std::string ba_file("ba.max");
    int max_grid_size = 32;
    int min_ba_size = 12800;
    {
	ParmParse pp;
	pp.query("ba_file", ba_file);
	pp.query("max_grid_size", max_grid_size);
	pp.query("min_ba_size", min_ba_size);
    }

    int nAtOnce = std::min(ParallelDescriptor::NProcs(), 32);

    // read boxarray    
    { 
	const int MyProc  = ParallelDescriptor::MyProc();
	const int NProcs  = ParallelDescriptor::NProcs();
	const int NSets   = (NProcs + (nAtOnce - 1)) / nAtOnce;
	const int MySet   = MyProc/nAtOnce;

	for (int iSet = 0; iSet < NSets; ++iSet)
	{
	    if (MySet == iSet)
	    {
                std::ifstream ifs(ba_file.c_str(), std::ios::in);
		ba.readFrom(ifs);

		const int iBuff     = 0;
		const int wakeUpPID = (MyProc + nAtOnce);
		const int tag       = (MyProc % nAtOnce);
		if (wakeUpPID < NProcs)
		    ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
	    }

	    if (MySet == (iSet + 1))
	    {
		//
		// Next set waits.
		//
		int iBuff;
		int waitForPID = (MyProc - nAtOnce);
		int tag        = (MyProc % nAtOnce);
		ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
	    }
	}
    }

    while (ba.size() < min_ba_size) {
	ba.refine(2);
	ba.maxSize(max_grid_size);
    };

    int nlevels;
    {
	int bsmin = 1024000;
	for (int i=0; i<ba.size(); ++i) {
	    for (int idim=0; idim<BL_SPACEDIM; ++idim) {
		bsmin = std::min(bsmin, ba[i].length(idim));
	    }
	}
	
	nlevels = 1;
	int L = bsmin;
	for (int lev=1; lev<10; ++lev) {
	    int Lold = L;
	    L /= 2;
	    if (L*2 != Lold) {
		break;
	    } else {
		nlevels++;
		if (L == 2) break;
	    }
	}

	if (ParallelDescriptor::IOProcessor()) {
	    std::cout << "min length = " << bsmin << std::endl;
	    std::cout << "num Pts    = " << ba.numPts() << std::endl;
	    std::cout << "num boxes  = " << ba.size() << std::endl;
	    std::cout << "num levels = " << nlevels << std::endl;
	}
    }

    ParallelDescriptor::Barrier();

    PArray<MultiFab> mfs(nlevels, PArrayManage);
    PArray<BoxArray> bas(nlevels, PArrayManage);
    bas.set(0, new BoxArray(ba));
    mfs.set(0, new MultiFab(ba, 1, 1));
    mfs[0].setVal(1.0);
    for (int lev=1; lev<nlevels; ++lev) {
	bas.set(lev, new BoxArray(bas[lev-1]));
	bas[lev].coarsen(2);
	mfs.set(lev, new MultiFab(bas[lev], 1, 1));
	mfs[lev].setVal(1.0);
    }

    Array<Real> points(nlevels);
    for (int lev=0; lev<nlevels; ++lev) {
	points[lev] = mfs[lev].norm1();
	if (ParallelDescriptor::IOProcessor()) {
	    std::cout << points[lev] << " points on level " << lev << std::endl; 
	}
    }

    int nrounds = 1000;
    {
	ParmParse pp;
	pp.query("nrounds", nrounds);
    }

    Real err = 0.0;

    ParallelDescriptor::Barrier();
    Real wt0 = ParallelDescriptor::second();

    for (int iround = 0; iround < nrounds; ++iround) {
	for (int c=0; c<2; ++c) {
	    for (int lev = 0; lev < nlevels; ++lev) {
		mfs[lev].FillBoundary_nowait(true);
		mfs[lev].FillBoundary_finish();
	    }
	    for (int lev = nlevels-1; lev >= 0; --lev) {
		mfs[lev].FillBoundary_nowait(true);
		mfs[lev].FillBoundary_finish();
	    }
	}
	Real e = double(iround+ParallelDescriptor::MyProc());
	ParallelDescriptor::ReduceRealMax(e);
	err += e;
    }
	    
    ParallelDescriptor::Barrier();
    Real wt1 = ParallelDescriptor::second();

    if (ParallelDescriptor::IOProcessor()) {
#ifdef BL_USE_UPCXX
        std::cout << "Using UPCXX" << std::endl;
#else
	if (ParallelDescriptor::MPIOneSided()) {
	    std::cout << "Using MPI Onesided" << std::endl;
	} else {
	    std::cout << "Using MPI" << std::endl;
	}
#endif
	std::cout << "----------------------------------------------" << std::endl;
	std::cout << "Fill Boundary Time: " << wt1-wt0 << std::endl;
	std::cout << "----------------------------------------------" << std::endl;
	std::cout << "ignore this line " << err << std::endl;
    }

    //
    // When MPI3 shared memory is used, the dtor of MultiFab calls MPI functions.
    // Because the scope of PArray < MultiFab > mfs is beyond the call to 
    // BoxLib::Finalize(), which in turn calls MPI_Finalize(), we destroy these
    // MultiFabs by hand now.
    //
    mfs.clear();

    BoxLib::Finalize();
}
