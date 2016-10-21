
#include <AmrAdv.H>
#include <AmrAdv_F.H>

void
AmrAdv::InitData ()
{
    if (restart_chkfile.empty())
    {
	InitFromScratch();

	if (plot_int > 0) {
	    WritePlotFile();
	}
    }
    else
    {
	InitFromCheckpoint();
    }
}

void
AmrAdv::InitFromScratch ()
{
    const Real time = 0.0;

    // define coarse level BoxArray and DistributionMap
    {
	finest_level = 0;

	const BoxArray& ba = MakeBaseGrids();
	DistributionMapping dm(ba, ParallelDescriptor::NProcs());

	MakeNewLevel(0, time, ba, dm);

	InitLevelData(0);
    }

    if (max_level > 0) // build fine levels
    {
	Array<BoxArray> new_grids(max_level+1);
	new_grids[0] = grids[0];
	do
	{
	    int new_finest;

	    // Add (at most) one level at a time.
	    MakeNewGrids(finest_level,time,new_finest,new_grids);

	    if (new_finest <= finest_level) break;
	    finest_level = new_finest;

	    DistributionMapping dm(new_grids[new_finest], ParallelDescriptor::NProcs());

	    MakeNewLevel(new_finest, time, new_grids[new_finest], dm);

	    InitLevelData(new_finest);
	}
	while (finest_level < max_level);

	// Iterate grids to ensure fine grids encompass all interesting junk.
	for (int it=0; it<4; ++it)  // try at most 4 times
	{
	    for (int i = 1; i <= finest_level; ++i) {
		new_grids[i] = grids[i];
	    }

	    int new_finest;
	    MakeNewGrids(0, time, new_finest, new_grids);
	    
	    if (new_finest < finest_level) break;
	    finest_level = new_finest;

	    bool grids_the_same = true;
	    for (int lev = 1; lev <= new_finest; ++lev) {
		if (new_grids[lev] != grids[lev]) {
		    grids_the_same = false;
		    DistributionMapping dm(new_grids[lev], ParallelDescriptor::NProcs());
		    MakeNewLevel(lev, time, new_grids[lev], dm);
		    InitLevelData(new_finest);
		}
	    }
	    if (grids_the_same) break;
	}
    }

    AverageDown();
}

void AmrAdv::InitLevelData (int lev)
{
    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real cur_time = t_new[lev];

    MultiFab& state = *phi_new[lev];

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

	initdata(lev, cur_time, ARLIM_3D(lo), ARLIM_3D(hi),
		 BL_TO_FORTRAN_3D(state[mfi]), ZFILL(dx),
		 ZFILL(prob_lo));
    }
}
