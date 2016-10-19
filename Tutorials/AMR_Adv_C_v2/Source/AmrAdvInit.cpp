
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

	SetBoxArray(0, MakeBaseGrids());
	SetDistributionMap(0, DistributionMapping(grids[0], ParallelDescriptor::NProcs()));

	MakeNewLevel(0, time);

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

	    SetBoxArray(new_finest, new_grids[new_finest]);
	    SetDistributionMap(new_finest, DistributionMapping(new_grids[new_finest],
							       ParallelDescriptor::NProcs()));

	    MakeNewLevel(new_finest, time);

	    InitLevelData(new_finest);
	}
	while (finest_level < max_level);
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
