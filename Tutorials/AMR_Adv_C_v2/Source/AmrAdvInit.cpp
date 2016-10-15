
#include <AmrAdv.H>
#include <AmrAdv_F.H>

void
AmrAdv::InitData ()
{
    if (restart_chkfile.empty()) {
	InitFromScratch();
    } else {
	InitFromCheckpoint();
    }
}

void
AmrAdv::InitFromScratch ()
{
    // define coarse level BoxArray and DistributionMap
    {
	SetBoxArray(0, MakeBaseGrids());
	SetDistributionMap(0, DistributionMapping(grids[0], ParallelDescriptor::NProcs()));

	int ncomp = 1;
	int nghost = 0;
	phi_new[0] = std::unique_ptr<MultiFab>(new MultiFab(grids[0], ncomp, nghost, dmap[0]));
	phi_old[0] = std::unique_ptr<MultiFab>(new MultiFab(grids[0], ncomp, nghost, dmap[0]));

	InitLevelData(0);
    }
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
