
#include <WarpX.H>

void
WarpX::InitData ()
{
    BL_PROFILE("WarpX::InitData()");

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
WarpX::InitFromScratch ()
{
    BL_ASSERT(max_level == 0);

    const Real time = 0.0;
    
    // define coarse level BoxArray and DistributionMap
    {
	finest_level = 0;

	const BoxArray& ba = MakeBaseGrids();
	DistributionMapping dm(ba, ParallelDescriptor::NProcs());

	MakeNewLevel(0, time, ba, dm);

	InitLevelData(0);
    }

    // if max_level > 0, define fine levels

    mypc->AllocData();
    mypc->InitData();
}

void
WarpX::InitLevelData (int lev)
{
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	current[lev][i]->setVal(0.0);
	Efield[lev][i]->setVal(0.0);
	Bfield[lev][i]->setVal(0.0);
    }
}
