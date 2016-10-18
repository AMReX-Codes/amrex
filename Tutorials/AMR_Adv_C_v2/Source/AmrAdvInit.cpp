
#include <AmrAdv.H>
#include <AmrAdv_F.H>
#include <PlotFileUtil.H>

std::string
AmrAdv::PlotFileName (int lev) const
{
    return BoxLib::Concatenate(plot_file, lev, 5);
}

Array<const MultiFab*>
AmrAdv::PlotFileMF () const
{
    Array<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
	r.push_back(phi_new[i].get());
    }
    return r;
}

std::vector<std::string>
AmrAdv::PlotFileVarNames () const
{
    return {"phi"};
}

void
AmrAdv::InitData ()
{
    if (restart_chkfile.empty())
    {
	InitFromScratch();

	if (plot_int > 0) {
	    const std::string& plotfilename = PlotFileName(0);
	    const auto& mf = PlotFileMF();
	    const auto& varnames = PlotFileVarNames();
	    
	    BoxLib::WriteMultiLevelPlotifle(plotfilename, finest_level+1, mf, varnames,
					    Geom(), t_new[0], istep, refRatio());
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
    const int ncomp = 1;
    const int nghost = 0;
    const Real time = 0.0;
    const Real old_time = -1.e200;

    // define coarse level BoxArray and DistributionMap
    {
	finest_level = 0;

	SetBoxArray(0, MakeBaseGrids());
	SetDistributionMap(0, DistributionMapping(grids[0], ParallelDescriptor::NProcs()));

	phi_new[0] = std::unique_ptr<MultiFab>(new MultiFab(grids[0], ncomp, nghost, dmap[0]));
	phi_old[0] = std::unique_ptr<MultiFab>(new MultiFab(grids[0], ncomp, nghost, dmap[0]));
	t_new[0] = time;
	t_old[0] = old_time;

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

	    phi_new[new_finest] = std::unique_ptr<MultiFab>
		(new MultiFab(grids[new_finest], ncomp, nghost, dmap[new_finest]));
	    phi_old[new_finest] = std::unique_ptr<MultiFab>
		(new MultiFab(grids[new_finest], ncomp, nghost, dmap[new_finest]));
	    t_new[new_finest] = time;
	    t_old[new_finest] = old_time;

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
