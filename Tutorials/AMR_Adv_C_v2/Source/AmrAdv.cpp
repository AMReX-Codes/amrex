
#include <ParallelDescriptor.H>
#include <ParmParse.H>
#include <MultiFabUtil.H>

#include <AmrAdv.H>

AmrAdv::AmrAdv ()
{
    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = maxLevel() + 1;

    istep.resize(nlevs_max, 0);
    isubstep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= maxLevel(); ++lev) {
	nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 0.0);
    dt_min.resize(nlevs_max, 0.0);

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);
}

AmrAdv::~AmrAdv ()
{

}

void
AmrAdv::ReadParameters ()
{
    ParmParse pp("adv");

    pp.query("max_step", max_step);
    pp.query("stop_time", stop_time);

    pp.query("cfl", cfl);

    pp.query("regrid_int", regrid_int);

    pp.query("restart", restart_chkfile);

    pp.query("check_file", check_file);
    pp.query("plot_file", plot_file);
    pp.query("check_int", check_int);
    pp.query("plot_int", plot_int);
}

void
AmrAdv::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
	BoxLib::average_down(*phi_new[lev+1], *phi_new[lev],
			     geom[lev+1], geom[lev],
			     0, phi_new[lev]->nComp(), refRatio(lev));
    }
}

