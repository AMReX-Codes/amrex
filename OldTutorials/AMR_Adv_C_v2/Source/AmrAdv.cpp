
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>

#include <AmrAdv.H>
#include <AmrAdvBC.H>

using namespace amrex;

AmrAdv::AmrAdv ()
{
    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = maxLevel() + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= maxLevel(); ++lev) {
	nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);

    flux_reg.resize(nlevs_max+1);
}

AmrAdv::~AmrAdv ()
{

}

void
AmrAdv::ReadParameters ()
{
    {
	ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
	pp.query("max_step", max_step);
	pp.query("stop_time", stop_time);
    }

    {
	ParmParse pp("amr"); // Traditionally, these have prefix, amr.

	pp.query("regrid_int", regrid_int);

	pp.query("check_file", check_file);
	pp.query("check_int", check_int);

	pp.query("plot_file", plot_file);
	pp.query("plot_int", plot_int);

	pp.query("restart", restart_chkfile);
    }

    {
	ParmParse pp("adv");
	
	pp.query("cfl", cfl);
	
	pp.query("do_reflux", do_reflux);
    }
}

void
AmrAdv::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev-1]->nComp();
    const int nghost = phi_new[lev-1]->nGrow();
    
    phi_new[lev].reset(new MultiFab(ba, dm, ncomp, nghost));
    phi_old[lev].reset(new MultiFab(ba, dm, ncomp, nghost));

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    FillCoarsePatch(lev, time, *phi_new[lev], 0, ncomp);
}


void
AmrAdv::RemakeLevel (int lev, Real time, const BoxArray& ba,
		     const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev]->nComp();
    const int nghost = phi_new[lev]->nGrow();

#if __cplusplus >= 201402L
    auto new_state = std::make_unique<MultiFab>(ba, dm, ncomp, nghost);
    auto old_state = std::make_unique<MultiFab>(ba, dm, ncomp, nghost);
#else
    std::unique_ptr<MultiFab> new_state(new MultiFab(ba, dm, ncomp, nghost));
    std::unique_ptr<MultiFab> old_state(new MultiFab(ba, dm, ncomp, nghost));
#endif

    FillPatch(lev, time, *new_state, 0, ncomp);

    std::swap(new_state, phi_new[lev]);
    std::swap(old_state, phi_old[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }    
}

void
AmrAdv::ClearLevel (int lev)
{
    phi_new[lev].reset(nullptr);
    phi_old[lev].reset(nullptr);
    flux_reg[lev].reset(nullptr);
}

void
AmrAdv::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
	amrex::average_down(*phi_new[lev+1], *phi_new[lev],
			     geom[lev+1], geom[lev],
			     0, phi_new[lev]->nComp(), refRatio(lev));
    }
}

void
AmrAdv::AverageDownTo (int crse_lev)
{
    amrex::average_down(*phi_new[crse_lev+1], *phi_new[crse_lev],
			 geom[crse_lev+1], geom[crse_lev],
			 0, phi_new[crse_lev]->nComp(), refRatio(crse_lev));
}

long
AmrAdv::CountCells (int lev)
{
    const int N = grids[lev].size();

    long cnt = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:cnt)
#endif
    for (int i = 0; i < N; ++i)
    {
        cnt += grids[lev][i].numPts();
    }

    return cnt;
}

void
AmrAdv::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    if (lev == 0)
    {
	Array<MultiFab*> smf;
	Array<Real> stime;
	GetData(0, time, smf, stime);

	AmrAdvPhysBC physbc;
	amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
				     geom[lev], physbc);
    }
    else
    {
	Array<MultiFab*> cmf, fmf;
	Array<Real> ctime, ftime;
	GetData(lev-1, time, cmf, ctime);
	GetData(lev  , time, fmf, ftime);

	AmrAdvPhysBC cphysbc, fphysbc;
	Interpolater* mapper = &cell_cons_interp;

	int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaryies
	int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
	Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

	amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
				   0, icomp, ncomp, geom[lev-1], geom[lev],
				   cphysbc, fphysbc, refRatio(lev-1),
				   mapper, bcs);
    }
}

void
AmrAdv::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    Array<MultiFab*> cmf;
    Array<Real> ctime;
    GetData(lev-1, time, cmf, ctime);
    
    if (cmf.size() != 1) {
	amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    AmrAdvPhysBC cphysbc, fphysbc;
    Interpolater* mapper = &cell_cons_interp;
    
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaryies
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
				 cphysbc, fphysbc, refRatio(lev-1),
				 mapper, bcs);
}

void
AmrAdv::GetData (int lev, Real time, Array<MultiFab*>& data, Array<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
	data.push_back(phi_new[lev].get());
	datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
	data.push_back(phi_old[lev].get());
	datatime.push_back(t_old[lev]);
    }
    else
    {
	data.push_back(phi_old[lev].get());
	data.push_back(phi_new[lev].get());
	datatime.push_back(t_old[lev]);
	datatime.push_back(t_new[lev]);
    }
}
