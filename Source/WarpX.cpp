
#include <limits>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>

#include <WarpX.H>
#include <WarpXConst.H>

#include <WarpX_f.H>

using namespace amrex;

long WarpX::current_deposition_algo = 3;
long WarpX::charge_deposition_algo = 0;
long WarpX::field_gathering_algo = 1;
long WarpX::particle_pusher_algo = 0;

long WarpX::nox = 1;
long WarpX::noy = 1;
long WarpX::noz = 1;

bool WarpX::use_laser         = false;

#if (BL_SPACEDIM == 3)
IntVect WarpX::Bx_nodal_flag(1,0,0);
IntVect WarpX::By_nodal_flag(0,1,0);
IntVect WarpX::Bz_nodal_flag(0,0,1);
#elif (BL_SPACEDIM == 2)
IntVect WarpX::Bx_nodal_flag(1,0);  // x is the first dimension to AMReX
IntVect WarpX::By_nodal_flag(0,0);  // y is the missing dimension to 2D AMReX
IntVect WarpX::Bz_nodal_flag(0,1);  // z is the second dimension to 2D AMReX
#endif

#if (BL_SPACEDIM == 3)
IntVect WarpX::Ex_nodal_flag(0,1,1);
IntVect WarpX::Ey_nodal_flag(1,0,1);
IntVect WarpX::Ez_nodal_flag(1,1,0);
#elif (BL_SPACEDIM == 2)
IntVect WarpX::Ex_nodal_flag(0,1);  // x is the first dimension to AMReX
IntVect WarpX::Ey_nodal_flag(1,1);  // y is the missing dimension to 2D AMReX
IntVect WarpX::Ez_nodal_flag(1,0);  // z is the second dimension to 2D AMReX
#endif

#if (BL_SPACEDIM == 3)
IntVect WarpX::jx_nodal_flag(0,1,1);
IntVect WarpX::jy_nodal_flag(1,0,1);
IntVect WarpX::jz_nodal_flag(1,1,0);
#elif (BL_SPACEDIM == 2)
IntVect WarpX::jx_nodal_flag(0,1);  // x is the first dimension to AMReX
IntVect WarpX::jy_nodal_flag(1,1);  // y is the missing dimension to 2D AMReX
IntVect WarpX::jz_nodal_flag(1,0);  // z is the second dimension to 2D AMReX
#endif

WarpX* WarpX::m_instance = nullptr;

WarpX&
WarpX::GetInstance ()
{
    if (!m_instance) {
	m_instance = new WarpX();
    }
    return *m_instance;
}

void
WarpX::ResetInstance ()
{
    delete m_instance;
    m_instance = nullptr;
}

WarpX::WarpX ()
{
    m_instance = this;

    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = maxLevel() + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
#if 0
    // no subcycling yet
    for (int lev = 1; lev <= maxLevel(); ++lev) {
	nsubsteps[lev] = MaxRefRatio(lev-1);
    }
#endif

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    // Particle Container
    mypc = std::unique_ptr<MultiParticleContainer> (new MultiParticleContainer(this));

    current.resize(nlevs_max);
    Efield.resize(nlevs_max);
    Bfield.resize(nlevs_max);
    cfbndry.resize(nlevs_max-1);
}

WarpX::~WarpX ()
{
}

void
WarpX::ReadParameters ()
{
    {
	ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
	pp.query("max_step", max_step);
	pp.query("stop_time", stop_time);
    }

    {
	ParmParse pp("amr"); // Traditionally, these have prefix, amr.

	pp.query("check_file", check_file);
	pp.query("check_int", check_int);

	pp.query("plot_file", plot_file);
	pp.query("plot_int", plot_int);

	pp.query("restart", restart_chkfile);
    }

    {
	ParmParse pp("warpx");

	pp.query("cfl", cfl);
	pp.query("verbose", verbose);
	pp.query("regrid_int", regrid_int);

        // PML
        if (Geometry::isAllPeriodic()) {
            do_pml = 0;  // no PML for all periodic boundaries
        } else {
            pp.query("do_pml", do_pml);
            pp.query("pml_ncell", pml_ncell);
        }

	pp.query("do_moving_window", do_moving_window);
	if (do_moving_window)
	{
	    std::string s;
	    pp.get("moving_window_dir", s);
	    if (s == "x" || s == "X") {
		moving_window_dir = 0;
	    }
#if (BL_SPACEDIM == 3)
	    else if (s == "y" || s == "Y") {
		moving_window_dir = 1;
	    }
#endif
	    else if (s == "z" || s == "Z") {
		moving_window_dir = BL_SPACEDIM-1;
	    }
	    else {
		const std::string msg = "Unknown moving_window_dir: "+s;
		amrex::Abort(msg.c_str());
	    }

	    moving_window_x = geom[0].ProbLo(moving_window_dir);

	    pp.get("moving_window_v", moving_window_v);
	    moving_window_v *= PhysConst::c;
	}

	pp.query("do_plasma_injection", do_plasma_injection);
	if (do_plasma_injection) {
	  pp.get("num_injected_species", num_injected_species);
	  injected_plasma_species.resize(num_injected_species);
	  pp.getarr("injected_plasma_species", injected_plasma_species,
		 0, num_injected_species);
	}

	pp.query("use_laser", use_laser);

        pp.query("plot_raw_fields", plot_raw_fields);
        if (ParallelDescriptor::NProcs() == 1) {
            plot_proc_number = false;
        }
        pp.query("plot_part_per_cell", plot_part_per_cell);
        pp.query("plot_part_per_grid", plot_part_per_grid);
        pp.query("plot_part_per_proc", plot_part_per_proc);
        pp.query("plot_proc_number"  , plot_proc_number);
        pp.query("plot_divb"         , plot_divb);
    }

    {
	ParmParse pp("interpolation");
	pp.query("nox", nox);
	pp.query("noy", noy);
	pp.query("noz", noz);
	if (nox != noy || nox != noz) {
	    amrex::Abort("warpx.nox, noy and noz must be equal");
	}
	if (nox < 1) {
	    amrex::Abort("warpx.nox must >= 1");
	}
    }

    {
	ParmParse pp("algo");
	pp.query("current_deposition", current_deposition_algo);
	pp.query("charge_deposition", charge_deposition_algo);
	pp.query("field_gathering", field_gathering_algo);
	pp.query("particle_pusher", particle_pusher_algo);
    }
}

// This is a virtual function.
void
WarpX::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& new_grids,
                                const DistributionMapping& new_dmap)
{
    AllocLevelData(lev, new_grids, new_dmap);
    InitLevelData(lev, time);
}

void
WarpX::ClearLevel (int lev)
{
    for (int i = 0; i < 3; ++i) {
	current[lev][i].reset();
	Efield [lev][i].reset();
	Bfield [lev][i].reset();
    }
    if (lev < max_level) {
        cfbndry[lev].reset();
    }
}

void
WarpX::AllocLevelData (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
    const int ng = WarpX::nox;  // need to update this

    // Create the MultiFabs for B
    Bfield[lev][0].reset( new MultiFab(amrex::convert(ba,Bx_nodal_flag),dm,1,ng));
    Bfield[lev][1].reset( new MultiFab(amrex::convert(ba,By_nodal_flag),dm,1,ng));
    Bfield[lev][2].reset( new MultiFab(amrex::convert(ba,Bz_nodal_flag),dm,1,ng));

    // Create the MultiFabs for E
    Efield[lev][0].reset( new MultiFab(amrex::convert(ba,Ex_nodal_flag),dm,1,ng));
    Efield[lev][1].reset( new MultiFab(amrex::convert(ba,Ey_nodal_flag),dm,1,ng));
    Efield[lev][2].reset( new MultiFab(amrex::convert(ba,Ez_nodal_flag),dm,1,ng));

    // Create the MultiFabs for the current
    current[lev][0].reset( new MultiFab(amrex::convert(ba,jx_nodal_flag),dm,1,ng));
    current[lev][1].reset( new MultiFab(amrex::convert(ba,jy_nodal_flag),dm,1,ng));
    current[lev][2].reset( new MultiFab(amrex::convert(ba,jz_nodal_flag),dm,1,ng));
}

std::array<Real,3>
WarpX::CellSize (int lev)
{
    const auto& gm = GetInstance().Geom(lev);
    const Real* dx = gm.CellSize();
#if (BL_SPACEDIM == 3)
    return { dx[0], dx[1], dx[2] };
#elif (BL_SPACEDIM == 2)
    return { dx[0], 1.0, dx[1] };
#else
    static_assert(BL_SPACEDIM != 1, "1D is not supported");
#endif
}

std::array<Real,3>
WarpX::LowerCorner(const Box& bx, int lev)
{
    const auto& gm = GetInstance().Geom(lev);
    const RealBox grid_box{bx, gm.CellSize(), gm.ProbLo()};
    const Real* xyzmin = grid_box.lo();
#if (BL_SPACEDIM == 3)
    return { xyzmin[0], xyzmin[1], xyzmin[2] };
#elif (BL_SPACEDIM == 2)
    return { xyzmin[0], -1.e100, xyzmin[1] };
#endif
}

std::array<Real,3>
WarpX::UpperCorner(const Box& bx, int lev)
{
    const auto& gm = GetInstance().Geom(lev);
    const RealBox grid_box{bx, gm.CellSize(), gm.ProbLo()};
    const Real* xyzmax = grid_box.hi();
#if (BL_SPACEDIM == 3)
    return { xyzmax[0], xyzmax[1], xyzmax[2] };
#elif (BL_SPACEDIM == 2)
    return { xyzmax[0], 1.e100, xyzmax[1] };
#endif
}

void
WarpX::ComputeDivB (MultiFab& divB, int dcomp, const Array<const MultiFab*>& B, const Real* dx)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(divB, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        WARPX_COMPUTE_DIVB(bx.loVect(), bx.hiVect(),
                           BL_TO_FORTRAN_N_ANYD(divB[mfi],dcomp),
                           D_DECL(BL_TO_FORTRAN_ANYD((*B[0])[mfi]),
                                  BL_TO_FORTRAN_ANYD((*B[1])[mfi]),
                                  BL_TO_FORTRAN_ANYD((*B[2])[mfi])),
                           dx);
    }
}
