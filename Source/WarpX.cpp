
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
#include <WarpX_f.H>
#include <WarpXConst.H>
#include <WarpXWrappers.h>
#include <WarpXUtil.H>

using namespace amrex;

Vector<Real> WarpX::B_external(3, 0.0);

int WarpX::do_moving_window = 0;
int WarpX::moving_window_dir = -1;

Real WarpX::gamma_boost = 1.;
Real WarpX::beta_boost = 0.;
Vector<int> WarpX::boost_direction = {0,0,0};

long WarpX::current_deposition_algo = 3;
long WarpX::charge_deposition_algo = 0;
long WarpX::field_gathering_algo = 1;
long WarpX::particle_pusher_algo = 0;
int WarpX::maxwell_fdtd_solver_id = 0;

long WarpX::nox = 1;
long WarpX::noy = 1;
long WarpX::noz = 1;

bool WarpX::use_laser         = false;
bool WarpX::use_filter        = false;
bool WarpX::serialize_ics     = false;
bool WarpX::refine_plasma     = false;

bool WarpX::do_boosted_frame_diagnostic = false;
int  WarpX::num_snapshots_lab = std::numeric_limits<int>::lowest();
Real WarpX::dt_snapshots_lab  = std::numeric_limits<Real>::lowest();
bool WarpX::do_boosted_frame_fields = true;
bool WarpX::do_boosted_frame_particles = true;

bool WarpX::do_dynamic_scheduling = true;

#if (AMREX_SPACEDIM == 3)
IntVect WarpX::Bx_nodal_flag(1,0,0);
IntVect WarpX::By_nodal_flag(0,1,0);
IntVect WarpX::Bz_nodal_flag(0,0,1);
#elif (AMREX_SPACEDIM == 2)
IntVect WarpX::Bx_nodal_flag(1,0);  // x is the first dimension to AMReX
IntVect WarpX::By_nodal_flag(0,0);  // y is the missing dimension to 2D AMReX
IntVect WarpX::Bz_nodal_flag(0,1);  // z is the second dimension to 2D AMReX
#endif

#if (AMREX_SPACEDIM == 3)
IntVect WarpX::Ex_nodal_flag(0,1,1);
IntVect WarpX::Ey_nodal_flag(1,0,1);
IntVect WarpX::Ez_nodal_flag(1,1,0);
#elif (AMREX_SPACEDIM == 2)
IntVect WarpX::Ex_nodal_flag(0,1);  // x is the first dimension to AMReX
IntVect WarpX::Ey_nodal_flag(1,1);  // y is the missing dimension to 2D AMReX
IntVect WarpX::Ez_nodal_flag(1,0);  // z is the second dimension to 2D AMReX
#endif

#if (AMREX_SPACEDIM == 3)
IntVect WarpX::jx_nodal_flag(0,1,1);
IntVect WarpX::jy_nodal_flag(1,0,1);
IntVect WarpX::jz_nodal_flag(1,1,0);
#elif (AMREX_SPACEDIM == 2)
IntVect WarpX::jx_nodal_flag(0,1);  // x is the first dimension to AMReX
IntVect WarpX::jy_nodal_flag(1,1);  // y is the missing dimension to 2D AMReX
IntVect WarpX::jz_nodal_flag(1,0);  // z is the second dimension to 2D AMReX
#endif

int WarpX::n_field_gather_buffer = 0;
int WarpX::n_current_deposition_buffer = -1;

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

    if (do_plasma_injection) {
        for (int i = 0; i < num_injected_species; ++i) {
            int ispecies = injected_plasma_species[i];
            WarpXParticleContainer& pc = mypc->GetParticleContainer(ispecies);
            auto& ppc = dynamic_cast<PhysicalParticleContainer&>(pc);
            ppc.injected = true;
        }
    }

    Efield_aux.resize(nlevs_max);
    Bfield_aux.resize(nlevs_max);

    F_fp.resize(nlevs_max);
    rho_fp.resize(nlevs_max);
    current_fp.resize(nlevs_max);
    Efield_fp.resize(nlevs_max);
    Bfield_fp.resize(nlevs_max);

    F_cp.resize(nlevs_max);
    rho_cp.resize(nlevs_max);
    current_cp.resize(nlevs_max);
    Efield_cp.resize(nlevs_max);
    Bfield_cp.resize(nlevs_max);

    Efield_cax.resize(nlevs_max);
    Bfield_cax.resize(nlevs_max);
    current_buffer_masks.resize(nlevs_max);
    gather_buffer_masks.resize(nlevs_max);
    current_buf.resize(nlevs_max);

    pml.resize(nlevs_max);

#ifdef WARPX_DO_ELECTROSTATIC
    masks.resize(nlevs_max);
    gather_masks.resize(nlevs_max);
#endif // WARPX_DO_ELECTROSTATIC

    costs.resize(nlevs_max);

#ifdef WARPX_USE_PSATD
    Efield_fp_fft.resize(nlevs_max);
    Bfield_fp_fft.resize(nlevs_max);
    current_fp_fft.resize(nlevs_max);
    rho_fp_fft.resize(nlevs_max);

    Efield_cp_fft.resize(nlevs_max);
    Bfield_cp_fft.resize(nlevs_max);
    current_cp_fft.resize(nlevs_max);
    rho_cp_fft.resize(nlevs_max);

    dataptr_fp_fft.resize(nlevs_max);
    dataptr_cp_fft.resize(nlevs_max);

    ba_valid_fp_fft.resize(nlevs_max);
    ba_valid_cp_fft.resize(nlevs_max);

    domain_fp_fft.resize(nlevs_max);
    domain_cp_fft.resize(nlevs_max);

    comm_fft.resize(nlevs_max,MPI_COMM_NULL);
    color_fft.resize(nlevs_max,-1);
#endif
}

WarpX::~WarpX ()
{
    int nlevs_max = maxLevel() +1;
    for (int lev = 0; lev < nlevs_max; ++lev) {
        ClearLevel(lev);
    }
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

        ReadBoostedFrameParameters(gamma_boost, beta_boost, boost_direction);

        pp.queryarr("B_external", B_external);

	pp.query("do_moving_window", do_moving_window);
	if (do_moving_window)
	{
	    std::string s;
	    pp.get("moving_window_dir", s);
	    if (s == "x" || s == "X") {
		moving_window_dir = 0;
	    }
#if (AMREX_SPACEDIM == 3)
	    else if (s == "y" || s == "Y") {
		moving_window_dir = 1;
	    }
#endif
	    else if (s == "z" || s == "Z") {
		moving_window_dir = AMREX_SPACEDIM-1;
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
          if (moving_window_v >= 0){
              // Inject particles continuously from the right end of the box
              current_injection_position = geom[0].ProbHi(moving_window_dir);
          } else {
              // Inject particles continuously from the left end of the box
              current_injection_position = geom[0].ProbLo(moving_window_dir);
          }
	}

        pp.query("do_boosted_frame_diagnostic", do_boosted_frame_diagnostic);
        if (do_boosted_frame_diagnostic) {

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(gamma_boost > 1.0,
                "gamma_boost must be > 1 to use the boosted frame diagnostic.");

            std::string s;
    	    pp.get("boost_direction", s);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( (s == "z" || s == "Z"),
                "The boosted frame diagnostic currently only works if the boost is in the z direction.");

            pp.get("num_snapshots_lab", num_snapshots_lab);
            pp.get("dt_snapshots_lab", dt_snapshots_lab);
            pp.get("gamma_boost", gamma_boost);

            pp.query("do_boosted_frame_fields", do_boosted_frame_fields);
            pp.query("do_boosted_frame_particles", do_boosted_frame_particles);


            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(do_moving_window,
                "The moving window should be on if using the boosted frame diagnostic.");

	        pp.get("moving_window_dir", s);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( (s == "z" || s == "Z"),
                "The boosted frame diagnostic currently only works if the moving window is in the z direction.");
        }

        pp.query("do_electrostatic", do_electrostatic);
        pp.query("n_buffer", n_buffer);
        pp.query("const_dt", const_dt);

	pp.query("use_laser", use_laser);
	pp.query("use_filter", use_filter);
	pp.query("serialize_ics", serialize_ics);
	pp.query("refine_plasma", refine_plasma);
        pp.query("do_dive_cleaning", do_dive_cleaning);
        pp.query("n_field_gather_buffer", n_field_gather_buffer);
        pp.query("n_current_deposition_buffer", n_current_deposition_buffer);

        pp.query("do_pml", do_pml);
        pp.query("pml_ncell", pml_ncell);
        pp.query("pml_delta", pml_delta);

        pp.query("plot_raw_fields", plot_raw_fields);
        pp.query("plot_raw_fields_guards", plot_raw_fields_guards);
        if (ParallelDescriptor::NProcs() == 1) {
            plot_proc_number = false;
        }
        pp.query("plot_part_per_cell", plot_part_per_cell);
        pp.query("plot_part_per_grid", plot_part_per_grid);
        pp.query("plot_part_per_proc", plot_part_per_proc);
        pp.query("plot_proc_number"  , plot_proc_number);
        pp.query("plot_dive"         , plot_dive);
        pp.query("plot_divb"         , plot_divb);

        if (maxLevel() > 0) {
            pp.query("plot_finepatch", plot_finepatch);
            pp.query("plot_crsepatch", plot_crsepatch);
        }

        {
            bool plotfile_min_max = true;
            pp.query("plotfile_min_max", plotfile_min_max);
            if (plotfile_min_max) {
                plotfile_headerversion = amrex::VisMF::Header::Version_v1;
            } else {
                plotfile_headerversion = amrex::VisMF::Header::NoFabHeader_v1;
            }
            pp.query("usesingleread", use_single_read);
            pp.query("usesinglewrite", use_single_write);
            ParmParse ppv("vismf");
            ppv.add("usesingleread", use_single_read);
            ppv.add("usesinglewrite", use_single_write);
            pp.query("mffile_nstreams", mffile_nstreams);
            VisMF::SetMFFileInStreams(mffile_nstreams);
            pp.query("field_io_nfiles", field_io_nfiles);
            VisMF::SetNOutFiles(field_io_nfiles);
            pp.query("particle_io_nfiles", particle_io_nfiles);
            ParmParse ppp("particles");
            ppp.add("particles_nfiles", particle_io_nfiles);
        }

        if (maxLevel() > 0) {
            Vector<Real> lo, hi;
            pp.getarr("fine_tag_lo", lo);
            pp.getarr("fine_tag_hi", hi);
            fine_tag_lo = RealVect{lo};
            fine_tag_hi = RealVect{hi};
        }

        pp.query("load_balance_int", load_balance_int);
        pp.query("load_balance_with_sfc", load_balance_with_sfc);
        pp.query("load_balance_knapsack_factor", load_balance_knapsack_factor);

        pp.query("do_dynamic_scheduling", do_dynamic_scheduling);
    }

    {
	ParmParse pp("interpolation");
	pp.query("nox", nox);
	pp.query("noy", noy);
	pp.query("noz", noz);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( nox == noy and nox == noz ,
	    "warpx.nox, noy and noz must be equal");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( nox >= 1, "warpx.nox must >= 1");
    }

    {
	ParmParse pp("algo");
	pp.query("current_deposition", current_deposition_algo);
	pp.query("charge_deposition", charge_deposition_algo);
	pp.query("field_gathering", field_gathering_algo);
	pp.query("particle_pusher", particle_pusher_algo);
	std::string s_solver = "";
	pp.query("maxwell_fdtd_solver", s_solver);
    std::transform(s_solver.begin(),
                   s_solver.end(),
                   s_solver.begin(),
                   ::tolower);
	// if maxwell_fdtd_solver is specified, set the value
	// of maxwell_fdtd_solver_id accordingly.
      // Otherwise keep the default value maxwell_fdtd_solver_id=0
      if (s_solver != "") {
        if (s_solver == "yee") {
            maxwell_fdtd_solver_id = 0;
        } else if (s_solver == "ckc") {
            maxwell_fdtd_solver_id = 1;
        } else {
            amrex::Abort("Unknown FDTD Solver type " + s_solver);
        }
      }
    }

#ifdef WARPX_USE_PSATD
    {
        ParmParse pp("psatd");
        pp.query("ngroups_fft", ngroups_fft);
        pp.query("fftw_plan_measure", fftw_plan_measure);
        pp.query("nox", nox_fft);
        pp.query("noy", noy_fft);
        pp.query("noz", noz_fft);
    }
#endif
}

// This is a virtual function.
void
WarpX::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& new_grids,
                                const DistributionMapping& new_dmap)
{
    AllocLevelData(lev, new_grids, new_dmap);
    InitLevelData(lev, time);

#ifdef WARPX_USE_PSATD
    AllocLevelDataFFT(lev);
    InitLevelDataFFT(lev, time);
#endif
}

void
WarpX::ClearLevel (int lev)
{
    for (int i = 0; i < 3; ++i) {
	Efield_aux[lev][i].reset();
	Bfield_aux[lev][i].reset();

	current_fp[lev][i].reset();
	Efield_fp [lev][i].reset();
	Bfield_fp [lev][i].reset();

	current_cp[lev][i].reset();
	Efield_cp [lev][i].reset();
	Bfield_cp [lev][i].reset();

	Efield_cax[lev][i].reset();
	Bfield_cax[lev][i].reset();
        current_buf[lev][i].reset();
    }

    current_buffer_masks[lev].reset();
    gather_buffer_masks[lev].reset();

    F_fp  [lev].reset();
    rho_fp[lev].reset();
    F_cp  [lev].reset();
    rho_cp[lev].reset();

    costs[lev].reset();

#ifdef WARPX_USE_PSATD
    for (int i = 0; i < 3; ++i) {
        Efield_fp_fft[lev][i].reset();
        Bfield_fp_fft[lev][i].reset();
        current_fp_fft[lev][i].reset();

        Efield_cp_fft[lev][i].reset();
        Bfield_cp_fft[lev][i].reset();
        current_cp_fft[lev][i].reset();
    }

    rho_fp_fft[lev].reset();
    rho_cp_fft[lev].reset();

    dataptr_fp_fft[lev].reset();
    dataptr_cp_fft[lev].reset();

    ba_valid_fp_fft[lev] = BoxArray();
    ba_valid_cp_fft[lev] = BoxArray();

    FreeFFT(lev);
#endif
}

void
WarpX::AllocLevelData (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
    // Ex, Ey, Ez, Bx, By, and Bz have the same number of ghost cells.
    // jx, jy, jz and rho have the same number of ghost cells.
    // E and B have the same number of ghost cells as j and rho if NCI filter is not used,
    // but different number of ghost cells in z-direction if NCI filter is used.
    int ngx = (WarpX::nox % 2) ? WarpX::nox+1 : WarpX::nox;  // Always even number
    int ngy = (WarpX::noy % 2) ? WarpX::noy+1 : WarpX::noy;  // Always even number
    int ngz_nonci = (WarpX::noz % 2) ? WarpX::noz+1 : WarpX::noz;  // Always even number
    int ngz;
    if (warpx_use_fdtd_nci_corr()) {
        int ng = WarpX::noz + (mypc->nstencilz_fdtd_nci_corr-1);
        ngz = (ng % 2) ? ng+1 : ng;
    } else {
        ngz = ngz_nonci;
    }

#if (AMREX_SPACEDIM == 3)
    IntVect ngE(ngx,ngy,ngz);
    IntVect ngJ(ngx,ngy,ngz_nonci);
#elif (AMREX_SPACEDIM == 2)
    IntVect ngE(ngx,ngz);
    IntVect ngJ(ngx,ngz_nonci);
#endif

    IntVect ngRho = ngJ+1; //One extra ghost cell, so that it's safe to deposit charge density
                           // after pushing particle. 

    if (mypc->nSpeciesDepositOnMainGrid() && n_current_deposition_buffer == 0) {
        n_current_deposition_buffer = 1;
    }

    if (n_current_deposition_buffer < 0) {
        n_current_deposition_buffer = ngJ.max();
    }

    int ngF = (do_moving_window) ? 2 : 0;
    // CKC solver requires one additional guard cell
    if (maxwell_fdtd_solver_id == 1) ngF = std::max( ngF, 1 );

    //
    // The fine patch
    //
    Bfield_fp[lev][0].reset( new MultiFab(amrex::convert(ba,Bx_nodal_flag),dm,1,ngE));
    Bfield_fp[lev][1].reset( new MultiFab(amrex::convert(ba,By_nodal_flag),dm,1,ngE));
    Bfield_fp[lev][2].reset( new MultiFab(amrex::convert(ba,Bz_nodal_flag),dm,1,ngE));

    Efield_fp[lev][0].reset( new MultiFab(amrex::convert(ba,Ex_nodal_flag),dm,1,ngE));
    Efield_fp[lev][1].reset( new MultiFab(amrex::convert(ba,Ey_nodal_flag),dm,1,ngE));
    Efield_fp[lev][2].reset( new MultiFab(amrex::convert(ba,Ez_nodal_flag),dm,1,ngE));

    current_fp[lev][0].reset( new MultiFab(amrex::convert(ba,jx_nodal_flag),dm,1,ngJ));
    current_fp[lev][1].reset( new MultiFab(amrex::convert(ba,jy_nodal_flag),dm,1,ngJ));
    current_fp[lev][2].reset( new MultiFab(amrex::convert(ba,jz_nodal_flag),dm,1,ngJ));

    if (do_dive_cleaning)
    {
        F_fp[lev].reset  (new MultiFab(amrex::convert(ba,IntVect::TheUnitVector()),dm,1, ngF));
        rho_fp[lev].reset(new MultiFab(amrex::convert(ba,IntVect::TheUnitVector()),dm,2,ngRho));
    }
#ifdef WARPX_USE_PSATD
    else
    {
        rho_fp[lev].reset(new MultiFab(amrex::convert(ba,IntVect::TheUnitVector()),dm,2,ngRho));
    }
#endif

    //
    // The Aux patch (i.e., the full solution)
    //
    if (lev == 0)
    {
        for (int idir = 0; idir < 3; ++idir) {
            Efield_aux[lev][idir].reset(new MultiFab(*Efield_fp[lev][idir], amrex::make_alias, 0, 1));
            Bfield_aux[lev][idir].reset(new MultiFab(*Bfield_fp[lev][idir], amrex::make_alias, 0, 1));
        }
    }
    else
    {
        Bfield_aux[lev][0].reset( new MultiFab(amrex::convert(ba,Bx_nodal_flag),dm,1,ngE));
        Bfield_aux[lev][1].reset( new MultiFab(amrex::convert(ba,By_nodal_flag),dm,1,ngE));
        Bfield_aux[lev][2].reset( new MultiFab(amrex::convert(ba,Bz_nodal_flag),dm,1,ngE));

        Efield_aux[lev][0].reset( new MultiFab(amrex::convert(ba,Ex_nodal_flag),dm,1,ngE));
        Efield_aux[lev][1].reset( new MultiFab(amrex::convert(ba,Ey_nodal_flag),dm,1,ngE));
        Efield_aux[lev][2].reset( new MultiFab(amrex::convert(ba,Ez_nodal_flag),dm,1,ngE));
    }

    //
    // The coarse patch
    //
    if (lev > 0)
    {
        BoxArray cba = ba;
        cba.coarsen(refRatio(lev-1));

        // Create the MultiFabs for B
        Bfield_cp[lev][0].reset( new MultiFab(amrex::convert(cba,Bx_nodal_flag),dm,1,ngE));
        Bfield_cp[lev][1].reset( new MultiFab(amrex::convert(cba,By_nodal_flag),dm,1,ngE));
        Bfield_cp[lev][2].reset( new MultiFab(amrex::convert(cba,Bz_nodal_flag),dm,1,ngE));

        // Create the MultiFabs for E
        Efield_cp[lev][0].reset( new MultiFab(amrex::convert(cba,Ex_nodal_flag),dm,1,ngE));
        Efield_cp[lev][1].reset( new MultiFab(amrex::convert(cba,Ey_nodal_flag),dm,1,ngE));
        Efield_cp[lev][2].reset( new MultiFab(amrex::convert(cba,Ez_nodal_flag),dm,1,ngE));

        // Create the MultiFabs for the current
        current_cp[lev][0].reset( new MultiFab(amrex::convert(cba,jx_nodal_flag),dm,1,ngJ));
        current_cp[lev][1].reset( new MultiFab(amrex::convert(cba,jy_nodal_flag),dm,1,ngJ));
        current_cp[lev][2].reset( new MultiFab(amrex::convert(cba,jz_nodal_flag),dm,1,ngJ));

        if (do_dive_cleaning)
        {
            F_cp[lev].reset  (new MultiFab(amrex::convert(cba,IntVect::TheUnitVector()),dm,1, ngF));
            rho_cp[lev].reset(new MultiFab(amrex::convert(cba,IntVect::TheUnitVector()),dm,2,ngRho));
        }
#ifdef WARPX_USE_PSATD
        else
        {
            rho_cp[lev].reset(new MultiFab(amrex::convert(cba,IntVect::TheUnitVector()),dm,2,ngRho));
        }
#endif
    }

    //
    // Copy of the coarse aux
    //
    if (lev > 0 && (n_field_gather_buffer > 0 || n_current_deposition_buffer > 0))
    {
        BoxArray cba = ba;
        cba.coarsen(refRatio(lev-1));

        if (n_field_gather_buffer > 0) {
            // Create the MultiFabs for B
            Bfield_cax[lev][0].reset( new MultiFab(amrex::convert(cba,Bx_nodal_flag),dm,1,ngE));
            Bfield_cax[lev][1].reset( new MultiFab(amrex::convert(cba,By_nodal_flag),dm,1,ngE));
            Bfield_cax[lev][2].reset( new MultiFab(amrex::convert(cba,Bz_nodal_flag),dm,1,ngE));

            // Create the MultiFabs for E
            Efield_cax[lev][0].reset( new MultiFab(amrex::convert(cba,Ex_nodal_flag),dm,1,ngE));
            Efield_cax[lev][1].reset( new MultiFab(amrex::convert(cba,Ey_nodal_flag),dm,1,ngE));
            Efield_cax[lev][2].reset( new MultiFab(amrex::convert(cba,Ez_nodal_flag),dm,1,ngE));

            gather_buffer_masks[lev].reset( new iMultiFab(ba, dm, 1, 0) );
        }

        if (n_current_deposition_buffer > 0) {
            current_buf[lev][0].reset( new MultiFab(amrex::convert(cba,jx_nodal_flag),dm,1,ngJ));
            current_buf[lev][1].reset( new MultiFab(amrex::convert(cba,jy_nodal_flag),dm,1,ngJ));
            current_buf[lev][2].reset( new MultiFab(amrex::convert(cba,jz_nodal_flag),dm,1,ngJ));
            current_buffer_masks[lev].reset( new iMultiFab(ba, dm, 1, 0) );
        }
    }

    if (load_balance_int > 0) {
        costs[lev].reset(new MultiFab(ba, dm, 1, 0));
    }
}

std::array<Real,3>
WarpX::CellSize (int lev)
{
    const auto& gm = GetInstance().Geom(lev);
    const Real* dx = gm.CellSize();
#if (AMREX_SPACEDIM == 3)
    return { dx[0], dx[1], dx[2] };
#elif (AMREX_SPACEDIM == 2)
    return { dx[0], 1.0, dx[1] };
#else
    static_assert(AMREX_SPACEDIM != 1, "1D is not supported");
#endif
}

amrex::RealBox
WarpX::getRealBox(const Box& bx, int lev)
{
    const auto& gm = GetInstance().Geom(lev);
    const RealBox grid_box{bx, gm.CellSize(), gm.ProbLo()};
    return( grid_box );
}

std::array<Real,3>
WarpX::LowerCorner(const Box& bx, int lev)
{
    const RealBox grid_box = getRealBox( bx, lev );
    const Real* xyzmin = grid_box.lo();
#if (AMREX_SPACEDIM == 3)
    return { xyzmin[0], xyzmin[1], xyzmin[2] };
#elif (AMREX_SPACEDIM == 2)
    return { xyzmin[0], -1.e100, xyzmin[1] };
#endif
}

std::array<Real,3>
WarpX::UpperCorner(const Box& bx, int lev)
{
    const RealBox grid_box = getRealBox( bx, lev );
    const Real* xyzmax = grid_box.hi();
#if (AMREX_SPACEDIM == 3)
    return { xyzmax[0], xyzmax[1], xyzmax[2] };
#elif (AMREX_SPACEDIM == 2)
    return { xyzmax[0], 1.e100, xyzmax[1] };
#endif
}

IntVect
WarpX::RefRatio (int lev)
{
    return GetInstance().refRatio(lev);
}

void
WarpX::ComputeDivB (MultiFab& divB, int dcomp,
                    const std::array<const MultiFab*, 3>& B,
                    const std::array<Real,3>& dx)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(divB, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        WRPX_COMPUTE_DIVB(bx.loVect(), bx.hiVect(),
                           BL_TO_FORTRAN_N_ANYD(divB[mfi],dcomp),
                           BL_TO_FORTRAN_ANYD((*B[0])[mfi]),
                           BL_TO_FORTRAN_ANYD((*B[1])[mfi]),
                           BL_TO_FORTRAN_ANYD((*B[2])[mfi]),
                           dx.data());
    }
}

void
WarpX::ComputeDivE (MultiFab& divE, int dcomp,
                    const std::array<const MultiFab*, 3>& E,
                    const std::array<Real,3>& dx)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(divE, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        WRPX_COMPUTE_DIVE(bx.loVect(), bx.hiVect(),
                           BL_TO_FORTRAN_N_ANYD(divE[mfi],dcomp),
                           BL_TO_FORTRAN_ANYD((*E[0])[mfi]),
                           BL_TO_FORTRAN_ANYD((*E[1])[mfi]),
                           BL_TO_FORTRAN_ANYD((*E[2])[mfi]),
                           dx.data());
    }
}

void
WarpX::applyFilter (MultiFab& dstmf, const MultiFab& srcmf)
{
    const int ncomp = srcmf.nComp();
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox tmpfab;
        for (MFIter mfi(dstmf,true); mfi.isValid(); ++mfi)
        {
            const auto& srcfab = srcmf[mfi];
            auto& dstfab = dstmf[mfi];
            const Box& tbx = mfi.growntilebox();
            const Box& gbx = amrex::grow(tbx,1);
            tmpfab.resize(gbx,ncomp);
            tmpfab.setVal(0.0, gbx, 0, ncomp);
            const Box& ibx = gbx & srcfab.box();
            tmpfab.copy(srcfab, ibx, 0, ibx, 0, ncomp);
            WRPX_FILTER(BL_TO_FORTRAN_BOX(tbx),
                        BL_TO_FORTRAN_ANYD(tmpfab),
                        BL_TO_FORTRAN_ANYD(dstfab),
                        ncomp);
        }
    }
}

void
WarpX::BuildBufferMasks ()
{
    int ngbuffer = std::max(n_field_gather_buffer, n_current_deposition_buffer);
    for (int lev = 1; lev <= maxLevel(); ++lev)
    {
        for (int ipass = 0; ipass < 2; ++ipass)
        {
            int ngbuffer = (ipass == 0) ? n_current_deposition_buffer : n_field_gather_buffer;
            iMultiFab* bmasks = (ipass == 0) ? current_buffer_masks[lev].get() : gather_buffer_masks[lev].get();
            if (bmasks)
            {
                iMultiFab tmp(bmasks->boxArray(), bmasks->DistributionMap(), 1, ngbuffer);
                const int covered = 1;
                const int notcovered = 0;
                const int physbnd = 1;
                const int interior = 1;
                const Box& dom = Geom(lev).Domain();
                const auto& period = Geom(lev).periodicity();
                tmp.BuildMask(dom, period, covered, notcovered, physbnd, interior);
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (MFIter mfi(*bmasks, true); mfi.isValid(); ++mfi)
                {
                    const Box& tbx = mfi.tilebox();
                    warpx_build_buffer_masks (BL_TO_FORTRAN_BOX(tbx),
                                              BL_TO_FORTRAN_ANYD((*bmasks)[mfi]),
                                              BL_TO_FORTRAN_ANYD(tmp[mfi]),
                                              &ngbuffer);
                }
            }
        }
    }
}

const iMultiFab*
WarpX::CurrentBufferMasks (int lev)
{
    return GetInstance().getCurrentBufferMasks(lev);
}

const iMultiFab*
WarpX::GatherBufferMasks (int lev)
{
    return GetInstance().getGatherBufferMasks(lev);
}
