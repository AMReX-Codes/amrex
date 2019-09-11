#include <WarpX.H>
#include <WarpX_f.H>
#include <WarpXConst.H>
#include <WarpXWrappers.h>
#include <WarpXUtil.H>
#include <WarpXAlgorithmSelection.H>
#include <WarpX_FDTD.H>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#ifdef BL_USE_SENSEI_INSITU
#   include <AMReX_AmrMeshInSituBridge.H>
#endif

#ifdef _OPENMP
#   include <omp.h>
#endif

#include <limits>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <numeric>

using namespace amrex;

Vector<Real> WarpX::B_external(3, 0.0);

int WarpX::do_moving_window = 0;
int WarpX::moving_window_dir = -1;
Real WarpX::moving_window_v = std::numeric_limits<amrex::Real>::max();

Real WarpX::gamma_boost = 1.;
Real WarpX::beta_boost = 0.;
Vector<int> WarpX::boost_direction = {0,0,0};
int WarpX::do_compute_max_step_from_zmax = 0;
Real WarpX::zmax_plasma_to_compute_max_step = 0.;

long WarpX::current_deposition_algo;
long WarpX::charge_deposition_algo;
long WarpX::field_gathering_algo;
long WarpX::particle_pusher_algo;
int WarpX::maxwell_fdtd_solver_id;

long WarpX::n_rz_azimuthal_modes = 1;
long WarpX::ncomps = 1;

long WarpX::nox = 1;
long WarpX::noy = 1;
long WarpX::noz = 1;

bool WarpX::use_fdtd_nci_corr = false;
int  WarpX::l_lower_order_in_v = true;

bool WarpX::use_filter        = false;
bool WarpX::serialize_ics     = false;
bool WarpX::refine_plasma     = false;

int WarpX::num_mirrors = 0;

int  WarpX::sort_int = -1;

bool WarpX::do_boosted_frame_diagnostic = false;
std::string WarpX::lab_data_directory = "lab_frame_data";
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

IntVect WarpX::filter_npass_each_dir(1);

int WarpX::n_field_gather_buffer = 0;
int WarpX::n_current_deposition_buffer = -1;

int WarpX::do_nodal = false;

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
    warpx_do_continuous_injection = mypc->doContinuousInjection();
    if (warpx_do_continuous_injection){
        if (moving_window_v >= 0){
            // Inject particles continuously from the right end of the box
            current_injection_position = geom[0].ProbHi(moving_window_dir);
        } else {
            // Inject particles continuously from the left end of the box
            current_injection_position = geom[0].ProbLo(moving_window_dir);
        }
    }
    do_boosted_frame_particles = mypc->doBoostedFrameDiags();

    Efield_aux.resize(nlevs_max);
    Bfield_aux.resize(nlevs_max);

    F_fp.resize(nlevs_max);
    rho_fp.resize(nlevs_max);
    current_fp.resize(nlevs_max);
    Efield_fp.resize(nlevs_max);
    Bfield_fp.resize(nlevs_max);

    current_store.resize(nlevs_max);

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
    charge_buf.resize(nlevs_max);

    pml.resize(nlevs_max);

#ifdef WARPX_DO_ELECTROSTATIC
    masks.resize(nlevs_max);
    gather_masks.resize(nlevs_max);
#endif // WARPX_DO_ELECTROSTATIC

    costs.resize(nlevs_max);

#ifdef WARPX_USE_PSATD
    spectral_solver_fp.resize(nlevs_max);
    spectral_solver_cp.resize(nlevs_max);
#endif
#ifdef WARPX_USE_PSATD_HYBRID
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

#ifdef BL_USE_SENSEI_INSITU
    insitu_bridge = nullptr;
#endif

    // NCI Godfrey filters can have different stencils
    // at different levels (the stencil depends on c*dt/dz)
    nci_godfrey_filter_exeybz.resize(nlevs_max);
    nci_godfrey_filter_bxbyez.resize(nlevs_max);
}

WarpX::~WarpX ()
{
    int nlevs_max = maxLevel() +1;
    for (int lev = 0; lev < nlevs_max; ++lev) {
        ClearLevel(lev);
    }

#ifdef BL_USE_SENSEI_INSITU
    delete insitu_bridge;
#endif
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
        pp.query("do_subcycling", do_subcycling);
        pp.query("override_sync_int", override_sync_int);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(do_subcycling != 1 || max_level <= 1,
                                     "Subcycling method 1 only works for 2 levels.");

    ReadBoostedFrameParameters(gamma_boost, beta_boost, boost_direction);

    // pp.query returns 1 if argument zmax_plasma_to_compute_max_step is
    // specified by the user, 0 otherwise.
    do_compute_max_step_from_zmax =
        pp.query("zmax_plasma_to_compute_max_step",
                  zmax_plasma_to_compute_max_step);

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

    pp.query("do_boosted_frame_diagnostic", do_boosted_frame_diagnostic);
    if (do_boosted_frame_diagnostic) {

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(gamma_boost > 1.0,
               "gamma_boost must be > 1 to use the boosted frame diagnostic.");

        pp.query("lab_data_directory", lab_data_directory);

        std::string s;
        pp.get("boost_direction", s);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( (s == "z" || s == "Z"),
               "The boosted frame diagnostic currently only works if the boost is in the z direction.");

        pp.get("num_snapshots_lab", num_snapshots_lab);

        // Read either dz_snapshots_lab or dt_snapshots_lab
        bool snapshot_interval_is_specified = 0;
        Real dz_snapshots_lab = 0;
        snapshot_interval_is_specified += pp.query("dt_snapshots_lab", dt_snapshots_lab);
        if ( pp.query("dz_snapshots_lab", dz_snapshots_lab) ){
            dt_snapshots_lab = dz_snapshots_lab/PhysConst::c;
            snapshot_interval_is_specified = 1;
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            snapshot_interval_is_specified,
            "When using back-transformed diagnostics, user should specify either dz_snapshots_lab or dt_snapshots_lab.");

        pp.get("gamma_boost", gamma_boost);

        pp.query("do_boosted_frame_fields", do_boosted_frame_fields);

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(do_moving_window,
               "The moving window should be on if using the boosted frame diagnostic.");

        pp.get("moving_window_dir", s);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( (s == "z" || s == "Z"),
               "The boosted frame diagnostic currently only works if the moving window is in the z direction.");
    }

    pp.query("do_electrostatic", do_electrostatic);
    pp.query("n_buffer", n_buffer);
    pp.query("const_dt", const_dt);

    // Read filter and fill IntVect filter_npass_each_dir with
    // proper size for AMREX_SPACEDIM
	pp.query("use_filter", use_filter);
    Vector<int> parse_filter_npass_each_dir(AMREX_SPACEDIM,1);
    pp.queryarr("filter_npass_each_dir", parse_filter_npass_each_dir);
    filter_npass_each_dir[0] = parse_filter_npass_each_dir[0];
    filter_npass_each_dir[1] = parse_filter_npass_each_dir[1];
#if (AMREX_SPACEDIM == 3)
    filter_npass_each_dir[2] = parse_filter_npass_each_dir[2];
#endif

    pp.query("num_mirrors", num_mirrors);
    if (num_mirrors>0){
        mirror_z.resize(num_mirrors);
        pp.getarr("mirror_z", mirror_z, 0, num_mirrors);
        mirror_z_width.resize(num_mirrors);
        pp.getarr("mirror_z_width", mirror_z_width, 0, num_mirrors);
        mirror_z_npoints.resize(num_mirrors);
        pp.getarr("mirror_z_npoints", mirror_z_npoints, 0, num_mirrors);
    }

	pp.query("serialize_ics", serialize_ics);
	pp.query("refine_plasma", refine_plasma);
        pp.query("do_dive_cleaning", do_dive_cleaning);
        pp.query("n_field_gather_buffer", n_field_gather_buffer);
        pp.query("n_current_deposition_buffer", n_current_deposition_buffer);
	pp.query("sort_int", sort_int);

        pp.query("do_pml", do_pml);
        pp.query("pml_ncell", pml_ncell);
        pp.query("pml_delta", pml_delta);
        pp.query("pml_has_particles", pml_has_particles);
        pp.query("do_pml_j_damping", do_pml_j_damping);
        pp.query("do_pml_in_domain", do_pml_in_domain);

        Vector<int> parse_do_pml_Lo(AMREX_SPACEDIM,1);
        pp.queryarr("do_pml_Lo", parse_do_pml_Lo);
        do_pml_Lo[0] = parse_do_pml_Lo[0];
        do_pml_Lo[1] = parse_do_pml_Lo[1];
#if (AMREX_SPACEDIM == 3)
        do_pml_Lo[2] = parse_do_pml_Lo[2];
#endif
        Vector<int> parse_do_pml_Hi(AMREX_SPACEDIM,1);
        pp.queryarr("do_pml_Hi", parse_do_pml_Hi);
        do_pml_Hi[0] = parse_do_pml_Hi[0];
        do_pml_Hi[1] = parse_do_pml_Hi[1];
#if (AMREX_SPACEDIM == 3)
        do_pml_Hi[2] = parse_do_pml_Hi[2];
#endif

        if ( (do_pml_j_damping==1)&&(do_pml_in_domain==0) ){
          amrex::Abort("J-damping can only be done when PML are inside simulation domain (do_pml_in_domain=1)");
        }

        pp.query("dump_openpmd", dump_openpmd);
        pp.query("openpmd_backend", openpmd_backend);
        pp.query("dump_plotfiles", dump_plotfiles);
        pp.query("plot_raw_fields", plot_raw_fields);
        pp.query("plot_raw_fields_guards", plot_raw_fields_guards);
        pp.query("plot_coarsening_ratio", plot_coarsening_ratio);
        bool user_fields_to_plot;
        user_fields_to_plot = pp.queryarr("fields_to_plot", fields_to_plot);
        if (not user_fields_to_plot){
            // If not specified, set default values
            fields_to_plot = {"Ex", "Ey", "Ez", "Bx", "By",
                              "Bz", "jx", "jy", "jz",
                              "part_per_cell"};
        }
        // set plot_rho to true of the users requests it, so that
        // rho is computed at each iteration.
        if (std::find(fields_to_plot.begin(), fields_to_plot.end(), "rho")
            != fields_to_plot.end()){
            plot_rho = true;
        }
        // Sanity check if user requests to plot F
        if (std::find(fields_to_plot.begin(), fields_to_plot.end(), "F")
            != fields_to_plot.end()){
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(do_dive_cleaning,
                "plot F only works if warpx.do_dive_cleaning = 1");
        }
        // If user requests to plot proc_number for a serial run,
        // delete proc_number from fields_to_plot
        if (ParallelDescriptor::NProcs() == 1){
            fields_to_plot.erase(std::remove(fields_to_plot.begin(),
                                             fields_to_plot.end(),
                                             "proc_number"),
                                 fields_to_plot.end());
        }

        // Check that the coarsening_ratio can divide the blocking factor
        for (int lev=0; lev<maxLevel(); lev++){
          for (int comp=0; comp<AMREX_SPACEDIM; comp++){
            if ( blockingFactor(lev)[comp] % plot_coarsening_ratio != 0 ){
              amrex::Abort("plot_coarsening_ratio should be an integer "
                           "divisor of the blocking factor.");
            }
          }
        }

        pp.query("plot_finepatch", plot_finepatch);
        if (maxLevel() > 0) {
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

        pp.query("do_nodal", do_nodal);
        if (do_nodal) {
            Bx_nodal_flag = IntVect::TheNodeVector();
            By_nodal_flag = IntVect::TheNodeVector();
            Bz_nodal_flag = IntVect::TheNodeVector();
            Ex_nodal_flag = IntVect::TheNodeVector();
            Ey_nodal_flag = IntVect::TheNodeVector();
            Ez_nodal_flag = IntVect::TheNodeVector();
            jx_nodal_flag = IntVect::TheNodeVector();
            jy_nodal_flag = IntVect::TheNodeVector();
            jz_nodal_flag = IntVect::TheNodeVector();
            // Use same shape factors in all directions, for gathering
            l_lower_order_in_v = false;
        }

        // Only needs to be set with WARPX_DIM_RZ, otherwise defaults to 1.
        pp.query("n_rz_azimuthal_modes", n_rz_azimuthal_modes);

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
        current_deposition_algo = GetAlgorithmInteger(pp, "current_deposition");
        charge_deposition_algo = GetAlgorithmInteger(pp, "charge_deposition");
        field_gathering_algo = GetAlgorithmInteger(pp, "field_gathering");
        particle_pusher_algo = GetAlgorithmInteger(pp, "particle_pusher");
        maxwell_fdtd_solver_id = GetAlgorithmInteger(pp, "maxwell_fdtd_solver");
    }

#ifdef WARPX_USE_PSATD
    {
        ParmParse pp("psatd");
        pp.query("hybrid_mpi_decomposition", fft_hybrid_mpi_decomposition);
        pp.query("ngroups_fft", ngroups_fft);
        pp.query("fftw_plan_measure", fftw_plan_measure);
        pp.query("nox", nox_fft);
        pp.query("noy", noy_fft);
        pp.query("noz", noz_fft);
    }
#endif

    {
        insitu_start = 0;
        insitu_int = 0;
        insitu_config = "";
        insitu_pin_mesh = 0;

        ParmParse pp("insitu");
        pp.query("int", insitu_int);
        pp.query("start", insitu_start);
        pp.query("config", insitu_config);
        pp.query("pin_mesh", insitu_pin_mesh);
    }

    // for slice generation //
    {
       ParmParse pp("slice");
       amrex::Vector<Real> slice_lo(AMREX_SPACEDIM);
       amrex::Vector<Real> slice_hi(AMREX_SPACEDIM);
       Vector<int> slice_crse_ratio(AMREX_SPACEDIM);
       // set default slice_crse_ratio //
       for (int idim=0; idim < AMREX_SPACEDIM; ++idim )
       {
          slice_crse_ratio[idim] = 1;
       }
       pp.queryarr("dom_lo",slice_lo,0,AMREX_SPACEDIM);
       pp.queryarr("dom_hi",slice_hi,0,AMREX_SPACEDIM);
       pp.queryarr("coarsening_ratio",slice_crse_ratio,0,AMREX_SPACEDIM);
       pp.query("plot_int",slice_plot_int);
       slice_realbox.setLo(slice_lo);
       slice_realbox.setHi(slice_hi);
       slice_cr_ratio = IntVect(AMREX_D_DECL(1,1,1));
       for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
       {
          if (slice_crse_ratio[idim] > 1 ) {
             slice_cr_ratio[idim] = slice_crse_ratio[idim];
          }
       }

    }

}

// This is a virtual function.
void
WarpX::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& new_grids,
                                const DistributionMapping& new_dmap)
{
    AllocLevelData(lev, new_grids, new_dmap);
    InitLevelData(lev, time);

#ifdef WARPX_USE_PSATD
    if (fft_hybrid_mpi_decomposition){
#ifdef WARPX_USE_PSATD_HYBRID
        AllocLevelDataFFT(lev);
        InitLevelDataFFT(lev, time);
#else
    amrex::Abort("The option `psatd.fft_hybrid_mpi_decomposition` does not work on GPU.");
#endif
    }
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

        current_store[lev][i].reset();

	current_cp[lev][i].reset();
	Efield_cp [lev][i].reset();
	Bfield_cp [lev][i].reset();

	Efield_cax[lev][i].reset();
	Bfield_cax[lev][i].reset();
        current_buf[lev][i].reset();
    }

    charge_buf[lev].reset();

    current_buffer_masks[lev].reset();
    gather_buffer_masks[lev].reset();

    F_fp  [lev].reset();
    rho_fp[lev].reset();
    F_cp  [lev].reset();
    rho_cp[lev].reset();

    costs[lev].reset();

#ifdef WARPX_USE_PSATD_HYBRID
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
    // When using subcycling, the particles on the fine level perform two pushes
    // before being redistributed ; therefore, we need one extra guard cell
    // (the particles may move by 2*c*dt)
    const int ngx_tmp = (maxLevel() > 0 && do_subcycling == 1) ? WarpX::nox+1 : WarpX::nox;
    const int ngy_tmp = (maxLevel() > 0 && do_subcycling == 1) ? WarpX::noy+1 : WarpX::noy;
    const int ngz_tmp = (maxLevel() > 0 && do_subcycling == 1) ? WarpX::noz+1 : WarpX::noz;

    // Ex, Ey, Ez, Bx, By, and Bz have the same number of ghost cells.
    // jx, jy, jz and rho have the same number of ghost cells.
    // E and B have the same number of ghost cells as j and rho if NCI filter is not used,
    // but different number of ghost cells in z-direction if NCI filter is used.
    // The number of cells should be even, in order to easily perform the
    // interpolation from coarse grid to fine grid.
    int ngx = (ngx_tmp % 2) ? ngx_tmp+1 : ngx_tmp;  // Always even number
    int ngy = (ngy_tmp % 2) ? ngy_tmp+1 : ngy_tmp;  // Always even number
    int ngz_nonci = (ngz_tmp % 2) ? ngz_tmp+1 : ngz_tmp;  // Always even number
    int ngz;
    if (WarpX::use_fdtd_nci_corr) {
        int ng = ngz_tmp + NCIGodfreyFilter::stencil_width;
        ngz = (ng % 2) ? ng+1 : ng;
    } else {
        ngz = ngz_nonci;
    }

    // J is only interpolated from fine to coarse (not coarse to fine)
    // and therefore does not need to be even.
    int ngJx = ngx_tmp;
    int ngJy = ngy_tmp;
    int ngJz = ngz_tmp;

    // When calling the moving window (with one level of refinement),  we shift
    // the fine grid by 2 cells ; therefore, we need at least 2 guard cells
    // on level 1. This may not be necessary for level 0.
    if (do_moving_window) {
        ngx = std::max(ngx,2);
        ngy = std::max(ngy,2);
        ngz = std::max(ngz,2);
        ngJx = std::max(ngJx,2);
        ngJy = std::max(ngJy,2);
        ngJz = std::max(ngJz,2);
    }

#if (AMREX_SPACEDIM == 3)
    IntVect ngE(ngx,ngy,ngz);
    IntVect ngJ(ngJx,ngJy,ngJz);
#elif (AMREX_SPACEDIM == 2)
    IntVect ngE(ngx,ngz);
    IntVect ngJ(ngJx,ngJz);
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

#ifdef WARPX_USE_PSATD
    if (fft_hybrid_mpi_decomposition == false){
        // All boxes should have the same number of guard cells
        // (to avoid temporary parallel copies)
        // Thus take the max of the required number of guards for each field
        // Also: the number of guard cell should be enough to contain
        // the stencil of the FFT solver. Here, this number (`ngFFT`)
        // is determined *empirically* to be the order of the solver
        // for nodal, and half the order of the solver for staggered.
        IntVect ngFFT;
        if (do_nodal) {
            ngFFT = IntVect(AMREX_D_DECL(nox_fft, noy_fft, noz_fft));
        } else {
            ngFFT = IntVect(AMREX_D_DECL(nox_fft/2, noy_fft/2, noz_fft/2));
        }
        for (int i_dim=0; i_dim<AMREX_SPACEDIM; i_dim++ ){
            int ng_required = ngFFT[i_dim];
            // Get the max
            ng_required = std::max( ng_required, ngE[i_dim] );
            ng_required = std::max( ng_required, ngJ[i_dim] );
            ng_required = std::max( ng_required, ngRho[i_dim] );
            ng_required = std::max( ng_required, ngF );
            // Set the guard cells to this max
            ngE[i_dim] = ng_required;
            ngJ[i_dim] = ng_required;
            ngRho[i_dim] = ng_required;
            ngF = ng_required;
        }
    }
#endif

    AllocLevelMFs(lev, ba, dm, ngE, ngJ, ngRho, ngF);
}

void
WarpX::AllocLevelMFs (int lev, const BoxArray& ba, const DistributionMapping& dm,
                      const IntVect& ngE, const IntVect& ngJ, const IntVect& ngRho, int ngF)
{

#if defined WARPX_DIM_RZ
    // With RZ multimode, there is a real and imaginary component
    // for each mode, except mode 0 which is purely real
    // Component 0 is mode 0.
    // Odd components are the real parts.
    // Even components are the imaginary parts.
    ncomps = n_rz_azimuthal_modes*2 - 1;
#endif

    //
    // The fine patch
    //
    Bfield_fp[lev][0].reset( new MultiFab(amrex::convert(ba,Bx_nodal_flag),dm,ncomps,ngE));
    Bfield_fp[lev][1].reset( new MultiFab(amrex::convert(ba,By_nodal_flag),dm,ncomps,ngE));
    Bfield_fp[lev][2].reset( new MultiFab(amrex::convert(ba,Bz_nodal_flag),dm,ncomps,ngE));

    Efield_fp[lev][0].reset( new MultiFab(amrex::convert(ba,Ex_nodal_flag),dm,ncomps,ngE));
    Efield_fp[lev][1].reset( new MultiFab(amrex::convert(ba,Ey_nodal_flag),dm,ncomps,ngE));
    Efield_fp[lev][2].reset( new MultiFab(amrex::convert(ba,Ez_nodal_flag),dm,ncomps,ngE));

    current_fp[lev][0].reset( new MultiFab(amrex::convert(ba,jx_nodal_flag),dm,ncomps,ngJ));
    current_fp[lev][1].reset( new MultiFab(amrex::convert(ba,jy_nodal_flag),dm,ncomps,ngJ));
    current_fp[lev][2].reset( new MultiFab(amrex::convert(ba,jz_nodal_flag),dm,ncomps,ngJ));

    if (do_dive_cleaning || plot_rho)
    {
        rho_fp[lev].reset(new MultiFab(amrex::convert(ba,IntVect::TheUnitVector()),dm,2*ncomps,ngRho));
    }

    if (do_subcycling == 1 && lev == 0)
    {
        current_store[lev][0].reset( new MultiFab(amrex::convert(ba,jx_nodal_flag),dm,ncomps,ngJ));
        current_store[lev][1].reset( new MultiFab(amrex::convert(ba,jy_nodal_flag),dm,ncomps,ngJ));
        current_store[lev][2].reset( new MultiFab(amrex::convert(ba,jz_nodal_flag),dm,ncomps,ngJ));
    }

    if (do_dive_cleaning)
    {
        F_fp[lev].reset  (new MultiFab(amrex::convert(ba,IntVect::TheUnitVector()),dm,ncomps, ngF));
    }
#ifdef WARPX_USE_PSATD
    else
    {
        rho_fp[lev].reset(new MultiFab(amrex::convert(ba,IntVect::TheUnitVector()),dm,2*ncomps,ngRho));
    }
    if (fft_hybrid_mpi_decomposition == false){
        // Allocate and initialize the spectral solver
        std::array<Real,3> dx = CellSize(lev);
#if (AMREX_SPACEDIM == 3)
        RealVect dx_vect(dx[0], dx[1], dx[2]);
#elif (AMREX_SPACEDIM == 2)
        RealVect dx_vect(dx[0], dx[2]);
#endif
        // Get the cell-centered box, with guard cells
        BoxArray realspace_ba = ba;  // Copy box
        realspace_ba.enclosedCells().grow(ngE); // cell-centered + guard cells
        // Define spectral solver
        spectral_solver_fp[lev].reset( new SpectralSolver( realspace_ba, dm,
            nox_fft, noy_fft, noz_fft, do_nodal, dx_vect, dt[lev] ) );
    }
#endif

    //
    // The Aux patch (i.e., the full solution)
    //
    if (lev == 0)
    {
        for (int idir = 0; idir < 3; ++idir) {
            Efield_aux[lev][idir].reset(new MultiFab(*Efield_fp[lev][idir], amrex::make_alias, 0, ncomps));
            Bfield_aux[lev][idir].reset(new MultiFab(*Bfield_fp[lev][idir], amrex::make_alias, 0, ncomps));
        }
    }
    else
    {
        Bfield_aux[lev][0].reset( new MultiFab(amrex::convert(ba,Bx_nodal_flag),dm,ncomps,ngE));
        Bfield_aux[lev][1].reset( new MultiFab(amrex::convert(ba,By_nodal_flag),dm,ncomps,ngE));
        Bfield_aux[lev][2].reset( new MultiFab(amrex::convert(ba,Bz_nodal_flag),dm,ncomps,ngE));

        Efield_aux[lev][0].reset( new MultiFab(amrex::convert(ba,Ex_nodal_flag),dm,ncomps,ngE));
        Efield_aux[lev][1].reset( new MultiFab(amrex::convert(ba,Ey_nodal_flag),dm,ncomps,ngE));
        Efield_aux[lev][2].reset( new MultiFab(amrex::convert(ba,Ez_nodal_flag),dm,ncomps,ngE));
    }

    //
    // The coarse patch
    //
    if (lev > 0)
    {
        BoxArray cba = ba;
        cba.coarsen(refRatio(lev-1));

        // Create the MultiFabs for B
        Bfield_cp[lev][0].reset( new MultiFab(amrex::convert(cba,Bx_nodal_flag),dm,ncomps,ngE));
        Bfield_cp[lev][1].reset( new MultiFab(amrex::convert(cba,By_nodal_flag),dm,ncomps,ngE));
        Bfield_cp[lev][2].reset( new MultiFab(amrex::convert(cba,Bz_nodal_flag),dm,ncomps,ngE));

        // Create the MultiFabs for E
        Efield_cp[lev][0].reset( new MultiFab(amrex::convert(cba,Ex_nodal_flag),dm,ncomps,ngE));
        Efield_cp[lev][1].reset( new MultiFab(amrex::convert(cba,Ey_nodal_flag),dm,ncomps,ngE));
        Efield_cp[lev][2].reset( new MultiFab(amrex::convert(cba,Ez_nodal_flag),dm,ncomps,ngE));

        // Create the MultiFabs for the current
        current_cp[lev][0].reset( new MultiFab(amrex::convert(cba,jx_nodal_flag),dm,ncomps,ngJ));
        current_cp[lev][1].reset( new MultiFab(amrex::convert(cba,jy_nodal_flag),dm,ncomps,ngJ));
        current_cp[lev][2].reset( new MultiFab(amrex::convert(cba,jz_nodal_flag),dm,ncomps,ngJ));

        if (do_dive_cleaning || plot_rho){
            rho_cp[lev].reset(new MultiFab(amrex::convert(cba,IntVect::TheUnitVector()),dm,2*ncomps,ngRho));
        }
        if (do_dive_cleaning)
        {
            F_cp[lev].reset  (new MultiFab(amrex::convert(cba,IntVect::TheUnitVector()),dm,ncomps, ngF));
        }
#ifdef WARPX_USE_PSATD
        else
        {
            rho_cp[lev].reset(new MultiFab(amrex::convert(cba,IntVect::TheUnitVector()),dm,2*ncomps,ngRho));
        }
        if (fft_hybrid_mpi_decomposition == false){
            // Allocate and initialize the spectral solver
            std::array<Real,3> cdx = CellSize(lev-1);
    #if (AMREX_SPACEDIM == 3)
            RealVect cdx_vect(cdx[0], cdx[1], cdx[2]);
    #elif (AMREX_SPACEDIM == 2)
            RealVect cdx_vect(cdx[0], cdx[2]);
    #endif
            // Get the cell-centered box, with guard cells
            BoxArray realspace_ba = cba;  // Copy box
            realspace_ba.enclosedCells().grow(ngE); // cell-centered + guard cells
            // Define spectral solver
            spectral_solver_cp[lev].reset( new SpectralSolver( realspace_ba, dm,
                nox_fft, noy_fft, noz_fft, do_nodal, cdx_vect, dt[lev] ) );
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
            Bfield_cax[lev][0].reset( new MultiFab(amrex::convert(cba,Bx_nodal_flag),dm,ncomps,ngE));
            Bfield_cax[lev][1].reset( new MultiFab(amrex::convert(cba,By_nodal_flag),dm,ncomps,ngE));
            Bfield_cax[lev][2].reset( new MultiFab(amrex::convert(cba,Bz_nodal_flag),dm,ncomps,ngE));

            // Create the MultiFabs for E
            Efield_cax[lev][0].reset( new MultiFab(amrex::convert(cba,Ex_nodal_flag),dm,ncomps,ngE));
            Efield_cax[lev][1].reset( new MultiFab(amrex::convert(cba,Ey_nodal_flag),dm,ncomps,ngE));
            Efield_cax[lev][2].reset( new MultiFab(amrex::convert(cba,Ez_nodal_flag),dm,ncomps,ngE));

            gather_buffer_masks[lev].reset( new iMultiFab(ba, dm, ncomps, 1) );
            // Gather buffer masks have 1 ghost cell, because of the fact
            // that particles may move by more than one cell when using subcycling.
        }

        if (n_current_deposition_buffer > 0) {
            current_buf[lev][0].reset( new MultiFab(amrex::convert(cba,jx_nodal_flag),dm,ncomps,ngJ));
            current_buf[lev][1].reset( new MultiFab(amrex::convert(cba,jy_nodal_flag),dm,ncomps,ngJ));
            current_buf[lev][2].reset( new MultiFab(amrex::convert(cba,jz_nodal_flag),dm,ncomps,ngJ));
            if (rho_cp[lev]) {
                charge_buf[lev].reset( new MultiFab(amrex::convert(cba,IntVect::TheUnitVector()),dm,2*ncomps,ngRho));
            }
            current_buffer_masks[lev].reset( new iMultiFab(ba, dm, ncomps, 1) );
            // Current buffer masks have 1 ghost cell, because of the fact
            // that particles may move by more than one cell when using subcycling.
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
WarpX::Evolve (int numsteps) {
    BL_PROFILE_REGION("WarpX::Evolve()");

#ifdef WARPX_DO_ELECTROSTATIC
    if (do_electrostatic) {
        EvolveES(numsteps);
    } else {
      EvolveEM(numsteps);
    }
#else
    EvolveEM(numsteps);
#endif // WARPX_DO_ELECTROSTATIC
}

void
WarpX::ComputeDivB (MultiFab& divB, int dcomp,
                    const std::array<const MultiFab*, 3>& B,
                    const std::array<Real,3>& dx)
{
    Real dxinv = 1./dx[0], dyinv = 1./dx[1], dzinv = 1./dx[2];

#ifdef WARPX_DIM_RZ
    const Real rmin = GetInstance().Geom(0).ProbLo(0);
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(divB, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& Bxfab = B[0]->array(mfi);
        auto const& Byfab = B[1]->array(mfi);
        auto const& Bzfab = B[2]->array(mfi);
        auto const& divBfab = divB.array(mfi);

        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            warpx_computedivb(i, j, k, dcomp, divBfab, Bxfab, Byfab, Bzfab, dxinv, dyinv, dzinv
#ifdef WARPX_DIM_RZ
                              ,rmin
#endif
                              );
        });
    }
}

void
WarpX::ComputeDivB (MultiFab& divB, int dcomp,
                    const std::array<const MultiFab*, 3>& B,
                    const std::array<Real,3>& dx, int ngrow)
{
    Real dxinv = 1./dx[0], dyinv = 1./dx[1], dzinv = 1./dx[2];

#ifdef WARPX_DIM_RZ
    const Real rmin = GetInstance().Geom(0).ProbLo(0);
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(divB, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox(ngrow);
        auto const& Bxfab = B[0]->array(mfi);
        auto const& Byfab = B[1]->array(mfi);
        auto const& Bzfab = B[2]->array(mfi);
        auto const& divBfab = divB.array(mfi);

        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            warpx_computedivb(i, j, k, dcomp, divBfab, Bxfab, Byfab, Bzfab, dxinv, dyinv, dzinv
#ifdef WARPX_DIM_RZ
                              ,rmin
#endif
                              );
        });
    }
}

void
WarpX::ComputeDivE (MultiFab& divE, int dcomp,
                    const std::array<const MultiFab*, 3>& E,
                    const std::array<Real,3>& dx)
{
    Real dxinv = 1./dx[0], dyinv = 1./dx[1], dzinv = 1./dx[2];

#ifdef WARPX_DIM_RZ
    const Real rmin = GetInstance().Geom(0).ProbLo(0);
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(divE, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& Exfab = E[0]->array(mfi);
        auto const& Eyfab = E[1]->array(mfi);
        auto const& Ezfab = E[2]->array(mfi);
        auto const& divEfab = divE.array(mfi);

        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            warpx_computedive(i, j, k, dcomp, divEfab, Exfab, Eyfab, Ezfab, dxinv, dyinv, dzinv
#ifdef WARPX_DIM_RZ
                              ,rmin
#endif
                              );
        });
    }
}

void
WarpX::ComputeDivE (MultiFab& divE, int dcomp,
                    const std::array<const MultiFab*, 3>& E,
                    const std::array<Real,3>& dx, int ngrow)
{
    Real dxinv = 1./dx[0], dyinv = 1./dx[1], dzinv = 1./dx[2];

#ifdef WARPX_DIM_RZ
    const Real rmin = GetInstance().Geom(0).ProbLo(0);
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(divE, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox(ngrow);
        auto const& Exfab = E[0]->array(mfi);
        auto const& Eyfab = E[1]->array(mfi);
        auto const& Ezfab = E[2]->array(mfi);
        auto const& divEfab = divE.array(mfi);

        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            warpx_computedive(i, j, k, dcomp, divEfab, Exfab, Eyfab, Ezfab, dxinv, dyinv, dzinv
#ifdef WARPX_DIM_RZ
                              ,rmin
#endif
                              );
        });
    }
}

void
WarpX::BuildBufferMasks ()
{
    for (int lev = 1; lev <= maxLevel(); ++lev)
    {
        for (int ipass = 0; ipass < 2; ++ipass)
        {
            int ngbuffer = (ipass == 0) ? n_current_deposition_buffer : n_field_gather_buffer;
            iMultiFab* bmasks = (ipass == 0) ? current_buffer_masks[lev].get() : gather_buffer_masks[lev].get();
            if (bmasks)
            {
                const int ngtmp = ngbuffer + bmasks->nGrow();
                iMultiFab tmp(bmasks->boxArray(), bmasks->DistributionMap(), 1, ngtmp);
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
                    const Box& tbx = mfi.growntilebox();
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

void
WarpX::StoreCurrent (int lev)
{
    for (int idim = 0; idim < 3; ++idim) {
        if (current_store[lev][idim]) {
            MultiFab::Copy(*current_store[lev][idim], *current_fp[lev][idim],
                           0, 0, 1, current_store[lev][idim]->nGrowVect());
        }
    }
}

void
WarpX::RestoreCurrent (int lev)
{
    for (int idim = 0; idim < 3; ++idim) {
        if (current_store[lev][idim]) {
            std::swap(current_fp[lev][idim], current_store[lev][idim]);
        }
    }
}

std::string
WarpX::Version ()
{
#ifdef WARPX_GIT_VERSION
    return std::string(WARPX_GIT_VERSION);
#else
    return std::string("Unknown");
#endif
}

std::string
WarpX::PicsarVersion ()
{
#ifdef WARPX_GIT_VERSION
    return std::string(PICSAR_GIT_VERSION);
#else
    return std::string("Unknown");
#endif
}
