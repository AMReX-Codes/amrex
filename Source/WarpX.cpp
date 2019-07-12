
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
#include <WarpXAlgorithmSelection.H>
#include <WarpX_FDTD.H>

#ifdef BL_USE_SENSEI_INSITU
#include <AMReX_AmrMeshInSituBridge.H>
#endif

using namespace amrex;

Vector<Real> WarpX::B_external(3, 0.0);

int WarpX::do_moving_window = 0;
int WarpX::moving_window_dir = -1;

Real WarpX::gamma_boost = 1.;
Real WarpX::beta_boost = 0.;
Vector<int> WarpX::boost_direction = {0,0,0};
int WarpX::do_compute_max_step_from_zmax = 0;
Real WarpX::zmax_plasma_to_compute_max_step = 0.;

long WarpX::use_picsar_deposition = 1;
long WarpX::current_deposition_algo;
long WarpX::charge_deposition_algo;
long WarpX::field_gathering_algo;
long WarpX::particle_pusher_algo;
int WarpX::maxwell_fdtd_solver_id;

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

    current_fp_owner_masks.resize(nlevs_max);
    current_cp_owner_masks.resize(nlevs_max);
    rho_fp_owner_masks.resize(nlevs_max);
    rho_cp_owner_masks.resize(nlevs_max);

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

    spectral_solver_fp.resize(nlevs_max);
    spectral_solver_cp.resize(nlevs_max);

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
        pp.get("dt_snapshots_lab", dt_snapshots_lab);
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

        pp.query("dump_openpmd", dump_openpmd);
        pp.query("dump_plotfiles", dump_plotfiles);
        pp.query("plot_raw_fields", plot_raw_fields);
        pp.query("plot_raw_fields_guards", plot_raw_fields_guards);
        if (ParallelDescriptor::NProcs() == 1) {
            plot_proc_number = false;
        }
        // Fields to dump into plotfiles
        pp.query("plot_E_field"      , plot_E_field);
        pp.query("plot_B_field"      , plot_B_field);
        pp.query("plot_J_field"      , plot_J_field);
        pp.query("plot_part_per_cell", plot_part_per_cell);
        pp.query("plot_part_per_grid", plot_part_per_grid);
        pp.query("plot_part_per_proc", plot_part_per_proc);
        pp.query("plot_proc_number"  , plot_proc_number);
        pp.query("plot_dive"         , plot_dive);
        pp.query("plot_divb"         , plot_divb);
        pp.query("plot_rho"          , plot_rho);
        pp.query("plot_F"            , plot_F);
        pp.query("plot_coarsening_ratio", plot_coarsening_ratio);

        // Check that the coarsening_ratio can divide the blocking factor
        for (int lev=0; lev<maxLevel(); lev++){
          for (int comp=0; comp<AMREX_SPACEDIM; comp++){
            if ( blockingFactor(lev)[comp] % plot_coarsening_ratio != 0 ){
              amrex::Abort("plot_coarsening_ratio should be an integer divisor of the blocking factor.");
            }
          }
        }

        if (plot_F){
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(do_dive_cleaning,
                "plot_F only works if warpx.do_dive_cleaning = 1");
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
        // If not in RZ mode, read use_picsar_deposition
        // In RZ mode, use_picsar_deposition is on, as the C++ version 
        // of the deposition does not support RZ
#ifndef WARPX_RZ
        pp.query("use_picsar_deposition", use_picsar_deposition);
#endif
        current_deposition_algo = GetAlgorithmInteger(pp, "current_deposition");
        if (!use_picsar_deposition){
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE( current_deposition_algo >= 2, 
                "if not use_picsar_deposition, cannot use Esirkepov deposition.");
        }
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
        // Override value
        if (fft_hybrid_mpi_decomposition==false) ngroups_fft=ParallelDescriptor::NProcs();
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

        current_store[lev][i].reset();

	current_cp[lev][i].reset();
	Efield_cp [lev][i].reset();
	Bfield_cp [lev][i].reset();

	Efield_cax[lev][i].reset();
	Bfield_cax[lev][i].reset();
        current_buf[lev][i].reset();

        current_fp_owner_masks[lev][i].reset();
        current_cp_owner_masks[lev][i].reset();
    }

    rho_fp_owner_masks[lev].reset();
    rho_cp_owner_masks[lev].reset();

    charge_buf[lev].reset();

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

    AllocLevelMFs(lev, ba, dm, ngE, ngJ, ngRho, ngF);
}

void
WarpX::AllocLevelMFs (int lev, const BoxArray& ba, const DistributionMapping& dm,
                      const IntVect& ngE, const IntVect& ngJ, const IntVect& ngRho, int ngF)
{
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

    const auto& period = Geom(lev).periodicity();
    current_fp_owner_masks[lev][0] = std::move(current_fp[lev][0]->OwnerMask(period));
    current_fp_owner_masks[lev][1] = std::move(current_fp[lev][1]->OwnerMask(period));
    current_fp_owner_masks[lev][2] = std::move(current_fp[lev][2]->OwnerMask(period));

    if (do_dive_cleaning || plot_rho)
    {
        rho_fp[lev].reset(new MultiFab(amrex::convert(ba,IntVect::TheUnitVector()),dm,2,ngRho));
        rho_fp_owner_masks[lev] = std::move(rho_fp[lev]->OwnerMask(period));
    }

    if (do_subcycling == 1 && lev == 0)
    {
        current_store[lev][0].reset( new MultiFab(amrex::convert(ba,jx_nodal_flag),dm,1,ngJ));
        current_store[lev][1].reset( new MultiFab(amrex::convert(ba,jy_nodal_flag),dm,1,ngJ));
        current_store[lev][2].reset( new MultiFab(amrex::convert(ba,jz_nodal_flag),dm,1,ngJ));
    }

    if (do_dive_cleaning)
    {
        F_fp[lev].reset  (new MultiFab(amrex::convert(ba,IntVect::TheUnitVector()),dm,1, ngF));
    }
#ifdef WARPX_USE_PSATD
    else
    {
        rho_fp[lev].reset(new MultiFab(amrex::convert(ba,IntVect::TheUnitVector()),dm,2,ngRho));
        rho_fp_owner_masks[lev] = std::move(rho_fp[lev]->OwnerMask(period));
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

        const auto& cperiod = Geom(lev-1).periodicity();
        current_cp_owner_masks[lev][0] = std::move(current_cp[lev][0]->OwnerMask(cperiod));
        current_cp_owner_masks[lev][1] = std::move(current_cp[lev][1]->OwnerMask(cperiod));
        current_cp_owner_masks[lev][2] = std::move(current_cp[lev][2]->OwnerMask(cperiod));

        if (do_dive_cleaning || plot_rho){
            rho_cp[lev].reset(new MultiFab(amrex::convert(cba,IntVect::TheUnitVector()),dm,2,ngRho));
            rho_cp_owner_masks[lev] = std::move(rho_cp[lev]->OwnerMask(cperiod));
        }
        if (do_dive_cleaning)
        {
            F_cp[lev].reset  (new MultiFab(amrex::convert(cba,IntVect::TheUnitVector()),dm,1, ngF));
        }
#ifdef WARPX_USE_PSATD
        else
        {
            rho_cp[lev].reset(new MultiFab(amrex::convert(cba,IntVect::TheUnitVector()),dm,2,ngRho));
            rho_cp_owner_masks[lev] = std::move(rho_cp[lev]->OwnerMask(cperiod));
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

            gather_buffer_masks[lev].reset( new iMultiFab(ba, dm, 1, 1) );
            // Gather buffer masks have 1 ghost cell, because of the fact
            // that particles may move by more than one cell when using subcycling.
        }

        if (n_current_deposition_buffer > 0) {
            current_buf[lev][0].reset( new MultiFab(amrex::convert(cba,jx_nodal_flag),dm,1,ngJ));
            current_buf[lev][1].reset( new MultiFab(amrex::convert(cba,jy_nodal_flag),dm,1,ngJ));
            current_buf[lev][2].reset( new MultiFab(amrex::convert(cba,jz_nodal_flag),dm,1,ngJ));
            if (do_dive_cleaning || plot_rho) {
                charge_buf[lev].reset( new MultiFab(amrex::convert(cba,IntVect::TheUnitVector()),dm,2,ngRho));
            }
            current_buffer_masks[lev].reset( new iMultiFab(ba, dm, 1, 1) );
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

#ifdef WARPX_RZ
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
#ifdef WARPX_RZ
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

#ifdef WARPX_RZ
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
#ifdef WARPX_RZ
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

#ifdef WARPX_RZ
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
#ifdef WARPX_RZ
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

#ifdef WARPX_RZ
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
#ifdef WARPX_RZ
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
