#include <iomanip>
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include <unistd.h>

using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::pair;
using std::string;

#include <AMReX_CONSTANTS.H>
#include <Nyx.H>
#include <Nyx_F.H>
#include <Derive_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_Particles_F.H>
#include <AMReX_Utility.H>
#include <AMReX_Print.H>

#if BL_USE_MPI
#include "MemInfo.H"
#endif

#ifdef GRAVITY
#include "Gravity.H"
#endif

#ifdef FORCING
#include "Forcing.H"
#endif

#ifdef REEBER
#include <ReeberAnalysis.H>
#endif

#ifdef GIMLET
#include <DoGimletAnalysis.H>
#include <postprocess_tau_fields.H>
#endif

#ifdef AGN
#include "agn_F.H"
#endif

using namespace amrex;
using namespace perilla;

extern "C" {
  int get_comp_urho();
  int get_comp_temp();
  int get_comp_e_int();
}

const int NyxHaloFinderSignal = 42;
const int GimletSignal = 55;

static int sum_interval = -1;
static Real fixed_dt    = -1.0;
static Real initial_dt  = -1.0;
static Real dt_cutoff   =  0;

int simd_width = 1;

int Nyx::strict_subcycling = 0;

Real Nyx::old_a      = -1.0;
Real Nyx::new_a      = -1.0;
Real Nyx::old_a_time = -1.0;
Real Nyx::new_a_time = -1.0;

Vector<Real> Nyx::plot_z_values;
Vector<Real> Nyx::analysis_z_values;

bool Nyx::dump_old = false;
int Nyx::verbose      = 0;

Real Nyx::cfl = 0.8;
Real Nyx::init_shrink = 1.0;
Real Nyx::change_max  = 1.1;

BCRec Nyx::phys_bc;
int Nyx::do_reflux = 1;
int Nyx::NUM_STATE = -1;
int Nyx::NUM_GROW  = -1;

int Nyx::nsteps_from_plotfile = -1;

ErrorList Nyx::err_list;

int Nyx::Density = -1;
int Nyx::Eden = -1;
int Nyx::Eint = -1;
int Nyx::Xmom = -1;
int Nyx::Ymom = -1;
int Nyx::Zmom = -1;

int Nyx::Temp_comp = -1;
int Nyx::  Ne_comp = -1;
int Nyx:: Zhi_comp = -1;

int Nyx::NumSpec  = 0;
int Nyx::NumAux   = 0;
int Nyx::NumAdv   = 0;

int Nyx::FirstSpec = -1;
int Nyx::FirstAux  = -1;
int Nyx::FirstAdv  = -1;

Real Nyx::small_dens = -1.e200;
Real Nyx::small_temp = -1.e200;
Real Nyx::gamma      =  0;

Real Nyx::comoving_OmB;
Real Nyx::comoving_OmM;
Real Nyx::comoving_h;

int Nyx::do_hydro = -1;
int Nyx::add_ext_src = 0;
int Nyx::heat_cool_type = 0;
int Nyx::strang_split = 1;
#ifdef SDC
int Nyx::sdc_split    = 0;
#endif

Real Nyx::average_gas_density = 0;
Real Nyx::average_dm_density = 0;
Real Nyx::average_neutr_density = 0;
Real Nyx::average_total_density = 0;

int         Nyx::inhomo_reion = 0;
std::string Nyx::inhomo_zhi_file = "";
int         Nyx::inhomo_grid = -1;

static int  slice_int    = -1;
std::string slice_file   = "slice_";
static int  slice_nfiles = 128;

// Real Nyx::ave_lev_vorticity[10];
// Real Nyx::std_lev_vorticity[10];

#ifdef GRAVITY
Gravity* Nyx::gravity  =  0;
int Nyx::do_grav       = -1;
#else
int Nyx::do_grav       =  0;
#endif

#ifdef FORCING
StochasticForcing* Nyx::forcing = 0;
int Nyx::do_forcing = -1;
#else
int Nyx::do_forcing =  0;
#endif

int Nyx::allow_untagging    = 0;
int Nyx::use_const_species  = 0;
int Nyx::normalize_species  = 0;
int Nyx::do_special_tagging = 0;
int Nyx::ppm_type           = 1;
int Nyx::ppm_reference      = 1;
int Nyx::corner_coupling    = 1;

int Nyx::use_colglaz        = 0;
int Nyx::version_2          = 0;

int Nyx::use_flattening     = 1;
int Nyx::ppm_flatten_before_integrals = 0;

Real Nyx:: h_species        = 0.0;
Real Nyx::he_species        = 0.0;

int Nyx::use_exact_gravity  = 0;

#ifdef AGN
Real Nyx::mass_halo_min     = 1.e10;
Real Nyx::mass_seed         = 1.e5;
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

int Nyx::write_parameters_in_plotfile = true;
int Nyx::print_fortran_warnings       = true;

// Do we use separate SPH particles to initialize
//  the density and momentum on the grid?
int  Nyx::init_with_sph_particles = 0;

// Do we write the particles in single (IEEE32)
//  or doublue (NATIVE) precision?
#ifdef BL_SINGLE_PRECISION_PARTICLES
std::string Nyx::particle_plotfile_format = "IEEE32";
#else
std::string Nyx::particle_plotfile_format = "NATIVE";
#endif

// this will be reset upon restart
Real         Nyx::previousCPUTimeUsed = 0.0;

Real         Nyx::startCPUTime = 0.0;

int reeber_int(0);
int gimlet_int(0);

#ifdef USE_PERILLA
void Nyx::initPerilla(Real time)
{
  MultiFab& S_new = get_new_data(State_Type);
  RG_S = new RegionGraph(S_new.IndexArray().size());
  RG_S->buildTileArray(S_new);

  ext_src_old= new MultiFab(grids, dmap, NUM_STATE, NUM_GROW);
  ext_src_old->setVal(0.);
  //if (add_ext_src)
      //get_old_source(prev_time, dt, *ext_src_old);

  MultiFab&  S        = get_new_data(State_Type);
  MultiFab&  D        = get_new_data(DiagEOS_Type);
  S_old_tmp= new MultiFab(S.boxArray(), S.DistributionMap(), NUM_STATE, NUM_GROW);
  D_old_tmp= new MultiFab(D.boxArray(), D.DistributionMap(), NUM_STATE, NUM_GROW);

  grav_vector= new MultiFab(grids, dmap, BL_SPACEDIM, NUM_GROW);
  grav_vector->setVal(0.);

  hydro_src= new MultiFab(grids, dmap, NUM_STATE, 0);
  hydro_src->setVal(0.);

  divu_cc= new MultiFab(grids, dmap, 1, 0);
  divu_cc->setVal(0.);

  //RG_Sborder = new RegionGraph(S_border.IndexArray().size());
  //RG_Sborder->buildTileArray(S_border);
}

void Nyx::finalizePerilla (Real time)
{
  if(RG_S) delete RG_S;
  //if(RG_Sborder) delete RG_Sborder;
  if(ext_src_old) delete ext_src_old;
//  if(S_old_tmp) delete S_old_tmp;
//  if(D_old_tmp) delete D_old_tmp;
  if(grav_vector) delete grav_vector;
  if(hydro_src) delete hydro_src;
  if(divu_cc) delete divu_cc;
}
#endif

// Note: Nyx::variableSetUp is in Nyx_setup.cpp
void
Nyx::variable_cleanup ()
{
#ifdef GRAVITY
    if (gravity != 0)
    {
        if (verbose > 1 && ParallelDescriptor::IOProcessor())
            std::cout << "Deleting gravity in variable_cleanup...\n";
        delete gravity;
        gravity = 0;
    }
#endif
#ifdef FORCING
    if (forcing != 0)
    {
        if (verbose > 1 && ParallelDescriptor::IOProcessor())
            std::cout << "Deleting forcing in variable_cleanup...\n";
        delete forcing;
        forcing = 0;
    }
#endif

    desc_lst.clear();
}

void
Nyx::read_params ()
{
    BL_PROFILE("Nyx::read_params()");
    static bool done = false;

    if (done) return;  // (caseywstark) when would this happen?

    done = true;  // ?

    ParmParse pp_nyx("nyx");

    pp_nyx.query("v", verbose);
    pp_nyx.get("init_shrink", init_shrink);
    pp_nyx.get("cfl", cfl);
    pp_nyx.query("change_max", change_max);
    pp_nyx.query("fixed_dt", fixed_dt);
    pp_nyx.query("initial_dt", initial_dt);
    pp_nyx.query("sum_interval", sum_interval);
    pp_nyx.query("do_reflux", do_reflux);
    do_reflux = (do_reflux ? 1 : 0);
    pp_nyx.get("dt_cutoff", dt_cutoff);

    pp_nyx.query("dump_old", dump_old);

    pp_nyx.query("small_dens", small_dens);
    pp_nyx.query("small_temp", small_temp);
    pp_nyx.query("gamma", gamma);

    pp_nyx.query("strict_subcycling",strict_subcycling);

#ifdef AMREX_USE_CVODE
    pp_nyx.query("simd_width", simd_width);
    if (simd_width < 1) amrex::Abort("simd_width must be a positive integer");
    set_simd_width(simd_width);

    if (verbose > 1) amrex::Print()
        << "SIMD width (# zones) for heating/cooling integration: "
        << simd_width << std::endl;
#endif

    // Get boundary conditions
    Vector<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    pp_nyx.getarr("lo_bc", lo_bc, 0, BL_SPACEDIM);
    pp_nyx.getarr("hi_bc", hi_bc, 0, BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        phys_bc.setLo(i, lo_bc[i]);
        phys_bc.setHi(i, hi_bc[i]);
    }

    //
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (Geometry::isAnyPeriodic())
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (Geometry::isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "Nyx::read_params:periodic in direction "
                              << dir
                              << " but low BC is not Interior" << std::endl;
                    amrex::Error();
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "Nyx::read_params:periodic in direction "
                              << dir
                              << " but high BC is not Interior" << std::endl;
                    amrex::Error();
                }
            }
        }
    }
    else
    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (lo_bc[dir] == Interior)
            {
                std::cerr << "Nyx::read_params:interior bc in direction "
                          << dir
                          << " but not periodic" << std::endl;
                amrex::Error();
            }
            if (hi_bc[dir] == Interior)
            {
                std::cerr << "Nyx::read_params:interior bc in direction "
                          << dir
                          << " but not periodic" << std::endl;
                amrex::Error();
            }
        }
    }

    pp_nyx.get("comoving_OmB", comoving_OmB);
    pp_nyx.get("comoving_OmM", comoving_OmM);
    pp_nyx.get("comoving_h", comoving_h);

    fort_set_omb(comoving_OmB);
    fort_set_omm(comoving_OmM);
    fort_set_hubble(comoving_h);

    pp_nyx.get("do_hydro", do_hydro);
#ifdef NO_HYDRO
    if (do_hydro == 1)
        amrex::Error("Cant have do_hydro == 1 when NO_HYDRO is true");
#endif

#ifdef NO_HYDRO
#ifndef GRAVITY
        amrex::Error("Dont know what to do with both hydro and gravity off");
#endif
#endif

    pp_nyx.query("add_ext_src", add_ext_src);
    pp_nyx.query("strang_split", strang_split);
#ifdef SDC
    pp_nyx.query("sdc_split", sdc_split);
    if (sdc_split == 1 && strang_split == 1)
        amrex::Error("Cant have strang_split == 1 and sdc_split == 1");
    if (sdc_split == 0 && strang_split == 0)
        amrex::Error("Cant have strang_split == 0 and sdc_split == 0");
    if (sdc_split != 1 && strang_split != 1)
        amrex::Error("Cant have strang_split != 1 and sdc_split != 1");
#else
    if (strang_split != 1)
        amrex::Error("Cant have strang_split != 1 with USE_SDC != TRUE");
#endif

#ifdef FORCING
    pp_nyx.get("do_forcing", do_forcing);
#ifdef NO_HYDRO
    if (do_forcing == 1)
        amrex::Error("Cant have do_forcing == 1 when NO_HYDRO is true ");
#endif
    if (do_forcing == 1 && add_ext_src == 0)
       amrex::Error("Nyx::must set add_ext_src to 1 if do_forcing = 1 ");
#else
    if (do_forcing == 1)
       amrex::Error("Nyx::you set do_forcing = 1 but forgot to set USE_FORCING = TRUE ");
#endif

    pp_nyx.query("heat_cool_type", heat_cool_type);
    if (heat_cool_type == 7)
    {
      amrex::Print() << "----- WARNING WARNING WARNING WARNING WARNING -----" << std::endl;
      amrex::Print() << "                                                   " << std::endl;
      amrex::Print() << "      SIMD CVODE is currently EXPERIMENTAL.        " << std::endl;
      amrex::Print() << "      Use at your own risk.                        " << std::endl;
      amrex::Print() << "                                                   " << std::endl;
      amrex::Print() << "----- WARNING WARNING WARNING WARNING WARNING -----" << std::endl;
      Vector<int> n_cell(BL_SPACEDIM);

      ParmParse pp_amr("amr");
      pp_amr.getarr("n_cell", n_cell, 0, BL_SPACEDIM);
      if (n_cell[0] % simd_width) {
        const std::string errmsg = "Currently the SIMD CVODE solver requires that n_cell[0] % simd_width = 0";
        amrex::Abort(errmsg);
      }
    }

    pp_nyx.query("use_exact_gravity", use_exact_gravity);

    pp_nyx.query("inhomo_reion", inhomo_reion);

    if (inhomo_reion) {
        pp_nyx.get("inhomo_zhi_file", inhomo_zhi_file);
        pp_nyx.get("inhomo_grid", inhomo_grid);
    }

#ifdef HEATCOOL
    if (heat_cool_type != 3 && heat_cool_type != 5 && heat_cool_type != 7)
       amrex::Error("Nyx:: nonzero heat_cool_type must equal 3 or 5 or 7");
    if (heat_cool_type == 0)
       amrex::Error("Nyx::contradiction -- HEATCOOL is defined but heat_cool_type == 0");

    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Integrating heating/cooling method with the following method: ";
      switch (heat_cool_type) {
        case 3:
          std::cout << "VODE";
          break;
        case 5:
          std::cout << "CVODE";
          break;
        case 7:
          std::cout << "SIMD CVODE";
          break;
      }
      std::cout << std::endl;
    }

#ifndef AMREX_USE_CVODE
    #ifndef AMREX_USE_SUNDIALS3
    if (heat_cool_type == 5 || heat_cool_type == 7)
        amrex::Error("Nyx:: cannot set heat_cool_type = 5 or 7 unless USE_CVODE=TRUE or USE_SUNDIALS3=TRUE");
    #endif
#else
    #ifdef SDC
    if (heat_cool_type == 7 && sdc_split == 1)
        amrex::Error("Nyx:: cannot set heat_cool_type = 7 with sdc_split = 1");
    #endif
#endif

#else
    if (heat_cool_type > 0)
       amrex::Error("Nyx::you set heat_cool_type > 0 but forgot to set USE_HEATCOOL = TRUE");
    if (inhomo_reion > 0)
       amrex::Error("Nyx::you set inhomo_reion > 0 but forgot to set USE_HEATCOOL = TRUE");
#endif

    pp_nyx.query("allow_untagging", allow_untagging);
    pp_nyx.query("use_const_species", use_const_species);
    pp_nyx.query("normalize_species", normalize_species);
    pp_nyx.query("ppm_type", ppm_type);
    pp_nyx.query("ppm_reference", ppm_reference);
    pp_nyx.query("ppm_flatten_before_integrals", ppm_flatten_before_integrals);
    pp_nyx.query("use_flattening", use_flattening);
    pp_nyx.query("use_colglaz", use_colglaz);
    pp_nyx.query("version_2", version_2);
    pp_nyx.query("corner_coupling", corner_coupling);

    if (do_hydro == 1)
    {
        if (do_hydro == 1 && use_const_species == 1)
        {
           pp_nyx.get("h_species" ,  h_species);
           pp_nyx.get("he_species", he_species);
           fort_set_xhydrogen(h_species);
           if (ParallelDescriptor::IOProcessor())
           {
               std::cout << "Nyx::setting species concentrations to "
                         << h_species << " and " << he_species
                         << " in the hydro and in the EOS " << std::endl;
           }
        }

        //
        if (use_colglaz == 1)
        {
           if (ppm_type == 0 && ParallelDescriptor::IOProcessor())
               std::cout << "Nyx::setting use_colglaz = 1 with ppm_type = 0 \n";
           if (ppm_type != 0)
               amrex::Error("Nyx::ppm_type must be 0 with use_colglaz = 1");
        }

        // ppm_flatten_before_integrals is only done for ppm_type != 0
        if (ppm_type == 0 && ppm_flatten_before_integrals > 0)
        {
            std::cerr << "ppm_flatten_before_integrals > 0 not implemented for ppm_type != 0 \n";
            amrex::Error();
        }

        if (version_2 > 0 && ppm_type == 0)
           amrex::Error("Nyx::version_2 only defined for ppm_type = 1");

        if (version_2 !=0 && version_2 != 1 && version_2 != 2 && version_2 != 3)
           amrex::Error("Nyx:: don't know what to do with version_2 flag");

        // Make sure ppm_type is set correctly.
        if (ppm_type != 0 && ppm_type != 1 && ppm_type != 2)
        {
           amrex::Error("Nyx::ppm_type must be 0, 1 or 2");
        }
    }

    // Make sure not to call refluxing if we're not actually doing any hydro.
    if (do_hydro == 0) do_reflux = 0;

#ifdef GRAVITY
    pp_nyx.get("do_grav", do_grav);
#endif

    read_particle_params();

    read_init_params();

    pp_nyx.query("write_parameter_file",write_parameters_in_plotfile);
    pp_nyx.query("print_fortran_warnings",print_fortran_warnings);

    read_comoving_params();

    if (pp_nyx.contains("plot_z_values"))
    {
      int num_z_values = pp_nyx.countval("plot_z_values");
      plot_z_values.resize(num_z_values);
      pp_nyx.queryarr("plot_z_values",plot_z_values,0,num_z_values);
    }

    if (pp_nyx.contains("analysis_z_values"))
    {
      int num_z_values = pp_nyx.countval("analysis_z_values");
      analysis_z_values.resize(num_z_values);
      pp_nyx.queryarr("analysis_z_values",analysis_z_values,0,num_z_values);
    }

    // How often do we want to write x,y,z 2-d slices of S_new
    pp_nyx.query("slice_int",    slice_int);
    pp_nyx.query("slice_file",   slice_file);
    pp_nyx.query("slice_nfiles", slice_nfiles);

    pp_nyx.query("gimlet_int", gimlet_int);

#ifdef AGN
    pp_nyx.query("mass_halo_min", mass_halo_min);
    pp_nyx.query("mass_seed", mass_seed);
#endif
}

Nyx::Nyx ()
{
    BL_PROFILE("Nyx::Nyx()");
#ifndef NO_HYDRO
    if (do_hydro == 1)
    {
        flux_reg = 0;
    }
#endif
    fine_mask = 0;
}

Nyx::Nyx (Amr&            papa,
          int             lev,
          const Geometry& level_geom,
          const BoxArray& bl,
          const DistributionMapping& dm,
          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time)
{
    BL_PROFILE("Nyx::Nyx(Amr)");
    build_metrics();
    fine_mask = 0;

#ifndef NO_HYDRO
    if (do_hydro == 1)
    {
        flux_reg = 0;
        if (level > 0 && do_reflux)
            flux_reg = new FluxRegister(grids, dmap, crse_ratio, level, NUM_STATE);
    }
#endif

#ifdef GRAVITY
    // Initialize to zero here in case we run with do_grav = false.
    MultiFab& new_grav_mf = get_new_data(Gravity_Type);
    new_grav_mf.setVal(0);

    if (do_grav)
    {
        // gravity is a static object, only alloc if not already there
        if (gravity == 0) {
          gravity = new Gravity(parent, parent->finestLevel(), &phys_bc, Density);
	}

        gravity->install_level(level, this);

        if (verbose && level == 0 && ParallelDescriptor::IOProcessor()) {
            std::cout << "Setting the gravity type to "
                      << gravity->get_gravity_type() << '\n';
	}
   }
#endif

#ifdef FORCING
    const Real* prob_lo = geom.ProbLo();
    const Real* prob_hi = geom.ProbHi();

    if (do_forcing)
    {
        // forcing is a static object, only alloc if not already there
        if (forcing == 0)
           forcing = new StochasticForcing();

        forcing->init(BL_SPACEDIM, prob_lo, prob_hi);
    }
#endif

    // Initialize the "a" variable
    if (level == 0 && time == 0.0 && old_a_time < 0.)
    {
       old_a_time = 0.0;
       new_a_time = 0.0;

       old_a = 1.0 / (1.0 + initial_z);
       new_a = old_a;
    }

#ifdef HEATCOOL
     // Initialize "this_z" in the atomic_rates_module
    if (heat_cool_type == 3 || heat_cool_type == 5 || heat_cool_type == 7)
         fort_interp_to_this_z(&initial_z);
#endif

#ifdef AGN
     // Initialize the uniform(0,1) random number generator.
     init_uniform01_rng();
#endif
}

Nyx::~Nyx ()
{
#ifndef NO_HYDRO
    if (do_hydro == 1)
        delete flux_reg;
#endif
    delete fine_mask;
}

void
Nyx::restart (Amr&     papa,
              istream& is,
              bool     b_read_special)
{
    BL_PROFILE("Nyx::restart()");
    AmrLevel::restart(papa, is, b_read_special);

    build_metrics();


    // get the elapsed CPU time to now;
    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
      // get elapsed CPU time
      std::ifstream CPUFile;
      std::string FullPathCPUFile = parent->theRestartFile();
      FullPathCPUFile += "/CPUtime";
      CPUFile.open(FullPathCPUFile.c_str(), std::ios::in);

      CPUFile >> previousCPUTimeUsed;
      CPUFile.close();

      std::cout << "read CPU time: " << previousCPUTimeUsed << "\n";
    }

#ifndef NO_HYDRO
    if (do_hydro == 1)
    {
        BL_ASSERT(flux_reg == 0);
        if (level > 0 && do_reflux)
            flux_reg = new FluxRegister(grids, dmap, crse_ratio, level, NUM_STATE);
    }
#endif

#ifdef GRAVITY
    if (do_grav && level == 0)
    {
        BL_ASSERT(gravity == 0);
        gravity = new Gravity(parent, parent->finestLevel(), &phys_bc, Density);
    }
#endif

#ifdef FORCING
    const Real* prob_lo = geom.ProbLo();
    const Real* prob_hi = geom.ProbHi();

    if (do_forcing)
    {
        // forcing is a static object, only alloc if not already there
        if (forcing == 0)
           forcing = new StochasticForcing();

        forcing->init(BL_SPACEDIM, prob_lo, prob_hi);
    }
#endif
}

void
Nyx::build_metrics ()
{
}

void
Nyx::setTimeLevel (Real time,
                   Real dt_old,
                   Real dt_new)
{
    if (ParallelDescriptor::IOProcessor()) {
       std::cout << "Setting the current time in the state data to "
                 << parent->cumTime() << std::endl;
    }
    AmrLevel::setTimeLevel(time, dt_old, dt_new);
}

void
Nyx::init (AmrLevel& old)
{
    BL_PROFILE("Nyx::init(old)");
    Nyx* old_level = (Nyx*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new = parent->dtLevel(level);
#ifdef NO_HYDRO
    Real cur_time = old_level->state[PhiGrav_Type].curTime();
    Real prev_time = old_level->state[PhiGrav_Type].prevTime();
#else
    Real cur_time = old_level->state[State_Type].curTime();
    Real prev_time = old_level->state[State_Type].prevTime();
#endif

    Real dt_old = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

#ifndef NO_HYDRO
    if (do_hydro == 1)
    {
        MultiFab& S_new = get_new_data(State_Type);
        MultiFab& D_new = get_new_data(DiagEOS_Type);

        for (FillPatchIterator
                 fpi(old, S_new, 0, cur_time,   State_Type, 0, NUM_STATE),
                dfpi(old, D_new, 0, cur_time, DiagEOS_Type, 0, D_new.nComp());
                fpi.isValid() && dfpi.isValid();
                ++fpi,++dfpi)
        {
            FArrayBox&  tmp =  fpi();
            FArrayBox& dtmp = dfpi();
            S_new[fpi].copy(tmp);
            D_new[fpi].copy(dtmp);
        }
    }
#endif

#ifdef GRAVITY
    MultiFab& Phi_new = get_new_data(PhiGrav_Type);
    for (FillPatchIterator fpi(old, Phi_new, 0, cur_time, PhiGrav_Type, 0, 1);
         fpi.isValid(); ++fpi)
    {
        Phi_new[fpi].copy(fpi());
    }
#endif

#ifdef SDC
    MultiFab& IR_new = get_new_data(SDC_IR_Type);
    for (FillPatchIterator fpi(old, IR_new, 0, cur_time, SDC_IR_Type, 0, 1);
         fpi.isValid(); ++fpi)
    {
        IR_new[fpi].copy(fpi());
    }
#endif
}

//
// This version inits the data on a new level that did not
// exist before regridding.
//
void
Nyx::init ()
{
    BL_PROFILE("Nyx::init()");
    Real dt        = parent->dtLevel(level);
#ifdef NO_HYDRO
    Real cur_time  = get_level(level-1).state[PhiGrav_Type].curTime();
    Real prev_time = get_level(level-1).state[PhiGrav_Type].prevTime();
#else
    Real cur_time  = get_level(level-1).state[State_Type].curTime();
    Real prev_time = get_level(level-1).state[State_Type].prevTime();
#endif
    Real dt_old    = (cur_time - prev_time) / (Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time, dt_old, dt);

#ifndef NO_HYDRO
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);
    FillCoarsePatch(S_new, 0, cur_time,   State_Type, 0, S_new.nComp());
    FillCoarsePatch(D_new, 0, cur_time, DiagEOS_Type, 0, D_new.nComp());
#endif

#ifdef GRAVITY
    MultiFab& Phi_new = get_new_data(PhiGrav_Type);
    FillCoarsePatch(Phi_new, 0, cur_time, PhiGrav_Type, 0, Phi_new.nComp());
#endif

    // We set dt to be large for this new level to avoid screwing up
    // computeNewDt.
    parent->setDtLevel(1.e100, level);
}

Real
Nyx::initial_time_step ()
{
    BL_PROFILE("Nyx::initial_time_step()");
    Real dummy_dt = 0;
    Real init_dt = 0;

    if (initial_dt > 0)
    {
        init_dt = initial_dt;
    }
    else
    {
        init_dt = init_shrink * est_time_step(dummy_dt);
    }

    bool dt_changed = false;
    if (level == 0 && plot_z_values.size() > 0)
        plot_z_est_time_step(init_dt,dt_changed);

    if (level == 0 && analysis_z_values.size() > 0)
        analysis_z_est_time_step(init_dt,dt_changed);

    return init_dt;
}

Real
Nyx::est_time_step (Real dt_old)
{
    BL_PROFILE("Nyx::est_time_step()");
    if (fixed_dt > 0)
        return fixed_dt;

    // This is just a dummy value to start with
    Real est_dt = 1.0e+200;

#ifndef NO_HYDRO
    const MultiFab& stateMF = get_new_data(State_Type);
#endif

#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif

#ifndef NO_HYDRO
    if (do_hydro)
    {
        Real a = get_comoving_a(cur_time);
        const Real* dx = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel reduction(min:est_dt)
#endif
	{
          Real dt = 1.e200;
	  for (MFIter mfi(stateMF,true); mfi.isValid(); ++mfi)
	    {
	      const Box& box = mfi.tilebox();

	      fort_estdt
                (BL_TO_FORTRAN(stateMF[mfi]), box.loVect(), box.hiVect(), dx,
                 &dt, &a);
	    }
          est_dt = std::min(est_dt, dt);
	}

        // If in comoving coordinates, then scale dt (based on u and c) by a
        est_dt *= a;

        ParallelDescriptor::ReduceRealMin(est_dt);
        est_dt *= cfl;
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "...estdt from hydro at level "
                      << level << ": "
                      << est_dt << '\n';
    }
#endif

#ifdef GRAVITY
    particle_est_time_step(est_dt);
#endif

    if (level==0)
        comoving_est_time_step(cur_time,est_dt);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Nyx::est_time_step at level "
                  << level
                  << ":  estdt = "
                  << est_dt << '\n';

    return est_dt;
}

void
Nyx::computeNewDt (int                   finest_level,
                   int                   sub_cycle,
                   Vector<int>&           n_cycle,
                   const Vector<IntVect>& ref_ratio,
                   Vector<Real>&          dt_min,
                   Vector<Real>&          dt_level,
                   Real                  stop_time,
                   int                   post_regrid_flag)
{
    BL_PROFILE("Nyx::computeNewDt()");
    //
    // We are at the start of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    int i;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        Nyx& adv_level = get_level(i);
        dt_min[i] = adv_level.est_time_step(dt_level[i]);
    }

    if (fixed_dt <= 0.0)
    {
        if (post_regrid_flag == 1)
        {
            //
            // Limit dt's by pre-regrid dt
            //
            for (i = 0; i <= finest_level; i++)
            {
                dt_min[i] = std::min(dt_min[i], dt_level[i]);
            }
            //
            // Find the minimum over all levels
            //
            for (i = 0; i <= finest_level; i++)
            {
                n_factor *= n_cycle[i];
                dt_0 = std::min(dt_0, n_factor * dt_min[i]);
            }
        }
        else
        {
            bool sub_unchanged=true;
            if ((parent->maxLevel() > 0) && (level == 0) &&
                (parent->subcyclingMode() == "Optimal") &&
                (parent->okToRegrid(level) || parent->levelSteps(0) == 0) )
            {
                int new_cycle[finest_level+1];
                for (i = 0; i <= finest_level; i++)
                    new_cycle[i] = n_cycle[i];
                // The max allowable dt
                Real dt_max[finest_level+1];
                for (i = 0; i <= finest_level; i++)
                {
                    dt_max[i] = dt_min[i];
                }
                // find the maximum number of cycles allowed.
                int cycle_max[finest_level+1];
                cycle_max[0] = 1;
                for (i = 1; i <= finest_level; i++)
                {
                    cycle_max[i] = parent->MaxRefRatio(i-1);
                }
                // estimate the amout of work to advance each level.
                Real est_work[finest_level+1];
                for (i = 0; i <= finest_level; i++)
                {
                    est_work[i] = parent->getLevel(i).estimateWork();
                }
                // this value will be used only if the subcycling pattern is changed.
                dt_0 = parent->computeOptimalSubcycling(finest_level+1, new_cycle, dt_max, est_work, cycle_max);
                for (i = 0; i <= finest_level; i++)
                {
                    if (n_cycle[i] != new_cycle[i])
                    {
                        sub_unchanged = false;
                        n_cycle[i] = new_cycle[i];
                    }
                }

            }

            if (sub_unchanged)
            //
            // Limit dt's by change_max * old dt
            //
            {
                for (i = 0; i <= finest_level; i++)
                {
                    if (verbose && ParallelDescriptor::IOProcessor())
                    {
                        if (dt_min[i] > change_max*dt_level[i])
                        {
                            std::cout << "Nyx::compute_new_dt : limiting dt at level "
                                      << i << '\n';
                            std::cout << " ... new dt computed: " << dt_min[i]
                                      << '\n';
                            std::cout << " ... but limiting to: "
                                      << change_max * dt_level[i] << " = " << change_max
                                      << " * " << dt_level[i] << '\n';
                        }
                    }

                    dt_min[i] = std::min(dt_min[i], change_max * dt_level[i]);
                }
                //
                // Find the minimum over all levels
                //
                for (i = 0; i <= finest_level; i++)
                {
                    n_factor *= n_cycle[i];
                    dt_0 = std::min(dt_0, n_factor * dt_min[i]);
                }
            }
            else
            {
                if (verbose && ParallelDescriptor::IOProcessor())
                {
                   std::cout << "Nyx: Changing subcycling pattern. New pattern:\n";
                   for (i = 1; i <= finest_level; i++)
                    std::cout << "   Lev / n_cycle: " << i << " " << n_cycle[i] << '\n';
                }
            }
        }
    }
    else
    {
        dt_0 = fixed_dt;
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001 * dt_0;
#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif
    if (stop_time >= 0.0)
    {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    // Shrink the time step if necessary in order to hit the next plot_z_value
    if (level == 0 && ( plot_z_values.size() > 0 || analysis_z_values.size() > 0 ) )
    {
        bool dt_changed_plot     = false;
        bool dt_changed_analysis = false;
   
        if (plot_z_values.size() > 0)
           plot_z_est_time_step(dt_0,dt_changed_plot);

        if (analysis_z_values.size() > 0)
           analysis_z_est_time_step(dt_0,dt_changed_analysis);

        // Update the value of a if we didn't change dt in the call to plot_z_est_time_step or analysis_z_est_time_step.
        // If we didn't change dt there, then we have already done the integration.
        // If we did    change dt there, then we need to re-integrate here.
        if (dt_changed_plot || dt_changed_analysis)
            integrate_comoving_a(cur_time,dt_0);
    }
    else 
    {
        integrate_comoving_a(cur_time,dt_0);
    }
     
    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }
}

void
Nyx::computeInitialDt (int                   finest_level,
                       int                   sub_cycle,
                       Vector<int>&           n_cycle,
                       const Vector<IntVect>& ref_ratio,
                       Vector<Real>&          dt_level,
                       Real                  stop_time)
{
    BL_PROFILE("Nyx::computeInitialDt()");
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    int i;
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    if (parent->subcyclingMode() == "Optimal")
    {
        int new_cycle[finest_level+1];
        for (i = 0; i <= finest_level; i++)
            new_cycle[i] = n_cycle[i];
        Real dt_max[finest_level+1];
        for (i = 0; i <= finest_level; i++)
        {
            dt_max[i] = get_level(i).initial_time_step();
        }
        // Find the maximum number of cycles allowed
        int cycle_max[finest_level+1];
        cycle_max[0] = 1;
        for (i = 1; i <= finest_level; i++)
        {
            cycle_max[i] = parent->MaxRefRatio(i-1);
        }
        // estimate the amout of work to advance each level.
        Real est_work[finest_level+1];
        for (i = 0; i <= finest_level; i++)
        {
            est_work[i] = parent->getLevel(i).estimateWork();
        }
        dt_0 = parent->computeOptimalSubcycling(finest_level+1, new_cycle, dt_max, est_work, cycle_max);
        for (i = 0; i <= finest_level; i++)
        {
            n_cycle[i] = new_cycle[i];
        }
        if (verbose && ParallelDescriptor::IOProcessor() && finest_level > 0)
        {
           std::cout << "Nyx: Initial subcycling pattern:\n";
           for (i = 0; i <= finest_level; i++)
               std::cout << "Level " << i << ": " << n_cycle[i] << '\n';
        }
    }
    else
    {
        for (i = 0; i <= finest_level; i++)
        {
            dt_level[i] = get_level(i).initial_time_step();
            n_factor *= n_cycle[i];
            dt_0 = std::min(dt_0, n_factor * dt_level[i]);
        }
    }
    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001 * dt_0;
#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif
    if (stop_time >= 0)
    {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }

    integrate_comoving_a(cur_time,dt_0);
}

bool
Nyx::writePlotNow ()
{
    BL_PROFILE("Nyx::writePlotNow()");
    if (level > 0)
        amrex::Error("Should only call writePlotNow at level 0!");

    bool found_one = false;

    if (plot_z_values.size() > 0)
    {

#ifdef NO_HYDRO
        Real prev_time = state[PhiGrav_Type].prevTime();
        Real  cur_time = state[PhiGrav_Type].curTime();
#else
        Real prev_time = state[State_Type].prevTime();
        Real  cur_time = state[State_Type].curTime();
#endif
        Real a_old = get_comoving_a(prev_time);
        Real z_old = (1. / a_old) - 1.;

        Real a_new = get_comoving_a( cur_time);
        Real z_new = (1. / a_new) - 1.;

        for (int i = 0; i < plot_z_values.size(); i++)
        {
            if (std::abs(z_new - plot_z_values[i]) < (0.01 * (z_old - z_new)) )
                found_one = true;
        }
    }

    if (found_one) {
        return true;
    } else {
        return false;
    }
}

bool
Nyx::doAnalysisNow ()
{
    BL_PROFILE("Nyx::doAnalysisNow()");
    if (level > 0)
        amrex::Error("Should only call doAnalysisNow at level 0!");

    bool found_one = false;

    if (analysis_z_values.size() > 0)
    {

#ifdef NO_HYDRO
        Real prev_time = state[PhiGrav_Type].prevTime();
        Real  cur_time = state[PhiGrav_Type].curTime();
#else
        Real prev_time = state[State_Type].prevTime();
        Real  cur_time = state[State_Type].curTime();
#endif
        Real a_old = get_comoving_a(prev_time);
        Real z_old = (1. / a_old) - 1.;

        Real a_new = get_comoving_a( cur_time);
        Real z_new = (1. / a_new) - 1.;

        for (int i = 0; i < analysis_z_values.size(); i++)
        {
            if (std::abs(z_new - analysis_z_values[i]) < (0.01 * (z_old - z_new)) )
                found_one = true;
        }
    }

    if (found_one) {
        return true;
    } else {
        return false;
    }
}

void
Nyx::do_energy_diagnostics ()
{
    // nothing to see here, folks
}

void
Nyx::post_timestep (int iteration)
{
if(wtid()==0){
    BL_PROFILE("Nyx::post_timestep()");
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();
    const int ncycle = parent->nCycle(level);

    //
    // Remove virtual particles at this level if we have any.
    //
    remove_virtual_particles();

    //
    // Remove Ghost particles on the final iteration
    //
    if (iteration == ncycle)
        remove_ghost_particles();

    //
    // Sync up if we're level 0 or if we have particles that may have moved
    // off the next finest level and need to be added to our own level.
    //
    if ((iteration < ncycle and level < finest_level) || level == 0)
    {
        for (int i = 0; i < theActiveParticles().size(); i++)
            theActiveParticles()[i]->Redistribute(level,
                                                  theActiveParticles()[i]->finestLevel(),
                                                  iteration);
    }

#ifndef NO_HYDRO
    if (do_reflux && level < finest_level)
    {
        MultiFab& S_new_crse = get_new_data(State_Type);
#ifdef GRAVITY
        MultiFab drho_and_drhoU;
#ifdef CGRAV
        if (do_grav &&
            (gravity->get_gravity_type() == "PoissonGrav"||gravity->get_gravity_type() == "CompositeGrav"))
#else
        if (do_grav && gravity->get_gravity_type() == "PoissonGrav")
#endif
        {
            // Define the update to rho and rhoU due to refluxing.
            drho_and_drhoU.define(grids, dmap, BL_SPACEDIM + 1, 0);
            MultiFab::Copy(drho_and_drhoU, S_new_crse, Density, 0,
                           BL_SPACEDIM + 1, 0);
            drho_and_drhoU.mult(-1.0);
        }
#endif // GRAVITY

        //We must reflux if the next finer level is subcycled relative to this level;
        //   otherwise the reflux was done as part of the multilevel advance
        if (parent->nCycle(level+1) != 1)
          reflux();

        // We need to do this before anything else because refluxing changes the
        // values of coarse cells underneath fine grids with the assumption
        // they'll be over-written by averaging down
        if (level < finest_level)
            average_down();

        // This needs to be done after any changes to the state from refluxing.
        enforce_nonnegative_species(S_new_crse);

#ifdef GRAVITY
#ifdef CGRAV
        if (do_grav &&
            (gravity->get_gravity_type() == "PoissonGrav"||gravity->get_gravity_type() == "CompositeGrav")
            && gravity->get_no_sync() == 0)
#else
        if (do_grav && gravity->get_gravity_type() == "PoissonGrav" && gravity->get_no_sync() == 0)
#endif
        {
            MultiFab::Add(drho_and_drhoU, S_new_crse, Density, 0, BL_SPACEDIM+1, 0);

            MultiFab dphi(grids, dmap, 1, 0);
            dphi.setVal(0);

            gravity->reflux_phi(level, dphi);

            // Compute (cross-level) gravity sync based on drho, dphi
            Vector<std::unique_ptr<MultiFab> > grad_delta_phi_cc(finest_level - level + 1);
            for (int lev = level; lev <= finest_level; lev++)
            {
                grad_delta_phi_cc[lev-level].reset(
                                      new MultiFab(get_level(lev).boxArray(),
						   get_level(lev).DistributionMap(),
						   BL_SPACEDIM, 0));
                grad_delta_phi_cc[lev-level]->setVal(0);
            }

            gravity->gravity_sync(level,finest_level,iteration,ncycle,drho_and_drhoU,dphi,
				  amrex::GetVecOfPtrs(grad_delta_phi_cc));
            dphi.clear();

            for (int lev = level; lev <= finest_level; lev++)
            {
                Real dt_lev = parent->dtLevel(lev);
                MultiFab&  S_new_lev = get_level(lev).get_new_data(State_Type);
                Real cur_time = state[State_Type].curTime();
                Real a_new = get_comoving_a(cur_time);

                const auto& ba = get_level(lev).boxArray();
                const auto& dm = get_level(lev).DistributionMap();
                MultiFab grad_phi_cc(ba, dm, BL_SPACEDIM, 0);
                gravity->get_new_grav_vector(lev, grad_phi_cc, cur_time);

#ifdef _OPENMP
#pragma omp parallel
#endif
		{
		  FArrayBox sync_src;
		  FArrayBox dstate;

		  for (MFIter mfi(S_new_lev,true); mfi.isValid(); ++mfi)
                  {
                    const Box& bx = mfi.tilebox();
                    dstate.resize(bx, BL_SPACEDIM + 1);
                    if (lev == level)
                    {
		      dstate.copy(drho_and_drhoU[mfi]);
                    }
                    else
                    {
		      dstate.setVal(0);
                    }

                    // Compute sync source
                    sync_src.resize(bx, BL_SPACEDIM+1);
                    int i = mfi.index();
                    fort_syncgsrc
                        (bx.loVect(), bx.hiVect(), BL_TO_FORTRAN(grad_phi_cc[i]),
                         BL_TO_FORTRAN((*grad_delta_phi_cc[lev-level])[i]),
                         BL_TO_FORTRAN(S_new_lev[i]), BL_TO_FORTRAN(dstate),
                         BL_TO_FORTRAN(sync_src), &a_new, dt_lev);

                    sync_src.mult(0.5 * dt_lev);
                    S_new_lev[mfi].plus(sync_src, 0, Xmom, BL_SPACEDIM);
                    S_new_lev[mfi].plus(sync_src, BL_SPACEDIM, Eden, 1);
		  }
		}
	    }
        }
#endif
    }
#endif // end ifndef NO_HYDRO

    if (level < finest_level)
        average_down();

    if (level == 0)
    {
        int nstep = parent->levelSteps(0);
#ifndef NO_HYDRO
        if ( (do_hydro == 1) && (sum_interval > 0) && (nstep % sum_interval == 0))
        {
            sum_integrated_quantities();
        }
#endif
        write_info();

#if BL_USE_MPI
        // Memory monitoring:
        MemInfo* mInfo = MemInfo::GetInstance();
        char info[32];
        snprintf(info, sizeof(info), "Step %4d", nstep);
        mInfo->LogSummary(info);
#endif
    }

#ifndef NO_HYDRO
    if (do_hydro)
    {
       MultiFab& S_new = get_new_data(State_Type);
       MultiFab& D_new = get_new_data(DiagEOS_Type);

       // First reset internal energy before call to compute_temp
       MultiFab reset_e_src(S_new.boxArray(), S_new.DistributionMap(), 1, NUM_GROW);
       reset_e_src.setVal(0.0);
       reset_internal_energy(S_new,D_new,reset_e_src);

       // Re-compute temperature after all the other updates.
       compute_new_temp(S_new,D_new);
    }
#endif
}
}

void
Nyx::post_restart ()
{
    BL_PROFILE("Nyx::post_restart()");
    if (level == 0)
        particle_post_restart(parent->theRestartFile());

    if (level == 0)
        comoving_a_post_restart(parent->theRestartFile());

    if (inhomo_reion) init_zhi();

#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif

    // Update the value of a only if restarting from chk00000
    //   (special case for which computeNewDt is *not* called from Amr::coarseTimeStep)
    if (level == 0 && cur_time == 0.0)
        integrate_comoving_a(cur_time,parent->dtLevel(0));

#ifdef TISF
     int blub = parent->finestLevel();
     fort_set_finest_level(&blub);
#endif

#ifdef GRAVITY

    if (do_grav)
    {
        if (level == 0)
        {
            for (int lev = 0; lev <= parent->finestLevel(); lev++)
            {
                AmrLevel& this_level = get_level(lev);
                gravity->install_level(lev, &this_level);
            }

            gravity->set_mass_offset(cur_time);

            if (
#ifdef CGRAV
            (gravity->get_gravity_type() == "PoissonGrav"||gravity->get_gravity_type() == "CompositeGrav")
#else
	    gravity->get_gravity_type() == "PoissonGrav"
#endif
)
            {
                // Do multilevel solve here.  We now store phi in the checkpoint file so we can use it
                //  at restart.
                int ngrow_for_solve = 1;
                int use_previous_phi_as_guess = 1;
                gravity->multilevel_solve_for_new_phi(0,parent->finestLevel(),ngrow_for_solve,use_previous_phi_as_guess);

#ifndef AGN
                if (do_dm_particles)
#endif
                {
                    for (int k = 0; k <= parent->finestLevel(); k++)
                    {
                        const auto& ba = get_level(k).boxArray();
                        const auto& dm = get_level(k).DistributionMap();
                        MultiFab grav_vec_new(ba, dm, BL_SPACEDIM, 0);
                        gravity->get_new_grav_vector(k, grav_vec_new, cur_time);
                    }
                }
            }
        }
    }
#endif

#ifdef FORCING
    if (do_forcing)
    {
        if (level == 0)
           forcing_post_restart(parent->theRestartFile());
    }
#endif

#ifndef NO_HYDRO
    if (level == 0)
    {
       // Need to compute this *before* regridding in case this is needed
       compute_average_density();
       set_small_values();
    }
#endif
}

#ifndef NO_HYDRO
void
Nyx::set_small_values ()
{
       if (do_hydro == 0) {
          return;
       }

       Real small_pres;

       const Real cur_time = state[State_Type].curTime();
       Real a = get_comoving_a(cur_time);

       Real average_temperature;
       compute_average_temperature(average_temperature);
       //
       // Get the number of species from the network model.
       //
       fort_get_num_spec(&NumSpec);
       fort_get_num_aux (&NumAux);

       fort_set_small_values
            (&average_gas_density, &average_temperature,
             &a,  &small_dens, &small_temp, &small_pres);

       if (verbose && ParallelDescriptor::IOProcessor())
       {
          std::cout << "... setting small_dens to " << small_dens << '\n';
          std::cout << "... setting small_temp to " << small_temp << '\n';
          std::cout << "... setting small_pres to " << small_pres << '\n';
       }
}
#endif

void
Nyx::postCoarseTimeStep (Real cumtime)
{
   BL_PROFILE("Nyx::postCoarseTimeStep()");

   AmrLevel::postCoarseTimeStep(cumtime);

#ifdef AGN
   halo_find(parent->dtLevel(level));
#endif 

#ifdef GIMLET
   LyA_statistics();
#endif

    //
    // postCoarseTimeStep() is only called by level 0.
    //
    if (Nyx::theDMPC() && particle_move_type == "Random")
        particle_move_random();

   int nstep = parent->levelSteps(0);

   if (slice_int > -1 && nstep%slice_int == 0)
   {
      BL_PROFILE("Nyx::postCoarseTimeStep: get_all_slice_data");

    if(slice_int != 2) {
      const Real* dx        = geom.CellSize();

      MultiFab& S_new = get_new_data(State_Type);
      MultiFab& D_new = get_new_data(DiagEOS_Type);

      Real x_coord = (geom.ProbLo()[0] + geom.ProbHi()[0]) / 2 + dx[0]/2;
      Real y_coord = (geom.ProbLo()[1] + geom.ProbHi()[1]) / 2 + dx[1]/2;
      Real z_coord = (geom.ProbLo()[2] + geom.ProbHi()[2]) / 2 + dx[2]/2;

      if (ParallelDescriptor::IOProcessor()) {
         std::cout << "Outputting slices at x = " << x_coord << "; y = " << y_coord << "; z = " << z_coord << std::endl;
      }

      const std::string& slicefilename = amrex::Concatenate(slice_file, nstep);
      UtilCreateCleanDirectory(slicefilename, true);

      int nfiles_current = amrex::VisMF::GetNOutFiles();
      amrex::VisMF::SetNOutFiles(slice_nfiles);

      // Slice state data
      std::unique_ptr<MultiFab> x_slice = amrex::get_slice_data(0, x_coord, S_new, geom, 0, S_new.nComp()-2);
      std::unique_ptr<MultiFab> y_slice = amrex::get_slice_data(1, y_coord, S_new, geom, 0, S_new.nComp()-2);
      std::unique_ptr<MultiFab> z_slice = amrex::get_slice_data(2, z_coord, S_new, geom, 0, S_new.nComp()-2);

      std::string xs = slicefilename + "/State_x";
      std::string ys = slicefilename + "/State_y";
      std::string zs = slicefilename + "/State_z";

      {
        BL_PROFILE("Nyx::postCoarseTimeStep: writeXSlice");
        amrex::VisMF::Write(*x_slice, xs);
      }
      {
        BL_PROFILE("Nyx::postCoarseTimeStep: writeYSlice");
        amrex::VisMF::Write(*y_slice, ys);
      }
      {
        BL_PROFILE("Nyx::postCoarseTimeStep: writeZSlice");
        amrex::VisMF::Write(*z_slice, zs);
      }

      // Slice diag_eos
      x_slice = amrex::get_slice_data(0, x_coord, D_new, geom, 0, D_new.nComp());
      y_slice = amrex::get_slice_data(1, y_coord, D_new, geom, 0, D_new.nComp());
      z_slice = amrex::get_slice_data(2, z_coord, D_new, geom, 0, D_new.nComp());

      xs = slicefilename + "/Diag_x";
      ys = slicefilename + "/Diag_y";
      zs = slicefilename + "/Diag_z";

      {
        BL_PROFILE("Nyx::postCoarseTimeStep: writeDiagSlices");
        amrex::VisMF::Write(*x_slice, xs);
        amrex::VisMF::Write(*y_slice, ys);
        amrex::VisMF::Write(*z_slice, zs);
      }

      amrex::VisMF::SetNOutFiles(nfiles_current);

      if (ParallelDescriptor::IOProcessor()) {
         std::cout << "Done with slices." << std::endl;
      }


    } else {

      MultiFab& S_new = get_new_data(State_Type);
      MultiFab& D_new = get_new_data(DiagEOS_Type);

      const std::string& slicefilename = amrex::Concatenate(slice_file, nstep);
      UtilCreateCleanDirectory(slicefilename, true);

      int nfiles_current = amrex::VisMF::GetNOutFiles();
      amrex::VisMF::SetNOutFiles(slice_nfiles);

      int maxBoxSize(64);
      amrex::Vector<std::string> SMFNames(3);
      SMFNames[0] = slicefilename + "/State_x";
      SMFNames[1] = slicefilename + "/State_y";
      SMFNames[2] = slicefilename + "/State_z";
      amrex::Vector<std::string> DMFNames(3);
      DMFNames[0] = slicefilename + "/Diag_x";
      DMFNames[1] = slicefilename + "/Diag_y";
      DMFNames[2] = slicefilename + "/Diag_z";

      for(int dir(0); dir < 3; ++dir) {
        Box sliceBox(geom.Domain());
        int dir_coord = geom.ProbLo()[dir] + (geom.Domain().length(dir) / 2);
        amrex::Print() << "Outputting slices at dir_coord[" << dir << "] = " << dir_coord << '\n';
        sliceBox.setSmall(dir, dir_coord);
        sliceBox.setBig(dir, dir_coord);
        BoxArray sliceBA(sliceBox);
        sliceBA.maxSize(maxBoxSize);
        DistributionMapping sliceDM(sliceBA);

        MultiFab SSliceMF(sliceBA, sliceDM, S_new.nComp()-2, 0);
        SSliceMF.copy(S_new, 0, 0, SSliceMF.nComp());
        amrex::VisMF::Write(SSliceMF, SMFNames[dir]);

        MultiFab DSliceMF(sliceBA, sliceDM, D_new.nComp(), 0);
        DSliceMF.copy(D_new, 0, 0, DSliceMF.nComp());
        amrex::VisMF::Write(DSliceMF, DMFNames[dir]);
      }

      amrex::VisMF::SetNOutFiles(nfiles_current);

      if (ParallelDescriptor::IOProcessor()) {
         std::cout << "Done with slices." << std::endl;
      }

    }

   }
}

void
Nyx::post_regrid (int lbase,
                  int new_finest)
{
    BL_PROFILE("Nyx::post_regrid()");
#ifndef NO_HYDRO
#ifdef TISF
     fort_set_finest_level(&new_finest);
#endif
#endif

    if (level == lbase) {
        particle_redistribute(lbase, false);
    }

#ifdef GRAVITY

    int which_level_being_advanced = parent->level_being_advanced();

    bool do_grav_solve_here;

    if (which_level_being_advanced >= 0)
    {
        do_grav_solve_here = (level == which_level_being_advanced) && (lbase == which_level_being_advanced);
    } else {
        do_grav_solve_here = (level == lbase);
    }

    // Only do solve here if we will be using it in the timestep right after without re-solving,
    //      or if this is called from somewhere other than Amr::timeStep
    const Real cur_time = state[PhiGrav_Type].curTime();
    if (do_grav && (cur_time > 0) && do_grav_solve_here)
    {
#ifdef CGRAV
        if (gravity->get_gravity_type() == "PoissonGrav" || gravity->get_gravity_type() == "CompositeGrav")
#else
        if (gravity->get_gravity_type() == "PoissonGrav")
#endif
        {
            int ngrow_for_solve = parent->levelCount(level) + 1;
            int use_previous_phi_as_guess = 1;
            gravity->multilevel_solve_for_new_phi(level, new_finest, ngrow_for_solve, use_previous_phi_as_guess);
        }
    }
#endif
    delete fine_mask;
    fine_mask = 0;
}

void
Nyx::post_init (Real stop_time)
{
    BL_PROFILE("Nyx::post_init()");
    if (level > 0) {
        return;
    }

    // If we restarted from a plotfile, we need to reset the level_steps counter
    if ( ! parent->theRestartPlotFile().empty()) {
        parent->setLevelSteps(0,nsteps_from_plotfile);
    }

    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level - 1; k >= 0; --k) {
        get_level(k).average_down();
    }

#ifdef GRAVITY
    if (do_grav)
    {
        const Real cur_time = state[PhiGrav_Type].curTime();
        if
#ifdef CGRAV
            (gravity->get_gravity_type() == "PoissonGrav" ||
             gravity->get_gravity_type() == "CompositeGrav")
#else
	    (gravity->get_gravity_type() == "PoissonGrav")
#endif
        {
            //
            // Calculate offset before first multilevel solve.
            //
            gravity->set_mass_offset(cur_time);

            //
            // Solve on full multilevel hierarchy
            //
            int ngrow_for_solve = 1;
            gravity->multilevel_solve_for_new_phi(0, finest_level, ngrow_for_solve);
        }

        // Make this call just to fill the initial state data.
        for (int k = 0; k <= finest_level; k++)
        {
            const auto& ba = get_level(k).boxArray();
            const auto& dm = get_level(k).DistributionMap();
            MultiFab grav_vec_new(ba, dm, BL_SPACEDIM, 0);
            gravity->get_new_grav_vector(k, grav_vec_new, cur_time);
        }
    }
#endif

#ifndef NO_HYDRO
    if ( (do_hydro == 1) && (sum_interval > 0) && (parent->levelSteps(0) % sum_interval == 0) )
    {
        sum_integrated_quantities();
    }
    else
    {
        // Even if we don't call `sum_integrated_quantities` we need to compute
        // average_density before regridding
        compute_average_density();
    }

    if (do_hydro == 1)
    {
        set_small_values();
    }
#endif

    write_info();
}

int
Nyx::okToContinue ()
{
    if (level > 0) {
        return 1;
    }

    int test = 1;
    if (parent->dtLevel(0) < dt_cutoff) {
        test = 0;
    }

    if ((test == 1) && (final_a > 0))
    {
#ifdef NO_HYDRO
        Real cur_time = state[PhiGrav_Type].curTime();
#else
        Real cur_time = state[State_Type].curTime();
#endif
        Real a = get_comoving_a(cur_time);
        if (a >= final_a) test = 0;
        if (verbose && ParallelDescriptor::IOProcessor())
        {
            if (test == 0) {
                std::cout << "...a " << a
                          << " is greater than or equal to final_a " << final_a
                          << '\n';
	    }
        }
    }
    return test;
}

#ifdef AUX_UPDATE
void
Nyx::advance_aux (Real time,
                  Real dt)
{
    BL_PROFILE("Nyx::advance_aux()");
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... special update for auxiliary variables \n";

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_old,true); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        FArrayBox& old_fab = S_old[mfi];
        FArrayBox& new_fab = S_new[mfi];
        fort_auxupdate
            (BL_TO_FORTRAN(old_fab), BL_TO_FORTRAN(new_fab), box.loVect(),
             box.hiVect(), &dt);
    }
}
#endif

#ifndef NO_HYDRO
void
Nyx::reflux ()
{
    BL_PROFILE("Nyx::reflux()");

    BL_ASSERT(level<parent->finestLevel());

    get_flux_reg(level+1).Reflux(get_new_data(State_Type), 1.0, 0, 0, NUM_STATE,
                                 geom);
}
#endif // NO_HYDRO

void
Nyx::average_down ()
{
    BL_PROFILE("Nyx::average_down()");
    if (level == parent->finestLevel()) return;

#ifndef NO_HYDRO
    // With State_Type we do DiagEOS_Type
    average_down(State_Type);
#endif

#ifdef GRAVITY
    average_down(PhiGrav_Type);
    average_down(Gravity_Type);
#endif
}

#ifndef NO_HYDRO
void
Nyx::enforce_nonnegative_species (MultiFab& S_new)
{
    BL_PROFILE("Nyx::enforce_nonnegative_species()");
//#ifdef _OPENMP
//#pragma omp parallel
//#endif
//    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    for (RGIter rgi(RG_S,true); rgi.isValid(); ++rgi)
    {
        const Box& bx = rgi.tilebox();
        int f = rgi.currentRegion;
        int mfi = S_new.IndexArray()[f];
        fort_enforce_nonnegative_species
	  (BL_TO_FORTRAN(S_new[mfi]), bx.loVect(), bx.hiVect(),
	   &print_fortran_warnings);
    }
}

void
Nyx::enforce_consistent_e (MultiFab& S)
{
    BL_PROFILE("Nyx::enforce_consistent_e()");
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        const int* lo = box.loVect();
        const int* hi = box.hiVect();
        fort_enforce_consistent_e
	  (lo, hi, BL_TO_FORTRAN(S[mfi]));
    }
}
#endif

void
Nyx::average_down (int state_index)
{
    BL_PROFILE("Nyx::average_down(si)");
#ifndef NO_HYDRO
    // We average DiagEOS_Type when average_down is called with State_Type
    if (state_index == DiagEOS_Type) return;
#endif

    if (level == parent->finestLevel()) return;

    Nyx& fine_lev = get_level(level+1);

    const Geometry& fgeom = fine_lev.geom;
    const Geometry& cgeom =          geom;

#ifndef NO_HYDRO
    if (state_index == State_Type)
    {
        MultiFab& S_crse =          get_new_data(State_Type);
        MultiFab& S_fine = fine_lev.get_new_data(State_Type);

        amrex::average_down(S_fine, S_crse,
                            fgeom, cgeom,
                            0, S_fine.nComp(), fine_ratio);

        MultiFab& D_crse =          get_new_data(DiagEOS_Type);
        MultiFab& D_fine = fine_lev.get_new_data(DiagEOS_Type);

        amrex::average_down(D_fine, D_crse,
                            fgeom, cgeom,
                            0, D_fine.nComp(), fine_ratio);
    }
    else
#endif
    {
      MultiFab& S_crse = get_new_data(state_index);
      MultiFab& S_fine = fine_lev.get_new_data(state_index);

      const int num_comps = S_fine.nComp();

      amrex::average_down(S_fine,S_crse,fgeom,cgeom,0,num_comps,fine_ratio);
    }
}

void
Nyx::errorEst (TagBoxArray& tags,
               int          clearval,
               int          tagval,
               Real         time,
               int          n_error_buf,
               int          ngrow)
{
    BL_PROFILE("Nyx::errorEst()");
    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    Real avg;

    for (int j = 0; j < err_list.size(); j++)
    {
        auto mf = derive(err_list[j].name(), time, err_list[j].nGrow());

        BL_ASSERT(!(mf == 0));

        if (err_list[j].errType() == ErrorRec::UseAverage)
        {
            if (err_list[j].name() == "density")
            {
                avg = average_gas_density;
            }
            else if (err_list[j].name() == "particle_mass_density")
            {
                avg = average_dm_density;
#ifdef NEUTRINO_PARTICLES
                avg += average_neutr_density;
#endif
            }
            else if (err_list[j].name() == "total_density")
            {
                avg = average_total_density;
            }
#if 0
            else if (err_list[j].name() == "magvort")
            {
                avg = std::fabs(ave_lev_vorticity[level]);
                stddev = std_lev_vorticity[level];
                thresh = avg + std::max(stddev,avg);
                //std::cout << "errorEst, level " << level << ": magvort avg " << avg << ", stddev " << stddev
                //        << ", max " << std::max(stddev,avg) << ", thresh " << thresh << std::endl;
                thresh = std::max(thresh, 500.0);
                //std::cout << "errorEst, level " << level << ": thresh cut " << thresh << std::endl;
                avg = thresh;
            }
#endif
            else
            {
                if (ParallelDescriptor::IOProcessor())
                    std::cout << "Dont know the average of this variable "
                              << err_list[j].name() << '\n';
                avg = 0;
            }
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            Vector<int> itags;

            for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
            {
		// FABs
		FArrayBox&  datfab  = (*mf)[mfi];
		TagBox&     tagfab  = tags[mfi];

		// Box in physical space
		int         idx     = mfi.index();
		RealBox     gridloc = RealBox(grids[idx],geom.CellSize(),geom.ProbLo());

		// tile box
		const Box&  tilebx  = mfi.tilebox();

		//fab box
		const Box&  datbox  = datfab.box();

		// We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
		// So we are going to get a temporary integer array.
		tagfab.get_itags(itags, tilebx);

		// data pointer and index space
		int*        tptr    = itags.dataPtr();
		const int*  tlo     = tilebx.loVect();
		const int*  thi     = tilebx.hiVect();
		//
		const int*  lo      = tlo;
		const int*  hi      = thi;
		//
		const Real* xlo     = gridloc.lo();
		//
		Real*       dat     = datfab.dataPtr();
		const int*  dlo     = datbox.loVect();
		const int*  dhi     = datbox.hiVect();
		const int   ncomp   = datfab.nComp();

                if (err_list[j].errType() == ErrorRec::Standard)
                {
                    err_list[j].errFunc()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
                                          &clearval, dat, ARLIM(dlo), ARLIM(dhi),
                                          lo, hi, &ncomp, domain_lo, domain_hi,
                                          dx, xlo, prob_lo, &time, &level);
                }
                else if (err_list[j].errType() == ErrorRec::UseAverage)
                {
                   err_list[j].errFunc2()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
                                          &clearval, dat, ARLIM(dlo), ARLIM(dhi),
                                          lo, hi, &ncomp, domain_lo, domain_hi,
                                          dx, &level, &avg);
                }

                //
                // Don't forget to set the tags in the TagBox.
                //
                if (allow_untagging == 1)
                {
                    tagfab.tags_and_untags(itags, tilebx);
                }
                else
                {
                   tagfab.tags(itags, tilebx);
                }
            }
        }
    }
}

std::unique_ptr<MultiFab>
Nyx::derive (const std::string& name,
             Real               time,
             int                ngrow)
{
    BL_PROFILE("Nyx::derive()");
    if (name == "Rank")
    {
	std::unique_ptr<MultiFab> derive_dat (new MultiFab(grids, dmap, 1, 0));
        for (MFIter mfi(*derive_dat); mfi.isValid(); ++mfi)
        {
           (*derive_dat)[mfi].setVal(ParallelDescriptor::MyProc());
        }
        return derive_dat;
    } else {
        return particle_derive(name, time, ngrow);
    }
}

void
Nyx::derive (const std::string& name,
             Real               time,
             MultiFab&          mf,
             int                dcomp)
{
    BL_PROFILE("Nyx::derive(mf)");
    AmrLevel::derive(name, time, mf, dcomp);
}

void
Nyx::network_init ()
{
    fort_network_init();
}

#ifndef NO_HYDRO
void
Nyx::reset_internal_energy (MultiFab& S_new, MultiFab& D_new, MultiFab& reset_e_src)
{
    BL_PROFILE("Nyx::reset_internal_energy()");
    // Synchronize (rho e) and (rho E) so they are consistent with each other

    const Real  cur_time = state[State_Type].curTime();
    Real        a        = get_comoving_a(cur_time);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        reset_internal_e
            (bx.loVect(), bx.hiVect(),
             BL_TO_FORTRAN(S_new[mfi]), BL_TO_FORTRAN(D_new[mfi]),
	     BL_TO_FORTRAN(reset_e_src[mfi]),
             &print_fortran_warnings, &a);
    }
}
#endif

#ifndef NO_HYDRO
void
Nyx::compute_new_temp (MultiFab& S_new, MultiFab& D_new)
{
    BL_PROFILE("Nyx::compute_new_temp()");

    Real cur_time  = state[State_Type].curTime();
    Real a        = get_comoving_a(cur_time);

#ifdef HEATCOOL 
    if (heat_cool_type == 3 || heat_cool_type == 5 || heat_cool_type == 7) 
    {
       const Real z = 1.0/a - 1.0;
       fort_interp_to_this_z(&z);
    }
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        if (heat_cool_type == 7) {
          fort_compute_temp_vec
              (bx.loVect(), bx.hiVect(),
              BL_TO_FORTRAN(S_new[mfi]),
              BL_TO_FORTRAN(D_new[mfi]), &a,
               &print_fortran_warnings);
        } else {
            fort_compute_temp
              (bx.loVect(), bx.hiVect(),
              BL_TO_FORTRAN(S_new[mfi]),
              BL_TO_FORTRAN(D_new[mfi]), &a,
               &print_fortran_warnings);
        }
    }

    // Compute the maximum temperature
    Real max_temp = D_new.norm0(Temp_comp);

    int imax = -1;
    int jmax = -1;
    int kmax = -1;

    Real den_maxt;

    // Find the cell which has the maximum temp -- but only if not the first
    // time step because in the first time step too many points have the same
    // value.
    Real prev_time   = state[State_Type].prevTime();
    if (prev_time > 0.0 && verbose > 0)
    {
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();

            fort_compute_max_temp_loc
                (bx.loVect(), bx.hiVect(),
                BL_TO_FORTRAN(S_new[mfi]),
                BL_TO_FORTRAN(D_new[mfi]),
                &max_temp,&den_maxt,&imax,&jmax,&kmax);
        }

        if (verbose > 1 && ParallelDescriptor::IOProcessor())
            if (imax > -1 && jmax > -1 && kmax > -1)
            {
              std::cout << "Maximum temp. at level " << level << " is " << max_temp
                        << " at density " << den_maxt
                        << " at (i,j,k) " << imax << " " << jmax << " " << kmax << std::endl;
            }
    }
}
#endif

#ifndef NO_HYDRO
void
Nyx::compute_rho_temp (Real& rho_T_avg, Real& T_avg, Real& Tinv_avg, Real& T_meanrho)
{
    BL_PROFILE("Nyx::compute_rho_temp()");
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);

    Real rho_T_sum=0.0,   T_sum=0.0, Tinv_sum=0.0, T_meanrho_sum=0.0;
    Real   rho_sum=0.0, vol_sum=0.0,    vol_mn_sum=0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:rho_T_sum, rho_sum, T_sum, Tinv_sum, T_meanrho_sum, vol_sum, vol_mn_sum)
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        fort_compute_rho_temp
            (bx.loVect(), bx.hiVect(), geom.CellSize(),
             BL_TO_FORTRAN(S_new[mfi]),
             BL_TO_FORTRAN(D_new[mfi]), &average_gas_density,
             &rho_T_sum, &T_sum, &Tinv_sum, &T_meanrho_sum, &rho_sum, &vol_sum, &vol_mn_sum);
    }
    Real sums[7] = {rho_T_sum, rho_sum, T_sum, Tinv_sum, T_meanrho_sum, vol_sum, vol_mn_sum};
    ParallelDescriptor::ReduceRealSum(sums,7);

    rho_T_avg = sums[0] / sums[1];  // density weighted T
        T_avg = sums[2] / sums[5];  // volume weighted T
     Tinv_avg = sums[3] / sums[1];  // 21cm T
    if (sums[6] > 0) {
       T_meanrho = sums[4] / sums[6];  // T at mean density
       T_meanrho = pow(10.0, T_meanrho);
    }
}
#endif

#ifndef NO_HYDRO
void
Nyx::compute_gas_fractions (Real T_cut, Real rho_cut,
                            Real& whim_mass_frac, Real& whim_vol_frac,
                            Real& hh_mass_frac,   Real& hh_vol_frac,
                            Real& igm_mass_frac,  Real& igm_vol_frac)
{
    BL_PROFILE("Nyx::compute_gas_fractions()");
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);

    Real whim_mass=0.0, whim_vol=0.0, hh_mass=0.0, hh_vol=0.0, igm_mass=0.0, igm_vol=0.0;
    Real mass_sum=0.0, vol_sum=0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:whim_mass, whim_vol, hh_mass, hh_vol, igm_mass, igm_vol, mass_sum, vol_sum)
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        fort_compute_gas_frac
            (bx.loVect(), bx.hiVect(), geom.CellSize(),
             BL_TO_FORTRAN(S_new[mfi]),
             BL_TO_FORTRAN(D_new[mfi]), &average_gas_density, &T_cut, &rho_cut,
             &whim_mass, &whim_vol, &hh_mass, &hh_vol, &igm_mass, &igm_vol, &mass_sum, &vol_sum);
    }
    Real sums[8] = {whim_mass, whim_vol, hh_mass, hh_vol, igm_mass, igm_vol, mass_sum, vol_sum};
    ParallelDescriptor::ReduceRealSum(sums,8);

    whim_mass_frac = sums[0] / sums[6];
    whim_vol_frac  = sums[1] / sums[7];
    hh_mass_frac   = sums[2] / sums[6];
    hh_vol_frac    = sums[3] / sums[7];
    igm_mass_frac  = sums[4] / sums[6];
    igm_vol_frac   = sums[5] / sums[7];
}
#endif

Real
Nyx::getCPUTime()
{

  int numCores = ParallelDescriptor::NProcs();
#ifdef _OPENMP
  numCores = numCores*omp_get_max_threads();
#endif

  Real T = numCores*(ParallelDescriptor::second() - startCPUTime) +
    previousCPUTimeUsed;

  return T;
}

void
Nyx::InitErrorList() {
    //err_list.clear(true);
    //err_list.add("FULLSTATE",1,ErrorRec::Special,FORT_DENERROR);
}


//static Box the_same_box (const Box& b) { return b; }

void
Nyx::InitDeriveList() {
}


void
Nyx::LevelDirectoryNames(const std::string &dir,
                         const std::string &secondDir,
                         std::string &LevelDir,
                         std::string &FullPath)
{
    LevelDir = amrex::Concatenate("Level_", level, 1);
    //
    // Now for the full pathname of that directory.
    //
    FullPath = dir;
    if( ! FullPath.empty() && FullPath.back() != '/') {
        FullPath += '/';
    }
    FullPath += secondDir;
    FullPath += "/";
    FullPath += LevelDir;
}


void
Nyx::CreateLevelDirectory (const std::string &dir)
{
    AmrLevel::CreateLevelDirectory(dir);  // ---- this sets levelDirectoryCreated = true

    std::string dm(dir + "/" + Nyx::retrieveDM());
    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(dm, 0755)) {
        amrex::CreateDirectoryFailed(dm);
      }
    }

    std::string LevelDir, FullPath;
    LevelDirectoryNames(dir, Nyx::retrieveDM(), LevelDir, FullPath);
    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
        amrex::CreateDirectoryFailed(FullPath);
      }
    }

#ifdef AGN
    std::string agn(dir + "/" + Nyx::retrieveAGN());
    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(agn, 0755)) {
        amrex::CreateDirectoryFailed(agn);
      }
    }

    LevelDirectoryNames(dir, Nyx::retrieveAGN(), LevelDir, FullPath);
    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
        amrex::CreateDirectoryFailed(FullPath);
      }
    }
#endif

    if(parent->UsingPrecreateDirectories()) {
      if(Nyx::theDMPC()) {
        Nyx::theDMPC()->SetLevelDirectoriesCreated(true);
      }
#ifdef AGN
      if(Nyx::theAPC()) {
        Nyx::theAPC()->SetLevelDirectoriesCreated(true);
      }
#endif
    }

}


