
// #include <unistd.h>

// #include <iomanip>

// #include <algorithm>
// #include <cstdio>
// #include <vector>
// #include <iostream>
// #include <string>
// #include <ctime>

// using std::cout;
// using std::cerr;
// using std::endl;
// using std::istream;
// using std::ostream;
// using std::pair;
// using std::string;

// #include <Utility.H>
// #include <CONSTANTS.H>
#include <SMC.H>
// #include <SMC_F.H>
// #include <VisMF.H>
// #include <TagBox.H>
#include <ParmParse.H>

#ifdef _OPENMP
#include <omp.h>
#endif

static Real fixed_dt     = -1.0;
static Real initial_dt   = -1.0;
static Real dt_cutoff    = 0.0;

int  SMC::ncons       = 0;
int  SMC::nprim       = 0;
int  SMC::verbose     = 0;
int  SMC::cfl_int     = 0;
Real SMC::cfl         = 0.1;
Real SMC::init_shrink = 0.5;
Real SMC::change_max  = 1.1;

// Note: SMC::variableSetUp is in SMC_setup.cpp

void
SMC::variableCleanUp () 
{
    desc_lst.clear();
}

void
SMC::read_params ()
{
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("smc");   

    pp.query("v",verbose);
    pp.query("cfl_int",cfl_int);
    pp.query("cfl",cfl);
    pp.query("init_shrink",init_shrink);
    pp.query("change_max",change_max);

    pp.query("fixed_dt",fixed_dt);
    pp.query("initial_dt",initial_dt);
    pp.query("dt_cutoff",dt_cutoff);
}

SMC::SMC (Amr&            papa,
	  int             lev,
	  const Geometry& level_geom,
	  const BoxArray& bl,
	  Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,time) 
{
    int ngrow = 4;
    prim.define(grids,nprim,ngrow,Fab_allocate);
}

SMC::~SMC () 
{
    ;
}

void
SMC::initData ()
{
    BL_PROFILE("SMC::initData()");
}

void
SMC::init (AmrLevel &old)
{
    BL_PROFILE("SMC::init(old)");
}

//
// This version inits the data on a new level that did not
// exist before regridding.
//
void
SMC::init ()
{
    BL_PROFILE("SMC::init()");
}

Real
SMC::initialTimeStep ()
{
    Real dummy_dt = 0.0;
    Real init_dt  = 0.0;

    if (initial_dt > 0.0) 
    {
       init_dt = initial_dt;
    } 
    else 
    {
       init_dt = init_shrink*estTimeStep(dummy_dt);
    }

    return init_dt;
}

Real
SMC::estTimeStep (Real dt_old)
{
    BL_PROFILE("SMC::estTimeStep()");

    if (fixed_dt > 0.0)
        return fixed_dt;

    // This is just a dummy value to start with 
    Real estdt  = 1.0e+200;

    const MultiFab& stateMF = get_new_data(State_Type);

    // xxxxxx


    if (verbose && ParallelDescriptor::IOProcessor())
	std::cout << "SMC::estTimeStep at level " << level << ":  estdt = " << estdt << '\n';

    return estdt;
}

void
SMC::computeNewDt (int                   finest_level,
		   int                   sub_cycle,
		   Array<int>&           n_cycle,
		   const Array<IntVect>& ref_ratio,
		   Array<Real>&          dt_min,
		   Array<Real>&          dt_level,
		   Real                  stop_time,
		   int                   post_regrid_flag)
{
    BL_PROFILE("SMC::computeNewDt()");

}

void
SMC::computeInitialDt (int                   finest_level,
		       int                   sub_cycle,
		       Array<int>&           n_cycle,
		       const Array<IntVect>& ref_ratio,
		       Array<Real>&          dt_level,
		       Real                  stop_time)
{
    BL_PROFILE("SMC::computeInitialDt()");

}

void
SMC::post_timestep (int iteration)
{
    BL_PROFILE("SMC::post_timestep()");

}

void
SMC::post_restart ()
{
    BL_PROFILE("SMC::post_restart()");

}

void
SMC::post_regrid (int lbase,
		  int new_finest)
{
}

void
SMC::post_init (Real stop_time)
{
}

int
SMC::okToContinue ()
{
    if (level > 0)
        return 1;

    int test = 1;
    if (parent->dtLevel(0) < dt_cutoff) test = 0;

    return test; 
}

