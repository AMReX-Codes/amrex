#ifndef _SMC_H_
#define _SMC_H_

#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include "AMReX_AbstractTask.H"
#include "AMReX_TaskGraph.H"
#include "RTS.H"
#include "AMReX_AsyncMFIter.H"

#ifdef _OPEMP
#include <omp.h>
#endif

class MyAction;
namespace amrex {

class SMC
{
    public:

	SMC ();
	~SMC ();

	void evolve ();

    private:

	static bool initialized;

	// data
	amrex::Geometry geom;
	amrex::MultiFab U;
	amrex::MultiFab Utmp;
	amrex::MultiFab Uprime;
	amrex::MultiFab Q;
	amrex::MultiFab mu;
	amrex::MultiFab xi;
	amrex::MultiFab lam;
	amrex::MultiFab Ddiag;

	amrex::Real t, dt, courno;
	amrex::Real dx[3];

	// number of components
	static int ncons;
	static int nprim;
	static int nspec;
	static int nplot;
	// number of ghost cells
	static int ngrow;

	// runtime parameters;
	static amrex::Vector<int> ncell;
	static amrex::Vector<int> max_grid_size;
	static int max_step;
	static amrex::Real stop_time;
	static amrex::Vector<amrex::Real> prob_lo;
	static amrex::Vector<amrex::Real> prob_hi;
	//
	static int  verbose;
	//
	static int  cfl_int;
	static amrex::Real cfl;
	static amrex::Real init_shrink;
	static amrex::Real fixed_dt;
	//
	static int  plot_int;
	//
	// timers
	amrex::Real wt_fb1;
	amrex::Real wt_fb2;
	amrex::Real wt_chem1;
	amrex::Real wt_chem2;
	amrex::Real wt_hypdiff;

	static void init_runtime ();
	static void init_stencil ();
	static void init_chemistry ();
	static void init_variables ();

	void build_multifabs ();
	void init_from_scratch ();

	void writePlotFile (int istep);

	void advance (int istep);
	void compute_dU0(FArrayBox& UU_fab, FArrayBox& Uprime_fab, FArrayBox& Q_fab, FArrayBox& mu_fab, FArrayBox& xi_fab, FArrayBox& lam_fab, FArrayBox& Ddiag_fab, Box, Box, Box);
	void compute_dU1(FArrayBox& UU_fab, FArrayBox& Uprime_fab, FArrayBox& Q_fab, FArrayBox& mu_fab, FArrayBox& xi_fab, FArrayBox& lam_fab, FArrayBox& Ddiag_fab, Box);
	void update_rk3(amrex::Real a, amrex::FArrayBox& Ua_fab, amrex::Real b, const amrex::FArrayBox& Ub_fab, amrex::Real c, const amrex::FArrayBox& Uc_fab, int stage, Box);
	void reset_density(amrex::FArrayBox& UU_fab, Box, Box, int ngrow);
   friend MyAction;
};

}//end namespace

#endif
