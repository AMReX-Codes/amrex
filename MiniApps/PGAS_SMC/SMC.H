#ifndef _SMC_H_
#define _SMC_H_

#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#ifdef _OPEMP
#include <omp.h>
#endif

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
    static int overlap;

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
    void compute_dUdt(amrex::MultiFab& UU, int istep=-1);
    void set_dt(int istep);
    void update_rk3(amrex::Real a, amrex::MultiFab& Ua, amrex::Real b, const amrex::MultiFab& Ub, amrex::Real c, const amrex::MultiFab& Uc);
    void reset_density(amrex::MultiFab& UU);
};

#endif
