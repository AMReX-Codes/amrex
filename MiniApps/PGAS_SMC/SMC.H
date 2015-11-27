#ifndef _SMC_H_
#define _SMC_H_

#include <Array.H>
#include <IntVect.H>
#include <Geometry.H>
#include <MultiFab.H>

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
    Geometry geom;
    MultiFab U;
    MultiFab Utmp;
    MultiFab Uprime;
    MultiFab Q;
    MultiFab mu;
    MultiFab xi;
    MultiFab lam;
    MultiFab Ddiag;

    Real t, dt, courno;
    Real dx[3];

    // number of components
    static int ncons;
    static int nprim;
    static int nspec;
    static int nplot;
    // number of ghost cells
    static int ngrow;

    // runtime parameters;
    static Array<int> ncell;
    static Array<int> max_grid_size;
    static int max_step;
    static Real stop_time;
    static Array<Real> prob_lo;
    static Array<Real> prob_hi;
    //
    static int  verbose;
    //
    static int  cfl_int;
    static Real cfl;
    static Real init_shrink;
    static Real fixed_dt;
    //
    static int  plot_int;
    //
    static int overlap;

    // timers
    Real wt_fb1;
    Real wt_fb2;
    Real wt_chem1;
    Real wt_chem2;
    Real wt_hypdiff;

    static void init_runtime ();
    static void init_stencil ();
    static void init_chemistry ();
    static void init_variables ();

    void build_multifabs ();
    void init_from_scratch ();

    void writePlotFile (int istep);

    void advance (int istep);
    void compute_dUdt(MultiFab& UU, int istep=-1);
    void set_dt(int istep);
    void update_rk3(Real a, MultiFab& Ua, Real b, const MultiFab& Ub, Real c, const MultiFab& Uc);
    void reset_density(MultiFab& UU);
};

#endif
