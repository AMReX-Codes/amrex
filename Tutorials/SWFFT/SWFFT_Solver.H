#ifndef MY_SWS_H_
#define MY_SWS_H_

#include <AMReX_MultiFab.H>

// We solve Lap(soln) = rhs
class SWFFT_Solver
{
public:
    SWFFT_Solver ();

    void solve ();

private:

    // runtime parameters
    int n_cell = 128;
    int max_grid_size = 64;
    int verbose = 2;

    // parameters in setting up the problem
    const amrex::Real a = 1.e-3;
    const amrex::Real b = 1.0;
    const amrex::Real sigma = 10.0;
    const amrex::Real w = 0.05;

    // data
    amrex::MultiFab rhs;
    amrex::MultiFab alpha;
    std::array<amrex::MultiFab,3> beta;
    amrex::MultiFab soln;
    amrex::MultiFab the_soln;
    amrex::Geometry geom;
    
    void init_rhs ();
    void comp_the_solution ();
};

#endif
