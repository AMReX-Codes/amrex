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
    int n_cell = 32;
    int max_grid_size = 32;
    int verbose = 2;

    // data
    amrex::MultiFab rhs;
    amrex::MultiFab soln;
    amrex::MultiFab the_soln;
    amrex::Geometry geom;
    
    void init_rhs ();
    void comp_the_solution ();
};

#endif
