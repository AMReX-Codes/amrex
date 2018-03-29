#ifndef MY_SWS_H_
#define MY_SWS_H_

#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

// We solve Lap(soln) = rhs
class SWFFT_Test
{
public:
    SWFFT_Test ();

    void solve ();

private:

    // Runtime parameters
    amrex::IntVect max_grid_size;   // Maximum allowable grid size in each coord direction
    amrex::IntVect n_cell;          // Number of cells             in each coord direction

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
