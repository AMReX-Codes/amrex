#ifndef MY_SWS_H_
#define MY_SWS_H_

#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

// We solve Lap(phi_dft) = phi_spatial
class SWFFT_Test
{
public:
    SWFFT_Test ();

    void computeFFT ();

private:

    // Runtime parameters
    amrex::IntVect max_grid_size;   // Maximum allowable grid size in each coord direction
    amrex::IntVect n_cell;          // Number of cells             in each coord direction

    int verbose = 2;
    int prob_type = 0;

    // data
    amrex::MultiFab phi_spatial;
    amrex::MultiFab phi_dft;
    amrex::Geometry geom;
    
    void init_phi_spatial ();
    void WritePlotFile(const int step=0, const amrex::Real time=0.0);
};

#endif
