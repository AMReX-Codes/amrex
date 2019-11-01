#ifndef GUARDCELLMANAGER_H_
#define GUARDCELLMANAGER_H_

#include <AMReX_IntVect.H>

class guardCellManager{

public:

    void Init(
        const bool do_subcycling,
        const bool do_fdtd_nci_corr,
        const bool do_nodal,
        const bool do_moving_window,
        const bool do_fft_mpi_dec,
        const bool aux_is_nodal,
        const int moving_window_dir,
        const int nox,
        const int nox_fft, const int noy_fft, const int noz_fft,
        const int nci_corr_stencil,
        const int maxwell_fdtd_solver_id,
        const int max_level);

    // Guard cells to initialize multifabs
    amrex::IntVect ngE = amrex::IntVect::TheZeroVector();
    amrex::IntVect ngJ = amrex::IntVect::TheZeroVector();
    amrex::IntVect ngRho = amrex::IntVect::TheZeroVector();
    amrex::IntVect ngExtra = amrex::IntVect::TheZeroVector();
    amrex::IntVect ngF = amrex::IntVect::TheZeroVector();
    int ngF_int = 0;
    
    //amrex::IntVect ng_alloc_EB = amrex::IntVect::TheZeroVector();
    //amrex::IntVect ng_alloc_J = amrex::IntVect::TheZeroVector();
    //amrex::IntVect ng_alloc_Rho = amrex::IntVect::TheZeroVector();
    //amrex::IntVect ng_alloc_F = amrex::IntVect::TheZeroVector();
    //int ng_alloc_F_int = 0;

    // Guard cells to exchange data
    amrex::IntVect ng_FieldSolver = amrex::IntVect::TheZeroVector();
    amrex::IntVect ng_FieldGather = amrex::IntVect::TheZeroVector();
    amrex::IntVect ng_Aux = amrex::IntVect::TheZeroVector();

    // Extra guard cells for fine level of E and B
    amrex::IntVect ng_Extra = amrex::IntVect::TheZeroVector();
};

#endif // GUARDCELLMANAGER_H_
