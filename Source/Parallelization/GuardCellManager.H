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
    amrex::IntVect ngExtra;
    amrex::IntVect ngE;
    amrex::IntVect ngJ;
    amrex::IntVect ngRho;
    amrex::IntVect ngF;
    int ngF_int;

    // Guard cells to exchange data
    amrex::IntVect ngB_FieldSolver;
    amrex::IntVect ngE_FieldSolver;
    amrex::IntVect ng_FieldGather;
    amrex::IntVect ngJ_CurrentDepo;
    amrex::IntVect ng_MovingWindow;
    amrex::IntVect ng_NCIFilter;
};

#endif // GUARDCELLMANAGER_H_
