#include "GuardCellManager.H"

using namespace amrex;

void
guardCellManager::Init(
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
    const int max_level)
{
    // When using subcycling, the particles on the fine level perform two pushes
    // before being redistributed ; therefore, we need one extra guard cell
    // (the particles may move by 2*c*dt)
    const int ngx_tmp = (max_level > 0 && do_subcycling == 1) ? nox+1 : nox;
    const int ngy_tmp = (max_level > 0 && do_subcycling == 1) ? nox+1 : nox;
    const int ngz_tmp = (max_level > 0 && do_subcycling == 1) ? nox+1 : nox;

    // Ex, Ey, Ez, Bx, By, and Bz have the same number of ghost cells.
    // jx, jy, jz and rho have the same number of ghost cells.
    // E and B have the same number of ghost cells as j and rho if NCI filter is not used,
    // but different number of ghost cells in z-direction if NCI filter is used.
    // The number of cells should be even, in order to easily perform the
    // interpolation from coarse grid to fine grid.
    int ngx = (ngx_tmp % 2) ? ngx_tmp+1 : ngx_tmp;  // Always even number
    int ngy = (ngy_tmp % 2) ? ngy_tmp+1 : ngy_tmp;  // Always even number
    int ngz_nonci = (ngz_tmp % 2) ? ngz_tmp+1 : ngz_tmp;  // Always even number
    int ngz;
    if (do_fdtd_nci_corr) {
        int ng = ngz_tmp + nci_corr_stencil;
        ngz = (ng % 2) ? ng+1 : ng;
    } else {
        ngz = ngz_nonci;
    }

    // J is only interpolated from fine to coarse (not coarse to fine)
    // and therefore does not need to be even.
    int ngJx = ngx_tmp;
    int ngJy = ngy_tmp;
    int ngJz = ngz_tmp;

    // When calling the moving window (with one level of refinement),  we shift
    // the fine grid by 2 cells ; therefore, we need at least 2 guard cells
    // on level 1. This may not be necessary for level 0.
    if (do_moving_window) {
        ngx = std::max(ngx,2);
        ngy = std::max(ngy,2);
        ngz = std::max(ngz,2);
        ngJx = std::max(ngJx,2);
        ngJy = std::max(ngJy,2);
        ngJz = std::max(ngJz,2);
    }

/*
#if (AMREX_SPACEDIM == 3)
    IntVect ngE(ngx,ngy,ngz);
    IntVect ngJ(ngJx,ngJy,ngJz);
#elif (AMREX_SPACEDIM == 2)
    IntVect ngE(ngx,ngz);
    IntVect ngJ(ngJx,ngJz);
#endif

    IntVect ngRho = ngJ+1; //One extra ghost cell, so that it's safe to deposit charge density
    // after pushing particle.
    int ngF = (do_moving_window) ? 2 : 0;
*/

#if (AMREX_SPACEDIM == 3)
    ngE = IntVect(ngx,ngy,ngz);
    ngJ = IntVect(ngJx,ngJy,ngJz);
#elif (AMREX_SPACEDIM == 2)
    ngE = IntVect(ngx,ngz);
    ngJ = IntVect(ngJx,ngJz);
#endif

    ngRho = ngJ+1; //One extra ghost cell, so that it's safe to deposit charge density
    // after pushing particle.
    ngF_int = (do_moving_window) ? 2 : 0;
    // CKC solver requires one additional guard cell
    if (maxwell_fdtd_solver_id == 1) ngF_int = std::max( ngF_int, 1 );
    ngF = IntVect(AMREX_D_DECL(ngF_int, ngF_int, ngF_int));

#ifdef WARPX_USE_PSATD
    if (do_fft_mpi_dec == false){
        // All boxes should have the same number of guard cells
        // (to avoid temporary parallel copies)
        // Thus take the max of the required number of guards for each field
        // Also: the number of guard cell should be enough to contain
        // the stencil of the FFT solver. Here, this number (`ngFFT`)
        // is determined *empirically* to be the order of the solver
        // for nodal, and half the order of the solver for staggered.
        IntVect ngFFT;
        if (do_nodal) {
            ngFFT = IntVect(AMREX_D_DECL(nox_fft, noy_fft, noz_fft));
        } else {
            ngFFT = IntVect(AMREX_D_DECL(nox_fft/2, noy_fft/2, noz_fft/2));
        }
        for (int i_dim=0; i_dim<AMREX_SPACEDIM; i_dim++ ){
            int ng_required = ngFFT[i_dim];
            // Get the max
            ng_required = std::max( ng_required, ngE[i_dim] );
            ng_required = std::max( ng_required, ngJ[i_dim] );
            ng_required = std::max( ng_required, ngRho[i_dim] );
            ng_required = std::max( ng_required, ngF[i_dim] );
            // Set the guard cells to this max
            ngE[i_dim] = ng_required;
            ngJ[i_dim] = ng_required;
            ngF[i_dim] = ng_required;
            ngRho[i_dim] = ng_required;
            ngF_int = ng_required;
        }
    }
    ngF = IntVect(AMREX_D_DECL(ngF_int, ngF_int, ngF_int));
#endif        

    ngExtra = IntVect(static_cast<int>(aux_is_nodal and !do_nodal));

    // Guard cells for field solver
    ngB_FieldSolver = IntVect(AMREX_D_DECL(1,1,1));
    ngE_FieldSolver = IntVect(AMREX_D_DECL(1,1,1));
    ng_MovingWindow = IntVect(AMREX_D_DECL(0,0,0)); // Multiplied by number of cells moved at each timestep
    ng_MovingWindow[moving_window_dir] = 1;
    int FGcell[4] = {0,1,1,2}; // Index is nox
    ng_FieldGather = IntVect(AMREX_D_DECL(FGcell[nox],FGcell[nox],FGcell[nox]));
    ngJ_CurrentDepo = ng_FieldGather;
    ng_NCIFilter = IntVect(AMREX_D_DECL(0,0,4));
}
