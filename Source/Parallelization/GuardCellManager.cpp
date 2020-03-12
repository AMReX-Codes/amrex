/* Copyright 2019-2020 Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "GuardCellManager.H"
#include "Filter/NCIGodfreyFilter.H"
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>

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
    const int max_level,
    const amrex::Array<amrex::Real,3> v_galilean,
    const bool safe_guard_cells)
{
    // When using subcycling, the particles on the fine level perform two pushes
    // before being redistributed ; therefore, we need one extra guard cell
    // (the particles may move by 2*c*dt)
    int ngx_tmp = (max_level > 0 && do_subcycling == 1) ? nox+1 : nox;
    int ngy_tmp = (max_level > 0 && do_subcycling == 1) ? nox+1 : nox;
    int ngz_tmp = (max_level > 0 && do_subcycling == 1) ? nox+1 : nox;

    if ((v_galilean[0]!=0) or
        (v_galilean[1]!=0) or
        (v_galilean[2]!=0)){
      // Add one guard cell in the case of the galilean algorithm
      ngx_tmp += 1;
      ngy_tmp += 1;
      ngz_tmp += 1;
    }

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

#if (AMREX_SPACEDIM == 3)
    ng_alloc_EB = IntVect(ngx,ngy,ngz);
    ng_alloc_J = IntVect(ngJx,ngJy,ngJz);
#elif (AMREX_SPACEDIM == 2)
    ng_alloc_EB = IntVect(ngx,ngz);
    ng_alloc_J = IntVect(ngJx,ngJz);
#endif

    ng_alloc_Rho = ng_alloc_J+1; //One extra ghost cell, so that it's safe to deposit charge density
    // after pushing particle.
    int ng_alloc_F_int = (do_moving_window) ? 2 : 0;
    // CKC solver requires one additional guard cell
    if (maxwell_fdtd_solver_id == 1) ng_alloc_F_int = std::max( ng_alloc_F_int, 1 );
    ng_alloc_F = IntVect(AMREX_D_DECL(ng_alloc_F_int, ng_alloc_F_int, ng_alloc_F_int));

#ifdef WARPX_USE_PSATD
    if (do_fft_mpi_dec == false){
        // All boxes should have the same number of guard cells
        // (to avoid temporary parallel copies)
        // Thus take the max of the required number of guards for each field
        // Also: the number of guard cell should be enough to contain
        // the stencil of the FFT solver. Here, this number (`ngFFT`)
        // is determined *empirically* to be the order of the solver
        // for nodal, and half the order of the solver for staggered.

        int ngFFt_x = do_nodal ? nox_fft : nox_fft/2.;
        int ngFFt_y = do_nodal ? noy_fft : noy_fft/2.;
        int ngFFt_z = do_nodal ? noz_fft : noz_fft/2.;

        ParmParse pp("psatd");
        pp.query("nx_guard", ngFFt_x);
        pp.query("ny_guard", ngFFt_y);
        pp.query("nz_guard", ngFFt_z);

        IntVect ngFFT = IntVect(AMREX_D_DECL(ngFFt_x, ngFFt_y, ngFFt_z));

        for (int i_dim=0; i_dim<AMREX_SPACEDIM; i_dim++ ){
            int ng_required = ngFFT[i_dim];
            // Get the max
            ng_required = std::max( ng_required, ng_alloc_EB[i_dim] );
            ng_required = std::max( ng_required, ng_alloc_J[i_dim] );
            ng_required = std::max( ng_required, ng_alloc_Rho[i_dim] );
            ng_required = std::max( ng_required, ng_alloc_F[i_dim] );
            // Set the guard cells to this max
            ng_alloc_EB[i_dim] = ng_required;
            ng_alloc_J[i_dim] = ng_required;
            ng_alloc_F[i_dim] = ng_required;
            ng_alloc_Rho[i_dim] = ng_required;
            ng_alloc_F_int = ng_required;
        }
    }
    ng_alloc_F = IntVect(AMREX_D_DECL(ng_alloc_F_int, ng_alloc_F_int, ng_alloc_F_int));
#endif

    ng_Extra = IntVect(static_cast<int>(aux_is_nodal and !do_nodal));

    // Compute number of cells required for Field Solver
#ifdef WARPX_USE_PSATD
    ng_FieldSolver = ng_alloc_EB;
    ng_FieldSolverF = ng_alloc_EB;
#else
    ng_FieldSolver = IntVect(AMREX_D_DECL(1,1,1));
    ng_FieldSolverF = IntVect(AMREX_D_DECL(1,1,1));
#endif

    if (safe_guard_cells){
        // Run in safe mode: exchange all allocated guard cells at each
        // call of FillBoundary
        ng_FieldSolver = ng_alloc_EB;
        ng_FieldSolverF = ng_alloc_F;
        ng_FieldGather = ng_alloc_EB;
        ng_UpdateAux = ng_alloc_EB;
        if (do_moving_window){
            ng_MovingWindow = ng_alloc_EB;
        }
    } else {

        ng_FieldSolver = ng_FieldSolver.min(ng_alloc_EB);

        // Compute number of cells required for Field Gather
        int FGcell[4] = {0,1,1,2}; // Index is nox
        IntVect ng_FieldGather_noNCI = IntVect(AMREX_D_DECL(FGcell[nox],FGcell[nox],FGcell[nox]));
        // Add one cell if momentum_conserving gather in a staggered-field simulation
        ng_FieldGather_noNCI += ng_Extra;
        // Not sure why, but need one extra guard cell when using MR
        if (max_level >= 1) ng_FieldGather_noNCI += ng_Extra;
        ng_FieldGather_noNCI = ng_FieldGather_noNCI.min(ng_alloc_EB);
        // If NCI filter, add guard cells in the z direction
        IntVect ng_NCIFilter = IntVect::TheZeroVector();
        if (do_fdtd_nci_corr)
            ng_NCIFilter[AMREX_SPACEDIM-1] = NCIGodfreyFilter::m_stencil_width;
        // Note: communications of guard cells for bilinear filter are handled
        // separately.
        ng_FieldGather = ng_FieldGather_noNCI + ng_NCIFilter;

        // Guard cells for auxiliary grid.
        // Not sure why there is a 2* here...
        ng_UpdateAux = 2*ng_FieldGather_noNCI + ng_NCIFilter;

        // Make sure we do not exchange more guard cells than allocated.
        ng_FieldGather = ng_FieldGather.min(ng_alloc_EB);
        ng_UpdateAux = ng_UpdateAux.min(ng_alloc_EB);
        ng_FieldSolverF = ng_FieldSolverF.min(ng_alloc_F);
        // Only FillBoundary(ng_FieldGather) is called between consecutive
        // field solves. So ng_FieldGather must have enough cells
        // for the field solve too.
        ng_FieldGather = ng_FieldGather.max(ng_FieldSolver);

        if (do_moving_window){
            ng_MovingWindow[moving_window_dir] = 1;
        }
    }
}
