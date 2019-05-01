#ifndef WARPX_SUM_GUARD_CELLS_H_
#define WARPX_SUM_GUARD_CELLS_H_

#include <AMReX_MultiFab.H>

/* \brief Sum the values of `mf`, where the different boxes overlap
 * (i.e. in the guard cells)
 *
 * This is typically called for the sources of the Maxwell equations (J/rho)
 * after deposition from the macroparticles.
 *
 *  - When WarpX is compiled with a finite-difference scheme: this only
 *    updates the *valid* cells of `mf`
 *  - When WarpX is compiled with a spectral scheme (WARPX_USE_PSATD): this
 *    updates both the *valid* cells and *guard* cells. (This is because a
 *    spectral solver requires the value of the sources over a large stencil.)
 */
void
WarpXSumGuardCells(amrex::MultiFab& mf, const amrex::Periodicity& period,
                   const int icomp=0, const int ncomp=1){
#ifdef WARPX_USE_PSATD
    // Update both valid cells and guard cells
    const int n_grow = mf.nGrow();
    // Create temporary MultiFab, for safe summation
    amrex::MultiFab tmp(mf.boxArray(), mf.DistributionMap(), ncomp, n_grow);
    amrex::MultiFab::Copy(tmp, mf, icomp, 0, ncomp, n_grow);
    mf.setVal(0.0, icomp, ncomp, 0);
    // Perform summation on overlapping grid points
    mf.copy(tmp, 0, icomp, ncomp, n_grow, n_grow, period, amrex::FabArrayBase::ADD);
#else
    // Update only the valid cells
    mf.SumBoundary(icomp, ncomp, period);
#endif
}

/* \brief Sum the values of `src` where the different boxes overlap
 * (i.e. in the guard cells) and copy them into `dst`
 *
 * This is typically called for the sources of the Maxwell equations (J/rho)
 * after deposition from the macroparticles + filtering.
 *
 *  - When WarpX is compiled with a finite-difference scheme: this only
 *    updates the *valid* cells of `dst`
 *  - When WarpX is compiled with a spectral scheme (WARPX_USE_PSATD): this
 *    updates both the *valid* cells and *guard* cells. (This is because a
 *    spectral solver requires the value of the sources over a large stencil.)
 *
 * Note: `i_comp` is the component where the results will be stored in `dst`;
 *       The component from which we copy in `src` is always 0.
 */
void
WarpXSumGuardCells(amrex::MultiFab& dst, amrex::MultiFab& src,
                   const amrex::Periodicity& period,
                   const int icomp=0, const int ncomp=1){
#ifdef WARPX_USE_PSATD
    // Update both valid cells and guard cells
    dst.ParallelCopy( src, 0, icomp, ncomp,
            src.nGrow(), dst.nGrow(), period, amrex::FabArrayBase::ADD );
#else
    // Update only the valid cells
    src.SumBoundary(icomp, ncomp, period);
    amrex::MultiFab::Copy( dst, src, 0, icomp, ncomp, 0 );
#endif
}

#endif // WARPX_SUM_GUARD_CELLS_H_
