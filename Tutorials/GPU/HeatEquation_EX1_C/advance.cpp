
#include "myfunc.H"
#include "myfunc_F.H"

using namespace amrex;

void advance (MultiFab& old_phi, MultiFab& new_phi,
	      std::array<MultiFab, AMREX_SPACEDIM>& flux,
	      Real dt, const Geometry& geom)
{
    // Fill the ghost cells of each grid from the other grids
    // includes periodic domain boundaries
    old_phi.FillBoundary(geom.periodicity());

    // Fill non-periodic physical boundaries
    fill_physbc(old_phi, geom);

    int Ncomp = old_phi.nComp();
    int ng_p = old_phi.nGrow();
    int ng_f = flux[0].nGrow();

    const Real* dx = geom.CellSizeF();

    //
    // Note that this simple example is not optimized.
    // The following two MFIter loops could be merged
    // and we do not have to use flux MultiFab.
    // 

    // Compute fluxes one grid at a time
    for ( MFIter mfi(old_phi); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        for (int idir = 1; idir <= AMREX_SPACEDIM; ++idir) {

            const Box& sbx = mfi.nodaltilebox(idir-1);

            FORT_LAUNCH(sbx, compute_flux,
                        BL_TO_FORTRAN_BOX(sbx),
                        BL_TO_FORTRAN_ANYD(old_phi[mfi]),
                        BL_TO_FORTRAN_ANYD(flux[idir-1][mfi]),
                        dx, idir);

        }

    }

    // Advance the solution one grid at a time
    for ( MFIter mfi(old_phi); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        FORT_LAUNCH(bx, update_phi,
                    BL_TO_FORTRAN_BOX(bx),
                    BL_TO_FORTRAN_ANYD(old_phi[mfi]),
                    BL_TO_FORTRAN_ANYD(new_phi[mfi]),
                    BL_TO_FORTRAN_ANYD(flux[0][mfi]),
                    BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
                    BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif
                    dx, dt);
    }

}
