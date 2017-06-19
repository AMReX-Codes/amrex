
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

    const Real* dx = geom.CellSize();

    //
    // Note that this simple example is not optimized.
    // The following two MFIter loops could be merged
    // and we do not have to use flux MultiFab.
    // 

    // Compute fluxes one grid at a time
    for ( MFIter mfi(old_phi); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        compute_flux(BL_TO_FORTRAN_BOX(bx),
                     BL_TO_FORTRAN_ANYD(old_phi[mfi]),
                     BL_TO_FORTRAN_ANYD(flux[0][mfi]),
                     BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
                     BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif
                     dx);
    }
    
    // Advance the solution one grid at a time
    for ( MFIter mfi(old_phi); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();
        
        update_phi(BL_TO_FORTRAN_BOX(bx),
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
