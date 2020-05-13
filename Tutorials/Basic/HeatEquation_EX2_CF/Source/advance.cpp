
#include "myfunc.H"
#include "myfunc_F.H"

#include <AMReX_BCUtil.H>

using namespace amrex;

void advance (MultiFab& phi_old,
              MultiFab& phi_new,
	      Array<MultiFab, AMREX_SPACEDIM>& flux,
	      Real dt,
              const Geometry& geom,
              const Vector<BCRec>& bc)
{
    // Fill the ghost cells of each grid from the other grids
    // includes periodic domain boundaries
    phi_old.FillBoundary(geom.periodicity());

    // Fill non-periodic physical boundaries
    FillDomainBoundary(phi_old, geom, bc);

    int Ncomp = phi_old.nComp();
    int ng_p = phi_old.nGrow();
    int ng_f = flux[0].nGrow();

    const Real* dx = geom.CellSize();

    //
    // Note that this simple example is not optimized.
    // The following two MFIter loops could be merged
    // and we do not have to use flux MultiFab.
    // 

    const Box& domain_bx = geom.Domain();

    // Compute fluxes one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        compute_flux(BL_TO_FORTRAN_BOX(bx),
                     BL_TO_FORTRAN_BOX(domain_bx),
                     BL_TO_FORTRAN_ANYD(phi_old[mfi]),
                     BL_TO_FORTRAN_ANYD(flux[0][mfi]),
                     BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
                     BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif
                     dx, bc[0].data());
    }
    
    // Advance the solution one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();
        
        update_phi(BL_TO_FORTRAN_BOX(bx),
                   BL_TO_FORTRAN_ANYD(phi_old[mfi]),
                   BL_TO_FORTRAN_ANYD(phi_new[mfi]),
                   BL_TO_FORTRAN_ANYD(flux[0][mfi]),
                   BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
                   BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif
                   dx, &dt);
    }
}
