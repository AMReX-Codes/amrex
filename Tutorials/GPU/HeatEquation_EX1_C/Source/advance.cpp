
#include "myfunc.H"
#include "myfunc_F.H"

AMREX_CUDA_GLOBAL
void compute_flux (Box bx, GeometryData geom, BaseFab<Real> &phi_old,
                   AMREX_D_DECL(BaseFab<Real> &fluxX, BaseFab<Real> &fluxY, BaseFab<Real> &fluxZ))
{
        compute_flux(BL_TO_FORTRAN_BOX(bx),
                     BL_TO_FORTRAN_BOX(geom.Domain()),
                     BL_TO_FORTRAN_ANYD(phi_old),
                     BL_TO_FORTRAN_ANYD(fluxX),
                     BL_TO_FORTRAN_ANYD(fluxY),
#if (AMREX_SPACEDIM == 3)   
                     BL_TO_FORTRAN_ANYD(fluxZ),
#endif
                     geom.CellSize());
}

AMREX_CUDA_GLOBAL
void update_phi (Box bx, GeometryData geom, BaseFab<Real> &phi_old, BaseFab<Real> &phi_new, 
                 AMREX_D_DECL(BaseFab<Real> &fluxX, BaseFab<Real> &fluxY, BaseFab<Real> &fluxZ), Real dt)
{
        update_phi(BL_TO_FORTRAN_BOX(bx),
                   BL_TO_FORTRAN_ANYD(phi_old),
                   BL_TO_FORTRAN_ANYD(phi_new),
                   BL_TO_FORTRAN_ANYD(fluxX),
                   BL_TO_FORTRAN_ANYD(fluxY),
#if (AMREX_SPACEDIM == 3)   
                   BL_TO_FORTRAN_ANYD(fluxZ),
#endif
                   geom.CellSize(), dt);
}



void advance (MultiFab& phi_old,
              MultiFab& phi_new,
	      std::array<MultiFab, AMREX_SPACEDIM>& flux,
	      Real dt,
              Geometry& geom)
{
    // Fill the ghost cells of each grid from the other grids
    // includes periodic domain boundaries
    phi_old.FillBoundary(geom.periodicity());

    int Ncomp = phi_old.nComp();
    int ng_p = phi_old.nGrow();
    int ng_f = flux[0].nGrow();

    //
    // Note that this simple example is not optimized.
    // The following two MFIter loops could be merged
    // and we do not have to use flux MultiFab.
    // 

    // =======================================================

    // Compute fluxes one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        AMREX_SIMPLE_LAUNCH(compute_flux, 1, 1, mfi.validbox(), 
                                  geom.data(), phi_old[mfi],
                                  AMREX_D_DECL(flux[0][mfi], flux[1][mfi], flux[2][mfi]));
    }
    syncDevice();

    // Advance the solution one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        AMREX_SIMPLE_LAUNCH(update_phi, 1, 1, mfi.validbox(),
                            geom.data(), phi_old[mfi], phi_new[mfi],
                            AMREX_D_DECL(flux[0][mfi], flux[1][mfi], flux[2][mfi]),
                            dt);
    }
    syncDevice();

}
