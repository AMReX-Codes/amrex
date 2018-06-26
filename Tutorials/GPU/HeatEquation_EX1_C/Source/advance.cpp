
#include "myfunc.H"
#include "myfunc_F.H"
#include "AMReX_Utility.H"

AMREX_CUDA_DEVICE
Box getCellBox(const Box & bx)
{
// Return intersection of the cell for this thread and the entire domain.
// If more threads are assigned than mesh cells in the domain, intersection will return an empty box.
// If box is empty, skip the work in the MFIter loop for that thread.
// If no CUDA, return the entire box.

#ifdef AMREX_USE_CUDA
     return (bx & Box(IntVect(AMREX_D_DECL(bx.smallEnd()[0] + (threadIdx.x) + blockDim.x*(blockIdx.x),
                                           bx.smallEnd()[1] + (threadIdx.y) + blockDim.y*(blockIdx.y),
                                           bx.smallEnd()[2] + (threadIdx.z) + blockDim.z*(blockIdx.z))),
                      IntVect(AMREX_D_DECL(bx.smallEnd()[0] + (threadIdx.x) + blockDim.x*(blockIdx.x),
                                           bx.smallEnd()[1] + (threadIdx.y) + blockDim.y*(blockIdx.y),
                                           bx.smallEnd()[2] + (threadIdx.z) + blockDim.z*(blockIdx.z)))));
#else
     return bx;
#endif
}

AMREX_CUDA_GLOBAL
void compute_flux (Box bx, GeometryData geom, BaseFab<Real> &phi_old,
                   AMREX_D_DECL(BaseFab<Real> &fluxX, BaseFab<Real> &fluxY, BaseFab<Real> &fluxZ))
{
     Box threadBox = getCellBox(bx);

     if (threadBox.ok())
     {

       compute_flux(BL_TO_FORTRAN_BOX(threadBox),
                    BL_TO_FORTRAN_BOX(geom.Domain()),
                    BL_TO_FORTRAN_ANYD(phi_old),
                    BL_TO_FORTRAN_ANYD(fluxX),
                    BL_TO_FORTRAN_ANYD(fluxY),
#if (AMREX_SPACEDIM == 3)   
                    BL_TO_FORTRAN_ANYD(fluxZ),
#endif
                    geom.CellSize());
     }
}

AMREX_CUDA_GLOBAL
void update_phi (Box bx, GeometryData geom, BaseFab<Real> &phi_old, BaseFab<Real> &phi_new, 
                 AMREX_D_DECL(BaseFab<Real> &fluxX, BaseFab<Real> &fluxY, BaseFab<Real> &fluxZ), Real dt)
{
     Box threadBox = getCellBox(bx);

     if (threadBox.ok())
     {
       update_phi(BL_TO_FORTRAN_BOX(threadBox),
                  BL_TO_FORTRAN_ANYD(phi_old),
                  BL_TO_FORTRAN_ANYD(phi_new),
                  BL_TO_FORTRAN_ANYD(fluxX),
                  BL_TO_FORTRAN_ANYD(fluxY),
#if (AMREX_SPACEDIM == 3)   
                  BL_TO_FORTRAN_ANYD(fluxZ),
#endif
                  geom.CellSize(), dt);
     }
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
        const Box& vbx = mfi.validbox();

        AMREX_BOX_LAUNCH(vbx,
                          compute_flux, vbx, 
                          geom.data(), phi_old[mfi],
                          AMREX_D_DECL(flux[0][mfi], flux[1][mfi], flux[2][mfi]));
    }
    syncDevice();

    // Advance the solution one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();

        AMREX_BOX_LAUNCH(vbx,
                          update_phi, vbx,
                          geom.data(), phi_old[mfi], phi_new[mfi],
                          AMREX_D_DECL(flux[0][mfi], flux[1][mfi], flux[2][mfi]),
                          dt);
    }
    syncDevice(); 

}
