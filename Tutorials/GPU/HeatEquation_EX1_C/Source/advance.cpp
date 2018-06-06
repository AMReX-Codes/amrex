
#include "myfunc.H"
#include "myfunc_F.H"

void advance (MultiFab& phi_old,
              MultiFab& phi_new,
	      std::array<MultiFab, AMREX_SPACEDIM>& flux,
	      Real dt,
              const Geometry& geom)
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

    Box* bx;
    Real* dx;
    Box* domain_bx;

    cudaMallocHost(&bx, size_t(sizeof(Box)));
    cudaMallocHost(&domain_bx, size_t(sizeof(Box)));
    cudaMallocHost(&dx, size_t(AMREX_SPACEDIM*sizeof(Real)));

    for (int i=0; i<AMREX_SPACEDIM; ++i)
    {
      dx[i] = geom.CellSize()[i];
    }
    *domain_bx = geom.Domain();

    cudaDeviceSynchronize();

    // Compute fluxes one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        amrex::Print() << "AMREX_SPACEDIM = " << AMREX_SPACEDIM << std::endl;
        *bx = mfi.validbox();

        compute_flux<<<1, 1>>>(BL_TO_FORTRAN_BOX(*bx),
                     BL_TO_FORTRAN_BOX(*domain_bx),
                     BL_TO_FORTRAN_ANYD(phi_old[mfi]),
                     BL_TO_FORTRAN_ANYD(flux[0][mfi]),
                     BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
                     BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif
                     dx);

    }

    cudaDeviceSynchronize();
 
    // Advance the solution one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        *bx = mfi.validbox();
        
        update_phi<<<1,1>>>
                  (BL_TO_FORTRAN_BOX(*bx),
                   BL_TO_FORTRAN_ANYD(phi_old[mfi]),
                   BL_TO_FORTRAN_ANYD(phi_new[mfi]),
                   BL_TO_FORTRAN_ANYD(flux[0][mfi]),
                   BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
                   BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif
                   dx, dt);

    }

    cudaDeviceSynchronize();
    cudaFreeHost(bx);
    cudaFreeHost(dx);
    cudaFreeHost(domain_bx);
}
