
#include "myfunc.H"
#include "myfunc_F.H"

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
/*
    Real* dx;
    Box* domain_bx;

    cudaMallocHost(&domain_bx, size_t(sizeof(Box)));
    cudaMallocHost(&dx, size_t(AMREX_SPACEDIM*sizeof(Real)));

    for (int i=0; i<AMREX_SPACEDIM; ++i)
    {
      dx[i] = geom.CellSize()[i];
    }
    *domain_bx = geom.Domain();
*/

    GeometryData* geomData = geom.dataPtr();

    // Compute fluxes one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        Box* bx;
        cudaMallocHost(&bx, size_t(sizeof(Box)));
        *bx = mfi.validbox();

        compute_flux<<<1, 1>>>(BL_TO_FORTRAN_BOX(*bx),
                     BL_TO_FORTRAN_BOX(geom.dataPtr()->Domain()),
                     BL_TO_FORTRAN_ANYD(phi_old[mfi]),
                     BL_TO_FORTRAN_ANYD(flux[0][mfi]),
                     BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
                     BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif
                     geom.dataPtr()->CellSize());

        cudaDeviceSynchronize();
        cudaFreeHost(bx);
    }

    // Advance the solution one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        Box* bx;
        cudaMallocHost(&bx, size_t(sizeof(Box)));
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
                   geom.dataPtr()->CellSize(), dt);

        cudaDeviceSynchronize();
        cudaFreeHost(bx);
    }
/*
    cudaFreeHost(dx);
    cudaFreeHost(domain_bx);
*/
}
