
#include "myfunc.H"
#include "myfunc_F.H"
#include <AMReX_CUDA_Utility.H>
#include <AMReX_Managed.H>
#include <AMReX_Device.H>

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
	const GeometryData& geomdata = geom.data();
	FArrayBox* phiOld = &(phi_old[mfi]);
	FArrayBox* fluxX = &(flux[0][mfi]);
	FArrayBox* fluxY = &(flux[1][mfi]);
	FArrayBox* fluxZ = &(flux[2][mfi]);

        AMREX_BOX_L_LAUNCH(vbx,
	[=] AMREX_CUDA_DEVICE ()
	{
             Box threadBox = getThreadBox(vbx);
             if (threadBox.ok())
             {
                compute_flux(BL_TO_FORTRAN_BOX(threadBox),
                             BL_TO_FORTRAN_BOX(geomdata.Domain()),
                             BL_TO_FORTRAN_ANYD(*phiOld),
                             BL_TO_FORTRAN_ANYD(*fluxX),
                             BL_TO_FORTRAN_ANYD(*fluxY),
#if (AMREX_SPACEDIM == 3)   
                             BL_TO_FORTRAN_ANYD(*fluxZ),
#endif
                             geomdata.CellSize());
             }
	});

    }
    Device::synchronize();

    // Advance the solution one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
	const GeometryData& geomdata = geom.data();
	FArrayBox* phiOld = &(phi_old[mfi]);
	FArrayBox* phiNew = &(phi_new[mfi]);
	FArrayBox* fluxX = &(flux[0][mfi]);
	FArrayBox* fluxY = &(flux[1][mfi]);
	FArrayBox* fluxZ = &(flux[2][mfi]);

        AMREX_BOX_L_LAUNCH(vbx, 
	[=] AMREX_CUDA_DEVICE ()
	{
            Box threadBox = getThreadBox(vbx);

            if (threadBox.ok())
            {
                update_phi(BL_TO_FORTRAN_BOX(threadBox),
                           BL_TO_FORTRAN_ANYD(*phiOld),
                           BL_TO_FORTRAN_ANYD(*phiNew),
                           BL_TO_FORTRAN_ANYD(*fluxX),
                           BL_TO_FORTRAN_ANYD(*fluxY),
#if (AMREX_SPACEDIM == 3)   
                           BL_TO_FORTRAN_ANYD(*fluxZ),
#endif
                           geomdata.CellSize(), dt);
            }
	});

    }
    Device::synchronize(); 

}
