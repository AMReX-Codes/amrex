
#include "myfunc.H"
#include "myfunc_F.H"
#include <AMReX_Gpu.H>

void advance (MultiFab& phi_old,
              MultiFab& phi_new,
              Array<MultiFab, AMREX_SPACEDIM>& flux,
              Real dt,
              Geometry const& geom)
{

    // Fill the ghost cells of each grid from the other grids
    // includes periodic domain boundaries
    phi_old.FillBoundary(geom.periodicity());
    //
    // Note that this simple example is not optimized.
    // The following two MFIter loops could be merged
    // and we do not have to use flux MultiFab.
    // 
    // =======================================================

    const Real dxinv = geom.InvCellSize(0);
    const Real dyinv = geom.InvCellSize(1);
#if (AMREX_SPACEDIM == 3)   
    const Real dzinv = geom.InvCellSize(2);
#endif
    int Ncomp = phi_old.nComp();
    int ng_p = phi_old.nGrow();
    int ng_f = flux[0].nGrow();

    // Compute fluxes one grid at a time


    for(MFIter mfi(phi_old, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& xbx = mfi.nodaltilebox(0);
        const Box& ybx = mfi.nodaltilebox(1);
#if (AMREX_SPACEDIM == 3)   
        const Box& zbx = mfi.nodaltilebox(2);
#endif
        FArrayBox* phi_fab  = phi_old.fabPtr(mfi);
        FArrayBox* fluxx = flux[0].fabPtr(mfi);
        FArrayBox* fluxy = flux[1].fabPtr(mfi);
#if (AMREX_SPACEDIM == 3)   
        FArrayBox* fluxz = flux[2].fabPtr(mfi);
#endif

        compute_flux_x_acc(BL_TO_FORTRAN_BOX(xbx),
                           BL_TO_FORTRAN_ANYD(*fluxx),
                           BL_TO_FORTRAN_ANYD(*phi_fab), dxinv);

        compute_flux_y_acc(BL_TO_FORTRAN_BOX(ybx),
                           BL_TO_FORTRAN_ANYD(*fluxy),
                           BL_TO_FORTRAN_ANYD(*phi_fab), dyinv);

#if (AMREX_SPACEDIM == 3)   
        compute_flux_z_acc(BL_TO_FORTRAN_BOX(zbx),
                           BL_TO_FORTRAN_ANYD(*fluxz),
                           BL_TO_FORTRAN_ANYD(*phi_fab), dzinv);
#endif
    }
    Gpu::Device::synchronize();

    for(MFIter mfi(phi_old, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const FArrayBox* pold_fab  = phi_old.fabPtr(mfi);
        const FArrayBox* fluxx = flux[0].fabPtr(mfi);
        const FArrayBox* fluxy = flux[1].fabPtr(mfi);
#if (AMREX_SPACEDIM == 3)   
        const FArrayBox* fluxz = flux[2].fabPtr(mfi);
#endif
        FArrayBox* pnew_fab = phi_new.fabPtr(mfi);

        update_phi_acc(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_ANYD(*fluxx),
                       BL_TO_FORTRAN_ANYD(*fluxy),
#if (AMREX_SPACEDIM == 3)   
                       BL_TO_FORTRAN_ANYD(*fluxz),
#endif
                       BL_TO_FORTRAN_ANYD(*pold_fab),
                       BL_TO_FORTRAN_ANYD(*pnew_fab),
                       dt,dxinv,dyinv
#if (AMREX_SPACEDIM == 3)   
                       ,dzinv
#endif
                       );
    }
    Gpu::Device::synchronize();

}
