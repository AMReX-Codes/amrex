
#include "myfunc.H"
#include "myfunc_F.H"

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
    const Real dzinv = geom.InvCellSize(2);
    int Ncomp = phi_old.nComp();
    int ng_p = phi_old.nGrow();
    int ng_f = flux[0].nGrow();

    // Compute fluxes one grid at a time


    for(MFIter mfi(phi_old, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& xbx = mfi.nodaltilebox(0);
        const Box& ybx = mfi.nodaltilebox(1);
        const Box& zbx = mfi.nodaltilebox(2);
        FArrayBox* phi_fab  = phi_old.fabPtr(mfi);
        FArrayBox* fluxx = flux[0].fabPtr(mfi);
        FArrayBox* fluxy = flux[1].fabPtr(mfi);
        FArrayBox* fluxz = flux[2].fabPtr(mfi);

        compute_flux_x(BL_TO_FORTRAN_BOX(xbx),
                           BL_TO_FORTRAN_ANYD(*fluxx),
                           BL_TO_FORTRAN_ANYD(*phi_fab), dxinv);

        compute_flux_y(BL_TO_FORTRAN_BOX(ybx),
                           BL_TO_FORTRAN_ANYD(*fluxy),
                           BL_TO_FORTRAN_ANYD(*phi_fab), dyinv);

        compute_flux_z(BL_TO_FORTRAN_BOX(zbx),
                           BL_TO_FORTRAN_ANYD(*fluxz),
                           BL_TO_FORTRAN_ANYD(*phi_fab), dzinv);
    }
    Gpu::Device::synchronize();

    for(MFIter mfi(phi_old, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const FArrayBox* pold_fab  = phi_old.fabPtr(mfi);
        const FArrayBox* fluxx = flux[0].fabPtr(mfi);
        const FArrayBox* fluxy = flux[1].fabPtr(mfi);
        const FArrayBox* fluxz = flux[2].fabPtr(mfi);
        FArrayBox* pnew_fab = phi_new.fabPtr(mfi);

        update_phi(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_ANYD(*fluxx),
                       BL_TO_FORTRAN_ANYD(*fluxy),
                       BL_TO_FORTRAN_ANYD(*fluxz),
                       BL_TO_FORTRAN_ANYD(*pold_fab),
                       BL_TO_FORTRAN_ANYD(*pnew_fab),
                       dt,dxinv,dyinv ,dzinv
                       );
    }
    Gpu::Device::synchronize();
}

void init_phi(MultiFab& phi_new, Geometry const& geom)
{
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
    GpuArray<Real,AMREX_SPACEDIM> prob_lo = geom.ProbLoArray();
    for(MFIter mfi(phi_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx   = mfi.tilebox();
        FArrayBox& fab  = phi_new[mfi];
        init_phi(BL_TO_FORTRAN_BOX(bx),
                     BL_TO_FORTRAN_ANYD(fab),
                     dx.data(),prob_lo.data());
    }
}
