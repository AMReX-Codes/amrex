
#include "myfunc.H"
#include "mykernel.H"

void advance (MultiFab& phi_old,
              MultiFab& phi_new,
              Array<MultiFab, AMREX_SPACEDIM> & flux,
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

#if !defined (AMREX_USE_ACC) && !defined (AMREX_OMP_OFFLOAD)
    // Compute fluxes one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& xbx = mfi.nodaltilebox(0);
        const Box& ybx = mfi.nodaltilebox(1);
        const Box& zbx = mfi.nodaltilebox(2);
        Array4<Real> fluxx = flux[0].array(mfi);
        Array4<Real> fluxy = flux[1].array(mfi);
        Array4<Real> fluxz = flux[2].array(mfi);
        const Array4<Real> phi = phi_old.array(mfi);

        AMREX_FOR_3D ( xbx, i, j, k,
        {
            compute_flux_x(i,j,k,fluxx,phi,dxinv);
        });

        AMREX_FOR_3D ( ybx, i, j, k,
        {
            compute_flux_y(i,j,k,fluxy,phi,dyinv);
        });

        AMREX_FOR_3D ( zbx, i, j, k,
        {
            compute_flux_z(i,j,k,fluxz,phi,dzinv);
        });
    }

    // Advance the solution one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
        const Array4<Real> fluxx = flux[0].array(mfi);
        const Array4<Real> fluxy = flux[1].array(mfi);
        const Array4<Real> fluxz = flux[2].array(mfi);
        const Array4<Real> phiOld = phi_old.array(mfi);
        Array4<Real>       phiNew = phi_old.array(mfi);

        AMREX_FOR_3D ( vbx, i, j, k,
        {
            update_phi(i,j,k,phiOld,phiNew,fluxx,fluxy,fluxz,dt,dxinv,dyinv,dzinv);
        });
    }

#elif defined (AMREX_USE_ACC) && !defined (AMREX_OMP_OFFLOAD)

    for(MFIter mfi(phi_old, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& xbx = mfi.nodaltilebox(0);
        const Box& ybx = mfi.nodaltilebox(1);
        const Box& zbx = mfi.nodaltilebox(2);
        FArrayBox* phi_fab  = phi_old.fabPtr(mfi);
        FArrayBox* fluxx = flux[0].fabPtr(mfi);
        FArrayBox* fluxy = flux[1].fabPtr(mfi);
        FArrayBox* fluxz = flux[2].fabPtr(mfi);

        compute_flux_x_acc(BL_TO_FORTRAN_BOX(xbx),
                           BL_TO_FORTRAN_ANYD(*fluxx),
                           BL_TO_FORTRAN_ANYD(*phi_fab), dxinv);

        compute_flux_y_acc(BL_TO_FORTRAN_BOX(ybx),
                           BL_TO_FORTRAN_ANYD(*fluxy),
                           BL_TO_FORTRAN_ANYD(*phi_fab), dyinv);

        compute_flux_z_acc(BL_TO_FORTRAN_BOX(zbx),
                           BL_TO_FORTRAN_ANYD(*fluxz),
                           BL_TO_FORTRAN_ANYD(*phi_fab), dzinv);
    }

    for(MFIter mfi(phi_old, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const FArrayBox* pold_fab  = phi_old.fabPtr(mfi);
        const FArrayBox* fluxx = flux[0].fabPtr(mfi);
        const FArrayBox* fluxy = flux[1].fabPtr(mfi);
        const FArrayBox* fluxz = flux[2].fabPtr(mfi);
        FArrayBox* pnew_fab = phi_new.fabPtr(mfi);

        update_phi_acc(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_ANYD(*fluxx),
                       BL_TO_FORTRAN_ANYD(*fluxy),
                       BL_TO_FORTRAN_ANYD(*fluxz),
                       BL_TO_FORTRAN_ANYD(*pold_fab),
                       BL_TO_FORTRAN_ANYD(*pnew_fab),
                       dt,dxinv,dyinv,dzinv);
    }

#elif !defined (AMREX_USE_ACC) && defined (AMREX_OMP_OFFLOAD)

    for(MFIter mfi(phi_old, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& xbx = mfi.nodaltilebox(0);
        const Box& ybx = mfi.nodaltilebox(1);
        const Box& zbx = mfi.nodaltilebox(2);
        FArrayBox* phi_fab  = phi_old.fabPtr(mfi);
        FArrayBox* fluxx = flux[0].fabPtr(mfi);
        FArrayBox* fluxy = flux[1].fabPtr(mfi);
        FArrayBox* fluxz = flux[2].fabPtr(mfi);

        compute_flux_x_omp(BL_TO_FORTRAN_BOX(xbx),
                           BL_TO_FORTRAN_ANYD(*fluxx),
                           BL_TO_FORTRAN_ANYD(*phi_fab), dxinv);

        compute_flux_y_omp(BL_TO_FORTRAN_BOX(ybx),
                           BL_TO_FORTRAN_ANYD(*fluxy),
                           BL_TO_FORTRAN_ANYD(*phi_fab), dyinv);

        compute_flux_z_omp(BL_TO_FORTRAN_BOX(zbx),
                           BL_TO_FORTRAN_ANYD(*fluxz),
                           BL_TO_FORTRAN_ANYD(*phi_fab), dzinv);
    }

    for(MFIter mfi(phi_old, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const FArrayBox* pold_fab  = phi_old.fabPtr(mfi);
        const FArrayBox* fluxx = flux[0].fabPtr(mfi);
        const FArrayBox* fluxy = flux[1].fabPtr(mfi);
        const FArrayBox* fluxz = flux[2].fabPtr(mfi);
        FArrayBox* pnew_fab = phi_new.fabPtr(mfi);

        update_phi_omp(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_ANYD(*fluxx),
                       BL_TO_FORTRAN_ANYD(*fluxy),
                       BL_TO_FORTRAN_ANYD(*fluxz),
                       BL_TO_FORTRAN_ANYD(*pold_fab),
                       BL_TO_FORTRAN_ANYD(*pnew_fab),
                       dt,dxinv,dyinv,dzinv);
    }

#endif
}
