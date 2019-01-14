
#include "myfunc.H"
#include "mykernel.H"

void advance (MultiFab& phi_old,
              MultiFab& phi_new,
	      Array<MultiFab, AMREX_SPACEDIM> const& flux,
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
}
