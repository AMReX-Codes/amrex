#include "myfunc.H"
#include "mykernel.H"
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>
#include <AMReX_BC_TYPES.H>

void advance (MultiFab& phi_old,
              MultiFab& phi_new,
	      Array<MultiFab, AMREX_SPACEDIM>& flux,
	      Real dt,
              Geometry const& geom,
	      Vector<BCRec> const& BoundaryCondition)
{
    
    // Fill the ghost cells of each grid from the other grids
    // includes periodic domain boundaries
    phi_old.FillBoundary(geom.periodicity());

    FillDomainBoundary(phi_old, geom, BoundaryCondition);

    const BCRec& bc = BoundaryCondition[0];
    
    //
    // Note that this simple example is not optimized.
    // The following two MFIter loops could be merged
    // and we do not have to use flux MultiFab.
    // 
    // =======================================================

    const Real dxinv = geom.InvCellSize(0);
    const Real dyinv = geom.InvCellSize(1);
    const Real dzinv = geom.InvCellSize(2);

    const Box& domain_bx = geom.Domain();
    const Dim3 dom_lo = lbound(domain_bx);
    const Dim3 dom_hi = ubound(domain_bx);

    // Compute fluxes one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& xbx = mfi.nodaltilebox(0);
	auto const& fluxx = flux[0].array(mfi);
        const Box& ybx = mfi.nodaltilebox(1);
	auto const& fluxy = flux[1].array(mfi);
#if (AMREX_SPACEDIM > 2)
        const Box& zbx = mfi.nodaltilebox(2);
	auto const& fluxz = flux[2].array(mfi);
#endif
	const Box& bx = mfi.validbox();
	const Dim3 lo = lbound(bx);
	const Dim3 hi = ubound(bx);
	
        auto const& phi = phi_old.array(mfi);

	amrex::ParallelFor(xbx,
	    [=] AMREX_GPU_DEVICE (int i, int j, int k)
	    {
		compute_flux_x(i,j,k,fluxx,phi,dxinv,
			       lo.x, hi.x, dom_lo.x, dom_hi.x, bc.lo(0), bc.hi(0));
	    });
	
	amrex::ParallelFor(ybx,
	    [=] AMREX_GPU_DEVICE (int i, int j, int k)
	    {
		compute_flux_y(i,j,k,fluxy,phi,dyinv,
			       lo.y, hi.y, dom_lo.y, dom_hi.y, bc.lo(1), bc.hi(1));
	    });
#if (AMREX_SPACEDIM > 2)
	amrex::ParallelFor(zbx,
	    [=] AMREX_GPU_DEVICE (int i, int j, int k)
	    {
		compute_flux_z(i,j,k,fluxz,phi,dzinv,
			       lo.z, hi.z, dom_lo.z, dom_hi.z, bc.lo(2), bc.hi(2));
	    });
#endif
    }

    // Advance the solution one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
        auto const& fluxx = flux[0].array(mfi);
        auto const& fluxy = flux[1].array(mfi);
        auto const& fluxz = flux[2].array(mfi);
        auto const& phiOld = phi_old.array(mfi);
        auto const& phiNew = phi_new.array(mfi);

        amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            update_phi(i,j,k,phiOld,phiNew,fluxx,fluxy,fluxz,dt,dxinv,dyinv,dzinv);
        });
    }
}

void init_phi(amrex::MultiFab& phi_new, amrex::Geometry const& geom){

    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
    GpuArray<Real,AMREX_SPACEDIM> prob_lo = geom.ProbLoArray();
    // =======================================
    // Initialize phi_new by calling a Fortran routine.
    // MFIter = MultiFab Iterator
    for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phiNew = phi_new.array(mfi);
        AMREX_FOR_3D ( vbx, i, j, k,
        {
            init_phi(i,j,k,phiNew,dx,prob_lo);
        });
    }
}
