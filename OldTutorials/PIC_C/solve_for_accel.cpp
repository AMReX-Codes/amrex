#include <iostream>
#include <memory>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BLFort.H>
#include <AMReX_MacBndry.H>
#include <AMReX_VisMF.H>

using namespace amrex;

void solve_with_mlmg  (const Vector<MultiFab*>& rhs, 
                       const Vector<MultiFab*>& phi, 
                       const Vector<Vector<MultiFab*> >& grad_phi_edge, 
                       const Vector<Geometry>& geom, int base_level, int finest_level,
                       Real tol, Real abs_tol);

void 
solve_for_accel(const Vector<MultiFab*>& rhs,
		const Vector<MultiFab*>& phi,
		const Vector<MultiFab*>& grad_phi, 
		const Vector<Geometry>& geom, int base_level, int finest_level, Real offset)
{
    Real tol     = 1.e-10;
    Real abs_tol = 1.e-14;

    Vector< Vector<std::unique_ptr<MultiFab> > > grad_phi_edge(rhs.size());

    for (int lev = base_level; lev <= finest_level ; lev++)
    {
	const DistributionMapping& dm = rhs[lev]->DistributionMap();
        grad_phi_edge[lev].resize(BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; ++n) {
	    BoxArray ba = rhs[lev]->boxArray();
	    grad_phi_edge[lev][n].reset(new MultiFab(ba.surroundingNodes(n), dm, 1, 1));
	}
    }

    Real     strt    = ParallelDescriptor::second();

    // ***************************************************
    // Make sure the RHS sums to 0 if fully periodic
    // ***************************************************
    for (int lev = base_level; lev <= finest_level; lev++) {
	Real n0 = rhs[lev]->norm0();
//	if (ParallelDescriptor::IOProcessor())
//	    std::cout << "Max of rhs in solve_for_phi before correction at level  " 
//                    << lev << " " << n0 << std::endl;
    }

    for (int lev = base_level; lev <= finest_level; lev++)
        rhs[lev]->plus(-offset, 0, 1, 0);

    for (int lev = base_level; lev <= finest_level; lev++) {
	Real n0 = rhs[lev]->norm0();
//	if (ParallelDescriptor::IOProcessor())
//	    std::cout << "Max of rhs in solve_for_phi  after correction at level  " 
//                    << lev << " " << n0 << std::endl;
    }

    // ***************************************************
    // Solve for phi and return both phi and grad_phi_edge
    // ***************************************************

    solve_with_mlmg  (rhs,phi,amrex::GetVecOfVecOfPtrs(grad_phi_edge),
                      geom,base_level,finest_level,tol,abs_tol);

    // Average edge-centered gradients to cell centers and fill the values in ghost cells.
    for (int lev = base_level; lev <= finest_level; lev++)
    {
        amrex::average_face_to_cellcenter(*grad_phi[lev],
					   amrex::GetVecOfConstPtrs(grad_phi_edge[lev]),
					   geom[lev]);
	grad_phi[lev]->FillBoundary(0,BL_SPACEDIM,geom[lev].periodicity());
    }

    // VisMF::Write(grad_phi[0],"GradPhi");

    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;
#if 0
#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "solve_for_phi() time = " << end << std::endl;
#ifdef BL_LAZY
        });
#endif
#endif
    }
}
