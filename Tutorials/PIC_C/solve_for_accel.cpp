#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MacBndry.H>
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#include <stencil_types.H>

void solve_with_f90  (PArray<MultiFab>& rhs, PArray<MultiFab>& phi, Array< PArray<MultiFab> >& grad_phi_edge, 
                      const PArray<Geometry>& geom, int base_level, Real tol, Real abs_tol);
#ifdef USEHPGMG
void solve_with_hpgmg(PArray<MultiFab>& rhs, Array< PArray<MultiFab> >& grad_phi_edge, const PArray<Geometry>& geom, Real tol, Real abs_tol);
#endif

void 
solve_for_accel(PArray<MultiFab>& rhs, PArray<MultiFab>& phi, PArray<MultiFab>& grad_phi, 
		const PArray<Geometry>& geom, int base_level)
{
 
    Real tol     = 1.e-10;
    Real abs_tol = 1.e-14;

    Array< PArray<MultiFab> > grad_phi_edge;
    grad_phi_edge.resize(rhs.size());

    for (int lev = 0; lev < rhs.size(); lev++)
    {
        grad_phi_edge[lev].resize(BL_SPACEDIM, PArrayManage);
        for (int n = 0; n < BL_SPACEDIM; ++n)
            grad_phi_edge[lev].set(n, new MultiFab(BoxArray(rhs[lev].boxArray()).surroundingNodes(n), 1, 1));
    }

    Real     strt    = ParallelDescriptor::second();

    // ***************************************************
    // Make sure the RHS sums to 0 if fully periodic
    // ***************************************************
    for (int lev = base_level; lev < rhs.size(); lev++) {
	Real n0 = rhs[lev].norm0();
	if (ParallelDescriptor::IOProcessor())
	    std::cout << "Max of rhs in solve_for_phi before correction at level  " << lev << " " << n0 << std::endl;
    }

    // This is a correction for fully periodic domains only
    if ( Geometry::isAllPeriodic() && (rhs[base_level].boxArray().numPts() == geom[base_level].Domain().numPts()) )
    {
	Real sum0 = rhs[base_level].sum(0);
	Real npts = rhs[base_level].boxArray().d_numPts();
        Real sum = sum0/npts;

	if (ParallelDescriptor::IOProcessor()) {
	    std::cout << "Sum                                  is " << sum0 << std::endl;
	    std::cout << "Npts                                 is " << npts << std::endl;
	    std::cout << "Sum of particle weights over level 0 is " << sum  << std::endl;
	}

        for (int lev = base_level; lev < rhs.size(); lev++)
            rhs[lev].plus(-sum, 0, 1, 0);

	for (int lev = base_level; lev < rhs.size(); lev++) {
	    Real n0 = rhs[lev].norm0();
	    if (ParallelDescriptor::IOProcessor())
		std::cout << "Max of rhs in solve_for_phi  after correction at level  " << lev << " " << n0 << std::endl;
	}
    }

    // ***************************************************
    // Solve for phi and return both phi and grad_phi_edge
    // ***************************************************

#ifdef USEHPGMG
   solve_with_hpgmg(rhs,phi,grad_phi_edge,geom,base_level,tol,abs_tol);
#else
   solve_with_f90  (rhs,phi,grad_phi_edge,geom,base_level,tol,abs_tol);
#endif

    // Average edge-centered gradients to cell centers.
    for (int lev = 0; lev < rhs.size(); lev++)
    {
        BoxLib::average_face_to_cellcenter(grad_phi[lev], grad_phi_edge[lev], geom[lev]);
        geom[lev].FillPeriodicBoundary(grad_phi[lev],true);  // wz: why only fill periodic boundary?
    }

    // VisMF::Write(grad_phi,"GradPhi");

    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "solve_for_phi() time = " << end << std::endl;
#ifdef BL_LAZY
        });
#endif
    }
}
