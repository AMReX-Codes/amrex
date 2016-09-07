#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MacBndry.H>
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#include <stencil_types.H>
#include <VisMF.H>

void solve_with_f90  (PArray<MultiFab>& rhs, PArray<MultiFab>& phi, Array< PArray<MultiFab> >& grad_phi_edge, 
                      const Array<Geometry>& geom, int base_level, int finest_level, Real tol, Real abs_tol);
void 
solve_for_accel(PArray<MultiFab>& rhs, PArray<MultiFab>& phi, PArray<MultiFab>& grad_phi, 
		const Array<Geometry>& geom, int base_level, int finest_level, Real offset)
{
 
    Real tol     = 1.e-10;
    Real abs_tol = 1.e-14;

    Array< PArray<MultiFab> > grad_phi_edge;
    grad_phi_edge.resize(rhs.size());

    for (int lev = base_level; lev <= finest_level ; lev++)
    {
        grad_phi_edge[lev].resize(BL_SPACEDIM, PArrayManage);
        for (int n = 0; n < BL_SPACEDIM; ++n)
            grad_phi_edge[lev].set(n, new MultiFab(BoxArray(rhs[lev].boxArray()).surroundingNodes(n), 1, 1));
    }

    // ***************************************************
    // Make sure the RHS sums to 0 if fully periodic
    // ***************************************************
    for (int lev = base_level; lev <= finest_level; lev++) {
	Real n0 = rhs[lev].norm0();
	if (ParallelDescriptor::IOProcessor())
	    std::cout << "Max of rhs in solve_for_phi before correction at level  " 
                      << lev << " " << n0 << std::endl;
    }

    for (int lev = base_level; lev <= finest_level; lev++)
        rhs[lev].plus(-offset, 0, 1, 0);

    for (int lev = base_level; lev <= finest_level; lev++) {
	Real n0 = rhs[lev].norm0();
	if (ParallelDescriptor::IOProcessor())
	    std::cout << "Max of rhs in solve_for_phi  after correction at level  " 
                      << lev << " " << n0 << std::endl;
    }

    // ***************************************************
    // Solve for phi and return both phi and grad_phi_edge
    // ***************************************************

   solve_with_f90  (rhs,phi,grad_phi_edge,geom,base_level,finest_level,tol,abs_tol);

    // Average edge-centered gradients to cell centers and fill the values in ghost cells.
    for (int lev = base_level; lev <= finest_level; lev++)
    {
        BoxLib::average_face_to_cellcenter(grad_phi[lev], grad_phi_edge[lev], geom[lev]);
	grad_phi[lev].FillBoundary(0,BL_SPACEDIM,geom[lev].periodicity());
    }

    // VisMF::Write(grad_phi[0],"GradPhi");
}
