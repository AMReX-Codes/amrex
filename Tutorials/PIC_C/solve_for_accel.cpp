#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MacBndry.H>
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#include <stencil_types.H>

#include "Particles.H"

void solve_with_f90  (MultiFab& rhs, PArray<MultiFab>& grad_phi_edge, const Geometry& geom, Real tol, Real abs_tol);
void solve_with_hpgmg(MultiFab& rhs, PArray<MultiFab>& grad_phi_edge, const Geometry& geom, Real tol, Real abs_tol);

void 
solve_for_accel(MultiFab& rhs, MultiFab& grad_phi, const Geometry& geom)
{
 
    Real tol     = 1.e-10;
    Real abs_tol = 1.e-12;

    int MAX_LEV = 10;
    Array< PArray<MultiFab> > grad_phi_edge;
    grad_phi_edge.resize(MAX_LEV);
    grad_phi_edge[0].resize(BL_SPACEDIM, PArrayManage);
    for (int n = 0; n < BL_SPACEDIM; ++n)
        grad_phi_edge[0].set(n, new MultiFab(BoxArray(rhs.boxArray()).surroundingNodes(n), 1, 1));

    Real     strt    = ParallelDescriptor::second();

#ifdef USEHPGMG
   solve_with_hpgmg(rhs,grad_phi_edge[0],geom,tol,abs_tol);
#else
   solve_with_f90  (rhs,grad_phi_edge[0],geom,tol,abs_tol);
#endif

    // Average edge-centered gradients to cell centers.
    BoxLib::average_face_to_cellcenter(grad_phi, grad_phi_edge[0], geom);
    geom.FillPeriodicBoundary(grad_phi,true);

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
