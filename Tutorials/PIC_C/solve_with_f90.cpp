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

void 
solve_with_f90(MultiFab& rhs, MultiFab& grad_phi, const Geometry& geom)
{
    Real     strt    = ParallelDescriptor::second();

    BCRec* phys_bc;
    int mg_bc[2*BL_SPACEDIM];

    int MAX_LEV = 10;
    Array< PArray<MultiFab> > grad_phi_edge;
    grad_phi_edge.resize(MAX_LEV);
    grad_phi_edge[0].resize(BL_SPACEDIM, PArrayManage);
    for (int n = 0; n < BL_SPACEDIM; ++n)
        grad_phi_edge[0].set(n, new MultiFab(BoxArray(rhs.boxArray()).surroundingNodes(n), 1, 1));

    // build (and initialize) a multifab for the solution on the box array with 1 component, 1 ghost cell
    MultiFab phi(BoxArray(rhs.boxArray()), 1, 1);
    phi.setVal(0.);

    std::cout << "Max of rhs in solve_for_phi before correction " << rhs.norm0() << std::endl;

    // This is a correction for fully periodic domains only
    if (Geometry::isAllPeriodic())
    {
        Real numpts = rhs.boxArray().d_numPts();
        Real sum = rhs.sum(0) / numpts;
        rhs.plus(-sum, 0, 1, 0);
    }

    std::cout << "Max of rhs in solve_for_phi  after correction " << rhs.norm0() << std::endl;

    // This tells the solver that we are using periodic bc's
    if (Geometry::isAllPeriodic())
    {
        for (int dir = 0; dir < BL_SPACEDIM; ++dir)
        {
            mg_bc[2*dir + 0] = 0;
            mg_bc[2*dir + 1] = 0;
        }
    }

    MacBndry bndry(rhs.boxArray(), 1, geom);
    //
    // Set Dirichlet boundary condition for phi in phi grow cells, use to
    // initialize bndry.
    //
    const int src_comp  = 0;
    const int dest_comp = 0;
    const int num_comp  = 1;

    // Need to set the boundary values here so they can get copied into "bndry"  
    {
        bndry.setBndryValues(phi, src_comp, dest_comp, num_comp, *phys_bc);
    }

    std::vector<BoxArray> bav(1);
    bav[0] = phi.boxArray();
    std::vector<DistributionMapping> dmv(1);
    dmv[0] = rhs.DistributionMap();
    std::vector<Geometry> fgeom(1);
    fgeom[0] = geom;

    Array< Array<Real> > xa(1);
    Array< Array<Real> > xb(1);

    xa[0].resize(BL_SPACEDIM);
    xb[0].resize(BL_SPACEDIM);

//  if (level == 0)
    {
        for (int i = 0; i < BL_SPACEDIM; ++i)
        {
            xa[0][i] = 0;
            xb[0][i] = 0;
        }
    }
#if 0
    else
    {
        const Real* dx_crse = parent->Geom(level-1).CellSize();
        for (int i = 0; i < BL_SPACEDIM; ++i)
        {
            xa[0][i] = 0.5 * dx_crse[i];
            xb[0][i] = 0.5 * dx_crse[i];
        }
    }
#endif

    MultiFab* phi_p[1] = {&phi};
    MultiFab* rhs_p[1] = {&rhs};

    const Real  tol     = 1.e-10;
    const Real  abs_tol = 0.;

    {
        int stencil_type = CC_CROSS_STENCIL;
        int always_use_bnorm = 1;
        int need_grad_phi = 1;
        Real final_resnorm;
        MGT_Solver mgt_solver(fgeom, mg_bc, bav, dmv, false, stencil_type, false, 0, 1, 2);
        mgt_solver.set_const_gravity_coeffs(xa, xb);
        mgt_solver.solve(phi_p, rhs_p, bndry, tol, abs_tol, always_use_bnorm, final_resnorm,need_grad_phi); 

        const Real* dx = geom.CellSize();
        mgt_solver.get_fluxes(0,grad_phi_edge[0],dx);
    }

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
