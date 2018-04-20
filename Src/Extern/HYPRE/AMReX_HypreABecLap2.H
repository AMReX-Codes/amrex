#ifndef AMREX_HYPREABECLAP2_H_
#define AMREX_HYPREABECLAP2_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MacBndry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_LO_BCTYPES.H>

#include "_hypre_sstruct_mv.h"
#include "_hypre_sstruct_ls.h"

namespace amrex
{

class HypreABecLap2
{
public:

    HypreABecLap2 (const BoxArray& grids,
                   const DistributionMapping& dmap,
                   const Geometry& geom_,
                   MPI_Comm comm_);
    ~HypreABecLap2 ();

    void setScalars(Real sa, Real sb);
    void setACoeffs(const MultiFab& alpha);
    void setBCoeffs(const std::array<const MultiFab*,BL_SPACEDIM>& beta);
    void setVerbose(int _verbose);
    void solve(MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
               int max_iter, LinOpBCType bc_type, Real bc_value);
    void solve(MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
               int max_iter, const BndryData& _bndry);

private :

    MPI_Comm comm;
    Geometry geom;

    int verbose;

    MultiFab acoefs;
    std::array<MultiFab,BL_SPACEDIM> bcoefs;
    Real scalar_a, scalar_b;
    BndryData bd;

    HYPRE_SStructGrid     hgrid;
    HYPRE_SStructStencil  stencil;
    HYPRE_SStructGraph    graph;
    HYPRE_SStructMatrix   A, A0;
    HYPRE_SStructVector   b;
    HYPRE_SStructVector   x;
    HYPRE_Solver          solver;

    Real rel_tol, abs_tol;
    int max_iter;

    void loadBndryData (LinOpBCType bc_type, Real bc_value);
    void loadMatrix ();
    void finalizeMatrix ();
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
    void finalizeVectors ();

    void setupSolver (Real rel_tol_, Real abs_tol_, int max_iter_);
    void solveDoIt ();
    void clearSolver ();

    void getSolution (MultiFab& soln);
};

}

#endif
