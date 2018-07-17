#ifndef AMREX_HYPREABECLAP3_H_
#define AMREX_HYPREABECLAP3_H_


#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MacBndry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_LO_BCTYPES.H>

#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"
#include <algorithm>

namespace amrex
{

class HypreABecLap3
{
public:

    HypreABecLap3 (const BoxArray& grids,
                   const DistributionMapping& dmap,
                   const Geometry& geom_,
                   MPI_Comm comm_);
    ~HypreABecLap3 ();

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

    iMultiFab GbInd;

    std::vector<int> CellsGIndex;
    std::vector<int> numCellsProc;
    std::vector<int> StartIndex;
    
    HYPRE_IJMatrix A;
    HYPRE_ParCSRMatrix par_A;
    HYPRE_IJVector b;
    HYPRE_ParVector par_b;
    HYPRE_IJVector x;
    HYPRE_ParVector par_x;
    HYPRE_Solver solver;

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
