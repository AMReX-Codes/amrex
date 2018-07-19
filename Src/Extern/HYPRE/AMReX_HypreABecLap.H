#ifndef AMREX_HypreABecLap_H_
#define AMREX_HypreABecLap_H_

#include <iomanip>
#include <iostream>
#include <memory>

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MacBndry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_LO_BCTYPES.H>

#include <HYPRE_struct_ls.h>

namespace amrex {

class HypreABecLap {

public:

    HypreABecLap (const BoxArray& grids,
                  const DistributionMapping& dmap,
                  const Geometry& geom_,
                  MPI_Comm comm_);
    ~HypreABecLap ();
    
    void setScalars (Real sa, Real sb);
    void setACoeffs (const MultiFab& alpha);
    void setBCoeffs (const Array<const MultiFab*,AMREX_SPACEDIM>& beta);
    void setVerbose (int _verbose);
    void solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
                int max_iter, const BndryData& bndry);
    
protected:

    static constexpr int stencil_size = 2*AMREX_SPACEDIM+1;
    
    MPI_Comm comm;
    
    int verbose;
    
    int is_periodic[AMREX_SPACEDIM];
    Real dx[AMREX_SPACEDIM];
    
    Geometry geom;
    
    std::unique_ptr<MultiFab> acoefs;
    Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> bcoefs;

    Real scalar_a, scalar_b;
    
    int maxiter;
    Real reltol, abstol;
    
    HYPRE_StructGrid grid;
    HYPRE_StructMatrix A;
    HYPRE_StructVector b, x;
    HYPRE_StructSolver solver, precond;
};

}

#endif

