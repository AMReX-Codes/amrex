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
                int max_iter, const BndryData& bndry, int max_bndry_order);
    
protected:

    static constexpr int stencil_size = 2*AMREX_SPACEDIM+1;
    
    MPI_Comm comm = MPI_COMM_NULL;
    
    int verbose = 0;
    
    int is_periodic[AMREX_SPACEDIM];
    
    Geometry geom;
    
    std::unique_ptr<MultiFab> acoefs;
    Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> bcoefs;

    Real scalar_a, scalar_b;
    
    HYPRE_StructGrid grid = NULL;
    HYPRE_StructMatrix A = NULL;
    HYPRE_StructVector b = NULL;
    HYPRE_StructVector x = NULL;
    HYPRE_StructSolver solver = NULL;
};

}

#endif

