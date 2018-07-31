#ifndef AMREX_HypreABecLap_H_
#define AMREX_HypreABecLap_H_

#include <iomanip>
#include <iostream>
#include <memory>

#include <AMReX_Hypre.H>
#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MacBndry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_LO_BCTYPES.H>

#include <HYPRE_struct_ls.h>

namespace amrex {

class HypreABecLap
    : public Hypre
{
public:

    HypreABecLap (const BoxArray& grids, const DistributionMapping& dmap,
                  const Geometry& geom_, MPI_Comm comm_);
    virtual ~HypreABecLap ();
    
    virtual void solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
                        int max_iter, const BndryData& bndry, int max_bndry_order) final;
    
private:

    HYPRE_StructGrid grid = NULL;
    HYPRE_StructMatrix A = NULL;
    HYPRE_StructVector b = NULL;
    HYPRE_StructVector x = NULL;
    HYPRE_StructSolver solver = NULL;

    void prepareSolver ();
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
    void getSolution (MultiFab& soln);
};

}

#endif

