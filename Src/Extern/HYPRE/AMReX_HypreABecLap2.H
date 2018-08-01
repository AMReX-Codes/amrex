#ifndef AMREX_HYPREABECLAP2_H_
#define AMREX_HYPREABECLAP2_H_

#include <AMReX_Hypre.H>

#include "_hypre_sstruct_mv.h"
#include "_hypre_sstruct_ls.h"

namespace amrex
{

class HypreABecLap2
    : public Hypre
{
public:

    HypreABecLap2 (const BoxArray& grids, const DistributionMapping& dmap,
                   const Geometry& geom_, MPI_Comm comm_);
    virtual ~HypreABecLap2 ();

    virtual void solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
                        int max_iter, const BndryData& bndry, int max_bndry_order) final;

private :

    HYPRE_SStructGrid     hgrid = NULL;
    HYPRE_SStructStencil  stencil = NULL;
    HYPRE_SStructGraph    graph = NULL;
    HYPRE_SStructMatrix   A = NULL;
    HYPRE_SStructVector   b = NULL;
    HYPRE_SStructVector   x = NULL;
    HYPRE_Solver          solver = NULL;

    void prepareSolver ();
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
    void getSolution (MultiFab& soln);
};

}

#endif
