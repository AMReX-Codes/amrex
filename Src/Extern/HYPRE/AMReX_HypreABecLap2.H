#ifndef AMREX_HYPREABECLAP2_H_
#define AMREX_HYPREABECLAP2_H_
#include <AMReX_Config.H>

#include <AMReX_Hypre.H>

#include "_hypre_sstruct_mv.h"
#include "_hypre_sstruct_ls.h"

namespace amrex
{

class HypreABecLap2 final
    : public Hypre
{
public:

    HypreABecLap2 (const BoxArray& grids, const DistributionMapping& dmap,
                   const Geometry& geom_, MPI_Comm comm_);
    ~HypreABecLap2 () override;

    HypreABecLap2 (HypreABecLap2 const&) = delete;
    HypreABecLap2 (HypreABecLap2 &&) = delete;
    HypreABecLap2& operator= (HypreABecLap2 const&) = delete;
    HypreABecLap2& operator= (HypreABecLap2 &&) = delete;

    void solve (MultiFab& soln, const MultiFab& rhs, Real reltol, Real abstol,
                int maxiter, const BndryData& bndry, int max_bndry_order) override;

private :

    HYPRE_SStructGrid     hgrid = nullptr;
    HYPRE_SStructStencil  stencil = nullptr;
    HYPRE_SStructGraph    graph = nullptr;
    HYPRE_SStructMatrix   A = nullptr;
    HYPRE_SStructVector   b = nullptr;
    HYPRE_SStructVector   x = nullptr;
    HYPRE_Solver          solver = nullptr;

    void getSolution (MultiFab& soln);

public: // for cuda
    void prepareSolver ();
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
};

}

#endif
