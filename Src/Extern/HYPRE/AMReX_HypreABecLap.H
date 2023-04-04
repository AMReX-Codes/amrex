#ifndef AMREX_HypreABecLap_H_
#define AMREX_HypreABecLap_H_
#include <AMReX_Config.H>

#include <AMReX_Hypre.H>
#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BndryData.H>
#include <AMReX_ParmParse.H>
#include <AMReX_LO_BCTYPES.H>

#include <HYPRE_struct_ls.h>

#include <memory>

namespace amrex {

class HypreABecLap final
    : public Hypre
{
public:

    HypreABecLap (const BoxArray& grids, const DistributionMapping& dmap,
                  const Geometry& geom_, MPI_Comm comm_);
    ~HypreABecLap () override;

    HypreABecLap (HypreABecLap const&) = delete;
    HypreABecLap (HypreABecLap &&) = delete;
    HypreABecLap& operator= (HypreABecLap const&) = delete;
    HypreABecLap& operator= (HypreABecLap &&) = delete;

    void solve (MultiFab& soln, const MultiFab& rhs, Real reltol, Real abstol,
                int maxiter, const BndryData& bndry, int max_bndry_order) override;

private:

    HYPRE_StructGrid grid = nullptr;
    HYPRE_StructMatrix A = nullptr;
    HYPRE_StructVector b = nullptr;
    HYPRE_StructVector x = nullptr;
    HYPRE_StructSolver solver = nullptr;

    void getSolution (MultiFab& soln);

public: // for cuda
    void prepareSolver ();
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
};

}

#endif
