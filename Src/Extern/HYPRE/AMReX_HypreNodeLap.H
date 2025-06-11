#ifndef AMREX_HYPRE_NODE_LAP_H_
#define AMREX_HYPRE_NODE_LAP_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_BndryData.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_HypreIJIface.H>

#include <algorithm>
#include <memory>
#include <type_traits>

namespace amrex {

class MLNodeLinOp;

class HypreNodeLap
{
public:
    HypreNodeLap (const BoxArray& grids_, const DistributionMapping& dmap_,
                  const Geometry& geom_, const FabFactory<FArrayBox>& factory_,
                  const iMultiFab& owner_mask_, const iMultiFab& dirichlet_mask_,
                  MPI_Comm comm_, MLNodeLinOp const* linop_, int verbose_,
                  std::string options_namespace_);
    ~HypreNodeLap ();

    HypreNodeLap (HypreNodeLap const&) = delete;
    HypreNodeLap (HypreNodeLap &&) = delete;
    HypreNodeLap& operator= (HypreNodeLap const&) = delete;
    HypreNodeLap& operator= (HypreNodeLap &&) = delete;

    using Int = HYPRE_Int;
    // To generate global ids for Hypre, we need to remove duplicates on
    // nodes shared by multiple Boxes.  So we need to use a type that
    // supports atomicAdd.  HYPRE_Int is either int or long long.  The
    // latter does not have native atomicAdd support in CUDA/HIP.
    using AtomicInt = std::conditional_t<sizeof(Int) == 4, unsigned int, unsigned long long>;

    void solve (MultiFab& soln, const MultiFab& rhs,
                Real rel_tol, Real abs_tol, int max_iter);

    void setHypreOptionsNamespace(const std::string& ns)
    { options_namespace = ns; }

private:

    BoxArray grids;
    DistributionMapping dmap;
    Geometry geom;
    FabFactory<FArrayBox> const* factory = nullptr;
    iMultiFab const* owner_mask = nullptr;
    iMultiFab const* dirichlet_mask = nullptr;
    MPI_Comm comm = MPI_COMM_NULL;
    MLNodeLinOp const* linop = nullptr;
    int verbose = 0;

    std::unique_ptr<HypreIJIface> hypre_ij;

    // Non-owning references to Hypre matrix, rhs, and solution data
    HYPRE_IJMatrix A = nullptr;
    HYPRE_IJVector b = nullptr;
    HYPRE_IJVector x = nullptr;

    LayoutData<Int> nnodes_grid;
    LayoutData<Gpu::DeviceVector<Int> > node_id_vec;
    FabArray<BaseFab<AtomicInt> > node_id;
    iMultiFab local_node_id;
    LayoutData<Int> id_offset;
    MultiFab tmpsoln;

    std::string options_namespace{"hypre"};

public: // for CUDA
    Int fill_local_node_id ();
#ifdef AMREX_USE_GPU
    Int fill_local_node_id_gpu ();
#endif
    Int fill_local_node_id_cpu ();
    void fill_global_node_id ();
    static void adjust_singular_matrix (Int const* ncols, Int const* cols, Int const* rows, Real* mat);
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
    void getSolution (MultiFab& soln);
};

}

#endif
