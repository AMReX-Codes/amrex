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

/**
 * \file AMReX_HypreNodeLap.H
 *
 * HYPRE IJ wrapper for nodal Laplacians used by MLNodeLinOp.
 */

namespace amrex {

class MLNodeLinOp;

/**
 * \brief IJ-based nodal Laplacian that mirrors MLNodeLinOp's discretization.
 */
class HypreNodeLap
{
public:
    /**
     * \brief Construct an IJ-based nodal Laplacian mirroring MLNodeLinOp.
     *
     * \param grids_           Cell-centered grids that define the nodal mesh.
     * \param dmap_            Distribution mapping for \p grids_.
     * \param geom_            Geometry describing the level.
     * \param factory_         Fab factory referenced when creating workspace.
     * \param owner_mask_      Ownership mask used to deduplicate shared nodes.
     * \param dirichlet_mask_  Mask describing Dirichlet nodes.
     * \param comm_            MPI communicator supplied to HYPRE.
     * \param linop_           MLNodeLinOp used as the reference discretization.
     * \param verbose_         Verbosity forwarded to HYPRE.
     * \param options_namespace_ ParmParse namespace for runtime options.
     */
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
    // To generate global ids for HYPRE, we need to remove duplicates on
    // nodes shared by multiple Boxes.  So we need to use a type that
    // supports atomicAdd.  HYPRE_Int is either int or long long.  The
    // latter does not have native atomicAdd support in CUDA/HIP.
    using AtomicInt = std::conditional_t<sizeof(Int) == 4, unsigned int, unsigned long long>;

    /**
     * \brief Solve the nodal Laplacian using HYPRE's IJ interface.
     *
     * \param soln    Solution MultiFab (nodal) updated on return.
     * \param rhs     RHS MultiFab (nodal).
     * \param rel_tol Relative tolerance.
     * \param abs_tol Absolute tolerance.
     * \param max_iter Maximum allowed iterations.
     */
    void solve (MultiFab& soln, const MultiFab& rhs,
                Real rel_tol, Real abs_tol, int max_iter);

    //! Override the ParmParse namespace used when configuring HYPRE (\p ns).
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

    // Non-owning references to HYPRE matrix, rhs, and solution data
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
    /**
     * \brief Assign local ids to each owned node.
     *
     * \return Total owned-node count on this rank.
     */
    Int fill_local_node_id ();
#ifdef AMREX_USE_GPU
    //! GPU-specific local-id fill routine; returns the owned-node count.
    Int fill_local_node_id_gpu ();
#endif
    //! CPU fallback local-id fill routine; returns the owned-node count.
    Int fill_local_node_id_cpu ();
    //! Convert local ids to globally unique ids recognized by HYPRE.
    void fill_global_node_id ();
    /**
     * \brief Apply the fix-up HYPRE expects for singular matrices (Dirichlet pinning).
     *
     * \param ncols Number of columns in each CSR row.
     * \param cols  Column indices.
     * \param rows  Row indices.
     * \param mat   Matrix entries updated in place.
     */
    static void adjust_singular_matrix (Int const* ncols, Int const* cols, Int const* rows, Real* mat);
    //! Copy AMReX solution guess \p soln and RHS \p rhs into HYPRE's IJ vectors prior to solving.
    void loadVectors (MultiFab& soln, const MultiFab& rhs);
    //! Copy the solved data \p soln out of HYPRE and back into AMReX storage.
    void getSolution (MultiFab& soln);
};

}

#endif
