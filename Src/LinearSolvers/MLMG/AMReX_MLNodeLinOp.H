#ifndef AMREX_ML_NODE_LINOP_H_H
#define AMREX_ML_NODE_LINOP_H_H
#include <AMReX_Config.H>

#include <AMReX_MLLinOp.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
#include <AMReX_HypreNodeLap.H>
#endif

namespace amrex {

/**
 * \file AMReX_MLNodeLinOp.H
 *
 * Base class for nodal linear operators.
 */

//! Coarsening strategy of nodal operators.
enum struct MLNodeLinOpCoarseningStrategy : int { Sigma, RAP };

template <typename MF>
class MLNodeLinOpT
    : public MLLinOpT<MF>
{
public:

    friend class HypreMLNodeLap;

    using FAB = typename MLLinOpT<MF>::FAB;
    using RT  = typename MLLinOpT<MF>::RT;
    using BCType    = typename MLLinOpT<MF>::BCType;
    using BCMode    = typename MLLinOpT<MF>::BCMode;
    using StateMode = typename MLLinOpT<MF>::StateMode;
    using Location  = typename MLLinOpT<MF>::Location;

    using CoarseningStrategy = MLNodeLinOpCoarseningStrategy;

    MLNodeLinOpT ();
    ~MLNodeLinOpT () override = default;

    MLNodeLinOpT (const MLNodeLinOpT&) = delete;
    MLNodeLinOpT (MLNodeLinOpT&&) = delete;
    MLNodeLinOpT& operator= (const MLNodeLinOpT&) = delete;
    MLNodeLinOpT& operator= (MLNodeLinOpT&&) = delete;

    /**
     * \brief Bind the nodal operator to an AMR hierarchy.
     *
     * \param a_geom    Per-level geometries.
     * \param a_grids   Per-level grids.
     * \param a_dmap    Distribution mappings.
     * \param a_info    LPInfo overrides.
     * \param a_factory Optional nodal factories per level.
     * \param a_eb_limit_coarsening True to forbid coarsening beyond where EB
     *                              geometry is available (-1 uses default).
     */
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FAB> const*>& a_factory = {},
                 int a_eb_limit_coarsening = -1);

    /**
     * \brief Control how many smoother passes are executed inside
     * each Fsmooth call.
     *
     * Used only by the host Gauss-Seidel sweeps of MLNodeLaplacian and
     * MLNodeABecLaplacian; ignored elsewhere.
     *
     * \param nsweeps Number of smoothing sweeps per Fsmooth invocation.
     */
    void setSmoothNumSweeps (int nsweeps) noexcept {
        m_smooth_num_sweeps = nsweeps;
        this->template spForward<MLNodeLinOpT>([&] (auto* sp) { sp->setSmoothNumSweeps(nsweeps); });
    }

    void setLevelBC (int /*amrlev*/, const MF* /*levelbcdata*/,
                     const MF* = nullptr, const MF* = nullptr,
                     const MF* = nullptr) final {}

    /**
     * \brief Apply the nodal operator with the requested BC/state modes.
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param out    Destination MultiFab updated with L(\p in).
     * \param in     Input MultiFab.
     * \param bc_mode Boundary-condition interpretation (solution vs correction).
     * \param s_mode  State interpretation (solution or correction).
     * \param bndry   Optional cached BC data.
     */
    void apply (int amrlev, int mglev, MF& out, MF& in, BCMode bc_mode,
                        StateMode s_mode, const MLMGBndryT<MF>* bndry=nullptr) const final;

    /**
     * \brief Invoke the nodal smoother on (\p amrlev,\p mglev).
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param sol    Correction updated in place.
     * \param rhs    Right-hand side for the smoother.
     * \param skip_fillboundary True to assume ghost nodes are already filled.
     * \param niter  Number of smoothing sweeps to perform.
     */
    void smooth (int amrlev, int mglev, MF& sol, const MF& rhs,
                 bool skip_fillboundary, int niter) const override;

    //! Residual helper that uses solution BCs and optional coarse data (operates on AMR level \p amrlev).
    void solutionResidual (int amrlev, MF& resid, MF& x, const MF& b,
                                   const MF* crse_bcdata=nullptr) override;
    //! Residual helper for correction solves (optionally homogeneous BCs) on (\p amrlev,\p mglev).
    void correctionResidual (int amrlev, int mglev, MF& resid, MF& x, const MF& b,
                                     BCMode bc_mode, const MF* crse_bcdata=nullptr) override;

    //! Return component-wise solvability offsets for (\p amrlev,\p mglev, \p rhs).
    Vector<RT> getSolvabilityOffset (int amrlev, int mglev,
                                               MF const& rhs) const override;
    //! Subtract offsets produced by getSolvabilityOffset() from \p rhs.
    void fixSolvabilityByOffset (int amrlev, int mglev, MF& rhs,
                                         Vector<RT> const& offset) const override;

    //! Prepare BC caches, masks, and metadata prior to solving.
    void prepareForSolve () override;

    //! Prepare auxiliary data used when the operator acts as a preconditioner.
    void preparePrecond () override;

    //! Zero Dirichlet nodes in \p mf before a correction or residual pass.
    void setDirichletNodesToZero (int amrlev, int mglev, MF& mf) const override;

    bool isSingular (int amrlev) const override
        { return (amrlev == 0) ? m_is_singular : false; }
    bool isBottomSingular () const override { return m_is_bottom_singular; }

    //! Dot-product helper.
    RT xdoty (int amrlev, int mglev, const MF& x, const MF& y, bool local) const final;

    //! Dot product used when the operator acts as a preconditioner.
    RT dotProductPrecond (Vector<MF const*> const& x,
                          Vector<MF const*> const& y) const final;

    //! L2 norm used when the operator acts as a preconditioner.
    RT norm2Precond (Vector<MF const*> const& x) const final;

    /**
     * \brief Apply physical BCs prior to smoothing/apply.
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param phi    MultiFab whose ghost values are set.
     * \param bc_mode Boundary-condition interpretation (solution vs correction).
     * \param state_mode State interpretation (solution vs correction).
     * \param skip_fillboundary True to skip FillBoundary calls.
     */
    virtual void applyBC (int amrlev, int mglev, MF& phi, BCMode bc_mode, StateMode state_mode,
                          bool skip_fillboundary=false) const;

    /**
     * \brief Low-level apply invoked by the base class.
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param out    Destination MultiFab receiving ``L(in)``.
     * \param in     Source MultiFab.
     */
    virtual void Fapply (int amrlev, int mglev, MF& out, const MF& in) const = 0;
    /**
     * \brief Low-level smoother invoked by the base class.
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param sol    Correction updated in place.
     * \param rhs    Right-hand side for the smoother.
     */
    virtual void Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs) const = 0;

    //! Synchronize nodal values across box/partition boundaries on (\p amrlev,\p mglev).
    void nodalSync (int amrlev, int mglev, MF& mf) const;

    /**
     * \brief Build a nodal mask that marks points owned by this rank (used for sync).
     *
     * \param ba   Grid layout on the current level.
     * \param dm   Distribution mapping describing MPI ownership.
     * \param geom Geometry describing periodicity/domain info.
     */
    static std::unique_ptr<iMultiFab> makeOwnerMask (const BoxArray& ba,
                                                     const DistributionMapping& dm,
                                                     const Geometry& geom);

    //! Build Dirichlet/dot-product masks and determine bottom singularity.
    void buildMasks ();

    // omask is either 0 or 1. 1 means the node is an unknown. 0 means it's known.
    /**
     * \brief Mark nodes that belong to overset (known) regions.
     *
     * \param amrlev AMR level index.
     * \param a_dmask Mask with 1 = unknown node, 0 = known/covered node.
     */
    void setOversetMask (int amrlev, const iMultiFab& a_dmask);

    //! Zero \p mf at the known nodes of the overset mask.
    void applyOverset (int amrlev, MF& mf) const override;

    virtual void fixUpResidualMask (int /*amrlev*/, iMultiFab& /*resmsk*/) { }

    //! Infinity norm helper operating on level \p amrlev (set \p local to true to skip MPI reduction).
    RT normInf (int amrlev, MF const& mf, bool local) const override;

    //! Nodal average-down for residuals between AMR levels (no-op here).
    void avgDownResAmr (int, MF&, MF const&) const final { }

    /**
     * \brief Prolong AMR-level data during FMG initialization.
     *
     * \param famrlev Fine AMR level index.
     * \param fine    Destination MultiFab.
     * \param crse    Source coarse MultiFab.
     * \param nghost  Grow cells to fill.
     */
    void interpolationAmr (int famrlev, MF& fine, const MF& crse,
                                   IntVect const& nghost) const override;

    /**
     * \brief Average nodal solutions down the AMR hierarchy and synchronize overlaps.
     *
     * \param sol Vector of nodal MultiFabs spanning all AMR levels.
     */
    void averageDownAndSync (Vector<MF>& sol) const override;

    /**
     * \brief Overwrite fine data with the prolongation of coarse data
     * (\p fine = prolong(\p crse)).
     *
     * \param amrlev AMR level index.
     * \param fmglev Fine MG level index.
     * \param fine   Destination fine MultiFab (overwritten in place).
     * \param crse   Source coarse MultiFab (its ghost cells may be modified).
     */
    void interpAssign (int amrlev, int fmglev, MF& fine, MF& crse) const override;

#if (AMREX_SPACEDIM > 1)
    [[nodiscard]] std::unique_ptr<MLAlgMG> makeAlgMG (int mglev) const override;
#endif

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    /**
     * \brief Build a Hypre nodal Laplacian using the current metadata.
     *
     * \param bottom_verbose Verbosity forwarded to Hypre.
     * \param options_namespace ParmParse namespace for Hypre options.
     */
    [[nodiscard]] std::unique_ptr<HypreNodeLap> makeHypreNodeLap(
        int bottom_verbose,
        const std::string& options_namespace) const override;

    virtual void fillIJMatrix (MFIter const& /*mfi*/,
                               Array4<HypreNodeLap::AtomicInt const> const& /*gid*/,
                               Array4<int const> const& /*lid*/,
                               HypreNodeLap::Int* /*ncols*/,
                               HypreNodeLap::Int* /*cols*/,
                               Real* /*mat*/) const
    {
        amrex::Abort("MLNodeLinOp::fillIJMatrix: how did we get here?");
    }
#endif

    /**
     * \brief Fill the matrix rows of one box for the algebraic solver.
     *
     * Same packed layout as the Hypre version: ncols[lid] entries per row,
     * diagonal first, neighbors without a row dropped. The MG level is on
     * AMR level 0. Global row ids are numeric_limits<Long>::max() and local
     * ids negative where there is no row.
     */
    virtual void fillAlgMatrix (int /*mglev*/, MFIter const& /*mfi*/,
                                Array4<Long const> const& /*gid*/,
                                Array4<int const> const& /*lid*/,
                                Long* /*ncols*/, Long* /*cols*/, Real* /*mat*/) const
    {
        amrex::Abort("MLNodeLinOp::fillAlgMatrix: not supported by this operator");
    }

    //! Fill the right-hand side of the rows of one box on MG level \p mglev.
    virtual void fillRHS (int /*mglev*/, MFIter const& /*mfi*/,
                          Array4<int const> const& /*lid*/,
                          Real* /*rhs*/,
                          Array4<Real const> const& /*bfab*/) const
    {
        amrex::Abort("MLNodeLinOp:fillRHS: how did we get here?");
    }

protected:

    void spCopySettings () override
    {
        MLLinOpT<MF>::spCopySettings();
        setSmoothNumSweeps(m_smooth_num_sweeps);
    }

    //! Convert a geometry array to RT.
    [[nodiscard]] static GpuArray<RT,AMREX_SPACEDIM>
    toRT (GpuArray<Real,AMREX_SPACEDIM> const& a) noexcept
    {
        return {AMREX_D_DECL(static_cast<RT>(a[0]), static_cast<RT>(a[1]), static_cast<RT>(a[2]))};
    }

    void resizeMultiGrid (int new_size) override;

    void buildMGHierarchy () override;

    std::unique_ptr<iMultiFab> m_owner_mask_top;  // ownership of nodes
    std::unique_ptr<iMultiFab> m_owner_mask_bottom;
    Vector<Vector<std::unique_ptr<iMultiFab> > > m_dirichlet_mask;  // dirichlet?
    Vector<std::unique_ptr<iMultiFab> > m_cc_fine_mask;          // cell-centered mask for cells covered by fine
    Vector<std::unique_ptr<iMultiFab> > m_nd_fine_mask;          // nodal mask: 0: this level node, 1: c/f boundary, 2: fine node
    Vector<std::unique_ptr<LayoutData<int> > > m_has_fine_bndry; // does this fab contain c/f boundary?
    MF m_bottom_dot_mask;
    MF m_coarse_dot_mask;
    Vector<MF> m_precond_weight_mask;

    Vector<std::unique_ptr<iMultiFab> > m_norm_fine_mask;

#ifdef AMREX_USE_EB
    CoarseningStrategy m_coarsening_strategy = CoarseningStrategy::RAP;
#else
    CoarseningStrategy m_coarsening_strategy = CoarseningStrategy::Sigma;
#endif

    bool m_masks_built = false;
    bool m_overset_dirichlet_mask = false;
    int m_smooth_num_sweeps = 2;
    mutable bool m_in_solution_mode = true;

private:
    bool m_is_singular = false;
    bool m_is_bottom_singular = false; // also when coarsening dropped all overset nodes
};

extern template class MLNodeLinOpT<MultiFab>;
#ifdef AMREX_MLMG_HAS_FLOAT
extern template class MLNodeLinOpT<fMultiFab>;
#endif

//! Double precision base class for nodal linear operators.
using MLNodeLinOp = MLNodeLinOpT<MultiFab>;

}

#endif
