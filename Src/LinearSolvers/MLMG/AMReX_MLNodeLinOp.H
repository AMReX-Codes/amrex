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

class MLNodeLinOp
    : public MLLinOp
{
public:

    enum struct CoarseningStrategy : int { Sigma, RAP };

    MLNodeLinOp ();
    ~MLNodeLinOp () override = default;

    MLNodeLinOp (const MLNodeLinOp&) = delete;
    MLNodeLinOp (MLNodeLinOp&&) = delete;
    MLNodeLinOp& operator= (const MLNodeLinOp&) = delete;
    MLNodeLinOp& operator= (MLNodeLinOp&&) = delete;

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
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {},
                 int a_eb_limit_coarsening = -1);

    /**
     * \brief Control how many smoother passes are executed inside
     * each Fsmooth call.
     *
     * \param nsweeps Number of smoothing sweeps per Fsmooth invocation.
     */
    void setSmoothNumSweeps (int nsweeps) noexcept {
        m_smooth_num_sweeps = nsweeps;
    }

    void setLevelBC (int /*amrlev*/, const MultiFab* /*levelbcdata*/,
                     const MultiFab* = nullptr, const MultiFab* = nullptr,
                     const MultiFab* = nullptr) final {}

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
    void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                        StateMode s_mode, const MLMGBndry* bndry=nullptr) const final;

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
    void smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                 bool skip_fillboundary, int niter) const override;

    //! Residual helper that uses solution BCs and optional coarse data (operates on AMR level \p amrlev).
    void solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                   const MultiFab* crse_bcdata=nullptr) override;
    //! Residual helper for correction solves (optionally homogeneous BCs) on (\p amrlev,\p mglev).
    void correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                     BCMode bc_mode, const MultiFab* crse_bcdata=nullptr) override;

    //! Return component-wise solvability offsets for (\p amrlev,\p mglev, \p rhs).
    Vector<Real> getSolvabilityOffset (int amrlev, int mglev,
                                               MultiFab const& rhs) const override;
    //! Subtract offsets produced by getSolvabilityOffset() from \p rhs.
    void fixSolvabilityByOffset (int amrlev, int mglev, MultiFab& rhs,
                                         Vector<Real> const& offset) const override;

    //! Prepare BC caches, masks, and metadata prior to solving.
    void prepareForSolve () override;

    //! Prepare auxiliary data used when the operator acts as a preconditioner.
    void preparePrecond () override;

    //! Zero Dirichlet nodes in \p mf before a correction or residual pass.
    void setDirichletNodesToZero (int amrlev, int mglev, MultiFab& mf) const override;

    bool isSingular (int amrlev) const override
        { return (amrlev == 0) ? m_is_bottom_singular : false; }
    bool isBottomSingular () const override { return m_is_bottom_singular; }

    //! Dot-product helper.
    Real xdoty (int amrlev, int mglev, const MultiFab& x, const MultiFab& y, bool local) const final;

    //! Dot product used when the operator acts as a preconditioner.
    Real dotProductPrecond (Vector<MultiFab const*> const& x,
                            Vector<MultiFab const*> const& y) const final;

    //! L2 norm used when the operator acts as a preconditioner.
    Real norm2Precond (Vector<MultiFab const*> const& x) const final;

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
    virtual void applyBC (int amrlev, int mglev, MultiFab& phi, BCMode bc_mode, StateMode state_mode,
                          bool skip_fillboundary=false) const;

    /**
     * \brief Low-level apply invoked by the base class.
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param out    Destination MultiFab receiving ``L(in)``.
     * \param in     Source MultiFab.
     */
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const = 0;
    /**
     * \brief Low-level smoother invoked by the base class.
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param sol    Correction updated in place.
     * \param rhs    Right-hand side for the smoother.
     */
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const = 0;

    //! Synchronize nodal values across box/partition boundaries on (\p amrlev,\p mglev).
    void nodalSync (int amrlev, int mglev, MultiFab& mf) const;

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

    virtual void fixUpResidualMask (int /*amrlev*/, iMultiFab& /*resmsk*/) { }

    //! Infinity norm helper operating on level \p amrlev (set \p local to true to skip MPI reduction).
    Real normInf (int amrlev, MultiFab const& mf, bool local) const override;

    //! Nodal average-down for residuals between AMR levels (no-op here).
    void avgDownResAmr (int, MultiFab&, MultiFab const&) const final { }

    /**
     * \brief Prolong AMR-level data during FMG initialization.
     *
     * \param famrlev Fine AMR level index.
     * \param fine    Destination MultiFab.
     * \param crse    Source coarse MultiFab.
     * \param nghost  Grow cells to fill.
     */
    void interpolationAmr (int famrlev, MultiFab& fine, const MultiFab& crse,
                                   IntVect const& nghost) const override;

    /**
     * \brief Average nodal solutions down the AMR hierarchy and synchronize overlaps.
     *
     * \param sol Vector of nodal MultiFabs spanning all AMR levels.
     */
    void averageDownAndSync (Vector<MultiFab>& sol) const override;

    /**
     * \brief Overwrite fine data with the prolongation of coarse data
     * (\p fine = prolong(\p crse)).
     *
     * \param amrlev AMR level index.
     * \param fmglev Fine MG level index.
     * \param fine   Destination fine MultiFab (overwritten in place).
     * \param crse   Source coarse MultiFab (its ghost cells may be modified).
     */
    void interpAssign (int amrlev, int fmglev, MultiFab& fine, MultiFab& crse) const override;

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

    virtual void fillRHS (MFIter const& /*mfi*/,
                          Array4<int const> const& /*lid*/,
                          Real* /*rhs*/,
                          Array4<Real const> const& /*bfab*/) const
    {
        amrex::Abort("MLNodeLinOp:fillRHS: how did we get here?");
    }
#endif

protected:

    void resizeMultiGrid (int new_size) override;

    std::unique_ptr<iMultiFab> m_owner_mask_top;  // ownership of nodes
    std::unique_ptr<iMultiFab> m_owner_mask_bottom;
    Vector<Vector<std::unique_ptr<iMultiFab> > > m_dirichlet_mask;  // dirichlet?
    Vector<std::unique_ptr<iMultiFab> > m_cc_fine_mask;          // cell-centered mask for cells covered by fine
    Vector<std::unique_ptr<iMultiFab> > m_nd_fine_mask;          // nodal mask: 0: this level node, 1: c/f boundary, 2: fine node
    Vector<std::unique_ptr<LayoutData<int> > > m_has_fine_bndry; // does this fab contain c/f boundary?
    MultiFab m_bottom_dot_mask;
    MultiFab m_coarse_dot_mask;
    Vector<MultiFab> m_precond_weight_mask;

    Vector<std::unique_ptr<iMultiFab> > m_norm_fine_mask;

#ifdef AMREX_USE_EB
    CoarseningStrategy m_coarsening_strategy = CoarseningStrategy::RAP;
#else
    CoarseningStrategy m_coarsening_strategy = CoarseningStrategy::Sigma;
#endif

    bool m_masks_built = false;
    bool m_overset_dirichlet_mask = false;
#ifdef AMREX_USE_GPU
    int m_smooth_num_sweeps = 4;
#else
    int m_smooth_num_sweeps = 2;
#endif
    mutable bool m_in_solution_mode = true;

private:
    bool m_is_bottom_singular = false;
};

}

#endif
