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

template <typename MF = amrex::MultiFab>
class MLNodeLinOpT
    : public MLLinOpT<MF>
{
public:

    using FAB        = typename MF::fab_type;
    using RT         = typename MF::value_type;

    using BCType     = LinOpBCType;
    using BCMode     = typename MLLinOpT<MF>::BCMode;
    using StateMode  = typename MLLinOpT<MF>::StateMode;
    using Location   = typename MLLinOpT<MF>::Location;

    enum struct CoarseningStrategy : int { Sigma, RAP };

    MLNodeLinOpT ();
    ~MLNodeLinOpT () override = default;

    MLNodeLinOpT (const MLNodeLinOpT<MF>&) = delete;
    MLNodeLinOpT (MLNodeLinOpT<MF>&&) = delete;
    MLNodeLinOpT<MF>& operator= (const MLNodeLinOpT<MF>&) = delete;
    MLNodeLinOpT<MF>& operator= (MLNodeLinOpT<MF>&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FAB> const*>& a_factory = {},
                 int a_eb_limit_coarsening = -1);

    void setSmoothNumSweeps (int nsweeps) noexcept {
        m_smooth_num_sweeps = nsweeps;
    }

    void setLevelBC (int /*amrlev*/, const MF* /*levelbcdata*/,
                     const MF* = nullptr, const MF* = nullptr,
                     const MF* = nullptr) final {}

    void apply (int amrlev, int mglev, MF& out, MF& in, BCMode bc_mode,
                        StateMode s_mode, const MLMGBndryT<MF>* bndry=nullptr) const final;

    void smooth (int amrlev, int mglev, MF& sol, const MF& rhs,
                 bool skip_fillboundary, int niter) const override;

    void solutionResidual (int amrlev, MF& resid, MF& x, const MF& b,
                                   const MF* crse_bcdata=nullptr) override;
    void correctionResidual (int amrlev, int mglev, MF& resid, MF& x, const MF& b,
                                     BCMode bc_mode, const MF* crse_bcdata=nullptr) override;

    Vector<RT> getSolvabilityOffset (int amrlev, int mglev,
                                               MF const& rhs) const override;
    void fixSolvabilityByOffset (int amrlev, int mglev, MF& rhs,
                                         Vector<RT> const& offset) const override;

    void prepareForSolve () override;

    void preparePrecond () override;

    void setDirichletNodesToZero (int amrlev, int mglev, MF& mf) const override;

    bool isSingular (int amrlev) const override
        { return (amrlev == 0) ? m_is_bottom_singular : false; }
    bool isBottomSingular () const override { return m_is_bottom_singular; }

    RT xdoty (int amrlev, int mglev, const MF& x, const MF& y, bool local) const final;

    RT dotProductPrecond (Vector<MF const*> const& x,
                            Vector<MF const*> const& y) const final;

    RT norm2Precond (Vector<MF const*> const& x) const final;

    virtual void applyBC (int amrlev, int mglev, MF& phi, BCMode bc_mode, StateMode state_mode,
                          bool skip_fillboundary=false) const;

    virtual void Fapply (int amrlev, int mglev, MF& out, const MF& in) const = 0;
    virtual void Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs) const = 0;

    void nodalSync (int amrlev, int mglev, MF& mf) const;

    static std::unique_ptr<iMultiFab> makeOwnerMask (const BoxArray& ba,
                                                     const DistributionMapping& dm,
                                                     const Geometry& geom);

    void buildMasks ();

    // omask is either 0 or 1. 1 means the node is an unknown. 0 means it's known.
    void setOversetMask (int amrlev, const iMultiFab& a_dmask);

    virtual void fixUpResidualMask (int /*amrlev*/, iMultiFab& /*resmsk*/) { }

    RT normInf (int amrlev, MF const& mf, bool local) const override;

    void avgDownResAmr (int, MF&, MF const&) const final { }

    void interpolationAmr (int famrlev, MF& fine, const MF& crse,
                                   IntVect const& nghost) const override;

    void averageDownAndSync (Vector<MF>& sol) const override;

    void interpAssign (int amrlev, int fmglev, MF& fine, MF& crse) const override;

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
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
        amrex::Abort("MLNodeLinOpT::fillIJMatrix: how did we get here?");
    }

    virtual void fillRHS (MFIter const& /*mfi*/,
                          Array4<int const> const& /*lid*/,
                          Real* /*rhs*/,
                          Array4<Real const> const& /*bfab*/) const
    {
        amrex::Abort("MLNodeLinOpT:fillRHS: how did we get here?");
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
#ifdef AMREX_USE_GPU
    int m_smooth_num_sweeps = 4;
#else
    int m_smooth_num_sweeps = 2;
#endif
    mutable bool m_in_solution_mode = true;

private:
    bool m_is_bottom_singular = false;
};

using MLNodeLinOp = MLNodeLinOpT<>;

}

#endif
