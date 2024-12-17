#ifndef AMREX_ML_NODE_LINOP_H_H
#define AMREX_ML_NODE_LINOP_H_H
#include <AMReX_Config.H>

#include <AMReX_MLLinOp.H>
#include <AMReX_iMultiFab.H>

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
#include <AMReX_HypreNodeLap.H>
#endif

namespace amrex {

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

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {},
                 int a_eb_limit_coarsening = -1);

    void setSmoothNumSweeps (int nsweeps) noexcept {
        m_smooth_num_sweeps = nsweeps;
    }

    void setLevelBC (int /*amrlev*/, const MultiFab* /*levelbcdata*/,
                     const MultiFab* = nullptr, const MultiFab* = nullptr,
                     const MultiFab* = nullptr) final {}

    void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                        StateMode s_mode, const MLMGBndry* bndry=nullptr) const final;

    void smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                         bool skip_fillboundary=false) const override;

    void solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                   const MultiFab* crse_bcdata=nullptr) override;
    void correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                     BCMode bc_mode, const MultiFab* crse_bcdata=nullptr) override;

    Vector<Real> getSolvabilityOffset (int amrlev, int mglev,
                                               MultiFab const& rhs) const override;
    void fixSolvabilityByOffset (int amrlev, int mglev, MultiFab& rhs,
                                         Vector<Real> const& offset) const override;

    void prepareForSolve () override;

    void preparePrecond () override;

    void setDirichletNodesToZero (int amrlev, int mglev, MultiFab& mf) const override;

    bool isSingular (int amrlev) const override
        { return (amrlev == 0) ? m_is_bottom_singular : false; }
    bool isBottomSingular () const override { return m_is_bottom_singular; }

    Real xdoty (int amrlev, int mglev, const MultiFab& x, const MultiFab& y, bool local) const final;

    Real dotProductPrecond (Vector<MultiFab const*> const& x,
                            Vector<MultiFab const*> const& y) const final;

    Real norm2Precond (Vector<MultiFab const*> const& x) const final;

    virtual void applyBC (int amrlev, int mglev, MultiFab& phi, BCMode bc_mode, StateMode state_mode,
                          bool skip_fillboundary=false) const;

    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const = 0;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const = 0;

    void nodalSync (int amrlev, int mglev, MultiFab& mf) const;

    static std::unique_ptr<iMultiFab> makeOwnerMask (const BoxArray& ba,
                                                     const DistributionMapping& dm,
                                                     const Geometry& geom);

    void buildMasks ();

    // omask is either 0 or 1. 1 means the node is an unknown. 0 means it's known.
    void setOversetMask (int amrlev, const iMultiFab& a_dmask);

    virtual void fixUpResidualMask (int /*amrlev*/, iMultiFab& /*resmsk*/) { }

    Real normInf (int amrlev, MultiFab const& mf, bool local) const override;

    void avgDownResAmr (int, MultiFab&, MultiFab const&) const final { }

    void interpolationAmr (int famrlev, MultiFab& fine, const MultiFab& crse,
                                   IntVect const& nghost) const override;

    void averageDownAndSync (Vector<MultiFab>& sol) const override;

    void interpAssign (int amrlev, int fmglev, MultiFab& fine, MultiFab& crse) const override;

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
