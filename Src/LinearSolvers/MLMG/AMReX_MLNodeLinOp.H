#ifndef AMREX_ML_NODE_LINOP_H_H
#define AMREX_ML_NODE_LINOP_H_H

#include <AMReX_MLLinOp.H>
#include <AMReX_iMultiFab.H>

#ifdef AMREX_USE_HYPRE
#include <AMReX_HypreNodeLap.H>
#endif

namespace amrex {

class MLNodeLinOp
    : public MLLinOp
{
public:

    friend class MLMG;
    friend class MLCGSolver;

    enum struct CoarseningStrategy : int { Sigma, RAP };

    MLNodeLinOp ();
    virtual ~MLNodeLinOp ();

    MLNodeLinOp (const MLNodeLinOp&) = delete;
    MLNodeLinOp (MLNodeLinOp&&) = delete;
    MLNodeLinOp& operator= (const MLNodeLinOp&) = delete;
    MLNodeLinOp& operator= (MLNodeLinOp&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    virtual void setLevelBC (int /*amrlev*/, const MultiFab* /*levelbcdata*/) final override {}

    virtual void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                        StateMode s_mode, const MLMGBndry* bndry=nullptr) const final override;

    virtual void smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                         bool skip_fillboundary=false) const final override;

    virtual void solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                   const MultiFab* crse_bcdata=nullptr) override;
    virtual void correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                     BCMode bc_mode, const MultiFab* crse_bcdata=nullptr) override;
    virtual void compFlux (int /*amrlev*/, const Array<MultiFab*,AMREX_SPACEDIM>& /*fluxes*/,
                           MultiFab& /*sol*/, Location /*loc*/) const final override {
        amrex::Abort("AMReX_MLNodeLinOp::compFlux::How did we get here?");
    }
    virtual void compGrad (int /*amrlev*/, const Array<MultiFab*,AMREX_SPACEDIM>& /*grad*/,
                           MultiFab& /*sol*/, Location /*loc*/) const final override {
        amrex::Abort("AMReX_MLNodeLinOp::compGrad::How did we get here?");
    }
    
    virtual void applyMetricTerm (int /*amrlev*/, int /*mglev*/, MultiFab& /*rhs*/) const final override {}
    virtual void unapplyMetricTerm (int /*amrlev*/, int /*mglev*/, MultiFab& /*rhs*/) const final override {}

    virtual void fillSolutionBC (int /*amrlev*/, MultiFab& /*sol*/,
                                 const MultiFab* /*crse_bcdata*/=nullptr) final override {
        amrex::Abort("AMReX_MLNodeLinOp::fillSolutionBC::How did we get here?");
    }

    virtual void applyInhomogNeumannTerm (int armlev, MultiFab& rhs) const override;

    virtual void prepareForSolve () override {}

    virtual bool isSingular (int amrlev) const override
        { return (amrlev == 0) ? m_is_bottom_singular : false; }
    virtual bool isBottomSingular () const override { return m_is_bottom_singular; }

    virtual Real xdoty (int amrlev, int mglev, const MultiFab& x, const MultiFab& y, bool local) const final override;

    virtual void applyBC (int amrlev, int mglev, MultiFab& phi, BCMode bc_mode, StateMode s_mode,
                          bool skip_fillboundary=false) const;

    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const = 0;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh) const = 0;

    virtual void nodalSync (int amrlev, int mglev, MultiFab& mf) const final override;

    virtual std::unique_ptr<MLLinOp> makeNLinOp (int /*grid_size*/) const final override {
        amrex::Abort("MLNodeLinOp::makeNLinOp: N-Solve not supported");
        return std::unique_ptr<MLLinOp>{};
    }

    static std::unique_ptr<iMultiFab> makeOwnerMask (const BoxArray& ba,
                                                     const DistributionMapping& dm,
                                                     const Geometry& geom);

    void buildMasks ();

    // omask is either 0 or 1. 1 means the node is an unknown. 0 means it's known.
    void setOversetMask (int amrlev, const iMultiFab& a_omask);

#ifdef AMREX_USE_HYPRE
    virtual std::unique_ptr<HypreNodeLap> makeHypreNodeLap (int bottom_verbose) const override;

    virtual void fillIJMatrix (MFIter const& mfi, Array4<HypreNodeLap::Int const> const& nid,
                               Array4<int const> const& owner,
                               Vector<HypreNodeLap::Int>& ncols, Vector<HypreNodeLap::Int>& rows,
                               Vector<HypreNodeLap::Int>& cols, Vector<Real>& mat) const
    {
        amrex::Abort("MLNodeLinOp: how did we get here?");
    }
#endif

protected:

    Vector<Vector<std::unique_ptr<iMultiFab> > > m_owner_mask;      // ownership of nodes
    Vector<Vector<std::unique_ptr<iMultiFab> > > m_dirichlet_mask;  // dirichlet?
    Vector<std::unique_ptr<iMultiFab> > m_cc_fine_mask;          // cell-centered mask for cells covered by fine
    Vector<std::unique_ptr<iMultiFab> > m_nd_fine_mask;          // nodal mask: 0: this level node, 1: c/f boundary, 2: fine node
    Vector<std::unique_ptr<LayoutData<int> > > m_has_fine_bndry; // does this fab contain c/f boundary?
    MultiFab m_bottom_dot_mask;
    MultiFab m_coarse_dot_mask;

#ifdef AMREX_USE_EB
    CoarseningStrategy m_coarsening_strategy = CoarseningStrategy::RAP;
#else
    CoarseningStrategy m_coarsening_strategy = CoarseningStrategy::Sigma;
#endif

    bool m_is_bottom_singular = false;
    bool m_masks_built = false;
    bool m_overset_dirichlet_mask = false;
};

}

#endif
