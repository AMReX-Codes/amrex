#ifndef AMREX_ML_NODE_LINOP_H_H
#define AMREX_ML_NODE_LINOP_H_H

#include <AMReX_MLLinOp.H>

namespace amrex {

class MLNodeLinOp
    : public MLLinOp
{
public:

    friend class MLMG;
    friend class MLCGSolver;

    MLNodeLinOp ();
    virtual ~MLNodeLinOp ();

    MLNodeLinOp (const MLNodeLinOp&) = delete;
    MLNodeLinOp (MLNodeLinOp&&) = delete;
    MLNodeLinOp& operator= (const MLNodeLinOp&) = delete;
    MLNodeLinOp& operator= (MLNodeLinOp&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo());

    virtual void setLevelBC (int amrlev, const MultiFab* levelbcdata) final {}

protected:

    virtual void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                        const MLMGBndry* bndry=nullptr) const final {
        amrex::Abort("TODO: apply");
    }
    virtual void smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                         bool skip_fillboundary=false) const final {
        amrex::Abort("TODO: smooth");
    }

    virtual void solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                   const MultiFab* crse_bcdata=nullptr) final {
        amrex::Abort("TODO: solutionResidual");
    }
    virtual void correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                     BCMode bc_mode, const MultiFab* crse_bcdata=nullptr) final {
        amrex::Abort("TODO: correctionResidual");
    }

    virtual void reflux (int crse_amrlev, MultiFab& res, const MultiFab& crse_sol, MultiFab& fine_sol) const final {
        amrex::Abort("TODO: reflux");
    }
    virtual void compFlux (int amrlev, const std::array<MultiFab*,AMREX_SPACEDIM>& fluxes, MultiFab& sol) const final {
        amrex::Abort("TODO: compFlux");
    }
    virtual void compGrad (int amrlev, const std::array<MultiFab*,AMREX_SPACEDIM>& grad, MultiFab& sol) const final {
        amrex::Abort("TODO: compGrad");
    }
    
    virtual void applyMetricTerm (int amrlev, int mglev, MultiFab& rhs) const final {}
    virtual void unapplyMetricTerm (int amrlev, int mglev, MultiFab& rhs) const final {}

    virtual void fillSolutionBC (int amrlev, MultiFab& sol, const MultiFab* crse_bcdata=nullptr) final {
        amrex::Abort("TODO: fillSolutionBC");
    }

    virtual void prepareForSolve () override {}
};

}

#endif
