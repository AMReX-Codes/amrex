#ifndef AMREX_ML_NODE_LINOP_H_H
#define AMREX_ML_NODE_LINOP_H_H

#include <AMReX_MLLinOp.H>
#include <AMReX_iMultiFab.H>

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
                        const MLMGBndry* bndry=nullptr) const final;

    virtual void smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                         bool skip_fillboundary=false) const final;

    virtual void solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                   const MultiFab* crse_bcdata=nullptr) final;
    virtual void correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                     BCMode bc_mode, const MultiFab* crse_bcdata=nullptr) final;
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

    virtual void applyBC (int amrlev, int mglev, MultiFab& phi, BCMode bc_mode,
                          bool skip_fillboundary=false) const = 0;
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const = 0;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh) const = 0;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only=0) const = 0;

    virtual iMultiFab const* getOwnerMask (int amrlev, int mglev) const final {
        return m_owner_mask[amrlev][mglev].get();
    }

    void nodalSync (int amrlev, int mglev, MultiFab& mf) const;

    static std::unique_ptr<iMultiFab> makeOwnerMask (const BoxArray& ba,
                                                     const DistributionMapping& dm,
                                                     const Geometry& geom);

    Vector<Vector<std::unique_ptr<iMultiFab> > > m_owner_mask;      // ownership of nodes
    Vector<Vector<std::unique_ptr<iMultiFab> > > m_dirichlet_mask;  // dirichlet?
    Vector<std::unique_ptr<iMultiFab> > m_crsefine_mask; // cell-cenetered mask for crse(0)/fine(1)
    Vector<std::unique_ptr<LayoutData<int> > > m_has_fine_bndry; // does this fab contains c/f boundary?
};

}

#endif
