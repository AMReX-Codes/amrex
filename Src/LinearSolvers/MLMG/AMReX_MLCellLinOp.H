#ifndef AMREX_ML_CELL_LINOP_H_H
#define AMREX_ML_CELL_LINOP_H_H

#include <AMReX_MLLinOp.H>

namespace amrex {

class MLCellLinOp
    : public MLLinOp
{
public:

    friend class MLMG;
    friend class MLCGSolver;

    MLCellLinOp ();
    virtual ~MLCellLinOp ();

    MLCellLinOp (const MLCellLinOp&) = delete;
    MLCellLinOp (MLCellLinOp&&) = delete;
    MLCellLinOp& operator= (const MLCellLinOp&) = delete;
    MLCellLinOp& operator= (MLCellLinOp&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    virtual void setLevelBC (int amrlev, const MultiFab* levelbcdata) final;

protected:

#if (AMREX_SPACEDIM != 3)
    struct MetricFactor {
        MetricFactor (const BoxArray& ba, const DistributionMapping& dm,
                      const Geometry& geom, bool null_metric);
        const Vector<Real>& cellCenters (const MFIter& mfi) const {
            return r_cellcenter[mfi];
        }
        const Vector<Real>& cellEdges (const MFIter& mfi) const {
            return r_celledge[mfi];
        }
        const Vector<Real>& invCellCenters (const MFIter& mfi) const {
            return inv_r_cellcenter[mfi];
        }
        const Vector<Real>& invCellEdges (const MFIter& mfi) const {
            return inv_r_celledge[mfi];
        }
        LayoutData<Vector<Real> > r_cellcenter;
        LayoutData<Vector<Real> > r_celledge;
        LayoutData<Vector<Real> > inv_r_cellcenter;
        LayoutData<Vector<Real> > inv_r_celledge;
    };
    Vector<Vector<std::unique_ptr<MetricFactor> > > m_metric_factor;
#endif

    Vector<std::unique_ptr<MLMGBndry> >   m_bndry_sol;
    Vector<std::unique_ptr<BndryRegister> > m_crse_sol_br;

    Vector<std::unique_ptr<MLMGBndry> > m_bndry_cor;
    Vector<std::unique_ptr<BndryRegister> > m_crse_cor_br;

    // In case of agglomeration, coarse MG grids on amr level 0 are
    // not simply coarsened from fine MG grids.  So we need to build
    // bcond and bcloc for each MG level.
    using RealTuple = std::array<Real,2*BL_SPACEDIM>;
    using BCTuple   = std::array<BoundCond,2*BL_SPACEDIM>;
    class BndryCondLoc
    {
    public:
        BndryCondLoc (const BoxArray& ba, const DistributionMapping& dm);
        void setLOBndryConds (const Geometry& geom, const Real* dx,
                              const std::array<BCType,AMREX_SPACEDIM>& lobc,
                              const std::array<BCType,AMREX_SPACEDIM>& hibc,
                              int ratio, const RealVect& a_loc);
        const BCTuple& bndryConds (const MFIter& mfi) const {
            return bcond[mfi];
        }
        const RealTuple& bndryLocs (const MFIter& mfi) const {
            return bcloc[mfi];
        }
    private:
        LayoutData<BCTuple>   bcond;
        LayoutData<RealTuple> bcloc;
    };
    Vector<Vector<std::unique_ptr<BndryCondLoc> > > m_bcondloc;

    // used to save interpolation coefficients of the first interior cells
    mutable Vector<Vector<BndryRegister> > m_undrrelxr;

    // boundary cell flags for covered, not_covered, outside_domain
    Vector<Vector<std::array<MultiMask,2*AMREX_SPACEDIM> > > m_maskvals;

    mutable Vector<YAFluxRegister> m_fluxreg;

    //
    // functions
    //

    void updateSolBC (int amrlev, const MultiFab& crse_bcdata) const;
    void updateCorBC (int amrlev, const MultiFab& crse_bcdata) const;

    void applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode,
                  const MLMGBndry* bndry=nullptr, bool skip_fillboundary=false) const;

    BoxArray makeNGrids (int grid_size) const;

    virtual void restriction (int, int, MultiFab& crse, MultiFab& fine) const final;

    virtual void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final;

    virtual void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs) final;

    virtual void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                        const MLMGBndry* bndry=nullptr) const final;
    virtual void smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                         bool skip_fillboundary=false) const final;

    virtual void solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                   const MultiFab* crse_bcdata=nullptr) final;

    virtual void correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                     BCMode bc_mode, const MultiFab* crse_bcdata=nullptr) final;

    // The assumption is crse_sol's boundary has been filled, but not fine_sol.
    virtual void reflux (int crse_amrlev,
                         MultiFab& res, const MultiFab& crse_sol, const MultiFab&,
                         MultiFab&, MultiFab& fine_sol, const MultiFab&) const final;
    virtual void compFlux (int amrlev, const std::array<MultiFab*,AMREX_SPACEDIM>& fluxes, MultiFab& sol) const final;
    virtual void compGrad (int amrlev, const std::array<MultiFab*,AMREX_SPACEDIM>& grad, MultiFab& sol) const final;

    virtual void applyMetricTerm (int amrlev, int mglev, MultiFab& rhs) const final;
    virtual void unapplyMetricTerm (int amrlev, int mglev, MultiFab& rhs) const final;
    virtual void fillSolutionBC (int amrlev, MultiFab& sol, const MultiFab* crse_bcdata=nullptr) final;

    virtual void prepareForSolve () override;

    virtual Real xdoty (int amrlev, int mglev, const MultiFab& x, const MultiFab& y, bool local) const final;

    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const = 0;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const = 0;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only=0) const = 0;
    virtual int getNComp() const;

private:

    void defineAuxData ();
    void defineBC ();

};

}

#endif
