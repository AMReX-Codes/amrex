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

    virtual void setLevelBC (int amrlev, const MultiFab* levelbcdata) final override;

    virtual bool needsUpdate () const override {
        return MLLinOp::needsUpdate();
    }
    virtual void update () override;

#ifdef AMREX_SOFT_PERF_COUNTERS
    struct Counters
    {
        size_t smooth_n = 0;
        size_t smooth_cell_n = 0;
        size_t apply_n = 0;
        size_t apply_cell_n = 0;
        size_t restrict_n = 0;
        size_t restrict_cell_n = 0;
        size_t interp_n = 0;
        size_t interp_cell_n = 0;

        void reset () {
            smooth_n = 0;
            smooth_cell_n = 0;
            apply_n = 0;
            apply_cell_n = 0;
            restrict_n = 0;
            restrict_cell_n = 0;
            interp_n = 0;
            interp_cell_n = 0;
        }

        void smooth (const MultiFab & mf) {
            smooth_n++;
            smooth_cell_n += mf.boxArray().numPts();
        }
        void apply (const MultiFab & mf) {
            apply_n++;
            apply_cell_n += mf.boxArray().numPts();
        }
        void restrict (const MultiFab & mf) {
            restrict_n++;
            restrict_cell_n += mf.boxArray().numPts();
        }
        void interpolate (const MultiFab & mf) {
            interp_n++;
            interp_cell_n += mf.boxArray().numPts();
        }

        void print () {
            Print() << ParallelContext::MyProcSub() << ": perf: MLCellLinOp::smooth()        called " << smooth_n
                    << " times on " << smooth_cell_n << " cells total" << std::endl;
            Print() << ParallelContext::MyProcSub() << ": perf: MLCellLinOp::apply()         called " << apply_n
                    << " times on " << apply_cell_n << " cells total" << std::endl;
            Print() << ParallelContext::MyProcSub() << ": perf: MLCellLinOp::restriction()   called " << restrict_n
                    << " times on " << restrict_cell_n << " cells total" << std::endl;
            Print() << ParallelContext::MyProcSub() << ": perf: MLCellLinOp::interpolation() called " << interp_n
                    << " times on " << interp_cell_n << " cells total" << std::endl;
        }
    };

    // software counters for performance analysis
    static Counters perf_counters;
#endif

    virtual bool isCrossStencil () const { return true; }

    void updateSolBC (int amrlev, const MultiFab& crse_bcdata) const;
    void updateCorBC (int amrlev, const MultiFab& crse_bcdata) const;

    virtual void applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode, StateMode s_mode,
                          const MLMGBndry* bndry=nullptr, bool skip_fillboundary=false) const;

    BoxArray makeNGrids (int grid_size) const;

    virtual void restriction (int, int, MultiFab& crse, MultiFab& fine) const override;

    virtual void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const override;

    virtual void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs) override;

    virtual void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                        StateMode s_mode, const MLMGBndry* bndry=nullptr) const override;
    virtual void smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                         bool skip_fillboundary=false) const final override;

    virtual void solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                   const MultiFab* crse_bcdata=nullptr) override;

    virtual void correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                     BCMode bc_mode, const MultiFab* crse_bcdata=nullptr) final override;

    // The assumption is crse_sol's boundary has been filled, but not fine_sol.
    virtual void reflux (int crse_amrlev,
                         MultiFab& res, const MultiFab& crse_sol, const MultiFab&,
                         MultiFab&, MultiFab& fine_sol, const MultiFab&) const final override;
    virtual void compFlux (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes,
                           MultiFab& sol, Location loc) const final override;
    virtual void compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                           MultiFab& sol, Location loc) const override;

    virtual void applyMetricTerm (int amrlev, int mglev, MultiFab& rhs) const final override;
    virtual void unapplyMetricTerm (int amrlev, int mglev, MultiFab& rhs) const final override;
    virtual void fillSolutionBC (int amrlev, MultiFab& sol, const MultiFab* crse_bcdata=nullptr) final override;

    virtual void prepareForSolve () override;

    virtual Real xdoty (int amrlev, int mglev, const MultiFab& x, const MultiFab& y, bool local) const final override;

    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const = 0;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const = 0;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, Location loc, const int face_only=0) const = 0;

protected:

#if (AMREX_SPACEDIM != 3)
    struct MetricFactor {
        MetricFactor (const BoxArray& ba, const DistributionMapping& dm,
                      const Geometry& geom, bool null_metric);
        const Vector<Real>& cellCenters (const MFIter& mfi) const noexcept {
            return r_cellcenter[mfi];
        }
        const Vector<Real>& cellEdges (const MFIter& mfi) const noexcept {
            return r_celledge[mfi];
        }
        const Vector<Real>& invCellCenters (const MFIter& mfi) const noexcept {
            return inv_r_cellcenter[mfi];
        }
        const Vector<Real>& invCellEdges (const MFIter& mfi) const noexcept {
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
    using RealTuple = Array<Real,2*BL_SPACEDIM>;
    using BCTuple   = Array<BoundCond,2*BL_SPACEDIM>;
    class BndryCondLoc
    {
    public:
        BndryCondLoc (const BoxArray& ba, const DistributionMapping& dm, int ncomp);
        void setLOBndryConds (const Geometry& geom, const Real* dx,
                              const Vector<Array<BCType,AMREX_SPACEDIM> >& lobc,
                              const Vector<Array<BCType,AMREX_SPACEDIM> >& hibc,
                              int ratio, const RealVect& a_loc);
        const Vector<BCTuple>& bndryConds (const MFIter& mfi) const noexcept {
            return bcond[mfi];
        }
        const Vector<RealTuple>& bndryLocs (const MFIter& mfi) const noexcept {
            return bcloc[mfi];
        }
        const BCTuple& bndryConds (const MFIter& mfi, int icomp) const noexcept {
            return bcond[mfi][icomp];
        }
        const RealTuple& bndryLocs (const MFIter& mfi, int icomp) const noexcept {
            return bcloc[mfi][icomp];
        }
    private:
        LayoutData<Vector<BCTuple> >   bcond;
        LayoutData<Vector<RealTuple> > bcloc;
        int m_ncomp;
    };
    Vector<Vector<std::unique_ptr<BndryCondLoc> > > m_bcondloc;

    // used to save interpolation coefficients of the first interior cells
    mutable Vector<Vector<BndryRegister> > m_undrrelxr;

    // boundary cell flags for covered, not_covered, outside_domain
    Vector<Vector<Array<MultiMask,2*AMREX_SPACEDIM> > > m_maskvals;

    mutable Vector<YAFluxRegister> m_fluxreg;

private:

    void defineAuxData ();
    void defineBC ();

};

}

#endif
