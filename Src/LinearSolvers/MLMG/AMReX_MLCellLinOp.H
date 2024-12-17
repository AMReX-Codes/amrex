#ifndef AMREX_ML_CELL_LINOP_H_
#define AMREX_ML_CELL_LINOP_H_
#include <AMReX_Config.H>

#include <AMReX_MLLinOp.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_YAFluxRegister.H>
#include <AMReX_MLLinOp_K.H>
#include <AMReX_MLMG_K.H>

#ifndef BL_NO_FORT
#include <AMReX_MLLinOp_F.H>
#endif

namespace amrex {

template <typename MF>
class MLCellLinOpT  // NOLINT(cppcoreguidelines-virtual-class-destructor)
    : public MLLinOpT<MF>
{
public:

    using FAB = typename MF::fab_type;
    using RT  = typename MF::value_type;

    using BCType = LinOpBCType;
    using BCMode    = typename MLLinOpT<MF>::BCMode;
    using StateMode = typename MLLinOpT<MF>::StateMode;
    using Location  = typename MLLinOpT<MF>::Location;

    MLCellLinOpT ();
    ~MLCellLinOpT () override = default;

    MLCellLinOpT (const MLCellLinOpT<MF>&) = delete;
    MLCellLinOpT (MLCellLinOpT<MF>&&) = delete;
    MLCellLinOpT<MF>& operator= (const MLCellLinOpT<MF>&) = delete;
    MLCellLinOpT<MF>& operator= (MLCellLinOpT<MF>&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FAB> const*>& a_factory = {});

    void setLevelBC (int amrlev, const MF* levelbcdata,
                             const MF* robinbc_a = nullptr,
                             const MF* robinbc_b = nullptr,
                             const MF* robinbc_f = nullptr) final;

    using MLLinOpT<MF>::setLevelBC;

    bool needsUpdate () const override {
        return MLLinOpT<MF>::needsUpdate();
    }
    void update () override;

    void setGaussSeidel (bool flag) noexcept { m_use_gauss_seidel = flag; }

    virtual bool isCrossStencil () const { return true; }
    virtual bool isTensorOp () const { return false; }

    void updateSolBC (int amrlev, const MF& crse_bcdata) const;
    void updateCorBC (int amrlev, const MF& crse_bcdata) const;

    virtual void applyBC (int amrlev, int mglev, MF& in, BCMode bc_mode, StateMode s_mode,
                          const MLMGBndryT<MF>* bndry=nullptr, bool skip_fillboundary=false) const;

    BoxArray makeNGrids (int grid_size) const;

    void restriction (int, int, MF& crse, MF& fine) const override;

    void interpolation (int amrlev, int fmglev, MF& fine, const MF& crse) const override;

    void interpAssign (int amrlev, int fmglev, MF& fine, MF& crse) const override;

    void interpolationAmr (int famrlev, MF& fine, const MF& crse,
                                   IntVect const& nghost) const override;

    void averageDownSolutionRHS (int camrlev, MF& crse_sol, MF& crse_rhs,
                                         const MF& fine_sol, const MF& fine_rhs) override;

    void apply (int amrlev, int mglev, MF& out, MF& in, BCMode bc_mode,
                        StateMode s_mode, const MLMGBndryT<MF>* bndry=nullptr) const override;
    void smooth (int amrlev, int mglev, MF& sol, const MF& rhs,
                         bool skip_fillboundary=false) const final;

    void solutionResidual (int amrlev, MF& resid, MF& x, const MF& b,
                                   const MF* crse_bcdata=nullptr) override;

    void prepareForFluxes (int amrlev, const MF* crse_bcdata = nullptr) override;

    void correctionResidual (int amrlev, int mglev, MF& resid, MF& x, const MF& b,
                                     BCMode bc_mode, const MF* crse_bcdata=nullptr) final;

    // The assumption is crse_sol's boundary has been filled, but not fine_sol.
    void reflux (int crse_amrlev,
                         MF& res, const MF& crse_sol, const MF&,
                         MF&, MF& fine_sol, const MF&) const final;
    void compFlux (int amrlev, const Array<MF*,AMREX_SPACEDIM>& fluxes,
                           MF& sol, Location loc) const override;
    void compGrad (int amrlev, const Array<MF*,AMREX_SPACEDIM>& grad,
                           MF& sol, Location loc) const override;

    void applyMetricTerm (int amrlev, int mglev, MF& rhs) const final;
    void unapplyMetricTerm (int amrlev, int mglev, MF& rhs) const final;
    Vector<RT> getSolvabilityOffset (int amrlev, int mglev,
                                             MF const& rhs) const override;
    void fixSolvabilityByOffset (int amrlev, int mglev, MF& rhs,
                                         Vector<RT> const& offset) const override;

    void prepareForSolve () override;

    RT xdoty (int amrlev, int mglev, const MF& x, const MF& y, bool local) const final;

    RT dotProductPrecond (Vector<MF const*> const& x, Vector<MF const*> const& y) const final;

    RT norm2Precond (Vector<MF const*> const& x) const final;

    virtual void Fapply (int amrlev, int mglev, MF& out, const MF& in) const = 0;
    virtual void Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs, int redblack) const = 0;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const Array<FAB*,AMREX_SPACEDIM>& flux,
                        const FAB& sol, Location loc, int face_only=0) const = 0;

    virtual void addInhomogNeumannFlux (int /*amrlev*/,
                                        const Array<MF*,AMREX_SPACEDIM>& /*grad*/,
                                        MF const& /*sol*/,
                                        bool /*mult_bcoef*/) const {}

    RT normInf (int amrlev, MF const& mf, bool local) const override;

    void averageDownAndSync (Vector<MF>& sol) const override;

    void avgDownResAmr (int clev, MF& cres, MF const& fres) const override;

    void beginPrecondBC () override;
    void endPrecondBC () override;

    struct BCTL {
        BoundCond type;
        RT location;
    };

    Vector<std::unique_ptr<MF> > m_robin_bcval;

#ifdef AMREX_USE_HYPRE
    void setInterpBndryHalfWidth (int w) { m_interpbndry_halfwidth = w; }
#endif

protected:

    bool m_has_metric_term = false;

    Vector<std::unique_ptr<MLMGBndryT<MF>> >   m_bndry_sol;
    Vector<std::unique_ptr<BndryRegisterT<MF>> > m_crse_sol_br;

    Vector<std::unique_ptr<MLMGBndryT<MF>> > m_bndry_cor;
    Vector<std::unique_ptr<BndryRegisterT<MF>> > m_crse_cor_br;

    Vector<std::unique_ptr<MLMGBndryT<MF>> > m_bndry_sol_zero;

    // In case of agglomeration, coarse MG grids on amr level 0 are
    // not simply coarsened from fine MG grids.  So we need to build
    // bcond and bcloc for each MG level.
    using RealTuple = Array<RT,2*BL_SPACEDIM>;
    using BCTuple   = Array<BoundCond,2*BL_SPACEDIM>;
    class BndryCondLoc
    {
    public:
        BndryCondLoc (const BoxArray& ba, const DistributionMapping& dm, int ncomp);

        void setLOBndryConds (const Geometry& geom, const Real* dx,
                              const Vector<Array<BCType,AMREX_SPACEDIM> >& lobc,
                              const Vector<Array<BCType,AMREX_SPACEDIM> >& hibc,
                              IntVect const& ratio, const RealVect& interior_bloc,
                              const Array<Real,AMREX_SPACEDIM>& domain_bloc_lo,
                              const Array<Real,AMREX_SPACEDIM>& domain_bloc_hi,
                              LinOpBCType crse_fine_bc_type);

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
        GpuArray<BCTL,2*AMREX_SPACEDIM> const* getBCTLPtr (const MFIter& mfi) const noexcept {
            return bctl[mfi];
        }
    private:
        LayoutData<Vector<BCTuple> >   bcond;
        LayoutData<Vector<RealTuple> > bcloc;
        LayoutData<GpuArray<BCTL,2*AMREX_SPACEDIM>*> bctl;
        Gpu::DeviceVector<GpuArray<BCTL,2*AMREX_SPACEDIM> > bctl_dv;
        int m_ncomp;
    };
    Vector<Vector<std::unique_ptr<BndryCondLoc> > > m_bcondloc;

    // used to save interpolation coefficients of the first interior cells
    mutable Vector<Vector<BndryRegisterT<MF>> > m_undrrelxr;

    // boundary cell flags for covered, not_covered, outside_domain
    Vector<Vector<Array<MultiMask,2*AMREX_SPACEDIM> > > m_maskvals;

    Vector<std::unique_ptr<iMultiFab> > m_norm_fine_mask;

    mutable Vector<YAFluxRegisterT<MF>> m_fluxreg;

    bool m_use_gauss_seidel = true; // use red-black Gauss-Seidel by default

private:

    void defineAuxData ();
    void defineBC ();

    void computeVolInv () const;
    mutable Vector<Vector<RT> > m_volinv; // used by solvability fix

    int m_interpbndry_halfwidth = 2;
};

template <typename T>
struct MLMGABCTag {
    Array4<T> fab;
    Array4<T const> bcval_lo;
    Array4<T const> bcval_hi;
    Array4<int const> mask_lo;
    Array4<int const> mask_hi;
    T bcloc_lo;
    T bcloc_hi;
    Box bx;
    BoundCond bctype_lo;
    BoundCond bctype_hi;
    int blen;
    int comp;
    int dir;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box() const noexcept { return bx; }
};

template <typename T>
struct MLMGPSTag {
    Array4<T> flo;
    Array4<T> fhi;
    Array4<int const> mlo;
    Array4<int const> mhi;
    T bcllo;
    T bclhi;
    Box bx;
    BoundCond bctlo;
    BoundCond bcthi;
    int blen;
    int comp;
    int dir;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box() const noexcept { return bx; }
};

#ifdef AMREX_USE_EB
template <typename T>
struct MLMGPSEBTag {
    Array4<T> flo;
    Array4<T> fhi;
    Array4<T const> ap;
    Array4<int const> mlo;
    Array4<int const> mhi;
    T bcllo;
    T bclhi;
    Box bx;
    BoundCond bctlo;
    BoundCond bcthi;
    int blen;
    int comp;
    int dir;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box() const noexcept { return bx; }
};
#endif

template <typename MF>
MLCellLinOpT<MF>::BndryCondLoc::BndryCondLoc (const BoxArray& ba,
                                              const DistributionMapping& dm,
                                              int ncomp)
    : bcond(ba, dm),
      bcloc(ba, dm),
      bctl(ba, dm),
      bctl_dv(bctl.local_size()*ncomp),
      m_ncomp(ncomp)
{
    auto* dp = bctl_dv.data();
    for (MFIter mfi(bcloc); mfi.isValid(); ++mfi) {
        bcond[mfi].resize(ncomp);
        bcloc[mfi].resize(ncomp);
        bctl[mfi] = dp;
        dp += ncomp;
    }
}

template <typename MF>
void
MLCellLinOpT<MF>::BndryCondLoc::
setLOBndryConds (const Geometry& geom, const Real* dx,
                 const Vector<Array<BCType,AMREX_SPACEDIM> >& lobc,
                 const Vector<Array<BCType,AMREX_SPACEDIM> >& hibc,
                 IntVect const& ratio, const RealVect& interior_bloc,
                 const Array<Real,AMREX_SPACEDIM>& domain_bloc_lo,
                 const Array<Real,AMREX_SPACEDIM>& domain_bloc_hi,
                 LinOpBCType crse_fine_bc_type)
{
    const Box& domain = geom.Domain();

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(bcloc); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        for (int icomp = 0; icomp < m_ncomp; ++icomp) {
            RealTuple & bloc  = bcloc[mfi][icomp];
            BCTuple   & bctag = bcond[mfi][icomp];
            MLMGBndryT<MF>::setBoxBC(bloc, bctag, bx, domain,
                                     lobc[icomp], hibc[icomp],
                                     dx, ratio, interior_bloc,
                                     domain_bloc_lo, domain_bloc_hi,
                                     geom.isPeriodicArray(),
                                     crse_fine_bc_type);
        }
    }

    Gpu::PinnedVector<GpuArray<BCTL,2*AMREX_SPACEDIM> > hv;
    hv.reserve(bctl_dv.size());
    for (MFIter mfi(bctl); mfi.isValid(); ++mfi)
    {
        for (int icomp = 0; icomp < m_ncomp; ++icomp) {
            GpuArray<BCTL,2*AMREX_SPACEDIM> tmp;
            for (int m = 0; m < 2*AMREX_SPACEDIM; ++m) {
                tmp[m].type = bcond[mfi][icomp][m];
                tmp[m].location = bcloc[mfi][icomp][m];
            }
            hv.push_back(std::move(tmp));
        }
    }
    Gpu::copyAsync(Gpu::hostToDevice, hv.begin(), hv.end(), bctl_dv.begin());
    Gpu::streamSynchronize();
}

template <typename MF>
MLCellLinOpT<MF>::MLCellLinOpT ()
{
    this->m_ixtype = IntVect::TheCellVector();
}

template <typename MF>
void
MLCellLinOpT<MF>::define (const Vector<Geometry>& a_geom,
                          const Vector<BoxArray>& a_grids,
                          const Vector<DistributionMapping>& a_dmap,
                          const LPInfo& a_info,
                          const Vector<FabFactory<FAB> const*>& a_factory)
{
    MLLinOpT<MF>::define(a_geom, a_grids, a_dmap, a_info, a_factory);
    defineAuxData();
    defineBC();
}

template <typename MF>
void
MLCellLinOpT<MF>::defineAuxData ()
{
    BL_PROFILE("MLCellLinOp::defineAuxData()");

    m_undrrelxr.resize(this->m_num_amr_levels);
    m_maskvals.resize(this->m_num_amr_levels);
    m_fluxreg.resize(this->m_num_amr_levels-1);
    m_norm_fine_mask.resize(this->m_num_amr_levels-1);

    const int ncomp = this->getNComp();

    for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev)
    {
        m_undrrelxr[amrlev].resize(this->m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < this->m_num_mg_levels[amrlev]; ++mglev)
        {
            m_undrrelxr[amrlev][mglev].define(this->m_grids[amrlev][mglev],
                                              this->m_dmap[amrlev][mglev],
                                              1, 0, 0, ncomp);
        }
    }

    for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev)
    {
        m_maskvals[amrlev].resize(this->m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < this->m_num_mg_levels[amrlev]; ++mglev)
        {
            for (OrientationIter oitr; oitr; ++oitr)
            {
                const Orientation face = oitr();
                const int ngrow = 1;
                const int extent = this->isCrossStencil() ? 0 : 1; // extend to corners
                m_maskvals[amrlev][mglev][face].define(this->m_grids[amrlev][mglev],
                                                       this->m_dmap[amrlev][mglev],
                                                       this->m_geom[amrlev][mglev],
                                                       face, 0, ngrow, extent, 1, true);
            }
        }
    }

    for (int amrlev = 0; amrlev < this->m_num_amr_levels-1; ++amrlev)
    {
        const IntVect ratio{this->m_amr_ref_ratio[amrlev]};
        m_fluxreg[amrlev].define(this->m_grids[amrlev+1][0],
                                 this->m_grids[amrlev][0],
                                 this->m_dmap[amrlev+1][0],
                                 this->m_dmap[amrlev][0],
                                 this->m_geom[amrlev+1][0],
                                 this->m_geom[amrlev][0],
                                 ratio, amrlev+1, ncomp);
        m_norm_fine_mask[amrlev] = std::make_unique<iMultiFab>
            (makeFineMask(this->m_grids[amrlev][0], this->m_dmap[amrlev][0],
                          this->m_grids[amrlev+1][0],
                          ratio, 1, 0));
    }

#if (AMREX_SPACEDIM != 3)
    m_has_metric_term = !this->m_geom[0][0].IsCartesian() && this->info.has_metric_term;
#endif
}

template <typename MF>
void
MLCellLinOpT<MF>::defineBC ()
{
    BL_PROFILE("MLCellLinOp::defineBC()");

    const int ncomp = this->getNComp();

    m_bndry_sol.resize(this->m_num_amr_levels);
    m_crse_sol_br.resize(this->m_num_amr_levels);

    m_bndry_cor.resize(this->m_num_amr_levels);
    m_crse_cor_br.resize(this->m_num_amr_levels);

    m_robin_bcval.resize(this->m_num_amr_levels);

    for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev)
    {
        m_bndry_sol[amrlev] = std::make_unique<MLMGBndryT<MF>>(this->m_grids[amrlev][0],
                                                               this->m_dmap[amrlev][0],
                                                               ncomp,
                                                               this->m_geom[amrlev][0]);
    }

    for (int amrlev = 1; amrlev < this->m_num_amr_levels; ++amrlev)
    {
        const int in_rad = 0;
        const int out_rad = 1;
        const int extent_rad = 2;
        const int crse_ratio = this->m_amr_ref_ratio[amrlev-1];
        BoxArray cba = this->m_grids[amrlev][0];
        cba.coarsen(crse_ratio);
        m_crse_sol_br[amrlev] = std::make_unique<BndryRegisterT<MF>>
            (cba, this->m_dmap[amrlev][0], in_rad, out_rad, extent_rad, ncomp);
    }

    for (int amrlev = 1; amrlev < this->m_num_amr_levels; ++amrlev)
    {
        const int in_rad = 0;
        const int out_rad = 1;
        const int extent_rad = 2;
        const int crse_ratio = this->m_amr_ref_ratio[amrlev-1];
        BoxArray cba = this->m_grids[amrlev][0];
        cba.coarsen(crse_ratio);
        m_crse_cor_br[amrlev] = std::make_unique<BndryRegisterT<MF>>
            (cba, this->m_dmap[amrlev][0], in_rad, out_rad, extent_rad, ncomp);
        m_crse_cor_br[amrlev]->setVal(RT(0.0));
    }

    // This has be to done after m_crse_cor_br is defined.
    for (int amrlev = 1; amrlev < this->m_num_amr_levels; ++amrlev)
    {
        m_bndry_cor[amrlev] = std::make_unique<MLMGBndryT<MF>>
            (this->m_grids[amrlev][0], this->m_dmap[amrlev][0], ncomp, this->m_geom[amrlev][0]);
        MF bc_data(this->m_grids[amrlev][0], this->m_dmap[amrlev][0], ncomp, 1);
        bc_data.setVal(0.0);

        m_bndry_cor[amrlev]->setBndryValues(*m_crse_cor_br[amrlev], 0, bc_data, 0, 0, ncomp,
                                            IntVect(this->m_amr_ref_ratio[amrlev-1]),
                                            InterpBndryDataT<MF>::IBD_max_order_DEF,
                                            m_interpbndry_halfwidth);

        Vector<Array<LinOpBCType,AMREX_SPACEDIM> > bclohi
            (ncomp,Array<LinOpBCType,AMREX_SPACEDIM>{{AMREX_D_DECL(BCType::Dirichlet,
                                                                   BCType::Dirichlet,
                                                                   BCType::Dirichlet)}});
        m_bndry_cor[amrlev]->setLOBndryConds(bclohi, bclohi, this->m_amr_ref_ratio[amrlev-1], RealVect{});
    }

    m_bcondloc.resize(this->m_num_amr_levels);
    for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev)
    {
        m_bcondloc[amrlev].resize(this->m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < this->m_num_mg_levels[amrlev]; ++mglev)
        {
            m_bcondloc[amrlev][mglev] = std::make_unique<BndryCondLoc>(this->m_grids[amrlev][mglev],
                                                                       this->m_dmap[amrlev][mglev],
                                                                       ncomp);
        }
    }
}

template <typename MF>
void
MLCellLinOpT<MF>::setLevelBC (int amrlev, const MF* a_levelbcdata, const MF* robinbc_a,
                              const MF* robinbc_b, const MF* robinbc_f)
{
    BL_PROFILE("MLCellLinOp::setLevelBC()");

    AMREX_ALWAYS_ASSERT(amrlev >= 0 && amrlev < this->m_num_amr_levels);

    const int ncomp = this->getNComp();

    MF zero;
    IntVect ng(1);
    if (this->hasHiddenDimension()) { ng[this->hiddenDirection()] = 0; }
    if (a_levelbcdata == nullptr) {
        zero.define(this->m_grids[amrlev][0], this->m_dmap[amrlev][0], ncomp, ng);
        zero.setVal(RT(0.0));
    } else {
        AMREX_ALWAYS_ASSERT(a_levelbcdata->nGrowVect().allGE(ng));
    }
    const MF& bcdata = (a_levelbcdata == nullptr) ? zero : *a_levelbcdata;

    IntVect br_ref_ratio(-1);

    if (amrlev == 0)
    {
        if (this->needsCoarseDataForBC())
        {
            // AMREX_ALWAYS_ASSERT(!this->hasHiddenDimension());
            if (this->hasHiddenDimension()) {
                int hidden_dir = this->hiddenDirection();
                AMREX_ALWAYS_ASSERT(this->m_coarse_data_crse_ratio[hidden_dir] == 1);
            }
            br_ref_ratio = this->m_coarse_data_crse_ratio.allGT(0) ? this->m_coarse_data_crse_ratio : IntVect(2);
            if (m_crse_sol_br[amrlev] == nullptr && br_ref_ratio.allGT(0))
            {
                const int in_rad = 0;
                const int out_rad = 1;
                const int extent_rad = 2;
                const IntVect crse_ratio = br_ref_ratio;
                BoxArray cba = this->m_grids[amrlev][0];
                cba.coarsen(crse_ratio);
                m_crse_sol_br[amrlev] = std::make_unique<BndryRegisterT<MF>>
                    (cba, this->m_dmap[amrlev][0], in_rad, out_rad, extent_rad, ncomp);
            }
            if (this->m_coarse_data_for_bc != nullptr) {
                AMREX_ALWAYS_ASSERT(this->m_coarse_data_crse_ratio.allGT(0));
                const Box& cbx = amrex::coarsen(this->m_geom[0][0].Domain(), this->m_coarse_data_crse_ratio);
                m_crse_sol_br[amrlev]->copyFrom(*this->m_coarse_data_for_bc, 0, 0, 0, ncomp,
                                                this->m_geom[0][0].periodicity(cbx));
            } else {
                m_crse_sol_br[amrlev]->setVal(RT(0.0));
            }
            m_bndry_sol[amrlev]->setBndryValues(*m_crse_sol_br[amrlev], 0,
                                                bcdata, 0, 0, ncomp, br_ref_ratio,
                                                InterpBndryDataT<MF>::IBD_max_order_DEF,
                                                m_interpbndry_halfwidth);
            br_ref_ratio = this->m_coarse_data_crse_ratio;
        }
        else
        {
            m_bndry_sol[amrlev]->setPhysBndryValues(bcdata,0,0,ncomp);
            br_ref_ratio = IntVect(1);
        }
    }
    else
    {
        m_bndry_sol[amrlev]->setPhysBndryValues(bcdata,0,0,ncomp);
        br_ref_ratio = IntVect(this->m_amr_ref_ratio[amrlev-1]);
    }

    auto crse_fine_bc_type = (amrlev == 0) ? this->m_coarse_fine_bc_type : LinOpBCType::Dirichlet;
    m_bndry_sol[amrlev]->setLOBndryConds(this->m_lobc, this->m_hibc, br_ref_ratio,
                                         this->m_coarse_bc_loc, crse_fine_bc_type);

    const Real* dx = this->m_geom[amrlev][0].CellSize();
    for (int mglev = 0; mglev < this->m_num_mg_levels[amrlev]; ++mglev)
    {
        m_bcondloc[amrlev][mglev]->setLOBndryConds(this->m_geom[amrlev][mglev], dx,
                                                   this->m_lobc, this->m_hibc,
                                                   br_ref_ratio, this->m_coarse_bc_loc,
                                                   this->m_domain_bloc_lo, this->m_domain_bloc_hi,
                                                   crse_fine_bc_type);
    }

    if (this->hasRobinBC()) {
        AMREX_ASSERT(robinbc_a != nullptr && robinbc_b != nullptr && robinbc_f != nullptr);
        m_robin_bcval[amrlev] = std::make_unique<MF>(this->m_grids[amrlev][0],
                                                     this->m_dmap[amrlev][0],
                                                     ncomp*3, 1);
        const Box& domain = this->m_geom[amrlev][0].Domain();
        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) { mfi_info.SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*m_robin_bcval[amrlev], mfi_info); mfi.isValid(); ++mfi) {
            Box const& vbx = mfi.validbox();
            Array4<RT const> const& ra = robinbc_a->const_array(mfi);
            Array4<RT const> const& rb = robinbc_b->const_array(mfi);
            Array4<RT const> const& rf = robinbc_f->const_array(mfi);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const Box& blo = amrex::adjCellLo(vbx, idim);
                const Box& bhi = amrex::adjCellHi(vbx, idim);
                bool outside_domain_lo = !(domain.contains(blo));
                bool outside_domain_hi = !(domain.contains(bhi));
                if ((!outside_domain_lo) && (!outside_domain_hi)) { continue; }
                for (int icomp = 0; icomp < ncomp; ++icomp) {
                    Array4<RT> const& rbc = (*m_robin_bcval[amrlev])[mfi].array(icomp*3);
                    if (this->m_lobc_orig[icomp][idim] == LinOpBCType::Robin && outside_domain_lo)
                    {
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(blo, i, j, k,
                        {
                            rbc(i,j,k,0) = ra(i,j,k,icomp);
                            rbc(i,j,k,1) = rb(i,j,k,icomp);
                            rbc(i,j,k,2) = rf(i,j,k,icomp);
                        });
                    }
                    if (this->m_hibc_orig[icomp][idim] == LinOpBCType::Robin && outside_domain_hi)
                    {
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bhi, i, j, k,
                        {
                            rbc(i,j,k,0) = ra(i,j,k,icomp);
                            rbc(i,j,k,1) = rb(i,j,k,icomp);
                            rbc(i,j,k,2) = rf(i,j,k,icomp);
                        });
                    }
                }
            }
        }
    }
}

template <typename MF>
void
MLCellLinOpT<MF>::update ()
{
    if (MLLinOpT<MF>::needsUpdate()) { MLLinOpT<MF>::update(); }
}

template <typename MF>
void
MLCellLinOpT<MF>::updateSolBC (int amrlev, const MF& crse_bcdata) const
{
    BL_PROFILE("MLCellLinOp::updateSolBC()");

    AMREX_ALWAYS_ASSERT(amrlev > 0);
    const int ncomp = this->getNComp();
    m_crse_sol_br[amrlev]->copyFrom(crse_bcdata, 0, 0, 0, ncomp,
                                    this->m_geom[amrlev-1][0].periodicity());
    m_bndry_sol[amrlev]->updateBndryValues(*m_crse_sol_br[amrlev], 0, 0, ncomp,
                                           IntVect(this->m_amr_ref_ratio[amrlev-1]),
                                           InterpBndryDataT<MF>::IBD_max_order_DEF,
                                           m_interpbndry_halfwidth);
}

template <typename MF>
void
MLCellLinOpT<MF>::updateCorBC (int amrlev, const MF& crse_bcdata) const
{
    BL_PROFILE("MLCellLinOp::updateCorBC()");
    AMREX_ALWAYS_ASSERT(amrlev > 0);
    const int ncomp = this->getNComp();
    m_crse_cor_br[amrlev]->copyFrom(crse_bcdata, 0, 0, 0, ncomp,
                                    this->m_geom[amrlev-1][0].periodicity());
    m_bndry_cor[amrlev]->updateBndryValues(*m_crse_cor_br[amrlev], 0, 0, ncomp,
                                           IntVect(this->m_amr_ref_ratio[amrlev-1]),
                                           InterpBndryDataT<MF>::IBD_max_order_DEF,
                                           m_interpbndry_halfwidth);
}

template <typename MF>
void
MLCellLinOpT<MF>::applyBC (int amrlev, int mglev, MF& in, BCMode bc_mode, StateMode,
                           const MLMGBndryT<MF>* bndry, bool skip_fillboundary) const
{
    BL_PROFILE("MLCellLinOp::applyBC()");
    // No coarsened boundary values, cannot apply inhomog at mglev>0.
    BL_ASSERT(mglev == 0 || bc_mode == BCMode::Homogeneous);
    BL_ASSERT(bndry != nullptr || bc_mode == BCMode::Homogeneous);

    const int ncomp = this->getNComp();
    const int cross = isCrossStencil();
    const int tensorop = isTensorOp();
    if (!skip_fillboundary) {
        in.FillBoundary(0, ncomp, this->m_geom[amrlev][mglev].periodicity(), cross);
    }

    int flagbc = bc_mode == BCMode::Inhomogeneous;
    const int imaxorder = this->maxorder;

    const Real* dxinv = this->m_geom[amrlev][mglev].InvCellSize();
    const RT dxi = static_cast<RT>(dxinv[0]);
    const RT dyi = (AMREX_SPACEDIM >= 2) ? static_cast<RT>(dxinv[1]) : RT(1.0);
    const RT dzi = (AMREX_SPACEDIM == 3) ? static_cast<RT>(dxinv[2]) : RT(1.0);

    const auto& maskvals = m_maskvals[amrlev][mglev];
    const auto& bcondloc = *m_bcondloc[amrlev][mglev];

    FAB foofab(Box::TheUnitBox(),ncomp);
    const auto& foo = foofab.const_array();

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) { mfi_info.SetDynamic(true); }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cross || tensorop || Gpu::notInLaunchRegion(),
                                     "non-cross stencil not support for gpu");

    const int hidden_direction = this->hiddenDirection();

#ifdef AMREX_USE_GPU
    if ((cross || tensorop) && Gpu::inLaunchRegion())
    {
        Vector<MLMGABCTag<RT>> tags;
        tags.reserve(in.local_size()*AMREX_SPACEDIM*ncomp);

        for (MFIter mfi(in); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            const auto& iofab = in.array(mfi);
            const auto & bdlv = bcondloc.bndryLocs(mfi);
            const auto & bdcv = bcondloc.bndryConds(mfi);

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (idim != hidden_direction) {
                    const Orientation olo(idim,Orientation::low);
                    const Orientation ohi(idim,Orientation::high);
                    const auto& bvlo = (bndry != nullptr) ?
                        bndry->bndryValues(olo).const_array(mfi) : foo;
                    const auto& bvhi = (bndry != nullptr) ?
                        bndry->bndryValues(ohi).const_array(mfi) : foo;
                    for (int icomp = 0; icomp < ncomp; ++icomp) {
                        tags.emplace_back(MLMGABCTag<RT>{iofab, bvlo, bvhi,
                                                 maskvals[olo].const_array(mfi),
                                                 maskvals[ohi].const_array(mfi),
                                                 bdlv[icomp][olo], bdlv[icomp][ohi],
                                                 amrex::adjCell(vbx,olo),
                                                 bdcv[icomp][olo], bdcv[icomp][ohi],
                                                 vbx.length(idim), icomp, idim});
                    }
                }
            }
        }

        ParallelFor(tags,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, MLMGABCTag<RT> const& tag) noexcept
        {
            if (tag.dir == 0)
            {
                mllinop_apply_bc_x(0, i, j, k, tag.blen, tag.fab,
                                   tag.mask_lo, tag.bctype_lo, tag.bcloc_lo, tag.bcval_lo,
                                   imaxorder, dxi, flagbc, tag.comp);
                mllinop_apply_bc_x(1, i+tag.blen+1, j, k, tag.blen, tag.fab,
                                   tag.mask_hi, tag.bctype_hi, tag.bcloc_hi, tag.bcval_hi,
                                   imaxorder, dxi, flagbc, tag.comp);
            }
#if (AMREX_SPACEDIM > 1)
            else
#if (AMREX_SPACEDIM > 2)
            if (tag.dir == 1)
#endif
            {
                mllinop_apply_bc_y(0, i, j, k, tag.blen, tag.fab,
                                   tag.mask_lo, tag.bctype_lo, tag.bcloc_lo, tag.bcval_lo,
                                   imaxorder, dyi, flagbc, tag.comp);
                mllinop_apply_bc_y(1, i, j+tag.blen+1, k, tag.blen, tag.fab,
                                   tag.mask_hi, tag.bctype_hi, tag.bcloc_hi, tag.bcval_hi,
                                   imaxorder, dyi, flagbc, tag.comp);
            }
#if (AMREX_SPACEDIM > 2)
            else {
                mllinop_apply_bc_z(0, i, j, k, tag.blen, tag.fab,
                                   tag.mask_lo, tag.bctype_lo, tag.bcloc_lo, tag.bcval_lo,
                                   imaxorder, dzi, flagbc, tag.comp);
                mllinop_apply_bc_z(1, i, j, k+tag.blen+1, tag.blen, tag.fab,
                                   tag.mask_hi, tag.bctype_hi, tag.bcloc_hi, tag.bcval_hi,
                                   imaxorder, dzi, flagbc, tag.comp);
            }
#endif
#endif
        });
    } else
#endif
    if (cross || tensorop)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(in, mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& vbx   = mfi.validbox();
            const auto& iofab = in.array(mfi);

            const auto & bdlv = bcondloc.bndryLocs(mfi);
            const auto & bdcv = bcondloc.bndryConds(mfi);

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                if (hidden_direction == idim) { continue; }
                const Orientation olo(idim,Orientation::low);
                const Orientation ohi(idim,Orientation::high);
                const Box blo = amrex::adjCellLo(vbx, idim);
                const Box bhi = amrex::adjCellHi(vbx, idim);
                const int blen = vbx.length(idim);
                const auto& mlo = maskvals[olo].array(mfi);
                const auto& mhi = maskvals[ohi].array(mfi);
                const auto& bvlo = (bndry != nullptr) ? bndry->bndryValues(olo).const_array(mfi) : foo;
                const auto& bvhi = (bndry != nullptr) ? bndry->bndryValues(ohi).const_array(mfi) : foo;
                for (int icomp = 0; icomp < ncomp; ++icomp) {
                    const BoundCond bctlo = bdcv[icomp][olo];
                    const BoundCond bcthi = bdcv[icomp][ohi];
                    const RT bcllo = bdlv[icomp][olo];
                    const RT bclhi = bdlv[icomp][ohi];
                    if (idim == 0) {
                        mllinop_apply_bc_x(0, blo, blen, iofab, mlo,
                                           bctlo, bcllo, bvlo,
                                           imaxorder, dxi, flagbc, icomp);
                        mllinop_apply_bc_x(1, bhi, blen, iofab, mhi,
                                           bcthi, bclhi, bvhi,
                                           imaxorder, dxi, flagbc, icomp);
                    } else if (idim == 1) {
                        mllinop_apply_bc_y(0, blo, blen, iofab, mlo,
                                           bctlo, bcllo, bvlo,
                                           imaxorder, dyi, flagbc, icomp);
                        mllinop_apply_bc_y(1, bhi, blen, iofab, mhi,
                                           bcthi, bclhi, bvhi,
                                           imaxorder, dyi, flagbc, icomp);
                    } else {
                        mllinop_apply_bc_z(0, blo, blen, iofab, mlo,
                                           bctlo, bcllo, bvlo,
                                           imaxorder, dzi, flagbc, icomp);
                        mllinop_apply_bc_z(1, bhi, blen, iofab, mhi,
                                           bcthi, bclhi, bvhi,
                                           imaxorder, dzi, flagbc, icomp);
                    }
                }
            }
        }
    }
    else
    {
#ifdef BL_NO_FORT
        amrex::Abort("amrex_mllinop_apply_bc not available when BL_NO_FORT=TRUE");
#else
        if constexpr (std::is_same_v<Real,RT>) {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            for (MFIter mfi(in, mfi_info); mfi.isValid(); ++mfi)
            {
                const Box& vbx   = mfi.validbox();

                const auto & bdlv = bcondloc.bndryLocs(mfi);
                const auto & bdcv = bcondloc.bndryConds(mfi);

                const RealTuple & bdl = bdlv[0];
                const BCTuple   & bdc = bdcv[0];

                for (OrientationIter oitr; oitr; ++oitr)
                {
                    const Orientation ori = oitr();

                    int cdr = ori;
                    RT  bcl = bdl[ori];
                    int bct = bdc[ori];

                    const auto& fsfab = (bndry != nullptr) ? bndry->bndryValues(ori)[mfi] : foofab;

                    const Mask& m = maskvals[ori][mfi];

                    amrex_mllinop_apply_bc(BL_TO_FORTRAN_BOX(vbx),
                                           BL_TO_FORTRAN_ANYD(in[mfi]),
                                           BL_TO_FORTRAN_ANYD(m),
                                           cdr, bct, bcl,
                                           BL_TO_FORTRAN_ANYD(fsfab),
                                           imaxorder, dxinv, flagbc, ncomp, cross);
                }
            }
        } else {
            amrex::Abort("Not supported");
        }
#endif
    }
}

template <typename MF>
BoxArray
MLCellLinOpT<MF>::makeNGrids (int grid_size) const
{
    const Box& dombx = this->m_geom[0].back().Domain();

    const BoxArray& old_ba = this->m_grids[0].back();
    const int N = old_ba.size();
    Vector<Box> bv;
    bv.reserve(N);
    for (int i = 0; i < N; ++i)
    {
        Box b = old_ba[i];
        b.coarsen(grid_size);
        b.refine(grid_size);
        IntVect sz = b.size();
        const IntVect nblks {AMREX_D_DECL(sz[0]/grid_size, sz[1]/grid_size, sz[2]/grid_size)};

        IntVect big = b.smallEnd() + grid_size - 1;
        b.setBig(big);

#if (AMREX_SPACEDIM == 3)
        for (int kk = 0; kk < nblks[2]; ++kk) {
#endif
#if (AMREX_SPACEDIM >= 2)
            for (int jj = 0; jj < nblks[1]; ++jj) {
#endif
                for (int ii = 0; ii < nblks[0]; ++ii)
                {
                    IntVect shft{AMREX_D_DECL(ii*grid_size,jj*grid_size,kk*grid_size)};
                    Box bb = amrex::shift(b,shft);
                    bb &= dombx;
                    bv.push_back(bb);
                }
#if (AMREX_SPACEDIM >= 2)
            }
#endif
#if (AMREX_SPACEDIM == 3)
        }
#endif
    }

    std::sort(bv.begin(), bv.end());
    bv.erase(std::unique(bv.begin(), bv.end()), bv.end());

    BoxList bl(std::move(bv));

    return BoxArray{std::move(bl)};
}

template <typename MF>
void
MLCellLinOpT<MF>::restriction (int amrlev, int cmglev, MF& crse, MF& fine) const
{
    const int ncomp = this->getNComp();
    IntVect ratio = (amrlev > 0) ? IntVect(2) : this->mg_coarsen_ratio_vec[cmglev-1];
    amrex::average_down(fine, crse, 0, ncomp, ratio);
}

template <typename MF>
void
MLCellLinOpT<MF>::interpolation (int amrlev, int fmglev, MF& fine, const MF& crse) const
{
    const int ncomp = this->getNComp();

    Dim3 ratio3 = {2,2,2};
    IntVect ratio = (amrlev > 0) ? IntVect(2) : this->mg_coarsen_ratio_vec[fmglev];
    AMREX_D_TERM(ratio3.x = ratio[0];,
                 ratio3.y = ratio[1];,
                 ratio3.z = ratio[2];);

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && fine.isFusingCandidate()) {
        auto const& finema = fine.arrays();
        auto const& crsema = crse.const_arrays();
        ParallelFor(fine, IntVect(0), ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            int ic = amrex::coarsen(i,ratio3.x);
            int jc = amrex::coarsen(j,ratio3.y);
            int kc = amrex::coarsen(k,ratio3.z);
            finema[box_no](i,j,k,n) += crsema[box_no](ic,jc,kc,n);
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(fine,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx    = mfi.tilebox();
            Array4<RT const> const& cfab = crse.const_array(mfi);
            Array4<RT> const& ffab = fine.array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                int ic = amrex::coarsen(i,ratio3.x);
                int jc = amrex::coarsen(j,ratio3.y);
                int kc = amrex::coarsen(k,ratio3.z);
                ffab(i,j,k,n) += cfab(ic,jc,kc,n);
            });
        }
    }
}

template <typename MF>
void
MLCellLinOpT<MF>::interpAssign (int amrlev, int fmglev, MF& fine, MF& crse) const
{
    const int ncomp = this->getNComp();

    const Geometry& crse_geom = this->Geom(amrlev,fmglev+1);
    const IntVect refratio = (amrlev > 0) ? IntVect(2) : this->mg_coarsen_ratio_vec[fmglev];
    const IntVect ng = crse.nGrowVect();

    MF cfine;
    const MF* cmf;

    if (amrex::isMFIterSafe(crse, fine))
    {
        crse.FillBoundary(crse_geom.periodicity());
        cmf = &crse;
    }
    else
    {
        BoxArray cba = fine.boxArray();
        cba.coarsen(refratio);
        cfine.define(cba, fine.DistributionMap(), ncomp, ng);
        cfine.setVal(RT(0.0));
        cfine.ParallelCopy(crse, 0, 0, ncomp, IntVect(0), ng, crse_geom.periodicity());
        cmf = & cfine;
    }

    bool isEB = fine.hasEBFabFactory();
    ignore_unused(isEB);

#ifdef AMREX_USE_EB
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(&(fine.Factory()));
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
#endif

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) { mfi_info.EnableTiling().SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& ff = fine.array(mfi);
        const auto& cc = cmf->array(mfi);
#ifdef AMREX_USE_EB
        bool call_lincc;
        if (isEB)
        {
            const auto& flag = (*flags)[mfi];
            if (flag.getType(amrex::grow(bx,1)) == FabType::regular) {
                call_lincc = true;
            } else {
                Array4<EBCellFlag const> const& flg = flag.const_array();
                AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE (bx, tbx,
                {
                    mlmg_eb_cc_interp_r<2>(tbx, ff, cc, flg, ncomp);
                });

                call_lincc = false;
            }
        }
        else
        {
            call_lincc = true;
        }
#else
        const bool call_lincc = true;
#endif
        if (call_lincc)
        {
#if (AMREX_SPACEDIM == 3)
            if (this->hasHiddenDimension()) {
                Box const& bx_2d = this->compactify(bx);
                auto const& ff_2d = this->compactify(ff);
                auto const& cc_2d = this->compactify(cc);
                AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE (bx_2d, tbx,
                {
                    TwoD::mlmg_lin_cc_interp_r2(tbx, ff_2d, cc_2d, ncomp);
                });
            } else
#endif
            {
                AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE (bx, tbx,
                {
                    mlmg_lin_cc_interp_r2(tbx, ff, cc, ncomp);
                });
            }
        }
    }
}

template <typename MF>
void
MLCellLinOpT<MF>::interpolationAmr (int famrlev, MF& fine, const MF& crse,
                                    IntVect const& /*nghost*/) const
{
    const int ncomp = this->getNComp();
    const int refratio = this->AMRRefRatio(famrlev-1);

#ifdef AMREX_USE_EB
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(this->Factory(famrlev));
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
#endif

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) { mfi_info.EnableTiling().SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine, mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& ff = fine.array(mfi);
        auto const& cc = crse.const_array(mfi);
#ifdef AMREX_USE_EB
        bool call_lincc;
        if (factory)
        {
            const auto& flag = (*flags)[mfi];
            if (flag.getType(amrex::grow(bx,1)) == FabType::regular) {
                call_lincc = true;
            } else {
                Array4<EBCellFlag const> const& flg = flag.const_array();
                switch(refratio) {
                case 2:
                {
                    AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE (bx, tbx,
                    {
                        mlmg_eb_cc_interp_r<2>(tbx, ff, cc, flg, ncomp);
                    });
                    break;
                }
                case 4:
                {
                    AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE (bx, tbx,
                    {
                        mlmg_eb_cc_interp_r<4>(tbx, ff, cc, flg, ncomp);
                    });
                    break;
                }
                default:
                    amrex::Abort("mlmg_eb_cc_interp: only refratio 2 and 4 are supported");
                }

                call_lincc = false;
            }
        }
        else
        {
            call_lincc = true;
        }
#else
        const bool call_lincc = true;
#endif
        if (call_lincc)
        {
            switch(refratio) {
            case 2:
            {
                AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE (bx, tbx,
                {
                    mlmg_lin_cc_interp_r2(tbx, ff, cc, ncomp);
                });
                break;
            }
            case 4:
            {
                AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE (bx, tbx,
                {
                    mlmg_lin_cc_interp_r4(tbx, ff, cc, ncomp);
                });
                break;
            }
            default:
                amrex::Abort("mlmg_lin_cc_interp: only refratio 2 and 4 are supported");
            }
        }
    }
}

template <typename MF>
void
MLCellLinOpT<MF>::averageDownSolutionRHS (int camrlev, MF& crse_sol, MF& crse_rhs,
                                          const MF& fine_sol, const MF& fine_rhs)
{
    const auto amrrr = this->AMRRefRatio(camrlev);
    const int ncomp = this->getNComp();
    amrex::average_down(fine_sol, crse_sol, 0, ncomp, amrrr);
    amrex::average_down(fine_rhs, crse_rhs, 0, ncomp, amrrr);
}

template <typename MF>
void
MLCellLinOpT<MF>::apply (int amrlev, int mglev, MF& out, MF& in, BCMode bc_mode,
                    StateMode s_mode, const MLMGBndryT<MF>* bndry) const
{
    BL_PROFILE("MLCellLinOp::apply()");
    applyBC(amrlev, mglev, in, bc_mode, s_mode, bndry);
    Fapply(amrlev, mglev, out, in);
}

template <typename MF>
void
MLCellLinOpT<MF>::smooth (int amrlev, int mglev, MF& sol, const MF& rhs,
                          bool skip_fillboundary) const
{
    BL_PROFILE("MLCellLinOp::smooth()");
    for (int redblack = 0; redblack < 2; ++redblack)
    {
        applyBC(amrlev, mglev, sol, BCMode::Homogeneous, StateMode::Solution,
                nullptr, skip_fillboundary);
        Fsmooth(amrlev, mglev, sol, rhs, redblack);
        skip_fillboundary = false;
    }
}

template <typename MF>
void
MLCellLinOpT<MF>::solutionResidual (int amrlev, MF& resid, MF& x, const MF& b,
                                    const MF* crse_bcdata)
{
    BL_PROFILE("MLCellLinOp::solutionResidual()");
    const int ncomp = this->getNComp();
    if (crse_bcdata != nullptr) {
        updateSolBC(amrlev, *crse_bcdata);
    }
    const int mglev = 0;
    apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, StateMode::Solution,
          m_bndry_sol[amrlev].get());

    AMREX_ASSERT(resid.nComp() == b.nComp());
    MF::Xpay(resid, RT(-1.0), b, 0, 0, ncomp, IntVect(0));
}

template <typename MF>
void
MLCellLinOpT<MF>::prepareForFluxes (int amrlev, const MF* crse_bcdata)
{
    if (crse_bcdata != nullptr) {
        updateSolBC(amrlev, *crse_bcdata);
    }
}

template <typename MF>
void
MLCellLinOpT<MF>::correctionResidual (int amrlev, int mglev, MF& resid, MF& x, const MF& b,
                                      BCMode bc_mode, const MF* crse_bcdata)
{
    BL_PROFILE("MLCellLinOp::correctionResidual()");
    const int ncomp = this->getNComp();
    if (bc_mode == BCMode::Inhomogeneous)
    {
        if (crse_bcdata)
        {
            AMREX_ASSERT(mglev == 0 && amrlev > 0);
            updateCorBC(amrlev, *crse_bcdata);
        }
        apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, StateMode::Correction,
              m_bndry_cor[amrlev].get());
    }
    else
    {
        AMREX_ASSERT(crse_bcdata == nullptr);
        apply(amrlev, mglev, resid, x, BCMode::Homogeneous, StateMode::Correction, nullptr);
    }

    MF::Xpay(resid, Real(-1.0), b, 0, 0, ncomp, IntVect(0));
}

template <typename MF>
void
MLCellLinOpT<MF>::reflux (int crse_amrlev, MF& res, const MF& crse_sol, const MF&,
                          MF&, MF& fine_sol, const MF&) const
{
    BL_PROFILE("MLCellLinOp::reflux()");

    auto& fluxreg = m_fluxreg[crse_amrlev];
    fluxreg.reset();

    const int ncomp = this->getNComp();

    const int fine_amrlev = crse_amrlev+1;

    Real dt = Real(1.0);
    const Real* crse_dx = this->m_geom[crse_amrlev][0].CellSize();
    const Real* fine_dx = this->m_geom[fine_amrlev][0].CellSize();

    const int mglev = 0;
    applyBC(fine_amrlev, mglev, fine_sol, BCMode::Inhomogeneous, StateMode::Solution,
            m_bndry_sol[fine_amrlev].get());

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) { mfi_info.EnableTiling().SetDynamic(true); }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        Array<FAB,AMREX_SPACEDIM> flux;
        Array<FAB*,AMREX_SPACEDIM> pflux {{ AMREX_D_DECL(flux.data(), flux.data()+1, flux.data()+2) }};
        Array<FAB const*,AMREX_SPACEDIM> cpflux {{ AMREX_D_DECL(flux.data(), flux.data()+1, flux.data()+2) }};

        for (MFIter mfi(crse_sol, mfi_info);  mfi.isValid(); ++mfi)
        {
            if (fluxreg.CrseHasWork(mfi))
            {
                const Box& tbx = mfi.tilebox();
                AMREX_D_TERM(flux[0].resize(amrex::surroundingNodes(tbx,0),ncomp);,
                             flux[1].resize(amrex::surroundingNodes(tbx,1),ncomp);,
                             flux[2].resize(amrex::surroundingNodes(tbx,2),ncomp););
                AMREX_D_TERM(Elixir elifx = flux[0].elixir();,
                             Elixir elify = flux[1].elixir();,
                             Elixir elifz = flux[2].elixir(););
                FFlux(crse_amrlev, mfi, pflux, crse_sol[mfi], Location::FaceCentroid);
                fluxreg.CrseAdd(mfi, cpflux, crse_dx, dt, RunOn::Gpu);
            }
        }

#ifdef AMREX_USE_OMP
#pragma omp barrier
#endif

        for (MFIter mfi(fine_sol, mfi_info);  mfi.isValid(); ++mfi)
        {
            if (fluxreg.FineHasWork(mfi))
            {
                const Box& tbx = mfi.tilebox();
                const int face_only = true;
                AMREX_D_TERM(flux[0].resize(amrex::surroundingNodes(tbx,0),ncomp);,
                             flux[1].resize(amrex::surroundingNodes(tbx,1),ncomp);,
                             flux[2].resize(amrex::surroundingNodes(tbx,2),ncomp););
                AMREX_D_TERM(Elixir elifx = flux[0].elixir();,
                             Elixir elify = flux[1].elixir();,
                             Elixir elifz = flux[2].elixir(););
                FFlux(fine_amrlev, mfi, pflux, fine_sol[mfi], Location::FaceCentroid, face_only);
                fluxreg.FineAdd(mfi, cpflux, fine_dx, dt, RunOn::Gpu);
            }
        }
    }

    fluxreg.Reflux(res);
}

template <typename MF>
void
MLCellLinOpT<MF>::compFlux (int amrlev, const Array<MF*,AMREX_SPACEDIM>& fluxes,
                            MF& sol, Location loc) const
{
    BL_PROFILE("MLCellLinOp::compFlux()");

    const int mglev = 0;
    const int ncomp = this->getNComp();
    applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution,
            m_bndry_sol[amrlev].get());

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) { mfi_info.EnableTiling().SetDynamic(true); }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        Array<FAB,AMREX_SPACEDIM> flux;
        Array<FAB*,AMREX_SPACEDIM> pflux {{ AMREX_D_DECL(flux.data(), flux.data()+1, flux.data()+2) }};
        for (MFIter mfi(sol, mfi_info);  mfi.isValid(); ++mfi)
        {
            const Box& tbx = mfi.tilebox();
            AMREX_D_TERM(flux[0].resize(amrex::surroundingNodes(tbx,0),ncomp);,
                         flux[1].resize(amrex::surroundingNodes(tbx,1),ncomp);,
                         flux[2].resize(amrex::surroundingNodes(tbx,2),ncomp););
            AMREX_D_TERM(Elixir elifx = flux[0].elixir();,
                         Elixir elify = flux[1].elixir();,
                         Elixir elifz = flux[2].elixir(););
            FFlux(amrlev, mfi, pflux, sol[mfi], loc);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const Box& nbx = mfi.nodaltilebox(idim);
                auto const& dst = fluxes[idim]->array(mfi);
                auto const& src =  pflux[idim]->const_array();
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D (nbx, ncomp, i, j, k, n,
                {
                    dst(i,j,k,n) = src(i,j,k,n);
                });
            }
        }
    }
}

template <typename MF>
void
MLCellLinOpT<MF>::compGrad (int amrlev, const Array<MF*,AMREX_SPACEDIM>& grad,
                            MF& sol, Location /*loc*/) const
{
    BL_PROFILE("MLCellLinOp::compGrad()");

    if (sol.nComp() > 1) {
        amrex::Abort("MLCellLinOp::compGrad called, but only works for single-component solves");
    }

    const int mglev = 0;
    applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, StateMode::Solution,
            m_bndry_sol[amrlev].get());

    const int ncomp = this->getNComp();

    AMREX_D_TERM(const RT dxi = static_cast<RT>(this->m_geom[amrlev][mglev].InvCellSize(0));,
                 const RT dyi = static_cast<RT>(this->m_geom[amrlev][mglev].InvCellSize(1));,
                 const RT dzi = static_cast<RT>(this->m_geom[amrlev][mglev].InvCellSize(2)););
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sol, TilingIfNotGPU());  mfi.isValid(); ++mfi)
    {
        AMREX_D_TERM(const Box& xbx = mfi.nodaltilebox(0);,
                     const Box& ybx = mfi.nodaltilebox(1);,
                     const Box& zbx = mfi.nodaltilebox(2););
        const auto& s = sol.array(mfi);
        AMREX_D_TERM(const auto& gx = grad[0]->array(mfi);,
                     const auto& gy = grad[1]->array(mfi);,
                     const auto& gz = grad[2]->array(mfi););

        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( xbx, ncomp, i, j, k, n,
        {
            gx(i,j,k,n) = dxi*(s(i,j,k,n) - s(i-1,j,k,n));
        });
#if (AMREX_SPACEDIM >= 2)
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( ybx, ncomp, i, j, k, n,
        {
            gy(i,j,k,n) = dyi*(s(i,j,k,n) - s(i,j-1,k,n));
        });
#endif
#if (AMREX_SPACEDIM == 3)
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( zbx, ncomp, i, j, k, n,
        {
            gz(i,j,k,n) = dzi*(s(i,j,k,n) - s(i,j,k-1,n));
        });
#endif
    }

    addInhomogNeumannFlux(amrlev, grad, sol, false);
}

template <typename MF>
void
MLCellLinOpT<MF>::applyMetricTerm (int amrlev, int mglev, MF& rhs) const
{
    amrex::ignore_unused(amrlev,mglev,rhs);
#if (AMREX_SPACEDIM != 3)
    if (!m_has_metric_term) { return; }

    const int ncomp = rhs.nComp();

    bool cc = rhs.ixType().cellCentered(0);

    const Geometry& geom = this->m_geom[amrlev][mglev];
    const RT dx = static_cast<RT>(geom.CellSize(0));
    const RT probxlo = static_cast<RT>(geom.ProbLo(0));

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        auto const& rhsarr = rhs.array(mfi);
#if (AMREX_SPACEDIM == 1)
        if (cc) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                RT rc = probxlo + (i+RT(0.5))*dx;
                rhsarr(i,j,k,n) *= rc*rc;
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                RT re = probxlo + i*dx;
                rhsarr(i,j,k,n) *= re*re;
            });
        }
#elif (AMREX_SPACEDIM == 2)
        if (cc) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                RT rc = probxlo + (i+RT(0.5))*dx;
                rhsarr(i,j,k,n) *= rc;
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                RT re = probxlo + i*dx;
                rhsarr(i,j,k,n) *= re;
            });
        }
#endif
    }
#endif
}

template <typename MF>
void
MLCellLinOpT<MF>::unapplyMetricTerm (int amrlev, int mglev, MF& rhs) const
{
    amrex::ignore_unused(amrlev,mglev,rhs);
#if (AMREX_SPACEDIM != 3)
    if (!m_has_metric_term) { return; }

    const int ncomp = rhs.nComp();

    bool cc = rhs.ixType().cellCentered(0);

    const Geometry& geom = this->m_geom[amrlev][mglev];
    const RT dx = static_cast<RT>(geom.CellSize(0));
    const RT probxlo = static_cast<RT>(geom.ProbLo(0));

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        auto const& rhsarr = rhs.array(mfi);
#if (AMREX_SPACEDIM == 1)
        if (cc) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                RT rcinv = RT(1.0)/(probxlo + (i+RT(0.5))*dx);
                rhsarr(i,j,k,n) *= rcinv*rcinv;
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                RT re = probxlo + i*dx;
                RT reinv = (re==RT(0.0)) ? RT(0.0) : RT(1.)/re;
                rhsarr(i,j,k,n) *= reinv*reinv;
            });
        }
#elif (AMREX_SPACEDIM == 2)
        if (cc) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                RT rcinv = RT(1.0)/(probxlo + (i+RT(0.5))*dx);
                rhsarr(i,j,k,n) *= rcinv;
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( tbx, ncomp, i, j, k, n,
            {
                RT re = probxlo + i*dx;
                RT reinv = (re==RT(0.0)) ? RT(0.0) : RT(1.)/re;
                rhsarr(i,j,k,n) *= reinv;
            });
        }
#endif
    }
#endif
}

template <typename MF>
auto
MLCellLinOpT<MF>::getSolvabilityOffset (int amrlev, int mglev, MF const& rhs) const
    -> Vector<RT>
{
    computeVolInv();

    const int ncomp = this->getNComp();
    Vector<RT> offset(ncomp);

#ifdef AMREX_USE_EB
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(this->Factory(amrlev,mglev));
    if (factory && !factory->isAllRegular())
    {
        if constexpr (std::is_same<MF,MultiFab>()) {
            const MultiFab& vfrac = factory->getVolFrac();
            for (int c = 0; c < ncomp; ++c) {
                offset[c] = amrex::Dot(rhs, c, vfrac, 0, 1, IntVect(0), true)
                    * m_volinv[amrlev][mglev];
            }
        } else {
            amrex::Abort("TODO: MLMG with EB only works with MultiFab");
        }
    }
    else
#endif
    {
        for (int c = 0; c < ncomp; ++c) {
            offset[c] = rhs.sum(c,IntVect(0),true) * m_volinv[amrlev][mglev];
        }
    }

    ParallelAllReduce::Sum(offset.data(), ncomp, ParallelContext::CommunicatorSub());

    return offset;
}

template <typename MF>
void
MLCellLinOpT<MF>::fixSolvabilityByOffset (int /*amrlev*/, int /*mglev*/, MF& rhs,
                                          Vector<RT> const& offset) const
{
    const int ncomp = this->getNComp();
    for (int c = 0; c < ncomp; ++c) {
        rhs.plus(-offset[c], c, 1);
    }
#ifdef AMREX_USE_EB
    if (!rhs.isAllRegular()) {
        if constexpr (std::is_same<MF,MultiFab>()) {
            amrex::EB_set_covered(rhs, 0, ncomp, 0, 0.0_rt);
        } else {
            amrex::Abort("amrex::EB_set_covered only works with MultiFab");
        }
    }
#endif
}

template <typename MF>
void
MLCellLinOpT<MF>::prepareForSolve ()
{
    BL_PROFILE("MLCellLinOp::prepareForSolve()");

    const int imaxorder = this->maxorder;
    const int ncomp = this->getNComp();
    const int hidden_direction = this->hiddenDirection();
    for (int amrlev = 0;  amrlev < this->m_num_amr_levels; ++amrlev)
    {
        for (int mglev = 0; mglev < this->m_num_mg_levels[amrlev]; ++mglev)
        {
            const auto& bcondloc = *m_bcondloc[amrlev][mglev];
            const auto& maskvals = m_maskvals[amrlev][mglev];

            const RT dxi = static_cast<RT>(this->m_geom[amrlev][mglev].InvCellSize(0));
            const RT dyi = static_cast<RT>((AMREX_SPACEDIM >= 2) ? this->m_geom[amrlev][mglev].InvCellSize(1) : Real(1.0));
            const RT dzi = static_cast<RT>((AMREX_SPACEDIM == 3) ? this->m_geom[amrlev][mglev].InvCellSize(2) : Real(1.0));

            auto& undrrelxr = this->m_undrrelxr[amrlev][mglev];
            MF foo(this->m_grids[amrlev][mglev], this->m_dmap[amrlev][mglev], ncomp, 0, MFInfo().SetAlloc(false));

#ifdef AMREX_USE_EB
            const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(this->m_factory[amrlev][mglev].get());
            const FabArray<EBCellFlagFab>* flags =
                (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
            auto area = (factory) ? factory->getAreaFrac()
                : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
            amrex::ignore_unused(area);
#endif

#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion()) {
#ifdef AMREX_USE_EB
                if (factory && !factory->isAllRegular()) {
#if defined(AMREX_USE_CUDA) && defined(_WIN32)
                    if           (!std::is_same<MF,MultiFab>()) {
#else
                    if constexpr (!std::is_same<MF,MultiFab>()) {
#endif
                        amrex::Abort("MLCellLinOp with EB only works with MultiFab");
                    } else {
                        Vector<MLMGPSEBTag<RT>> tags;
                        tags.reserve(foo.local_size()*AMREX_SPACEDIM*ncomp);

                        for (MFIter mfi(foo); mfi.isValid(); ++mfi)
                        {
                            const Box& vbx = mfi.validbox();

                            const auto & bdlv = bcondloc.bndryLocs(mfi);
                            const auto & bdcv = bcondloc.bndryConds(mfi);

                            auto fabtyp = (flags) ? (*flags)[mfi].getType(vbx) : FabType::regular;

                            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                            {
                                if (idim != hidden_direction && fabtyp != FabType::covered) {
                                    const Orientation olo(idim,Orientation::low);
                                    const Orientation ohi(idim,Orientation::high);
                                    auto const& ap = (fabtyp == FabType::singlevalued)
                                        ? area[idim]->const_array(mfi) : Array4<Real const>{};
                                    for (int icomp = 0; icomp < ncomp; ++icomp) {
                                        tags.emplace_back(MLMGPSEBTag<RT>{undrrelxr[olo].array(mfi),
                                                                          undrrelxr[ohi].array(mfi),
                                                                          ap,
                                                                          maskvals[olo].const_array(mfi),
                                                                          maskvals[ohi].const_array(mfi),
                                                                          bdlv[icomp][olo], bdlv[icomp][ohi],
                                                                          amrex::adjCell(vbx,olo),
                                                                          bdcv[icomp][olo], bdcv[icomp][ohi],
                                                                          vbx.length(idim), icomp, idim});
                                    }
                                }
                            }
                        }

                        ParallelFor(tags,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k, MLMGPSEBTag<RT> const& tag) noexcept
                        {
                            if (tag.ap) {
                                if (tag.dir == 0)
                                {
                                    mllinop_comp_interp_coef0_x_eb
                                        (0, i           , j, k, tag.blen, tag.flo, tag.mlo, tag.ap,
                                         tag.bctlo, tag.bcllo, imaxorder, dxi, tag.comp);
                                    mllinop_comp_interp_coef0_x_eb
                                        (1, i+tag.blen+1, j, k, tag.blen, tag.fhi, tag.mhi, tag.ap,
                                         tag.bcthi, tag.bclhi, imaxorder, dxi, tag.comp);
                                }
#if (AMREX_SPACEDIM > 1)
                                else
#if (AMREX_SPACEDIM > 2)
                                if (tag.dir == 1)
#endif
                                {
                                    mllinop_comp_interp_coef0_y_eb
                                        (0, i, j           , k, tag.blen, tag.flo, tag.mlo, tag.ap,
                                         tag.bctlo, tag.bcllo, imaxorder, dyi, tag.comp);
                                    mllinop_comp_interp_coef0_y_eb
                                        (1, i, j+tag.blen+1, k, tag.blen, tag.fhi, tag.mhi, tag.ap,
                                         tag.bcthi, tag.bclhi, imaxorder, dyi, tag.comp);
                                }
#if (AMREX_SPACEDIM > 2)
                                else {
                                    mllinop_comp_interp_coef0_z_eb
                                        (0, i, j, k           , tag.blen, tag.flo, tag.mlo, tag.ap,
                                         tag.bctlo, tag.bcllo, imaxorder, dzi, tag.comp);
                                    mllinop_comp_interp_coef0_z_eb
                                        (1, i, j, k+tag.blen+1, tag.blen, tag.fhi, tag.mhi, tag.ap,
                                         tag.bcthi, tag.bclhi, imaxorder, dzi, tag.comp);
                                }
#endif
#endif
                            } else {
                                if (tag.dir == 0)
                                {
                                    mllinop_comp_interp_coef0_x
                                        (0, i           , j, k, tag.blen, tag.flo, tag.mlo,
                                         tag.bctlo, tag.bcllo, imaxorder, dxi, tag.comp);
                                    mllinop_comp_interp_coef0_x
                                        (1, i+tag.blen+1, j, k, tag.blen, tag.fhi, tag.mhi,
                                         tag.bcthi, tag.bclhi, imaxorder, dxi, tag.comp);
                                }
#if (AMREX_SPACEDIM > 1)
                                else
#if (AMREX_SPACEDIM > 2)
                                if (tag.dir == 1)
#endif
                                {
                                    mllinop_comp_interp_coef0_y
                                        (0, i, j           , k, tag.blen, tag.flo, tag.mlo,
                                         tag.bctlo, tag.bcllo, imaxorder, dyi, tag.comp);
                                    mllinop_comp_interp_coef0_y
                                        (1, i, j+tag.blen+1, k, tag.blen, tag.fhi, tag.mhi,
                                         tag.bcthi, tag.bclhi, imaxorder, dyi, tag.comp);
                                }
#if (AMREX_SPACEDIM > 2)
                                else {
                                    mllinop_comp_interp_coef0_z
                                        (0, i, j, k           , tag.blen, tag.flo, tag.mlo,
                                         tag.bctlo, tag.bcllo, imaxorder, dzi, tag.comp);
                                    mllinop_comp_interp_coef0_z
                                        (1, i, j, k+tag.blen+1, tag.blen, tag.fhi, tag.mhi,
                                         tag.bcthi, tag.bclhi, imaxorder, dzi, tag.comp);
                                }
#endif
#endif
                            }
                        });
                    }
                } else
#endif
                {
                    Vector<MLMGPSTag<RT>> tags;
                    tags.reserve(foo.local_size()*AMREX_SPACEDIM*ncomp);

                    for (MFIter mfi(foo); mfi.isValid(); ++mfi)
                    {
                        const Box& vbx = mfi.validbox();

                        const auto & bdlv = bcondloc.bndryLocs(mfi);
                        const auto & bdcv = bcondloc.bndryConds(mfi);

                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                        {
                            if (idim != hidden_direction) {
                                const Orientation olo(idim,Orientation::low);
                                const Orientation ohi(idim,Orientation::high);
                                for (int icomp = 0; icomp < ncomp; ++icomp) {
                                    tags.emplace_back(MLMGPSTag<RT>{undrrelxr[olo].array(mfi),
                                                            undrrelxr[ohi].array(mfi),
                                                            maskvals[olo].const_array(mfi),
                                                            maskvals[ohi].const_array(mfi),
                                                            bdlv[icomp][olo], bdlv[icomp][ohi],
                                                            amrex::adjCell(vbx,olo),
                                                            bdcv[icomp][olo], bdcv[icomp][ohi],
                                                            vbx.length(idim), icomp, idim});
                                }
                            }
                        }
                    }

                    ParallelFor(tags,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, MLMGPSTag<RT> const& tag) noexcept
                    {
                        if (tag.dir == 0)
                        {
                            mllinop_comp_interp_coef0_x
                                (0, i           , j, k, tag.blen, tag.flo, tag.mlo,
                                 tag.bctlo, tag.bcllo, imaxorder, dxi, tag.comp);
                            mllinop_comp_interp_coef0_x
                                (1, i+tag.blen+1, j, k, tag.blen, tag.fhi, tag.mhi,
                                 tag.bcthi, tag.bclhi, imaxorder, dxi, tag.comp);
                        }
#if (AMREX_SPACEDIM > 1)
                        else
#if (AMREX_SPACEDIM > 2)
                        if (tag.dir == 1)
#endif
                        {
                            mllinop_comp_interp_coef0_y
                                (0, i, j           , k, tag.blen, tag.flo, tag.mlo,
                                 tag.bctlo, tag.bcllo, imaxorder, dyi, tag.comp);
                            mllinop_comp_interp_coef0_y
                                (1, i, j+tag.blen+1, k, tag.blen, tag.fhi, tag.mhi,
                                 tag.bcthi, tag.bclhi, imaxorder, dyi, tag.comp);
                        }
#if (AMREX_SPACEDIM > 2)
                        else {
                            mllinop_comp_interp_coef0_z
                                (0, i, j, k           , tag.blen, tag.flo, tag.mlo,
                                 tag.bctlo, tag.bcllo, imaxorder, dzi, tag.comp);
                            mllinop_comp_interp_coef0_z
                                (1, i, j, k+tag.blen+1, tag.blen, tag.fhi, tag.mhi,
                                 tag.bcthi, tag.bclhi, imaxorder, dzi, tag.comp);
                        }
#endif
#endif
                    });
                }
            } else
#endif
            {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
                for (MFIter mfi(foo, MFItInfo{}.SetDynamic(true)); mfi.isValid(); ++mfi)
                {
                    const Box& vbx = mfi.validbox();

                    const auto & bdlv = bcondloc.bndryLocs(mfi);
                    const auto & bdcv = bcondloc.bndryConds(mfi);

#ifdef AMREX_USE_EB
                    auto fabtyp = (flags) ? (*flags)[mfi].getType(vbx) : FabType::regular;
#endif
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                    {
                        if (idim == hidden_direction) { continue; }
                        const Orientation olo(idim,Orientation::low);
                        const Orientation ohi(idim,Orientation::high);
                        const Box blo = amrex::adjCellLo(vbx, idim);
                        const Box bhi = amrex::adjCellHi(vbx, idim);
                        const int blen = vbx.length(idim);
                        const auto& mlo = maskvals[olo].array(mfi);
                        const auto& mhi = maskvals[ohi].array(mfi);
                        const auto& flo = undrrelxr[olo].array(mfi);
                        const auto& fhi = undrrelxr[ohi].array(mfi);
                        for (int icomp = 0; icomp < ncomp; ++icomp) {
                            const BoundCond bctlo = bdcv[icomp][olo];
                            const BoundCond bcthi = bdcv[icomp][ohi];
                            const auto bcllo = bdlv[icomp][olo];
                            const auto bclhi = bdlv[icomp][ohi];
#ifdef AMREX_USE_EB
                            if (fabtyp == FabType::singlevalued) {
                                if constexpr (!std::is_same<MF,MultiFab>()) {
                                    amrex::Abort("MLCellLinOp with EB only works with MultiFab");
                                } else {
                                    auto const& ap = area[idim]->const_array(mfi);
                                    if (idim == 0) {
                                        mllinop_comp_interp_coef0_x_eb
                                            (0, blo, blen, flo, mlo, ap, bctlo, bcllo,
                                             imaxorder, dxi, icomp);
                                        mllinop_comp_interp_coef0_x_eb
                                            (1, bhi, blen, fhi, mhi, ap, bcthi, bclhi,
                                             imaxorder, dxi, icomp);
                                    } else if (idim == 1) {
                                        mllinop_comp_interp_coef0_y_eb
                                            (0, blo, blen, flo, mlo, ap, bctlo, bcllo,
                                             imaxorder, dyi, icomp);
                                        mllinop_comp_interp_coef0_y_eb
                                            (1, bhi, blen, fhi, mhi, ap, bcthi, bclhi,
                                             imaxorder, dyi, icomp);
                                    } else {
                                        mllinop_comp_interp_coef0_z_eb
                                            (0, blo, blen, flo, mlo, ap, bctlo, bcllo,
                                             imaxorder, dzi, icomp);
                                        mllinop_comp_interp_coef0_z_eb
                                            (1, bhi, blen, fhi, mhi, ap, bcthi, bclhi,
                                             imaxorder, dzi, icomp);
                                    }
                                }
                            } else if (fabtyp == FabType::regular)
#endif
                            {
                                if (idim == 0) {
                                    mllinop_comp_interp_coef0_x
                                        (0, blo, blen, flo, mlo, bctlo, bcllo,
                                         imaxorder, dxi, icomp);
                                    mllinop_comp_interp_coef0_x
                                        (1, bhi, blen, fhi, mhi, bcthi, bclhi,
                                         imaxorder, dxi, icomp);
                                } else if (idim == 1) {
                                    mllinop_comp_interp_coef0_y
                                        (0, blo, blen, flo, mlo, bctlo, bcllo,
                                         imaxorder, dyi, icomp);
                                    mllinop_comp_interp_coef0_y
                                        (1, bhi, blen, fhi, mhi, bcthi, bclhi,
                                         imaxorder, dyi, icomp);
                                } else {
                                    mllinop_comp_interp_coef0_z
                                        (0, blo, blen, flo, mlo, bctlo, bcllo,
                                         imaxorder, dzi, icomp);
                                    mllinop_comp_interp_coef0_z
                                        (1, bhi, blen, fhi, mhi, bcthi, bclhi,
                                         imaxorder, dzi, icomp);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template <typename MF>
auto
MLCellLinOpT<MF>::xdoty (int /*amrlev*/, int /*mglev*/, const MF& x, const MF& y, bool local) const
    -> RT
{
    const int ncomp = this->getNComp();
    const IntVect nghost(0);
    RT result = amrex::Dot(x,0,y,0,ncomp,nghost,true);
    if (!local) {
        ParallelAllReduce::Sum(result, ParallelContext::CommunicatorSub());
    }
    return result;
}

template <typename MF>
auto
MLCellLinOpT<MF>::dotProductPrecond (Vector<MF const*> const& x, Vector<MF const*> const& y) const -> RT
{
    const int ncomp = this->getNComp();
    const IntVect nghost(0);
    RT result = 0;
    for (int ilev = 0; ilev < this->NAMRLevels()-1; ++ilev) {
        result += amrex::Dot(*m_norm_fine_mask[ilev], *x[ilev], 0, *y[ilev], 0, ncomp, nghost, true);
    }
    result += amrex::Dot(*x[this->NAMRLevels()-1], 0,
                         *y[this->NAMRLevels()-1], 0, ncomp, nghost, true);
    ParallelAllReduce::Sum(result, ParallelContext::CommunicatorSub());
    return result;
}

template <typename MF>
auto
MLCellLinOpT<MF>::norm2Precond (Vector<MF const*> const& x) const -> RT
{
    const int ncomp = this->getNComp();
    const IntVect nghost(0);
    RT result = 0;
    for (int ilev = 0; ilev < this->NAMRLevels()-1; ++ilev) {
        result += amrex::Dot(*m_norm_fine_mask[ilev], *x[ilev], 0, ncomp, nghost, true);
    }
    result += amrex::Dot(*x[this->NAMRLevels()-1], 0, ncomp, nghost, true);
    ParallelAllReduce::Sum(result, ParallelContext::CommunicatorSub());
    return std::sqrt(result);
}

template <typename MF>
void
MLCellLinOpT<MF>::computeVolInv () const
{
    if (!m_volinv.empty()) { return; }

    m_volinv.resize(this->m_num_amr_levels);
    for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev) {
        m_volinv[amrlev].resize(this->NMGLevels(amrlev));
    }

    // We don't need to compute for every level

    auto f = [&] (int amrlev, int mglev) {
#ifdef AMREX_USE_EB
        const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(this->Factory(amrlev,mglev));
        if (factory && !factory->isAllRegular())
        {
            if constexpr (std::is_same<MF,MultiFab>()) {
                const auto& vfrac = factory->getVolFrac();
                m_volinv[amrlev][mglev] = vfrac.sum(0,true);
            } else {
                amrex::Abort("MLCellLinOp with EB only works with MultiFab");
            }
        }
        else
#endif
        {
            auto const npts = (this->m_coarse_fine_bc_type == LinOpBCType::Dirichlet)
                ? this->compactify(this->Geom(amrlev,mglev).Domain()).d_numPts()
                : this->m_grids[amrlev][mglev].d_numPts();
            AMREX_ASSERT(npts > 0.);
            m_volinv[amrlev][mglev] = RT(1.0 / npts);
        }
    };

    // amrlev = 0, mglev = 0
    f(0,0);

    int mgbottom = this->NMGLevels(0)-1;
    f(0,mgbottom);

#ifdef AMREX_USE_EB
    RT temp1, temp2;
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(this->Factory(0,0));
    if (factory && !factory->isAllRegular())
    {
        ParallelAllReduce::Sum<RT>({m_volinv[0][0], m_volinv[0][mgbottom]},
                                   ParallelContext::CommunicatorSub());
        temp1 = RT(1.0)/m_volinv[0][0];
        temp2 = RT(1.0)/m_volinv[0][mgbottom];
    }
    else
    {
        temp1 = m_volinv[0][0];
        temp2 = m_volinv[0][mgbottom];
    }
    m_volinv[0][0] = temp1;
    m_volinv[0][mgbottom] = temp2;
#endif
}

template <typename MF>
auto
MLCellLinOpT<MF>::normInf (int amrlev, MF const& mf, bool local) const -> RT
{
    const int ncomp = this->getNComp();
    const int finest_level = this->NAMRLevels() - 1;
    RT norm = RT(0.0);
#ifdef AMREX_USE_EB
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(this->Factory(amrlev));
    if (factory && !factory->isAllRegular()) {
#if defined(AMREX_USE_CUDA) && defined(_WIN32)
        if           (!std::is_same<MF,MultiFab>()) {
#else
        if constexpr (!std::is_same<MF,MultiFab>()) {
#endif
            amrex::Abort("MLCellLinOpT with EB only works with MultiFab");
        } else {
            const MultiFab& vfrac = factory->getVolFrac();
            if (amrlev == finest_level) {
#ifdef AMREX_USE_GPU
                if (Gpu::inLaunchRegion()) {
                    auto const& ma = mf.const_arrays();
                    auto const& vfrac_ma = vfrac.const_arrays();
                    norm = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{},
                                     mf, IntVect(0), ncomp,
                                     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n)
                                         -> GpuTuple<Real>
                                     {
                                         return std::abs(ma[box_no](i,j,k,n)
                                                                 *vfrac_ma[box_no](i,j,k));
                                     });
                } else
#endif
                {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:norm)
#endif
                    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
                        Box const& bx = mfi.tilebox();
                        auto const& fab = mf.const_array(mfi);
                        auto const& v = vfrac.const_array(mfi);
                        AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
                        {
                            norm = std::max(norm, std::abs(fab(i,j,k,n)*v(i,j,k)));
                        });
                    }
                }
            } else {
#ifdef AMREX_USE_GPU
                if (Gpu::inLaunchRegion()) {
                    auto const& ma = mf.const_arrays();
                    auto const& mask_ma = m_norm_fine_mask[amrlev]->const_arrays();
                    auto const& vfrac_ma = vfrac.const_arrays();
                    norm = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{},
                                     mf, IntVect(0), ncomp,
                                     [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n)
                                         -> GpuTuple<Real>
                                     {
                                         if (mask_ma[box_no](i,j,k)) {
                                             return std::abs(ma[box_no](i,j,k,n)
                                                                     *vfrac_ma[box_no](i,j,k));
                                         } else {
                                             return Real(0.0);
                                         }
                                     });
                } else
#endif
                {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:norm)
#endif
                    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
                        Box const& bx = mfi.tilebox();
                        auto const& fab = mf.const_array(mfi);
                        auto const& mask = m_norm_fine_mask[amrlev]->const_array(mfi);
                        auto const& v = vfrac.const_array(mfi);
                        AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
                        {
                            if (mask(i,j,k)) {
                                norm = std::max(norm, std::abs(fab(i,j,k,n)*v(i,j,k)));
                            }
                        });
                    }
                }
            }
        }
    } else
#endif
    {
        if (amrlev == finest_level) {
            norm = mf.norminf(0, ncomp, IntVect(0), true);
        } else {
            norm = mf.norminf(*m_norm_fine_mask[amrlev], 0, ncomp, IntVect(0), true);
        }
    }

    if (!local) { ParallelAllReduce::Max(norm, ParallelContext::CommunicatorSub()); }
    return norm;
}

template <typename MF>
void
MLCellLinOpT<MF>::averageDownAndSync (Vector<MF>& sol) const
{
    int ncomp = this->getNComp();
    for (int falev = this->NAMRLevels()-1; falev > 0; --falev)
    {
#ifdef AMREX_USE_EB
        if (!sol[falev].isAllRegular()) {
            if constexpr (std::is_same<MF,MultiFab>()) {
                amrex::EB_average_down(sol[falev], sol[falev-1], 0, ncomp, this->AMRRefRatio(falev-1));
            } else {
                amrex::Abort("EB_average_down only works with MultiFab");
            }
        } else
#endif
        {
            amrex::average_down(sol[falev], sol[falev-1], 0, ncomp, this->AMRRefRatio(falev-1));
        }
    }
}

template <typename MF>
void
MLCellLinOpT<MF>::avgDownResAmr (int clev, MF& cres, MF const& fres) const
{
#ifdef AMREX_USE_EB
    if (!fres.isAllRegular()) {
        if constexpr (std::is_same<MF,MultiFab>()) {
            amrex::EB_average_down(fres, cres, 0, this->getNComp(),
                                   this->AMRRefRatio(clev));
        } else {
            amrex::Abort("EB_average_down only works with MultiFab");
        }
    } else
#endif
    {
        amrex::average_down(fres, cres, 0, this->getNComp(),
                            this->AMRRefRatio(clev));
    }
}

template <typename MF>
void
MLCellLinOpT<MF>::beginPrecondBC ()
{
    this->m_precond_mode = true;

    if (m_bndry_sol_zero.empty()) {
        m_bndry_sol_zero.resize(m_bndry_sol.size());
        const int ncomp = this->getNComp();
        for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev) {
            m_bndry_sol_zero[amrlev] = std::make_unique<MLMGBndryT<MF>>
                (this->m_grids[amrlev][0],
                 this->m_dmap[amrlev][0],
                 ncomp,
                 this->m_geom[amrlev][0]);
        }
        std::swap(m_bndry_sol, m_bndry_sol_zero);
        MF const* coarse_data_for_bc_save = this->m_coarse_data_for_bc;
        this->m_coarse_data_for_bc = nullptr;
        for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev) {
            this->setLevelBC(amrlev, nullptr);
        }
        this->m_coarse_data_for_bc = coarse_data_for_bc_save;
    } else {
        std::swap(m_bndry_sol, m_bndry_sol_zero);
    }
}

template <typename MF>
void
MLCellLinOpT<MF>::endPrecondBC ()
{
    this->m_precond_mode = false;
    std::swap(m_bndry_sol, m_bndry_sol_zero);
}

extern template class MLCellLinOpT<MultiFab>;

using MLCellLinOp = MLCellLinOpT<MultiFab>;

}

#endif
