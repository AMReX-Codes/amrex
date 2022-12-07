#ifndef AMREX_MLALAPLACIAN_H_
#define AMREX_MLALAPLACIAN_H_
#include <AMReX_Config.H>

#include <AMReX_MLCellABecLap.H>
#include <AMReX_MLALap_K.H>
#include <AMReX_MultiFabUtil.H>

namespace amrex {

template <typename MF>
class MLALaplacianT
    : public MLCellABecLapT<MF>
{
public:

    using FAB = typename MF::fab_type;
    using RT  = typename MF::value_type;

    using BCType = LinOpBCType;
    using Location  = typename MLLinOpT<MF>::Location;

    MLALaplacianT () {}
    MLALaplacianT (const Vector<Geometry>& a_geom,
                   const Vector<BoxArray>& a_grids,
                   const Vector<DistributionMapping>& a_dmap,
                   const LPInfo& a_info = LPInfo(),
                   const Vector<FabFactory<FAB> const*>& a_factory = {},
                   const int a_ncomp = 1);
    virtual ~MLALaplacianT ();

    MLALaplacianT (const MLALaplacianT<MF>&) = delete;
    MLALaplacianT (MLALaplacianT<MF>&&) = delete;
    MLALaplacianT<MF>& operator= (const MLALaplacianT<MF>&) = delete;
    MLALaplacianT<MF>& operator= (MLALaplacianT<MF>&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FAB> const*>& a_factory = {});

    void setScalars (RT a, RT b) noexcept;
    void setACoeffs (int amrlev, const MF& alpha);

    virtual int getNComp () const override { return m_ncomp; }

    virtual bool needsUpdate () const override {
        return (m_needs_update || MLCellABecLapT<MF>::needsUpdate());
    }
    virtual void update () override;

    virtual void prepareForSolve () final override;
    virtual bool isSingular (int amrlev) const final override { return m_is_singular[amrlev]; }
    virtual bool isBottomSingular () const final override { return m_is_singular[0]; }
    virtual void Fapply (int amrlev, int mglev, MF& out, const MF& in) const final override;
    virtual void Fsmooth (int amrlev, int mglev, MF& sol, const MF& rsh, int redblack) const final override;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const Array<FAB*,AMREX_SPACEDIM>& flux,
                        const FAB& sol, Location /* loc */,
                        const int face_only=0) const final override;

    virtual void normalize (int marlve, int mglev, MF& mf) const final override;

    virtual RT getAScalar () const final override { return m_a_scalar; }
    virtual RT getBScalar () const final override { return m_b_scalar; }
    virtual MF const* getACoeffs (int amrlev, int mglev) const final override
        { return &(m_a_coeffs[amrlev][mglev]); }
    virtual Array<MF const*,AMREX_SPACEDIM> getBCoeffs (int /*amrlev*/, int /*mglev*/) const final override
        { return {{ AMREX_D_DECL(nullptr,nullptr,nullptr)}}; }

    virtual std::unique_ptr<MLLinOpT<MF>> makeNLinOp (int /*grid_size*/) const final override {
        amrex::Abort("MLALaplacian::makeNLinOp: Not implemented");
        return std::unique_ptr<MLLinOpT<MF>>{};
    }

    void averageDownCoeffsSameAmrLevel (int amrlev, Vector<MF>& a);
    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);

private:

    bool m_needs_update = true;

    RT m_a_scalar = std::numeric_limits<RT>::quiet_NaN();
    RT m_b_scalar = std::numeric_limits<RT>::quiet_NaN();
    Vector<Vector<MF> > m_a_coeffs;

    Vector<int> m_is_singular;

    int m_ncomp = 1;

    void updateSingularFlag ();
};

template <typename MF>
MLALaplacianT<MF>::MLALaplacianT (const Vector<Geometry>& a_geom,
                                  const Vector<BoxArray>& a_grids,
                                  const Vector<DistributionMapping>& a_dmap,
                                  const LPInfo& a_info,
                                  const Vector<FabFactory<FAB> const*>& a_factory,
                                  const int a_ncomp)
    : m_ncomp(a_ncomp)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

template <typename MF>
void
MLALaplacianT<MF>::define (const Vector<Geometry>& a_geom,
                           const Vector<BoxArray>& a_grids,
                           const Vector<DistributionMapping>& a_dmap,
                           const LPInfo& a_info,
                           const Vector<FabFactory<FAB> const*>& a_factory)
{
    BL_PROFILE("MLALaplacian::define()");

    MLCellABecLapT<MF>::define(a_geom, a_grids, a_dmap, a_info, a_factory);

    const int ncomp = this->getNComp();

    m_a_coeffs.resize(this->m_num_amr_levels);
    for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev)
    {
        m_a_coeffs[amrlev].resize(this->m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < this->m_num_mg_levels[amrlev]; ++mglev)
        {
            m_a_coeffs[amrlev][mglev].define(this->m_grids[amrlev][mglev],
                                             this->m_dmap[amrlev][mglev], ncomp, 0);
        }
    }
}

template <typename MF>
MLALaplacianT<MF>::~MLALaplacianT ()
{}

template <typename MF>
void
MLALaplacianT<MF>::setScalars (RT a, RT b) noexcept
{
    m_a_scalar = a;
    m_b_scalar = b;
    if (a == RT(0.0))
    {
        for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev)
        {
            m_a_coeffs[amrlev][0].setVal(RT(0.0));
        }
    }
}

template <typename MF>
void
MLALaplacianT<MF>::setACoeffs (int amrlev, const MF& alpha)
{
    const int ncomp = this->getNComp();
    m_a_coeffs[amrlev][0].LocalCopy(alpha, 0, 0, ncomp, IntVect(0));
    m_needs_update = true;
}

template <typename MF>
void
MLALaplacianT<MF>::averageDownCoeffs ()
{
    BL_PROFILE("MLALaplacian::averageDownCoeffs()");

    for (int amrlev = this->m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
        auto& fine_a_coeffs = m_a_coeffs[amrlev];

        averageDownCoeffsSameAmrLevel(amrlev, fine_a_coeffs);
        averageDownCoeffsToCoarseAmrLevel(amrlev);
    }

    averageDownCoeffsSameAmrLevel(0, m_a_coeffs[0]);
}

template <typename MF>
void
MLALaplacianT<MF>::averageDownCoeffsSameAmrLevel (int amrlev, Vector<MF>& a)
{
    const int ncomp = this->getNComp();
    const int nmglevs = a.size();
    for (int mglev = 1; mglev < nmglevs; ++mglev)
    {
        if (m_a_scalar == RT(0.0))
        {
            a[mglev].setVal(RT(0.0));
        }
        else
        {
            AMREX_ASSERT(amrlev == 0 || !this->hasHiddenDimension());
            IntVect ratio = (amrlev > 0) ? IntVect(this->mg_coarsen_ratio) : this->mg_coarsen_ratio_vec[mglev-1];
            amrex::average_down(a[mglev-1], a[mglev], 0, ncomp, ratio);
        }
    }
}

template <typename MF>
void
MLALaplacianT<MF>::averageDownCoeffsToCoarseAmrLevel (int flev)
{
    const int ncomp = this->getNComp();
    auto& fine_a_coeffs = m_a_coeffs[flev  ].back();
    auto& crse_a_coeffs = m_a_coeffs[flev-1].front();

    if (m_a_scalar != RT(0.0)) {
        // We coarsen from the back of flev to the front of flev-1.
        // So we use this->mg_coarsen_ratio.
        amrex::average_down(fine_a_coeffs, crse_a_coeffs, 0, ncomp, this->mg_coarsen_ratio);
    }
}

template <typename MF>
void
MLALaplacianT<MF>::updateSingularFlag ()
{
    m_is_singular.clear();
    m_is_singular.resize(this->m_num_amr_levels, false);
    auto itlo = std::find(this->m_lobc[0].begin(), this->m_lobc[0].end(), BCType::Dirichlet);
    auto ithi = std::find(this->m_hibc[0].begin(), this->m_hibc[0].end(), BCType::Dirichlet);
    if (itlo == this->m_lobc[0].end() && ithi == this->m_hibc[0].end())
    {  // No Dirichlet
        for (int alev = 0; alev < this->m_num_amr_levels; ++alev)
        {
            if (this->m_domain_covered[alev])
            {
                if (m_a_scalar == RT(0.0))
                {
                    m_is_singular[alev] = true;
                }
                else
                {
                    // We are only testing component 0 here, assuming the others
                    // are similar.
                    RT asum = m_a_coeffs[alev].back().sum(0,IntVect(0));
                    RT amax = m_a_coeffs[alev].back().norminf(0,1,IntVect(0));
                    m_is_singular[alev] = (asum <= amax * RT(1.e-12));
                }
            }
        }
    }
}

template <typename MF>
void
MLALaplacianT<MF>::prepareForSolve ()
{
    BL_PROFILE("MLALaplacian::prepareForSolve()");
    MLCellABecLapT<MF>::prepareForSolve();
    averageDownCoeffs();
    updateSingularFlag();
    m_needs_update = false;
}

template <typename MF>
void
MLALaplacianT<MF>::update ()
{
    if (MLCellABecLapT<MF>::needsUpdate()) MLCellABecLapT<MF>::update();
    averageDownCoeffs();
    updateSingularFlag();
    m_needs_update = false;
}

template <typename MF>
void
MLALaplacianT<MF>::Fapply (int amrlev, int mglev, MF& out, const MF& in) const
{
    BL_PROFILE("MLALaplacian::Fapply()");

    const int ncomp = this->getNComp();

    const MF& acoef = m_a_coeffs[amrlev][mglev];

    const GpuArray<RT,AMREX_SPACEDIM> dxinv
        {AMREX_D_DECL(RT(this->m_geom[amrlev][mglev].InvCellSize(0)),
                      RT(this->m_geom[amrlev][mglev].InvCellSize(1)),
                      RT(this->m_geom[amrlev][mglev].InvCellSize(2)))};
#if (AMREX_SPACEDIM < 3)
    const RT dx = RT(this->m_geom[amrlev][mglev].CellSize(0));
    const RT probxlo = RT(this->m_geom[amrlev][mglev].ProbLo(0));
#endif

#if (AMREX_SPACEDIM == 3)
    GpuArray<RT,2> dhinv {this->get_d0(dxinv[0], dxinv[1], dxinv[2]),
                          this->get_d1(dxinv[0], dxinv[1], dxinv[2])};
#endif

    const RT ascalar = m_a_scalar;
    const RT bscalar = m_b_scalar;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& xfab = in.array(mfi);
        const auto& yfab = out.array(mfi);
        const auto& afab = acoef.array(mfi);

#if (AMREX_SPACEDIM != 3)
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlalap_adotx_m(tbx, yfab, xfab, afab, dxinv, ascalar, bscalar, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlalap_adotx(tbx, yfab, xfab, afab, dxinv, ascalar, bscalar, ncomp);
            });
        }
#else
        if (this->hasHiddenDimension()) {
            Box const& bx2d = this->compactify(bx);
            const auto& xfab2d = this->compactify(xfab);
            const auto& yfab2d = this->compactify(yfab);
            const auto& afab2d = this->compactify(afab);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx2d, tbx2d,
            {
                TwoD::mlalap_adotx(tbx2d, yfab2d, xfab2d, afab2d, dhinv, ascalar, bscalar, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlalap_adotx(tbx, yfab, xfab, afab, dxinv, ascalar, bscalar, ncomp);
            });
        }
#endif
    }
}

template <typename MF>
void
MLALaplacianT<MF>::normalize (int amrlev, int mglev, MF& mf) const
{
    BL_PROFILE("MLALaplacian::normalize()");

    const int ncomp = this->getNComp();

    const MF& acoef = m_a_coeffs[amrlev][mglev];

    const GpuArray<RT,AMREX_SPACEDIM> dxinv
        {AMREX_D_DECL(RT(this->m_geom[amrlev][mglev].InvCellSize(0)),
                      RT(this->m_geom[amrlev][mglev].InvCellSize(1)),
                      RT(this->m_geom[amrlev][mglev].InvCellSize(2)))};
#if (AMREX_SPACEDIM < 3)
    const RT dx = RT(this->m_geom[amrlev][mglev].CellSize(0));
    const RT probxlo = RT(this->m_geom[amrlev][mglev].ProbLo(0));
#endif

#if (AMREX_SPACEDIM == 3)
    GpuArray<RT,2> dhinv {this->get_d0(dxinv[0], dxinv[1], dxinv[2]),
                          this->get_d1(dxinv[0], dxinv[1], dxinv[2])};
#endif

    const RT ascalar = m_a_scalar;
    const RT bscalar = m_b_scalar;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& fab = mf.array(mfi);
        const auto& afab = acoef.array(mfi);

#if (AMREX_SPACEDIM != 3)
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlalap_normalize_m(tbx, fab, afab, dxinv, ascalar, bscalar, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlalap_normalize(tbx, fab, afab, dxinv, ascalar, bscalar, ncomp);
            });
        }
#else
        if (this->hasHiddenDimension()) {
            Box const& bx2d = this->compactify(bx);
            const auto&  fab2d = this->compactify(fab);
            const auto& afab2d = this->compactify(afab);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx2d, tbx2d,
            {
                TwoD::mlalap_normalize(tbx2d, fab2d, afab2d, dhinv, ascalar, bscalar, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlalap_normalize(tbx, fab, afab, dxinv, ascalar, bscalar, ncomp);
            });
        }
#endif
    }
}

template <typename MF>
void
MLALaplacianT<MF>::Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs, int redblack) const
{
    BL_PROFILE("MLALaplacian::Fsmooth()");

    const int ncomp = this->getNComp();

    const MF& acoef = m_a_coeffs[amrlev][mglev];
    const auto& undrrelxr = this->m_undrrelxr[amrlev][mglev];
    const auto& maskvals  = this->m_maskvals [amrlev][mglev];

    OrientationIter oitr;

    const auto& f0 = undrrelxr[oitr()]; ++oitr;
    const auto& f1 = undrrelxr[oitr()]; ++oitr;
#if (AMREX_SPACEDIM > 1)
    const auto& f2 = undrrelxr[oitr()]; ++oitr;
    const auto& f3 = undrrelxr[oitr()]; ++oitr;
#if (AMREX_SPACEDIM > 2)
    const auto& f4 = undrrelxr[oitr()]; ++oitr;
    const auto& f5 = undrrelxr[oitr()]; ++oitr;
#endif
#endif

    const MultiMask& mm0 = maskvals[0];
    const MultiMask& mm1 = maskvals[1];
#if (AMREX_SPACEDIM > 1)
    const MultiMask& mm2 = maskvals[2];
    const MultiMask& mm3 = maskvals[3];
#if (AMREX_SPACEDIM > 2)
    const MultiMask& mm4 = maskvals[4];
    const MultiMask& mm5 = maskvals[5];
#endif
#endif

    const Real* dxinv = this->m_geom[amrlev][mglev].InvCellSize();
    AMREX_D_TERM(const RT dhx = m_b_scalar*RT(dxinv[0]*dxinv[0]);,
                 const RT dhy = m_b_scalar*RT(dxinv[1]*dxinv[1]);,
                 const RT dhz = m_b_scalar*RT(dxinv[2]*dxinv[2]););

#if (AMREX_SPACEDIM == 3)
    RT dh0 = this->get_d0(dhx, dhy, dhz);
    RT dh1 = this->get_d1(dhx, dhy, dhz);
#endif

#if (AMREX_SPACEDIM < 3)
    const RT dx = RT(this->m_geom[amrlev][mglev].CellSize(0));
    const RT probxlo = RT(this->m_geom[amrlev][mglev].ProbLo(0));
#endif

    const RT alpha = m_a_scalar;

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sol,mfi_info); mfi.isValid(); ++mfi)
    {
        const auto& m0 = mm0.array(mfi);
        const auto& m1 = mm1.array(mfi);
#if (AMREX_SPACEDIM > 1)
        const auto& m2 = mm2.array(mfi);
        const auto& m3 = mm3.array(mfi);
#if (AMREX_SPACEDIM > 2)
        const auto& m4 = mm4.array(mfi);
        const auto& m5 = mm5.array(mfi);
#endif
#endif

        const Box& tbx = mfi.tilebox();
        const Box& vbx = mfi.validbox();
        const auto& solnfab = sol.array(mfi);
        const auto& rhsfab  = rhs.array(mfi);
        const auto& afab    = acoef.array(mfi);

        const auto& f0fab = f0.array(mfi);
        const auto& f1fab = f1.array(mfi);
#if (AMREX_SPACEDIM > 1)
        const auto& f2fab = f2.array(mfi);
        const auto& f3fab = f3.array(mfi);
#if (AMREX_SPACEDIM > 2)
        const auto& f4fab = f4.array(mfi);
        const auto& f5fab = f5.array(mfi);
#endif
#endif

#if (AMREX_SPACEDIM == 1)
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
            {
                mlalap_gsrb_m(thread_box, solnfab, rhsfab, alpha, dhx,
                              afab,
                              f0fab, m0,
                              f1fab, m1,
                              vbx, redblack,
                              dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
            {
                mlalap_gsrb(thread_box, solnfab, rhsfab, alpha, dhx,
                            afab,
                            f0fab, m0,
                            f1fab, m1,
                            vbx, redblack, ncomp);
            });
        }

#endif

#if (AMREX_SPACEDIM == 2)
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
            {
                mlalap_gsrb_m(thread_box, solnfab, rhsfab, alpha, dhx, dhy,
                              afab,
                              f0fab, m0,
                              f1fab, m1,
                              f2fab, m2,
                              f3fab, m3,
                              vbx, redblack,
                              dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
            {
                mlalap_gsrb(thread_box, solnfab, rhsfab, alpha, dhx, dhy,
                            afab,
                            f0fab, m0,
                            f1fab, m1,
                            f2fab, m2,
                            f3fab, m3,
                            vbx, redblack, ncomp);
            });
        }
#endif

#if (AMREX_SPACEDIM == 3)
        if (this->hasHiddenDimension()) {
            Box const& tbx_2d = this->compactify(tbx);
            Box const& vbx_2d = this->compactify(vbx);
            const auto& solnfab_2d = this->compactify(solnfab);
            const auto& rhsfab_2d = this->compactify(rhsfab);
            const auto& afab_2d = this->compactify(afab);
            const auto& f0fab_2d = this->compactify(this->get_d0(f0fab,f1fab,f2fab));
            const auto& f1fab_2d = this->compactify(this->get_d1(f0fab,f1fab,f2fab));
            const auto& f2fab_2d = this->compactify(this->get_d0(f3fab,f4fab,f5fab));
            const auto& f3fab_2d = this->compactify(this->get_d1(f3fab,f4fab,f5fab));
            const auto& m0_2d = this->compactify(this->get_d0(m0,m1,m2));
            const auto& m1_2d = this->compactify(this->get_d1(m0,m1,m2));
            const auto& m2_2d = this->compactify(this->get_d0(m3,m4,m5));
            const auto& m3_2d = this->compactify(this->get_d1(m3,m4,m5));
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx_2d, thread_box,
            {
                TwoD::mlalap_gsrb(thread_box, solnfab_2d, rhsfab_2d, alpha, dh0, dh1,
                                  afab_2d,
                                  f0fab_2d, m0_2d,
                                  f1fab_2d, m1_2d,
                                  f2fab_2d, m2_2d,
                                  f3fab_2d, m3_2d,
                                  vbx_2d, redblack, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
            {
                mlalap_gsrb(thread_box, solnfab, rhsfab, alpha, dhx, dhy, dhz,
                            afab,
                            f0fab, m0,
                            f1fab, m1,
                            f2fab, m2,
                            f3fab, m3,
                            f4fab, m4,
                            f5fab, m5,
                            vbx, redblack, ncomp);
            });
        }
#endif
    }
}

template <typename MF>
void
MLALaplacianT<MF>::FFlux (int amrlev, const MFIter& mfi,
                          const Array<FAB*,AMREX_SPACEDIM>& flux,
                          const FAB& sol, Location, const int face_only) const
{
    AMREX_ASSERT(!this->hasHiddenDimension());

    BL_PROFILE("MLALaplacian::FFlux()");

    const int ncomp = this->getNComp();
    const int mglev = 0;
    const Box& box = mfi.tilebox();
    const Real* dxinv = this->m_geom[amrlev][mglev].InvCellSize();

    AMREX_D_TERM(const auto& fxarr = flux[0]->array();,
                 const auto& fyarr = flux[1]->array();,
                 const auto& fzarr = flux[2]->array(););
    const auto& solarr = sol.array();

#if (AMREX_SPACEDIM != 3)
    const RT dx = RT(this->m_geom[amrlev][mglev].CellSize(0));
    const RT probxlo = RT(this->m_geom[amrlev][mglev].ProbLo(0));
#endif

#if (AMREX_SPACEDIM == 3)
    if (face_only) {
        RT fac = m_b_scalar * RT(dxinv[0]);
        Box blo = amrex::bdryLo(box, 0);
        int blen = box.length(0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlalap_flux_xface(tbox, fxarr, solarr, fac, blen, ncomp);
        });
        fac = m_b_scalar * RT(dxinv[1]);
        blo = amrex::bdryLo(box, 1);
        blen = box.length(1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlalap_flux_yface(tbox, fyarr, solarr, fac, blen, ncomp);
        });
        fac = m_b_scalar * RT(dxinv[2]);
        blo = amrex::bdryLo(box, 2);
        blen = box.length(2);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlalap_flux_zface(tbox, fzarr, solarr, fac, blen, ncomp);
        });
    } else {
        RT fac = m_b_scalar * RT(dxinv[0]);
        Box bflux = amrex::surroundingNodes(box, 0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlalap_flux_x(tbox, fxarr, solarr, fac, ncomp);
        });
        fac = m_b_scalar * RT(dxinv[1]);
        bflux = amrex::surroundingNodes(box, 1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlalap_flux_y(tbox, fyarr, solarr, fac, ncomp);
        });
        fac = m_b_scalar * RT(dxinv[2]);
        bflux = amrex::surroundingNodes(box, 2);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlalap_flux_z(tbox, fzarr, solarr, fac, ncomp);
        });
    }
#elif (AMREX_SPACEDIM == 2)
    if (face_only) {
        RT fac = m_b_scalar * RT(dxinv[0]);
        Box blo = amrex::bdryLo(box, 0);
        int blen = box.length(0);
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlalap_flux_xface_m(tbox, fxarr, solarr, fac, blen, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlalap_flux_xface(tbox, fxarr, solarr, fac, blen, ncomp);
            });
        }
        fac = m_b_scalar * RT(dxinv[1]);
        blo = amrex::bdryLo(box, 1);
        blen = box.length(1);
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlalap_flux_yface_m(tbox, fyarr, solarr, fac, blen, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlalap_flux_yface(tbox, fyarr, solarr, fac, blen, ncomp);
            });
        }
    } else {
        RT fac = m_b_scalar * RT(dxinv[0]);
        Box bflux = amrex::surroundingNodes(box, 0);
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlalap_flux_x_m(tbox, fxarr, solarr, fac, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlalap_flux_x(tbox, fxarr, solarr, fac, ncomp);
            });
        }
        fac = m_b_scalar * RT(dxinv[1]);
        bflux = amrex::surroundingNodes(box, 1);
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlalap_flux_y_m(tbox, fyarr, solarr, fac, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlalap_flux_y(tbox, fyarr, solarr, fac, ncomp);
            });
        }
    }
#else
    if (face_only) {
        RT fac = m_b_scalar * RT(dxinv[0]);
        Box blo = amrex::bdryLo(box, 0);
        int blen = box.length(0);
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlalap_flux_xface_m(tbox, fxarr, solarr, fac, blen, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlalap_flux_xface(tbox, fxarr, solarr, fac, blen, ncomp);
            });
        }
    } else {
        RT fac = m_b_scalar * RT(dxinv[0]);
        Box bflux = amrex::surroundingNodes(box, 0);
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlalap_flux_x_m(tbox, fxarr, solarr, fac, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlalap_flux_x(tbox, fxarr, solarr, fac, ncomp);
            });
        }
    }
#endif
}

extern template class MLALaplacianT<MultiFab>;

using MLALaplacian = MLALaplacianT<MultiFab>;

}

#endif
