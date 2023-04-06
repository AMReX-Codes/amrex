#ifndef AMREX_MLPOISSON_H_
#define AMREX_MLPOISSON_H_
#include <AMReX_Config.H>

#include <AMReX_MLCellABecLap.H>
#include <AMReX_MLPoisson_K.H>
#include <AMReX_MLALaplacian.H>

namespace amrex {

// del dot grad phi

template <typename MF>
class MLPoissonT
    : public MLCellABecLapT<MF>
{
public:

    using FAB = typename MF::fab_type;
    using RT  = typename MF::value_type;

    using BCType = LinOpBCType;
    using Location  = typename MLLinOpT<MF>::Location;

    MLPoissonT () = default;
    MLPoissonT (const Vector<Geometry>& a_geom,
                const Vector<BoxArray>& a_grids,
                const Vector<DistributionMapping>& a_dmap,
                const LPInfo& a_info = LPInfo(),
                const Vector<FabFactory<FAB> const*>& a_factory = {});
    MLPoissonT (const Vector<Geometry>& a_geom,
                const Vector<BoxArray>& a_grids,
                const Vector<DistributionMapping>& a_dmap,
                const Vector<iMultiFab const*>& a_overset_mask, // 1: unknown, 0: known
                const LPInfo& a_info = LPInfo(),
                const Vector<FabFactory<FAB> const*>& a_factory = {});
    ~MLPoissonT () override;

    MLPoissonT (const MLPoissonT<MF>&) = delete;
    MLPoissonT (MLPoissonT<MF>&&) = delete;
    MLPoissonT<MF>& operator= (const MLPoissonT<MF>&) = delete;
    MLPoissonT<MF>& operator= (MLPoissonT<MF>&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FAB> const*>& a_factory = {});

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const Vector<iMultiFab const*>& a_overset_mask,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FAB> const*>& a_factory = {});

    void prepareForSolve () final;
    [[nodiscard]] bool isSingular (int amrlev) const final { return m_is_singular[amrlev]; }
    [[nodiscard]] bool isBottomSingular () const final { return m_is_singular[0]; }
    void Fapply (int amrlev, int mglev, MF& out, const MF& in) const final;
    void Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs, int redblack) const final;
    void FFlux (int amrlev, const MFIter& mfi,
                        const Array<FAB*,AMREX_SPACEDIM>& flux,
                        const FAB& sol, Location loc, int face_only=0) const final;

    void normalize (int amrlev, int mglev, MF& mf) const final;

    [[nodiscard]] RT getAScalar () const final { return RT(0.0); }
    [[nodiscard]] RT getBScalar () const final { return RT(-1.0); }
    [[nodiscard]] MF const* getACoeffs (int /*amrlev*/, int /*mglev*/) const final { return nullptr; }
    [[nodiscard]] Array<MF const*,AMREX_SPACEDIM> getBCoeffs (int /*amrlev*/, int /*mglev*/) const final
        { return {{ AMREX_D_DECL(nullptr,nullptr,nullptr)}}; }

    [[nodiscard]] std::unique_ptr<MLLinOpT<MF>> makeNLinOp (int grid_size) const final;

    [[nodiscard]] bool supportNSolve () const final;

    void copyNSolveSolution (MF& dst, MF const& src) const final;

    //! Compute dphi/dn on domain faces after the solver has converged.
    void get_dpdn_on_domain_faces (Array<MF*,AMREX_SPACEDIM> const& dpdn,
                                   MF const& phi);

private:

    Vector<int> m_is_singular;
};

template <typename MF>
MLPoissonT<MF>::MLPoissonT (const Vector<Geometry>& a_geom,
                            const Vector<BoxArray>& a_grids,
                            const Vector<DistributionMapping>& a_dmap,
                            const LPInfo& a_info,
                            const Vector<FabFactory<FAB> const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

template <typename MF>
MLPoissonT<MF>::MLPoissonT (const Vector<Geometry>& a_geom,
                            const Vector<BoxArray>& a_grids,
                            const Vector<DistributionMapping>& a_dmap,
                            const Vector<iMultiFab const*>& a_overset_mask,
                            const LPInfo& a_info,
                            const Vector<FabFactory<FAB> const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_overset_mask, a_info, a_factory);
}

template <typename MF>
void
MLPoissonT<MF>::define (const Vector<Geometry>& a_geom,
                        const Vector<BoxArray>& a_grids,
                        const Vector<DistributionMapping>& a_dmap,
                        const LPInfo& a_info,
                        const Vector<FabFactory<FAB> const*>& a_factory)
{
    BL_PROFILE("MLPoisson::define()");
    MLCellABecLapT<MF>::define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

template <typename MF>
void
MLPoissonT<MF>::define (const Vector<Geometry>& a_geom,
                        const Vector<BoxArray>& a_grids,
                        const Vector<DistributionMapping>& a_dmap,
                        const Vector<iMultiFab const*>& a_overset_mask,
                        const LPInfo& a_info,
                        const Vector<FabFactory<FAB> const*>& a_factory)
{
    BL_PROFILE("MLPoisson::define(overset)");
    MLCellABecLapT<MF>::define(a_geom, a_grids, a_dmap, a_overset_mask, a_info, a_factory);
}

template <typename MF>
MLPoissonT<MF>::~MLPoissonT () = default;

template <typename MF>
void
MLPoissonT<MF>::prepareForSolve ()
{
    BL_PROFILE("MLPoisson::prepareForSolve()");

    MLCellABecLapT<MF>::prepareForSolve();

    m_is_singular.clear();
    m_is_singular.resize(this->m_num_amr_levels, false);
    auto itlo = std::find(this->m_lobc[0].begin(), this->m_lobc[0].end(), BCType::Dirichlet);
    auto ithi = std::find(this->m_hibc[0].begin(), this->m_hibc[0].end(), BCType::Dirichlet);
    if (itlo == this->m_lobc[0].end() && ithi == this->m_hibc[0].end())
    {  // No Dirichlet
        for (int alev = 0; alev < this->m_num_amr_levels; ++alev)
        {
            // For now this assumes that overset regions are treated as Dirichlet bc's
            if (this->m_domain_covered[alev] && !this->m_overset_mask[alev][0])
            {
                m_is_singular[alev] = true;
            }
        }
    }
}

template <typename MF>
void
MLPoissonT<MF>::Fapply (int amrlev, int mglev, MF& out, const MF& in) const
{
    BL_PROFILE("MLPoisson::Fapply()");

    const Real* dxinv = this->m_geom[amrlev][mglev].InvCellSize();

    AMREX_D_TERM(const RT dhx = RT(dxinv[0]*dxinv[0]);,
                 const RT dhy = RT(dxinv[1]*dxinv[1]);,
                 const RT dhz = RT(dxinv[2]*dxinv[2]););

#if (AMREX_SPACEDIM == 3)
    RT dh0 = this->get_d0(dhx, dhy, dhz);
    RT dh1 = this->get_d1(dhx, dhy, dhz);
#endif

#if (AMREX_SPACEDIM < 3)
    const RT dx = RT(this->m_geom[amrlev][mglev].CellSize(0));
    const RT probxlo = RT(this->m_geom[amrlev][mglev].ProbLo(0));
#endif

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && out.isFusingCandidate() && !this->hasHiddenDimension()) {
        auto const& xma = in.const_arrays();
        auto const& yma = out.arrays();
        if (this->m_overset_mask[amrlev][mglev]) {
            AMREX_ASSERT(!this->m_has_metric_term);
            const auto& osmma = this->m_overset_mask[amrlev][mglev]->const_arrays();
            ParallelFor(out,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                amrex::ignore_unused(j,k);
                mlpoisson_adotx_os(AMREX_D_DECL(i,j,k), yma[box_no], xma[box_no], osmma[box_no],
                                   AMREX_D_DECL(dhx,dhy,dhz));
            });
        } else {
#if (AMREX_SPACEDIM < 3)
            if (this->m_has_metric_term) {
                ParallelFor(out,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                {
                    amrex::ignore_unused(j,k);
                    mlpoisson_adotx_m(AMREX_D_DECL(i,j,k), yma[box_no], xma[box_no],
                                      AMREX_D_DECL(dhx,dhy,dhz), dx, probxlo);
                });
            } else
#endif
            {
                ParallelFor(out,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                {
                    amrex::ignore_unused(j,k);
                    mlpoisson_adotx(AMREX_D_DECL(i,j,k), yma[box_no], xma[box_no],
                                    AMREX_D_DECL(dhx,dhy,dhz));
                });
            }
        }
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto& xfab = in.array(mfi);
            const auto& yfab = out.array(mfi);

            if (this->m_overset_mask[amrlev][mglev]) {
                AMREX_ASSERT(!this->m_has_metric_term);
                const auto& osm = this->m_overset_mask[amrlev][mglev]->const_array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D( bx, i, j, k,
                {
                    amrex::ignore_unused(j,k);
                    mlpoisson_adotx_os(AMREX_D_DECL(i,j,k), yfab, xfab, osm,
                                       AMREX_D_DECL(dhx,dhy,dhz));
                });
            } else {
#if (AMREX_SPACEDIM == 3)
                if (this->hasHiddenDimension()) {
                    Box const& bx2d = this->compactify(bx);
                    const auto& xfab2d = this->compactify(xfab);
                    const auto& yfab2d = this->compactify(yfab);
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx2d, i, j, k,
                    {
                        amrex::ignore_unused(k);
                        TwoD::mlpoisson_adotx(i, j, yfab2d, xfab2d, dh0, dh1);
                    });
                } else {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                    {
                        mlpoisson_adotx(i, j, k, yfab, xfab, dhx, dhy, dhz);
                    });
                }
#elif (AMREX_SPACEDIM == 2)
                if (this->m_has_metric_term) {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                    {
                        amrex::ignore_unused(k);
                        mlpoisson_adotx_m(i, j, yfab, xfab, dhx, dhy, dx, probxlo);
                    });
                } else {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                    {
                        amrex::ignore_unused(k);
                        mlpoisson_adotx(i, j, yfab, xfab, dhx, dhy);
                    });
                }
#elif (AMREX_SPACEDIM == 1)
                if (this->m_has_metric_term) {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                    {
                        amrex::ignore_unused(j,k);
                        mlpoisson_adotx_m(i, yfab, xfab, dhx, dx, probxlo);
                    });
                } else {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                    {
                        amrex::ignore_unused(j,k);
                        mlpoisson_adotx(i, yfab, xfab, dhx);
                    });
                }
#endif
            }
        }
    }
}

template <typename MF>
void
MLPoissonT<MF>::normalize (int amrlev, int mglev, MF& mf) const
{
    amrex::ignore_unused(amrlev,mglev,mf);
#if (AMREX_SPACEDIM != 3)
    BL_PROFILE("MLPoisson::normalize()");

    if (!this->m_has_metric_term) return;

    const Real* dxinv = this->m_geom[amrlev][mglev].InvCellSize();
    AMREX_D_TERM(const RT dhx = RT(dxinv[0]*dxinv[0]);,
                 const RT dhy = RT(dxinv[1]*dxinv[1]);,
                 const RT dhz = RT(dxinv[2]*dxinv[2]););
    const RT dx = RT(this->m_geom[amrlev][mglev].CellSize(0));
    const RT probxlo = RT(this->m_geom[amrlev][mglev].ProbLo(0));

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && mf.isFusingCandidate()) {
        auto const& ma = mf.arrays();
        ParallelFor(mf,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            mlpoisson_normalize(i,j,k, ma[box_no], AMREX_D_DECL(dhx,dhy,dhz), dx, probxlo);
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto& fab = mf.array(mfi);

#if (AMREX_SPACEDIM == 2)
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlpoisson_normalize(i,j,k, fab, dhx, dhy, dx, probxlo);
            });
#else
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlpoisson_normalize(i,j,k, fab, dhx, dx, probxlo);
            });
#endif
        }
    }
#endif
}

template <typename MF>
void
MLPoissonT<MF>::Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs, int redblack) const
{
    BL_PROFILE("MLPoisson::Fsmooth()");

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
    AMREX_D_TERM(const RT dhx = RT(dxinv[0]*dxinv[0]);,
                 const RT dhy = RT(dxinv[1]*dxinv[1]);,
                 const RT dhz = RT(dxinv[2]*dxinv[2]););

#if (AMREX_SPACEDIM == 3)
    RT dh0 = RT(this->get_d0(dhx, dhy, dhz));
    RT dh1 = RT(this->get_d1(dhx, dhy, dhz));
#endif

#if (AMREX_SPACEDIM < 3)
    const RT dx = RT(this->m_geom[amrlev][mglev].CellSize(0));
    const RT probxlo = RT(this->m_geom[amrlev][mglev].ProbLo(0));
#endif

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && sol.isFusingCandidate()
        && ! this->hasHiddenDimension())
    {
        const auto& m0ma = mm0.const_arrays();
        const auto& m1ma = mm1.const_arrays();
#if (AMREX_SPACEDIM > 1)
        const auto& m2ma = mm2.const_arrays();
        const auto& m3ma = mm3.const_arrays();
#if (AMREX_SPACEDIM > 2)
        const auto& m4ma = mm4.const_arrays();
        const auto& m5ma = mm5.const_arrays();
#endif
#endif

        const auto& solnma = sol.arrays();
        const auto& rhsma  = rhs.const_arrays();

        AMREX_ALWAYS_ASSERT(rhs.nGrowVect() == 0);

        const auto& f0ma = f0.const_arrays();
        const auto& f1ma = f1.const_arrays();
#if (AMREX_SPACEDIM > 1)
        const auto& f2ma = f2.const_arrays();
        const auto& f3ma = f3.const_arrays();
#if (AMREX_SPACEDIM > 2)
        const auto& f4ma = f4.const_arrays();
        const auto& f5ma = f5.const_arrays();
#endif
#endif

#if (AMREX_SPACEDIM == 1)
        if (this->m_overset_mask[amrlev][mglev]) {
            AMREX_ASSERT(!this->m_has_metric_term);
            const auto& osmma = this->m_overset_mask[amrlev][mglev]->const_arrays();
            ParallelFor(sol,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                Box vbx(rhsma[box_no]);
                mlpoisson_gsrb_os(i, j, k, solnma[box_no], rhsma[box_no],
                                  osmma[box_no], dhx,
                                  f0ma[box_no], m0ma[box_no],
                                  f1ma[box_no], m1ma[box_no],
                                  vbx, redblack);
            });
        } else if (this->m_has_metric_term) {
            ParallelFor(sol,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                Box vbx(rhsma[box_no]);
                mlpoisson_gsrb_m(i, j, k, solnma[box_no], rhsma[box_no], dhx,
                                 f0ma[box_no], m0ma[box_no],
                                 f1ma[box_no], m1ma[box_no],
                                 vbx, redblack,
                                 dx, probxlo);
            });
        } else {
            ParallelFor(sol,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                Box vbx(rhsma[box_no]);
                mlpoisson_gsrb(i, j, k, solnma[box_no], rhsma[box_no], dhx,
                               f0ma[box_no], m0ma[box_no],
                               f1ma[box_no], m1ma[box_no],
                               vbx, redblack);
            });
        }
#endif

#if (AMREX_SPACEDIM == 2)
        if (this->m_overset_mask[amrlev][mglev]) {
            AMREX_ASSERT(!this->m_has_metric_term);
            const auto& osmma = this->m_overset_mask[amrlev][mglev]->const_arrays();
            ParallelFor(sol,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                Box vbx(rhsma[box_no]);
                mlpoisson_gsrb_os(i, j, k, solnma[box_no], rhsma[box_no],
                                  osmma[box_no], dhx, dhy,
                                  f0ma[box_no], m0ma[box_no],
                                  f1ma[box_no], m1ma[box_no],
                                  f2ma[box_no], m2ma[box_no],
                                  f3ma[box_no], m3ma[box_no],
                                  vbx, redblack);
            });
        } else if (this->m_has_metric_term) {
            ParallelFor(sol,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                Box vbx(rhsma[box_no]);
                mlpoisson_gsrb_m(i, j, k, solnma[box_no], rhsma[box_no], dhx, dhy,
                                 f0ma[box_no], m0ma[box_no],
                                 f1ma[box_no], m1ma[box_no],
                                 f2ma[box_no], m2ma[box_no],
                                 f3ma[box_no], m3ma[box_no],
                                 vbx, redblack,
                                 dx, probxlo);
            });
        } else {
            ParallelFor(sol,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                Box vbx(rhsma[box_no]);
                mlpoisson_gsrb(i, j, k, solnma[box_no], rhsma[box_no], dhx, dhy,
                               f0ma[box_no], m0ma[box_no],
                               f1ma[box_no], m1ma[box_no],
                               f2ma[box_no], m2ma[box_no],
                               f3ma[box_no], m3ma[box_no],
                               vbx, redblack);
            });
        }
#endif

#if (AMREX_SPACEDIM == 3)
        if (this->m_overset_mask[amrlev][mglev]) {
            AMREX_ASSERT(!this->m_has_metric_term);
            const auto& osmma = this->m_overset_mask[amrlev][mglev]->const_arrays();
            ParallelFor(sol,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                Box vbx(rhsma[box_no]);
                mlpoisson_gsrb_os(i, j, k, solnma[box_no], rhsma[box_no],
                                  osmma[box_no], dhx, dhy, dhz,
                                  f0ma[box_no], m0ma[box_no],
                                  f1ma[box_no], m1ma[box_no],
                                  f2ma[box_no], m2ma[box_no],
                                  f3ma[box_no], m3ma[box_no],
                                  f4ma[box_no], m4ma[box_no],
                                  f5ma[box_no], m5ma[box_no],
                                  vbx, redblack);
            });
        } else {
            ParallelFor(sol,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                Box vbx(rhsma[box_no]);
                mlpoisson_gsrb(i, j, k, solnma[box_no], rhsma[box_no], dhx, dhy, dhz,
                               f0ma[box_no], m0ma[box_no],
                               f1ma[box_no], m1ma[box_no],
                               f2ma[box_no], m2ma[box_no],
                               f3ma[box_no], m3ma[box_no],
                               f4ma[box_no], m4ma[box_no],
                               f5ma[box_no], m5ma[box_no],
                               vbx, redblack);
            });
        }
#endif
    } else
#endif
    {
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
            if (this->m_overset_mask[amrlev][mglev]) {
                AMREX_ASSERT(!this->m_has_metric_term);
                const auto& osm = this->m_overset_mask[amrlev][mglev]->const_array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( tbx, i, j, k,
                {
                    mlpoisson_gsrb_os(i, j, k, solnfab, rhsfab, osm, dhx,
                                      f0fab, m0,
                                      f1fab, m1,
                                      vbx, redblack);
                });
            } else if (this->m_has_metric_term) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( tbx, i, j, k,
                {
                    mlpoisson_gsrb_m(i, j, k, solnfab, rhsfab, dhx,
                                     f0fab, m0,
                                     f1fab, m1,
                                     vbx, redblack,
                                     dx, probxlo);
                });
            } else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( tbx, i, j, k,
                {
                    mlpoisson_gsrb(i, j, k, solnfab, rhsfab, dhx,
                                   f0fab, m0,
                                   f1fab, m1,
                                   vbx, redblack);
                });
            }
#endif

#if (AMREX_SPACEDIM == 2)
            if (this->m_overset_mask[amrlev][mglev]) {
                AMREX_ASSERT(!this->m_has_metric_term);
                const auto& osm = this->m_overset_mask[amrlev][mglev]->const_array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( tbx, i, j, k,
                {
                    mlpoisson_gsrb_os(i, j, k, solnfab, rhsfab, osm, dhx, dhy,
                                      f0fab, m0,
                                      f1fab, m1,
                                      f2fab, m2,
                                      f3fab, m3,
                                      vbx, redblack);
                });
            } else if (this->m_has_metric_term) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( tbx, i, j, k,
                {
                    mlpoisson_gsrb_m(i, j, k, solnfab, rhsfab, dhx, dhy,
                                     f0fab, m0,
                                     f1fab, m1,
                                     f2fab, m2,
                                     f3fab, m3,
                                     vbx, redblack,
                                     dx, probxlo);
                });
            } else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( tbx, i, j, k,
                {
                    mlpoisson_gsrb(i, j, k, solnfab, rhsfab, dhx, dhy,
                                   f0fab, m0,
                                   f1fab, m1,
                                   f2fab, m2,
                                   f3fab, m3,
                                   vbx, redblack);
                });
            }

#endif

#if (AMREX_SPACEDIM == 3)
            if (this->m_overset_mask[amrlev][mglev]) {
                AMREX_ASSERT(!this->m_has_metric_term);
                const auto& osm = this->m_overset_mask[amrlev][mglev]->const_array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( tbx, i, j, k,
                {
                    mlpoisson_gsrb_os(i, j, k, solnfab, rhsfab, osm, dhx, dhy, dhz,
                                      f0fab, m0,
                                      f1fab, m1,
                                      f2fab, m2,
                                      f3fab, m3,
                                      f4fab, m4,
                                      f5fab, m5,
                                      vbx, redblack);
                });
            } else if (this->hasHiddenDimension()) {
                Box const& tbx_2d = this->compactify(tbx);
                Box const& vbx_2d = this->compactify(vbx);
                const auto& solnfab_2d = this->compactify(solnfab);
                const auto& rhsfab_2d = this->compactify(rhsfab);
                const auto& f0fab_2d = this->compactify(this->get_d0(f0fab,f1fab,f2fab));
                const auto& f1fab_2d = this->compactify(this->get_d1(f0fab,f1fab,f2fab));
                const auto& f2fab_2d = this->compactify(this->get_d0(f3fab,f4fab,f5fab));
                const auto& f3fab_2d = this->compactify(this->get_d1(f3fab,f4fab,f5fab));
                const auto& m0_2d = this->compactify(this->get_d0(m0,m1,m2));
                const auto& m1_2d = this->compactify(this->get_d1(m0,m1,m2));
                const auto& m2_2d = this->compactify(this->get_d0(m3,m4,m5));
                const auto& m3_2d = this->compactify(this->get_d1(m3,m4,m5));
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( tbx_2d, i, j, k,
                {
                    TwoD::mlpoisson_gsrb(i, j, k, solnfab_2d, rhsfab_2d, dh0, dh1,
                                         f0fab_2d, m0_2d,
                                         f1fab_2d, m1_2d,
                                         f2fab_2d, m2_2d,
                                         f3fab_2d, m3_2d,
                                         vbx_2d, redblack);
                });
            } else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( tbx, i, j, k,
                {
                    mlpoisson_gsrb(i, j, k, solnfab, rhsfab, dhx, dhy, dhz,
                                   f0fab, m0,
                                   f1fab, m1,
                                   f2fab, m2,
                                   f3fab, m3,
                                   f4fab, m4,
                                   f5fab, m5,
                                   vbx, redblack);
                });
            }
#endif
        }
    }
}

template <typename MF>
void
MLPoissonT<MF>::FFlux (int amrlev, const MFIter& mfi,
                       const Array<FAB*,AMREX_SPACEDIM>& flux,
                       const FAB& sol, Location, const int face_only) const
{
    AMREX_ASSERT(!this->hasHiddenDimension());

    BL_PROFILE("MLPoisson::FFlux()");

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
        RT fac = RT(dxinv[0]);
        Box blo = amrex::bdryLo(box, 0);
        int blen = box.length(0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlpoisson_flux_xface(tbox, fxarr, solarr, fac, blen);
        });
        fac = RT(dxinv[1]);
        blo = amrex::bdryLo(box, 1);
        blen = box.length(1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlpoisson_flux_yface(tbox, fyarr, solarr, fac, blen);
        });
        fac = RT(dxinv[2]);
        blo = amrex::bdryLo(box, 2);
        blen = box.length(2);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlpoisson_flux_zface(tbox, fzarr, solarr, fac, blen);
        });
    } else {
        RT fac = RT(dxinv[0]);
        Box bflux = amrex::surroundingNodes(box, 0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlpoisson_flux_x(tbox, fxarr, solarr, fac);
        });
        fac = RT(dxinv[1]);
        bflux = amrex::surroundingNodes(box, 1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlpoisson_flux_y(tbox, fyarr, solarr, fac);
        });
        fac = RT(dxinv[2]);
        bflux = amrex::surroundingNodes(box, 2);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlpoisson_flux_z(tbox, fzarr, solarr, fac);
        });
    }
#elif (AMREX_SPACEDIM == 2)
    if (face_only) {
        RT fac = RT(dxinv[0]);
        Box blo = amrex::bdryLo(box, 0);
        int blen = box.length(0);
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlpoisson_flux_xface_m(tbox, fxarr, solarr, fac, blen, dx, probxlo);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlpoisson_flux_xface(tbox, fxarr, solarr, fac, blen);
            });
        }
        fac = RT(dxinv[1]);
        blo = amrex::bdryLo(box, 1);
        blen = box.length(1);
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlpoisson_flux_yface_m(tbox, fyarr, solarr, fac, blen, dx, probxlo);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlpoisson_flux_yface(tbox, fyarr, solarr, fac, blen);
            });
        }
    } else {
        RT fac = RT(dxinv[0]);
        Box bflux = amrex::surroundingNodes(box, 0);
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlpoisson_flux_x_m(tbox, fxarr, solarr, fac, dx, probxlo);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlpoisson_flux_x(tbox, fxarr, solarr, fac);
            });
        }
        fac = RT(dxinv[1]);
        bflux = amrex::surroundingNodes(box, 1);
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlpoisson_flux_y_m(tbox, fyarr, solarr, fac, dx, probxlo);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlpoisson_flux_y(tbox, fyarr, solarr, fac);
            });
        }
    }
#else
    if (face_only) {
        RT fac = RT(dxinv[0]);
        Box blo = amrex::bdryLo(box, 0);
        int blen = box.length(0);
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlpoisson_flux_xface_m(tbox, fxarr, solarr, fac, blen, dx, probxlo);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlpoisson_flux_xface(tbox, fxarr, solarr, fac, blen);
            });
        }
    } else {
        RT fac = RT(dxinv[0]);
        Box bflux = amrex::surroundingNodes(box, 0);
        if (this->m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlpoisson_flux_x_m(tbox, fxarr, solarr, fac, dx, probxlo);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlpoisson_flux_x(tbox, fxarr, solarr, fac);
            });
        }
    }
#endif
}

template <typename MF>
bool
MLPoissonT<MF>::supportNSolve () const
{
    bool support = true;
    if (this->m_domain_covered[0]) support = false;
    if (this->doAgglomeration()) support = false;
    if (AMREX_SPACEDIM != 3) support = false;
    return support;
}

template <typename MF>
std::unique_ptr<MLLinOpT<MF>>
MLPoissonT<MF>::makeNLinOp (int grid_size) const
{
    const Geometry& geom = this->m_geom[0].back();
    const BoxArray& ba = this->makeNGrids(grid_size);

    DistributionMapping dm;
    {
        const std::vector<std::vector<int> >& sfc = DistributionMapping::makeSFC(ba);
        Vector<int> pmap(ba.size());
        AMREX_ALWAYS_ASSERT(ParallelContext::CommunicatorSub() == ParallelDescriptor::Communicator());
        const int nprocs = ParallelDescriptor::NProcs();
        for (int iproc = 0; iproc < nprocs; ++iproc) {
            for (int ibox : sfc[iproc]) {
                pmap[ibox] = iproc;
            }
        }
        dm.define(std::move(pmap));
    }

    LPInfo minfo{};
    minfo.has_metric_term = this->info.has_metric_term;

    std::unique_ptr<MLLinOpT<MF>> r{new MLALaplacianT<MF>({geom}, {ba}, {dm}, minfo)};
    auto nop = dynamic_cast<MLALaplacianT<MF>*>(r.get());
    if (!nop) {
        return nullptr;
    }

    nop->m_parent = this;

    nop->setMaxOrder(this->maxorder);
    nop->setVerbose(this->verbose);

    nop->setDomainBC(this->m_lobc, this->m_hibc);

    if (this->needsCoarseDataForBC())
    {
        const Real* dx0 = this->m_geom[0][0].CellSize();
        const Real fac = Real(0.5)*this->m_coarse_data_crse_ratio;
        RealVect cbloc {AMREX_D_DECL(dx0[0]*fac, dx0[1]*fac, dx0[2]*fac)};
        nop->setCoarseFineBCLocation(cbloc);
    }

    nop->setScalars(1.0, -1.0);

    const Real* dxinv = geom.InvCellSize();
    RT dxscale = RT(dxinv[0]);
#if (AMREX_SPACEDIM >= 2)
    dxscale = std::max(dxscale,RT(dxinv[1]));
#endif
#if (AMREX_SPACEDIM == 3)
    dxscale = std::max(dxscale,RT(dxinv[2]));
#endif

    MF alpha(ba, dm, 1, 0);
    alpha.setVal(RT(1.e30)*dxscale*dxscale);

    MF foo(this->m_grids[0].back(), this->m_dmap[0].back(), 1, 0, MFInfo().SetAlloc(false));
    const FabArrayBase::CPC& cpc = alpha.getCPC(IntVect(0),foo,IntVect(0),Periodicity::NonPeriodic());
    alpha.setVal(RT(0.0), cpc, 0, 1);

    nop->setACoeffs(0, alpha);

    return r;
}

template <typename MF>
void
MLPoissonT<MF>::copyNSolveSolution (MF& dst, MF const& src) const
{
    dst.ParallelCopy(src);
}

template <typename MF>
void
MLPoissonT<MF>::get_dpdn_on_domain_faces (Array<MF*,AMREX_SPACEDIM> const& dpdn,
                                          MF const& phi)
{
    BL_PROFILE("MLPoisson::dpdn_faces()");

    // We do not need to call applyBC because this function is used by the
    // OpenBC solver after solver has converged.  That means the BC has been
    // filled to check the residual.

    Box const& domain0 = this->m_geom[0][0].Domain();
    AMREX_D_TERM(const RT dxi = RT(this->m_geom[0][0].InvCellSize(0));,
                 const RT dyi = RT(this->m_geom[0][0].InvCellSize(1));,
                 const RT dzi = RT(this->m_geom[0][0].InvCellSize(2));)

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi);  mfi.isValid(); ++mfi)
    {
        Box const& vbx = mfi.validbox();
        for (OrientationIter oit; oit.isValid(); ++oit) {
            Orientation face = oit();
            if (vbx[face] == domain0[face]) {
                int dir = face.coordDir();
                auto const& p = phi.const_array(mfi);
                auto const& gp = dpdn[dir]->array(mfi);
                Box const& b2d = amrex::bdryNode(vbx,face);
                if (dir == 0) {
                    // because it's dphi/dn, not dphi/dx.
                    RT fac = dxi * (face.isLow() ? RT(-1.0) : RT(1.));
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(b2d, i, j, k,
                    {
                        gp(i,j,k) = fac * (p(i,j,k) - p(i-1,j,k));
                    });
                }
#if (AMREX_SPACEDIM > 1)
                else if (dir == 1) {
                    RT fac = dyi * (face.isLow() ? RT(-1.0) : RT(1.));
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(b2d, i, j, k,
                    {
                        gp(i,j,k) = fac * (p(i,j,k) - p(i,j-1,k));
                    });
                }
#if (AMREX_SPACEDIM > 2)
                else {
                    RT fac = dzi * (face.isLow() ? RT(-1.0) : RT(1.));
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(b2d, i, j, k,
                    {
                        gp(i,j,k) = fac * (p(i,j,k) - p(i,j,k-1));
                    });
                }
#endif
#endif
            }
        }
    }
}

extern template class MLPoissonT<MultiFab>;

using MLPoisson = MLPoissonT<MultiFab>;

}

#endif
