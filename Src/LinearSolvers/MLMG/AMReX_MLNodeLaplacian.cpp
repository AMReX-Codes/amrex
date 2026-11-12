#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLNodeLap_K.H>
#include <AMReX_MultiFabUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_algoim.H>
#endif

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

namespace amrex {

template <typename MF>
MLNodeLaplacianT<MF>::MLNodeLaplacianT (const Vector<Geometry>& a_geom,
                                        const Vector<BoxArray>& a_grids,
                                        const Vector<DistributionMapping>& a_dmap,
                                        const LPInfo& a_info,
                                        const Vector<FabFactory<FAB> const*>& a_factory,
                                        RT  a_const_sigma)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory, a_const_sigma);
}

#ifdef AMREX_USE_EB
template <typename MF>
MLNodeLaplacianT<MF>::MLNodeLaplacianT (const Vector<Geometry>& a_geom,
                                        const Vector<BoxArray>& a_grids,
                                        const Vector<DistributionMapping>& a_dmap,
                                        const LPInfo& a_info,
                                        const Vector<EBFArrayBoxFactory const*>& a_factory,
                                        RT  a_const_sigma)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory, a_const_sigma);
}
#endif

template <typename MF>
void
MLNodeLaplacianT<MF>::define (const Vector<Geometry>& a_geom,
                              const Vector<BoxArray>& a_grids,
                              const Vector<DistributionMapping>& a_dmap,
                              const LPInfo& a_info,
                              const Vector<FabFactory<FAB> const*>& a_factory,
                              RT  a_const_sigma)
{
    BL_PROFILE("MLNodeLaplacian::define()");

    // This makes sure grids are cell-centered;
    Vector<BoxArray> cc_grids = a_grids;
    for (auto& ba : cc_grids) {
        ba.enclosedCells();
    }

    MLNodeLinOpT<MF>::define(a_geom, cc_grids, a_dmap, this->dpInfo(a_info), a_factory);

    m_const_sigma = a_const_sigma;
    m_sigma.resize(this->m_num_amr_levels);
    for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev)
    {
        m_sigma[amrlev].resize(this->m_num_mg_levels[amrlev]);
        const int mglev = 0;
        const int idim = 0;
#ifdef AMREX_USE_EB
        bool allocate_sigma_mfs = true;
#else
        bool allocate_sigma_mfs = m_const_sigma == RT(0.0);
#endif
        if (allocate_sigma_mfs) {
            m_sigma[amrlev][mglev][idim] = std::make_unique<MF>(this->m_grids[amrlev][mglev],
                                                                this->m_dmap[amrlev][mglev],
                                                                1, 1, MFInfo(),
                                                                *this->m_factory[amrlev][0]);
            m_sigma[amrlev][mglev][idim]->setVal(m_const_sigma);
#ifdef AMREX_USE_EB
            m_sigma[amrlev][mglev][idim]->setDomainBndry(0.0, this->m_geom[amrlev][mglev]);
#endif
        }
    }

#ifdef AMREX_USE_EB
#if (AMREX_SPACEDIM == 2)
    const int ncomp_i = 5;
#else
    const int ncomp_i = algoim::numIntgs;
#endif
    m_integral.resize(this->m_num_amr_levels);
    m_surface_integral.resize(this->m_num_amr_levels);
    m_eb_vel_dot_n.resize(this->m_num_amr_levels);
    if constexpr (std::is_same_v<MF,MultiFab>) { // EB needs double precision
        for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev)
        {
            m_integral[amrlev] = std::make_unique<MF>(this->m_grids[amrlev][0],
                                                      this->m_dmap[amrlev][0],
                                                      ncomp_i, 1, MFInfo(),
                                                      *this->m_factory[amrlev][0]);
        }
    } else {
        amrex::ignore_unused(ncomp_i);
    }
#endif

#if (AMREX_SPACEDIM == 2)
    m_is_rz = this->m_geom[0][0].IsRZ();
#endif

    // Mixed precision needs sigma coarsening, and neither it nor the single
    // precision operator has EB.
    if (this->info.mixed_precision || !std::is_same_v<MF,MultiFab>) {
        this->m_coarsening_strategy = CoarseningStrategy::Sigma;
    }

    // The copy's sigma is scaled into the single precision range.
    m_sp_sigma_scale = this->spScaleFor(std::abs(a_const_sigma));
    this->makeSPLinOp(a_info, [&] (auto const& sp_info) {
        return std::make_unique<MLNodeLaplacianT<fMultiFab>>
            (this->spGeom(a_geom), a_grids, a_dmap, sp_info, Vector<FabFactory<BaseFab<float>> const*>{},
             static_cast<float>(m_sp_sigma_scale * a_const_sigma));
    });
    spUpdateScale();

    m_needs_update = true;
}

#ifdef AMREX_USE_EB
template <typename MF>
void
MLNodeLaplacianT<MF>::define (const Vector<Geometry>& a_geom,
                              const Vector<BoxArray>& a_grids,
                              const Vector<DistributionMapping>& a_dmap,
                              const LPInfo& a_info,
                              const Vector<EBFArrayBoxFactory const*>& a_factory,
                              RT  a_const_sigma)
{
    if constexpr (std::is_same_v<MF,MultiFab>) {
        Vector<FabFactory<FAB> const*> _factory;
        for (const auto *x : a_factory) {
            _factory.push_back(static_cast<FabFactory<FAB> const*>(x));
        }
        define(a_geom, a_grids, a_dmap, a_info, _factory, a_const_sigma);
    } else {
        amrex::ignore_unused(a_geom, a_grids, a_dmap, a_info, a_factory, a_const_sigma);
        amrex::Abort("MLNodeLaplacian: EB requires double precision");
    }
}
#endif

template <typename MF>
void
MLNodeLaplacianT<MF>::resizeMultiGrid (int new_size)
{
    if (!   m_sigma.empty()) {
        if (m_sigma[0].size() > new_size) {
            m_sigma[0].resize(new_size);
        }
    }

    if (!   m_stencil.empty()) {
        if (m_stencil[0].size() > new_size) {
            m_stencil[0].resize(new_size);
        }
    }

    if (!   m_s0_norm0.empty()) {
        if (m_s0_norm0[0].size() > new_size) {
            m_s0_norm0[0].resize(new_size);
        }
    }

    MLNodeLinOpT<MF>::resizeMultiGrid(new_size);
}

template <typename MF>
void
MLNodeLaplacianT<MF>::unimposeNeumannBC (int amrlev, MF& rhs) const
{
    if (this->m_coarsening_strategy == CoarseningStrategy::RAP) {
        const Box& nddom = amrex::surroundingNodes(this->Geom(amrlev).Domain());
        const auto lobc = this->LoBC();
        const auto hibc = this->HiBC();

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) { mfi_info.EnableTiling().SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(rhs,mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<RT> const& rhsarr = rhs.array(mfi);
            mlndlap_unimpose_neumann_bc(bx, rhsarr, nddom, lobc, hibc);
        }
    }
}

template <typename MF>
auto
MLNodeLaplacianT<MF>::getSolvabilityOffset (int amrlev, int mglev, MF const& rhs) const -> Vector<RT>
{
    amrex::ignore_unused(amrlev);
    AMREX_ASSERT(amrlev==0 && (mglev+1==this->m_num_mg_levels[0] || mglev==0));
    AMREX_ASSERT(this->getNComp() == 1);

    if (this->m_coarsening_strategy == CoarseningStrategy::RAP) {
#ifdef AMREX_USE_EB
        const auto *factory = ebFactory(amrlev);
        if (mglev == 0 && factory && !factory->isAllRegular()) {
            const MF& vfrac = *ebVolFrac(factory);
            const auto& vfrac_ma = vfrac.const_arrays();

            Box dom = this->Geom(amrlev,mglev).Domain();
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (this->m_lobc[0][idim] != LinOpBCType::Neumann &&
                    this->m_lobc[0][idim] != LinOpBCType::inflow)
                {
                    dom.growLo(idim, 10);
                }
                if (this->m_hibc[0][idim] != LinOpBCType::Neumann &&
                    this->m_hibc[0][idim] != LinOpBCType::inflow)
                {
                    dom.growHi(idim, 10);
                }
            }

            const auto& mask = (mglev+1 == this->m_num_mg_levels[0]) ? this->m_bottom_dot_mask : this->m_coarse_dot_mask;
            const auto& mask_ma = mask.const_arrays();
            const auto& rhs_ma = rhs.const_arrays();
            auto r = ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<RT,RT>{},
                               rhs, IntVect(0),
                               [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                               -> GpuTuple<RT,RT>
                               {
                                   RT scale = RT(0.0);
#if (AMREX_SPACEDIM == 3)
                                   int const koff = 1;
                                   RT const fac = RT(0.125);
#else
                                   int const koff = 0;
                                   RT const fac = RT(0.25);
#endif
                                   for (int kc = k-koff; kc <= k; ++kc) {
                                   for (int jc = j-1   ; jc <= j; ++jc) {
                                   for (int ic = i-1   ; ic <= i; ++ic) {
                                       if (dom.contains(ic,jc,kc)) {
                                           scale += vfrac_ma[box_no](ic,jc,kc) * fac;
                                       }
                                   }}}
                                   return { mask_ma[box_no](i,j,k) * rhs_ma[box_no](i,j,k),
                                            mask_ma[box_no](i,j,k) * scale };
                               });

            RT s1 = amrex::get<0>(r);
            RT s2 = amrex::get<1>(r);
            ParallelAllReduce::Sum<RT>({s1,s2}, ParallelContext::CommunicatorSub());
            return {s1/s2};
        } else
#endif
        {
            Box nddom = amrex::surroundingNodes(this->Geom(amrlev,mglev).Domain());
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (this->m_lobc[0][idim] != LinOpBCType::Neumann &&
                    this->m_lobc[0][idim] != LinOpBCType::inflow)
                {
                    nddom.growLo(idim, 10); // so that the test in ParReduce will fail
                }
                if (this->m_hibc[0][idim] != LinOpBCType::Neumann &&
                    this->m_hibc[0][idim] != LinOpBCType::inflow)
                {
                    nddom.growHi(idim, 10);
                }
            }

            const auto& mask = (mglev+1 == this->m_num_mg_levels[0]) ? this->m_bottom_dot_mask : this->m_coarse_dot_mask;
            const auto& mask_ma = mask.const_arrays();
            const auto& rhs_ma = rhs.const_arrays();
            auto r = ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<RT,RT>{},
                               rhs, IntVect(0),
                               [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                                   -> GpuTuple<RT,RT>
                               {
                                   RT scale = RT(1.0);
                                   if (i == nddom.smallEnd(0) ||
                                       i == nddom.bigEnd(0)) {
                                       scale *= RT(0.5);
                                   }
#if (AMREX_SPACEDIM >= 2)
                                   if (j == nddom.smallEnd(1) ||
                                       j == nddom.bigEnd(1)) {
                                       scale *= RT(0.5);
                                   }
#endif
#if (AMREX_SPACEDIM == 3)
                                   if (k == nddom.smallEnd(2) ||
                                       k == nddom.bigEnd(2)) {
                                       scale *= RT(0.5);
                                   }
#endif
                                   return { mask_ma[box_no](i,j,k) * rhs_ma[box_no](i,j,k),
                                            mask_ma[box_no](i,j,k) * scale };
                               });

            RT s1 = amrex::get<0>(r);
            RT s2 = amrex::get<1>(r);
            ParallelAllReduce::Sum<RT>({s1,s2}, ParallelContext::CommunicatorSub());
            return {s1/s2};
        }
    } else {
        return MLNodeLinOpT<MF>::getSolvabilityOffset(amrlev, mglev, rhs);
    }
}

template <typename MF>
void
MLNodeLaplacianT<MF>::fixSolvabilityByOffset (int amrlev, int mglev, MF& rhs,
                                              Vector<RT> const& a_offset) const
{
    RT offset = a_offset[0];

    if (this->m_coarsening_strategy == CoarseningStrategy::RAP) {
#ifdef AMREX_USE_EB
        const auto *factory = ebFactory(amrlev);
        if (mglev == 0 && factory && !factory->isAllRegular()) {
            const MF& vfrac = *ebVolFrac(factory);
            const auto& vfrac_ma = vfrac.const_arrays();

            Box dom = this->Geom(amrlev,mglev).Domain();
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (this->m_lobc[0][idim] != LinOpBCType::Neumann &&
                    this->m_lobc[0][idim] != LinOpBCType::inflow)
                {
                    dom.growLo(idim, 10);
                }
                if (this->m_hibc[0][idim] != LinOpBCType::Neumann &&
                    this->m_hibc[0][idim] != LinOpBCType::inflow)
                {
                    dom.growHi(idim, 10);
                }
            }

            auto const& rhs_ma = rhs.arrays();
            ParallelFor(rhs, IntVect(0),
                        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                        {
                            RT scale = RT(0.0);
#if (AMREX_SPACEDIM == 3)
                            int const koff = 1;
                            RT const fac = RT(0.125);
#else
                            int const koff = 0;
                            RT const fac = RT(0.25);
#endif
                            for (int kc = k-koff; kc <= k; ++kc) {
                            for (int jc = j-1   ; jc <= j; ++jc) {
                            for (int ic = i-1   ; ic <= i; ++ic) {
                                if (dom.contains(ic,jc,kc)) {
                                    scale += vfrac_ma[box_no](ic,jc,kc) * fac;
                                }
                            }}}
                            rhs_ma[box_no](i,j,k) -= offset * scale;
                        });
        } else
#endif
        {
            Box nddom = amrex::surroundingNodes(this->Geom(amrlev,mglev).Domain());
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (this->m_lobc[0][idim] != LinOpBCType::Neumann &&
                    this->m_lobc[0][idim] != LinOpBCType::inflow)
                {
                    nddom.growLo(idim, 10); // so that the test in ParReduce will fail
                }
                if (this->m_hibc[0][idim] != LinOpBCType::Neumann &&
                    this->m_hibc[0][idim] != LinOpBCType::inflow)
                {
                    nddom.growHi(idim, 10);
                }
            }

            auto const& rhs_ma = rhs.arrays();
            ParallelFor(rhs, IntVect(0),
                        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                        {
                            RT scale = RT(1.0);
                            if (i == nddom.smallEnd(0) ||
                                i == nddom.bigEnd(0)) {
                                scale *= RT(0.5);
                            }
#if (AMREX_SPACEDIM >= 2)
                            if (j == nddom.smallEnd(1) ||
                                j == nddom.bigEnd(1)) {
                                scale *= RT(0.5);
                            }
#endif
#if (AMREX_SPACEDIM == 3)
                            if (k == nddom.smallEnd(2) ||
                                k == nddom.bigEnd(2)) {
                                scale *= RT(0.5);
                            }
#endif
                            rhs_ma[box_no](i,j,k) -= offset * scale;
                        });
        }
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else {
        rhs.plus(-offset, 0, 1);
    }
}

template <typename MF>
void
MLNodeLaplacianT<MF>::setSigma (int amrlev, const MF& a_sigma)
{
    AMREX_ALWAYS_ASSERT(m_sigma[amrlev][0][0]);

    // If we are going to use sigma with AMREX_SPACEDIM components but have only allocated sigma with idim=0 before,
    //    we need to allocate sigma for the other directions here
    if (a_sigma.nComp() > 1)
    {
        AMREX_ALWAYS_ASSERT(a_sigma.nComp() == AMREX_SPACEDIM);
        for (int idim = 1; idim < AMREX_SPACEDIM; idim++) {
            m_sigma[amrlev][0][idim] = std::make_unique<MF>(this->m_grids[amrlev][0],
                                                            this->m_dmap[amrlev][0],
                                                            1, 1, MFInfo());
        }
        setMapped(true);

        for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
            amrex::Copy(*m_sigma[amrlev][0][idim], a_sigma, idim, 0, 1, 0);
        }

    } else {
        amrex::Copy(*m_sigma[amrlev][0][0], a_sigma, 0, 0, 1, 0);
    }

    m_needs_update = true;
}

template <typename MF>
void
MLNodeLaplacianT<MF>::spForwardSigma ()
{
    if constexpr (MLLinOpT<MF>::mixed_precision_capable) {
        // One scale for all AMR levels, from the largest sigma.
        RT smax = RT(0.0);
        for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev) {
            for (auto const& sg : m_sigma[amrlev][0]) {
                if (sg) { smax = std::max(smax, sg->norminf(0, 1, IntVect(0), true)); }
            }
        }
        ParallelAllReduce::Max(smax, ParallelContext::CommunicatorSub());
        m_sp_sigma_scale = this->spScaleFor(smax);
        spUpdateScale();

        int ncomp = 1;
        if (m_use_mapped) { ncomp = AMREX_SPACEDIM; }
        RT const s = m_sp_sigma_scale;
        for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev) {
            fMultiFab fsigma(this->m_grids[amrlev][0], this->m_dmap[amrlev][0], ncomp, 0);
            auto const& fa = fsigma.arrays();
            for (int n = 0; n < ncomp; ++n) {
                auto const& src = m_sigma[amrlev][0][n] ? m_sigma[amrlev][0][n] : m_sigma[amrlev][0][0];
                auto const& sa = src->const_arrays();
                ParallelFor(fsigma, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
                {
                    fa[b](i,j,k,n) = static_cast<float>(s * sa[b](i,j,k));
                });
            }
            Gpu::streamSynchronize();
            this->template spForward<MLNodeLaplacianT>([&] (auto* sp) { sp->setSigma(amrlev, fsigma); });
        }
    }
}

template <typename MF>
void
MLNodeLaplacianT<MF>::FillBoundaryCoeff (MF& sigma, const Geometry& geom)
{
    BL_PROFILE("MLNodeLaplacian::FillBoundaryCoeff()");

    sigma.FillBoundary(geom.periodicity());

    if (this->m_coarsening_strategy == CoarseningStrategy::Sigma)
    {
        const Box& domain = geom.Domain();
        const auto lobc = this->LoBC();
        const auto hibc = this->HiBC();

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) { mfi_info.SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(sigma, mfi_info); mfi.isValid(); ++mfi)
        {
            Array4<RT> const& sfab = sigma.array(mfi);
            mlndlap_fillbc_cc<RT>(mfi.validbox(),sfab,domain,lobc,hibc);
        }
    }
}

template <typename MF>
void
MLNodeLaplacianT<MF>::fixUpResidualMask (int amrlev, iMultiFab& resmsk)
{
    if (!this->m_masks_built) { this->buildMasks(); }

    const iMultiFab& cfmask = *this->m_nd_fine_mask[amrlev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(resmsk,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<int> const& rmsk = resmsk.array(mfi);
        Array4<int const> const& fmsk = cfmask.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
        {
            if (fmsk(i,j,k) == nodelap_detail::crse_fine_node) { rmsk(i,j,k) = 1; }
        });
    }
}

template <typename MF>
void
MLNodeLaplacianT<MF>::prepareForSolve ()
{
    BL_PROFILE("MLNodeLaplacian::prepareForSolve()");

    // The single precision copy would form its RAP coarse operators in single
    // precision, which loses small couplings and slows convergence.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!this->m_sp_linop || this->m_coarsening_strategy != CoarseningStrategy::RAP,
        "MLNodeLaplacian: mixed precision needs CoarseningStrategy::Sigma (setCoarseningStrategy)");

    if (this->m_sp_linop && m_sigma[0][0][0]) { spForwardSigma(); }

    MLNodeLinOpT<MF>::prepareForSolve();

    this->buildMasks();

    averageDownCoeffs();

#ifdef AMREX_USE_EB
    buildIntegral();
    if (m_build_surface_integral) { buildSurfaceIntegral(); }
#endif

    buildStencil();

    m_needs_update = false;
}

template <typename MF>
void
MLNodeLaplacianT<MF>::update ()
{
    BL_PROFILE("MLNodeLaplacian::update()");
    // prepareForSolve rebuilds the sigma-dependent masks, averaged coefficients,
    // and stencil.  The EB integrals are not rebuilt: buildIntegral() and
    // buildSurfaceIntegral() self-guard (m_integral_built / m_surface_integral_built)
    // and are no-ops once they have already been built.
    prepareForSolve();
}

template <typename MF>
void
MLNodeLaplacianT<MF>::restriction (int amrlev, int cmglev, MF& crse, MF& fine) const
{
    BL_PROFILE("MLNodeLaplacian::restriction()");

    this->applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

    IntVect const ratio = (amrlev > 0) ? IntVect(2) : this->mg_coarsen_ratio_vec[cmglev-1];

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MF cfine;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), ratio);
        cfine.define(ba, fine.DistributionMap(), 1, 0, MFInfo().SetArena(The_Async_Arena()));
    }

    MF* pcrse = (need_parallel_copy) ? &cfine : &crse;
    const iMultiFab& dmsk = *this->m_dirichlet_mask[amrlev][cmglev-1];

    const auto& stencil = m_stencil[amrlev][cmglev-1];

    bool regular_coarsening = true;
#if (AMREX_SPACEDIM == 1)
    int idir = 0;
#else
    int idir = 2;
    if (amrlev == 0) {
        regular_coarsening = this->mg_coarsen_ratio_vec[cmglev-1] == this->mg_coarsen_ratio;
        if (ratio[1] == 1) {
            idir = 1;
        } else if (ratio[0] == 1) {
            idir = 0;
        }
    }
#endif

#ifdef AMREX_USE_GPU
    auto pcrse_ma = pcrse->arrays();
    auto fine_ma = fine.const_arrays();
    auto msk_ma = dmsk.const_arrays();

    if (Gpu::inLaunchRegion()) {
        if (this->m_coarsening_strategy == CoarseningStrategy::Sigma)
        {
            if (regular_coarsening)
            {
                ParallelFor(*pcrse, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
                {
                    mlndlap_restriction(i,j,k,pcrse_ma[box_no],fine_ma[box_no],msk_ma[box_no]);
                });
            }
            else
            {
                ParallelFor(*pcrse, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
                {
                    mlndlap_semi_restriction(i,j,k,pcrse_ma[box_no],fine_ma[box_no],msk_ma[box_no],idir);
                });
            }
        }
        else
        {
            auto st_ma = stencil->const_arrays();
            ParallelFor(*pcrse, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                mlndlap_restriction_rap(i,j,k,pcrse_ma[box_no],fine_ma[box_no],st_ma[box_no],msk_ma[box_no]);
            });
        }
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(*pcrse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<RT> cfab = pcrse->array(mfi);
            Array4<RT const> const& ffab = fine.const_array(mfi);
            Array4<int const> const& mfab = dmsk.const_array(mfi);
            if (this->m_coarsening_strategy == CoarseningStrategy::Sigma)
            {
                if (regular_coarsening)
                {
                    amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                    {
                        mlndlap_restriction(i,j,k,cfab,ffab,mfab);
                    });
                }
                else
                {
                    amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                    {
                        mlndlap_semi_restriction(i,j,k,cfab,ffab,mfab,idir);
                    });
                }
            }
            else
            {
                Array4<RT const> const& stfab = stencil->const_array(mfi);
                amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                {
                    mlndlap_restriction_rap(i,j,k,cfab,ffab,stfab,mfab);
                });
            }
        }
    }

    if (need_parallel_copy) {
        crse.ParallelCopy(cfine);
    }
}

template <typename MF>
void
MLNodeLaplacianT<MF>::interpolation (int amrlev, int fmglev, MF& fine, const MF& crse) const
{
    BL_PROFILE("MLNodeLaplacian::interpolation()");

    const auto& sigma = m_sigma[amrlev][fmglev];
    const auto& stencil = m_stencil[amrlev][fmglev];

    IntVect const ratio = (amrlev > 0) ? IntVect(2) : this->mg_coarsen_ratio_vec[fmglev];

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MF cfine;
    const MF* cmf = &crse;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), ratio);
        cfine.define(ba, fine.DistributionMap(), 1, 0, MFInfo().SetArena(The_Async_Arena()));
        cfine.ParallelCopy(crse);
        cmf = &cfine;
    }

    const iMultiFab& dmsk = *this->m_dirichlet_mask[amrlev][fmglev];

    bool regular_coarsening = true;
#if (AMREX_SPACEDIM == 1)
    int idir = 0;
#else
    int idir = 2;
    if (amrlev == 0) {
        regular_coarsening = this->mg_coarsen_ratio_vec[fmglev] == this->mg_coarsen_ratio;
        if (ratio[1] == 1) {
            idir = 1;
        } else if (ratio[0] == 1) {
            idir = 0;
        }
    }
#endif
    if (sigma[0] == nullptr) {
        AMREX_ALWAYS_ASSERT(regular_coarsening);
    }

#ifdef AMREX_USE_GPU
    auto fine_ma = fine.arrays();
    auto crse_ma = cmf->const_arrays();
    auto msk_ma = dmsk.const_arrays();

    if (Gpu::inLaunchRegion()) {
        if (this->m_coarsening_strategy == CoarseningStrategy::RAP)
        {
            auto sten_ma = stencil->const_arrays();
            ParallelFor(fine, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                mlndlap_interpadd_rap(i, j, k, fine_ma[box_no], crse_ma[box_no], sten_ma[box_no], msk_ma[box_no]);
            });
        }
        else if (sigma[0] == nullptr)
        {
            ParallelFor(fine, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                mlndlap_interpadd_c(i, j, k, fine_ma[box_no], crse_ma[box_no], msk_ma[box_no]);
            });
        }
        else if ( (m_use_harmonic_average && fmglev > 0) ||
                   m_use_mapped )
        {
            AMREX_D_TERM(MultiArray4<RT const> const& sx_ma = sigma[0]->const_arrays();,
                         MultiArray4<RT const> const& sy_ma = sigma[1]->const_arrays();,
                         MultiArray4<RT const> const& sz_ma = sigma[2]->const_arrays(););
            ParallelFor(fine, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                mlndlap_interpadd_ha(i, j, k, fine_ma[box_no], crse_ma[box_no], AMREX_D_DECL(sx_ma[box_no], sy_ma[box_no], sz_ma[box_no]), msk_ma[box_no]);
            });
        }
        else if (regular_coarsening)
        {
            auto sig_ma = sigma[0]->const_arrays();
            ParallelFor(fine, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                mlndlap_interpadd_aa(i, j, k, fine_ma[box_no], crse_ma[box_no], sig_ma[box_no], msk_ma[box_no]);
            });
        }
        else
        {
            auto sig_ma = sigma[0]->const_arrays();
            ParallelFor(fine, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                mlndlap_semi_interpadd_aa(i, j, k, fine_ma[box_no], crse_ma[box_no], sig_ma[box_no], msk_ma[box_no], idir);
            });
        }
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(fine, true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<RT> const& ffab = fine.array(mfi);
            Array4<RT const> const& cfab = cmf->const_array(mfi);
            Array4<int const> const& mfab = dmsk.const_array(mfi);
            if (this->m_coarsening_strategy == CoarseningStrategy::RAP)
            {
                Array4<RT const> const& stfab = stencil->const_array(mfi);
                amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                {
                    mlndlap_interpadd_rap(i,j,k,ffab,cfab,stfab,mfab);
                });
            }
            else if (sigma[0] == nullptr)
            {
                amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                {
                    mlndlap_interpadd_c(i,j,k,ffab,cfab,mfab);
                });
            }
            else if ( (m_use_harmonic_average && fmglev > 0) ||
                       m_use_mapped )
            {
                AMREX_D_TERM(Array4<RT const> const& sxfab = sigma[0]->const_array(mfi);,
                             Array4<RT const> const& syfab = sigma[1]->const_array(mfi);,
                             Array4<RT const> const& szfab = sigma[2]->const_array(mfi);)
                amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                {
                    mlndlap_interpadd_ha(i,j,k,ffab,cfab,AMREX_D_DECL(sxfab,syfab,szfab),mfab);
                });
            }
            else
            {
                Array4<RT const> const& sfab = sigma[0]->const_array(mfi);
                if (regular_coarsening)
                {
                    amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                    {
                        mlndlap_interpadd_aa(i,j,k,ffab,cfab,sfab,mfab);
                    });
                }
                else
                {
                    amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                    {
                        mlndlap_semi_interpadd_aa(i,j,k,ffab,cfab,sfab,mfab,idir);
                    });
                }
            }
        }
    }
}

template <typename MF>
void
MLNodeLaplacianT<MF>::averageDownSolutionRHS (int camrlev, MF& crse_sol, MF& crse_rhs,
                                              const MF& fine_sol, const MF& fine_rhs)
{
    const auto& amrrr = this->AMRRefRatio(camrlev);
    amrex::average_down(fine_sol, crse_sol, 0, 1, amrrr);

    if (this->isSingular(0))
    {
        MF frhs(fine_rhs.boxArray(), fine_rhs.DistributionMap(), 1, amrrr-1,
                MFInfo().SetArena(The_Async_Arena()));
        amrex::Copy(frhs, fine_rhs, 0, 0, 1, 0);
        restrictInteriorNodes(camrlev, crse_rhs, frhs);
    }
}

template <typename MF>
void
MLNodeLaplacianT<MF>::interpCoarseFineSolution (int famrlev, MF& fine_sol,
                                                const MF& crse_sol) const
{
    const int ncomp = this->getNComp();

    MF crse_copy = this->makeCoarseAmr(famrlev, IntVect(0), MFInfo().SetArena(The_Async_Arena()));
    crse_copy.setVal(0.0);
    crse_copy.ParallelCopy(crse_sol, 0, 0, ncomp, IntVect(0), IntVect(0),
                           this->Geom(famrlev-1,0).periodicity());
    MF fine_interp(fine_sol.boxArray(), fine_sol.DistributionMap(), ncomp, 0,
                   MFInfo().SetArena(The_Async_Arena()));
    this->interpolationAmr(famrlev, fine_interp, crse_copy, IntVect(0));

    const Box nddom = amrex::surroundingNodes(this->Geom(famrlev,0).Domain());
    const auto lobc = this->LoBC();
    const auto hibc = this->HiBC();
    GpuArray<int,AMREX_SPACEDIM> dlo{}, dhi{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        dlo[d] = (lobc[d] == LinOpBCType::Dirichlet) ? nddom.smallEnd(d) : nddom.smallEnd(d)-1;
        dhi[d] = (hibc[d] == LinOpBCType::Dirichlet) ? nddom.bigEnd(d) : nddom.bigEnd(d)+1;
    }
    auto const& f = fine_sol.arrays();
    auto const& t = fine_interp.const_arrays();
    auto const& dm = this->m_dirichlet_mask[famrlev][0]->const_arrays();
    amrex::ParallelFor(fine_sol, IntVect(0), ncomp,
    [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n) noexcept
    {
        // Not unknowns, nor nodes the user's overset mask marks known
        if (dm[bno](i,j,k) == 0 || dm[bno](i,j,k) == nodelap_detail::overset_node) { return; }
        const IntVect iv(AMREX_D_DECL(i,j,k));
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (iv[d] == dlo[d] || iv[d] == dhi[d]) { return; } // Dirichlet domain node
        }
        f[bno](i,j,k,n) = t[bno](i,j,k,n);
    });
    Gpu::streamSynchronize();
}

template <typename MF>
void
MLNodeLaplacianT<MF>::restrictInteriorNodes (int camrlev, MF& crhs, MF& a_frhs) const
{
    const BoxArray& fba = a_frhs.boxArray();
    const DistributionMapping& fdm = a_frhs.DistributionMap();
    const int amrrr = this->AMRRefRatio(camrlev);

    MF* frhs = nullptr;
    std::unique_ptr<MF> mf;
    if (a_frhs.nGrowVect().allGE(amrrr-1))
    {
        frhs = &a_frhs;
    }
    else
    {
        mf = std::make_unique<MF>(fba, fdm, 1, amrrr-1);
        frhs = mf.get();
        amrex::Copy(*frhs, a_frhs, 0, 0, 1, 0);
    }

    const Geometry& cgeom = this->m_geom[camrlev  ][0];
    const Geometry& fgeom = this->m_geom[camrlev+1][0];

    const Box& f_nd_domain = amrex::surroundingNodes(fgeom.Domain());

    const auto lobc = this->LoBC();
    const auto hibc = this->HiBC();

    const iMultiFab& fdmsk = *this->m_dirichlet_mask[camrlev+1][0];
    const auto& stencil    =  m_nosigma_stencil[camrlev+1];

    MF cfine(amrex::coarsen(fba, amrrr), fdm, 1, 0, MFInfo().SetArena(The_Async_Arena()));

    frhs->setBndry(0.0);
    frhs->FillBoundary(fgeom.periodicity());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cfine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<RT> const& cfab = cfine.array(mfi);
        Array4<RT const> const& ffab = frhs->const_array(mfi);
        Array4<int const> const& mfab = fdmsk.const_array(mfi);
        if (this->m_coarsening_strategy == CoarseningStrategy::Sigma) {
            if (amrrr == 2) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_restriction<2>(i,j,k,cfab,ffab,mfab,f_nd_domain,lobc,hibc);
                });
            } else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_restriction<4>(i,j,k,cfab,ffab,mfab,f_nd_domain,lobc,hibc);
                });
            }
        } else {
            Array4<RT const> const& stfab = stencil->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_restriction_rap(i,j,k,cfab,ffab,stfab,mfab);
            });
        }
    }

    MF tmp_crhs(crhs.boxArray(), crhs.DistributionMap(), 1, 0, MFInfo().SetArena(The_Async_Arena()));
    tmp_crhs.setVal(0.0);
    tmp_crhs.ParallelCopy(cfine, cgeom.periodicity());

    const iMultiFab& c_nd_mask = *this->m_nd_fine_mask[camrlev];
    const auto& has_fine_bndry = *this->m_has_fine_bndry[camrlev];

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) { mfi_info.EnableTiling().SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(crhs, mfi_info); mfi.isValid(); ++mfi)
    {
        if (has_fine_bndry[mfi])
        {
            const Box& bx = mfi.tilebox();
            Array4<RT> const& dfab = crhs.array(mfi);
            Array4<RT const> const& sfab = tmp_crhs.const_array(mfi);
            Array4<int const> const& mfab = c_nd_mask.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                if (mfab(i,j,k) == nodelap_detail::fine_node) { dfab(i,j,k) = sfab(i,j,k); }
            });
        }
    }
}

template <typename MF>
void
MLNodeLaplacianT<MF>::normalize (int amrlev, int mglev, MF& mf) const
{
    BL_PROFILE("MLNodeLaplacian::normalize()");

    if (m_sigma[0][0][0] == nullptr) { return; }

    const auto& sigma = m_sigma[amrlev][mglev];
    const auto& stencil = m_stencil[amrlev][mglev];
    const auto dxinv = this->toRT(this->m_geom[amrlev][mglev].InvCellSizeArray());
    const iMultiFab& dmsk = *this->m_dirichlet_mask[amrlev][mglev];
    const RT s0_norm0 = m_s0_norm0[amrlev][mglev];

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && mf.isFusingCandidate()) {
        const auto& ma = mf.arrays();
        const auto& dmsk_ma = dmsk.const_arrays();

        if (this->m_coarsening_strategy == CoarseningStrategy::RAP)
        {
            const auto& sten_ma = stencil->const_arrays();
            ParallelFor(mf,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                mlndlap_normalize_sten(i,j,k,ma[box_no],sten_ma[box_no],dmsk_ma[box_no],s0_norm0);
            });
        }
        else if ( (m_use_harmonic_average && mglev > 0) ||
                   m_use_mapped )
        {
            AMREX_D_TERM(const auto& sx_ma = sigma[0]->const_arrays();,
                         const auto& sy_ma = sigma[1]->const_arrays();,
                         const auto& sz_ma = sigma[2]->const_arrays(););
            ParallelFor(mf,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                mlndlap_normalize_ha(i,j,k,ma[box_no],AMREX_D_DECL(sx_ma[box_no],sy_ma[box_no],sz_ma[box_no]),dmsk_ma[box_no],dxinv);
            });
        }
        else
        {
            const auto& sx_ma = sigma[0]->const_arrays();
            ParallelFor(mf,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                mlndlap_normalize_aa(i,j,k,ma[box_no],sx_ma[box_no],dmsk_ma[box_no],dxinv);
            });
        }
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<RT> const& arr = mf.array(mfi);
            Array4<int const> const& dmskarr = dmsk.const_array(mfi);
            if (this->m_coarsening_strategy == CoarseningStrategy::RAP)
            {
                Array4<RT const> const& stenarr = stencil->const_array(mfi);

                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_normalize_sten(i,j,k,arr,stenarr,dmskarr,s0_norm0);
                });
            }
            else if ( (m_use_harmonic_average && mglev > 0) ||
                       m_use_mapped )
            {
                AMREX_D_TERM(Array4<RT const> const& sxarr = sigma[0]->const_array(mfi);,
                             Array4<RT const> const& syarr = sigma[1]->const_array(mfi);,
                             Array4<RT const> const& szarr = sigma[2]->const_array(mfi););

                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_normalize_ha(i,j,k,arr,AMREX_D_DECL(sxarr,syarr,szarr),dmskarr,dxinv);
                });
            }
            else
            {
                Array4<RT const> const& sarr = sigma[0]->const_array(mfi);

                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_normalize_aa(i,j,k,arr,sarr,dmskarr,dxinv);
                });
            }
        }
    }
}

template <typename MF>
void
MLNodeLaplacianT<MF>::checkPoint (std::string const& file_name) const
{
    if (ParallelContext::IOProcessorSub())
    {
        UtilCreateCleanDirectory(file_name, false);
        {
            std::string HeaderFileName(file_name+"/Header");
            std::ofstream HeaderFile;
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) {
                FileOpenFailed(HeaderFileName);
            }

            HeaderFile.precision(17);

            // MLLinOp stuff
            HeaderFile << "verbose = " << this->verbose << "\n"
                       << "nlevs = " << this->NAMRLevels() << "\n"
                       << "do_agglomeration = " << this->info.do_agglomeration << "\n"
                       << "do_consolidation = " << this->info.do_consolidation << "\n"
                       << "agg_grid_size = " << this->info.agg_grid_size << "\n"
                       << "con_grid_size = " << this->info.con_grid_size << "\n"
                       << "has_metric_term = " << this->info.has_metric_term << "\n"
                       << "max_coarsening_level = " << this->info.max_coarsening_level << "\n";
#if (AMREX_SPACEDIM == 1)
            HeaderFile << "lobc = " << static_cast<int>(this->m_lobc[0][0]) << "\n";
#elif (AMREX_SPACEDIM == 2)
            HeaderFile << "lobc = " << static_cast<int>(this->m_lobc[0][0])
                       << " "       << static_cast<int>(this->m_lobc[0][1]) << "\n";
#else
            HeaderFile << "lobc = " << static_cast<int>(this->m_lobc[0][0])
                       << " "       << static_cast<int>(this->m_lobc[0][1])
                       << " "       << static_cast<int>(this->m_lobc[0][2]) << "\n";
#endif
#if (AMREX_SPACEDIM == 1)
            HeaderFile << "hibc = " << static_cast<int>(this->m_hibc[0][0]) << "\n";
#elif (AMREX_SPACEDIM == 2)
            HeaderFile << "hibc = " << static_cast<int>(this->m_hibc[0][0])
                       << " "       << static_cast<int>(this->m_hibc[0][1]) << "\n";
#else
            HeaderFile << "hibc = " << static_cast<int>(this->m_hibc[0][0])
                       << " "       << static_cast<int>(this->m_hibc[0][1])
                       << " "       << static_cast<int>(this->m_hibc[0][2]) << "\n";
#endif
            // m_coarse_data_for_bc: not used
            HeaderFile << "maxorder = " << this->getMaxOrder() << "\n";

            // MLNodeLaplacian stuff
            HeaderFile << "is_rz = " << m_is_rz << "\n";
            HeaderFile << "m_const_sigma = " << m_const_sigma << "\n";
            HeaderFile << "use_gauss_seidel = " << m_use_gauss_seidel << "\n";
            HeaderFile << "use_harmonic_average = " << m_use_harmonic_average << "\n";
            HeaderFile << "coarsen_strategy = " << static_cast<int>(this->m_coarsening_strategy) << "\n";
            // No level bc multifab
        }

        for (int ilev = 0; ilev < this->NAMRLevels(); ++ilev)
        {
            UtilCreateCleanDirectory(file_name+"/Level_"+std::to_string(ilev), false);
            std::string HeaderFileName(file_name+"/Level_"+std::to_string(ilev)+"/Header");
            std::ofstream HeaderFile;
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) {
                FileOpenFailed(HeaderFileName);
            }

            HeaderFile.precision(17);

            HeaderFile << this->Geom(ilev) << "\n";
            this->m_grids[ilev][0].writeOn(HeaderFile);  HeaderFile << "\n";
        }
    }

    ParallelContext::BarrierSub();

    if constexpr (std::is_same_v<MF,MultiFab>) {
        for (int ilev = 0; ilev < this->NAMRLevels(); ++ilev)
        {
            if (m_sigma[ilev][0][0]) {
                VisMF::Write(*m_sigma[ilev][0][0], file_name+"/Level_"+std::to_string(ilev)+"/sigma");
            }
        }
    } else {
        amrex::Abort("MLNodeLaplacian::checkPoint: double precision only");
    }
}

#ifdef AMREX_USE_EB
template <typename MF>
void
MLNodeLaplacianT<MF>::setEBInflowVelocity (int amrlev, const MF& eb_vel)
{
    if constexpr (std::is_same_v<MF,MultiFab>) {
        const int mglev = 0;
        if (m_eb_vel_dot_n[amrlev] == nullptr) {
            m_eb_vel_dot_n[amrlev] = std::make_unique<MF>(
                    this->m_grids[amrlev][mglev], this->m_dmap[amrlev][mglev],
                    1, 1, MFInfo(), *this->m_factory[amrlev][mglev]);
        }

        m_eb_vel_dot_n[amrlev]->setVal(0.0);

        const auto *ebfactory = dynamic_cast<EBFArrayBoxFactory const*>(this->m_factory[amrlev][mglev].get());

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) { mfi_info.EnableTiling().SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*m_eb_vel_dot_n[amrlev], mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto& flagfab = ebfactory->getMultiEBCellFlagFab()[mfi];

            if (flagfab.getType(bx) == FabType::singlevalued) {
                Array4<RT> const& eb_vel_dot_n = m_eb_vel_dot_n[amrlev]->array(mfi);
                Array4<RT const> const& ebvelin = eb_vel.const_array(mfi);
                Array4<RT const> const& bnorm = ebfactory->getBndryNormal().const_array(mfi);

                ParallelFor(bx, [eb_vel_dot_n,ebvelin,bnorm]
                 AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    for(int n = 0; n < AMREX_SPACEDIM; ++n)
                    {
                        eb_vel_dot_n(i,j,k) += ebvelin(i,j,k,n)*bnorm(i,j,k,n);
                    }
                });
            }
        }

        m_eb_vel_dot_n[amrlev]->FillBoundary(this->m_geom[amrlev][mglev].periodicity());

#if (AMREX_SPACEDIM == 2)
        const int ncomp_si = 3;
#else
        const int ncomp_si = algoim::numSurfIntgs;
#endif
        if (m_surface_integral[amrlev] == nullptr) {
            m_surface_integral[amrlev] = std::make_unique<MF>(this->m_grids[amrlev][0],
                                                            this->m_dmap[amrlev][0],
                                                            ncomp_si, 1, MFInfo(),
                                                            *this->m_factory[amrlev][0]);
            m_surface_integral_built = false;
        }
        // Turn on flag for building surface integrals
        m_build_surface_integral = true;
    } else {
        amrex::ignore_unused(amrlev, eb_vel);
        amrex::Abort("MLNodeLaplacian::setEBInflowVelocity: EB requires double precision");
    }
}
#endif

template class MLNodeLaplacianT<MultiFab>;
#ifdef AMREX_MLMG_HAS_FLOAT
template class MLNodeLaplacianT<fMultiFab>;
#endif

/// \cond DOXYGEN_IGNORE
namespace detail {
// Many members are defined in the other AMReX_MLNodeLaplacian_*.cpp files.
// With link-time optimization GCC makes the references to them from this
// file weak, and a linker does not take an object file out of a static
// library for a weak reference.  These are strong references to each file.
// Nothing calls this function.
void mlndlap_link_parts ()
{
    mlndlap_anchor_sync();
    mlndlap_anchor_sten();
    mlndlap_anchor_misc();
    mlndlap_anchor_matrix();
#ifdef AMREX_USE_EB
    mlndlap_anchor_eb();
#endif
}
}
/// \endcond

}
