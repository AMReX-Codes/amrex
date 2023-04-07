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

MLNodeLaplacian::MLNodeLaplacian (const Vector<Geometry>& a_geom,
                                  const Vector<BoxArray>& a_grids,
                                  const Vector<DistributionMapping>& a_dmap,
                                  const LPInfo& a_info,
                                  const Vector<FabFactory<FArrayBox> const*>& a_factory,
                                  Real  a_const_sigma)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory, a_const_sigma);
}

#ifdef AMREX_USE_EB
MLNodeLaplacian::MLNodeLaplacian (const Vector<Geometry>& a_geom,
                                  const Vector<BoxArray>& a_grids,
                                  const Vector<DistributionMapping>& a_dmap,
                                  const LPInfo& a_info,
                                  const Vector<EBFArrayBoxFactory const*>& a_factory,
                                  Real  a_const_sigma)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory, a_const_sigma);
}
#endif

void
MLNodeLaplacian::define (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info,
                         const Vector<FabFactory<FArrayBox> const*>& a_factory,
                         Real  a_const_sigma)
{
    BL_PROFILE("MLNodeLaplacian::define()");

    // This makes sure grids are cell-centered;
    Vector<BoxArray> cc_grids = a_grids;
    for (auto& ba : cc_grids) {
        ba.enclosedCells();
    }

    MLNodeLinOp::define(a_geom, cc_grids, a_dmap, a_info, a_factory);

    m_const_sigma = a_const_sigma;
    m_sigma.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_sigma[amrlev].resize(m_num_mg_levels[amrlev]);
        const int mglev = 0;
        const int idim = 0;
#ifdef AMREX_USE_EB
        bool allocate_sigma_mfs = true;
#else
        bool allocate_sigma_mfs = m_const_sigma == Real(0.0);
#endif
        if (allocate_sigma_mfs) {
            m_sigma[amrlev][mglev][idim] = std::make_unique<MultiFab>(m_grids[amrlev][mglev],
                                                                      m_dmap[amrlev][mglev],
                                                                      1, 1, MFInfo(),
                                                                      *m_factory[amrlev][0]);
            m_sigma[amrlev][mglev][idim]->setVal(m_const_sigma);
#ifdef AMREX_USE_EB
            m_sigma[amrlev][mglev][idim]->setDomainBndry(0.0, m_geom[amrlev][mglev]);
#endif
        }
    }

#ifdef AMREX_USE_EB
#if (AMREX_SPACEDIM == 2)
    const int ncomp_i = 5;
#else
    const int ncomp_i = algoim::numIntgs;
#endif
    m_integral.resize(m_num_amr_levels);
    m_surface_integral.resize(m_num_amr_levels);
    m_eb_vel_dot_n.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_integral[amrlev] = std::make_unique<MultiFab>(m_grids[amrlev][0],
                                                        m_dmap[amrlev][0],
                                                        ncomp_i, 1, MFInfo(),
                                                        *m_factory[amrlev][0]);
    }
#endif

#if (AMREX_SPACEDIM == 2)
    m_is_rz = m_geom[0][0].IsRZ();
#endif
}

#ifdef AMREX_USE_EB
void
MLNodeLaplacian::define (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info,
                         const Vector<EBFArrayBoxFactory const*>& a_factory,
                         Real  a_const_sigma)
{
    Vector<FabFactory<FArrayBox> const*> _factory;
    for (const auto *x : a_factory) {
        _factory.push_back(static_cast<FabFactory<FArrayBox> const*>(x));
    }
    define(a_geom, a_grids, a_dmap, a_info, _factory, a_const_sigma);
}
#endif

void
MLNodeLaplacian::resizeMultiGrid (int new_size)
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

    MLNodeLinOp::resizeMultiGrid(new_size);
}

void
MLNodeLaplacian::unimposeNeumannBC (int amrlev, MultiFab& rhs) const
{
    if (m_coarsening_strategy == CoarseningStrategy::RAP) {
        const Box& nddom = amrex::surroundingNodes(Geom(amrlev).Domain());
        const auto lobc = LoBC();
        const auto hibc = HiBC();

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(rhs,mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& rhsarr = rhs.array(mfi);
            mlndlap_unimpose_neumann_bc(bx, rhsarr, nddom, lobc, hibc);
        }
    }
}

Vector<Real>
MLNodeLaplacian::getSolvabilityOffset (int amrlev, int mglev, MultiFab const& rhs) const
{
    amrex::ignore_unused(amrlev);
    AMREX_ASSERT(amrlev==0 && (mglev+1==m_num_mg_levels[0] || mglev==0));
    AMREX_ASSERT(getNComp() == 1);

    if (m_coarsening_strategy == CoarseningStrategy::RAP) {
#ifdef AMREX_USE_EB
        const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        if (mglev == 0 && factory && !factory->isAllRegular()) {
            const MultiFab& vfrac = factory->getVolFrac();
            const auto& vfrac_ma = vfrac.const_arrays();

            Box dom = Geom(amrlev,mglev).Domain();
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (m_lobc[0][idim] != LinOpBCType::Neumann &&
                    m_lobc[0][idim] != LinOpBCType::inflow)
                {
                    dom.growLo(idim, 10);
                }
                if (m_hibc[0][idim] != LinOpBCType::Neumann &&
                    m_hibc[0][idim] != LinOpBCType::inflow)
                {
                    dom.growHi(idim, 10);
                }
            }

            const auto& mask = (mglev+1 == m_num_mg_levels[0]) ? m_bottom_dot_mask : m_coarse_dot_mask;
            const auto& mask_ma = mask.const_arrays();
            const auto& rhs_ma = rhs.const_arrays();
            auto r = ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Real,Real>{},
                               rhs, IntVect(0),
                               [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                               -> GpuTuple<Real,Real>
                               {
                                   Real scale = 0.0_rt;
#if (AMREX_SPACEDIM == 3)
                                   int const koff = 1;
                                   Real const fac = 0.125_rt;
#else
                                   int const koff = 0;
                                   Real const fac = 0.25_rt;
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

            Real s1 = amrex::get<0>(r);
            Real s2 = amrex::get<1>(r);
            ParallelAllReduce::Sum<Real>({s1,s2}, ParallelContext::CommunicatorSub());
            return {s1/s2};
        } else
#endif
        {
            Box nddom = amrex::surroundingNodes(Geom(amrlev,mglev).Domain());
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (m_lobc[0][idim] != LinOpBCType::Neumann &&
                    m_lobc[0][idim] != LinOpBCType::inflow)
                {
                    nddom.growLo(idim, 10); // so that the test in ParReduce will faill
                }
                if (m_hibc[0][idim] != LinOpBCType::Neumann &&
                    m_hibc[0][idim] != LinOpBCType::inflow)
                {
                    nddom.growHi(idim, 10);
                }
            }

            const auto& mask = (mglev+1 == m_num_mg_levels[0]) ? m_bottom_dot_mask : m_coarse_dot_mask;
            const auto& mask_ma = mask.const_arrays();
            const auto& rhs_ma = rhs.const_arrays();
            auto r = ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Real,Real>{},
                               rhs, IntVect(0),
                               [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                                   -> GpuTuple<Real,Real>
                               {
                                   Real scale = 1.0_rt;
                                   if (i == nddom.smallEnd(0) ||
                                       i == nddom.bigEnd(0)) {
                                       scale *= 0.5_rt;
                                   }
#if (AMREX_SPACEDIM >= 2)
                                   if (j == nddom.smallEnd(1) ||
                                       j == nddom.bigEnd(1)) {
                                       scale *= 0.5_rt;
                                   }
#endif
#if (AMREX_SPACEDIM == 3)
                                   if (k == nddom.smallEnd(2) ||
                                       k == nddom.bigEnd(2)) {
                                       scale *= 0.5_rt;
                                   }
#endif
                                   return { mask_ma[box_no](i,j,k) * rhs_ma[box_no](i,j,k),
                                            mask_ma[box_no](i,j,k) * scale };
                               });

            Real s1 = amrex::get<0>(r);
            Real s2 = amrex::get<1>(r);
            ParallelAllReduce::Sum<Real>({s1,s2}, ParallelContext::CommunicatorSub());
            return {s1/s2};
        }
    } else {
        return MLNodeLinOp::getSolvabilityOffset(amrlev, mglev, rhs);
    }
}

void
MLNodeLaplacian::fixSolvabilityByOffset (int amrlev, int mglev, MultiFab& rhs,
                                         Vector<Real> const& a_offset) const
{
    Real offset = a_offset[0];

    if (m_coarsening_strategy == CoarseningStrategy::RAP) {
#ifdef AMREX_USE_EB
        const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        if (mglev == 0 && factory && !factory->isAllRegular()) {
            const MultiFab& vfrac = factory->getVolFrac();
            const auto& vfrac_ma = vfrac.const_arrays();

            Box dom = Geom(amrlev,mglev).Domain();
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (m_lobc[0][idim] != LinOpBCType::Neumann &&
                    m_lobc[0][idim] != LinOpBCType::inflow)
                {
                    dom.growLo(idim, 10);
                }
                if (m_hibc[0][idim] != LinOpBCType::Neumann &&
                    m_hibc[0][idim] != LinOpBCType::inflow)
                {
                    dom.growHi(idim, 10);
                }
            }

            auto const& rhs_ma = rhs.arrays();
            ParallelFor(rhs, IntVect(0),
                        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                        {
                            Real scale = 0.0_rt;
#if (AMREX_SPACEDIM == 3)
                            int const koff = 1;
                            Real const fac = 0.125_rt;
#else
                            int const koff = 0;
                            Real const fac = 0.25_rt;
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
            Box nddom = amrex::surroundingNodes(Geom(amrlev,mglev).Domain());
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (m_lobc[0][idim] != LinOpBCType::Neumann &&
                    m_lobc[0][idim] != LinOpBCType::inflow)
                {
                    nddom.growLo(idim, 10); // so that the test in ParReduce will faill
                }
                if (m_hibc[0][idim] != LinOpBCType::Neumann &&
                    m_hibc[0][idim] != LinOpBCType::inflow)
                {
                    nddom.growHi(idim, 10);
                }
            }

            auto const& rhs_ma = rhs.arrays();
            ParallelFor(rhs, IntVect(0),
                        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                        {
                            Real scale = 1.0_rt;
                            if (i == nddom.smallEnd(0) ||
                                i == nddom.bigEnd(0)) {
                                scale *= 0.5_rt;
                            }
#if (AMREX_SPACEDIM >= 2)
                            if (j == nddom.smallEnd(1) ||
                                j == nddom.bigEnd(1)) {
                                scale *= 0.5_rt;
                            }
#endif
#if (AMREX_SPACEDIM == 3)
                            if (k == nddom.smallEnd(2) ||
                                k == nddom.bigEnd(2)) {
                                scale *= 0.5_rt;
                            }
#endif
                            rhs_ma[box_no](i,j,k) -= offset * scale;
                        });
        }
        Gpu::streamSynchronize();
    } else {
        rhs.plus(-offset, 0, 1);
    }
}

void
MLNodeLaplacian::setSigma (int amrlev, const MultiFab& a_sigma)
{
    AMREX_ALWAYS_ASSERT(m_sigma[amrlev][0][0]);

    // If we are going to use sigma with AMREX_SPACEDIM components but have only allocated sigma with idim=0 before,
    //    we need to allocate sigma for the other directions here
    if (a_sigma.nComp() > 1)
    {
        AMREX_ALWAYS_ASSERT(a_sigma.nComp() == AMREX_SPACEDIM);
        for (int idim = 1; idim < AMREX_SPACEDIM; idim++)
            m_sigma[amrlev][0][idim] = std::make_unique<MultiFab>(m_grids[amrlev][0],
                                                                  m_dmap[amrlev][0],
                                                                  1, 1, MFInfo());
        setMapped(true);

        for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
            MultiFab::Copy(*m_sigma[amrlev][0][idim], a_sigma, idim, 0, 1, 0);

    } else {
        MultiFab::Copy(*m_sigma[amrlev][0][0], a_sigma, 0, 0, 1, 0);
    }
}

void
MLNodeLaplacian::FillBoundaryCoeff (MultiFab& sigma, const Geometry& geom)
{
    BL_PROFILE("MLNodeLaplacian::FillBoundaryCoeff()");

    sigma.FillBoundary(geom.periodicity());

    if (m_coarsening_strategy == CoarseningStrategy::Sigma)
    {
        const Box& domain = geom.Domain();
        const auto lobc = LoBC();
        const auto hibc = HiBC();

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(sigma, mfi_info); mfi.isValid(); ++mfi)
        {
            Array4<Real> const& sfab = sigma.array(mfi);
            mlndlap_fillbc_cc<Real>(mfi.validbox(),sfab,domain,lobc,hibc);
        }
    }
}

void
MLNodeLaplacian::fixUpResidualMask (int amrlev, iMultiFab& resmsk)
{
    if (!m_masks_built) buildMasks();

    const iMultiFab& cfmask = *m_nd_fine_mask[amrlev];

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
            if (fmsk(i,j,k) == crse_fine_node) rmsk(i,j,k) = 1;
        });
    }
}

void
MLNodeLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLNodeLaplacian::prepareForSolve()");

    MLNodeLinOp::prepareForSolve();

    buildMasks();

    averageDownCoeffs();

#ifdef AMREX_USE_EB
    buildIntegral();
    if (m_build_surface_integral) buildSurfaceIntegral();
#endif

    buildStencil();
}

void
MLNodeLaplacian::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
    BL_PROFILE("MLNodeLaplacian::restriction()");

    applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
        cfine.define(ba, fine.DistributionMap(), 1, 0);
    }

    MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][cmglev-1];

    const auto& stencil = m_stencil[amrlev][cmglev-1];

    bool regular_coarsening = true;
#if (AMREX_SPACEDIM == 1)
    int idir = 0;
#else
    int idir = 2;
    if (cmglev > 0) {
        regular_coarsening = mg_coarsen_ratio_vec[cmglev-1] == mg_coarsen_ratio;
        IntVect ratio = mg_coarsen_ratio_vec[cmglev-1];
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
        if (m_coarsening_strategy == CoarseningStrategy::Sigma)
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
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(*pcrse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> cfab = pcrse->array(mfi);
            Array4<Real const> const& ffab = fine.const_array(mfi);
            Array4<int const> const& mfab = dmsk.const_array(mfi);
            if (m_coarsening_strategy == CoarseningStrategy::Sigma)
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
                Array4<Real const> const& stfab = stencil->const_array(mfi);
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

void
MLNodeLaplacian::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
    BL_PROFILE("MLNodeLaplacian::interpolation()");

    const auto& sigma = m_sigma[amrlev][fmglev];
    const auto& stencil = m_stencil[amrlev][fmglev];

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    const MultiFab* cmf = &crse;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
        cfine.define(ba, fine.DistributionMap(), 1, 0);
        cfine.ParallelCopy(crse);
        cmf = &cfine;
    }

    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][fmglev];

    bool regular_coarsening = true;
#if (AMREX_SPACEDIM == 1)
    int idir = 0;
#else
    int idir = 2;
    if (fmglev > 0) {
        regular_coarsening = mg_coarsen_ratio_vec[fmglev] == mg_coarsen_ratio;
        IntVect ratio = mg_coarsen_ratio_vec[fmglev];
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
        if (m_coarsening_strategy == CoarseningStrategy::RAP)
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
            AMREX_D_TERM(MultiArray4<Real const> const& sx_ma = sigma[0]->const_arrays();,
                         MultiArray4<Real const> const& sy_ma = sigma[1]->const_arrays();,
                         MultiArray4<Real const> const& sz_ma = sigma[2]->const_arrays(););
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
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(fine, true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<Real> const& ffab = fine.array(mfi);
            Array4<Real const> const& cfab = cmf->const_array(mfi);
            Array4<int const> const& mfab = dmsk.const_array(mfi);
            if (m_coarsening_strategy == CoarseningStrategy::RAP)
            {
                Array4<Real const> const& stfab = stencil->const_array(mfi);
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
                AMREX_D_TERM(Array4<Real const> const& sxfab = sigma[0]->const_array(mfi);,
                             Array4<Real const> const& syfab = sigma[1]->const_array(mfi);,
                             Array4<Real const> const& szfab = sigma[2]->const_array(mfi);)
                amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                {
                    mlndlap_interpadd_ha(i,j,k,ffab,cfab,AMREX_D_DECL(sxfab,syfab,szfab),mfab);
                });
            }
            else
            {
                Array4<Real const> const& sfab = sigma[0]->const_array(mfi);
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

void
MLNodeLaplacian::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs)
{
    const auto& amrrr = AMRRefRatio(camrlev);
    amrex::average_down(fine_sol, crse_sol, 0, 1, amrrr);

    if (isSingular(0))
    {
        MultiFab frhs(fine_rhs.boxArray(), fine_rhs.DistributionMap(), 1, amrrr-1);
        MultiFab::Copy(frhs, fine_rhs, 0, 0, 1, 0);
        restrictInteriorNodes(camrlev, crse_rhs, frhs);
    }
}

void
MLNodeLaplacian::restrictInteriorNodes (int camrlev, MultiFab& crhs, MultiFab& a_frhs) const
{
    const BoxArray& fba = a_frhs.boxArray();
    const DistributionMapping& fdm = a_frhs.DistributionMap();
    const int amrrr = AMRRefRatio(camrlev);

    MultiFab* frhs = nullptr;
    std::unique_ptr<MultiFab> mf;
    if (a_frhs.nGrowVect().allGE(IntVect(amrrr-1)))
    {
        frhs = &a_frhs;
    }
    else
    {
        mf = std::make_unique<MultiFab>(fba, fdm, 1, amrrr-1);
        frhs = mf.get();
        MultiFab::Copy(*frhs, a_frhs, 0, 0, 1, 0);
    }

    const Geometry& cgeom = m_geom[camrlev  ][0];
    const Geometry& fgeom = m_geom[camrlev+1][0];

    const Box& f_nd_domain = amrex::surroundingNodes(fgeom.Domain());

    const auto lobc = LoBC();
    const auto hibc = HiBC();

    const iMultiFab& fdmsk = *m_dirichlet_mask[camrlev+1][0];
    const auto& stencil    =  m_nosigma_stencil[camrlev+1];

    MultiFab cfine(amrex::coarsen(fba, amrrr), fdm, 1, 0);

    frhs->setBndry(0.0);
    frhs->FillBoundary(fgeom.periodicity());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cfine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& cfab = cfine.array(mfi);
        Array4<Real const> const& ffab = frhs->const_array(mfi);
        Array4<int const> const& mfab = fdmsk.const_array(mfi);
        if (m_coarsening_strategy == CoarseningStrategy::Sigma) {
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
            Array4<Real const> const& stfab = stencil->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_restriction_rap(i,j,k,cfab,ffab,stfab,mfab);
            });
        }
    }

    MultiFab tmp_crhs(crhs.boxArray(), crhs.DistributionMap(), 1, 0);
    tmp_crhs.setVal(0.0);
    tmp_crhs.ParallelCopy(cfine, cgeom.periodicity());

    const iMultiFab& c_nd_mask = *m_nd_fine_mask[camrlev];
    const auto& has_fine_bndry = *m_has_fine_bndry[camrlev];

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(crhs, mfi_info); mfi.isValid(); ++mfi)
    {
        if (has_fine_bndry[mfi])
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& dfab = crhs.array(mfi);
            Array4<Real const> const& sfab = tmp_crhs.const_array(mfi);
            Array4<int const> const& mfab = c_nd_mask.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                if (mfab(i,j,k) == fine_node) dfab(i,j,k) = sfab(i,j,k);
            });
        }
    }
}

void
MLNodeLaplacian::normalize (int amrlev, int mglev, MultiFab& mf) const
{
    BL_PROFILE("MLNodeLaplacian::normalize()");

    if (m_sigma[0][0][0] == nullptr) return;

    const auto& sigma = m_sigma[amrlev][mglev];
    const auto& stencil = m_stencil[amrlev][mglev];
    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];
    const Real s0_norm0 = m_s0_norm0[amrlev][mglev];

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && mf.isFusingCandidate()) {
        const auto& ma = mf.arrays();
        const auto& dmsk_ma = dmsk.const_arrays();

        if (m_coarsening_strategy == CoarseningStrategy::RAP)
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
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& arr = mf.array(mfi);
            Array4<int const> const& dmskarr = dmsk.const_array(mfi);
            if (m_coarsening_strategy == CoarseningStrategy::RAP)
            {
                Array4<Real const> const& stenarr = stencil->const_array(mfi);

                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_normalize_sten(i,j,k,arr,stenarr,dmskarr,s0_norm0);
                });
            }
            else if ( (m_use_harmonic_average && mglev > 0) ||
                       m_use_mapped )
            {
                AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                             Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                             Array4<Real const> const& szarr = sigma[2]->const_array(mfi););

                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_normalize_ha(i,j,k,arr,AMREX_D_DECL(sxarr,syarr,szarr),dmskarr,dxinv);
                });
            }
            else
            {
                Array4<Real const> const& sarr = sigma[0]->const_array(mfi);

                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_normalize_aa(i,j,k,arr,sarr,dmskarr,dxinv);
                });
            }
        }
    }
}

void
MLNodeLaplacian::checkPoint (std::string const& file_name) const
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

            // MLLinop stuff
            HeaderFile << "verbose = " << verbose << "\n"
                       << "nlevs = " << NAMRLevels() << "\n"
                       << "do_agglomeration = " << info.do_agglomeration << "\n"
                       << "do_consolidation = " << info.do_consolidation << "\n"
                       << "agg_grid_size = " << info.agg_grid_size << "\n"
                       << "con_grid_size = " << info.con_grid_size << "\n"
                       << "has_metric_term = " << info.has_metric_term << "\n"
                       << "max_coarsening_level = " << info.max_coarsening_level << "\n";
#if (AMREX_SPACEDIM == 1)
            HeaderFile << "lobc = " << static_cast<int>(m_lobc[0][0]) << "\n";
#elif (AMREX_SPACEDIM == 2)
            HeaderFile << "lobc = " << static_cast<int>(m_lobc[0][0])
                       << " "       << static_cast<int>(m_lobc[0][1]) << "\n";
#else
            HeaderFile << "lobc = " << static_cast<int>(m_lobc[0][0])
                       << " "       << static_cast<int>(m_lobc[0][1])
                       << " "       << static_cast<int>(m_lobc[0][2]) << "\n";
#endif
#if (AMREX_SPACEDIM == 1)
            HeaderFile << "hibc = " << static_cast<int>(m_hibc[0][0]) << "\n";
#elif (AMREX_SPACEDIM == 2)
            HeaderFile << "hibc = " << static_cast<int>(m_hibc[0][0])
                       << " "       << static_cast<int>(m_hibc[0][1]) << "\n";
#else
            HeaderFile << "hibc = " << static_cast<int>(m_hibc[0][0])
                       << " "       << static_cast<int>(m_hibc[0][1])
                       << " "       << static_cast<int>(m_hibc[0][2]) << "\n";
#endif
            // m_coarse_data_for_bc: not used
            HeaderFile << "maxorder = " << getMaxOrder() << "\n";

            // MLNodeLaplacian stuff
            HeaderFile << "is_rz = " << m_is_rz << "\n";
            HeaderFile << "m_const_sigma = " << m_const_sigma << "\n";
            HeaderFile << "use_gauss_seidel = " << m_use_gauss_seidel << "\n";
            HeaderFile << "use_harmonic_average = " << m_use_harmonic_average << "\n";
            HeaderFile << "coarsen_strategy = " << static_cast<int>(m_coarsening_strategy) << "\n";
            // No level bc multifab
        }

        for (int ilev = 0; ilev < NAMRLevels(); ++ilev)
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

            HeaderFile << Geom(ilev) << "\n";
            m_grids[ilev][0].writeOn(HeaderFile);  HeaderFile << "\n";
        }
    }

    ParallelContext::BarrierSub();

    for (int ilev = 0; ilev < NAMRLevels(); ++ilev)
    {
        if (m_sigma[ilev][0][0]) {
            VisMF::Write(*m_sigma[ilev][0][0], file_name+"/Level_"+std::to_string(ilev)+"/sigma");
        }
    }
}

#ifdef AMREX_USE_EB
void
MLNodeLaplacian::setEBInflowVelocity (int amrlev, const MultiFab& eb_vel)
{
    const int mglev = 0;
    if (m_eb_vel_dot_n[amrlev] == nullptr) {
        m_eb_vel_dot_n[amrlev] = std::make_unique<MultiFab>(
                m_grids[amrlev][mglev], m_dmap[amrlev][mglev],
                1, 1, MFInfo(), *m_factory[amrlev][mglev]);
    }

    m_eb_vel_dot_n[amrlev]->setVal(0.0);

    const auto *ebfactory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_eb_vel_dot_n[amrlev], mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& flagfab = ebfactory->getMultiEBCellFlagFab()[mfi];

        if (flagfab.getType(bx) == FabType::singlevalued) {
            Array4<Real> const& eb_vel_dot_n = m_eb_vel_dot_n[amrlev]->array(mfi);
            Array4<Real const> const& ebvelin = eb_vel.const_array(mfi);
            Array4<Real const> const& bnorm = ebfactory->getBndryNormal().const_array(mfi);

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

    m_eb_vel_dot_n[amrlev]->FillBoundary(m_geom[amrlev][mglev].periodicity());

#if (AMREX_SPACEDIM == 2)
    const int ncomp_si = 3;
#else
    const int ncomp_si = algoim::numSurfIntgs;
#endif
    m_surface_integral[amrlev] = std::make_unique<MultiFab>(m_grids[amrlev][0],
                                                    m_dmap[amrlev][0],
                                                    ncomp_si, 1, MFInfo(),
                                                    *m_factory[amrlev][0]);
    // Turn on flag for building surface integrals
    m_build_surface_integral = true;
}
#endif

}
