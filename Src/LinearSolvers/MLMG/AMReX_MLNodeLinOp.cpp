
#include <AMReX_MLNodeLinOp.H>
#include <AMReX_MLNodeLinOp_K.H>
#include <AMReX_MLMG_K.H>
#include <AMReX_MultiFabUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

namespace amrex {

MLNodeLinOp::MLNodeLinOp ()
{
    m_ixtype = IntVect::TheNodeVector();
}

void
MLNodeLinOp::define (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info,
                     const Vector<FabFactory<FArrayBox> const*>& a_factory,
                     int a_eb_limit_coarsening)
{
    bool eb_limit_coarsening;
    if (a_eb_limit_coarsening < 0) { // default
#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
        eb_limit_coarsening = true;
#else
        eb_limit_coarsening = false;
#endif
    } else {
        eb_limit_coarsening = a_eb_limit_coarsening;
    }

    MLLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory, eb_limit_coarsening);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!hasHiddenDimension(),
                                     "Nodal solver cannot have any hidden dimensions");

    m_dirichlet_mask.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        m_dirichlet_mask[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_dirichlet_mask[amrlev][mglev] = std::make_unique<iMultiFab>
                (amrex::convert(m_grids[amrlev][mglev],IntVect::TheNodeVector()),
                 m_dmap[amrlev][mglev], 1, 0);
            m_dirichlet_mask[amrlev][mglev]->setVal(0); // non-Dirichlet by default
        }
    }

    m_owner_mask_top = makeOwnerMask(m_grids[0][0],
                                      m_dmap[0][0],
                                      m_geom[0][0]);
    if (m_num_mg_levels[0] == 1) {
        m_owner_mask_bottom = std::make_unique<iMultiFab>(*m_owner_mask_top, amrex::make_alias, 0,
                                                          m_owner_mask_top->nComp());
    } else {
        m_owner_mask_bottom = makeOwnerMask(m_grids[0][m_num_mg_levels[0]-1],
                                             m_dmap[0][m_num_mg_levels[0]-1],
                                             m_geom[0][m_num_mg_levels[0]-1]);
    }

    m_cc_fine_mask.resize(m_num_amr_levels);
    m_nd_fine_mask.resize(m_num_amr_levels);
    m_has_fine_bndry.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        if (amrlev < m_num_amr_levels-1)
        {
            m_nd_fine_mask[amrlev] = std::make_unique<iMultiFab>
                (amrex::convert(m_grids[amrlev][0],IntVect::TheNodeVector()),
                 m_dmap[amrlev][0], 1, 0);
            m_cc_fine_mask[amrlev] = std::make_unique<iMultiFab>
                (m_grids[amrlev][0], m_dmap[amrlev][0], 1, 1);
        } else {
            m_cc_fine_mask[amrlev] = std::make_unique<iMultiFab>
                (m_grids[amrlev][0], m_dmap[amrlev][0], 1, 1, MFInfo().SetAlloc(false));
        }
        m_has_fine_bndry[amrlev] = std::make_unique<LayoutData<int> >(m_grids[amrlev][0],
                                                                      m_dmap[amrlev][0]);
    }

    m_norm_fine_mask.resize(m_num_amr_levels-1);
    for (int amrlev = 0; amrlev < m_num_amr_levels-1; ++amrlev) {
        m_norm_fine_mask[amrlev] = std::make_unique<iMultiFab>
            (makeFineMask(amrex::convert(m_grids[amrlev][0], IntVect(1)), m_dmap[amrlev][0],
                          amrex::convert(m_grids[amrlev+1][0], IntVect(1)),
                          IntVect(m_amr_ref_ratio[amrlev]), 1, 0));
    }
}

void
MLNodeLinOp::prepareForSolve ()
{
    for (int amrlev = 0; amrlev < m_num_amr_levels-1; ++amrlev) {
        fixUpResidualMask(amrlev, *m_norm_fine_mask[amrlev]);
    }
}

std::unique_ptr<iMultiFab>
MLNodeLinOp::makeOwnerMask (const BoxArray& a_ba, const DistributionMapping& dm,
                            const Geometry& geom)
{
    const BoxArray& ba = amrex::convert(a_ba, IntVect::TheNodeVector());
    MultiFab foo(ba,dm,1,0, MFInfo().SetAlloc(false));
    return foo.OwnerMask(geom.periodicity());
}

void
MLNodeLinOp::nodalSync (int amrlev, int mglev, MultiFab& mf) const
{
    mf.OverrideSync(m_geom[amrlev][mglev].periodicity());
}

void
MLNodeLinOp::solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                               const MultiFab* /*crse_bcdata*/)
{
    const int mglev = 0;
    const int ncomp = b.nComp();
    apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, StateMode::Solution);

    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][0];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(resid, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& res = resid.array(mfi);
        Array4<Real const> const& bb = b.const_array(mfi);
        Array4<int const> const& dd = dmsk.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
        {
            if (dd(i,j,k)) {
                res(i,j,k,n) = 0.0;
            } else {
                res(i,j,k,n) = bb(i,j,k,n) - res(i,j,k,n);
            }
        });
    }
}

void
MLNodeLinOp::correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                 BCMode /*bc_mode*/, const MultiFab* /*crse_bcdata*/)
{
    apply(amrlev, mglev, resid, x, BCMode::Homogeneous, StateMode::Correction);
    int ncomp = b.nComp();
    MultiFab::Xpay(resid, -1.0, b, 0, 0, ncomp, 0);
}

void
MLNodeLinOp::apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                    StateMode s_mode, const MLMGBndry*) const
{
    applyBC(amrlev, mglev, in, bc_mode, s_mode);
    Fapply(amrlev, mglev, out, in);
}

void
MLNodeLinOp::smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                     bool skip_fillboundary) const
{
    if (!skip_fillboundary) {
        applyBC(amrlev, mglev, sol, BCMode::Homogeneous, StateMode::Correction);
    }
    Fsmooth(amrlev, mglev, sol, rhs);
}

Real
MLNodeLinOp::xdoty (int amrlev, int mglev, const MultiFab& x, const MultiFab& y, bool local) const
{
    amrex::ignore_unused(amrlev);
    AMREX_ASSERT(amrlev==0);
    AMREX_ASSERT(mglev+1==m_num_mg_levels[0] || mglev==0);
    const auto& mask = (mglev+1 == m_num_mg_levels[0]) ? m_bottom_dot_mask : m_coarse_dot_mask;
    const int ncomp = y.nComp();
    const IntVect nghost(0);
    Real result = amrex::Dot(mask, x, 0, y, 0, ncomp, nghost, true);
    if (!local) {
        ParallelAllReduce::Sum(result, ParallelContext::CommunicatorSub());
    }
    return result;
}

Real
MLNodeLinOp::dotProductPrecond (Vector<MultiFab const*> const& x,
                                Vector<MultiFab const*> const& y) const
{
    Real result = 0;
    const int ncomp = x[0]->nComp();
    for (int ilev = 0; ilev < NAMRLevels(); ++ilev) {
        result += amrex::Dot(m_precond_weight_mask[ilev],
                             *x[ilev],0,*y[ilev],0,ncomp,IntVect(0),true);
    }
    ParallelAllReduce::Sum(result, ParallelContext::CommunicatorSub());
    return result;
}

Real
MLNodeLinOp::norm2Precond (Vector<MultiFab const*> const& x) const
{
    Real result = 0;
    const int ncomp = x[0]->nComp();
    for (int ilev = 0; ilev < NAMRLevels(); ++ilev) {
        result += amrex::Dot(m_precond_weight_mask[ilev],
                             *x[ilev],0,ncomp,IntVect(0),true);
    }
    ParallelAllReduce::Sum(result, ParallelContext::CommunicatorSub());
    return std::sqrt(result);
}

Vector<Real>
MLNodeLinOp::getSolvabilityOffset (int amrlev, int mglev, MultiFab const& rhs) const
{
    amrex::ignore_unused(amrlev);
    AMREX_ASSERT(amrlev==0 && (mglev+1==m_num_mg_levels[0] || mglev==0));
    AMREX_ASSERT(getNComp() == 1);

    const auto& mask = (mglev+1 == m_num_mg_levels[0]) ? m_bottom_dot_mask : m_coarse_dot_mask;
    const auto& mask_ma = mask.const_arrays();
    const auto& rhs_ma = rhs.const_arrays();
    auto r = ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Real,Real>{},
                       rhs, IntVect(0),
                       [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                           -> GpuTuple<Real,Real>
                       {
                           return { mask_ma[box_no](i,j,k) * rhs_ma[box_no](i,j,k),
                                    mask_ma[box_no](i,j,k) };
                       });

    Real s1 = amrex::get<0>(r);
    Real s2 = amrex::get<1>(r);
    ParallelAllReduce::Sum<Real>({s1,s2}, ParallelContext::CommunicatorSub());
    return {s1/s2};
}

void
MLNodeLinOp::fixSolvabilityByOffset (int /*amrlev*/, int /*mglev*/, MultiFab& rhs,
                                     Vector<Real> const& offset) const
{
    rhs.plus(-offset[0], 0, 1);
}

namespace {

void MLNodeLinOp_set_dot_mask (MultiFab& dot_mask, iMultiFab const& omask, Geometry const& geom,
                               GpuArray<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                               GpuArray<LinOpBCType,AMREX_SPACEDIM> const& hibc,
                               MLNodeLinOp::CoarseningStrategy strategy)
{
    Box nddomain = amrex::surroundingNodes(geom.Domain());

    if (strategy != MLNodeLinOp::CoarseningStrategy::Sigma) {
        nddomain.grow(1000); // hack to avoid masks being modified at Neumann boundary
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dot_mask,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& dfab = dot_mask.array(mfi);
        Array4<int const> const& sfab = omask.const_array(mfi);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
        {
            mlndlap_set_dot_mask(tbx, dfab, sfab, nddomain, lobc, hibc);
        });
    }
}

}

void
MLNodeLinOp::buildMasks ()
{
    if (m_masks_built) { return; }

    BL_PROFILE("MLNodeLinOp::buildMasks()");

    m_masks_built = true;

    m_is_bottom_singular = false;
    auto itlo = std::find(m_lobc[0].begin(), m_lobc[0].end(), BCType::Dirichlet); // NOLINT
    auto ithi = std::find(m_hibc[0].begin(), m_hibc[0].end(), BCType::Dirichlet); // NOLINT
    if (itlo == m_lobc[0].end() && ithi == m_hibc[0].end())
    {  // No Dirichlet
        m_is_bottom_singular = (m_domain_covered[0] && !m_overset_dirichlet_mask);
    }

    const auto lobc = LoBC();
    const auto hibc = HiBC();

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            const Geometry& geom = m_geom[amrlev][mglev];
            const auto& period = geom.periodicity();
            const Box& ccdomain = geom.Domain();
            const Box& nddomain = amrex::surroundingNodes(ccdomain);

            auto& dmask = *m_dirichlet_mask[amrlev][mglev];

            iMultiFab ccm(m_grids[amrlev][mglev],m_dmap[amrlev][mglev],1,1);
            ccm.BuildMask(ccdomain,period,0,1,2,0);

            MFItInfo mfi_info;
            if (Gpu::notInLaunchRegion()) { mfi_info.SetDynamic(true); }

            if (m_overset_dirichlet_mask && mglev > 0) {
                const auto& dmask_fine = *m_dirichlet_mask[amrlev][mglev-1];
                amrex::average_down_nodal(dmask_fine, dmask, IntVect(2));
            }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(dmask, mfi_info); mfi.isValid(); ++mfi)
            {
                const Box& ndbx = mfi.validbox();
                Array4<int> const& mskarr = dmask.array(mfi);
                Array4<int const> const& ccarr = ccm.const_array(mfi);
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( ndbx, tbx,
                {
                    mlndlap_set_dirichlet_mask(tbx, mskarr, ccarr, nddomain, lobc, hibc);
                });
            }
        }
    }

    for (int amrlev = 0; amrlev < m_num_amr_levels-1; ++amrlev)
    {
        iMultiFab& cc_mask = *m_cc_fine_mask[amrlev];
        iMultiFab& nd_mask = *m_nd_fine_mask[amrlev];
        LayoutData<int>& has_cf = *m_has_fine_bndry[amrlev];
        const Box& ccdom = m_geom[amrlev][0].Domain();

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMRRefRatio(amrlev) == 2 || AMRRefRatio(amrlev) == 4,
                                         "ref_ratio != 2 and 4 not supported");

        cc_mask = amrex::makeFineMask(cc_mask, *m_cc_fine_mask[amrlev+1], cc_mask.nGrowVect(),
                                      IntVect(AMRRefRatio(amrlev)), m_geom[amrlev][0].periodicity(),
                                      0, 1, has_cf); // coarse: 0, fine: 1

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(cc_mask); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            Array4<int> const& fab = cc_mask.array(mfi);
            mlndlap_fillbc_cc<int>(bx,fab,ccdom,lobc,hibc);
        }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(nd_mask,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<int> const& nmsk = nd_mask.array(mfi);
            Array4<int const> const& cmsk = cc_mask.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
            {
                mlndlap_set_nodal_mask(i,j,k,nmsk,cmsk);
            });
        }
    }

    auto& has_cf = *m_has_fine_bndry[m_num_amr_levels-1];
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(has_cf); mfi.isValid(); ++mfi)
    {
        has_cf[mfi] = 0;
    }

    {
        int amrlev = 0;
        int mglev = m_num_mg_levels[amrlev]-1;
        const Geometry& geom = m_geom[amrlev][mglev];
        const iMultiFab& omask = *m_owner_mask_bottom;
        m_bottom_dot_mask.define(omask.boxArray(), omask.DistributionMap(), 1, 0);
        MLNodeLinOp_set_dot_mask(m_bottom_dot_mask, omask, geom, lobc, hibc, m_coarsening_strategy);
    }

    if (isBottomSingular())
    {
        int amrlev = 0;
        int mglev = 0;
        const Geometry& geom = m_geom[amrlev][mglev];
        const iMultiFab& omask = *m_owner_mask_top;
        m_coarse_dot_mask.define(omask.boxArray(), omask.DistributionMap(), 1, 0);
        MLNodeLinOp_set_dot_mask(m_coarse_dot_mask, omask, geom, lobc, hibc, m_coarsening_strategy);
    }
}

void
MLNodeLinOp::preparePrecond ()
{
    if (m_precond_weight_mask.empty()) {
        m_precond_weight_mask.resize(m_num_amr_levels);
        for (int ilev = 0; ilev < m_num_amr_levels; ++ilev) {
            m_precond_weight_mask[ilev].define(amrex::convert(m_grids[ilev][0],IntVect(1)),
                                               m_dmap[ilev][0], 1, 0);
            auto omask = makeOwnerMask(m_grids[ilev][0],
                                       m_dmap[ilev][0],
                                       m_geom[ilev][0]);
            const auto lobc = LoBC();
            const auto hibc = HiBC();
            Box nddomain = amrex::surroundingNodes(m_geom[ilev][0].Domain());
            if (m_coarsening_strategy != MLNodeLinOp::CoarseningStrategy::Sigma) {
                nddomain.grow(1000); // hack to avoid masks being modified at Neumann boundary
            }
            if (ilev < m_num_amr_levels-1) {
                auto const& fmask = *m_nd_fine_mask[ilev];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(m_precond_weight_mask[ilev],TilingIfNotGPU());
                     mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    Array4<Real> const& dfab = m_precond_weight_mask[ilev].array(mfi);
                    Array4<int const> const& sfab = omask->const_array(mfi);
                    Array4<int const> const& ffab = fmask.const_array(mfi);
                    AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                    {
                        mlndlap_set_dot_mask(tbx, dfab, sfab, ffab, nddomain, lobc, hibc);
                    });
                }
            } else {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(m_precond_weight_mask[ilev],TilingIfNotGPU());
                     mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    Array4<Real> const& dfab = m_precond_weight_mask[ilev].array(mfi);
                    Array4<int const> const& sfab = omask->const_array(mfi);
                    AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                    {
                        mlndlap_set_dot_mask(tbx, dfab, sfab, nddomain, lobc, hibc);
                    });
                }
            }
        }
    }
}

void
MLNodeLinOp::setDirichletNodesToZero (int amrlev, int mglev, MultiFab& mf) const
{
    auto const& maskma = m_dirichlet_mask[amrlev][mglev]->const_arrays();
    auto const& ma = mf.arrays();
    const int ncomp = getNComp();
    ParallelFor(mf, IntVect(0), ncomp,
    [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n)
    {
        if (maskma[bno](i,j,k)) { ma[bno](i,j,k,n) = RT(0.0); }
    });
    Gpu::streamSynchronize();
#ifdef AMREX_USE_EB
    EB_set_covered(mf, 0, ncomp, 0, RT(0.0));
#endif
}

void
MLNodeLinOp::setOversetMask (int amrlev, const iMultiFab& a_dmask)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*m_dirichlet_mask[amrlev][0], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Array4<int const> const& omsk = a_dmask.const_array(mfi);
        Array4<int> const& dmsk = m_dirichlet_mask[amrlev][0]->array(mfi);
        Box const& bx = mfi.tilebox();
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
        {
            dmsk(i,j,k) = 1 - omsk(i,j,k);
        });
    }
    m_overset_dirichlet_mask = true;
}

void
MLNodeLinOp::applyBC (int amrlev, int mglev, MultiFab& phi, BCMode/* bc_mode*/,
                      StateMode state_mode, bool skip_fillboundary) const
{
    BL_PROFILE("MLNodeLinOp::applyBC()");

    m_in_solution_mode = state_mode == StateMode::Solution;

    const Geometry& geom = m_geom[amrlev][mglev];
    const Box& nd_domain = amrex::surroundingNodes(geom.Domain());

    if (!skip_fillboundary) {
        phi.FillBoundary(geom.periodicity());
    }

    if (m_coarsening_strategy == CoarseningStrategy::Sigma)
    {
        const auto lobc = LoBC();
        const auto hibc = HiBC();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(phi); mfi.isValid(); ++mfi)
        {
            Array4<Real> const& fab = phi.array(mfi);
            mlndlap_applybc(mfi.validbox(),fab,nd_domain,lobc,hibc);
        }
    }
}

void
MLNodeLinOp::resizeMultiGrid (int new_size)
{
    if (new_size <= 0 || new_size >= m_num_mg_levels[0]) { return; }

    if (m_dirichlet_mask[0].size() > new_size) {
        m_dirichlet_mask[0].resize(new_size);
    }

    if (m_masks_built)
    {
        const auto lobc = LoBC();
        const auto hibc = HiBC();
        int amrlev = 0;
        int mglev = new_size-1;
        if (mglev == 0) {
            m_owner_mask_bottom = std::make_unique<iMultiFab>(*m_owner_mask_top, amrex::make_alias, 0,
                                                              m_owner_mask_top->nComp());
        } else {
            m_owner_mask_bottom = makeOwnerMask(m_grids[0][mglev],
                                                 m_dmap[0][mglev],
                                                 m_geom[0][mglev]);
        }
        const Geometry& geom = m_geom[amrlev][mglev];
        const iMultiFab& omask = *m_owner_mask_bottom;
        m_bottom_dot_mask = MultiFab();
        m_bottom_dot_mask.define(omask.boxArray(), omask.DistributionMap(), 1, 0);
        MLNodeLinOp_set_dot_mask(m_bottom_dot_mask, omask, geom, lobc, hibc, m_coarsening_strategy);
    }

    MLLinOp::resizeMultiGrid(new_size);
}

Real
MLNodeLinOp::normInf (int amrlev, MultiFab const& mf, bool local) const
{
    const int ncomp = this->getNComp();
    const int finest_level = NAMRLevels() - 1;
    if (amrlev == finest_level) {
        return mf.norminf(0, ncomp, IntVect(0), local);
    } else {
        return mf.norminf(*m_norm_fine_mask[amrlev], 0, ncomp, IntVect(0), local);
    }
}

void
MLNodeLinOp::interpolationAmr (int famrlev, MultiFab& fine, const MultiFab& crse,
                               IntVect const& nghost) const
{
    const int ncomp = getNComp();
    const int refratio = AMRRefRatio(famrlev-1);

    AMREX_ALWAYS_ASSERT(refratio == 2 || refratio == 4);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box fbx = mfi.tilebox();
        fbx.grow(nghost);
        Array4<Real> const& ffab = fine.array(mfi);
        Array4<Real const> const& cfab = crse.const_array(mfi);

        if (refratio == 2) {
            AMREX_HOST_DEVICE_FOR_4D ( fbx, ncomp, i, j, k, n,
            {
                mlmg_lin_nd_interp_r2(i,j,k,n,ffab,cfab);
            });
        } else {
            AMREX_HOST_DEVICE_FOR_4D ( fbx, ncomp, i, j, k, n,
            {
                mlmg_lin_nd_interp_r4(i,j,k,n,ffab,cfab);
            });
        }
    }
}

void
MLNodeLinOp::averageDownAndSync (Vector<MultiFab>& sol) const
{
    const int ncomp = getNComp();
    const int finest_amr_lev = NAMRLevels() - 1;

    nodalSync(finest_amr_lev, 0, sol[finest_amr_lev]);

    for (int falev = finest_amr_lev; falev > 0; --falev)
    {
        const auto& fmf = sol[falev  ];
        auto&       cmf = sol[falev-1];

        auto rr = AMRRefRatio(falev-1);
        MultiFab tmpmf(amrex::coarsen(fmf.boxArray(), rr), fmf.DistributionMap(), ncomp, 0);
        amrex::average_down(fmf, tmpmf, 0, ncomp, rr);
        cmf.ParallelCopy(tmpmf, 0, 0, ncomp);
        nodalSync(falev-1, 0, cmf);
    }
}

void
MLNodeLinOp::interpAssign (int amrlev, int fmglev, MultiFab& fine, MultiFab& crse) const
{
    const int ncomp = getNComp();

    const Geometry& crse_geom = Geom(amrlev,fmglev+1);
    const IntVect refratio = (amrlev > 0) ? IntVect(2) : mg_coarsen_ratio_vec[fmglev];
    AMREX_ALWAYS_ASSERT(refratio == 2);

    MultiFab cfine;
    const MultiFab* cmf;

    if (amrex::isMFIterSafe(crse, fine))
    {
        crse.FillBoundary(crse_geom.periodicity());
        cmf = &crse;
    }
    else
    {
        BoxArray cba = fine.boxArray();
        cba.coarsen(refratio);
        cfine.define(cba, fine.DistributionMap(), ncomp, 0);
        cfine.ParallelCopy(crse, 0, 0, ncomp, 0, 0, crse_geom.periodicity());
        cmf = & cfine;
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& fbx = mfi.tilebox();
        Array4<Real> const& ffab = fine.array(mfi);
        Array4<Real const> const& cfab = cmf->const_array(mfi);

        AMREX_HOST_DEVICE_FOR_4D ( fbx, ncomp, i, j, k, n,
        {
            mlmg_lin_nd_interp_r2(i,j,k,n,ffab,cfab);
        });
    }
}

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
std::unique_ptr<HypreNodeLap>
MLNodeLinOp::makeHypreNodeLap (int bottom_verbose, const std::string& options_namespace) const
{
    const BoxArray& ba = m_grids[0].back();
    const DistributionMapping& dm = m_dmap[0].back();
    const Geometry& geom = m_geom[0].back();
    const auto& factory = *(m_factory[0].back());
    const auto& owner_mask = *m_owner_mask_bottom;
    const auto& dirichlet_mask = *(m_dirichlet_mask[0].back());
    MPI_Comm comm = BottomCommunicator();

    return std::make_unique<amrex::HypreNodeLap>(ba, dm, geom, factory, owner_mask, dirichlet_mask,
                                                 comm, this, bottom_verbose, options_namespace);
}
#endif

}
