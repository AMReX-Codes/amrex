
#include <AMReX_MLNodeLinOp.H>
#include <AMReX_MLNodeLap_K.H>
#include <AMReX_MultiFabUtil.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

MLNodeLinOp::MLNodeLinOp ()
{
    m_ixtype = IntVect::TheNodeVector();
}

MLNodeLinOp::~MLNodeLinOp () {}

void
MLNodeLinOp::define (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info,
                     const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
#ifdef AMREX_USE_HYPRE
    bool eb_limit_coarsening = true;
#else
    bool eb_limit_coarsening = false;
#endif
    MLLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory, eb_limit_coarsening);

    m_owner_mask.resize(m_num_amr_levels);
    m_dirichlet_mask.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        m_owner_mask[amrlev].resize(m_num_mg_levels[amrlev]);
        m_dirichlet_mask[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_owner_mask[amrlev][mglev] = makeOwnerMask(m_grids[amrlev][mglev],
                                                        m_dmap[amrlev][mglev],
                                                        m_geom[amrlev][mglev]);
            m_dirichlet_mask[amrlev][mglev].reset
                (new iMultiFab(amrex::convert(m_grids[amrlev][mglev],IntVect::TheNodeVector()),
                               m_dmap[amrlev][mglev], 1, 0));
            m_dirichlet_mask[amrlev][mglev]->setVal(0); // non-Dirichlet by default
        }
    }


    m_cc_fine_mask.resize(m_num_amr_levels);
    m_nd_fine_mask.resize(m_num_amr_levels);
    m_has_fine_bndry.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        if (amrlev < m_num_amr_levels-1)
        {
            m_nd_fine_mask[amrlev].reset(new iMultiFab(amrex::convert(m_grids[amrlev][0],IntVect::TheNodeVector()),
                                                       m_dmap[amrlev][0], 1, 0));
            m_cc_fine_mask[amrlev].reset(new iMultiFab(m_grids[amrlev][0], m_dmap[amrlev][0], 1, 1));
        } else {
            m_cc_fine_mask[amrlev].reset(new iMultiFab(m_grids[amrlev][0], m_dmap[amrlev][0], 1, 1,
                                                       MFInfo().SetAlloc(false)));
        }
        m_has_fine_bndry[amrlev].reset(new LayoutData<int>(m_grids[amrlev][0], m_dmap[amrlev][0]));
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
    mf.OverrideSync(*m_owner_mask[amrlev][mglev], m_geom[amrlev][mglev].periodicity());
}

void
MLNodeLinOp::solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                               const MultiFab* /*crse_bcdata*/)
{
    const int mglev = 0;
    const int ncomp = b.nComp();
    apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, StateMode::Solution);

    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][0];
#ifdef _OPENMP
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
        applyBC(amrlev, mglev, sol, BCMode::Homogeneous, StateMode::Solution);
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
    const int nghost = 0;
    MultiFab tmp(x.boxArray(), x.DistributionMap(), ncomp, 0);
    MultiFab::Copy(tmp, x, 0, 0, ncomp, nghost);
    for (int i = 0; i < ncomp; i++) {
        MultiFab::Multiply(tmp, mask, 0, i, 1, nghost);
    }
    Real result = MultiFab::Dot(tmp,0,y,0,ncomp,nghost,true);
    if (!local) {
        ParallelAllReduce::Sum(result, ParallelContext::CommunicatorSub());
    }
    return result;
}

void
MLNodeLinOp::applyInhomogNeumannTerm (int amrlev, MultiFab& rhs) const
{
    amrex::ignore_unused(amrlev);
    int ncomp = rhs.nComp();
    for (int n = 0; n < ncomp; ++n)
    {
        auto itlo = std::find(m_lo_inhomog_neumann[n].begin(),
                              m_lo_inhomog_neumann[n].end(),   1);
        auto ithi = std::find(m_hi_inhomog_neumann[n].begin(),
                              m_hi_inhomog_neumann[n].end(),   1);
        if (itlo != m_lo_inhomog_neumann[n].end() or
            ithi != m_hi_inhomog_neumann[n].end())
        {
            amrex::Abort("Inhomogeneous Neumann not supported for nodal solver");
        }
    }
}

namespace {

void MLNodeLinOp_set_dot_mask (MultiFab& dot_mask, iMultiFab const& omask, Geometry const& geom,
                               GpuArray<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                               GpuArray<LinOpBCType,AMREX_SPACEDIM> const& hibc,
                               MLNodeLinOp::CoarseningStrategy strategy)
{
    Box nddomain = amrex::surroundingNodes(geom.Domain());

    if (strategy != MLNodeLinOp::CoarseningStrategy::Sigma) {
        nddomain.grow(1000); // hack to avoid masks being modified at Neuman boundary
    }

#ifdef _OPENMP
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
    if (m_masks_built) return;

    BL_PROFILE("MLNodeLinOp::buildMasks()");

    m_masks_built = true;

    m_is_bottom_singular = false;
    auto itlo = std::find(m_lobc[0].begin(), m_lobc[0].end(), BCType::Dirichlet);
    auto ithi = std::find(m_hibc[0].begin(), m_hibc[0].end(), BCType::Dirichlet);
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
            if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);

            if (m_overset_dirichlet_mask and mglev > 0) {
                const auto& dmask_fine = *m_dirichlet_mask[amrlev][mglev-1];
                amrex::average_down_nodal(dmask_fine, dmask, IntVect(2));
            }
#ifdef _OPENMP
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

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMRRefRatio(amrlev) == 2, "ref_ratio != 0 not supported");

        cc_mask = amrex::makeFineMask(cc_mask, *m_cc_fine_mask[amrlev+1], cc_mask.nGrowVect(),
                                      IntVect(AMRRefRatio(amrlev)), m_geom[amrlev][0].periodicity(),
                                      0, 1, has_cf); // coarse: 0, fine: 1

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(cc_mask); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            Array4<int> const& fab = cc_mask.array(mfi);
            mlndlap_fillbc_cc<int>(bx,fab,ccdom,lobc,hibc);
        }

#ifdef _OPENMP
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
#ifdef _OPENMP
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
        const iMultiFab& omask = *m_owner_mask[amrlev][mglev];
        m_bottom_dot_mask.define(omask.boxArray(), omask.DistributionMap(), 1, 0);
        MLNodeLinOp_set_dot_mask(m_bottom_dot_mask, omask, geom, lobc, hibc, m_coarsening_strategy);
    }

    if (m_is_bottom_singular)
    {
        int amrlev = 0;
        int mglev = 0;
        const Geometry& geom = m_geom[amrlev][mglev];
        const iMultiFab& omask = *m_owner_mask[amrlev][mglev];
        m_coarse_dot_mask.define(omask.boxArray(), omask.DistributionMap(), 1, 0);
        MLNodeLinOp_set_dot_mask(m_coarse_dot_mask, omask, geom, lobc, hibc, m_coarsening_strategy);
    }
}

void
MLNodeLinOp::setOversetMask (int amrlev, const iMultiFab& a_dmask)
{
#ifdef _OPENMP
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
MLNodeLinOp::applyBC (int amrlev, int mglev, MultiFab& phi, BCMode/* bc_mode*/, StateMode,
                      bool skip_fillboundary) const
{
    BL_PROFILE("MLNodeLinOp::applyBC()");

    const Geometry& geom = m_geom[amrlev][mglev];
    const Box& nd_domain = amrex::surroundingNodes(geom.Domain());

    if (!skip_fillboundary) {
        phi.FillBoundary(geom.periodicity());
    }

    if (m_coarsening_strategy == CoarseningStrategy::Sigma)
    {
        const auto lobc = LoBC();
        const auto hibc = HiBC();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(phi); mfi.isValid(); ++mfi)
        {
            Array4<Real> const& fab = phi.array(mfi);
            mlndlap_applybc(mfi.validbox(),fab,nd_domain,lobc,hibc);
        }
    }
}

#ifdef AMREX_USE_HYPRE
std::unique_ptr<HypreNodeLap>
MLNodeLinOp::makeHypreNodeLap (int bottom_verbose) const
{
    const BoxArray& ba = m_grids[0].back();
    const DistributionMapping& dm = m_dmap[0].back();
    const Geometry& geom = m_geom[0].back();
    const auto& factory = *(m_factory[0].back());
    const auto& owner_mask = *(m_owner_mask[0].back());
    const auto& dirichlet_mask = *(m_dirichlet_mask[0].back());
    MPI_Comm comm = BottomCommunicator();

    std::unique_ptr<HypreNodeLap> hypre_solver
        (new amrex::HypreNodeLap(ba, dm, geom, factory, owner_mask, dirichlet_mask,
                                 comm, this, bottom_verbose));

    return hypre_solver;
}
#endif

}

