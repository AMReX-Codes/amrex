
#include <AMReX_MLNodeLinOp.H>
#include <AMReX_MLNodeLap_F.H>
#include <AMReX_MLNodeLap_K.H>

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
    MLLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);

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
        }
    }


    m_cc_fine_mask.resize(m_num_amr_levels);
    m_nd_fine_mask.resize(m_num_amr_levels);
    m_has_fine_bndry.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        if (amrlev < m_num_amr_levels-1)
        {
            m_cc_fine_mask[amrlev].reset(new iMultiFab(m_grids[amrlev][0], m_dmap[amrlev][0], 1, 1));
            m_nd_fine_mask[amrlev].reset(new iMultiFab(amrex::convert(m_grids[amrlev][0],IntVect::TheNodeVector()),
                                                       m_dmap[amrlev][0], 1, 0));
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
                               const MultiFab* crse_bcdata)
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
                                 BCMode bc_mode, const MultiFab* crse_bcdata)
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
        ParallelAllReduce::Sum(result, Communicator(amrlev, mglev));
    }
    return result;
}

void
MLNodeLinOp::applyInhomogNeumannTerm (int amrlev, MultiFab& rhs) const
{
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
        m_is_bottom_singular = m_domain_covered[0];
    }

    const auto lobc = m_lobc[0];
    const auto hibc = m_hibc[0];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        std::vector< std::pair<int,Box> > isects;
        IArrayBox ccfab;

        for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
        {
            for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
            {
                const Geometry& geom = m_geom[amrlev][mglev];
                const auto& period = geom.periodicity();
                const Box& ccdomain = geom.Domain();
                const Box& nddomain = amrex::surroundingNodes(ccdomain);
                const std::vector<IntVect>& pshifts = period.shiftIntVect();

                Box ccdomain_p = ccdomain;
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    if (geom.isPeriodic(idim)) {
                        ccdomain_p.grow(idim, 1);
                    }
                }

                {
                    auto& dmask = *m_dirichlet_mask[amrlev][mglev];
                    const BoxArray& ccba = m_grids[amrlev][mglev];

                    MFItInfo mfi_info;
                    if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
                    for (MFIter mfi(dmask, mfi_info); mfi.isValid(); ++mfi)
                    {
                        const Box& ndbx = mfi.validbox();
                        const Box& ccbx = amrex::enclosedCells(ndbx);
                        const Box& ccbxg1 = amrex::grow(ccbx,1);
                        Array4<int> const& mskarr = dmask.array(mfi);

                        ccfab.resize(ccbxg1);
                        Elixir cceli = ccfab.elixir();
                        Array4<int> const& ccarr = ccfab.array();

                        if (ccdomain_p.contains(ccarr)) {
                            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(ccbxg1, i, j, k,
                            {
                                ccarr(i,j,k) = 1;
                            });
                        } else {
                            AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                            {
                                if (ccdomain_p.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                                    ccarr(i,j,k) = 1;
                                } else {
                                    ccarr(i,j,k) = 2;
                                }
                            });
                        }

                        for (const auto& iv : pshifts)
                        {
                            ccba.intersections(ccbxg1+iv, isects);
                            for (const auto& is : isects)
                            {
                                Box b = is.second-iv;
                                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(b, i, j, k,
                                {
                                    ccarr(i,j,k) = 0;
                                });
                            }
                        }

                        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( ndbx, tbx,
                        {
                            mlndlap_set_dirichlet_mask(tbx, mskarr, ccarr, nddomain,
                                                       lobc, hibc);
                        });
                    }
                }
            }
        }
    }

    for (int amrlev = 0; amrlev < m_num_amr_levels-1; ++amrlev)
    {
        iMultiFab& cc_mask = *m_cc_fine_mask[amrlev];
        iMultiFab& nd_mask = *m_nd_fine_mask[amrlev];
        LayoutData<int>& has_cf = *m_has_fine_bndry[amrlev];
        const BoxArray& fba = m_grids[amrlev+1][0];
        const BoxArray& cfba = amrex::coarsen(fba, AMRRefRatio(amrlev));

        const Box& ccdom = m_geom[amrlev][0].Domain();

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMRRefRatio(amrlev) == 2, "ref_ratio != 0 not supported");

        cc_mask.setVal(0);  // coarse by default

        const std::vector<IntVect>& pshifts = m_geom[amrlev][0].periodicity().shiftIntVect();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            std::vector< std::pair<int,Box> > isects;

            for (MFIter mfi(cc_mask); mfi.isValid(); ++mfi)
            {
                has_cf[mfi] = 0;
                const Box& bx = mfi.fabbox();
                Array4<int> const& fab = cc_mask.array(mfi);
                for (const auto& iv : pshifts)
                {
                    cfba.intersections(bx+iv, isects);
                    for (const auto& is : isects)
                    {
                        Box const& b = is.second-iv;
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(b,i,j,k,
                        {
                            fab(i,j,k) = 1;
                        });
                    }
                    if (!isects.empty()) has_cf[mfi] = 1;
                }

                mlndlap_fillbc_cc<int>(mfi.validbox(),fab,ccdom,m_lobc[0],m_hibc[0]);
            }
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
        const iMultiFab& omask = *m_owner_mask[amrlev][mglev];
        m_bottom_dot_mask.define(omask.boxArray(), omask.DistributionMap(), 1, 0);

        const Geometry& geom = m_geom[amrlev][mglev];
        Box nddomain = amrex::surroundingNodes(geom.Domain());

        if (m_coarsening_strategy != CoarseningStrategy::Sigma) {
            nddomain.grow(1000); // hack to avoid masks being modified at Neuman boundary
        }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(m_bottom_dot_mask,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& dfab = m_bottom_dot_mask.array(mfi);
            Array4<int const> const& sfab = omask.const_array(mfi);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlndlap_set_dot_mask(tbx, dfab, sfab, nddomain, lobc, hibc);
            });
        }
    }

    if (m_is_bottom_singular)
    {
        int amrlev = 0;
        int mglev = 0;
        const iMultiFab& omask = *m_owner_mask[amrlev][mglev];
        m_coarse_dot_mask.define(omask.boxArray(), omask.DistributionMap(), 1, 0);

        const Geometry& geom = m_geom[amrlev][mglev];
        Box nddomain = amrex::surroundingNodes(geom.Domain());

        if (m_coarsening_strategy != CoarseningStrategy::Sigma) {
            nddomain.grow(1000); // hack to avoid masks being modified at Neuman boundary
        }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(m_coarse_dot_mask,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& dfab = m_coarse_dot_mask.array(mfi);
            Array4<int const> const& sfab = omask.const_array(mfi);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlndlap_set_dot_mask(tbx, dfab, sfab, nddomain, lobc, hibc);
            });
        }
    }
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
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(phi); mfi.isValid(); ++mfi)
        {
            Array4<Real> const& fab = phi.array(mfi);
            mlndlap_applybc(mfi.validbox(),fab,nd_domain,m_lobc[0],m_hibc[0]);
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

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(NMGLevels(0) == 1,
                                     "MLNodeLaplacian: To use hypre, max_coarsening_level must be 0");

    std::unique_ptr<HypreNodeLap> hypre_solver
        (new amrex::HypreNodeLap(ba, dm, geom, factory, owner_mask, dirichlet_mask,
                                 comm, this, bottom_verbose));

    return hypre_solver;
}
#endif

}

