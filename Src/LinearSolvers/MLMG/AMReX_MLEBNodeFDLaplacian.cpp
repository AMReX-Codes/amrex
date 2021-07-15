#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_MLEBNodeFDLap_K.H>
#include <AMReX_MLNodeLap_K.H>
#include <AMReX_MLNodeTensorLap_K.H>

namespace amrex {

MLEBNodeFDLaplacian::MLEBNodeFDLaplacian (
    const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const LPInfo& a_info,
    const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

MLEBNodeFDLaplacian::~MLEBNodeFDLaplacian ()
{}

void
MLEBNodeFDLaplacian::setSigma (Array<Real,AMREX_SPACEDIM> const& a_sigma) noexcept
{
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        m_sigma[i] = a_sigma[i];
    }
}

void
MLEBNodeFDLaplacian::setEBDirichlet (Real a_phi_eb)
{
    m_s_phi_eb = a_phi_eb;
}

void
MLEBNodeFDLaplacian::define (const Vector<Geometry>& a_geom,
                             const Vector<BoxArray>& a_grids,
                             const Vector<DistributionMapping>& a_dmap,
                             const LPInfo& a_info,
                             const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    static_assert(AMREX_SPACEDIM > 1, "MLEBNodeFDLaplacian: 1D not supported");

    BL_PROFILE("MLEBNodeFDLaplacian::define()");

    // This makes sure grids are cell-centered;
    Vector<BoxArray> cc_grids = a_grids;
    for (auto& ba : cc_grids) {
        ba.enclosedCells();
    }

    if (a_grids.size() > 1) {
        amrex::Abort("MLEBNodeFDLaplacian: multi-level not supported");
    }

    Vector<FabFactory<FArrayBox> const*> _factory;
    for (auto x : a_factory) {
        _factory.push_back(static_cast<FabFactory<FArrayBox> const*>(x));
    }

    int eb_limit_coarsening = false;
    MLNodeLinOp::define(a_geom, cc_grids, a_dmap, a_info, _factory, eb_limit_coarsening);
}

std::unique_ptr<FabFactory<FArrayBox> >
MLEBNodeFDLaplacian::makeFactory (int amrlev, int mglev) const
{
    if (amrlev == 0 && mglev > 0) {
        return std::make_unique<FArrayBoxFactory>();
    } else {
        return makeEBFabFactory(m_geom[amrlev][mglev],
                                m_grids[amrlev][mglev],
                                m_dmap[amrlev][mglev],
                                {1,1,1}, EBSupport::full);
    }
}

void
MLEBNodeFDLaplacian::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
    BL_PROFILE("MLEBNodeFDLaplacian::restriction()");

    applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
        cfine.define(ba, fine.DistributionMap(), 1, 0);
    }

    MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][cmglev-1];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*pcrse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> cfab = pcrse->array(mfi);
        Array4<Real const> const& ffab = fine.const_array(mfi);
        Array4<int const> const& mfab = dmsk.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
        {
            mlndlap_restriction(i,j,k,cfab,ffab,mfab);
        });
    }

    if (need_parallel_copy) {
        crse.ParallelCopy(cfine);
    }
}

void
MLEBNodeFDLaplacian::interpolation (int amrlev, int fmglev, MultiFab& fine,
                                    const MultiFab& crse) const
{
    BL_PROFILE("MLEBNodeFDLaplacian::interpolation()");

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

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const& ffab = fine.array(mfi);
        Array4<Real const> const& cfab = cmf->const_array(mfi);
        Array4<int const> const& mfab = dmsk.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
        {
            mlndtslap_interpadd(i,j,k,ffab,cfab,mfab);
        });
    }
}

void
MLEBNodeFDLaplacian::averageDownSolutionRHS (int /*camrlev*/, MultiFab& /*crse_sol*/,
                                             MultiFab& /*crse_rhs*/,
                                             const MultiFab& /*fine_sol*/,
                                             const MultiFab& /*fine_rhs*/)
{
    amrex::Abort("MLEBNodeFDLaplacian::averageDownSolutionRHS: todo");
}

void
MLEBNodeFDLaplacian::reflux (int /*crse_amrlev*/, MultiFab& /*res*/,
                             const MultiFab& /*crse_sol*/, const MultiFab& /*crse_rhs*/,
                             MultiFab& /*fine_res*/, MultiFab& /*fine_sol*/,
                             const MultiFab& /*fine_rhs*/) const
{
    amrex::Abort("MLEBNodeFDLaplacian::reflux: TODO");
}

void
MLEBNodeFDLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLEBNodeFDLaplacian::prepareForSolve()");

    MLNodeLinOp::prepareForSolve();

    buildMasks();

    // Set covered nodes to Dirichlet
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev) {
            auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
            if (factory) {
                auto const& levset = factory->getLevelSet();
                auto& dmask = *m_dirichlet_mask[amrlev][mglev];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(dmask,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const Box& ndbx = mfi.tilebox();
                    Array4<int> const& mskarr = dmask.array(mfi);
                    Array4<Real const> const lstarr = levset.const_array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(ndbx, i, j, k,
                    {
                        if (lstarr(i,j,k) >= Real(0.0)) {
                            mskarr(i,j,k) = -1;
                        }
                    });
                }
            }
        }
    }

    m_acoef.clear();
    m_acoef.emplace_back(amrex::convert(m_grids[0][0],IntVect(1)),
                         m_dmap[0][0], 1, 0);
    const auto dxinv = m_geom[0][0].InvCellSizeArray();
    AMREX_D_TERM(const Real bcx = m_sigma[0]*dxinv[0]*dxinv[0];,
                 const Real bcy = m_sigma[1]*dxinv[1]*dxinv[1];,
                 const Real bcz = m_sigma[2]*dxinv[2]*dxinv[2];)
    m_a_huge = 1.e10 * AMREX_D_TERM(std::abs(bcx),+std::abs(bcy),+std::abs(bcz));
    Real ahuge = m_a_huge;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(m_acoef[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        auto const& acf = m_acoef[0].array(mfi);
        auto const& msk = m_dirichlet_mask[0][0]->const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
        {
            acf(i,j,k) = msk(i,j,k) ? ahuge : 0.0;
        });
    }

    for (int mglev = 1; mglev < m_num_mg_levels[0]; ++mglev) {
        m_acoef.emplace_back(amrex::convert(m_grids[0][mglev],IntVect(1)),
                             m_dmap[0][mglev], 1, 0);
        auto const& fine = m_acoef[mglev-1];
        auto      & crse = m_acoef[mglev];

        bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
        MultiFab cfine;
        if (need_parallel_copy) {
            const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
            cfine.define(ba, fine.DistributionMap(), 1, 0);
        }

        MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*pcrse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> cfab = pcrse->array(mfi);
            Array4<Real const> const& ffab = fine.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                cfab(i,j,k) = ffab(2*i,2*j,2*k);
            });
        }

        if (need_parallel_copy) {
            crse.ParallelCopy(cfine);
        }
    }
}

void
MLEBNodeFDLaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLEBNodeFDLaplacian::Fapply()");

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    AMREX_D_TERM(const Real bx = m_sigma[0]*dxinv[0]*dxinv[0];,
                 const Real by = m_sigma[1]*dxinv[1]*dxinv[1];,
                 const Real bz = m_sigma[2]*dxinv[2]*dxinv[2];)
    const auto phieb = (m_in_solution_mode) ? m_s_phi_eb : Real(0.0);

    auto const& dmask = *m_dirichlet_mask[amrlev][mglev];

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    if (factory) {
        auto const& edgecent = factory->getEdgeCent();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(out,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();
            Array4<Real const> const& xarr = in.const_array(mfi);
            Array4<Real> const& yarr = out.array(mfi);
            Array4<int const> const& dmarr = dmask.const_array(mfi);
            bool cutfab = edgecent[0]->ok(mfi);
            AMREX_D_TERM(Array4<Real const> const& ecx
                             = cutfab ? edgecent[0]->const_array(mfi) : Array4<Real const>{};,
                         Array4<Real const> const& ecy
                             = cutfab ? edgecent[1]->const_array(mfi) : Array4<Real const>{};,
                         Array4<Real const> const& ecz
                             = cutfab ? edgecent[2]->const_array(mfi) : Array4<Real const>{};)
            if (phieb == std::numeric_limits<Real>::lowest()) {
                auto const& phiebarr = m_phi_eb[amrlev].const_array(mfi);
                AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                {
                    mlebndfdlap_adotx_eb(i,j,k,yarr,xarr,dmarr,AMREX_D_DECL(ecx,ecy,ecz),
                                         phiebarr, AMREX_D_DECL(bx,by,bz));
                });
            } else {
                AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                {
                    mlebndfdlap_adotx_eb(i,j,k,yarr,xarr,dmarr,AMREX_D_DECL(ecx,ecy,ecz),
                                         phieb, AMREX_D_DECL(bx,by,bz));
                });
            }
        }
    } else {
        AMREX_ALWAYS_ASSERT(amrlev == 0);
        auto const& acoef = m_acoef[mglev];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(out,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();
            Array4<Real const> const& xarr = in.const_array(mfi);
            Array4<Real> const& yarr = out.array(mfi);
            Array4<int const> const& dmarr = dmask.const_array(mfi);
            Array4<Real const> const& acarr = acoef.const_array(mfi);
            AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
            {
                mlebndfdlap_adotx(i,j,k,yarr,xarr,dmarr,acarr,AMREX_D_DECL(bx,by,bz));
            });
        }
    }
}

void
MLEBNodeFDLaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const
{
    BL_PROFILE("MLEBNodeFDLaplacian::Fsmooth()");

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    AMREX_D_TERM(const Real bx = m_sigma[0]*dxinv[0]*dxinv[0];,
                 const Real by = m_sigma[1]*dxinv[1]*dxinv[1];,
                 const Real bz = m_sigma[2]*dxinv[2]*dxinv[2];)

    auto const& dmask = *m_dirichlet_mask[amrlev][mglev];

    for (int redblack = 0; redblack < 4; ++redblack) {
        if (redblack > 0) {
            applyBC(amrlev, mglev, sol, BCMode::Homogeneous, StateMode::Correction);
        }

        auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
        if (factory) {
            auto const& edgecent = factory->getEdgeCent();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(sol,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& box = mfi.tilebox();
                Array4<Real> const& solarr = sol.array(mfi);
                Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                Array4<int const> const& dmskarr = dmask.const_array(mfi);
                bool cutfab = edgecent[0]->ok(mfi);
                AMREX_D_TERM(Array4<Real const> const& ecx
                                 = cutfab ? edgecent[0]->const_array(mfi) : Array4<Real const>{};,
                             Array4<Real const> const& ecy
                                 = cutfab ? edgecent[1]->const_array(mfi) : Array4<Real const>{};,
                             Array4<Real const> const& ecz
                                 = cutfab ? edgecent[2]->const_array(mfi) : Array4<Real const>{};)

                AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                {
                    mlebndfdlap_gsrb_eb(i,j,k,solarr,rhsarr,dmskarr,AMREX_D_DECL(ecx,ecy,ecz),
                                        AMREX_D_DECL(bx,by,bz), redblack);
                });
            }
        } else {
            AMREX_ALWAYS_ASSERT(amrlev == 0);
            auto const& acoef = m_acoef[mglev];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(sol,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& box = mfi.tilebox();
                Array4<Real> const& solarr = sol.array(mfi);
                Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                Array4<int const> const& dmskarr = dmask.const_array(mfi);
                Array4<Real const> const& acarr = acoef.const_array(mfi);
                AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                {
                    mlebndfdlap_gsrb(i,j,k,solarr,rhsarr,dmskarr,acarr,
                                     AMREX_D_DECL(bx,by,bz), redblack);
                });
            }
        }
    }

    nodalSync(amrlev, mglev, sol);
}

void
MLEBNodeFDLaplacian::normalize (int amrlev, int mglev, MultiFab& mf) const
{
    if (amrlev == 0 && mglev > 0) {
        Real ahugeinv = Real(1.0) / m_a_huge;
        auto const& acoef = m_acoef[mglev];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();
            Array4<Real> const& fab = mf.array(mfi);
            Array4<Real const> const& acarr = acoef.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(box, i, j, k,
            {
                if (acarr(i,j,k) > Real(0.0)) {
                    fab(i,j,k) *= ahugeinv;
                }
            });
        }
    }
}

void
MLEBNodeFDLaplacian::fixUpResidualMask (int /*amrlev*/, iMultiFab& /*resmsk*/)
{
    amrex::Abort("MLEBNodeFDLaplacian::fixUpResidualMask: TODO");
}

void
MLEBNodeFDLaplacian::compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                               MultiFab& sol, Location /*loc*/) const
{
    BL_PROFILE("MLEBNodeFDLaplacian::compGrad()");

    AMREX_ASSERT(AMREX_D_TERM(grad[0]->ixType() == IndexType(IntVect(AMREX_D_DECL(0,1,1))),
                           && grad[1]->ixType() == IndexType(IntVect(AMREX_D_DECL(1,0,1))),
                           && grad[2]->ixType() == IndexType(IntVect(AMREX_D_DECL(1,1,0)))));
    const int mglev = 0;
    AMREX_D_TERM(const auto dxi = m_geom[amrlev][mglev].InvCellSize(0);,
                 const auto dyi = m_geom[amrlev][mglev].InvCellSize(1);,
                 const auto dzi = m_geom[amrlev][mglev].InvCellSize(2);)
    const auto phieb = m_s_phi_eb;

    auto const& dmask = *m_dirichlet_mask[amrlev][mglev];

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    AMREX_ASSERT(factory);
    auto const& edgecent = factory->getEdgeCent();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*grad[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        AMREX_D_TERM(const Box& xbox = mfi.tilebox(IntVect(AMREX_D_DECL(0,1,1)));,
                     const Box& ybox = mfi.tilebox(IntVect(AMREX_D_DECL(1,0,1)));,
                     const Box& zbox = mfi.tilebox(IntVect(AMREX_D_DECL(1,1,0)));)
        Array4<Real const> const& p = sol.const_array(mfi);
        AMREX_D_TERM(Array4<Real> const& gpx = grad[0]->array(mfi);,
                     Array4<Real> const& gpy = grad[1]->array(mfi);,
                     Array4<Real> const& gpz = grad[2]->array(mfi);)
        Array4<int const> const& dmarr = dmask.const_array(mfi);
        bool cutfab = edgecent[0]->ok(mfi);
        AMREX_D_TERM(Array4<Real const> const& ecx
                         = cutfab ? edgecent[0]->const_array(mfi) : Array4<Real const>{};,
                     Array4<Real const> const& ecy
                         = cutfab ? edgecent[1]->const_array(mfi) : Array4<Real const>{};,
                     Array4<Real const> const& ecz
                         = cutfab ? edgecent[2]->const_array(mfi) : Array4<Real const>{};)
        if (phieb == std::numeric_limits<Real>::lowest()) {
            auto const& phiebarr = m_phi_eb[amrlev].const_array(mfi);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM(
                xbox, txbox,
                {
                    mlebndfdlap_grad_x(txbox, gpx, p, dmarr, ecx, phiebarr, dxi);
                }
                , ybox, tybox,
                {
                    mlebndfdlap_grad_y(tybox, gpy, p, dmarr, ecy, phiebarr, dyi);
                }
                , zbox, tzbox,
                {
                    mlebndfdlap_grad_z(tzbox, gpz, p, dmarr, ecz, phiebarr, dzi);
                });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM(
                xbox, txbox,
                {
                    mlebndfdlap_grad_x(txbox, gpx, p, dmarr, ecx, phieb, dxi);
                }
                , ybox, tybox,
                {
                    mlebndfdlap_grad_y(tybox, gpy, p, dmarr, ecy, phieb, dyi);
                }
                , zbox, tzbox,
                {
                    mlebndfdlap_grad_z(tzbox, gpz, p, dmarr, ecz, phieb, dzi);
                });
        }
    }
}

#if defined(AMREX_USE_HYPRE)
void
MLEBNodeFDLaplacian::fillIJMatrix (MFIter const& mfi,
                                   Array4<HypreNodeLap::AtomicInt const> const& gid,
                                   Array4<int const> const& lid,
                                   HypreNodeLap::Int* const ncols,
                                   HypreNodeLap::Int* const cols,
                                   Real* const mat) const
{
    amrex::Abort("MLEBNodeFDLaplacian::fillIJMatrix: todo");
}

void
MLEBNodeFDLaplacian::fillRHS (MFIter const& mfi, Array4<int const> const& lid,
                              Real* const rhs, Array4<Real const> const& bfab) const
{
    amrex::Abort("MLEBNodeFDLaplacian::fillRHS: todo");
}
#endif

}
