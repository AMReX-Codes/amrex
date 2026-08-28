#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_MLEBNodeFDLap_K.H>
#include <AMReX_MLNodeLinOp_K.H>
#include <AMReX_MLNodeTensorLap_K.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_SingleBoxCGSolver.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBMultiFabUtil.H>
#endif

namespace amrex {

#ifdef AMREX_USE_EB
MLEBNodeFDLaplacian::MLEBNodeFDLaplacian (
    const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const LPInfo& a_info,
    const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}
#endif

MLEBNodeFDLaplacian::MLEBNodeFDLaplacian (
    const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const LPInfo& a_info)
{
    define(a_geom, a_grids, a_dmap, a_info);
}

void
MLEBNodeFDLaplacian::setSigma (Array<Real,AMREX_SPACEDIM> const& a_sigma) noexcept
{
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        m_sigma[i] = a_sigma[i];
    }
}

void
MLEBNodeFDLaplacian::setSigma (int amrlev, MultiFab const& a_sigma)
{
    m_needs_update = true;
    m_has_sigma_mf = true;
    m_sigma_mf[amrlev][0] = std::make_unique<MultiFab>
        (this->m_grids[amrlev][0], this->m_dmap[amrlev][0], 1, 1, MFInfo{},
         *(this->m_factory[amrlev][0]));
    MultiFab::Copy(*m_sigma_mf[amrlev][0], a_sigma, 0, 0, 1, 0);
#ifdef AMREX_USE_EB
    amrex::EB_set_covered(*m_sigma_mf[amrlev][0], Real(0.0));
#endif
}

void
MLEBNodeFDLaplacian::setRZ (bool flag) // NOLINT
{
#if (AMREX_SPACEDIM == 2)
    m_rz = flag;
#else
    amrex::ignore_unused(flag, m_rz);
#endif
}

void
MLEBNodeFDLaplacian::setAlpha (Real a_alpha) // NOLINT
{
#if (AMREX_SPACEDIM == 2)
    m_rz_alpha = a_alpha;
#else
    amrex::ignore_unused(a_alpha);
#endif
}

#ifdef AMREX_USE_EB

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
    for (const auto *x : a_factory) {
        _factory.push_back(static_cast<FabFactory<FArrayBox> const*>(x));
    }

    int eb_limit_coarsening = true;
    m_coarsening_strategy = CoarseningStrategy::Sigma; // This will fill nodes outside Neumann BC
    MLNodeLinOp::define(a_geom, cc_grids, a_dmap, a_info, _factory, eb_limit_coarsening);

    m_sigma_mf.resize(this->m_num_amr_levels);
    for (int ilev = 0; ilev < this->m_num_amr_levels; ++ilev) {
        m_sigma_mf[ilev].resize(this->m_num_mg_levels[ilev]);
    }
}

#endif

void
MLEBNodeFDLaplacian::define (const Vector<Geometry>& a_geom,
                             const Vector<BoxArray>& a_grids,
                             const Vector<DistributionMapping>& a_dmap,
                             const LPInfo& a_info)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMREX_SPACEDIM>1, "MLEBNodeFDLaplacian: 1D not supported");

    BL_PROFILE("MLEBNodeFDLaplacian::define()");

    // This makes sure grids are cell-centered;
    Vector<BoxArray> cc_grids = a_grids;
    for (auto& ba : cc_grids) {
        ba.enclosedCells();
    }

    if (a_grids.size() > 1) {
        amrex::Abort("MLEBNodeFDLaplacian: multi-level not supported");
    }

    m_coarsening_strategy = CoarseningStrategy::Sigma; // This will fill nodes outside Neumann BC
    MLNodeLinOp::define(a_geom, cc_grids, a_dmap, a_info);

    m_sigma_mf.resize(this->m_num_amr_levels);
    for (int ilev = 0; ilev < this->m_num_amr_levels; ++ilev) {
        m_sigma_mf[ilev].resize(this->m_num_mg_levels[ilev]);
    }
}

#ifdef AMREX_USE_EB
std::unique_ptr<FabFactory<FArrayBox> >
MLEBNodeFDLaplacian::makeFactory (int amrlev, int mglev) const
{
    if (EB2::TopIndexSpaceIfPresent()) {
        auto const* fact0 = Factory(amrlev,0) ? Factory(amrlev, 0) : Factory(0,0);
        return makeEBFabFactory(static_cast<EBFArrayBoxFactory const*>(fact0)->getEBIndexSpace(),
                                m_geom[amrlev][mglev],
                                m_grids[amrlev][mglev],
                                m_dmap[amrlev][mglev],
                                {1,1,1}, EBSupport::full);
    } else {
        return MLNodeLinOp::makeFactory(amrlev, mglev);
    }
}
#endif

void
MLEBNodeFDLaplacian::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
    BL_PROFILE("MLEBNodeFDLaplacian::restriction()");

    applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

    IntVect const ratio = (amrlev > 0) ? IntVect(2) : mg_coarsen_ratio_vec[cmglev-1];
    int semicoarsening_dir = info.semicoarsening_direction;

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), ratio);
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
        if (ratio == 2) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_restriction(i,j,k,cfab,ffab,mfab);
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_semi_restriction(i,j,k,cfab,ffab,mfab, semicoarsening_dir);
            });
        }
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

    IntVect const ratio = (amrlev > 0) ? IntVect(2) : mg_coarsen_ratio_vec[fmglev];
    int semicoarsening_dir = info.semicoarsening_direction;

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    const MultiFab* cmf = &crse;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), ratio);
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
        if (ratio == 2) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndtslap_interpadd(i,j,k,ffab,cfab,mfab);
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndtslap_semi_interpadd(i,j,k,ffab,cfab,mfab,semicoarsening_dir);
            });
        }
    }
}

void
MLEBNodeFDLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLEBNodeFDLaplacian::prepareForSolve()");

    MLNodeLinOp::prepareForSolve();

    buildMasks();

#ifdef AMREX_USE_EB
    // Set covered nodes to Dirichlet, but with a negative value.
    // compGrad relies on the negative value to detect EB.
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev) {
            const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
            if (factory && !factory->isAllRegular()) {
                auto const& levset_mf = factory->getLevelSet();
                auto const& levset_ar = levset_mf.const_arrays();
                auto& dmask_mf = *m_dirichlet_mask[amrlev][mglev];
                auto const& dmask_ar = dmask_mf.arrays();
                amrex::ParallelFor(dmask_mf,
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                {
                    if (levset_ar[box_no](i,j,k) >= Real(0.0)) {
                        dmask_ar[box_no](i,j,k) = -1;
                    }
                });
            }
        }
    }
#endif

    {
        int amrlev = 0;
        int mglev = m_num_mg_levels[amrlev]-1;
        auto const& dotmasks = m_bottom_dot_mask.arrays();
        auto const& dirmasks = m_dirichlet_mask[amrlev][mglev]->const_arrays();
        amrex::ParallelFor(m_bottom_dot_mask,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            if (dirmasks[box_no](i,j,k)) {
                dotmasks[box_no](i,j,k) = Real(0.);
            }
        });
    }

    AMREX_ASSERT(!isBottomSingular());

    Gpu::streamSynchronize();

#if (AMREX_SPACEDIM == 2)
    if (m_rz) {
        if (m_geom[0][0].ProbLo(0) == 0._rt) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_lobc[0][0] == BCType::Neumann,
                                             "The lo-x BC must be Neumann for 2d RZ");
        }
        if (m_sigma[0] == 0._rt) {
            m_sigma[0] = 1._rt; // For backward compatibility
        }
    }
#endif

    if (m_has_sigma_mf) {
        update_sigma();
    }
}

#ifdef AMREX_USE_EB
bool
MLEBNodeFDLaplacian::scaleRHS (int amrlev, MultiFab* rhs) const
{
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());

    if (!factory) {return false; }

    if (rhs) {
        auto const& dmask = *m_dirichlet_mask[amrlev][0];
        auto const& edgecent = factory->getEdgeCent();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();
            Array4<Real> const& rhsarr = rhs->array(mfi);
            Array4<int const> const& dmarr = dmask.const_array(mfi);
            bool cutfab = edgecent[0]->ok(mfi);
            if (cutfab) {
                AMREX_D_TERM(Array4<Real const> const& ecx = edgecent[0]->const_array(mfi);,
                             Array4<Real const> const& ecy = edgecent[1]->const_array(mfi);,
                             Array4<Real const> const& ecz = edgecent[2]->const_array(mfi));
                AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                {
                    mlebndfdlap_scale_rhs(i,j,k,rhsarr,dmarr,AMREX_D_DECL(ecx,ecy,ecz));
                });
            }
        }
    }

    return true;
}
#endif

void
MLEBNodeFDLaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLEBNodeFDLaplacian::Fapply()");

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
#if (AMREX_SPACEDIM == 2)
    const auto sig0 = m_sigma[0];
    const auto dx0 = m_geom[amrlev][mglev].CellSize(0);
    const auto dx1 = m_geom[amrlev][mglev].CellSize(1)/std::sqrt(m_sigma[1]);
    const auto xlo = m_geom[amrlev][mglev].ProbLo(0);
    const auto alpha = m_rz_alpha;
#endif
    AMREX_D_TERM(const Real bx = m_sigma[0]*dxinv[0]*dxinv[0];,
                 const Real by = m_sigma[1]*dxinv[1]*dxinv[1];,
                 const Real bz = m_sigma[2]*dxinv[2]*dxinv[2];)

    auto const& dmask = *m_dirichlet_mask[amrlev][mglev];

#ifdef AMREX_USE_EB
    const auto phieb = (m_in_solution_mode && !this->m_precond_mode) ? m_s_phi_eb : Real(0.0);
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    Array<const MultiCutFab*,AMREX_SPACEDIM> edgecent {AMREX_D_DECL(nullptr,nullptr,nullptr)};
    if (factory) {
        edgecent = factory->getEdgeCent();
    }
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        Array4<Real const> const& xarr = in.const_array(mfi);
        Array4<Real> const& yarr = out.array(mfi);
        Array4<int const> const& dmarr = dmask.const_array(mfi);
#ifdef AMREX_USE_EB
        bool cutfab = edgecent[0] && edgecent[0]->ok(mfi);
        if (cutfab && factory) { // clang-tidy is not that smart
            AMREX_D_TERM(Array4<Real const> const& ecx = edgecent[0]->const_array(mfi);,
                         Array4<Real const> const& ecy = edgecent[1]->const_array(mfi);,
                         Array4<Real const> const& ecz = edgecent[2]->const_array(mfi));
            auto const& levset = factory->getLevelSet().const_array(mfi);
            if (phieb == std::numeric_limits<Real>::lowest()) {
                auto const& phiebarr = m_phi_eb[amrlev].const_array(mfi);
#if (AMREX_SPACEDIM == 2)
                if (m_rz) {
                    if (m_has_sigma_mf) {
                        auto const& sigarr = m_sigma_mf[amrlev][mglev]->const_array(mfi);
                        auto const& vfrc = factory->getVolFrac().const_array(mfi);
                        AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                        {
                            mlebndfdlap_sig_adotx_rz_eb(i,j,k,yarr,xarr,levset,dmarr,ecx,ecy,
                                                        sigarr, vfrc, phiebarr, dx0, dx1, xlo, alpha);
                        });
                    } else {
                        AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                        {
                            mlebndfdlap_adotx_rz_eb(i,j,k,yarr,xarr,levset,dmarr,ecx,ecy,
                                                    phiebarr, sig0, dx0, dx1, xlo, alpha);
                        });
                    }
                } else
#endif
                if (m_has_sigma_mf) {
                    auto const& sigarr = m_sigma_mf[amrlev][mglev]->const_array(mfi);
                    auto const& vfrc = factory->getVolFrac().const_array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_sig_adotx_eb(i,j,k,yarr,xarr,levset,dmarr,AMREX_D_DECL(ecx,ecy,ecz),
                                                 sigarr, vfrc, phiebarr, AMREX_D_DECL(bx,by,bz));
                    });
                } else {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_adotx_eb(i,j,k,yarr,xarr,levset,dmarr,AMREX_D_DECL(ecx,ecy,ecz),
                                             phiebarr, AMREX_D_DECL(bx,by,bz));
                    });
                }
            } else {
#if (AMREX_SPACEDIM == 2)
                if (m_rz) {
                    if (m_has_sigma_mf) {
                        auto const& sigarr = m_sigma_mf[amrlev][mglev]->const_array(mfi);
                        auto const& vfrc = factory->getVolFrac().const_array(mfi);
                        AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                        {
                            mlebndfdlap_sig_adotx_rz_eb(i,j,k,yarr,xarr,levset,dmarr,ecx,ecy,
                                                        sigarr, vfrc, phieb, dx0, dx1, xlo, alpha);
                        });
                    } else {
                        AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                        {
                            mlebndfdlap_adotx_rz_eb(i,j,k,yarr,xarr,levset,dmarr,ecx,ecy,
                                                    phieb, sig0, dx0, dx1, xlo, alpha);
                        });
                    }
                } else
#endif
                if (m_has_sigma_mf) {
                    auto const& sigarr = m_sigma_mf[amrlev][mglev]->const_array(mfi);
                    auto const& vfrc = factory->getVolFrac().const_array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_sig_adotx_eb(i,j,k,yarr,xarr,levset,dmarr,AMREX_D_DECL(ecx,ecy,ecz),
                                                 sigarr, vfrc, phieb, AMREX_D_DECL(bx,by,bz));
                    });
                } else {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_adotx_eb(i,j,k,yarr,xarr,levset,dmarr,AMREX_D_DECL(ecx,ecy,ecz),
                                             phieb, AMREX_D_DECL(bx,by,bz));
                    });
                }
            }
        } else
#endif // AMREX_USE_EB
        {
#if (AMREX_SPACEDIM == 2)
            if (m_rz) {
                if (m_has_sigma_mf) {
                    auto const& sigarr = m_sigma_mf[amrlev][mglev]->const_array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_sig_adotx_rz(i,j,k,yarr,xarr,dmarr,sigarr,dx0,dx1,xlo,alpha);
                    });
                } else {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_adotx_rz(i,j,k,yarr,xarr,dmarr,sig0,dx0,dx1,xlo,alpha);
                    });
                }
            } else
#endif
            if (m_has_sigma_mf) {
                auto const& sigarr = m_sigma_mf[amrlev][mglev]->const_array(mfi);
                AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                {
                    mlebndfdlap_sig_adotx(i,j,k,yarr,xarr,dmarr,sigarr,AMREX_D_DECL(bx,by,bz));
                });
            } else {
                AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                {
                    mlebndfdlap_adotx(i,j,k,yarr,xarr,dmarr,AMREX_D_DECL(bx,by,bz));
                });
            }
        }
    }
}

void
MLEBNodeFDLaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const
{
    BL_PROFILE("MLEBNodeFDLaplacian::Fsmooth()");

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
#if (AMREX_SPACEDIM == 2)
    const auto sig0 = m_sigma[0];
    const auto dx0 = m_geom[amrlev][mglev].CellSize(0);
    const auto dx1 = m_geom[amrlev][mglev].CellSize(1)/std::sqrt(m_sigma[1]);
    const auto xlo = m_geom[amrlev][mglev].ProbLo(0);
    const auto alpha = m_rz_alpha;
#endif
    AMREX_D_TERM(const Real bx = m_sigma[0]*dxinv[0]*dxinv[0];,
                 const Real by = m_sigma[1]*dxinv[1]*dxinv[1];,
                 const Real bz = m_sigma[2]*dxinv[2]*dxinv[2];)

    auto const& dmask = *m_dirichlet_mask[amrlev][mglev];

    for (int redblack = 0; redblack < 2; ++redblack) {
        if (redblack > 0) {
            applyBC(amrlev, mglev, sol, BCMode::Homogeneous, StateMode::Correction);
        }

#ifdef AMREX_USE_EB
        const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
        Array<const MultiCutFab*,AMREX_SPACEDIM> edgecent {AMREX_D_DECL(nullptr,nullptr,nullptr)};
        if (factory) {
            edgecent = factory->getEdgeCent();
        }
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(sol,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();
            Array4<Real> const& solarr = sol.array(mfi);
            Array4<Real const> const& rhsarr = rhs.const_array(mfi);
            Array4<int const> const& dmskarr = dmask.const_array(mfi);
#ifdef AMREX_USE_EB
            bool cutfab = edgecent[0] && edgecent[0]->ok(mfi);
            if (cutfab && factory) { // clang-tidy is not that smart
                AMREX_D_TERM(Array4<Real const> const& ecx = edgecent[0]->const_array(mfi);,
                             Array4<Real const> const& ecy = edgecent[1]->const_array(mfi);,
                             Array4<Real const> const& ecz = edgecent[2]->const_array(mfi));
                auto const& levset = factory->getLevelSet().const_array(mfi);
#if (AMREX_SPACEDIM == 2)
                if (m_rz) {
                    if (m_has_sigma_mf) {
                        auto const& sigarr = m_sigma_mf[amrlev][mglev]->const_array(mfi);
                        auto const& vfrc = factory->getVolFrac().const_array(mfi);
                        AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                        {
                            mlebndfdlap_sig_gsrb_rz_eb(i,j,k,solarr,rhsarr,levset,dmskarr,ecx,ecy,
                                                       sigarr, vfrc, dx0, dx1, xlo, redblack, alpha);
                        });
                    } else {
                        AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                        {
                            mlebndfdlap_gsrb_rz_eb(i,j,k,solarr,rhsarr,levset,dmskarr,ecx,ecy,
                                                   sig0, dx0, dx1, xlo, redblack, alpha);
                        });
                    }
                } else
#endif
                if (m_has_sigma_mf) {
                    auto const& sigarr = m_sigma_mf[amrlev][mglev]->const_array(mfi);
                    auto const& vfrc = factory->getVolFrac().const_array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_sig_gsrb_eb(i,j,k,solarr,rhsarr,levset,dmskarr,AMREX_D_DECL(ecx,ecy,ecz),
                                                sigarr, vfrc, AMREX_D_DECL(bx,by,bz), redblack);
                    });
                } else {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_gsrb_eb(i,j,k,solarr,rhsarr,levset,dmskarr,AMREX_D_DECL(ecx,ecy,ecz),
                                            AMREX_D_DECL(bx,by,bz), redblack);
                    });
                }
            } else
#endif // AMREX_USE_EB
            {
#if (AMREX_SPACEDIM == 2)
                if (m_rz) {
                    if (m_has_sigma_mf) {
                        auto const& sigarr = m_sigma_mf[amrlev][mglev]->const_array(mfi);
                        AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                        {
                            mlebndfdlap_sig_gsrb_rz(i,j,k,solarr,rhsarr,dmskarr,sigarr,
                                                    dx0, dx1, xlo, redblack, alpha);
                        });
                    } else {
                        AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                        {
                            mlebndfdlap_gsrb_rz(i,j,k,solarr,rhsarr,dmskarr,
                                                sig0, dx0, dx1, xlo, redblack, alpha);
                        });
                    }
                } else
#endif
                if (m_has_sigma_mf) {
                    auto const& sigarr = m_sigma_mf[amrlev][mglev]->const_array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_sig_gsrb(i,j,k,solarr,rhsarr,dmskarr,sigarr,
                                             AMREX_D_DECL(bx,by,bz), redblack);
                    });
                } else {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_gsrb(i,j,k,solarr,rhsarr,dmskarr,
                                         AMREX_D_DECL(bx,by,bz), redblack);
                    });
                }
            }
        }
    }

    nodalSync(amrlev, mglev, sol);
}

void
MLEBNodeFDLaplacian::normalize (int amrlev, int mglev, MultiFab& mf) const
{
    amrex::ignore_unused(amrlev, mglev, mf);
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
    if (amrex::isMFIterSafe(*grad[0], sol)) {
        this->compGrad_doit(amrlev, grad, sol);
    } else {
        Array<MultiFab,AMREX_SPACEDIM> grad_tmp;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            grad_tmp[idim].define(amrex::convert(sol.boxArray(), grad[idim]->ixType()),
                                  sol.DistributionMap(), 1, 0,
                                  MFInfo{}.SetArena(The_Async_Arena()));
        }
        this->compGrad_doit(amrlev, GetArrOfPtrs(grad_tmp), sol);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            grad[idim]->ParallelCopy(grad_tmp[idim], 0, 0, 1);
        }
    }
}

void
MLEBNodeFDLaplacian::compGrad_doit (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                                    MultiFab& sol) const
{
    AMREX_ASSERT(AMREX_D_TERM(grad[0]->ixType() == IndexType(IntVect(AMREX_D_DECL(0,1,1))),
                           && grad[1]->ixType() == IndexType(IntVect(AMREX_D_DECL(1,0,1))),
                           && grad[2]->ixType() == IndexType(IntVect(AMREX_D_DECL(1,1,0)))));

    const int mglev = 0;
    AMREX_D_TERM(const auto dxi = m_geom[amrlev][mglev].InvCellSize(0);,
                 const auto dyi = m_geom[amrlev][mglev].InvCellSize(1);,
                 const auto dzi = m_geom[amrlev][mglev].InvCellSize(2);)

#ifdef AMREX_USE_EB
    auto const& dmask = *m_dirichlet_mask[amrlev][mglev];
    const auto phieb = m_s_phi_eb;
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    Array<const MultiCutFab*,AMREX_SPACEDIM> edgecent {AMREX_D_DECL(nullptr,nullptr,nullptr)};
    if (factory) {
        edgecent = factory->getEdgeCent();
    }
#endif

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
#ifdef AMREX_USE_EB
        if (factory && !factory->isAllRegular()) {
            Array4<int const> const& dmarr = dmask.const_array(mfi);
            bool cutfab = edgecent[0] && edgecent[0]->ok(mfi);
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
        } else
#endif
        {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM(
                    xbox, txbox,
                    {
                        mlebndfdlap_grad_x(txbox, gpx, p, dxi);
                    }
                    , ybox, tybox,
                    {
                        mlebndfdlap_grad_y(tybox, gpy, p, dyi);
                    }
                    , zbox, tzbox,
                    {
                        mlebndfdlap_grad_z(tzbox, gpz, p, dzi);
                    });
        }
    }
}

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
void
MLEBNodeFDLaplacian::fillIJMatrix (MFIter const& /*mfi*/,
                                   Array4<HypreNodeLap::AtomicInt const> const& /*gid*/,
                                   Array4<int const> const& /*lid*/,
                                   HypreNodeLap::Int* /*ncols*/,
                                   HypreNodeLap::Int* /*cols*/,
                                   Real* /*mat*/) const
{
    amrex::Abort("MLEBNodeFDLaplacian::fillIJMatrix: todo");
}

void
MLEBNodeFDLaplacian::fillRHS (MFIter const& /*mfi*/, Array4<int const> const& /*lid*/,
                              Real* /*rhs*/, Array4<Real const> const& /*bfab*/) const
{
    amrex::Abort("MLEBNodeFDLaplacian::fillRHS: todo");
}
#endif

void
MLEBNodeFDLaplacian::postSolve (Vector<MultiFab*> const& sol) const
{
#ifdef AMREX_USE_EB
    if (this->m_precond_mode) { return; }
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        const auto phieb = m_s_phi_eb;
        const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        if (!factory || factory->isAllRegular()) { return; }
        auto const& levset_mf = factory->getLevelSet();
        auto const& levset_ar = levset_mf.const_arrays();
        MultiFab& mf = *sol[amrlev];
        auto const& sol_ar = mf.arrays();
        if (phieb == std::numeric_limits<Real>::lowest()) {
            auto const& phieb_ar = m_phi_eb[amrlev].const_arrays();
            amrex::ParallelFor(mf, IntVect(1),
            [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
            {
                if (levset_ar[bi](i,j,k) >= Real(0.0)) {
                    sol_ar[bi](i,j,k) = phieb_ar[bi](i,j,k);
                }
            });
        } else {
            amrex::ParallelFor(mf, IntVect(1),
            [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k) noexcept
            {
                if (levset_ar[bi](i,j,k) >= Real(0.0)) {
                    sol_ar[bi](i,j,k) = phieb;
                }
            });
        }
    }
#else
    amrex::ignore_unused(sol);
#endif
}

void
MLEBNodeFDLaplacian::update ()
{
    if (MLNodeLinOp::needsUpdate()) {
        MLNodeLinOp::update();
    }

    if (m_needs_update && m_has_sigma_mf) {
        update_sigma();
    }
    m_needs_update = false;
}

void
MLEBNodeFDLaplacian::update_sigma ()
{
    AMREX_D_TERM(m_sigma[0] = Real(1.0);,
                 m_sigma[1] = Real(1.0);,
                 m_sigma[2] = Real(1.0));
    AMREX_ALWAYS_ASSERT(this->m_num_amr_levels == 1);
    for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev) {
        for (int mglev = 1; mglev < this->m_num_mg_levels[amrlev]; ++mglev) {
            if (m_sigma_mf[amrlev][mglev] == nullptr) {
                m_sigma_mf[amrlev][mglev] = std::make_unique<MultiFab>
                    (this->m_grids[amrlev][mglev], this->m_dmap[amrlev][mglev], 1, 1,
                     MFInfo{}, *(this->m_factory[amrlev][mglev]));
            }
            IntVect const ratio = (amrlev > 0) ? IntVect (2)
                : this->mg_coarsen_ratio_vec[mglev-1];
#ifdef AMREX_USE_EB
            amrex::EB_average_down
#else
            amrex::average_down
#endif
                (*m_sigma_mf[amrlev][mglev-1],
                 *m_sigma_mf[amrlev][mglev], 0, 1, ratio);
        }

        for (int mglev = 0; mglev < this->m_num_mg_levels[amrlev]; ++mglev) {
            auto const& geom = this->m_geom[amrlev][mglev];
            auto& sigma = *m_sigma_mf[amrlev][mglev];
            sigma.FillBoundary(geom.periodicity());

            const Box& domain = geom.Domain();
            const auto lobc = LoBC();
            const auto hibc = HiBC();

            MFItInfo mfi_info;
            if (Gpu::notInLaunchRegion()) { mfi_info.SetDynamic(true); }
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
}

namespace {
    struct LPBase
    {
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real xdoty (IntVect const& iv, int, Real vx, Real vy) const
        {
            return dotmsk(iv)*vx*vy;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void normalize (IntVect const&, int, Real&) const {}

        // The neighbor index returned by lowerNeighbor/upperNeighbor is for
        // the solution data only, which lives in a Box without ghost cells.
        // EB data (level set and edge centroids) does have ghost cells.
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        int lowerNeighbor (int i, int idim) const
        {
            // There is only a single Box. Thus a box boundary node is
            // either Dirichlet (including domain Dirichlet boundary or
            // coarse/fine Dirichlet boundary) or on the periodic or Neumann
            // domain boundary. Because this is called on non-Dirichlet
            // nodes only, if a boundary node (e.g., i == dlo[idim]) gets
            // here, it's either periodic or Neumann.
            return (i == dlo[idim])
                ? (is_periodic[idim] ? dhi[idim]-1 : i+1)
                : i-1;
        }

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        int upperNeighbor (int i, int idim) const
        {
            return (i == dhi[idim])
                ? (is_periodic[idim] ? dlo[idim]+1 : i-1)
                : i+1;
        }

        Array4<Real const> dotmsk;
        Array4<int const> dirmsk;
        GpuArray<Real,AMREX_SPACEDIM> beta;
        GpuArray<int,AMREX_SPACEDIM> dlo, dhi;
        GpuArray<bool,AMREX_SPACEDIM> is_periodic;
    };

    struct LP
        : public LPBase
    {
        LP (LPBase const& a_lpbase)
            : LPBase(a_lpbase)
            {}

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real apply (IntVect const& iv, int, Array4<Real> const& xa, int n) const
        {
            int const i = iv[0];
            int const j = iv[1];
#if (AMREX_SPACEDIM == 3)
            int const k = iv[2];
#else
            int const k = 0;
#endif
            if (dirmsk(i,j,k)) {
                return Real(0.0);
            }

            Real const xc = xa(i,j,k,n);
            int const im = lowerNeighbor(i,0);
            int const ip = upperNeighbor(i,0);
            int const jm = lowerNeighbor(j,1);
            int const jp = upperNeighbor(j,1);
            Real y = beta[0] * (xa(im,j,k,n) + xa(ip,j,k,n) - Real(2.0)*xc)
                +    beta[1] * (xa(i,jm,k,n) + xa(i,jp,k,n) - Real(2.0)*xc);
#if (AMREX_SPACEDIM == 3)
            int const km = lowerNeighbor(k,2);
            int const kp = upperNeighbor(k,2);
            y += beta[2] * (xa(i,j,km,n) + xa(i,j,kp,n) - Real(2.0)*xc);
#endif
            return y;
        }
    };

    struct LPSigma
        : public LPBase
    {
        LPSigma (LPBase const& a_lpbase, Array4<Real const> const& a_sigma)
            : LPBase(a_lpbase), sigma(a_sigma)
            {}

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real apply (IntVect const& iv, int, Array4<Real> const& xa, int n) const
        {
            int const i = iv[0];
            int const j = iv[1];
#if (AMREX_SPACEDIM == 3)
            int const k = iv[2];
#else
            int const k = 0;
#endif
            if (dirmsk(i,j,k)) {
                return Real(0.0);
            }

            Real const xc = xa(i,j,k,n);
            int const im = lowerNeighbor(i,0);
            int const ip = upperNeighbor(i,0);
            int const jm = lowerNeighbor(j,1);
            int const jp = upperNeighbor(j,1);
#if (AMREX_SPACEDIM == 2)
            Real const sigxm = Real(0.5)*(sigma(i-1,j-1,k) + sigma(i-1,j,k));
            Real const sigxp = Real(0.5)*(sigma(i  ,j-1,k) + sigma(i  ,j,k));
            Real y = beta[0] * (sigxm*(xa(im,j,k,n)-xc) + sigxp*(xa(ip,j,k,n)-xc));

            Real const sigym = Real(0.5)*(sigma(i-1,j-1,k) + sigma(i,j-1,k));
            Real const sigyp = Real(0.5)*(sigma(i-1,j  ,k) + sigma(i,j  ,k));
            y += beta[1] * (sigym*(xa(i,jm,k,n)-xc) + sigyp*(xa(i,jp,k,n)-xc));
#else
            Real const sigxm = Real(0.25)*(sigma(i-1,j-1,k-1) + sigma(i-1,j,k-1)
                                               + sigma(i-1,j-1,k) + sigma(i-1,j,k));
            Real const sigxp = Real(0.25)*(sigma(i,j-1,k-1) + sigma(i,j,k-1)
                                               + sigma(i,j-1,k) + sigma(i,j,k));
            Real y = beta[0] * (sigxm*(xa(im,j,k,n)-xc) + sigxp*(xa(ip,j,k,n)-xc));

            Real const sigym = Real(0.25)*(sigma(i-1,j-1,k-1) + sigma(i,j-1,k-1)
                                               + sigma(i-1,j-1,k) + sigma(i,j-1,k));
            Real const sigyp = Real(0.25)*(sigma(i-1,j,k-1) + sigma(i,j,k-1)
                                               + sigma(i-1,j,k) + sigma(i,j,k));
            y += beta[1] * (sigym*(xa(i,jm,k,n)-xc) + sigyp*(xa(i,jp,k,n)-xc));

            int const km = lowerNeighbor(k,2);
            int const kp = upperNeighbor(k,2);
            Real const sigzm = Real(0.25)*(sigma(i-1,j-1,k-1) + sigma(i,j-1,k-1)
                                               + sigma(i-1,j,k-1) + sigma(i,j,k-1));
            Real const sigzp = Real(0.25)*(sigma(i-1,j-1,k) + sigma(i,j-1,k)
                                               + sigma(i-1,j,k) + sigma(i,j,k));
            y += beta[2] * (sigzm*(xa(i,j,km,n)-xc) + sigzp*(xa(i,j,kp,n)-xc));
#endif
            return y;
        }

        Array4<Real const> sigma;
    };

#if (AMREX_SPACEDIM == 2)
    struct RZConstSigma
    {
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real radialEdgeSigma (int, int, int, Real) const
        {
            return sigr;
        }

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        static Real axialEdgeSigma (int, int, int, Real)
        {
            return Real(1.0);
        }

        Real sigr;
    };

    struct RZVarSigma
    {
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real radialEdgeSigma (int i, int j, int k, Real) const
        {
            return Real(0.5)*(sigma(i,j-1,k)+sigma(i,j,k));
        }

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real axialEdgeSigma (int i, int j, int k, Real r) const
        {
            Real const rp = r + Real(0.5)*dr;
            Real const rm = amrex::Math::abs(r-Real(0.5)*dr);
            return (sigma(i-1,j,k)*rm + sigma(i,j,k)*rp) / (rm+rp);
        }

        Array4<Real const> sigma;
        Real dr;
    };

#ifdef AMREX_USE_EB
    struct RZVarSigmaEB
    {
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real radialEdgeSigma (int i, int j, int k, Real) const
        {
            Real const vfm = vfrac(i,j-1,k);
            Real const vfp = vfrac(i,j,k);
            return (sigma(i,j-1,k)*vfm + sigma(i,j,k)*vfp) / (vfm+vfp);
        }

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real axialEdgeSigma (int i, int j, int k, Real r) const
        {
            Real const rp = r + Real(0.5)*dr;
            Real const rm = amrex::Math::abs(r-Real(0.5)*dr);
            Real const wm = vfrac(i-1,j,k)*rm;
            Real const wp = vfrac(i,j,k)*rp;
            return (sigma(i-1,j,k)*wm + sigma(i,j,k)*wp) / (wm+wp);
        }

        Array4<Real const> sigma;
        Array4<Real const> vfrac;
        Real dr;
    };
#endif

    template <bool UseEB>
    struct RZEBData {};

#ifdef AMREX_USE_EB
    template <>
    struct RZEBData<true>
    {
        Array4<Real const> levset;
        GpuArray<Array4<Real const>,AMREX_SPACEDIM> edgecent;
    };
#endif

    template <bool UseEB, typename SigmaPolicy>
    struct LPRZ
        : public LPBase, public SigmaPolicy, public RZEBData<UseEB>
    {
        LPRZ (LPBase const& a_lpbase, SigmaPolicy const& a_sigma_policy,
              RZEBData<UseEB> const& a_eb_data, Real a_dr, Real a_dz,
              Real a_rlo, Real a_alpha)
            : LPBase(a_lpbase), SigmaPolicy(a_sigma_policy), RZEBData<UseEB>(a_eb_data),
              dr(a_dr), dz(a_dz), rlo(a_rlo), alpha(a_alpha)
            {}

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real apply (IntVect const& iv, int, Array4<Real> const& xa, int n) const
        {
            int const i = iv[0];
            int const j = iv[1];
            int const k = 0;
            Real const r = rlo + Real(i)*dr;
            if (dirmsk(i,j,k) || (r == Real(0.0) && alpha != Real(0.0))) {
                return Real(0.0);
            }

            int const im = lowerNeighbor(i,0);
            int const ip = upperNeighbor(i,0);
            int const jm = lowerNeighbor(j,1);
            int const jp = upperNeighbor(j,1);
            Real const xc = xa(i,j,k,n);
            Real out;
            Real scale = Real(1.0);

            if (r == Real(0.0)) {
                Real const sigp = this->radialEdgeSigma(i,j,k,r);
                if constexpr (UseEB) {
                    if (this->levset(i+1,j,k) >= Real(0.0)) {
                        Real const hp = (this->edgecent[0](i,j,k) == Real(1.0))
                            ? Real(1.0) : Real(1.0)+Real(2.0)*this->edgecent[0](i,j,k);
                        out = -Real(4.0)*sigp*xc/(dr*dr*hp*hp);
                        scale = hp;
                    } else {
                        out = Real(4.0)*sigp*(xa(ip,j,k,n)-xc)/(dr*dr);
                    }
                } else {
                    out = Real(4.0)*sigp*(xa(ip,j,k,n)-xc)/(dr*dr);
                }
            } else {
                Real hp = Real(1.0);
                Real hm = Real(1.0);
                if constexpr (UseEB) {
                    hp = (this->edgecent[0](i,j,k) == Real(1.0))
                        ? Real(1.0) : Real(1.0)+Real(2.0)*this->edgecent[0](i,j,k);
                    hm = (this->edgecent[0](i-1,j,k) == Real(1.0))
                        ? Real(1.0) : Real(1.0)-Real(2.0)*this->edgecent[0](i-1,j,k);
                }

                Real const sigp = this->radialEdgeSigma(i,j,k,r);
                Real const sigm = this->radialEdgeSigma(i-1,j,k,r);
                Real tmp;
                if constexpr (UseEB) {
                    tmp = (this->levset(i+1,j,k) < Real(0.0))
                        ? sigp*(xa(ip,j,k,n)-xc)*(r+Real(0.5)*dr)
                        : -sigp*xc/hp*(r+Real(0.5)*hp*dr);
                    tmp += (this->levset(i-1,j,k) < Real(0.0))
                        ? sigm*(xa(im,j,k,n)-xc)*(r-Real(0.5)*dr)
                        : -sigm*xc/hm*(r-Real(0.5)*hm*dr);
                } else {
                    tmp = sigp*(xa(ip,j,k,n)-xc)*(r+Real(0.5)*dr)
                        + sigm*(xa(im,j,k,n)-xc)*(r-Real(0.5)*dr);
                }
                out = tmp*Real(2.0)/((hp+hm)*r*dr*dr);
                scale = amrex::min(hm,hp);
            }

            Real hp = Real(1.0);
            Real hm = Real(1.0);
            if constexpr (UseEB) {
                hp = (this->edgecent[1](i,j,k) == Real(1.0))
                    ? Real(1.0) : Real(1.0)+Real(2.0)*this->edgecent[1](i,j,k);
                hm = (this->edgecent[1](i,j-1,k) == Real(1.0))
                    ? Real(1.0) : Real(1.0)-Real(2.0)*this->edgecent[1](i,j-1,k);
            }

            Real const sigp = this->axialEdgeSigma(i,j,k,r);
            Real const sigm = this->axialEdgeSigma(i,j-1,k,r);
            Real tmp;
            if constexpr (UseEB) {
                tmp = (this->levset(i,j+1,k) < Real(0.0))
                    ? sigp*(xa(i,jp,k,n)-xc) : -sigp*xc/hp;
                tmp += (this->levset(i,j-1,k) < Real(0.0))
                    ? sigm*(xa(i,jm,k,n)-xc) : -sigm*xc/hm;
            } else {
                tmp = sigp*(xa(i,jp,k,n)-xc) + sigm*(xa(i,jm,k,n)-xc);
            }
            out += tmp*Real(2.0)/((hp+hm)*dz*dz);
            scale = amrex::min(scale,hm,hp);

            if (r != Real(0.0)) {
                out -= alpha*xc/(r*r);
            }
            return out*scale;
        }

        Real dr;
        Real dz;
        Real rlo;
        Real alpha;
    };
#endif

#ifdef AMREX_USE_EB
    struct EBConstSigma
    {
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        static Real edgeSigmaX (int, int, int)
        {
            return Real(1.0);
        }

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        static Real edgeSigmaY (int, int, int)
        {
            return Real(1.0);
        }

#if (AMREX_SPACEDIM == 3)
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        static Real edgeSigmaZ (int, int, int)
        {
            return Real(1.0);
        }
#endif
    };

    struct EBVarSigma
    {
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real edgeSigmaX (int i, int j, int k) const
        {
#if (AMREX_SPACEDIM == 2)
            Real const vf0 = vfrac(i,j-1,k);
            Real const vf1 = vfrac(i,j,k);
            return (sigma(i,j-1,k)*vf0 + sigma(i,j,k)*vf1) / (vf0+vf1);
#else
            Real const vf0 = vfrac(i,j-1,k-1);
            Real const vf1 = vfrac(i,j  ,k-1);
            Real const vf2 = vfrac(i,j-1,k  );
            Real const vf3 = vfrac(i,j  ,k  );
            return (sigma(i,j-1,k-1)*vf0 + sigma(i,j,k-1)*vf1
                    + sigma(i,j-1,k)*vf2 + sigma(i,j,k)*vf3)
                / (vf0+vf1+vf2+vf3);
#endif
        }

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real edgeSigmaY (int i, int j, int k) const
        {
#if (AMREX_SPACEDIM == 2)
            Real const vf0 = vfrac(i-1,j,k);
            Real const vf1 = vfrac(i,j,k);
            return (sigma(i-1,j,k)*vf0 + sigma(i,j,k)*vf1) / (vf0+vf1);
#else
            Real const vf0 = vfrac(i-1,j,k-1);
            Real const vf1 = vfrac(i  ,j,k-1);
            Real const vf2 = vfrac(i-1,j,k  );
            Real const vf3 = vfrac(i  ,j,k  );
            return (sigma(i-1,j,k-1)*vf0 + sigma(i,j,k-1)*vf1
                    + sigma(i-1,j,k)*vf2 + sigma(i,j,k)*vf3)
                / (vf0+vf1+vf2+vf3);
#endif
        }

#if (AMREX_SPACEDIM == 3)
        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real edgeSigmaZ (int i, int j, int k) const
        {
            Real const vf0 = vfrac(i-1,j-1,k);
            Real const vf1 = vfrac(i  ,j-1,k);
            Real const vf2 = vfrac(i-1,j  ,k);
            Real const vf3 = vfrac(i  ,j  ,k);
            return (sigma(i-1,j-1,k)*vf0 + sigma(i,j-1,k)*vf1
                    + sigma(i-1,j,k)*vf2 + sigma(i,j,k)*vf3)
                / (vf0+vf1+vf2+vf3);
        }
#endif

        Array4<Real const> sigma;
        Array4<Real const> vfrac;
    };

    template <typename SigmaPolicy>
    struct LPEB
        : public LPBase, public SigmaPolicy
    {
        LPEB (LPBase const& a_lpbase, Array4<Real const> const& a_levset,
              GpuArray<Array4<Real const>,AMREX_SPACEDIM> const& a_edgecent,
              SigmaPolicy const& a_sigma_policy)
            : LPBase(a_lpbase), SigmaPolicy(a_sigma_policy),
              levset(a_levset), edgecent(a_edgecent)
            {}

        [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        Real apply (IntVect const& iv, int, Array4<Real> const& xa, int n) const
        {
            int const i = iv[0];
            int const j = iv[1];
#if (AMREX_SPACEDIM == 3)
            int const k = iv[2];
#else
            int const k = 0;
#endif
            if (dirmsk(i,j,k)) {
                return Real(0.0);
            }

            Real const xc = xa(i,j,k,n);
            int const im = lowerNeighbor(i,0);
            int const ip = upperNeighbor(i,0);
            Real const hpx = (edgecent[0](i,j,k) == Real(1.0))
                ? Real(1.0) : Real(1.0)+Real(2.0)*edgecent[0](i,j,k);
            Real const hmx = (edgecent[0](i-1,j,k) == Real(1.0))
                ? Real(1.0) : Real(1.0)-Real(2.0)*edgecent[0](i-1,j,k);
            Real const sigxp = this->edgeSigmaX(i,j,k);
            Real const sigxm = this->edgeSigmaX(i-1,j,k);
            Real tmp = (levset(i+1,j,k) < Real(0.0))
                ? sigxp*(xa(ip,j,k,n)-xc) : -sigxp*xc/hpx;
            tmp += (levset(i-1,j,k) < Real(0.0))
                ? sigxm*(xa(im,j,k,n)-xc) : -sigxm*xc/hmx;
            Real y = beta[0]*tmp*Real(2.0)/(hpx+hmx);
            Real scale = amrex::min(hmx,hpx);

            int const jm = lowerNeighbor(j,1);
            int const jp = upperNeighbor(j,1);
            Real const hpy = (edgecent[1](i,j,k) == Real(1.0))
                ? Real(1.0) : Real(1.0)+Real(2.0)*edgecent[1](i,j,k);
            Real const hmy = (edgecent[1](i,j-1,k) == Real(1.0))
                ? Real(1.0) : Real(1.0)-Real(2.0)*edgecent[1](i,j-1,k);
            Real const sigyp = this->edgeSigmaY(i,j,k);
            Real const sigym = this->edgeSigmaY(i,j-1,k);
            tmp = (levset(i,j+1,k) < Real(0.0))
                ? sigyp*(xa(i,jp,k,n)-xc) : -sigyp*xc/hpy;
            tmp += (levset(i,j-1,k) < Real(0.0))
                ? sigym*(xa(i,jm,k,n)-xc) : -sigym*xc/hmy;
            y += beta[1]*tmp*Real(2.0)/(hpy+hmy);
            scale = amrex::min(scale,hmy,hpy);

#if (AMREX_SPACEDIM == 3)
            int const km = lowerNeighbor(k,2);
            int const kp = upperNeighbor(k,2);
            Real const hpz = (edgecent[2](i,j,k) == Real(1.0))
                ? Real(1.0) : Real(1.0)+Real(2.0)*edgecent[2](i,j,k);
            Real const hmz = (edgecent[2](i,j,k-1) == Real(1.0))
                ? Real(1.0) : Real(1.0)-Real(2.0)*edgecent[2](i,j,k-1);
            Real const sigzp = this->edgeSigmaZ(i,j,k);
            Real const sigzm = this->edgeSigmaZ(i,j,k-1);
            tmp = (levset(i,j,k+1) < Real(0.0))
                ? sigzp*(xa(i,j,kp,n)-xc) : -sigzp*xc/hpz;
            tmp += (levset(i,j,k-1) < Real(0.0))
                ? sigzm*(xa(i,j,km,n)-xc) : -sigzm*xc/hmz;
            y += beta[2]*tmp*Real(2.0)/(hpz+hmz);
            scale = amrex::min(scale,hmz,hpz);
#endif

            return y*scale;
        }

        Array4<Real const> levset;
        GpuArray<Array4<Real const>,AMREX_SPACEDIM> edgecent;
    };
#endif
}

void
MLEBNodeFDLaplacian::customBottomSolve (MLMGT<MultiFab>* mlmg, MultiFab& x, const MultiFab& b,
                                        Real eps_rel, Real eps_abs, int maxiter)
{
    amrex::ignore_unused(eps_rel, eps_abs, maxiter);

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    bool use_custom_solver = (x.size() == 1);
    if (use_custom_solver)
    {
        int const amrlev = 0;
        int const mglev = NMGLevels(0) - 1;

#ifdef AMREX_USE_EB
        auto const* eb_factory = dynamic_cast<EBFArrayBoxFactory const*>
            (m_factory[amrlev][mglev].get());
        bool const use_eb = eb_factory && !eb_factory->isAllRegular();
#endif

        int ret = 0;
        if (ParallelDescriptor::MyProc() == x.DistributionMap()[0])
        {
            auto const& geom = m_geom[amrlev][mglev];
            const auto dxinv = geom.InvCellSizeArray();
#if (AMREX_SPACEDIM == 2)
            const auto sig0 = m_sigma[0];
            const auto dx0 = geom.CellSize(0);
            const auto dx1 = geom.CellSize(1)/std::sqrt(m_sigma[1]);
            const auto xlo = geom.ProbLo(0);
            const auto alpha = m_rz_alpha;
#endif
            AMREX_D_TERM(const Real bx = m_sigma[0]*dxinv[0]*dxinv[0];,
                         const Real by = m_sigma[1]*dxinv[1]*dxinv[1];,
                         const Real bz = m_sigma[2]*dxinv[2]*dxinv[2];)

            auto const& dotmsk = m_bottom_dot_mask[0].const_array();
            auto const& dirmsk = (*m_dirichlet_mask[amrlev][mglev])[0].const_array();
            Box box = x.boxArray()[0];
            Box dbox = amrex::convert(geom.Domain(),IntVect(1));
            LPBase lpbase{dotmsk,dirmsk,
                          GpuArray<Real,AMREX_SPACEDIM>{AMREX_D_DECL(bx,by,bz)},
                          GpuArray<int,AMREX_SPACEDIM>
                              {AMREX_D_DECL(dbox.smallEnd(0),
                                            dbox.smallEnd(1),
                                            dbox.smallEnd(2))},
                          GpuArray<int,AMREX_SPACEDIM>
                              {AMREX_D_DECL(dbox.bigEnd(0),
                                            dbox.bigEnd(1),
                                            dbox.bigEnd(2))},
                          GpuArray<bool,AMREX_SPACEDIM>{AMREX_D_DECL(geom.isPeriodic(0),
                                                                     geom.isPeriodic(1),
                                                                     geom.isPeriodic(2))}};
#if (AMREX_SPACEDIM == 2)
            if (m_rz) {
#ifdef AMREX_USE_EB
                if (use_eb) {
                    auto const& edgecent = eb_factory->getEdgeCent();
                    AMREX_ALWAYS_ASSERT(edgecent[0] && edgecent[0]->ok(0));
                    GpuArray<Array4<Real const>,AMREX_SPACEDIM> const ec
                        {(*edgecent[0])[0].const_array(), (*edgecent[1])[0].const_array()};
                    auto const& levset = eb_factory->getLevelSet()[0].const_array();
                    if (m_has_sigma_mf) {
                        auto const& sigma = (*m_sigma_mf[amrlev][mglev])[0].const_array();
                        auto const& vfrac = eb_factory->getVolFrac()[0].const_array();
                        LPRZ<true,RZVarSigmaEB> lp
                            (lpbase, RZVarSigmaEB{sigma,vfrac,dx0},
                             RZEBData<true>{levset,ec}, dx0, dx1, xlo, alpha);
                        ret = bicgstab_solve(box, x[0], b[0], lp,
                                             eps_rel, eps_abs, maxiter);
                    } else {
                        LPRZ<true,RZConstSigma> lp
                            (lpbase, RZConstSigma{sig0}, RZEBData<true>{levset,ec},
                             dx0, dx1, xlo, alpha);
                        ret = bicgstab_solve(box, x[0], b[0], lp,
                                             eps_rel, eps_abs, maxiter);
                    }
                } else
#endif
                {
                    if (m_has_sigma_mf) {
                        auto const& sigma = (*m_sigma_mf[amrlev][mglev])[0].const_array();
                        LPRZ<false,RZVarSigma> lp
                            (lpbase, RZVarSigma{sigma,dx0}, RZEBData<false>{},
                             dx0, dx1, xlo, alpha);
                        ret = bicgstab_solve(box, x[0], b[0], lp,
                                             eps_rel, eps_abs, maxiter);
                    } else {
                        LPRZ<false,RZConstSigma> lp
                            (lpbase, RZConstSigma{sig0}, RZEBData<false>{},
                             dx0, dx1, xlo, alpha);
                        ret = bicgstab_solve(box, x[0], b[0], lp,
                                             eps_rel, eps_abs, maxiter);
                    }
                }
            } else
#endif
#ifdef AMREX_USE_EB
            if (use_eb) {
                auto const& edgecent = eb_factory->getEdgeCent();
                AMREX_ALWAYS_ASSERT(edgecent[0] && edgecent[0]->ok(0));
                GpuArray<Array4<Real const>,AMREX_SPACEDIM> const ec
                    {AMREX_D_DECL((*edgecent[0])[0].const_array(),
                                  (*edgecent[1])[0].const_array(),
                                  (*edgecent[2])[0].const_array())};
                auto const& levset = eb_factory->getLevelSet()[0].const_array();
                if (m_has_sigma_mf) {
                    auto const& sigma = (*m_sigma_mf[amrlev][mglev])[0].const_array();
                    auto const& vfrac = eb_factory->getVolFrac()[0].const_array();
                    LPEB<EBVarSigma> lp(lpbase, levset, ec, EBVarSigma{sigma,vfrac});
                    ret = bicgstab_solve(box, x[0], b[0], lp, eps_rel, eps_abs, maxiter);
                } else {
                    LPEB<EBConstSigma> lp(lpbase, levset, ec, EBConstSigma{});
                    ret = bicgstab_solve(box, x[0], b[0], lp, eps_rel, eps_abs, maxiter);
                }
            } else
#endif
            if (m_has_sigma_mf) {
                auto const& sigma = (*m_sigma_mf[amrlev][mglev])[0].const_array();
                LPSigma lp(lpbase, sigma);
                ret = bicgstab_solve(box, x[0], b[0], lp, eps_rel, eps_abs, maxiter);
            } else {
                LP lp(lpbase);
                ret = bicgstab_solve(box, x[0], b[0], lp, eps_rel, eps_abs, maxiter);
            }
        }

        if (ParallelContext::NProcsSub() > 1) {
            int root = ParallelContext::global_to_local_rank(x.DistributionMap()[0]);
            ParallelDescriptor::Bcast(&ret, 1, root, ParallelContext::CommunicatorSub());
        }

        mlmg->postCG(ret);
    } else
#endif
    {
        int ret = mlmg->bottomSolveWithCG(x, b, MLCGSolverT<MultiFab>::Type::BiCGStab);
        mlmg->postCG(ret);
    }
}

}
