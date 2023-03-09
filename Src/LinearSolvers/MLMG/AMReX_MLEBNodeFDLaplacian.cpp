#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_MLEBNodeFDLap_K.H>
#include <AMReX_MLNodeLap_K.H>
#include <AMReX_MLNodeTensorLap_K.H>
#include <AMReX_MultiFabUtil.H>

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
}

#ifdef AMREX_USE_EB
std::unique_ptr<FabFactory<FArrayBox> >
MLEBNodeFDLaplacian::makeFactory (int amrlev, int mglev) const
{
    return makeEBFabFactory(m_geom[amrlev][mglev],
                            m_grids[amrlev][mglev],
                            m_dmap[amrlev][mglev],
                            {1,1,1}, EBSupport::full);
}
#endif

void
MLEBNodeFDLaplacian::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
    BL_PROFILE("MLEBNodeFDLaplacian::restriction()");

    applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

    IntVect const ratio = mg_coarsen_ratio_vec[cmglev-1];
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

    IntVect const ratio = mg_coarsen_ratio_vec[fmglev];
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

#ifdef AMREX_USE_EB
    // Set covered nodes to Dirichlet, but with a negative value.
    // compGrad relies on the negative value to detect EB.
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev) {
            const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
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

    if (m_is_bottom_singular)
    {
        int amrlev = 0;
        int mglev = 0;
        auto const& dotmasks = m_coarse_dot_mask.arrays();
        auto const& dirmasks = m_dirichlet_mask[amrlev][mglev]->const_arrays();
        amrex::ParallelFor(m_coarse_dot_mask,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            if (dirmasks[box_no](i,j,k)) {
                dotmasks[box_no](i,j,k) = Real(0.);
            }
        });
    }

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
}

#ifdef AMREX_USE_EB
void
MLEBNodeFDLaplacian::scaleRHS (int amrlev, MultiFab& rhs) const
{
    auto const& dmask = *m_dirichlet_mask[amrlev][0];
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
    auto const& edgecent = factory->getEdgeCent();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        Array4<Real> const& rhsarr = rhs.array(mfi);
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
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(alpha == 0._rt, "alpha != 0 not implemented yet");
#endif
    AMREX_D_TERM(const Real bx = m_sigma[0]*dxinv[0]*dxinv[0];,
                 const Real by = m_sigma[1]*dxinv[1]*dxinv[1];,
                 const Real bz = m_sigma[2]*dxinv[2]*dxinv[2];)

    auto const& dmask = *m_dirichlet_mask[amrlev][mglev];

#ifdef AMREX_USE_EB
    const auto phieb = (m_in_solution_mode) ? m_s_phi_eb : Real(0.0);
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    auto const& edgecent = factory->getEdgeCent();
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
        bool cutfab = edgecent[0]->ok(mfi);
        if (cutfab) {
            AMREX_D_TERM(Array4<Real const> const& ecx = edgecent[0]->const_array(mfi);,
                         Array4<Real const> const& ecy = edgecent[1]->const_array(mfi);,
                         Array4<Real const> const& ecz = edgecent[2]->const_array(mfi));
            if (phieb == std::numeric_limits<Real>::lowest()) {
                auto const& phiebarr = m_phi_eb[amrlev].const_array(mfi);
#if (AMREX_SPACEDIM == 2)
                if (m_rz) {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_adotx_rz_eb(i,j,k,yarr,xarr,dmarr,ecx,ecy,
                                                phiebarr, sig0, dx0, dx1, xlo);
                    });
                } else
#endif
                {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_adotx_eb(i,j,k,yarr,xarr,dmarr,AMREX_D_DECL(ecx,ecy,ecz),
                                             phiebarr, AMREX_D_DECL(bx,by,bz));
                    });
                }
            } else {
#if (AMREX_SPACEDIM == 2)
                if (m_rz) {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_adotx_rz_eb(i,j,k,yarr,xarr,dmarr,ecx,ecy,
                                                phieb, sig0, dx0, dx1, xlo);
                    });
                } else
#endif
                {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_adotx_eb(i,j,k,yarr,xarr,dmarr,AMREX_D_DECL(ecx,ecy,ecz),
                                             phieb, AMREX_D_DECL(bx,by,bz));
                    });
                }
            }
        } else
#endif // AMREX_USE_EB
        {
#if (AMREX_SPACEDIM == 2)
            if (m_rz) {
                AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                {
                    mlebndfdlap_adotx_rz(i,j,k,yarr,xarr,dmarr,sig0,dx0,dx1,xlo);
                });
            } else
#endif
            {
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
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(alpha == 0._rt, "alpha != 0 not implemented yet");
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
        auto const& edgecent = factory->getEdgeCent();
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
            bool cutfab = edgecent[0]->ok(mfi);
            if (cutfab) {
                AMREX_D_TERM(Array4<Real const> const& ecx = edgecent[0]->const_array(mfi);,
                             Array4<Real const> const& ecy = edgecent[1]->const_array(mfi);,
                             Array4<Real const> const& ecz = edgecent[2]->const_array(mfi));
#if (AMREX_SPACEDIM == 2)
                if (m_rz) {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_gsrb_rz_eb(i,j,k,solarr,rhsarr,dmskarr,ecx,ecy,
                                               sig0, dx0, dx1, xlo, redblack);
                    });
                } else
#endif
                {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_gsrb_eb(i,j,k,solarr,rhsarr,dmskarr,AMREX_D_DECL(ecx,ecy,ecz),
                                            AMREX_D_DECL(bx,by,bz), redblack);
                    });
                }
            } else
#endif // AMREX_USE_EB
            {
#if (AMREX_SPACEDIM == 2)
                if (m_rz) {
                    AMREX_HOST_DEVICE_FOR_3D(box, i, j, k,
                    {
                        mlebndfdlap_gsrb_rz(i,j,k,solarr,rhsarr,dmskarr,
                                            sig0, dx0, dx1, xlo, redblack);
                    });
                } else
#endif
                {
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
    AMREX_ASSERT(factory);
    auto const& edgecent = factory->getEdgeCent();
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
#else
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
#endif
    }
}

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
void
MLEBNodeFDLaplacian::fillIJMatrix (MFIter const& /*mfi*/,
                                   Array4<HypreNodeLap::AtomicInt const> const& /*gid*/,
                                   Array4<int const> const& /*lid*/,
                                   HypreNodeLap::Int* const /*ncols*/,
                                   HypreNodeLap::Int* const /*cols*/,
                                   Real* const /*mat*/) const
{
    amrex::Abort("MLEBNodeFDLaplacian::fillIJMatrix: todo");
}

void
MLEBNodeFDLaplacian::fillRHS (MFIter const& /*mfi*/, Array4<int const> const& /*lid*/,
                              Real* const /*rhs*/, Array4<Real const> const& /*bfab*/) const
{
    amrex::Abort("MLEBNodeFDLaplacian::fillRHS: todo");
}
#endif

void
MLEBNodeFDLaplacian::postSolve (Vector<MultiFab>& sol) const
{
#ifdef AMREX_USE_EB
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        const auto phieb = m_s_phi_eb;
        const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        auto const& levset_mf = factory->getLevelSet();
        auto const& levset_ar = levset_mf.const_arrays();
        MultiFab& mf = sol[amrlev];
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

}
