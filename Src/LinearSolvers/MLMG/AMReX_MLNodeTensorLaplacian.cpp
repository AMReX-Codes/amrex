#include <AMReX_MLNodeTensorLaplacian.H>
#include <AMReX_MLNodeLinOp_K.H>
#include <AMReX_MLNodeTensorLap_K.H>
#include <AMReX_MultiFabUtil.H>

namespace amrex {

template <typename MF>
MLNodeTensorLaplacianT<MF>::MLNodeTensorLaplacianT (const Vector<Geometry>& a_geom,
                                              const Vector<BoxArray>& a_grids,
                                              const Vector<DistributionMapping>& a_dmap,
                                              const LPInfo& a_info)
{
    define(a_geom, a_grids, a_dmap, a_info);
}

template <typename MF>
void
MLNodeTensorLaplacianT<MF>::setSigma (Array<Real,nelems> const& a_sigma) noexcept
{
    for (int i = 0; i < nelems; ++i) { m_sigma[i] = a_sigma[i]; }
}

template <typename MF>
void
MLNodeTensorLaplacianT<MF>::setBeta (Array<Real,AMREX_SPACEDIM> const& a_beta) noexcept // NOLINT(readability-convert-member-functions-to-static)
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(a_beta);
#elif (AMREX_SPACEDIM == 2)
    m_sigma[0] = Real(1.) - a_beta[0]*a_beta[0];
    m_sigma[1] =          - a_beta[0]*a_beta[1];
    m_sigma[2] = Real(1.) - a_beta[1]*a_beta[1];
#elif (AMREX_SPACEDIM == 3)
    m_sigma[0] = Real(1.) - a_beta[0]*a_beta[0];
    m_sigma[1] =          - a_beta[0]*a_beta[1];
    m_sigma[2] =          - a_beta[0]*a_beta[2];
    m_sigma[3] = Real(1.) - a_beta[1]*a_beta[1];
    m_sigma[4] =          - a_beta[1]*a_beta[2];
    m_sigma[5] = Real(1.) - a_beta[2]*a_beta[2];
#endif
}

template <typename MF>
GpuArray<Real,MLNodeTensorLaplacianT<MF>::nelems>
MLNodeTensorLaplacianT<MF>::scaledSigma (int amrlev, int mglev) const noexcept
{
    auto s = m_sigma;
    auto const& dxinv = this->m_geom[amrlev][mglev].InvCellSizeArray();
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(dxinv);
#elif (AMREX_SPACEDIM == 2)
    s[0] *= dxinv[0]*dxinv[0];
    s[1] *= dxinv[0]*dxinv[1];
    s[2] *= dxinv[1]*dxinv[1];
#elif (AMREX_SPACEDIM == 3)
    s[0] *= dxinv[0]*dxinv[0];
    s[1] *= dxinv[0]*dxinv[1];
    s[2] *= dxinv[0]*dxinv[2];
    s[3] *= dxinv[1]*dxinv[1];
    s[4] *= dxinv[1]*dxinv[2];
    s[5] *= dxinv[2]*dxinv[2];
#endif
    return s;
}

template <typename MF>
void
MLNodeTensorLaplacianT<MF>::define (const Vector<Geometry>& a_geom,
                               const Vector<BoxArray>& a_grids,
                               const Vector<DistributionMapping>& a_dmap,
                               const LPInfo& a_info)
{
    BL_PROFILE("MLNodeTensorLaplacian::define()");

    // This makes sure grids are cell-centered;
    Vector<BoxArray> cc_grids = a_grids;
    for (auto& ba : cc_grids) {
        ba.enclosedCells();
    }

    this->m_coarsening_strategy = CoarseningStrategy::Sigma; // This will fill nodes outside Neumann BC
    MLNodeLinOpT<MF>::define(a_geom, cc_grids, a_dmap, a_info);
}

template <typename MF>
void
MLNodeTensorLaplacianT<MF>::restriction (int amrlev, int cmglev, MF& crse, MF& fine) const
{
    BL_PROFILE("MLNodeTensorLaplacian::restriction()");

    this->applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

    IntVect const ratio = (amrlev > 0) ? IntVect(2) : this->mg_coarsen_ratio_vec[cmglev-1];
    int semicoarsening_dir = this->info.semicoarsening_direction;

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MF cfine;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), ratio);
        cfine.define(ba, fine.DistributionMap(), 1, 0);
    }

    MF* pcrse = (need_parallel_copy) ? &cfine : &crse;
    const iMultiFab& dmsk = *this->m_dirichlet_mask[amrlev][cmglev-1];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*pcrse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<RT> cfab = pcrse->array(mfi);
        Array4<RT const> const& ffab = fine.const_array(mfi);
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

template <typename MF>
void
MLNodeTensorLaplacianT<MF>::interpolation (int amrlev, int fmglev, MF& fine,
                                      const MF& crse) const
{
    BL_PROFILE("MLNodeTensorLaplacian::interpolation()");

    IntVect const ratio = (amrlev > 0) ? IntVect(2) : this->mg_coarsen_ratio_vec[fmglev];
    int semicoarsening_dir = this->info.semicoarsening_direction;

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MF cfine;
    const MF* cmf = &crse;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), ratio);
        cfine.define(ba, fine.DistributionMap(), 1, 0);
        cfine.ParallelCopy(crse);
        cmf = &cfine;
    }

    const iMultiFab& dmsk = *this->m_dirichlet_mask[amrlev][fmglev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<RT> const& ffab = fine.array(mfi);
        Array4<RT const> const& cfab = cmf->const_array(mfi);
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

template <typename MF>
void
MLNodeTensorLaplacianT<MF>::averageDownSolutionRHS (int camrlev, MF& crse_sol, MF& /*crse_rhs*/,
                                               const MF& fine_sol, const MF& /*fine_rhs*/)
{
    const auto& amrrr = this->AMRRefRatio(camrlev);
    amrex::average_down(fine_sol, crse_sol, 0, 1, amrrr);

    if (this->isSingular(0))
    {
        amrex::Abort("MLNodeTensorLaplacian::averageDownSolutionRHS: TODO");
    }
}

template <typename MF>
void
MLNodeTensorLaplacianT<MF>::reflux (int /*crse_amrlev*/,
                               MF& /*res*/, const MF& /*crse_sol*/, const MF& /*crse_rhs*/,
                               MF& /*fine_res*/, MF& /*fine_sol*/, const MF& /*fine_rhs*/) const
{
    amrex::Abort("MLNodeTensorLaplacian::reflux: TODO");
}

template <typename MF>
void
MLNodeTensorLaplacianT<MF>::prepareForSolve ()
{
    BL_PROFILE("MLNodeTensorLaplacian::prepareForSolve()");

    MLNodeLinOpT<MF>::prepareForSolve();

    this->buildMasks();
}

template <typename MF>
void
MLNodeTensorLaplacianT<MF>::Fapply (int amrlev, int mglev, MF& out, const MF& in) const
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(amrlev, mglev, out, in);
#else
    BL_PROFILE("MLNodeTensorLaplacian::Fapply()");

    auto const& s = scaledSigma(amrlev, mglev);

    auto const& in_a = in.const_arrays();
    auto const& out_a = out.arrays();
    auto const& dmsk_a = this->m_dirichlet_mask[amrlev][mglev]->const_arrays();

    amrex::ParallelFor(out,
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
    {
        mlndtslap_adotx(i,j,k, out_a[box_no], in_a[box_no], dmsk_a[box_no], s);
    });
    if (!Gpu::inNoSyncRegion()) {
        Gpu::streamSynchronize();
    }
#endif
}

template <typename MF>
void
MLNodeTensorLaplacianT<MF>::smooth (int amrlev, int mglev, MF& sol, const MF& rhs,
                               bool skip_fillboundary, int niter) const
{
    BL_PROFILE("MLNodeTensorLaplacian::smooth()");
    for (int i = 0; i < niter; ++i) {
        for (int redblack = 0; redblack < 4; ++redblack) {
            if (!skip_fillboundary) {
                this->applyBC(amrlev, mglev, sol, BCMode::Homogeneous, StateMode::Correction);
            }
            m_redblack = redblack;
            Fsmooth(amrlev, mglev, sol, rhs);
            skip_fillboundary = false;
        }
        this->nodalSync(amrlev, mglev, sol);
    }
}

template <typename MF>
void
MLNodeTensorLaplacianT<MF>::Fsmooth (int amrlev, int mglev, MF& sol, const MF& rhs) const
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(amrlev, mglev, sol, rhs);
#else
    BL_PROFILE("MLNodeTensorLaplacian::Fsmooth()");

    auto const& s = scaledSigma(amrlev, mglev);

    auto const& sol_a = sol.arrays();
    auto const& rhs_a = rhs.const_arrays();
    auto const& dmsk_a = this->m_dirichlet_mask[amrlev][mglev]->const_arrays();
    int redblack = m_redblack;

    amrex::ParallelFor(sol,
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
    {
        if ((i+j+k+redblack) % 2 == 0) {
            mlndtslap_gauss_seidel(i, j, k, sol_a[box_no], rhs_a[box_no], dmsk_a[box_no], s);
        }
    });
    if (!Gpu::inNoSyncRegion()) {
        Gpu::streamSynchronize();
    }
#endif
}

template <typename MF>
void
MLNodeTensorLaplacianT<MF>::normalize (int amrlev, int mglev, MF& mf) const
{
    amrex::ignore_unused(amrlev,mglev,mf);
}

template <typename MF>
void
MLNodeTensorLaplacianT<MF>::fixUpResidualMask (int /*amrlev*/, iMultiFab& /*resmsk*/)
{
    amrex::Abort("MLNodeTensorLaplacian::fixUpResidualMask: TODO");
}

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
template <typename MF>
void
MLNodeTensorLaplacianT<MF>::fillIJMatrix (MFIter const& mfi,
                                     Array4<HypreNodeLap::AtomicInt const> const& gid,
                                     Array4<int const> const& lid,
                                     HypreNodeLap::Int* ncols,
                                     HypreNodeLap::Int* cols,
                                     Real* mat) const
{
    const int amrlev = 0;
    const int mglev = this->NMGLevels(amrlev)-1;
    auto const& s = scaledSigma(amrlev, mglev);

    const Box& ndbx = mfi.validbox();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE
            (static_cast<Long>(ndbx.numPts())*AMREX_D_TERM(3,*3,*3) <
             static_cast<Long>(std::numeric_limits<int>::max()),
             "The Box is too big.  We could use Long here, but it would much slower.");
        const int nmax = ndbx.numPts() * AMREX_D_TERM(3,*3,*3);
        const auto ndlo = amrex::lbound(ndbx);
        const auto ndlen = amrex::length(ndbx);
        amrex::Scan::PrefixSum<int>
            (nmax,
             [=] AMREX_GPU_DEVICE (int offset) noexcept
             {
                 Dim3 node = nodelap_detail::GetNode()(ndlo, ndlen, offset);
                 Dim3 node2 = nodelap_detail::GetNode2()(offset, node);
                 return (lid(node.x,node.y,node.z) >= 0 &&
                         gid(node2.x,node2.y,node2.z)
                         < std::numeric_limits<HypreNodeLap::AtomicInt>::max());
             },
             [=] AMREX_GPU_DEVICE (int offset, int ps) noexcept
             {
                 Dim3 node = nodelap_detail::GetNode()(ndlo, ndlen, offset);
                 mlndtslap_fill_ijmatrix_gpu(ps, node.x, node.y, node.z, offset,
                                             ndbx, gid, lid, ncols, cols, mat, s);
             },
             amrex::Scan::Type::exclusive);
    } else
#endif
    {
        mlndtslap_fill_ijmatrix_cpu(ndbx, gid, lid, ncols, cols, mat, s);
    }
}

template <typename MF>
void
MLNodeTensorLaplacianT<MF>::fillRHS (MFIter const& mfi, Array4<int const> const& lid,
                                Real* rhs, Array4<Real const> const& bfab) const
{
    const Box& bx = mfi.validbox();
    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
    {
        if (lid(i,j,k) >= 0) {
            rhs[lid(i,j,k)] = bfab(i,j,k);
        }
    });
}
#endif

template class MLNodeTensorLaplacianT<MultiFab>;
template class MLNodeTensorLaplacianT<fMultiFab>;

}
