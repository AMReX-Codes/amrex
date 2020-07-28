#include <AMReX_MLNodeTensorLaplacian.H>
#include <AMReX_MLNodeLap_K.H>
#include <AMReX_MLNodeTensorLap_K.H>
#include <AMReX_MultiFabUtil.H>

namespace amrex {

MLNodeTensorLaplacian::MLNodeTensorLaplacian (const Vector<Geometry>& a_geom,
                                              const Vector<BoxArray>& a_grids,
                                              const Vector<DistributionMapping>& a_dmap,
                                              const LPInfo& a_info)
{
    define(a_geom, a_grids, a_dmap, a_info);
}

MLNodeTensorLaplacian::~MLNodeTensorLaplacian ()
{}

void
MLNodeTensorLaplacian::setSigma (Array<Real,nelems> const& a_sigma) noexcept
{
    for (int i = 0; i < nelems; ++i) m_sigma[i] = a_sigma[i];
}

void
MLNodeTensorLaplacian::setBeta (Array<Real,AMREX_SPACEDIM> const& a_beta) noexcept
{
#if (AMREX_SPACEDIM == 2)
    m_sigma[0] = 1. - a_beta[0]*a_beta[0];
    m_sigma[1] =    - a_beta[0]*a_beta[1];
    m_sigma[2] = 1. - a_beta[1]*a_beta[1];
#elif (AMREX_SPACEDIM == 3)
    m_sigma[0] = 1. - a_beta[0]*a_beta[0];
    m_sigma[1] =    - a_beta[0]*a_beta[1];
    m_sigma[2] =    - a_beta[0]*a_beta[2];
    m_sigma[3] = 1. - a_beta[1]*a_beta[1];
    m_sigma[4] =    - a_beta[1]*a_beta[2];
    m_sigma[5] = 1. - a_beta[2]*a_beta[2];
#endif
}

void
MLNodeTensorLaplacian::define (const Vector<Geometry>& a_geom,
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

    MLNodeLinOp::define(a_geom, cc_grids, a_dmap, a_info);
}

void
MLNodeTensorLaplacian::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
    BL_PROFILE("MLNodeTensorLaplacian::restriction()");

    applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
        cfine.define(ba, fine.DistributionMap(), 1, 0);
    }

    MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][cmglev-1];

#ifdef _OPENMP
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
MLNodeTensorLaplacian::interpolation (int amrlev, int fmglev, MultiFab& fine,
                                      const MultiFab& crse) const
{
    BL_PROFILE("MLNodeTensorLaplacian::interpolation()");

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

#ifdef _OPENMP
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
MLNodeTensorLaplacian::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& /*crse_rhs*/,
                                               const MultiFab& fine_sol, const MultiFab& /*fine_rhs*/)
{
    const auto& amrrr = AMRRefRatio(camrlev);
    amrex::average_down(fine_sol, crse_sol, 0, 1, amrrr);

    if (isSingular(0))
    {
        amrex::Abort("MLNodeTensorLaplacian::averageDownSolutionRHS: TODO");
    }
}

void
MLNodeTensorLaplacian::reflux (int /*crse_amrlev*/,
                               MultiFab& /*res*/, const MultiFab& /*crse_sol*/, const MultiFab& /*crse_rhs*/,
                               MultiFab& /*fine_res*/, MultiFab& /*fine_sol*/, const MultiFab& /*fine_rhs*/) const
{
    amrex::Abort("MLNodeTensorLaplacian::reflux: TODO");
}

void
MLNodeTensorLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLNodeTensorLaplacian::prepareForSolve()");

    MLNodeLinOp::prepareForSolve();

    buildMasks();
}

void
MLNodeTensorLaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLNodeTensorLaplacian::Fapply()");

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const auto s = m_sigma;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real const> const& xarr = in.const_array(mfi);
        Array4<Real> const& yarr = out.array(mfi);

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
        {
            mlndtslap_adotx(tbx,yarr,xarr,s,dxinv);
        });
    }
}

void
MLNodeTensorLaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const
{
    BL_PROFILE("MLNodeTensorLaplacian::Fsmooth()");

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];
    const auto s = m_sigma;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sol); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        Array4<Real> const& solarr = sol.array(mfi);
        Array4<Real const> const& rhsarr = rhs.const_array(mfi);
        Array4<int const> const& dmskarr = dmsk.const_array(mfi);

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
        {
            mlndtslap_gauss_seidel(tbx, solarr, rhsarr, dmskarr, s, dxinv);
        });
    }

    nodalSync(amrlev, mglev, sol);
}

void
MLNodeTensorLaplacian::normalize (int amrlev, int mglev, MultiFab& mf) const
{
    amrex::ignore_unused(amrlev,mglev,mf);
    return;

#if 0
    BL_PROFILE("MLNodeTensorLaplacian::normalize()");

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];
    const auto s = m_sigma;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& arr = mf.array(mfi);
        Array4<int const> const& dmskarr = dmsk.const_array(mfi);

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
        {
            mlndtslap_normalize(tbx, arr, dmskarr, s, dxinv);
        });
    }
#endif
}

void
MLNodeTensorLaplacian::fixUpResidualMask (int /*amrlev*/, iMultiFab& /*resmsk*/)
{
    amrex::Abort("MLNodeTensorLaplacian::fixUpResidualMask: TODO");
}

#ifdef AMREX_USE_HYPRE
void
MLNodeTensorLaplacian::fillIJMatrix (MFIter const& mfi, Array4<HypreNodeLap::Int const> const& nid,
                                     Array4<int const> const& owner,
                                     Vector<HypreNodeLap::Int>& ncols, Vector<HypreNodeLap::Int>& rows,
                                     Vector<HypreNodeLap::Int>& cols, Vector<Real>& mat) const
{
    const int amrlev = 0;
    const int mglev = NMGLevels(amrlev)-1;
    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    Array4<int const> const& dmsk = m_dirichlet_mask[amrlev][mglev]->const_array(mfi);
    const auto s = m_sigma;
    mlndtslap_fill_ijmatrix(mfi.validbox(),nid,owner,ncols,rows,cols,mat,dmsk,s,dxinv);
}
#endif

}
