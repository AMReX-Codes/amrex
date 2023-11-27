#include <AMReX_MLNodeABecLaplacian.H>
#include <AMReX_MLNodeLap_K.H>
#include <AMReX_MLNodeABecLap_K.H>

namespace amrex {

MLNodeABecLaplacian::MLNodeABecLaplacian (const Vector<Geometry>& a_geom,
                                          const Vector<BoxArray>& a_grids,
                                          const Vector<DistributionMapping>& a_dmap,
                                          const LPInfo& a_info,
                                          const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

void
MLNodeABecLaplacian::define (const Vector<Geometry>& a_geom,
                             const Vector<BoxArray>& a_grids,
                             const Vector<DistributionMapping>& a_dmap,
                             const LPInfo& a_info,
                             const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
#ifdef AMREX_USE_EB
    amrex::Abort("MLNodeABecLaplacian does not support EB");
#endif

    BL_PROFILE("MLNodeABecLaplacian::define()");

    // This makes sure grids are cell-centered;
    Vector<BoxArray> cc_grids = a_grids;
    for (auto& ba : cc_grids) {
        ba.enclosedCells();
    }

    MLNodeLinOp::define(a_geom, cc_grids, a_dmap, a_info, a_factory);

    const int ncomp = getNComp();

    m_a_coeffs.resize(m_num_amr_levels);
    m_b_coeffs.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        m_a_coeffs[amrlev].resize(m_num_mg_levels[amrlev]);
        m_b_coeffs[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev) {
            m_a_coeffs[amrlev][mglev].define
                (amrex::convert(m_grids[amrlev][mglev], IntVect::TheNodeVector()),
                 m_dmap[amrlev][mglev], ncomp, 0);
            m_b_coeffs[amrlev][mglev].define
                (m_grids[amrlev][mglev], m_dmap[amrlev][mglev], ncomp, 1);
        }
    }
}

void
MLNodeABecLaplacian::setACoeffs (int amrlev, Real a_acoef)
{
    m_a_coeffs[amrlev][0].setVal(a_acoef);
    m_needs_update = true;
}

void
MLNodeABecLaplacian::setACoeffs (int amrlev, const MultiFab& a_acoef)
{
    const int ncomp = getNComp();
    m_a_coeffs[amrlev][0].LocalCopy(a_acoef, 0, 0, ncomp, IntVect(0));
    m_needs_update = true;
}

void
MLNodeABecLaplacian::setBCoeffs (int amrlev, Real a_bcoef)
{
    m_b_coeffs[amrlev][0].setVal(a_bcoef);
    m_needs_update = true;
}

void
MLNodeABecLaplacian::setBCoeffs (int amrlev, const MultiFab& a_bcoef)
{
    const int ncomp = getNComp();
    m_b_coeffs[amrlev][0].LocalCopy(a_bcoef, 0, 0, ncomp, IntVect(0));
    m_needs_update = true;
}

void
MLNodeABecLaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLNodeLaplacian::Fapply()");

    AMREX_ALWAYS_ASSERT(getNComp() == 1);

    auto const alpha = m_a_scalar;
    auto const beta  = m_b_scalar;
    const auto dxinvarr = m_geom[amrlev][mglev].InvCellSizeArray();

    auto const& acoef_ma = m_a_coeffs[amrlev][mglev].const_arrays();
    auto const& bcoef_ma = m_b_coeffs[amrlev][mglev].const_arrays();
    auto const& dmskarr_ma = m_dirichlet_mask[amrlev][mglev]->const_arrays();

    auto const& xarr_ma = in.const_arrays();
    auto const& yarr_ma = out.arrays();

    ParallelFor(out, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
    {
        auto lap = mlndlap_adotx_aa(i,j,k,xarr_ma[box_no],bcoef_ma[box_no],dmskarr_ma[box_no],
#if (AMREX_SPACEDIM == 2)
                                    false,
#endif
                                    dxinvarr);
        yarr_ma[box_no](i,j,k) = (dmskarr_ma[box_no](i,j,k)) ? Real(0.0)
            : alpha*acoef_ma[box_no](i,j,k)*xarr_ma[box_no](i,j,k) - beta*lap;
    });
    Gpu::streamSynchronize();
}

void
MLNodeABecLaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const
{
    BL_PROFILE("MLNodeABecLaplacian::Fsmooth()");

    auto const alpha = m_a_scalar;
    auto const beta  = m_b_scalar;
    const auto dxinvarr = m_geom[amrlev][mglev].InvCellSizeArray();

    auto const& acoef = m_a_coeffs[amrlev][mglev];
    auto const& bcoef = m_b_coeffs[amrlev][mglev];
    auto const& dmsk  = *(m_dirichlet_mask[amrlev][mglev]);

#ifdef AMREX_USE_GPU

    auto const& acoef_ma = acoef.const_arrays();
    auto const& bcoef_ma = bcoef.const_arrays();
    auto const& dmskarr_ma = dmsk.const_arrays();
    auto const& solarr_ma = sol.arrays();
    auto const& rhsarr_ma = rhs.const_arrays();

    for (int ns = 0; ns < m_smooth_num_sweeps; ++ns) {
        ParallelFor(sol, [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            auto lap = mlndlap_adotx_aa(i,j,k,solarr_ma[box_no],bcoef_ma[box_no],dmskarr_ma[box_no],
#if (AMREX_SPACEDIM == 2)
                                        false,
#endif
                                        dxinvarr);
            mlndabeclap_jacobi_aa(i,j,k, solarr_ma[box_no], lap, rhsarr_ma[box_no], alpha, beta,
                                  acoef_ma[box_no], bcoef_ma[box_no],
                                  dmskarr_ma[box_no], dxinvarr);
        });
        Gpu::streamSynchronize();
        if (m_smooth_num_sweeps > 1) { nodalSync(amrlev, mglev, sol); }
    }
#else

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(sol); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real const> const& aarr = acoef.array(mfi);
        Array4<Real const> const& barr = bcoef.array(mfi);
        Array4<Real> const& solarr = sol.array(mfi);
        Array4<Real const> const& rhsarr = rhs.const_array(mfi);
        Array4<int const> const& dmskarr = dmsk.const_array(mfi);
        for (int ns = 0; ns < m_smooth_num_sweeps; ++ns) {
            mlndabeclap_gauss_seidel_aa(bx, solarr, rhsarr, alpha, beta,
                                        aarr, barr, dmskarr, dxinvarr);
        }
    }
    nodalSync(amrlev, mglev, sol);
#endif
}

void
MLNodeABecLaplacian::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
    BL_PROFILE("MLNodeABecLaplacian::restriction()");

    applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
        cfine.define(ba, fine.DistributionMap(), 1, 0);
    }

    MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;

    auto pcrse_ma = pcrse->arrays();
    auto fine_ma = fine.const_arrays();
    auto msk_ma = m_dirichlet_mask[amrlev][cmglev-1]->const_arrays();

    ParallelFor(*pcrse, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
    {
        mlndlap_restriction(i,j,k,pcrse_ma[box_no],fine_ma[box_no],msk_ma[box_no]);
    });
    Gpu::streamSynchronize();

    if (need_parallel_copy) {
        crse.ParallelCopy(cfine);
    }
}

void
MLNodeABecLaplacian::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
    BL_PROFILE("MLNodeABecLaplacian::interpolation()");

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    const MultiFab* cmf = &crse;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
        cfine.define(ba, fine.DistributionMap(), 1, 0);
        cfine.ParallelCopy(crse);
        cmf = &cfine;
    }

    auto const& fine_ma = fine.arrays();
    auto const& crse_ma = cmf->const_arrays();
    auto const& msk_ma = m_dirichlet_mask[amrlev][fmglev]->const_arrays();
    auto const& sig_ma = m_b_coeffs[amrlev][fmglev].const_arrays();

    ParallelFor(fine, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
    {
        mlndlap_interpadd_aa(i, j, k, fine_ma[box_no], crse_ma[box_no],
                             sig_ma[box_no], msk_ma[box_no]);
    });
    Gpu::streamSynchronize();
}

void
MLNodeABecLaplacian::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                             const MultiFab& fine_sol, const MultiFab& fine_rhs)
{
    amrex::ignore_unused(camrlev,crse_sol,crse_rhs,fine_sol,fine_rhs);
    amrex::Abort("MLNodeABecLaplacian::averageDownSolutionRHS TODO");
}

void
MLNodeABecLaplacian::reflux (int crse_amrlev,
                             MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                             MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const
{
    amrex::ignore_unused(crse_amrlev,res,crse_sol,crse_rhs,fine_res,fine_sol,fine_rhs);
    amrex::Abort("MLNodeABecLaplacian::reflux TODO");
}

void
MLNodeABecLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLNodeABecLaplacian::prepareForSolve()");

    MLNodeLinOp::prepareForSolve();

    buildMasks();

    averageDownCoeffs();

    m_needs_update = false;
}

void
MLNodeABecLaplacian::update ()
{
    BL_PROFILE("MLNodeABecLaplacian::prepareForSolve()");
    averageDownCoeffs();
    m_needs_update = false;
}

void
MLNodeABecLaplacian::fixUpResidualMask (int amrlev, iMultiFab& resmsk)
{
    if (!m_masks_built) { buildMasks(); }

    auto const& fmsk = m_nd_fine_mask[amrlev]->const_arrays();
    auto const& rmsk = resmsk.arrays();

    amrex::ParallelFor(resmsk,
    [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
    {
        if (fmsk[bno](i,j,k) == crse_fine_node) { rmsk[bno](i,j,k) = 1; }
    });
    Gpu::streamSynchronize();
}

void
MLNodeABecLaplacian::averageDownCoeffs ()
{
    BL_PROFILE("MLNodeABecLaplacian::averageDownCoeffs()");

    for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev) {
        averageDownCoeffsSameAmrLevel(amrlev);
        averageDownCoeffsToCoarseAmrLevel(amrlev);
    }

    averageDownCoeffsSameAmrLevel(0);

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev) {
            m_b_coeffs[amrlev][mglev].FillBoundary(m_geom[amrlev][mglev].periodicity());

            const Box& domain = m_geom[amrlev][mglev].Domain();
            const auto lobc = LoBC();
            const auto hibc = HiBC();

            MFItInfo mfi_info;
            if (Gpu::notInLaunchRegion()) { mfi_info.SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(m_b_coeffs[amrlev][mglev], mfi_info); mfi.isValid(); ++mfi)
            {
                Array4<Real> const& sfab = m_b_coeffs[amrlev][mglev].array(mfi);
                mlndlap_fillbc_cc<Real>(mfi.validbox(),sfab,domain,lobc,hibc);
            }
        }
    }
}

void
MLNodeABecLaplacian::averageDownCoeffsToCoarseAmrLevel (int flev)
{
    const int mglev = 0;
    const int ncomp = getNComp();
    // xxxxx TODO: There is a potential issue of the coarse data not consistent
    // across periodic boundaries.
    amrex::average_down_nodal(m_a_coeffs[flev  ][mglev],
                              m_a_coeffs[flev-1][mglev],
                              IntVect(m_amr_ref_ratio[flev-1]));
    amrex::average_down(m_b_coeffs[flev  ][mglev],
                        m_b_coeffs[flev-1][mglev], 0, ncomp,
                        m_amr_ref_ratio[flev-1]);
}

void
MLNodeABecLaplacian::averageDownCoeffsSameAmrLevel (int amrlev)
{
    const int ncomp = getNComp();
    for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev) {
        IntVect ratio(mg_coarsen_ratio);
        amrex::average_down_nodal(m_a_coeffs[amrlev][mglev-1],
                                  m_a_coeffs[amrlev][mglev  ], ratio);
        amrex::average_down(m_b_coeffs[amrlev][mglev-1],
                            m_b_coeffs[amrlev][mglev  ], 0, ncomp, ratio);
    }
}

}
