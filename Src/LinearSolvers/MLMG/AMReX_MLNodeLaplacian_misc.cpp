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

void
MLNodeLaplacian::averageDownCoeffs ()
{
    BL_PROFILE("MLNodeLaplacian::averageDownCoeffs()");

    if (m_sigma[0][0][0] == nullptr) return;

    if (m_coarsening_strategy == CoarseningStrategy::Sigma)
    {
        for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
        {
            for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
            {
#if (AMREX_SPACEDIM == 1)
                int ndims = 1;
#else
                int ndims = (m_use_harmonic_average || m_use_mapped) ? AMREX_SPACEDIM : 1;
#endif
                for (int idim = 0; idim < ndims; ++idim)
                {
                    if (m_sigma[amrlev][mglev][idim] == nullptr) {
                        if (m_use_harmonic_average && mglev == 0) {
                            m_sigma[amrlev][mglev][idim] = std::make_unique<MultiFab>
                                (*m_sigma[amrlev][mglev][0], amrex::make_alias, 0, 1);
                        } else {
                            m_sigma[amrlev][mglev][idim] = std::make_unique<MultiFab>
                                (m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1);
                            m_sigma[amrlev][mglev][idim]->setVal(0.0);
                        }
                    }
                }
            }
        }
    }

    for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
        averageDownCoeffsSameAmrLevel(amrlev);
        averageDownCoeffsToCoarseAmrLevel(amrlev);
    }

    averageDownCoeffsSameAmrLevel(0);

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        if (m_use_harmonic_average || m_use_mapped) {
            int mglev = 0;
            FillBoundaryCoeff(*m_sigma[amrlev][mglev][0], m_geom[amrlev][mglev]);
            int starting_mglev = m_use_harmonic_average ? 1 : 0;
            for (mglev = starting_mglev; mglev < m_num_mg_levels[amrlev]; ++mglev)
            {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    if (m_sigma[amrlev][mglev][idim]) {
                        FillBoundaryCoeff(*m_sigma[amrlev][mglev][idim], m_geom[amrlev][mglev]);
                    }
                }
            }
        } else {
            int idim = 0;
            for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
            {
                if (m_sigma[amrlev][mglev][idim]) {
                    FillBoundaryCoeff(*m_sigma[amrlev][mglev][idim], m_geom[amrlev][mglev]);
                }
            }
        }
    }
}

void
MLNodeLaplacian::averageDownCoeffsToCoarseAmrLevel (int flev)
{
    if (m_sigma[0][0][0] == nullptr) return;

    const int mglev = 0;
    const int idim = 0;  // other dimensions are just aliases
#ifdef AMREX_USE_EB
    amrex::EB_average_down(*m_sigma[flev][mglev][idim], *m_sigma[flev-1][mglev][idim], 0, 1,
                           m_amr_ref_ratio[flev-1]);
#else
    amrex::average_down(*m_sigma[flev][mglev][idim], *m_sigma[flev-1][mglev][idim], 0, 1,
                        m_amr_ref_ratio[flev-1]);
#endif
}

void
MLNodeLaplacian::averageDownCoeffsSameAmrLevel (int amrlev)
{
    if (m_sigma[0][0][0] == nullptr) return;

    if (m_coarsening_strategy != CoarseningStrategy::Sigma) return;

#if (AMREX_SPACEDIM == 1)
    const int nsigma = 1;
#else
    const int nsigma = (m_use_harmonic_average || m_use_mapped) ? AMREX_SPACEDIM : 1;
#endif

    for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
    {
        int idir = 2;
        bool regular_coarsening = mg_coarsen_ratio_vec[mglev-1] == mg_coarsen_ratio;
        IntVect ratio = mg_coarsen_ratio_vec[mglev-1];
        if (ratio[1] == 1) {
            idir = 1;
        } else if (ratio[0] == 1) {
            idir = 0;
        }
        for (int idim = 0; idim < nsigma; ++idim)
        {
            const MultiFab& fine = *m_sigma[amrlev][mglev-1][idim];
            MultiFab& crse = *m_sigma[amrlev][mglev][idim];
            bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
            MultiFab cfine;
            if (need_parallel_copy) {
                const BoxArray& ba = amrex::coarsen(fine.boxArray(), ratio);
                cfine.define(ba, fine.DistributionMap(), 1, 0);
            }

            MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;

            if (regular_coarsening) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*pcrse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    Array4<Real> const& cfab = pcrse->array(mfi);
                    Array4<Real const> const& ffab = fine.const_array(mfi);
                    if (idim == 0) {
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
                        {
                            mlndlap_avgdown_coeff_x(i,j,k,cfab,ffab);
                        });
                    }
#if (AMREX_SPACEDIM >= 2)
                    else if (idim == 1) {
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
                        {
                            mlndlap_avgdown_coeff_y(i,j,k,cfab,ffab);
                        });
                    }
#if (AMREX_SPACEDIM == 3)
                    else {
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
                        {
                            mlndlap_avgdown_coeff_z(i,j,k,cfab,ffab);
                        });
                    }
#endif
#endif
                }
            } else {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*pcrse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    Array4<Real> const& cfab = pcrse->array(mfi);
                    Array4<Real const> const& ffab = fine.const_array(mfi);
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
                    {
                        mlndlap_semi_avgdown_coeff(i,j,k,cfab,ffab,idir);
                    });
                }
            }
            if (need_parallel_copy) {
                crse.ParallelCopy(cfine);
            }
        }
    }
}

void
MLNodeLaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLNodeLaplacian::Fapply()");

    const auto& sigma = m_sigma[amrlev][mglev];
    const auto& stencil = m_stencil[amrlev][mglev];
    const auto dxinvarr = m_geom[amrlev][mglev].InvCellSizeArray();
#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];


#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto xarr_ma = in.const_arrays();
        auto yarr_ma = out.arrays();
        auto dmskarr_ma = dmsk.const_arrays();

        if (m_coarsening_strategy == CoarseningStrategy::RAP)
        {
            auto stenarr_ma = stencil->const_arrays();
            ParallelFor(out, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
                yarr_ma[box_no](i,j,k) = mlndlap_adotx_sten(i,j,k,xarr_ma[box_no],stenarr_ma[box_no],dmskarr_ma[box_no]);
            });
        }
        else if (sigma[0] == nullptr)
        {
            Real const_sigma = m_const_sigma;
            ParallelFor(out, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
#if (AMREX_SPACEDIM == 2)
                yarr_ma[box_no](i,j,k) = mlndlap_adotx_c(i,j,k,xarr_ma[box_no],const_sigma,dmskarr_ma[box_no], is_rz, dxinvarr);
#else
                yarr_ma[box_no](i,j,k) = mlndlap_adotx_c(i,j,k,xarr_ma[box_no],const_sigma,dmskarr_ma[box_no], dxinvarr);
#endif
            });
        }
        else if ( (m_use_harmonic_average && mglev > 0) ||
                   m_use_mapped )
        {
            AMREX_D_TERM(MultiArray4<Real const> const& sxarr_ma = sigma[0]->const_arrays();,
                         MultiArray4<Real const> const& syarr_ma = sigma[1]->const_arrays();,
                         MultiArray4<Real const> const& szarr_ma = sigma[2]->const_arrays(););
            ParallelFor(out, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
#if (AMREX_SPACEDIM == 2)
                yarr_ma[box_no](i,j,k) = mlndlap_adotx_ha(i,j,k,xarr_ma[box_no],AMREX_D_DECL(sxarr_ma[box_no],syarr_ma[box_no],szarr_ma[box_no]), dmskarr_ma[box_no],
                                               is_rz, dxinvarr);
#else
                yarr_ma[box_no](i,j,k) = mlndlap_adotx_ha(i,j,k,xarr_ma[box_no],AMREX_D_DECL(sxarr_ma[box_no],syarr_ma[box_no],szarr_ma[box_no]), dmskarr_ma[box_no],
                                               dxinvarr);
#endif
            });
        }
        else
        {
            auto sarr_ma = sigma[0]->const_arrays();
            ParallelFor(out, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
            {
#if (AMREX_SPACEDIM == 2)
                yarr_ma[box_no](i,j,k) = mlndlap_adotx_aa(i,j,k,xarr_ma[box_no],sarr_ma[box_no],dmskarr_ma[box_no], is_rz, dxinvarr);
#else
                yarr_ma[box_no](i,j,k) = mlndlap_adotx_aa(i,j,k,xarr_ma[box_no],sarr_ma[box_no],dmskarr_ma[box_no], dxinvarr);
#endif
            });
        }
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(out,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real const> const& xarr = in.const_array(mfi);
            Array4<Real> const& yarr = out.array(mfi);
            Array4<int const> const& dmskarr = dmsk.const_array(mfi);

            if (m_coarsening_strategy == CoarseningStrategy::RAP)
            {
                Array4<Real const> const& stenarr = stencil->const_array(mfi);
                amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                {
                    yarr(i,j,k) = mlndlap_adotx_sten(i,j,k,xarr,stenarr,dmskarr);
                });
            }
            else if (sigma[0] == nullptr)
            {
                Real const_sigma = m_const_sigma;
#if (AMREX_SPACEDIM == 2)
                amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                {
                    yarr(i,j,k) = mlndlap_adotx_c(i,j,k,xarr,const_sigma,dmskarr, is_rz, dxinvarr);
                });
#else
                amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                {
                    yarr(i,j,k) = mlndlap_adotx_c(i,j,k,xarr,const_sigma,dmskarr, dxinvarr);
                });
#endif
            }
            else if ( (m_use_harmonic_average && mglev > 0) ||
                       m_use_mapped )
            {
                AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                             Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                             Array4<Real const> const& szarr = sigma[2]->const_array(mfi););
#if (AMREX_SPACEDIM == 2)
                amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                {
                    yarr(i,j,k) = mlndlap_adotx_ha(i,j,k,xarr,AMREX_D_DECL(sxarr,syarr,szarr), dmskarr,
                                                   is_rz, dxinvarr);
                });
#else
                amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                {
                    yarr(i,j,k) = mlndlap_adotx_ha(i,j,k,xarr,AMREX_D_DECL(sxarr,syarr,szarr), dmskarr,
                                                   dxinvarr);
                });
#endif
            }
            else
            {
                Array4<Real const> const& sarr = sigma[0]->const_array(mfi);
#if (AMREX_SPACEDIM == 2)
                amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                {
                    yarr(i,j,k) = mlndlap_adotx_aa(i,j,k,xarr,sarr,dmskarr, is_rz, dxinvarr);
                });
#else
                amrex::LoopConcurrentOnCpu(bx, [&] (int i, int j, int k) noexcept
                {
                    yarr(i,j,k) = mlndlap_adotx_aa(i,j,k,xarr,sarr,dmskarr, dxinvarr);
                });
#endif
           }
        }
    }
}

void
MLNodeLaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const
{
    BL_PROFILE("MLNodeLaplacian::Fsmooth()");

    const auto& sigma = m_sigma[amrlev][mglev];
    const auto& stencil = m_stencil[amrlev][mglev];
    const auto dxinvarr = m_geom[amrlev][mglev].InvCellSizeArray();
#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        auto solarr_ma = sol.arrays();
        auto rhsarr_ma = rhs.const_arrays();
        auto dmskarr_ma = dmsk.const_arrays();
        if (m_coarsening_strategy == CoarseningStrategy::RAP)
        {
            auto starr_ma = stencil->const_arrays();
            for (int ns = 0; ns < m_smooth_num_sweeps; ++ns)
            {
                ParallelFor(sol, [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                {
                    Real Ax = mlndlap_adotx_sten(i,j,k,solarr_ma[box_no],starr_ma[box_no],dmskarr_ma[box_no]);
                    mlndlap_jacobi_sten(i,j,k,solarr_ma[box_no],Ax,rhsarr_ma[box_no],starr_ma[box_no],dmskarr_ma[box_no]);
                });
            }
        }
        else if (sigma[0] == nullptr)
        {
            for (int ns = 0; ns < m_smooth_num_sweeps; ++ns)
            {
                Real const_sigma = m_const_sigma;
                ParallelFor(sol, [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                {
                    Real Ax = mlndlap_adotx_c(i,j,k,solarr_ma[box_no],const_sigma,dmskarr_ma[box_no],
#if (AMREX_SPACEDIM == 2)
                                              is_rz,
#endif
                                              dxinvarr);
                    mlndlap_jacobi_c(i,j,k, solarr_ma[box_no], Ax, rhsarr_ma[box_no], const_sigma,
                                     dmskarr_ma[box_no], dxinvarr);
                });
            }
        }
        else if ((m_use_harmonic_average && mglev > 0) || m_use_mapped)
        {
            AMREX_D_TERM(MultiArray4<Real const> const& sxarr_ma = sigma[0]->const_arrays();,
                         MultiArray4<Real const> const& syarr_ma = sigma[1]->const_arrays();,
                         MultiArray4<Real const> const& szarr_ma = sigma[2]->const_arrays(););
            for (int ns = 0; ns < m_smooth_num_sweeps; ++ns)
            {
                ParallelFor(sol, [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                {
                    Real Ax = mlndlap_adotx_ha(i,j,k,solarr_ma[box_no],AMREX_D_DECL(sxarr_ma[box_no],syarr_ma[box_no],szarr_ma[box_no]), dmskarr_ma[box_no],
#if (AMREX_SPACEDIM == 2)
                                               is_rz,
#endif
                                               dxinvarr);
                    mlndlap_jacobi_ha(i,j,k, solarr_ma[box_no], Ax, rhsarr_ma[box_no], AMREX_D_DECL(sxarr_ma[box_no],syarr_ma[box_no],szarr_ma[box_no]),
                                      dmskarr_ma[box_no], dxinvarr);
                });
            }
        }
        else
        {
            auto sarr_ma = sigma[0]->const_arrays();
            for (int ns = 0; ns < m_smooth_num_sweeps; ++ns)
            {
                ParallelFor(sol, [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                {
                    Real Ax = mlndlap_adotx_aa(i,j,k,solarr_ma[box_no],sarr_ma[box_no],dmskarr_ma[box_no],
#if (AMREX_SPACEDIM == 2)
                                               is_rz,
#endif
                                               dxinvarr);
                    mlndlap_jacobi_aa(i,j,k, solarr_ma[box_no], Ax, rhsarr_ma[box_no], sarr_ma[box_no],
                                      dmskarr_ma[box_no], dxinvarr);
                });
            }
        }

        Gpu::streamSynchronize();
        if (m_smooth_num_sweeps > 1) nodalSync(amrlev, mglev, sol);
    }
    else // cpu
#endif
    {
        bool regular_coarsening = true;
        if (amrlev == 0 && mglev > 0)
        {
            regular_coarsening = mg_coarsen_ratio_vec[mglev-1] == mg_coarsen_ratio;
        }
        if (sigma[0] == nullptr) {
            AMREX_ALWAYS_ASSERT(regular_coarsening);
        }

        if (m_use_gauss_seidel)
        {
            if (m_coarsening_strategy == CoarseningStrategy::RAP)
            {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
                for (MFIter mfi(sol); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.validbox();
                    Array4<Real> const& solarr = sol.array(mfi);
                    Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                    Array4<Real const> const& starr = stencil->const_array(mfi);
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

                    for (int ns = 0; ns < m_smooth_num_sweeps; ++ns) {
                        mlndlap_gauss_seidel_sten(bx,solarr,rhsarr,starr,dmskarr);
                    }
                }
            }
            else if (sigma[0] == nullptr)
            {
                Real const_sigma = m_const_sigma;
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
                for (MFIter mfi(sol); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.validbox();
                    Array4<Real> const& solarr = sol.array(mfi);
                    Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

                    for (int ns = 0; ns < m_smooth_num_sweeps; ++ns) {
                        mlndlap_gauss_seidel_c(bx, solarr, rhsarr,
                                               const_sigma, dmskarr, dxinvarr
#if (AMREX_SPACEDIM == 2)
                                               ,is_rz
#endif
                            );
                    }
                }
            }
            else if ( (m_use_harmonic_average && mglev > 0) || m_use_mapped )
            {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
                for (MFIter mfi(sol); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.validbox();
                    AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                                 Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                                 Array4<Real const> const& szarr = sigma[2]->const_array(mfi););
                    Array4<Real> const& solarr = sol.array(mfi);
                    Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

                    for (int ns = 0; ns < m_smooth_num_sweeps; ++ns) {
                        mlndlap_gauss_seidel_ha(bx, solarr, rhsarr,
                                                AMREX_D_DECL(sxarr,syarr,szarr),
                                                dmskarr, dxinvarr
#if (AMREX_SPACEDIM == 2)
                                                ,is_rz
#endif
                            );
                    }
                }
            }
            else
            {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
                for (MFIter mfi(sol); mfi.isValid(); ++mfi)
                {

                    const Box& bx = mfi.validbox();
                    Array4<Real const> const& sarr = sigma[0]->const_array(mfi);
                    Array4<Real> const& solarr = sol.array(mfi);
                    Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

                    if ( regular_coarsening )
                    {
                        for (int ns = 0; ns < m_smooth_num_sweeps; ++ns) {
                            mlndlap_gauss_seidel_aa(bx, solarr, rhsarr,
                                                    sarr, dmskarr, dxinvarr
#if (AMREX_SPACEDIM == 2)
                                                   ,is_rz
#endif
                                 );
                        }
                    } else {
                        for (int ns = 0; ns < m_smooth_num_sweeps; ++ns) {
                            mlndlap_gauss_seidel_with_line_solve_aa(bx, solarr, rhsarr,
                                                                    sarr, dmskarr, dxinvarr
#if (AMREX_SPACEDIM == 2)
                                                                   ,is_rz
#endif
                                );
                        }
                    }
                }
            }

            nodalSync(amrlev, mglev, sol);
        }
        else
        {
            MultiFab Ax(sol.boxArray(), sol.DistributionMap(), 1, 0);
            Fapply(amrlev, mglev, Ax, sol);

            if (m_coarsening_strategy == CoarseningStrategy::RAP)
            {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
                for (MFIter mfi(sol,true); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    Array4<Real> const& solarr = sol.array(mfi);
                    Array4<Real const> const& Axarr = Ax.const_array(mfi);
                    Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                    Array4<Real const> const& stenarr = stencil->const_array(mfi);
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

                    mlndlap_jacobi_sten(bx,solarr,Axarr,rhsarr,stenarr,dmskarr);
                }
            }
            else if (sigma[0] == nullptr)
            {
                Real const_sigma = m_const_sigma;
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
                for (MFIter mfi(sol,true); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    Array4<Real> const& solarr = sol.array(mfi);
                    Array4<Real const> const& Axarr = Ax.const_array(mfi);
                    Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

                    mlndlap_jacobi_c (bx, solarr, Axarr, rhsarr, const_sigma,
                                      dmskarr, dxinvarr);
                }
            }
            else if ( (m_use_harmonic_average && mglev > 0) || m_use_mapped )
            { // NOLINT(bugprone-branch-clone)
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
                for (MFIter mfi(sol,true); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                                 Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                                 Array4<Real const> const& szarr = sigma[2]->const_array(mfi););
                    Array4<Real> const& solarr = sol.array(mfi);
                    Array4<Real const> const& Axarr = Ax.const_array(mfi);
                    Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

                    mlndlap_jacobi_ha (bx, solarr, Axarr, rhsarr, AMREX_D_DECL(sxarr,syarr,szarr),
                                       dmskarr, dxinvarr);
                }
            }
            else
            {
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
                for (MFIter mfi(sol,true); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    Array4<Real const> const& sarr = sigma[0]->const_array(mfi);
                    Array4<Real> const& solarr = sol.array(mfi);
                    Array4<Real const> const& Axarr = Ax.const_array(mfi);
                    Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

                    mlndlap_jacobi_aa (bx, solarr, Axarr, rhsarr, sarr,
                                       dmskarr, dxinvarr);
                }
            }
        }
    }
}

void
MLNodeLaplacian::updateVelocity (const Vector<MultiFab*>& vel, const Vector<MultiFab const*>& sol) const
{
#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        const auto& sigma = m_sigma[amrlev][0][0];
        const auto dxinv = m_geom[amrlev][0].InvCellSizeArray();
#ifdef AMREX_USE_EB
        const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
        const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
        const MultiFab* intg = m_integral[amrlev].get();
#endif
        for (MFIter mfi(*vel[amrlev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& varr = vel[amrlev]->array(mfi);
            Array4<Real const> const& solarr = sol[amrlev]->const_array(mfi);
#ifdef AMREX_USE_EB
            bool regular = !factory;
            if (factory)
            {
                Array4<Real const> const& sigmaarr = sigma->const_array(mfi);
                auto type = (*flags)[mfi].getType(bx);
                Array4<Real const> const& vfracarr = vfrac->const_array(mfi);
                Array4<Real const> const& intgarr = intg->const_array(mfi);
                if (type == FabType::covered)
                {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, AMREX_SPACEDIM, i, j, k, n,
                    {
                        varr(i,j,k,n) = 0.0;
                    });
                }
                else if (type == FabType::singlevalued)
                {
                    AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                    {
                        mlndlap_mknewu_eb(i,j,k, varr, solarr, sigmaarr, vfracarr, intgarr, dxinv);
                    });
                }
                else
                {
                    regular = true;
                }
            }
            if (regular)
#endif
            {
                if (sigma) {
                    Array4<Real const> const& sigmaarr = sigma->const_array(mfi);
#if (AMREX_SPACEDIM == 2)
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                    {
                        mlndlap_mknewu(i,j,k,varr,solarr,sigmaarr,dxinv,is_rz);
                    });
#else
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                    {
                        mlndlap_mknewu(i,j,k,varr,solarr,sigmaarr,dxinv);
                    });
#endif
                } else {
                    Real const_sigma = m_const_sigma;
#if (AMREX_SPACEDIM == 2)
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                    {
                        mlndlap_mknewu_c(i,j,k,varr,solarr,const_sigma,dxinv,is_rz);
                    });
#else
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                    {
                        mlndlap_mknewu_c(i,j,k,varr,solarr,const_sigma,dxinv);
                    });
#endif
                }
            }
        }
    }
}

void
MLNodeLaplacian::compGrad (int amrlev, MultiFab& grad, MultiFab& sol) const
{
#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

    Real sigma = Real(-1.0);

    AMREX_ASSERT(grad.nComp() >= AMREX_SPACEDIM);

    const auto dxinv = m_geom[amrlev][0].InvCellSizeArray();
#ifdef AMREX_USE_EB
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
    const MultiFab* intg = m_integral[amrlev].get();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(grad, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& garr = grad.array(mfi);
        Array4<Real const> const& solarr = sol.const_array(mfi);

        AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, AMREX_SPACEDIM, i, j, k, n,
        {
            garr(i,j,k,n) = 0.0;
        });

#ifdef AMREX_USE_EB
        bool regular = !factory;
        if (factory)
        {
            auto type = (*flags)[mfi].getType(bx);
            Array4<Real const> const& vfracarr = vfrac->const_array(mfi);
            Array4<Real const> const& intgarr = intg->const_array(mfi);
            if (type == FabType::covered)
            { }
            else if (type == FabType::singlevalued)
            {
              AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
              {
                  mlndlap_mknewu_eb_c(i,j,k, garr, solarr, sigma, vfracarr, intgarr, dxinv);
              });
            }
            else
            {
                regular = true;
            }
        }
        if (regular)
#endif
        {

#if (AMREX_SPACEDIM == 2)
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
            {
                mlndlap_mknewu_c(i,j,k,garr,solarr,sigma,dxinv,is_rz);
            });
#else
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
            {
                mlndlap_mknewu_c(i,j,k,garr,solarr,sigma,dxinv);
            });
#endif
        }
    }
}

void
MLNodeLaplacian::getFluxes (const Vector<MultiFab*> & a_flux, const Vector<MultiFab*>& a_sol) const
{
#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

    AMREX_ASSERT(a_flux[0]->nComp() >= AMREX_SPACEDIM);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        const auto& sigma = m_sigma[amrlev][0][0];
        const auto dxinv = m_geom[amrlev][0].InvCellSizeArray();
#ifdef AMREX_USE_EB
        const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
        const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
        const MultiFab* intg = m_integral[amrlev].get();
#endif

        // Initialize to zero because we only want -(sigma * grad(phi))

        for (MFIter mfi(*a_flux[amrlev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& farr = a_flux[amrlev]->array(mfi);
            Array4<Real const> const& solarr = a_sol[amrlev]->const_array(mfi);

            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, AMREX_SPACEDIM, i, j, k, n,
            {
                farr(i,j,k,n) = 0.0;
            });

#ifdef AMREX_USE_EB
            bool regular = !factory;
            if (factory)
            {
                Array4<Real const> const& sigmaarr = sigma->array(mfi);
                auto type = (*flags)[mfi].getType(bx);
                Array4<Real const> const& vfracarr = vfrac->const_array(mfi);
                Array4<Real const> const& intgarr = intg->const_array(mfi);
                if (type == FabType::covered)
                { }
                else if (type == FabType::singlevalued)
                {
                    AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                    {
                        mlndlap_mknewu_eb(i,j,k, farr, solarr, sigmaarr, vfracarr, intgarr, dxinv);
                    });
                }
                else
                {
                    regular = true;
                }
            }
            if (regular)
#endif
            {
                if (sigma) {
                    Array4<Real const> const& sigmaarr = sigma->array(mfi);
#if (AMREX_SPACEDIM == 2)
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                    {
                        mlndlap_mknewu(i,j,k,farr,solarr,sigmaarr,dxinv,is_rz);
                    });
#else
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                    {
                        mlndlap_mknewu(i,j,k,farr,solarr,sigmaarr,dxinv);
                    });
#endif
                } else {
                    Real const_sigma = m_const_sigma;
#if (AMREX_SPACEDIM == 2)
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                    {
                        mlndlap_mknewu_c(i,j,k,farr,solarr,const_sigma,dxinv,is_rz);
                    });
#else
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                    {
                        mlndlap_mknewu_c(i,j,k,farr,solarr,const_sigma,dxinv);
                    });
#endif
                }
            }
        }
    }
}

void
MLNodeLaplacian::compDivergence (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel)
{
    compRHS(rhs, vel, Vector<const MultiFab*>(), Vector<MultiFab*>());
}

void
MLNodeLaplacian::compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel,  // NOLINT(readability-convert-member-functions-to-static)
                          const Vector<const MultiFab*>& rhnd,
                          const Vector<MultiFab*>& a_rhcc)
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(rhs,vel,rhnd,a_rhcc);
#else
    //
    // Note that div vel we copmute on a coarse/fine nodes is not a
    // composite divergence.  It has been restricted so that it is suitable
    // as RHS for our geometric mulitgrid solver with a MG hirerachy
    // including multiple AMR levels.
    //
    // Also note that even for RAP, we do doubling at Nuemann boundary,
    // because unimposeNeumannBC will be called on rhs for RAP.
    //

    BL_PROFILE("MLNodeLaplacian::compRHS()");

    if (!m_masks_built) buildMasks();

#ifdef AMREX_USE_EB
    if (!m_integral_built) buildIntegral();
    if (m_build_surface_integral && !m_surface_integral_built) buildSurfaceIntegral();
#endif

#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

    const auto lobc = LoBC();
    const auto hibc = HiBC();

    bool has_inflow = false;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        has_inflow = has_inflow || (lobc[idim] == LinOpBCType::inflow ||
                                    hibc[idim] == LinOpBCType::inflow);
    }

    Vector<std::unique_ptr<MultiFab> > rhcc(m_num_amr_levels);
    Vector<std::unique_ptr<MultiFab> > rhs_cc(m_num_amr_levels);

    for (int ilev = 0; ilev < m_num_amr_levels; ++ilev)
    {
        const Geometry& geom = m_geom[ilev][0];
        AMREX_ASSERT(vel[ilev]->nComp() >= AMREX_SPACEDIM);
        AMREX_ASSERT(vel[ilev]->nGrow() >= 1);

        if (has_inflow) { // Zero out transverse velocity so that it's not seen.
            Box domain = geom.Domain();
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (lobc[idim] != LinOpBCType::inflow) {
                    domain.growLo(idim,1);
                }
                if (hibc[idim] != LinOpBCType::inflow) {
                    domain.growHi(idim,1);
                }
            }
            const auto dlo = domain.smallEnd();
            const auto dhi = domain.bigEnd();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*vel[ilev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox(1);
                Array4<Real> const& vfab = vel[ilev]->array(mfi);
                if ( ! domain.contains(bx) ) {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                    {
                        IntVect cell(AMREX_D_DECL(i,j,k));
                        for (int in = 0; in < AMREX_SPACEDIM; ++in) {
                            for (int it = 0; it < AMREX_SPACEDIM; ++it) {
                                if (it != in) {
                                    if (cell[in] < dlo[in] || cell[in] > dhi[in]) {
                                        vfab(i,j,k,it) = Real(0.0);
                                    }
                                }
                            }
                        }
                    });
                }
            }
        }

        vel[ilev]->FillBoundary(0, AMREX_SPACEDIM, IntVect(1), geom.periodicity());

        if (ilev < a_rhcc.size() && a_rhcc[ilev])
        {
            rhcc[ilev] = std::make_unique<MultiFab>(a_rhcc[ilev]->boxArray(),
                                                    a_rhcc[ilev]->DistributionMap(), 1, 1);
            rhcc[ilev]->setVal(0.0);
            MultiFab::Copy(*rhcc[ilev], *a_rhcc[ilev], 0, 0, 1, 0);
            rhcc[ilev]->FillBoundary(geom.periodicity());

            rhs_cc[ilev] = std::make_unique<MultiFab>(rhs[ilev]->boxArray(),
                                                      rhs[ilev]->DistributionMap(), 1, 0);
        }

        const auto dxinvarr = geom.InvCellSizeArray();
        const Box& nddom = amrex::surroundingNodes(geom.Domain());

        const iMultiFab& dmsk = *m_dirichlet_mask[ilev][0];

#ifdef AMREX_USE_EB
        const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[ilev][0].get());
        const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
        const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
        const MultiCutFab* barea = (factory) ? &(factory->getBndryArea()) : nullptr;
        const MultiFab* intg = m_integral[ilev].get();
        const MultiFab* sintg = m_surface_integral[ilev].get();

        AMREX_ALWAYS_ASSERT(ilev == m_num_amr_levels-1 || AMRRefRatio(ilev) == 2
                            || factory == nullptr || factory->isAllRegular());
#endif

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*rhs[ilev],mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& rhsarr = rhs[ilev]->array(mfi);
            Array4<Real const> const& velarr = vel[ilev]->const_array(mfi);
            Array4<int const> const& dmskarr = dmsk.const_array(mfi);

#ifdef AMREX_USE_EB
            bool regular = !factory;
            FabType typ = FabType::regular;
            if (factory)
            {
                const auto& flag = (*flags)[mfi];
                const auto& ccbx = amrex::grow(amrex::enclosedCells(bx),1);
                typ = flag.getType(ccbx);
                if (typ == FabType::covered)
                {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                    {
                        rhsarr(i,j,k) = 0.0;
                    });
                }
                else if (typ == FabType::singlevalued)
                {
                    Array4<Real const> const& vfracarr = vfrac->const_array(mfi);
                    Array4<Real const> const& intgarr = intg->const_array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                    {
                        mlndlap_divu_eb(i,j,k,rhsarr,velarr,vfracarr,intgarr,dmskarr,dxinvarr,nddom,lobc,hibc);
                    });

                    if (m_eb_vel_dot_n[ilev]) {
                        Array4<Real const> const& eb_vel_dot_n = m_eb_vel_dot_n[ilev]->const_array(mfi);
                        Array4<Real const> const& bareaarr = barea->const_array(mfi);
                        Array4<Real const> const& sintgarr = sintg->const_array(mfi);

                        AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                        {
                            add_eb_flow_contribution(i,j,k,rhsarr,dmskarr,
                                dxinvarr,bareaarr,sintgarr,eb_vel_dot_n);
                        });
                    }
                }
                else
                {
                    regular = true;
                }
            }
            if (regular)
#endif
            {
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                {
                    mlndlap_divu(i,j,k,rhsarr,velarr,dmskarr,dxinvarr,
                                 nddom,lobc,hibc,is_rz);
                });
#else
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                {
                    mlndlap_divu(i,j,k,rhsarr,velarr,dmskarr,dxinvarr,
                                 nddom,lobc,hibc);
                });
#endif
            }

            mlndlap_impose_neumann_bc(bx, rhsarr, nddom, lobc, hibc);

            if (rhcc[ilev])
            {
                Array4<Real> const& rhs_cc_a = rhs_cc[ilev]->array(mfi);
                Array4<Real const> const& rhccarr = rhcc[ilev]->const_array(mfi);
#ifdef AMREX_USE_EB
                if (typ == FabType::singlevalued)
                {
                    Array4<Real const> const& vfracarr = vfrac->const_array(mfi);
                    Array4<Real const> const& intgarr = intg->const_array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                    {
                        rhs_cc_a(i,j,k) = mlndlap_rhcc_eb(i,j,k,rhccarr,vfracarr,intgarr,dmskarr);
                    });
                }
                else if (typ == FabType::covered)
                {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                    {
                        rhs_cc_a(i,j,k) = 0.0;
                    });
                }
                else
#endif
                {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                    {
                        rhs_cc_a(i,j,k) = mlndlap_rhcc(i, j, k, rhccarr, dmskarr);
                    });
                }

                mlndlap_impose_neumann_bc(bx, rhs_cc_a, nddom, lobc, hibc);
            }
        }
    }

    Vector<std::unique_ptr<MultiFab> > frhs(m_num_amr_levels);

    for (int ilev = 0; ilev < m_num_amr_levels-1; ++ilev)
    {
        const int amrrr = AMRRefRatio(ilev);
        const Geometry& fgeom = m_geom[ilev+1][0];
        AMREX_ALWAYS_ASSERT(amrrr == 2 || amrrr == 4);

        frhs[ilev] = std::make_unique<MultiFab>(amrex::coarsen(rhs[ilev+1]->boxArray(),amrrr),
                                                rhs[ilev+1]->DistributionMap(), 1, 0);

        const Box& ccfdom = fgeom.Domain();
        const auto fdxinv = fgeom.InvCellSizeArray();
        const iMultiFab& fdmsk = *m_dirichlet_mask[ilev+1][0];

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*frhs[ilev],mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& cbx = mfi.tilebox();
            const Box& fvbx = amrex::refine(mfi.validbox(),amrrr);
            const Box& cc_fvbx = amrex::enclosedCells(fvbx);

            Box bx_vel = cc_fvbx;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                if (m_lobc[0][idim] == LinOpBCType::inflow)
                {
                    if (bx_vel.smallEnd(idim) == ccfdom.smallEnd(idim)) {
                        bx_vel.growLo(idim, 1);
                    }
                }
                if (m_hibc[0][idim] == LinOpBCType::inflow)
                {
                    if (bx_vel.bigEnd(idim) == ccfdom.bigEnd(idim)) {
                        bx_vel.growHi(idim, 1);
                    }
                }
            }

            Array4<Real> const& rhsarr = frhs[ilev]->array(mfi);
            Array4<Real const> const& velarr = vel[ilev+1]->const_array(mfi);
            Array4<Real const> const& rhsarr_fine = rhs[ilev+1]->const_array(mfi);
            Array4<int const> const& mskarr = fdmsk.const_array(mfi);
#if (AMREX_SPACEDIM == 2)
            if (amrrr == 2) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(cbx, i, j, k,
                {
                    mlndlap_divu_fine_contrib<2>(i,j,k,fvbx,bx_vel,rhsarr,velarr,rhsarr_fine,
                                                 mskarr,fdxinv,is_rz);
                });
            } else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(cbx, i, j, k,
                {
                    mlndlap_divu_fine_contrib<4>(i,j,k,fvbx,bx_vel,rhsarr,velarr,rhsarr_fine,
                                                 mskarr,fdxinv,is_rz);
                });
            }
#elif (AMREX_SPACEDIM == 3)
            if (amrrr == 2) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(cbx, i, j, k,
                {
                    mlndlap_divu_fine_contrib<2>(i,j,k,fvbx,bx_vel,rhsarr,velarr,rhsarr_fine,
                                                 mskarr,fdxinv);
                });
            } else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(cbx, i, j, k,
                {
                    mlndlap_divu_fine_contrib<4>(i,j,k,fvbx,bx_vel,rhsarr,velarr,rhsarr_fine,
                                                 mskarr,fdxinv);
                });
            }
#endif

            if (rhcc[ilev+1])
            {
                // xxxxx TODO: incorrect if cut cells are too close to coarse/fine boundary
                Array4<Real const> const& rhccarr = rhcc[ilev+1]->const_array(mfi);
                if (amrrr == 2) {
                    AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                    {
                        mlndlap_rhcc_fine_contrib<2>(i,j,k,cc_fvbx,rhsarr,rhccarr,mskarr);
                    });
                } else {
                    AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                    {
                        mlndlap_rhcc_fine_contrib<4>(i,j,k,cc_fvbx,rhsarr,rhccarr,mskarr);
                    });
                }
            }
        }
    }

    for (int ilev = 0; ilev < m_num_amr_levels; ++ilev)
    {
        if (rhs_cc[ilev]) {
            MultiFab::Add(*rhs[ilev], *rhs_cc[ilev], 0, 0, 1, 0);
        }
    }

    for (int ilev = m_num_amr_levels-2; ilev >= 0; --ilev)
    {
        const Geometry& cgeom = m_geom[ilev][0];

        MultiFab crhs(rhs[ilev]->boxArray(), rhs[ilev]->DistributionMap(), 1, 0);
        crhs.setVal(0.0);
        crhs.ParallelAdd(*frhs[ilev], cgeom.periodicity());

        const Box& cccdom = cgeom.Domain();
        const Box& cccdom_p = cgeom.growPeriodicDomain(1);
        Box cccdom_pi = cccdom_p;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (lobc[idim] == LinOpBCType::inflow) {
                cccdom_pi.growLo(idim,1);
            }
            if (hibc[idim] == LinOpBCType::inflow) {
                cccdom_pi.growHi(idim,1);
            }
        }
        const Box& cnddom = amrex::surroundingNodes(cccdom);
        const auto cdxinv = cgeom.InvCellSizeArray();
        const iMultiFab& cdmsk = *m_dirichlet_mask[ilev][0];
        const iMultiFab& c_nd_mask = *m_nd_fine_mask[ilev];
        const iMultiFab& c_cc_mask = *m_cc_fine_mask[ilev];
        const auto& has_fine_bndry = *m_has_fine_bndry[ilev];

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*rhs[ilev],mfi_info); mfi.isValid(); ++mfi)
        {
            if (has_fine_bndry[mfi])
            {
                const Box& bx = mfi.tilebox();

                Array4<Real> const& rhsarr = rhs[ilev]->array(mfi);
                Array4<Real const> const& velarr = vel[ilev]->const_array(mfi);
                Array4<Real const> const& crhsarr = crhs.const_array(mfi);
                Array4<int const> const& cdmskarr = cdmsk.const_array(mfi);
                Array4<int const> const& ndmskarr = c_nd_mask.const_array(mfi);
                Array4<int const> const& ccmskarr = c_cc_mask.const_array(mfi);

                Array4<Real const> rhccarr = (rhcc[ilev])
                    ? rhcc[ilev]->const_array(mfi) : Array4<Real const>{};
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                {
                    mlndlap_divu_cf_contrib(i,j,k,rhsarr,velarr,crhsarr,rhccarr,
                                            cdmskarr,ndmskarr,ccmskarr,
                                            is_rz,
                                            cdxinv,cccdom_p,cccdom_pi,cnddom,lobc,hibc);
                });
#elif (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                {
                    mlndlap_divu_cf_contrib(i,j,k,rhsarr,velarr,crhsarr,rhccarr,
                                            cdmskarr,ndmskarr,ccmskarr,
                                            cdxinv,cccdom_p,cccdom_pi,cnddom,lobc,hibc);
                });
#endif
            }
        }
    }

    for (int ilev = 0; ilev < m_num_amr_levels; ++ilev)
    {
        if (ilev < rhnd.size() && rhnd[ilev]) {
            MultiFab::Add(*rhs[ilev], *rhnd[ilev], 0, 0, 1, 0);
        }
    }

#ifdef AMREX_USE_EB
    for (int ilev = 0; ilev < m_num_amr_levels; ++ilev) {
        amrex::EB_set_covered(*rhs[ilev], 0.0);
    }
#endif
#endif
}

}
