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
MLNodeLaplacian::buildStencil ()
{
    m_stencil.resize(m_num_amr_levels);
    m_nosigma_stencil.resize(m_num_amr_levels);
    m_s0_norm0.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_stencil[amrlev].resize(m_num_mg_levels[amrlev]);
        m_s0_norm0[amrlev].resize(m_num_mg_levels[amrlev],0.0);
    }

    if (m_coarsening_strategy != CoarseningStrategy::RAP) return;

    const int ncomp_s = (AMREX_SPACEDIM == 2) ? 5 : 9;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMREX_SPACEDIM != 1,
                                     "MLNodeLaplacian::buildStencil: 1d not supported");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_geom[0][0].IsRZ(),
                                     "MLNodeLaplacian::buildStencil: cylindrical not supported for RAP");

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        AMREX_ALWAYS_ASSERT(amrlev == m_num_amr_levels-1 || AMRRefRatio(amrlev) == 2);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            const int nghost = (0 == amrlev && mglev+1 == m_num_mg_levels[amrlev]) ? 1 : 4;
            m_stencil[amrlev][mglev] = std::make_unique<MultiFab>
                (amrex::convert(m_grids[amrlev][mglev], IntVect::TheNodeVector()),
                 m_dmap[amrlev][mglev], ncomp_s, nghost);
            m_stencil[amrlev][mglev]->setVal(0.0);
        }

        if (amrlev > 0) {
            m_nosigma_stencil[amrlev] = std::make_unique<MultiFab>
                (amrex::convert(m_grids[amrlev][0], IntVect::TheNodeVector()),
                 m_dmap[amrlev][0], ncomp_s, 4);
            m_nosigma_stencil[amrlev]->setVal(0.0);
        }

        {
            const Geometry& geom = m_geom[amrlev][0];
            const auto dxinvarr = geom.InvCellSizeArray();

#ifdef AMREX_USE_EB
            const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
            const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
            const MultiFab* intg = m_integral[amrlev].get();
            const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
#endif

            MFItInfo mfi_info;
            if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                FArrayBox sgfab;
                FArrayBox cnfab;

                for (MFIter mfi(*m_stencil[amrlev][0],mfi_info); mfi.isValid(); ++mfi)
                {
                    Box vbx = mfi.validbox();
                    AMREX_D_TERM(vbx.growLo(0,1);, vbx.growLo(1,1);, vbx.growLo(2,1));
                    Box bx = mfi.growntilebox(1);
                    bx &= vbx;
                    const Box& ccbx = amrex::enclosedCells(bx);
                    const Box& ccbxg1 = amrex::grow(ccbx,1);
                    const FArrayBox& sgfab_orig = (*m_sigma[amrlev][0][0])[mfi];
                    Array4<Real const> const& sgarr_orig = sgfab_orig.const_array();

                    Array4<Real> const& starr = m_stencil[amrlev][0]->array(mfi);
                    Array4<Real> const ns_starr = m_nosigma_stencil[amrlev] ?
                        m_nosigma_stencil[amrlev]->array(mfi) : Array4<Real>{};
#ifdef AMREX_USE_EB
                    Array4<Real const> const& intgarr = intg->const_array(mfi);

                    const int ncomp_c = (AMREX_SPACEDIM == 2) ? 6 : 27;
                    bool regular = !factory;
                    if (factory)
                    {
                        Array4<EBCellFlag const> const& flagarr = flags->const_array(mfi);
                        Array4<Real const> const& vfracarr = vfrac->const_array(mfi);
                        const auto& flag = (*flags)[mfi];
                        const auto& typ = flag.getType(ccbxg1);
                        if (typ == FabType::covered)
                        {
                            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx,ncomp_s,i,j,k,n,
                            {
                                starr(i,j,k,n) = 0.0;
                            });

                            if (ns_starr) {
                                AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx,ncomp_s,i,j,k,n,
                                {
                                    ns_starr(i,j,k,n) = 0.0;
                                });
                            }
                        }
                        else if (typ == FabType::singlevalued)
                        {
                            const Box& btmp = ccbxg1 & sgfab_orig.box();

                            cnfab.resize(ccbxg1, ncomp_c);
                            Elixir cneli = cnfab.elixir();
                            Array4<Real> const& cnarr = cnfab.array();

                            sgfab.resize(ccbxg1);
                            Elixir sgeli = sgfab.elixir();
                            Array4<Real> const& sgarr = sgfab.array();

                            AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                            {
                                if (btmp.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                                    mlndlap_set_connection(i,j,k,cnarr,intgarr,vfracarr,flagarr);
                                    sgarr(i,j,k) = sgarr_orig(i,j,k);
                                } else {
                                    for (int n = 0; n < ncomp_c; ++n) {
                                        cnarr(i,j,k,n) = 0.0;
                                    }
                                    sgarr(i,j,k) = 0.0;
                                }
                            });

                            AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                            {
                                mlndlap_set_stencil_eb(i, j, k, starr, sgarr, cnarr, dxinvarr);
                            });

                            if (ns_starr) {
                                AMREX_HOST_DEVICE_FOR_3D(btmp, i, j, k,
                                {
                                    sgarr(i,j,k) = Real(1.0);
                                });
                                AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                                {
                                    mlndlap_set_stencil_eb(i,j,k, ns_starr, sgarr, cnarr, dxinvarr);
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
                        const Box& btmp = ccbxg1 & sgfab_orig.box();

                        sgfab.resize(ccbxg1);
                        Elixir sgeli = sgfab.elixir();
                        Array4<Real> const& sgarr = sgfab.array();

                        AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                        {
                            if (btmp.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                                sgarr(i,j,k) = sgarr_orig(i,j,k);
                            } else {
                                sgarr(i,j,k) = 0.0;
                            }
                        });

                        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                        {
                            mlndlap_set_stencil(tbx,starr,sgarr,dxinvarr);
                        });

                        if (ns_starr) {
                            AMREX_HOST_DEVICE_FOR_3D(btmp, i, j, k,
                            {
                                sgarr(i,j,k) = Real(1.0);
                            });

                            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                            {
                                mlndlap_set_stencil(tbx, ns_starr, sgarr, dxinvarr);
                            });
                        }
                    }
                }
            }

            // set_stencil_s0 has to be in a separate MFIter from set_stencil
            // because it uses other cells' data.
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*m_stencil[amrlev][0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                Array4<Real> const& starr = m_stencil[amrlev][0]->array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_set_stencil_s0(i,j,k,starr);
                });

                if (m_nosigma_stencil[amrlev]) {
                    Array4<Real> const& ns_starr = m_nosigma_stencil[amrlev]->array(mfi);
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                    {
                        mlndlap_set_stencil_s0(i,j,k,ns_starr);
                    });
                }
            }

            m_stencil[amrlev][0]->FillBoundary(geom.periodicity());

            if (m_nosigma_stencil[amrlev]) {
                m_nosigma_stencil[amrlev]->FillBoundary(geom.periodicity());
            }
        }

        for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            const MultiFab& fine = *m_stencil[amrlev][mglev-1];
            MultiFab& crse = *m_stencil[amrlev][mglev];
            bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
            MultiFab cfine;
            if (need_parallel_copy) {
                const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
                cfine.define(ba, fine.DistributionMap(), fine.nComp(), 1);
                cfine.setVal(0.0);
            }

            MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*pcrse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box vbx = mfi.validbox();
                AMREX_D_TERM(vbx.growLo(0,1);, vbx.growLo(1,1);, vbx.growLo(2,1));
                Box bx = mfi.growntilebox(1);
                bx &= vbx;
                Array4<Real> const& csten = pcrse->array(mfi);
                Array4<Real const> const& fsten = fine.const_array(mfi);

                AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                {
                    mlndlap_stencil_rap(i,j,k,csten,fsten);
                });
            }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*pcrse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                Array4<Real> const& starr = pcrse->array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_set_stencil_s0(i,j,k,starr);
                });
            }

            if (need_parallel_copy) {
                crse.ParallelCopy(cfine);
            }

            m_stencil[amrlev][mglev]->FillBoundary(m_geom[amrlev][mglev].periodicity());
        }
    }

#ifdef AMREX_USE_EB
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev) {
#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion() && m_stencil[amrlev][mglev]->isFusingCandidate()) {
                auto const& stma = m_stencil[amrlev][mglev]->const_arrays();
                auto const& dmskma = m_dirichlet_mask[amrlev][mglev]->arrays();
                ParallelFor(*m_stencil[amrlev][mglev],
                [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                {
                    if (stma[box_no](i,j,k,0) == Real(0.0)) {
                        dmskma[box_no](i,j,k) = 1;
                    }
                });
                // We only need to sync once at the end of this function.
            } else
#endif
            {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*m_stencil[amrlev][mglev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    Box const& bx = mfi.tilebox();
                    Array4<Real const> const& starr = m_stencil[amrlev][mglev]->const_array(mfi);
                    Array4<int> const& dmskarr = m_dirichlet_mask[amrlev][mglev]->array(mfi);
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                    {
                        if (starr(i,j,k,0) == Real(0.0)) {
                            dmskarr(i,j,k) = 1;
                        }
                    });
                }
            }
        }
    }
#endif

#ifdef AMREX_USE_EB
    int max_eb_level = 0;
    for (int mglev = m_num_mg_levels[0]-1; mglev > 0; --mglev) {
        int mlo = m_dirichlet_mask[0][mglev]->min(0);
        if (!mlo) {
            // This level is good because not every nodes are Dirichlet.
            max_eb_level = mglev;
            break;
        }
    }
    if (max_eb_level+1 < m_num_mg_levels[0]) {
        resizeMultiGrid(max_eb_level+1);
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

    // This is only needed at the bottom.
    m_s0_norm0[0].back() = m_stencil[0].back()->norm0(0,0) * m_normalization_threshold;
}

}
