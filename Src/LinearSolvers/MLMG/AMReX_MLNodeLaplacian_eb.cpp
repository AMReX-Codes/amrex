#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLNodeLap_K.H>

#ifdef AMREX_USE_EB
#include <AMReX_algoim.H>
#endif

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

namespace amrex {

template <typename MF>
void
MLNodeLaplacianT<MF>::buildIntegral ()
{
    if constexpr (std::is_same_v<MF,MultiFab>) {
        if (m_integral_built) { return; }

        BL_PROFILE("MLNodeLaplacian::buildIntegral()");

        m_integral_built = true;

#if (AMREX_SPACEDIM == 2)
        for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev)
        {
            MF* intg = m_integral[amrlev].get();

            const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(this->m_factory[amrlev][0].get());
            if (factory)
            {
                const int ncomp = intg->nComp();
                const auto& flags = factory->getMultiEBCellFlagFab();
                const auto& vfrac = factory->getVolFrac();
                const auto& area = factory->getAreaFrac();
                const auto& bcent = factory->getBndryCent();

                MFItInfo mfi_info;
                if (Gpu::notInLaunchRegion()) { mfi_info.EnableTiling().SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*intg,mfi_info); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.growntilebox();
                    Array4<RT> const& garr = intg->array(mfi);
                    const auto& flag = flags[mfi];
                    auto typ = flag.getType(bx);

                    if (typ == FabType::covered) {
                        AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                        {
                            garr(i,j,k,n) = 0.0;
                        });
                    } else if (typ == FabType::regular) {
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                        {
                            mlndlap_set_integral(i,j,k,garr);
                        });
                    } else {
                        Array4<EBCellFlag const> const& flagarr = flags.const_array(mfi);
                        Array4<RT const> const& vfracarr = vfrac.const_array(mfi);
                        Array4<RT const> const& axarr = area[0]->const_array(mfi);
                        Array4<RT const> const& ayarr = area[1]->const_array(mfi);
                        Array4<RT const> const& bcarr = bcent.const_array(mfi);
                        AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                        {
                            mlndlap_set_integral_eb(i,j,k,garr,flagarr,vfracarr,axarr,ayarr,bcarr);
                        });
                    }
                }
            }
        }
#else
        for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev)
        {
            if (dynamic_cast<EBFArrayBoxFactory const*>(this->m_factory[amrlev][0].get())) {
                amrex::algoim::compute_integrals(*m_integral[amrlev]);
            }
        }
#endif
    }
    // The single precision copy has no EB factory, so it has no integrals.
}

template <typename MF>
void
MLNodeLaplacianT<MF>::buildSurfaceIntegral ()
{
    if constexpr (std::is_same_v<MF,MultiFab>) {
        if (m_surface_integral_built) { return; }

        BL_PROFILE("MLNodeLaplacian::buildSurfaceIntegral()");

        m_surface_integral_built = true;

#if (AMREX_SPACEDIM == 2)
        for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev)
        {
            MF* sintg = m_surface_integral[amrlev].get();

            const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(this->m_factory[amrlev][0].get());
            if (factory && sintg)
            {
                const int ncomp = sintg->nComp();
                const auto& flags = factory->getMultiEBCellFlagFab();
                const auto& bcent = factory->getBndryCent();
                const auto& barea = factory->getBndryArea();
                const auto& bnorm = factory->getBndryNormal();

                MFItInfo mfi_info;
                if (Gpu::notInLaunchRegion()) { mfi_info.EnableTiling().SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*sintg,mfi_info); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.growntilebox();
                    Array4<RT> const& garr = sintg->array(mfi);
                    const auto& flag = flags[mfi];
                    auto typ = flag.getType(bx);

                    if (typ == FabType::covered) {
                        AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                        {
                            garr(i,j,k,n) = 0.0;
                        });
                    } else if (typ == FabType::regular) {
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                        {
                            mlndlap_set_surface_integral(i,j,k,garr);
                        });
                    } else {
                        Array4<EBCellFlag const> const& flagarr = flags.const_array(mfi);
                        Array4<RT const> const& bcarr = bcent.const_array(mfi);
                        Array4<RT const> const& baarr = barea.const_array(mfi);
                        Array4<RT const> const& bnarr = bnorm.const_array(mfi);
                        AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                        {
                            mlndlap_set_surface_integral_eb(i,j,k,garr,flagarr,bcarr,baarr,bnarr);
                        });
                    }
                }
            }
        }
#else
        for (int amrlev = 0; amrlev < this->m_num_amr_levels; ++amrlev)
        {
            if (m_surface_integral[amrlev] &&
                dynamic_cast<EBFArrayBoxFactory const*>(this->m_factory[amrlev][0].get())) {
                amrex::algoim::compute_surface_integrals(*m_surface_integral[amrlev]);
            }
        }
#endif
    }
    // The single precision copy has no EB factory, so it has no integrals.
}

template void MLNodeLaplacianT<MultiFab>::buildIntegral ();
template void MLNodeLaplacianT<MultiFab>::buildSurfaceIntegral ();
#ifdef AMREX_MLMG_HAS_FLOAT
template void MLNodeLaplacianT<fMultiFab>::buildIntegral ();
template void MLNodeLaplacianT<fMultiFab>::buildSurfaceIntegral ();
#endif

}

namespace amrex::detail { void mlndlap_anchor_eb () {} } // see mlndlap_link_parts
