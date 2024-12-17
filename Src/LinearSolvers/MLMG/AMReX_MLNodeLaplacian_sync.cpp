#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLNodeLap_K.H>
#include <AMReX_MultiFabUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#endif

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

namespace amrex {

void
MLNodeLaplacian::compSyncResidualCoarse (MultiFab& sync_resid, const MultiFab& a_phi,
                                         const MultiFab& vold, const MultiFab* rhcc,
                                         const BoxArray& fine_grids, const IntVect& ref_ratio)
{
    BL_PROFILE("MLNodeLaplacian::SyncResCrse()");

#if (AMREX_SPACEDIM == 1)
    amrex::Abort("MLNodeLaplacian::compSyncResidualCoarse: 1D not supported");
#endif

    sync_resid.setVal(0.0);

    const Geometry& geom = m_geom[0][0];
    const DistributionMapping& dmap = m_dmap[0][0];
    const BoxArray& ccba = m_grids[0][0];
    const BoxArray& ndba = amrex::convert(ccba, IntVect::TheNodeVector());
    const BoxArray& ccfba = amrex::convert(fine_grids, IntVect::TheZeroVector());

    const auto lobc = LoBC();
    const auto hibc = HiBC();

    // cell-center, 1: coarse; 0: covered by fine
    const int owner = 1;
    const int nonowner = 0;
    iMultiFab crse_cc_mask = amrex::makeFineMask(ccba, dmap, IntVect(1), ccfba, ref_ratio,
                                                 geom.periodicity(), owner, nonowner);

    const Box& ccdom = geom.Domain();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(crse_cc_mask); mfi.isValid(); ++mfi)
    {
        Array4<int> const& fab = crse_cc_mask.array(mfi);
        mlndlap_fillbc_cc<int>(mfi.validbox(),fab,ccdom,lobc,hibc);
    }

    MultiFab phi(ndba, dmap, 1, 1);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Box& gbx = mfi.growntilebox();
        Array4<Real> const& fab = phi.array(mfi);
        Array4<Real const> const& a_fab = a_phi.const_array(mfi);
        Array4<int const> const& msk = crse_cc_mask.const_array(mfi);
        AMREX_HOST_DEVICE_FOR_3D(gbx, i, j, k,
        {
            if (bx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                fab(i,j,k) = a_fab(i,j,k);
                mlndlap_zero_fine(i,j,k,fab,msk,nonowner);
            } else {
                fab(i,j,k) = 0.0;
            }
        });
    }

    const auto& nddom = amrex::surroundingNodes(ccdom);

    const auto dxinv = geom.InvCellSizeArray();

#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

    const auto& sigma_orig = m_sigma[0][0][0];
    const iMultiFab& dmsk = *m_dirichlet_mask[0][0];

#ifdef AMREX_USE_EB
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[0][0].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    const MultiFab* intg = m_integral[0].get();
    const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
#endif

    bool neumann_doubling = true; // yes even for RAP, because unimposeNeumannBC will be called on rhs

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) { mfi_info.EnableTiling().SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox rhs, u;
#ifdef AMREX_USE_EB
        FArrayBox sten, cn;
#endif
        for (MFIter mfi(sync_resid,mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            auto typ = FabType::regular;
#ifdef AMREX_USE_EB
            if (factory) {
                typ = (*flags)[mfi].getType(amrex::enclosedCells(bx));
            }
#endif
            if (typ != FabType::covered)
            {
                const Box& bxg1 = amrex::grow(bx,1);
                const Box& ccbxg1 = amrex::enclosedCells(bxg1);
                Array4<int const> const& cccmsk = crse_cc_mask.const_array(mfi);

                bool has_fine;
                if (Gpu::inLaunchRegion()) {
                    AMREX_ASSERT(ccbxg1 == crse_cc_mask[mfi].box());
                    has_fine = Reduce::AnyOf(ccbxg1,
                                             [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> bool
                    {
                        return cccmsk(i,j,k) == nonowner;
                    });
                } else {
                    has_fine = mlndlap_any_fine_sync_cells(ccbxg1,cccmsk,nonowner);
                }

                if (has_fine)
                {
                    const Box& ccvbx = amrex::enclosedCells(mfi.validbox());

                    u.resize(ccbxg1, AMREX_SPACEDIM);
                    Elixir ueli = u.elixir();
                    Array4<Real> const& uarr = u.array();

                    Box b = ccbxg1 & ccvbx;
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                    {
                        if (m_lobc[0][idim] == LinOpBCType::inflow)
                        {
                            if (b.smallEnd(idim) == ccdom.smallEnd(idim)) {
                                b.growLo(idim, 1);
                            }
                        }
                        if (m_hibc[0][idim] == LinOpBCType::inflow)
                        {
                            if (b.bigEnd(idim) == ccdom.bigEnd(idim)) {
                                b.growHi(idim, 1);
                            }
                        }
                    }

                    Array4<Real const> const& voarr = vold.const_array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                    {
                        if (b.contains(IntVect(AMREX_D_DECL(i,j,k))) && cccmsk(i,j,k)){
                            AMREX_D_TERM(uarr(i,j,k,0) = voarr(i,j,k,0);,
                                         uarr(i,j,k,1) = voarr(i,j,k,1);,
                                         uarr(i,j,k,2) = voarr(i,j,k,2););
                        } else {
                            AMREX_D_TERM(uarr(i,j,k,0) = 0.0;,
                                         uarr(i,j,k,1) = 0.0;,
                                         uarr(i,j,k,2) = 0.0;);
                        }
                    });

                    rhs.resize(bx);
                    Elixir rhseli = rhs.elixir();
                    Array4<Real> const& rhsarr = rhs.array();
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

#ifdef AMREX_USE_EB
                    if (typ == FabType::singlevalued)
                    {
                        Array4<Real const> const& vfracarr = vfrac->const_array(mfi);
                        Array4<Real const> const& intgarr = intg->const_array(mfi);
                        AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                        {
                            mlndlap_divu_eb(i,j,k,rhsarr,uarr,vfracarr,intgarr,dmskarr,dxinv,nddom,lobc,hibc);
                        });
                    }
                    else
#endif
                    {
#if (AMREX_SPACEDIM == 2)
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                        {
                            mlndlap_divu(i,j,k,rhsarr,uarr,dmskarr,dxinv,nddom,lobc,hibc,is_rz);
                        });
#else
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                        {
                            mlndlap_divu(i,j,k,rhsarr,uarr,dmskarr,dxinv,nddom,lobc,hibc);
                        });
#endif
                    }

                    if (rhcc)
                    {
                        Array4<Real> rhccarr = uarr;
                        Array4<Real const> const& rhccarr_orig = rhcc->const_array(mfi);
                        const Box& b2 = ccbxg1 & ccvbx;
                        AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                        {
                            if (b2.contains(IntVect(AMREX_D_DECL(i,j,k))) && cccmsk(i,j,k)){
                                rhccarr(i,j,k) = rhccarr_orig(i,j,k);
                            } else {
                                rhccarr(i,j,k) = 0.0;
                            }
                        });

#ifdef AMREX_USE_EB
                        if (typ == FabType::singlevalued)
                        {
                            Array4<Real const> const& vfracarr = vfrac->const_array(mfi);
                            Array4<Real const> const& intgarr = intg->const_array(mfi);
                            AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                            {
                                Real rhs2 = mlndlap_rhcc_eb(i,j,k,rhccarr,vfracarr,intgarr,dmskarr);
                                rhsarr(i,j,k) += rhs2;
                            });
                        }
                        else
#endif
                        {
                            AMREX_HOST_DEVICE_FOR_3D (bx, i, j, k,
                            {
                                Real rhs2 = mlndlap_rhcc(i,j,k,rhccarr,dmskarr);
                                rhsarr(i,j,k) += rhs2;
                            });
                        }
                    }

                    Array4<Real> const& sync_resid_a = sync_resid.array(mfi);
                    Array4<Real const> const& phiarr = phi.const_array(mfi);
#ifdef AMREX_USE_EB
                    if (typ == FabType::singlevalued)
                    {
                        Array4<Real const> const& sigmaarr_orig = sigma_orig->const_array(mfi);

                        Box stbx = bx;
                        AMREX_D_TERM(stbx.growLo(0,1);, stbx.growLo(1,1);, stbx.growLo(2,1));
                        Box const& sgbx = amrex::grow(amrex::enclosedCells(stbx),1);

                        constexpr int ncomp_s = (AMREX_SPACEDIM == 2) ? 5 : 9;
                        sten.resize(stbx,ncomp_s);
                        Elixir steneli = sten.elixir();
                        Array4<Real> const& stenarr = sten.array();

                        constexpr int ncomp_c = (AMREX_SPACEDIM == 2) ? 6 : 27;
                        cn.resize(sgbx,ncomp_c+1);
                        Elixir cneli = cn.elixir();
                        Array4<Real> const& cnarr = cn.array();
                        Array4<Real> const& sgarr = cn.array(ncomp_c);

                        Array4<EBCellFlag const> const& flagarr = flags->const_array(mfi);
                        Array4<Real const> const& vfracarr = vfrac->const_array(mfi);
                        Array4<Real const> const& intgarr = intg->const_array(mfi);

                        const Box& ibx = sgbx & amrex::enclosedCells(mfi.validbox());
                        AMREX_HOST_DEVICE_FOR_3D(sgbx, i, j, k,
                        {
                            if (ibx.contains(IntVect(AMREX_D_DECL(i,j,k))) && cccmsk(i,j,k)) {
                                mlndlap_set_connection(i,j,k,cnarr,intgarr,vfracarr,flagarr);
                                sgarr(i,j,k) = sigmaarr_orig(i,j,k);
                            } else {
                                for (int n = 0; n < ncomp_c; ++n) {
                                    cnarr(i,j,k,n) = 0.0;
                                }
                                sgarr(i,j,k) = 0.0;
                            }
                        });

                        AMREX_HOST_DEVICE_FOR_3D(stbx, i, j, k,
                        {
                            mlndlap_set_stencil_eb(i, j, k, stenarr, sgarr, cnarr, dxinv);
                        });

                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                        {
                            mlndlap_set_stencil_s0(i, j, k, stenarr);
                            sync_resid_a(i,j,k) = mlndlap_adotx_sten(i, j, k, phiarr, stenarr, dmskarr);
                            mlndlap_crse_resid(i, j, k, sync_resid_a, rhsarr, cccmsk, nddom, lobc, hibc, neumann_doubling);
                        });
                    }
                    else
#endif
                    {
                        Array4<Real> sigmaarr = uarr;
                        const Box& ibx = ccbxg1 & amrex::enclosedCells(mfi.validbox());
                        if (sigma_orig) {
                            Array4<Real const> const& sigmaarr_orig = sigma_orig->const_array(mfi);
                            AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                            {
                                if (ibx.contains(IntVect(AMREX_D_DECL(i,j,k))) && cccmsk(i,j,k)) {
                                    sigmaarr(i,j,k) = sigmaarr_orig(i,j,k);
                                } else {
                                    sigmaarr(i,j,k) = 0.0;
                                }
                            });
                        } else {
                            Real const_sigma = m_const_sigma;
                            AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                            {
                                if (ibx.contains(IntVect(AMREX_D_DECL(i,j,k))) && cccmsk(i,j,k)) {
                                    sigmaarr(i,j,k) = const_sigma;
                                } else {
                                    sigmaarr(i,j,k) = 0.0;
                                }
                            });
                        }

#if (AMREX_SPACEDIM == 2)
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                        {
                            sync_resid_a(i,j,k) = mlndlap_adotx_aa(i, j, k, phiarr, sigmaarr, dmskarr, is_rz, dxinv);
                            mlndlap_crse_resid(i, j, k, sync_resid_a, rhsarr, cccmsk, nddom, lobc, hibc, neumann_doubling);
                        });
#else
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                        {
                            sync_resid_a(i,j,k) = mlndlap_adotx_aa(i, j, k, phiarr, sigmaarr, dmskarr, dxinv);
                            mlndlap_crse_resid(i, j, k, sync_resid_a, rhsarr, cccmsk, nddom, lobc, hibc, neumann_doubling);
                        });
#endif
                    }
                }
            }
        }
    }
}

void
MLNodeLaplacian::compSyncResidualFine (MultiFab& sync_resid, const MultiFab& phi, const MultiFab& vold,
                                       const MultiFab* rhcc)
{
    BL_PROFILE("MLNodeLaplacian::SyncResFine()");

#if (AMREX_SPACEDIM == 1)
    amrex::Abort("MLNodeLaplacian::compSyncResidualFine: 1D not supported");
#endif

    const auto& sigma_orig = m_sigma[0][0][0];
    const iMultiFab& dmsk = *m_dirichlet_mask[0][0];

    const auto lobc = LoBC();
    const auto hibc = HiBC();

#ifdef AMREX_USE_EB
    const auto *factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[0][0].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    const MultiFab* intg = m_integral[0].get();
    const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
#endif

    const Geometry& geom = m_geom[0][0];
    const Box& ccdom = geom.Domain();
    const Box& nddom = amrex::surroundingNodes(ccdom);
    const auto dxinv = geom.InvCellSizeArray();
#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) { mfi_info.EnableTiling().SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox rhs, u;
        IArrayBox tmpmask;
#ifdef AMREX_USE_EB
        FArrayBox sten, cn;
#endif
        for (MFIter mfi(sync_resid,mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            auto typ = FabType::regular;
#ifdef AMREX_USE_EB
            if (factory) {
                typ = (*flags)[mfi].getType(amrex::enclosedCells(bx));
            }
#endif
            if (typ != FabType::covered)
            {
                const Box& gbx = mfi.growntilebox();
                const Box& vbx = mfi.validbox();
                const Box& ccvbx = amrex::enclosedCells(vbx);
                const Box& bxg1 = amrex::grow(bx,1);
                const Box& ccbxg1 = amrex::enclosedCells(bxg1);

                u.resize(ccbxg1, AMREX_SPACEDIM);
                Elixir ueli = u.elixir();
                Array4<Real> const& uarr = u.array();

                Box ovlp = ccvbx & ccbxg1;
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                {
                    if (m_lobc[0][idim] == LinOpBCType::inflow)
                    {
                        if (ovlp.smallEnd(idim) == ccdom.smallEnd(idim)) {
                            ovlp.growLo(idim, 1);
                        }
                    }
                    if (m_hibc[0][idim] == LinOpBCType::inflow)
                    {
                        if (ovlp.bigEnd(idim) == ccdom.bigEnd(idim)) {
                            ovlp.growHi(idim, 1);
                        }
                    }
                }

                Array4<Real const> const& voarr = vold.const_array(mfi);
                AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                {
                    if (ovlp.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                        AMREX_D_TERM(uarr(i,j,k,0) = voarr(i,j,k,0);,
                                     uarr(i,j,k,1) = voarr(i,j,k,1);,
                                     uarr(i,j,k,2) = voarr(i,j,k,2););
                    } else {
                        AMREX_D_TERM(uarr(i,j,k,0) = 0.0;,
                                     uarr(i,j,k,1) = 0.0;,
                                     uarr(i,j,k,2) = 0.0;);
                    }
                });

                tmpmask.resize(bx);
                Elixir tmeli = tmpmask.elixir();
                Array4<int> const& tmpmaskarr = tmpmask.array();
                Array4<int const> const& dmskarr = dmsk.const_array(mfi);
                AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                {
                    tmpmaskarr(i,j,k) = 1-dmskarr(i,j,k);
                });

                rhs.resize(bx);
                Elixir rhseli = rhs.elixir();
                Array4<Real> const& rhsarr = rhs.array();

#ifdef AMREX_USE_EB
                if (typ == FabType::singlevalued)
                {
                    Array4<Real const> const& vfracarr = vfrac->const_array(mfi);
                    Array4<Real const> const& intgarr = intg->const_array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                    {
                        mlndlap_divu_eb(i,j,k,rhsarr,uarr,vfracarr,intgarr,tmpmaskarr,dxinv,nddom,lobc,hibc);
                    });
                }
                else
#endif
                {
#if (AMREX_SPACEDIM == 2)
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                    {
                        mlndlap_divu(i,j,k,rhsarr,uarr,tmpmaskarr,dxinv,nddom,lobc,hibc,is_rz);
                    });
#else
                    AMREX_HOST_DEVICE_PARALLEL_FOR_3D (bx, i, j, k,
                    {
                        mlndlap_divu(i,j,k,rhsarr,uarr,tmpmaskarr,dxinv,nddom,lobc,hibc);
                    });
#endif
                }

                if (rhcc)
                {
                    Array4<Real> rhccarr = uarr;
                    Array4<Real const> const& rhccarr_orig = rhcc->const_array(mfi);
                    const Box& ovlp3 = ccvbx & ccbxg1;
                    AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                    {
                        if (ovlp3.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                            rhccarr(i,j,k) = rhccarr_orig(i,j,k);
                        } else {
                            rhccarr(i,j,k) = 0.0;
                        }
                    });
#ifdef AMREX_USE_EB
                    if (typ == FabType::singlevalued)
                    {
                        Array4<Real const> const& vfracarr = vfrac->const_array(mfi);
                        Array4<Real const> const& intgarr = intg->const_array(mfi);
                        AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                        {
                            Real rhs2 = mlndlap_rhcc_eb(i,j,k,rhccarr,vfracarr,intgarr,tmpmaskarr);
                            rhsarr(i,j,k) += rhs2;
                        });
                    }
                    else
#endif
                    {
                        AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                        {
                            Real rhs2 = mlndlap_rhcc(i,j,k,rhccarr,tmpmaskarr);
                            rhsarr(i,j,k) += rhs2;
                        });
                    }
                }

                Array4<Real> const& sync_resid_a = sync_resid.array(mfi);
                Array4<Real const> const& phiarr = phi.const_array(mfi);
#ifdef AMREX_USE_EB
                if (typ == FabType::singlevalued)
                {
                    Array4<Real const> const& sigmaarr_orig = sigma_orig->const_array(mfi);

                    Box stbx = bx;
                    AMREX_D_TERM(stbx.growLo(0,1);, stbx.growLo(1,1);, stbx.growLo(2,1));
                    Box const& sgbx = amrex::grow(amrex::enclosedCells(stbx),1);

                    constexpr int ncomp_s = (AMREX_SPACEDIM == 2) ? 5 : 9;
                    sten.resize(stbx,ncomp_s);
                    Elixir steneli = sten.elixir();
                    Array4<Real> const& stenarr = sten.array();

                    constexpr int ncomp_c = (AMREX_SPACEDIM == 2) ? 6 : 27;
                    cn.resize(sgbx,ncomp_c+1);
                    Elixir cneli = cn.elixir();
                    Array4<Real> const& cnarr = cn.array();
                    Array4<Real> const& sgarr = cn.array(ncomp_c);

                    Array4<EBCellFlag const> const& flagarr = flags->const_array(mfi);
                    Array4<Real const> const& vfracarr = vfrac->const_array(mfi);
                    Array4<Real const> const& intgarr = intg->const_array(mfi);

                    const Box& ibx = sgbx & amrex::enclosedCells(mfi.validbox());
                    AMREX_HOST_DEVICE_FOR_3D(sgbx, i, j, k,
                    {
                        if (ibx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                            mlndlap_set_connection(i,j,k,cnarr,intgarr,vfracarr,flagarr);
                            sgarr(i,j,k) = sigmaarr_orig(i,j,k);
                        } else {
                            for (int n = 0; n < ncomp_c; ++n) {
                                cnarr(i,j,k,n) = 0.0;
                            }
                            sgarr(i,j,k) = 0.0;
                        }
                    });

                    AMREX_HOST_DEVICE_FOR_3D(stbx, i, j, k,
                    {
                        mlndlap_set_stencil_eb(i, j, k, stenarr, sgarr, cnarr, dxinv);
                    });

                    AMREX_HOST_DEVICE_FOR_3D(gbx, i, j, k,
                    {
                        if (bx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                            mlndlap_set_stencil_s0(i, j, k, stenarr);
                            sync_resid_a(i,j,k) = mlndlap_adotx_sten(i,j,k, phiarr, stenarr, tmpmaskarr);
                            sync_resid_a(i,j,k) = rhsarr(i,j,k) - sync_resid_a(i,j,k);
                        } else {
                            sync_resid_a(i,j,k) = 0.0;
                        }
                    });
                }
                else
#endif
                {
                    Array4<Real> sigmaarr = uarr;
                    const Box& ovlp2 = ccvbx & ccbxg1;
                    if (sigma_orig) {
                        Array4<Real const> const& sigmaarr_orig = sigma_orig->const_array(mfi);
                        AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                        {
                            if (ovlp2.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                                sigmaarr(i,j,k) = sigmaarr_orig(i,j,k);
                            } else {
                                sigmaarr(i,j,k) = 0.0;
                            }
                        });
                    } else {
                        Real const_sigma = m_const_sigma;
                        AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                        {
                            if (ovlp2.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                                sigmaarr(i,j,k) = const_sigma;
                            } else {
                                sigmaarr(i,j,k) = 0.0;
                            }
                        });
                    }

#if (AMREX_SPACEDIM == 2)
                    AMREX_HOST_DEVICE_FOR_3D(gbx, i, j, k,
                    {
                        if (bx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                            sync_resid_a(i,j,k) = mlndlap_adotx_aa(i,j,k, phiarr, sigmaarr, tmpmaskarr,
                                                                   is_rz, dxinv);
                            sync_resid_a(i,j,k) = rhsarr(i,j,k) - sync_resid_a(i,j,k);
                        } else {
                            sync_resid_a(i,j,k) = 0.0;
                        }
                    });
#else
                    AMREX_HOST_DEVICE_FOR_3D(gbx, i, j, k,
                    {
                        if (bx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                            sync_resid_a(i,j,k) = mlndlap_adotx_aa(i,j,k, phiarr, sigmaarr, tmpmaskarr,
                                                                   dxinv);
                            sync_resid_a(i,j,k) = rhsarr(i,j,k) - sync_resid_a(i,j,k);
                        } else {
                            sync_resid_a(i,j,k) = 0.0;
                        }
                    });
#endif
                }
            }

            // Do not impose neumann bc here because how SyncRegister works.
        }
    }
}

void
MLNodeLaplacian::reflux (int crse_amrlev,
                         MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                         MultiFab& a_fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(crse_amrlev,res,crse_sol,crse_rhs,a_fine_res,fine_sol,fine_rhs);
    amrex::Abort("MLNodeLaplacian::reflux: 1D not supported");
#else
    //
    //  Note that the residue we copmute on a coarse/fine node is not a
    //  composite divergence.  It has been restricted so that it is suitable
    //  as RHS for our geometric multigrid solver with a MG hirerachy
    //  including multiple AMR levels.
    //

    BL_PROFILE("MLNodeLaplacian::reflux()");

    const int amrrr = AMRRefRatio(crse_amrlev);
    AMREX_ALWAYS_ASSERT(amrrr == 2 || m_coarsening_strategy == CoarseningStrategy::Sigma);

    const Geometry& cgeom = m_geom[crse_amrlev  ][0];
    const Geometry& fgeom = m_geom[crse_amrlev+1][0];
    const auto cdxinv = cgeom.InvCellSizeArray();
    const auto fdxinv = fgeom.InvCellSizeArray();
    const Box& c_cc_domain = cgeom.Domain();
    const Box& c_cc_domain_p = cgeom.growPeriodicDomain(1);
    const Box& c_nd_domain = amrex::surroundingNodes(c_cc_domain);
    const Box& f_nd_domain = amrex::surroundingNodes(fgeom.Domain());

    const auto lobc = LoBC();
    const auto hibc = HiBC();

    bool neumann_doubling = false;
    if (m_coarsening_strategy == CoarseningStrategy::Sigma) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            neumann_doubling = neumann_doubling || (lobc[idim] == LinOpBCType::inflow  ||
                                                    lobc[idim] == LinOpBCType::Neumann ||
                                                    hibc[idim] == LinOpBCType::inflow  ||
                                                    hibc[idim] == LinOpBCType::Neumann);
        }
    }

#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

    const BoxArray& fba = fine_sol.boxArray();
    const DistributionMapping& fdm = fine_sol.DistributionMap();

    const iMultiFab& fdmsk = *m_dirichlet_mask[crse_amrlev+1][0];
    const auto& stencil    =  m_nosigma_stencil[crse_amrlev+1];

    MultiFab fine_res_for_coarse(amrex::coarsen(fba, amrrr), fdm, 1, 0);

    std::unique_ptr<MultiFab> tmp_fine_res;
    if (!a_fine_res.nGrowVect().allGE(amrrr-1)) {
        tmp_fine_res = std::make_unique<MultiFab>(a_fine_res.boxArray(),
                                                  a_fine_res.DistributionMap(), 1, amrrr-1);
        MultiFab::Copy(*tmp_fine_res, a_fine_res, 0, 0, 1, 0);
    }
    MultiFab& fine_res = (tmp_fine_res) ? *tmp_fine_res :  a_fine_res;

    fine_res.FillBoundary(fgeom.periodicity());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine_res_for_coarse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& cfab = fine_res_for_coarse.array(mfi);
        Array4<Real const> const& ffab = fine_res.const_array(mfi);
        Array4<int const> const& mfab = fdmsk.const_array(mfi);
        if (m_coarsening_strategy == CoarseningStrategy::Sigma) {
            if (amrrr == 2) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_restriction<2>(i,j,k,cfab,ffab,mfab,f_nd_domain,lobc,hibc);
                });
            } else {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_restriction<4>(i,j,k,cfab,ffab,mfab,f_nd_domain,lobc,hibc);
                });
            }
        } else {
            Array4<Real const> const& stfab = stencil->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_restriction_rap(i,j,k,cfab,ffab,stfab,mfab);
            });
        }
    }
    res.ParallelCopy(fine_res_for_coarse, cgeom.periodicity());

    MultiFab fine_contrib(amrex::coarsen(fba, amrrr), fdm, 1, 0);

    const auto& fsigma = m_sigma[crse_amrlev+1][0][0];

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) { mfi_info.EnableTiling().SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine_contrib,mfi_info); mfi.isValid(); ++mfi)
    {
        const Box& cbx = mfi.tilebox();
        const Box& fvbx = amrex::refine(mfi.validbox(),amrrr);
        const Box& cc_fvbx = amrex::enclosedCells(fvbx);

        Array4<Real> const& farr = fine_contrib.array(mfi);
        Array4<Real const> const& resarr = fine_res.const_array(mfi);
        Array4<Real const> const& rhsarr = fine_rhs.const_array(mfi);
        Array4<Real const> const& solarr = fine_sol.const_array(mfi);
        Array4<int const> const& marr = fdmsk.const_array(mfi);

        if (fsigma) {
            Array4<Real const> const& sigarr = fsigma->const_array(mfi);
#if (AMREX_SPACEDIM == 2)
            if (amrrr == 2) {
                AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                {
                    mlndlap_Ax_fine_contrib<2>(i,j,k,fvbx,cc_fvbx,farr,resarr,rhsarr,solarr,
                                               sigarr,marr,is_rz,fdxinv);
                });
            } else {
                AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                {
                    mlndlap_Ax_fine_contrib<4>(i,j,k,fvbx,cc_fvbx,farr,resarr,rhsarr,solarr,
                                               sigarr,marr,is_rz,fdxinv);
                });
            }
#elif (AMREX_SPACEDIM == 3)
            if (amrrr == 2) {
                AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                {
                    mlndlap_Ax_fine_contrib<2>(i,j,k,fvbx,cc_fvbx,farr,resarr,rhsarr,solarr,
                                               sigarr,marr,fdxinv);
                });
            } else {
                AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                {
                    mlndlap_Ax_fine_contrib<4>(i,j,k,fvbx,cc_fvbx,farr,resarr,rhsarr,solarr,
                                               sigarr,marr,fdxinv);
                });
            }
#endif
        } else {
            Real const_sigma = m_const_sigma;
#if (AMREX_SPACEDIM == 2)
            if (amrrr == 2) {
                AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                {
                    mlndlap_Ax_fine_contrib_cs<2>(i,j,k,fvbx,cc_fvbx,farr,resarr,rhsarr,solarr,
                                                  const_sigma,marr,is_rz,fdxinv);
                });
            } else {
                AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                {
                    mlndlap_Ax_fine_contrib_cs<4>(i,j,k,fvbx,cc_fvbx,farr,resarr,rhsarr,solarr,
                                                  const_sigma,marr,is_rz,fdxinv);
                });
            }
#elif (AMREX_SPACEDIM == 3)
            if (amrrr == 2) {
                AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                {
                    mlndlap_Ax_fine_contrib_cs<2>(i,j,k,fvbx,cc_fvbx,farr,resarr,rhsarr,solarr,
                                                  const_sigma,marr,fdxinv);
                });
            } else {
                AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                {
                    mlndlap_Ax_fine_contrib_cs<4>(i,j,k,fvbx,cc_fvbx,farr,resarr,rhsarr,solarr,
                                                  const_sigma,marr,fdxinv);
                });
            }
#endif
        }
    }

    MultiFab fine_contrib_on_crse(crse_sol.boxArray(), crse_sol.DistributionMap(), 1, 0);
    fine_contrib_on_crse.setVal(0.0);
    fine_contrib_on_crse.ParallelAdd(fine_contrib, cgeom.periodicity());

    const iMultiFab& cdmsk = *m_dirichlet_mask[crse_amrlev][0];
    const auto& nd_mask     = m_nd_fine_mask[crse_amrlev];
    const auto& cc_mask     = m_cc_fine_mask[crse_amrlev];
    const auto& has_fine_bndry = m_has_fine_bndry[crse_amrlev];

    const auto& csigma = m_sigma[crse_amrlev][0][0];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(res,mfi_info); mfi.isValid(); ++mfi)
    {
        if ((*has_fine_bndry)[mfi])
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& resarr = res.array(mfi);
            Array4<Real const> const& csolarr = crse_sol.const_array(mfi);
            Array4<Real const> const& crhsarr = crse_rhs.const_array(mfi);
            Array4<int const> const& cdmskarr = cdmsk.const_array(mfi);
            Array4<int const> const& ndmskarr = nd_mask->const_array(mfi);
            Array4<int const> const& ccmskarr = cc_mask->const_array(mfi);
            Array4<Real const> const& fcocarr = fine_contrib_on_crse.const_array(mfi);

            if (csigma) {
                Array4<Real const> const& csigarr = csigma->const_array(mfi);
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                {
                    mlndlap_res_cf_contrib(i,j,k,resarr,csolarr,crhsarr,csigarr,
                                           cdmskarr,ndmskarr,ccmskarr,fcocarr,
                                           cdxinv,c_cc_domain_p,c_nd_domain,
                                           is_rz,
                                           lobc,hibc, neumann_doubling);
                });
#elif (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                {
                    mlndlap_res_cf_contrib(i,j,k,resarr,csolarr,crhsarr,csigarr,
                                           cdmskarr,ndmskarr,ccmskarr,fcocarr,
                                           cdxinv,c_cc_domain_p,c_nd_domain,
                                           lobc,hibc, neumann_doubling);
                });
#endif
            } else {
                Real const_sigma = m_const_sigma;
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                {
                    mlndlap_res_cf_contrib_cs(i,j,k,resarr,csolarr,crhsarr,const_sigma,
                                              cdmskarr,ndmskarr,ccmskarr,fcocarr,
                                              cdxinv,c_cc_domain_p,c_nd_domain,
                                              is_rz,
                                              lobc,hibc, neumann_doubling);
                });
#elif (AMREX_SPACEDIM == 3)
                AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                {
                    mlndlap_res_cf_contrib_cs(i,j,k,resarr,csolarr,crhsarr,const_sigma,
                                              cdmskarr,ndmskarr,ccmskarr,fcocarr,
                                              cdxinv,c_cc_domain_p,c_nd_domain,
                                              lobc,hibc, neumann_doubling);
                });
#endif
            }
        }
    }
#ifdef AMREX_USE_EB
    // Make sure to zero out the residual on any nodes completely surrounded by covered cells
    amrex::EB_set_covered(res,0.0);
#endif
#endif
}

}
