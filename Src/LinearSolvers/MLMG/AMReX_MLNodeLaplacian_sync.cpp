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

template <typename MF>
void
MLNodeLaplacianT<MF>::compSyncResidualCoarse (MF& sync_resid, const MF& a_phi,
                                              const MF& vold, const MF* rhcc,
                                              const BoxArray& fine_grids, const IntVect& ref_ratio)
{
    BL_PROFILE("MLNodeLaplacian::SyncResCrse()");

#if (AMREX_SPACEDIM == 1)
    amrex::Abort("MLNodeLaplacian::compSyncResidualCoarse: 1D not supported");
#endif

    sync_resid.setVal(0.0);

    const Geometry& geom = this->m_geom[0][0];
    const DistributionMapping& dmap = this->m_dmap[0][0];
    const BoxArray& ccba = this->m_grids[0][0];
    const BoxArray& ndba = amrex::convert(ccba, IntVect::TheNodeVector());
    const BoxArray& ccfba = amrex::convert(fine_grids, IntVect::TheZeroVector());

    const auto lobc = this->LoBC();
    const auto hibc = this->HiBC();

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

    MF phi(ndba, dmap, 1, 1);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(phi,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Box& gbx = mfi.growntilebox();
        Array4<RT> const& fab = phi.array(mfi);
        Array4<RT const> const& a_fab = a_phi.const_array(mfi);
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

    const auto dxinv = this->toRT(geom.InvCellSizeArray());

#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

    const auto& sigma_orig = m_sigma[0][0][0];
    // Anisotropic (mapped) path: capture per-axis sigmas.
    const bool aniso_sigma = m_use_mapped && (sigma_orig != nullptr)
        AMREX_D_TERM(, && (m_sigma[0][0][1] != nullptr),
                     && (m_sigma[0][0][2] != nullptr));
    AMREX_D_TERM(,
                 const MF* sigma_orig_y = aniso_sigma ? m_sigma[0][0][1].get() : nullptr;,
                 const MF* sigma_orig_z = aniso_sigma ? m_sigma[0][0][2].get() : nullptr;)
    const iMultiFab& dmsk = *this->m_dirichlet_mask[0][0];

#ifdef AMREX_USE_EB
    const auto *factory = ebFactory(0);
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    const MF* intg = m_integral[0].get();
    const MF* vfrac = ebVolFrac(factory);
#endif

    bool neumann_doubling = true; // yes even for RAP, because unimposeNeumannBC will be called on rhs

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) { mfi_info.EnableTiling().SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FAB rhs, u;
        // Per-axis sigma scratch for the anisotropic path.
        FAB sig_y_local;
#if (AMREX_SPACEDIM == 3)
        FAB sig_z_local;
#endif
#ifdef AMREX_USE_EB
        FAB sten, cn;
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

                    u.resize(ccbxg1, AMREX_SPACEDIM, The_Async_Arena());
                    Array4<RT> const& uarr = u.array();

                    Box b = ccbxg1 & ccvbx;
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                    {
                        if (this->m_lobc[0][idim] == LinOpBCType::inflow)
                        {
                            if (b.smallEnd(idim) == ccdom.smallEnd(idim)) {
                                b.growLo(idim, 1);
                            }
                        }
                        if (this->m_hibc[0][idim] == LinOpBCType::inflow)
                        {
                            if (b.bigEnd(idim) == ccdom.bigEnd(idim)) {
                                b.growHi(idim, 1);
                            }
                        }
                    }

                    Array4<RT const> const& voarr = vold.const_array(mfi);
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

                    rhs.resize(bx, 1, The_Async_Arena());
                    Array4<RT> const& rhsarr = rhs.array();
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

#ifdef AMREX_USE_EB
                    if (typ == FabType::singlevalued)
                    {
                        Array4<RT const> const& vfracarr = vfrac->const_array(mfi);
                        Array4<RT const> const& intgarr = intg->const_array(mfi);
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
                        Array4<RT> rhccarr = uarr;
                        Array4<RT const> const& rhccarr_orig = rhcc->const_array(mfi);
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
                            Array4<RT const> const& vfracarr = vfrac->const_array(mfi);
                            Array4<RT const> const& intgarr = intg->const_array(mfi);
                            AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                            {
                                RT rhs2 = mlndlap_rhcc_eb<RT>(i,j,k,rhccarr,vfracarr,intgarr,dmskarr);
                                rhsarr(i,j,k) += rhs2;
                            });
                        }
                        else
#endif
                        {
                            AMREX_HOST_DEVICE_FOR_3D (bx, i, j, k,
                            {
                                RT rhs2 = mlndlap_rhcc<RT>(i,j,k,rhccarr,dmskarr);
                                rhsarr(i,j,k) += rhs2;
                            });
                        }
                    }

                    Array4<RT> const& sync_resid_a = sync_resid.array(mfi);
                    Array4<RT const> const& phiarr = phi.const_array(mfi);
#ifdef AMREX_USE_EB
                    if (typ == FabType::singlevalued)
                    {
                        Array4<RT const> const& sigmaarr_orig = sigma_orig->const_array(mfi);

                        Box stbx = bx;
                        AMREX_D_TERM(stbx.growLo(0,1);, stbx.growLo(1,1);, stbx.growLo(2,1));
                        Box const& sgbx = amrex::grow(amrex::enclosedCells(stbx),1);

                        constexpr int ncomp_s = (AMREX_SPACEDIM == 2) ? 5 : 9;
                        sten.resize(stbx,ncomp_s, The_Async_Arena());
                        Array4<RT> const& stenarr = sten.array();

                        constexpr int ncomp_c = (AMREX_SPACEDIM == 2) ? 6 : 27;
                        cn.resize(sgbx,ncomp_c+1, The_Async_Arena());
                        Array4<RT> const& cnarr = cn.array();
                        Array4<RT> const& sgarr = cn.array(ncomp_c);

                        Array4<EBCellFlag const> const& flagarr = flags->const_array(mfi);
                        Array4<RT const> const& vfracarr = vfrac->const_array(mfi);
                        Array4<RT const> const& intgarr = intg->const_array(mfi);

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
                        Array4<RT> sigmaarr = uarr;
                        const Box& ibx = ccbxg1 & amrex::enclosedCells(mfi.validbox());

                        // Per-axis sigma scratch when mapped.
                        if (aniso_sigma) {
                            AMREX_D_TERM(,
                                         sig_y_local.resize(ccbxg1, 1, The_Async_Arena());,
                                         sig_z_local.resize(ccbxg1, 1, The_Async_Arena());)
                        }
                        AMREX_D_TERM(,
                                     Array4<RT> syarr_w = aniso_sigma ? sig_y_local.array() : Array4<RT>{};,
                                     Array4<RT> szarr_w = aniso_sigma ? sig_z_local.array() : Array4<RT>{};)

                        if (sigma_orig) {
                            Array4<RT const> const& sigmaarr_orig = sigma_orig->const_array(mfi);
                            if (aniso_sigma) {
                                AMREX_D_TERM(,
                                             Array4<RT const> syarr_orig
                                                 = sigma_orig_y->const_array(mfi);,
                                             Array4<RT const> szarr_orig
                                                 = sigma_orig_z->const_array(mfi);)
                                AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                                {
                                    if (ibx.contains(IntVect(AMREX_D_DECL(i,j,k))) && cccmsk(i,j,k)) {
                                        AMREX_D_TERM(sigmaarr(i,j,k) = sigmaarr_orig(i,j,k);,
                                                     syarr_w(i,j,k) = syarr_orig(i,j,k);,
                                                     szarr_w(i,j,k) = szarr_orig(i,j,k);)
                                    } else {
                                        AMREX_D_TERM(sigmaarr(i,j,k) = 0.0;,
                                                     syarr_w(i,j,k) = 0.0;,
                                                     szarr_w(i,j,k) = 0.0;)
                                    }
                                });
                            } else {
                                AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                                {
                                    if (ibx.contains(IntVect(AMREX_D_DECL(i,j,k))) && cccmsk(i,j,k)) {
                                        sigmaarr(i,j,k) = sigmaarr_orig(i,j,k);
                                    } else {
                                        sigmaarr(i,j,k) = 0.0;
                                    }
                                });
                            }
                        } else {
                            // Const sigma: anisotropic n/a here.
                            RT const_sigma = m_const_sigma;
                            AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                            {
                                if (ibx.contains(IntVect(AMREX_D_DECL(i,j,k))) && cccmsk(i,j,k)) {
                                    sigmaarr(i,j,k) = const_sigma;
                                } else {
                                    sigmaarr(i,j,k) = 0.0;
                                }
                            });
                        }

                        if (aniso_sigma) {
                            Array4<RT const> sxarr_c = sigmaarr;
#if (AMREX_SPACEDIM >= 2)
                            Array4<RT const> syarr_c = syarr_w;
#endif
#if (AMREX_SPACEDIM == 3)
                            Array4<RT const> szarr_c = szarr_w;
#endif
#if (AMREX_SPACEDIM == 1)
                            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                            {
                                sync_resid_a(i,j,k) = mlndlap_adotx_ha(i,j,k, phiarr,
                                                                       sxarr_c, dmskarr, dxinv);
                                mlndlap_crse_resid(i,j,k, sync_resid_a, rhsarr, cccmsk, nddom, lobc, hibc, neumann_doubling);
                            });
#elif (AMREX_SPACEDIM == 2)
                            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                            {
                                sync_resid_a(i,j,k) = mlndlap_adotx_ha(i,j,k, phiarr,
                                                                       sxarr_c, syarr_c, dmskarr,
                                                                       is_rz, dxinv);
                                mlndlap_crse_resid(i,j,k, sync_resid_a, rhsarr, cccmsk, nddom, lobc, hibc, neumann_doubling);
                            });
#else
                            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                            {
                                sync_resid_a(i,j,k) = mlndlap_adotx_ha(i,j,k, phiarr,
                                                                       sxarr_c, syarr_c, szarr_c,
                                                                       dmskarr, dxinv);
                                mlndlap_crse_resid(i,j,k, sync_resid_a, rhsarr, cccmsk, nddom, lobc, hibc, neumann_doubling);
                            });
#endif
                        } else {
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
}

template <typename MF>
void
MLNodeLaplacianT<MF>::compSyncResidualFine (MF& sync_resid, const MF& phi, const MF& vold,
                                            const MF* rhcc)
{
    BL_PROFILE("MLNodeLaplacian::SyncResFine()");

#if (AMREX_SPACEDIM == 1)
    amrex::Abort("MLNodeLaplacian::compSyncResidualFine: 1D not supported");
#endif

    const auto& sigma_orig = m_sigma[0][0][0];
    // Per-axis sigma for the anisotropic (mapped) path.
    const bool aniso_sigma = m_use_mapped && (sigma_orig != nullptr)
        AMREX_D_TERM(, && (m_sigma[0][0][1] != nullptr),
                     && (m_sigma[0][0][2] != nullptr));
    AMREX_D_TERM(,
                 const MF* sigma_orig_y = aniso_sigma ? m_sigma[0][0][1].get() : nullptr;,
                 const MF* sigma_orig_z = aniso_sigma ? m_sigma[0][0][2].get() : nullptr;)
    const iMultiFab& dmsk = *this->m_dirichlet_mask[0][0];

    const auto lobc = this->LoBC();
    const auto hibc = this->HiBC();

#ifdef AMREX_USE_EB
    const auto *factory = ebFactory(0);
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    const MF* intg = m_integral[0].get();
    const MF* vfrac = ebVolFrac(factory);
#endif

    const Geometry& geom = this->m_geom[0][0];
    const Box& ccdom = geom.Domain();
    const Box& nddom = amrex::surroundingNodes(ccdom);
    const auto dxinv = this->toRT(geom.InvCellSizeArray());
#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) { mfi_info.EnableTiling().SetDynamic(true); }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FAB rhs, u;
        IArrayBox tmpmask;
        // Per-axis sigma scratch for the anisotropic path.
        FAB sig_y_local;
#if (AMREX_SPACEDIM == 3)
        FAB sig_z_local;
#endif
#ifdef AMREX_USE_EB
        FAB sten, cn;
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

                u.resize(ccbxg1, AMREX_SPACEDIM, The_Async_Arena());
                Array4<RT> const& uarr = u.array();

                Box ovlp = ccvbx & ccbxg1;
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                {
                    if (this->m_lobc[0][idim] == LinOpBCType::inflow)
                    {
                        if (ovlp.smallEnd(idim) == ccdom.smallEnd(idim)) {
                            ovlp.growLo(idim, 1);
                        }
                    }
                    if (this->m_hibc[0][idim] == LinOpBCType::inflow)
                    {
                        if (ovlp.bigEnd(idim) == ccdom.bigEnd(idim)) {
                            ovlp.growHi(idim, 1);
                        }
                    }
                }

                Array4<RT const> const& voarr = vold.const_array(mfi);
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

                tmpmask.resize(bx, 1, The_Async_Arena());
                Array4<int> const& tmpmaskarr = tmpmask.array();
                Array4<int const> const& dmskarr = dmsk.const_array(mfi);
                AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                {
                    tmpmaskarr(i,j,k) = dmskarr(i,j,k) ? 0 : 1;
                });

                rhs.resize(bx, 1, The_Async_Arena());
                Array4<RT> const& rhsarr = rhs.array();

#ifdef AMREX_USE_EB
                if (typ == FabType::singlevalued)
                {
                    Array4<RT const> const& vfracarr = vfrac->const_array(mfi);
                    Array4<RT const> const& intgarr = intg->const_array(mfi);
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
                    Array4<RT> rhccarr = uarr;
                    Array4<RT const> const& rhccarr_orig = rhcc->const_array(mfi);
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
                        Array4<RT const> const& vfracarr = vfrac->const_array(mfi);
                        Array4<RT const> const& intgarr = intg->const_array(mfi);
                        AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                        {
                            RT rhs2 = mlndlap_rhcc_eb<RT>(i,j,k,rhccarr,vfracarr,intgarr,tmpmaskarr);
                            rhsarr(i,j,k) += rhs2;
                        });
                    }
                    else
#endif
                    {
                        AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                        {
                            RT rhs2 = mlndlap_rhcc<RT>(i,j,k,rhccarr,tmpmaskarr);
                            rhsarr(i,j,k) += rhs2;
                        });
                    }
                }

                Array4<RT> const& sync_resid_a = sync_resid.array(mfi);
                Array4<RT const> const& phiarr = phi.const_array(mfi);
#ifdef AMREX_USE_EB
                if (typ == FabType::singlevalued)
                {
                    Array4<RT const> const& sigmaarr_orig = sigma_orig->const_array(mfi);

                    Box stbx = bx;
                    AMREX_D_TERM(stbx.growLo(0,1);, stbx.growLo(1,1);, stbx.growLo(2,1));
                    Box const& sgbx = amrex::grow(amrex::enclosedCells(stbx),1);

                    constexpr int ncomp_s = (AMREX_SPACEDIM == 2) ? 5 : 9;
                    sten.resize(stbx,ncomp_s, The_Async_Arena());
                    Array4<RT> const& stenarr = sten.array();

                    constexpr int ncomp_c = (AMREX_SPACEDIM == 2) ? 6 : 27;
                    cn.resize(sgbx,ncomp_c+1, The_Async_Arena());
                    Array4<RT> const& cnarr = cn.array();
                    Array4<RT> const& sgarr = cn.array(ncomp_c);

                    Array4<EBCellFlag const> const& flagarr = flags->const_array(mfi);
                    Array4<RT const> const& vfracarr = vfrac->const_array(mfi);
                    Array4<RT const> const& intgarr = intg->const_array(mfi);

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
                    Array4<RT> sigmaarr = uarr;
                    const Box& ovlp2 = ccvbx & ccbxg1;

                    if (aniso_sigma) {
                        AMREX_D_TERM(,
                                     sig_y_local.resize(ccbxg1, 1, The_Async_Arena());,
                                     sig_z_local.resize(ccbxg1, 1, The_Async_Arena());)
                    }
                    AMREX_D_TERM(,
                                 Array4<RT> syarr_w = aniso_sigma ? sig_y_local.array() : Array4<RT>{};,
                                 Array4<RT> szarr_w = aniso_sigma ? sig_z_local.array() : Array4<RT>{};)

                    if (sigma_orig) {
                        Array4<RT const> const& sigmaarr_orig = sigma_orig->const_array(mfi);
                        if (aniso_sigma) {
                            AMREX_D_TERM(,
                                         Array4<RT const> syarr_orig
                                             = sigma_orig_y->const_array(mfi);,
                                         Array4<RT const> szarr_orig
                                             = sigma_orig_z->const_array(mfi);)
                            AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                            {
                                if (ovlp2.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                                    AMREX_D_TERM(sigmaarr(i,j,k) = sigmaarr_orig(i,j,k);,
                                                 syarr_w(i,j,k) = syarr_orig(i,j,k);,
                                                 szarr_w(i,j,k) = szarr_orig(i,j,k);)
                                } else {
                                    AMREX_D_TERM(sigmaarr(i,j,k) = 0.0;,
                                                 syarr_w(i,j,k) = 0.0;,
                                                 szarr_w(i,j,k) = 0.0;)
                                }
                            });
                        } else {
                            AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                            {
                                if (ovlp2.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                                    sigmaarr(i,j,k) = sigmaarr_orig(i,j,k);
                                } else {
                                    sigmaarr(i,j,k) = 0.0;
                                }
                            });
                        }
                    } else {
                        RT const_sigma = m_const_sigma;
                        AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                        {
                            if (ovlp2.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                                sigmaarr(i,j,k) = const_sigma;
                            } else {
                                sigmaarr(i,j,k) = 0.0;
                            }
                        });
                    }

                    if (aniso_sigma) {
                        Array4<RT const> sxarr_c = sigmaarr;
#if (AMREX_SPACEDIM >= 2)
                        Array4<RT const> syarr_c = syarr_w;
#endif
#if (AMREX_SPACEDIM == 3)
                        Array4<RT const> szarr_c = szarr_w;
#endif
#if (AMREX_SPACEDIM == 1)
                        AMREX_HOST_DEVICE_FOR_3D(gbx, i, j, k,
                        {
                            if (bx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                                sync_resid_a(i,j,k) = mlndlap_adotx_ha(i,j,k, phiarr,
                                                                       sxarr_c, tmpmaskarr, dxinv);
                                sync_resid_a(i,j,k) = rhsarr(i,j,k) - sync_resid_a(i,j,k);
                            } else {
                                sync_resid_a(i,j,k) = 0.0;
                            }
                        });
#elif (AMREX_SPACEDIM == 2)
                        AMREX_HOST_DEVICE_FOR_3D(gbx, i, j, k,
                        {
                            if (bx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                                sync_resid_a(i,j,k) = mlndlap_adotx_ha(i,j,k, phiarr,
                                                                       sxarr_c, syarr_c, tmpmaskarr,
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
                                sync_resid_a(i,j,k) = mlndlap_adotx_ha(i,j,k, phiarr,
                                                                       sxarr_c, syarr_c, szarr_c,
                                                                       tmpmaskarr, dxinv);
                                sync_resid_a(i,j,k) = rhsarr(i,j,k) - sync_resid_a(i,j,k);
                            } else {
                                sync_resid_a(i,j,k) = 0.0;
                            }
                        });
#endif
                    } else {
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
            }

            // Do not impose neumann bc here because how SyncRegister works.
        }
    }
}

template <typename MF>
void
MLNodeLaplacianT<MF>::reflux (int crse_amrlev,
                              MF& res, const MF& crse_sol, const MF& crse_rhs,
                              MF& a_fine_res, MF& fine_sol, const MF& fine_rhs) const
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

    const int amrrr = this->AMRRefRatio(crse_amrlev);
    AMREX_ALWAYS_ASSERT(amrrr == 2 || this->m_coarsening_strategy == CoarseningStrategy::Sigma);

    const Geometry& cgeom = this->m_geom[crse_amrlev  ][0];
    const Geometry& fgeom = this->m_geom[crse_amrlev+1][0];
    const auto cdxinv = this->toRT(cgeom.InvCellSizeArray());
    const auto fdxinv = this->toRT(fgeom.InvCellSizeArray());
    const Box& c_cc_domain = cgeom.Domain();
    const Box& c_cc_domain_p = cgeom.growPeriodicDomain(1);
    const Box& c_nd_domain = amrex::surroundingNodes(c_cc_domain);
    const Box& f_nd_domain = amrex::surroundingNodes(fgeom.Domain());

    const auto lobc = this->LoBC();
    const auto hibc = this->HiBC();

    bool neumann_doubling = false;
    if (this->m_coarsening_strategy == CoarseningStrategy::Sigma) {
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

    const iMultiFab& fdmsk = *this->m_dirichlet_mask[crse_amrlev+1][0];
    const auto& stencil    =  m_nosigma_stencil[crse_amrlev+1];

    MF fine_res_for_coarse(amrex::coarsen(fba, amrrr), fdm, 1, 0,
                           MFInfo().SetArena(The_Async_Arena()));

    std::unique_ptr<MF> tmp_fine_res;
    if (!a_fine_res.nGrowVect().allGE(amrrr-1)) {
        tmp_fine_res = std::make_unique<MF>(a_fine_res.boxArray(),
                                            a_fine_res.DistributionMap(), 1, amrrr-1,
                                            MFInfo().SetArena(The_Async_Arena()));
        amrex::Copy(*tmp_fine_res, a_fine_res, 0, 0, 1, 0);
    }
    MF& fine_res = (tmp_fine_res) ? *tmp_fine_res :  a_fine_res;

    fine_res.FillBoundary(fgeom.periodicity());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine_res_for_coarse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<RT> const& cfab = fine_res_for_coarse.array(mfi);
        Array4<RT const> const& ffab = fine_res.const_array(mfi);
        Array4<int const> const& mfab = fdmsk.const_array(mfi);
        if (this->m_coarsening_strategy == CoarseningStrategy::Sigma) {
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
            Array4<RT const> const& stfab = stencil->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_restriction_rap(i,j,k,cfab,ffab,stfab,mfab);
            });
        }
    }
    res.ParallelCopy(fine_res_for_coarse, cgeom.periodicity());

    MF fine_contrib(amrex::coarsen(fba, amrrr), fdm, 1, 0, MFInfo().SetArena(The_Async_Arena()));

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

        Array4<RT> const& farr = fine_contrib.array(mfi);
        Array4<RT const> const& resarr = fine_res.const_array(mfi);
        Array4<RT const> const& rhsarr = fine_rhs.const_array(mfi);
        Array4<RT const> const& solarr = fine_sol.const_array(mfi);
        Array4<int const> const& marr = fdmsk.const_array(mfi);

        if (fsigma) {
            Array4<RT const> const& sigarr = fsigma->const_array(mfi);
            const bool aniso_sigma = m_use_mapped && (fsigma != nullptr)
                AMREX_D_TERM(, && (m_sigma[crse_amrlev+1][0][1] != nullptr),
                             && (m_sigma[crse_amrlev+1][0][2] != nullptr));
            AMREX_D_TERM(,
                         Array4<RT const> syarr = aniso_sigma
                             ? m_sigma[crse_amrlev+1][0][1]->const_array(mfi) : Array4<RT const>{};,
                         Array4<RT const> szarr = aniso_sigma
                             ? m_sigma[crse_amrlev+1][0][2]->const_array(mfi) : Array4<RT const>{};)
#if (AMREX_SPACEDIM == 2)
            if (aniso_sigma) {
                if (amrrr == 2) {
                    AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                    {
                        mlndlap_Ax_fine_contrib_ha<2>(i,j,k,fvbx,cc_fvbx,farr,resarr,rhsarr,solarr,
                                                      sigarr,syarr,marr,is_rz,fdxinv);
                    });
                } else {
                    AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                    {
                        mlndlap_Ax_fine_contrib_ha<4>(i,j,k,fvbx,cc_fvbx,farr,resarr,rhsarr,solarr,
                                                      sigarr,syarr,marr,is_rz,fdxinv);
                    });
                }
            } else if (amrrr == 2) {
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
            if (aniso_sigma) {
                if (amrrr == 2) {
                    AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                    {
                        mlndlap_Ax_fine_contrib_ha<2>(i,j,k,fvbx,cc_fvbx,farr,resarr,rhsarr,solarr,
                                                      sigarr,syarr,szarr,marr,fdxinv);
                    });
                } else {
                    AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                    {
                        mlndlap_Ax_fine_contrib_ha<4>(i,j,k,fvbx,cc_fvbx,farr,resarr,rhsarr,solarr,
                                                      sigarr,syarr,szarr,marr,fdxinv);
                    });
                }
            } else if (amrrr == 2) {
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
            RT const_sigma = m_const_sigma;
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

    MF fine_contrib_on_crse(crse_sol.boxArray(), crse_sol.DistributionMap(), 1, 0,
                            MFInfo().SetArena(The_Async_Arena()));
    fine_contrib_on_crse.setVal(0.0);
    fine_contrib_on_crse.ParallelAdd(fine_contrib, cgeom.periodicity());

    const iMultiFab& cdmsk = *this->m_dirichlet_mask[crse_amrlev][0];
    const auto& nd_mask     = this->m_nd_fine_mask[crse_amrlev];
    const auto& cc_mask     = this->m_cc_fine_mask[crse_amrlev];
    const auto& has_fine_bndry = this->m_has_fine_bndry[crse_amrlev];

    const auto& csigma = m_sigma[crse_amrlev][0][0];
    const bool aniso_csigma = m_use_mapped && (csigma != nullptr)
        AMREX_D_TERM(, && (m_sigma[crse_amrlev][0][1] != nullptr),
                     && (m_sigma[crse_amrlev][0][2] != nullptr));

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(res,mfi_info); mfi.isValid(); ++mfi)
    {
        if ((*has_fine_bndry)[mfi])
        {
            const Box& bx = mfi.tilebox();
            Array4<RT> const& resarr = res.array(mfi);
            Array4<RT const> const& csolarr = crse_sol.const_array(mfi);
            Array4<RT const> const& crhsarr = crse_rhs.const_array(mfi);
            Array4<int const> const& cdmskarr = cdmsk.const_array(mfi);
            Array4<int const> const& ndmskarr = nd_mask->const_array(mfi);
            Array4<int const> const& ccmskarr = cc_mask->const_array(mfi);
            Array4<RT const> const& fcocarr = fine_contrib_on_crse.const_array(mfi);

            if (csigma) {
                Array4<RT const> const& csigarr = csigma->const_array(mfi);
                AMREX_D_TERM(,
                             Array4<RT const> csigyarr = aniso_csigma
                                 ? m_sigma[crse_amrlev][0][1]->const_array(mfi) : Array4<RT const>{};,
                             Array4<RT const> csigzarr = aniso_csigma
                                 ? m_sigma[crse_amrlev][0][2]->const_array(mfi) : Array4<RT const>{};)
#if (AMREX_SPACEDIM == 2)
                if (aniso_csigma) {
                    AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                    {
                        mlndlap_res_cf_contrib_ha(i,j,k,resarr,csolarr,crhsarr,
                                                  csigarr,csigyarr,
                                                  cdmskarr,ndmskarr,ccmskarr,fcocarr,
                                                  cdxinv,c_cc_domain_p,c_nd_domain,
                                                  is_rz,
                                                  lobc,hibc, neumann_doubling);
                    });
                } else {
                    AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                    {
                        mlndlap_res_cf_contrib(i,j,k,resarr,csolarr,crhsarr,csigarr,
                                               cdmskarr,ndmskarr,ccmskarr,fcocarr,
                                               cdxinv,c_cc_domain_p,c_nd_domain,
                                               is_rz,
                                               lobc,hibc, neumann_doubling);
                    });
                }
#elif (AMREX_SPACEDIM == 3)
                if (aniso_csigma) {
                    AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                    {
                        mlndlap_res_cf_contrib_ha(i,j,k,resarr,csolarr,crhsarr,
                                                  csigarr,csigyarr,csigzarr,
                                                  cdmskarr,ndmskarr,ccmskarr,fcocarr,
                                                  cdxinv,c_cc_domain_p,c_nd_domain,
                                                  lobc,hibc, neumann_doubling);
                    });
                } else {
                    AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                    {
                        mlndlap_res_cf_contrib(i,j,k,resarr,csolarr,crhsarr,csigarr,
                                               cdmskarr,ndmskarr,ccmskarr,fcocarr,
                                               cdxinv,c_cc_domain_p,c_nd_domain,
                                               lobc,hibc, neumann_doubling);
                    });
                }
#endif
            } else {
                RT const_sigma = m_const_sigma;
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
    if constexpr (std::is_same_v<MF,MultiFab>) {
        amrex::EB_set_covered(res,0.0);
    }
#endif
#endif
}

template void MLNodeLaplacianT<MultiFab>::compSyncResidualCoarse (MultiFab&, const MultiFab&, const MultiFab&,
                                                                  const MultiFab*, const BoxArray&,
                                                                  const IntVect&);
template void MLNodeLaplacianT<MultiFab>::compSyncResidualFine (MultiFab&, const MultiFab&, const MultiFab&,
                                                                const MultiFab*);
template void MLNodeLaplacianT<MultiFab>::reflux (int, MultiFab&, const MultiFab&, const MultiFab&,
                                                  MultiFab&, MultiFab&, const MultiFab&) const;
#ifdef AMREX_MLMG_HAS_FLOAT
template void MLNodeLaplacianT<fMultiFab>::compSyncResidualCoarse (fMultiFab&, const fMultiFab&, const fMultiFab&,
                                                                   const fMultiFab*, const BoxArray&,
                                                                   const IntVect&);
template void MLNodeLaplacianT<fMultiFab>::compSyncResidualFine (fMultiFab&, const fMultiFab&, const fMultiFab&,
                                                                 const fMultiFab*);
template void MLNodeLaplacianT<fMultiFab>::reflux (int, fMultiFab&, const fMultiFab&, const fMultiFab&,
                                                   fMultiFab&, fMultiFab&, const fMultiFab&) const;
#endif

}

namespace amrex::detail { void mlndlap_anchor_sync () {} } // see mlndlap_link_parts
