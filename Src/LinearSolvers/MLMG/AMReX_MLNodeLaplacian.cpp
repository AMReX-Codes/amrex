#include <limits>

#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLNodeLap_K.H>
#include <AMReX_MultiFabUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_algoim.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

MLNodeLaplacian::MLNodeLaplacian (const Vector<Geometry>& a_geom,
                                  const Vector<BoxArray>& a_grids,
                                  const Vector<DistributionMapping>& a_dmap,
                                  const LPInfo& a_info,
                                  const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

#ifdef AMREX_USE_EB
MLNodeLaplacian::MLNodeLaplacian (const Vector<Geometry>& a_geom,
                                  const Vector<BoxArray>& a_grids,
                                  const Vector<DistributionMapping>& a_dmap,
                                  const LPInfo& a_info,
                                  const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}
#endif

MLNodeLaplacian::~MLNodeLaplacian ()
{}

void
MLNodeLaplacian::define (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info,
                         const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("MLNodeLaplacian::define()");

    // This makes sure grids are cell-centered;
    Vector<BoxArray> cc_grids = a_grids;
    for (auto& ba : cc_grids) {
        ba.enclosedCells();
    }

    MLNodeLinOp::define(a_geom, cc_grids, a_dmap, a_info, a_factory);

    m_sigma.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_sigma[amrlev].resize(m_num_mg_levels[amrlev]);
        const int mglev = 0;
        const int idim = 0;
        m_sigma[amrlev][mglev][idim].reset
            (new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1,
                          MFInfo(), *m_factory[amrlev][0]));
        m_sigma[amrlev][mglev][idim]->setVal(0.0);
    }

#ifdef AMREX_USE_EB
#if (AMREX_SPACEDIM == 2)
    const int ncomp_i = 5;
#else
    const int ncomp_i = algoim::numIntgs;
#endif
    m_integral.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
#ifdef AMREX_USE_EB
        m_integral[amrlev].reset
            (new MultiFab(m_grids[amrlev][0], m_dmap[amrlev][0], ncomp_i, 1,
                          MFInfo(), *m_factory[amrlev][0]));
#else
        m_integral[amrlev].reset
            (new MultiFab(m_grids[amrlev][0], m_dmap[amrlev][0], ncomp_i, 1));
#endif
    }
#endif

#if (AMREX_SPACEDIM == 2)
    m_is_rz = m_geom[0][0].IsRZ();
#endif
}

#ifdef AMREX_USE_EB
void
MLNodeLaplacian::define (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info,
                         const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    Vector<FabFactory<FArrayBox> const*> _factory;
    for (auto x : a_factory) {
        _factory.push_back(static_cast<FabFactory<FArrayBox> const*>(x));
    }
    define(a_geom, a_grids, a_dmap, a_info, _factory);
}
#endif

void
MLNodeLaplacian::unimposeNeumannBC (int amrlev, MultiFab& rhs) const
{
    if (m_coarsening_strategy == CoarseningStrategy::RAP) {
        const Box& nddom = amrex::surroundingNodes(Geom(amrlev).Domain());
        const auto lobc = LoBC();
        const auto hibc = HiBC();

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(rhs,mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& rhsarr = rhs.array(mfi);
            mlndlap_unimpose_neumann_bc(bx, rhsarr, nddom, lobc, hibc);
        }
    }
}

void
MLNodeLaplacian::setSigma (int amrlev, const MultiFab& a_sigma)
{
    MultiFab::Copy(*m_sigma[amrlev][0][0], a_sigma, 0, 0, 1, 0);
}

void
MLNodeLaplacian::compDivergence (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel)
{
    compRHS(rhs, vel, Vector<const MultiFab*>(), Vector<MultiFab*>());
}

void
MLNodeLaplacian::compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel,
                          const Vector<const MultiFab*>& rhnd,
                          const Vector<MultiFab*>& a_rhcc)
{
    BL_PROFILE("MLNodeLaplacian::compRHS()");

    if (!m_masks_built) buildMasks();

#ifdef AMREX_USE_EB
    if (!m_integral_built) buildIntegral();
#endif

#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

    const auto lobc = LoBC();
    const auto hibc = HiBC();

    Vector<std::unique_ptr<MultiFab> > rhcc(m_num_amr_levels);
    Vector<std::unique_ptr<MultiFab> > rhs_cc(m_num_amr_levels);

    for (int ilev = 0; ilev < m_num_amr_levels; ++ilev)
    {
        const Geometry& geom = m_geom[ilev][0];
        AMREX_ASSERT(vel[ilev]->nComp() >= AMREX_SPACEDIM);
        AMREX_ASSERT(vel[ilev]->nGrow() >= 1);
        vel[ilev]->FillBoundary(0, AMREX_SPACEDIM, IntVect(1), geom.periodicity());

        if (ilev < a_rhcc.size() && a_rhcc[ilev])
        {
            rhcc[ilev].reset(new MultiFab(a_rhcc[ilev]->boxArray(),
                                          a_rhcc[ilev]->DistributionMap(), 1, 1));
            rhcc[ilev]->setVal(0.0);
            MultiFab::Copy(*rhcc[ilev], *a_rhcc[ilev], 0, 0, 1, 0);
            rhcc[ilev]->FillBoundary(geom.periodicity());

            rhs_cc[ilev].reset(new MultiFab(rhs[ilev]->boxArray(),
                                            rhs[ilev]->DistributionMap(), 1, 0));
        }

        const auto dxinvarr = geom.InvCellSizeArray();
        const Box& nddom = amrex::surroundingNodes(geom.Domain());

        const iMultiFab& dmsk = *m_dirichlet_mask[ilev][0];

#ifdef AMREX_USE_EB
        auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[ilev][0].get());
        const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
        const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
        const MultiFab* intg = m_integral[ilev].get();
#endif

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
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
        const Geometry& cgeom = m_geom[ilev  ][0];
        const Geometry& fgeom = m_geom[ilev+1][0];

        frhs[ilev].reset(new MultiFab(amrex::coarsen(rhs[ilev+1]->boxArray(),2),
                                      rhs[ilev+1]->DistributionMap(), 1, 0));
        frhs[ilev]->setVal(0.0);

        const Box& cccdom = cgeom.Domain();
        const Box& nddom = amrex::surroundingNodes(fgeom.Domain());
        const auto fdxinv = fgeom.InvCellSizeArray();
        const iMultiFab& fdmsk = *m_dirichlet_mask[ilev+1][0];

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            FArrayBox vfab, rfab, rhccfab;
            for (MFIter mfi(*frhs[ilev],mfi_info); mfi.isValid(); ++mfi)
            {
                const Box& cvbx = mfi.validbox();
                const Box& fvbx = amrex::refine(cvbx,2);
                const Box& cbx = mfi.tilebox();
                const Box& fbx = amrex::refine(cbx,2);

                const Box& cc_fbx = amrex::enclosedCells(fbx);
                const Box& cc_fvbx = amrex::enclosedCells(fvbx);

                const Box& bx_vel = amrex::grow(cc_fbx,2) & amrex::grow(cc_fvbx,1);
                Box b = bx_vel & cc_fvbx;
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                {
                    if (m_lobc[0][idim] == LinOpBCType::inflow)
                    {
                        if (b.smallEnd(idim) == cccdom.smallEnd(idim)) {
                            b.growLo(idim, 1);
                        }
                    }
                    if (m_hibc[0][idim] == LinOpBCType::inflow)
                    {
                        if (b.bigEnd(idim) == cccdom.bigEnd(idim)) {
                            b.growHi(idim, 1);
                        }
                    }
                }

                vfab.resize(bx_vel, AMREX_SPACEDIM);
                Elixir veli = vfab.elixir();
                Array4<Real> const& varr = vfab.array();

                const Box& bx_rhs = amrex::grow(fbx,1);
                const Box& b2 = bx_rhs & amrex::grow(fvbx,-1);
                rfab.resize(bx_rhs);
                Elixir reli = rfab.elixir();
                Array4<Real> const& rarr = rfab.array();

                Array4<Real const> const& varr_orig = vel[ilev+1]->const_array(mfi);
                AMREX_HOST_DEVICE_FOR_4D(bx_vel, AMREX_SPACEDIM, i, j, k, n,
                {
                    if (b.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                        varr(i,j,k,n) = varr_orig(i,j,k,n);
                    } else {
                        varr(i,j,k,n) = 0.0;
                    }
                });

                Array4<Real const> const& rarr_orig = rhs[ilev+1]->const_array(mfi);
#if (AMREX_SPACEDIM == 2)
                AMREX_HOST_DEVICE_FOR_3D(bx_rhs, i, j, k,
                {
                    if (b2.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                        rarr(i,j,k) = rarr_orig(i,j,k);
                    } else {
                        rarr(i,j,k) = 0.0;
                    }
                    mlndlap_divu_compute_fine_contrib(i,j,k,fvbx,rarr,varr,fdxinv,
                                                      nddom,lobc,hibc,is_rz);
                });
#else
                AMREX_HOST_DEVICE_FOR_3D(bx_rhs, i, j, k,
                {
                    if (b2.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                        rarr(i,j,k) = rarr_orig(i,j,k);
                    } else {
                        rarr(i,j,k) = 0.0;
                    }
                    mlndlap_divu_compute_fine_contrib(i,j,k,fvbx,rarr,varr,fdxinv,
                                                      nddom,lobc,hibc);
                });
#endif

                Array4<Real> const& rhsarr = frhs[ilev]->array(mfi);
                Array4<int const> const& mskarr = fdmsk.const_array(mfi);
                AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                {
                    mlndlap_divu_add_fine_contrib(i,j,k,fvbx,rhsarr,rarr,mskarr);
                });

                if (rhcc[ilev+1])
                {
                    const Box& bx_rhcc = amrex::grow(cc_fbx,2);
                    const Box& b3 = bx_rhcc & cc_fvbx;

                    rhccfab.resize(bx_rhcc);
                    Elixir rhcceli = rhccfab.elixir();
                    Array4<Real> const& rhccarr = rhccfab.array();

                    Array4<Real const> const& rhccarr_orig = rhcc[ilev+1]->const_array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(bx_rhcc, i, j, k,
                    {
                        if (b3.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                            rhccarr(i,j,k) = rhccarr_orig(i,j,k);
                        } else {
                            rhccarr(i,j,k) = 0.0;
                        }
                    });

                    AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
                    {
                        mlndlap_rhcc_fine_contrib(i,j,k,fvbx,rhsarr,rhccarr,mskarr);
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
        const Box& cnddom = amrex::surroundingNodes(cccdom);
        const auto cdxinv = cgeom.InvCellSizeArray();
        const iMultiFab& cdmsk = *m_dirichlet_mask[ilev][0];
        const iMultiFab& c_nd_mask = *m_nd_fine_mask[ilev];
        const iMultiFab& c_cc_mask = *m_cc_fine_mask[ilev];
        const auto& has_fine_bndry = *m_has_fine_bndry[ilev];

        bool neumann_doubling = true; // yes even for RAP, because unimposeNeumannBC will be called on rhs

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*rhs[ilev]); mfi.isValid(); ++mfi)
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
                                            cdxinv,cnddom,lobc,hibc, neumann_doubling);
                });
#else
                AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                {
                    mlndlap_divu_cf_contrib(i,j,k,rhsarr,velarr,crhsarr,rhccarr,
                                            cdmskarr,ndmskarr,ccmskarr,
                                            cdxinv,cnddom,lobc,hibc, neumann_doubling);
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
}

void
MLNodeLaplacian::updateVelocity (const Vector<MultiFab*>& vel, const Vector<MultiFab const*>& sol) const
{
#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        const auto& sigma = *m_sigma[amrlev][0][0];
        const auto dxinv = m_geom[amrlev][0].InvCellSizeArray();
#ifdef AMREX_USE_EB
        auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
        const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
        const MultiFab* intg = m_integral[amrlev].get();
#endif
        for (MFIter mfi(*vel[amrlev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& varr = vel[amrlev]->array(mfi);
            Array4<Real const> const& solarr = sol[amrlev]->const_array(mfi);
            Array4<Real const> const& sigmaarr = sigma.const_array(mfi);
#ifdef AMREX_USE_EB
            bool regular = !factory;
            if (factory)
            {
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
            }
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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        const auto& sigma = *m_sigma[amrlev][0][0];
        const auto dxinv = m_geom[amrlev][0].InvCellSizeArray();
#ifdef AMREX_USE_EB
        auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
        const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
        const MultiFab* intg = m_integral[amrlev].get();
#endif

        // Initialize to zero because we only want -(sigma * grad(phi))

        for (MFIter mfi(sigma, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& farr = a_flux[amrlev]->array(mfi);
            Array4<Real const> const& solarr = a_sol[amrlev]->const_array(mfi);
            Array4<Real const> const& sigmaarr = sigma.array(mfi);

            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, AMREX_SPACEDIM, i, j, k, n,
            {
                farr(i,j,k,n) = 0.0;
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
            }
        }
    }
}

void
MLNodeLaplacian::averageDownCoeffs ()
{
    BL_PROFILE("MLNodeLaplacian::averageDownCoeffs()");

    if (m_coarsening_strategy == CoarseningStrategy::Sigma)
    {
        for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
        {
            for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
            {
                int ndims = (m_use_harmonic_average) ? AMREX_SPACEDIM : 1;
                for (int idim = 0; idim < ndims; ++idim)
                {
                    if (m_sigma[amrlev][mglev][idim] == nullptr) {
                        if (mglev == 0) {
                            m_sigma[amrlev][mglev][idim].reset
                                (new MultiFab(*m_sigma[amrlev][mglev][0], amrex::make_alias, 0, 1));
                        } else {
                            m_sigma[amrlev][mglev][idim].reset
                                (new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1));
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
        if (m_use_harmonic_average) {
            int mglev = 0;
            FillBoundaryCoeff(*m_sigma[amrlev][mglev][0], m_geom[amrlev][mglev]);
            for (mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
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
    if (m_coarsening_strategy != CoarseningStrategy::Sigma) return;

    const int nsigma = (m_use_harmonic_average) ? AMREX_SPACEDIM : 1;

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
#ifdef _OPENMP
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
                    } else if (idim == 1) {
#if (AMREX_SPACEDIM >= 2)
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
                        {
                            mlndlap_avgdown_coeff_y(i,j,k,cfab,ffab);
                        });
#endif
                    } else {
#if (AMREX_SPACEDIM == 3)
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
                        {
                            mlndlap_avgdown_coeff_z(i,j,k,cfab,ffab);
                        });
#endif
                    }
                }
            } else {
#ifdef _OPENMP
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
MLNodeLaplacian::FillBoundaryCoeff (MultiFab& sigma, const Geometry& geom)
{
    BL_PROFILE("MLNodeLaplacian::FillBoundaryCoeff()");

    sigma.FillBoundary(geom.periodicity());

    if (m_coarsening_strategy == CoarseningStrategy::Sigma)
    {
        const Box& domain = geom.Domain();
        const auto lobc = LoBC();
        const auto hibc = HiBC();

        MFItInfo mfi_info;
        if (Gpu::notInLaunchRegion()) mfi_info.SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(sigma, mfi_info); mfi.isValid(); ++mfi)
        {
            Array4<Real> const& sfab = sigma.array(mfi);
            mlndlap_fillbc_cc<Real>(mfi.validbox(),sfab,domain,lobc,hibc);
        }
    }
}

void
MLNodeLaplacian::buildStencil ()
{
    m_stencil.resize(m_num_amr_levels);
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
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            const int nghost = (0 == amrlev && mglev+1 == m_num_mg_levels[amrlev]) ? 1 : 4;
            m_stencil[amrlev][mglev].reset
                (new MultiFab(amrex::convert(m_grids[amrlev][mglev],
                                             IntVect::TheNodeVector()),
                              m_dmap[amrlev][mglev], ncomp_s, nghost));
            m_stencil[amrlev][mglev]->setVal(0.0);
        }

        {
            const Geometry& geom = m_geom[amrlev][0];
            const auto dxinvarr = geom.InvCellSizeArray();

#ifdef AMREX_USE_EB
            auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
            const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
            const MultiFab* intg = m_integral[amrlev].get();
            const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
#endif

            MFItInfo mfi_info;
            if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
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
                    }
                }
            }

            // set_stencil_s0 has to be in a separate MFIter from set_stencil
            // because it uses other cells' data.
#ifdef _OPENMP
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
            }

            m_stencil[amrlev][0]->FillBoundary(geom.periodicity());
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

#ifdef _OPENMP
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
#ifdef AMREX_USE_DPCPP
                // xxxxx DPCPP todo: this kernel hangs at JIT compilation
#ifndef AMREX_DPCPP_STENCIL_RAP_ON_GPU
                Gpu::LaunchSafeGuard lsg(false);
#endif
#endif
                AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                {
                    mlndlap_stencil_rap(i,j,k,csten,fsten);
                });
            }

#ifdef _OPENMP
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


    // This is only needed at the bottom.
    m_s0_norm0[0].back() = m_stencil[0].back()->norm0(0,0) * m_normalization_threshold;
}

void
MLNodeLaplacian::fixUpResidualMask (int amrlev, iMultiFab& resmsk)
{
    if (!m_masks_built) buildMasks();

    const iMultiFab& cfmask = *m_nd_fine_mask[amrlev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(resmsk,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<int> const& rmsk = resmsk.array(mfi);
        Array4<int const> const& fmsk = cfmask.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
        {
            if (fmsk(i,j,k) == crse_fine_node) rmsk(i,j,k) = 1;
        });
    }
}

void
MLNodeLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLNodeLaplacian::prepareForSolve()");

    MLNodeLinOp::prepareForSolve();

    buildMasks();

    averageDownCoeffs();

#ifdef AMREX_USE_EB
    buildIntegral();
#endif

    buildStencil();
}

void
MLNodeLaplacian::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
    BL_PROFILE("MLNodeLaplacian::restriction()");

    applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
        cfine.define(ba, fine.DistributionMap(), 1, 0);
    }

    MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][cmglev-1];

    const auto& stencil = m_stencil[amrlev][cmglev-1];

    bool regular_coarsening = true; int idir = 2;
    if (cmglev > 0) {
        regular_coarsening = mg_coarsen_ratio_vec[cmglev-1] == mg_coarsen_ratio;
        IntVect ratio = mg_coarsen_ratio_vec[cmglev-1];
        if (ratio[1] == 1) {
            idir = 1;
        } else if (ratio[0] == 1) {
            idir = 0;
        }
    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*pcrse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> cfab = pcrse->array(mfi);
        Array4<Real const> const& ffab = fine.const_array(mfi);
        Array4<int const> const& mfab = dmsk.const_array(mfi);
        if (m_coarsening_strategy == CoarseningStrategy::Sigma)
        {
	    if (regular_coarsening) 
	    {
            	AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            	{
               	    mlndlap_restriction(i,j,k,cfab,ffab,mfab);
            	});
	    }
	    else
	    {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_semi_restriction(i,j,k,cfab,ffab,mfab,idir);
                });
	    }
        }
        else
        {
            Array4<Real const> const& stfab = stencil->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_restriction_rap(i,j,k,cfab,ffab,stfab,mfab);
            });
        }
    }

    if (need_parallel_copy) {
        crse.ParallelCopy(cfine);
    }
}

void
MLNodeLaplacian::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
    BL_PROFILE("MLNodeLaplacian::interpolation()");

    const auto& sigma = m_sigma[amrlev][fmglev];
    const auto& stencil = m_stencil[amrlev][fmglev];

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

    bool regular_coarsening = true; int idir = 2;
    if (fmglev > 0) {
        regular_coarsening = mg_coarsen_ratio_vec[fmglev] == mg_coarsen_ratio;
        IntVect ratio = mg_coarsen_ratio_vec[fmglev];
	if (ratio[1] == 1) {
	    idir = 1;
        } else if (ratio[0] == 1) {
 	    idir = 0;
        }
    }


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const& ffab = fine.array(mfi);
        Array4<Real const> const& cfab = cmf->const_array(mfi);
        Array4<int const> const& mfab = dmsk.const_array(mfi);
        if (m_coarsening_strategy == CoarseningStrategy::RAP)
        {
            Array4<Real const> const& stfab = stencil->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_interpadd_rap(i,j,k,ffab,cfab,stfab,mfab);
            });
        }
        else if (m_use_harmonic_average && fmglev > 0)
        {
            AMREX_D_TERM(Array4<Real const> const& sxfab = sigma[0]->const_array(mfi);,
                         Array4<Real const> const& syfab = sigma[1]->const_array(mfi);,
                         Array4<Real const> const& szfab = sigma[2]->const_array(mfi););
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_interpadd_ha(i,j,k,ffab,cfab,AMREX_D_DECL(sxfab,syfab,szfab),mfab);
            });
        }
        else
        { 
	    Array4<Real const> const& sfab = sigma[0]->const_array(mfi);
	    if (regular_coarsening)
	    {
            	AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            	{
                    mlndlap_interpadd_aa(i,j,k,ffab,cfab,sfab,mfab);
                });
	    } 
	    else
	    {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    mlndlap_semi_interpadd_aa(i,j,k,ffab,cfab,sfab,mfab,idir);
                });
	    } 
        }
    }
}

void
MLNodeLaplacian::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs)
{
    const auto& amrrr = AMRRefRatio(camrlev);
    amrex::average_down(fine_sol, crse_sol, 0, 1, amrrr);

    if (isSingular(0))
    {
        MultiFab frhs(fine_rhs.boxArray(), fine_rhs.DistributionMap(), 1, 1);
        MultiFab::Copy(frhs, fine_rhs, 0, 0, 1, 0);
        restrictInteriorNodes(camrlev, crse_rhs, frhs);
    }
}

void
MLNodeLaplacian::restrictInteriorNodes (int camrlev, MultiFab& crhs, MultiFab& a_frhs) const
{
    const BoxArray& fba = a_frhs.boxArray();
    const DistributionMapping& fdm = a_frhs.DistributionMap();

    MultiFab* frhs = nullptr;
    std::unique_ptr<MultiFab> mf;
    if (a_frhs.nGrow() == 1)
    {
        frhs = &a_frhs;
    }
    else
    {
        mf.reset(new MultiFab(fba, fdm, 1, 1));
        frhs = mf.get();
        MultiFab::Copy(*frhs, a_frhs, 0, 0, 1, 0);
    }

    const Geometry& cgeom = m_geom[camrlev  ][0];

    const iMultiFab& fdmsk = *m_dirichlet_mask[camrlev+1][0];
    const auto& stencil    =  m_stencil[camrlev+1][0];

    MultiFab cfine(amrex::coarsen(fba, 2), fdm, 1, 0);

    frhs->setBndry(0.0);

    applyBC(camrlev+1, 0, *frhs, BCMode::Inhomogeneous, StateMode::Solution);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cfine, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& cfab = cfine.array(mfi);
        Array4<Real const> const& ffab = frhs->const_array(mfi);
        Array4<int const> const& mfab = fdmsk.const_array(mfi);
        if (m_coarsening_strategy == CoarseningStrategy::Sigma) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_restriction(i,j,k,cfab,ffab,mfab);
            });
        } else {
            Array4<Real const> const& stfab = stencil->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_restriction_rap(i,j,k,cfab,ffab,stfab,mfab);
            });
        }
    }

    MultiFab tmp_crhs(crhs.boxArray(), crhs.DistributionMap(), 1, 0);
    tmp_crhs.setVal(0.0);
    tmp_crhs.ParallelCopy(cfine, cgeom.periodicity());

    const iMultiFab& c_nd_mask = *m_nd_fine_mask[camrlev];
    const auto& has_fine_bndry = *m_has_fine_bndry[camrlev];

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(crhs, mfi_info); mfi.isValid(); ++mfi)
    {
        if (has_fine_bndry[mfi])
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const& dfab = crhs.array(mfi);
            Array4<Real const> const& sfab = tmp_crhs.const_array(mfi);
            Array4<int const> const& mfab = c_nd_mask.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                if (mfab(i,j,k) == fine_node) dfab(i,j,k) = sfab(i,j,k);
            });
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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
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
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                mlndlap_adotx_sten(i,j,k,yarr,xarr,stenarr,dmskarr);
            });
        }
        else if (m_use_harmonic_average && mglev > 0)
        {
            AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                         Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                         Array4<Real const> const& szarr = sigma[2]->const_array(mfi););
#if (AMREX_SPACEDIM == 2)
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                mlndlap_adotx_ha(i,j,k,yarr,xarr,AMREX_D_DECL(sxarr,syarr,szarr), dmskarr,
                                 is_rz, dxinvarr);
            });
#else
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                mlndlap_adotx_ha(i,j,k,yarr,xarr,AMREX_D_DECL(sxarr,syarr,szarr), dmskarr,
                                 dxinvarr);
            });
#endif
        }
        else
        {
            Array4<Real const> const& sarr = sigma[0]->const_array(mfi);
#if (AMREX_SPACEDIM == 2)
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                mlndlap_adotx_aa(i,j,k,yarr,xarr,sarr,dmskarr, is_rz, dxinvarr);
            });
#else
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                mlndlap_adotx_aa(i,j,k,yarr,xarr,sarr,dmskarr, dxinvarr);
            });
#endif
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
        constexpr int nsweeps = 4;
        MultiFab Ax(sol.boxArray(), sol.DistributionMap(), 1, 0);

        if (m_coarsening_strategy == CoarseningStrategy::RAP)
        {
            for (MFIter mfi(sol); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                Array4<Real> const& solarr = sol.array(mfi);
                Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                Array4<Real const> const& starr = stencil->const_array(mfi);
                Array4<int const> const& dmskarr = dmsk.const_array(mfi);
                Array4<Real> const& Axarr = Ax.array(mfi);

                for (int ns = 0; ns < nsweeps; ++ns) {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        mlndlap_adotx_sten(i,j,k,Axarr,solarr,starr,dmskarr);
                    });
                    AMREX_LAUNCH_DEVICE_LAMBDA ( bx, tbx,
                    {
                        mlndlap_jacobi_sten(tbx,solarr,Axarr,rhsarr,starr,dmskarr);
                    });
                }
            }
        }
        else if (m_use_harmonic_average && mglev > 0)
        {
            for (MFIter mfi(sol); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                             Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                             Array4<Real const> const& szarr = sigma[2]->const_array(mfi););
                Array4<Real> const& solarr = sol.array(mfi);
                Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                Array4<int const> const& dmskarr = dmsk.const_array(mfi);
                Array4<Real> const& Axarr = Ax.array(mfi);

                for (int ns = 0; ns < nsweeps; ++ns) {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        mlndlap_adotx_ha(i,j,k,Axarr,solarr,AMREX_D_DECL(sxarr,syarr,szarr), dmskarr,
#if (AMREX_SPACEDIM == 2)
                                         is_rz,
#endif
                                         dxinvarr);
                    });
                    AMREX_LAUNCH_DEVICE_LAMBDA ( bx, tbx,
                    {
                        mlndlap_jacobi_ha (tbx, solarr, Axarr, rhsarr, AMREX_D_DECL(sxarr,syarr,szarr),
                                           dmskarr, dxinvarr);
                    });
                }
            }
        }
        else
        {
            for (MFIter mfi(sol); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                Array4<Real const> const& sarr = sigma[0]->const_array(mfi);
                Array4<Real> const& solarr = sol.array(mfi);
                Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                Array4<int const> const& dmskarr = dmsk.const_array(mfi);
                Array4<Real> const& Axarr = Ax.array(mfi);

                for (int ns = 0; ns < nsweeps; ++ns) {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        mlndlap_adotx_aa(i,j,k,Axarr,solarr,sarr,dmskarr,
#if (AMREX_SPACEDIM == 2)
                                         is_rz,
#endif
                                         dxinvarr);
                    });
                    AMREX_LAUNCH_DEVICE_LAMBDA ( bx, tbx,
                    {
                        mlndlap_jacobi_aa (tbx, solarr, Axarr, rhsarr, sarr,
                                           dmskarr, dxinvarr);
                    });
                }
            }
        }

        if (nsweeps > 1) nodalSync(amrlev, mglev, sol);
    }
    else // cpu
#endif
    {
	bool regular_coarsening = true;
	if (amrlev == 0 and mglev > 0) 
    	{
            regular_coarsening = mg_coarsen_ratio_vec[mglev-1] == mg_coarsen_ratio;
        }

        constexpr int nsweeps = 2;
        if (m_use_gauss_seidel)
        {
            if (m_coarsening_strategy == CoarseningStrategy::RAP)
            {
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (MFIter mfi(sol); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.validbox();
                    Array4<Real> const& solarr = sol.array(mfi);
                    Array4<Real const> const& rhsarr = rhs.const_array(mfi);
                    Array4<Real const> const& starr = stencil->const_array(mfi);
                    Array4<int const> const& dmskarr = dmsk.const_array(mfi);

                    for (int ns = 0; ns < nsweeps; ++ns) {
                        mlndlap_gauss_seidel_sten(bx,solarr,rhsarr,starr,dmskarr);
                    }
                }
            }
            else if (m_use_harmonic_average && mglev > 0)
            {
#ifdef _OPENMP
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

                    for (int ns = 0; ns < nsweeps; ++ns) {
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
#ifdef _OPENMP
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
                        for (int ns = 0; ns < nsweeps; ++ns) {
                            mlndlap_gauss_seidel_aa(bx, solarr, rhsarr,
                                                    sarr, dmskarr, dxinvarr
#if (AMREX_SPACEDIM == 2)
                                                   ,is_rz
#endif
                                 );
                        }
		    } else {
			for (int ns = 0; ns < nsweeps; ++ns) {
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
#ifdef _OPENMP
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
            else if (m_use_harmonic_average && mglev > 0)
            {
#ifdef _OPENMP
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
#ifdef _OPENMP
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
MLNodeLaplacian::normalize (int amrlev, int mglev, MultiFab& mf) const
{
    BL_PROFILE("MLNodeLaplacian::normalize()");

    const auto& sigma = m_sigma[amrlev][mglev];
    const auto& stencil = m_stencil[amrlev][mglev];
    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];
    const Real s0_norm0 = m_s0_norm0[amrlev][mglev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& arr = mf.array(mfi);
        Array4<int const> const& dmskarr = dmsk.const_array(mfi);
        if (m_coarsening_strategy == CoarseningStrategy::RAP)
        {
            Array4<Real const> const& stenarr = stencil->const_array(mfi);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlndlap_normalize_sten(tbx,arr,stenarr,dmskarr,s0_norm0);
            });
        }
        else if (m_use_harmonic_average && mglev > 0)
        {
            AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                         Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                         Array4<Real const> const& szarr = sigma[2]->const_array(mfi););

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlndlap_normalize_ha(tbx,arr,AMREX_D_DECL(sxarr,syarr,szarr),dmskarr,dxinv);
            });
        }
        else
        {
            Array4<Real const> const& sarr = sigma[0]->const_array(mfi);

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlndlap_normalize_aa(tbx,arr,sarr,dmskarr,dxinv);
            });
        }
    }
}

void
MLNodeLaplacian::compSyncResidualCoarse (MultiFab& sync_resid, const MultiFab& a_phi,
                                         const MultiFab& vold, const MultiFab* rhcc,
                                         const BoxArray& fine_grids, const IntVect& ref_ratio)
{
    BL_PROFILE("MLNodeLaplacian::SyncResCrse()");

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
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(crse_cc_mask); mfi.isValid(); ++mfi)
    {
        Array4<int> const& fab = crse_cc_mask.array(mfi);
        mlndlap_fillbc_cc<int>(mfi.validbox(),fab,ccdom,lobc,hibc);
    }

    MultiFab phi(ndba, dmap, 1, 1);
#ifdef _OPENMP
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

    const MultiFab& sigma_orig = *m_sigma[0][0][0];
    const iMultiFab& dmsk = *m_dirichlet_mask[0][0];

#ifdef AMREX_USE_EB
    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[0][0].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    const MultiFab* intg = m_integral[0].get();
    const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
#endif

    bool neumann_doubling = true; // yes even for RAP, because unimposeNeumannBC will be called on rhs

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
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
		        if (b.contains(IntVect(AMREX_D_DECL(i,j,k))) and cccmsk(i,j,k)){
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
 			    if (b2.contains(IntVect(AMREX_D_DECL(i,j,k))) and cccmsk(i,j,k)){
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
                    Array4<Real const> const& sigmaarr_orig = sigma_orig.const_array(mfi);
#ifdef AMREX_USE_EB
                    if (typ == FabType::singlevalued)
                    {
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
                            if (ibx.contains(IntVect(AMREX_D_DECL(i,j,k))) and cccmsk(i,j,k)) {
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
                            mlndlap_adotx_sten(i, j, k, sync_resid_a, phiarr, stenarr, dmskarr);
                            mlndlap_crse_resid(i, j, k, sync_resid_a, rhsarr, cccmsk, nddom, lobc, hibc, neumann_doubling);
                        });
                    }
                    else
#endif
                    {
                        Array4<Real> sigmaarr = uarr;
                        const Box& ibx = ccbxg1 & amrex::enclosedCells(mfi.validbox());
                        AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                        {
                            if (ibx.contains(IntVect(AMREX_D_DECL(i,j,k))) and cccmsk(i,j,k)) {
                                sigmaarr(i,j,k) = sigmaarr_orig(i,j,k);
                            } else {
                                sigmaarr(i,j,k) = 0.0;
                            }
                        });

#if (AMREX_SPACEDIM == 2)
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                        {
                            mlndlap_adotx_aa(i, j, k, sync_resid_a, phiarr, sigmaarr, dmskarr, is_rz, dxinv);
                            mlndlap_crse_resid(i, j, k, sync_resid_a, rhsarr, cccmsk, nddom, lobc, hibc, neumann_doubling);
                        });
#else
                        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                        {
                            mlndlap_adotx_aa(i, j, k, sync_resid_a, phiarr, sigmaarr, dmskarr, dxinv);
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

    const MultiFab& sigma_orig = *m_sigma[0][0][0];
    const iMultiFab& dmsk = *m_dirichlet_mask[0][0];

    const auto lobc = LoBC();
    const auto hibc = HiBC();

#ifdef AMREX_USE_EB
    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[0][0].get());
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
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
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
                Array4<Real const> const& sigmaarr_orig = sigma_orig.const_array(mfi);
#ifdef AMREX_USE_EB
                if (typ == FabType::singlevalued)
                {
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
                            mlndlap_adotx_sten(i,j,k, sync_resid_a, phiarr, stenarr, tmpmaskarr);
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
                    AMREX_HOST_DEVICE_FOR_3D(ccbxg1, i, j, k,
                    {
                        if (ovlp2.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                            sigmaarr(i,j,k) = sigmaarr_orig(i,j,k);
                        } else {
                            sigmaarr(i,j,k) = 0.0;
                        }
                    });

#if (AMREX_SPACEDIM == 2)
                    AMREX_HOST_DEVICE_FOR_3D(gbx, i, j, k,
                    {
                        if (bx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                            mlndlap_adotx_aa(i,j,k, sync_resid_a, phiarr, sigmaarr, tmpmaskarr,
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
                            mlndlap_adotx_aa(i,j,k, sync_resid_a, phiarr, sigmaarr, tmpmaskarr,
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
                         MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const
{
    BL_PROFILE("MLNodeLaplacian::reflux()");

    const Geometry& cgeom = m_geom[crse_amrlev  ][0];
    const Geometry& fgeom = m_geom[crse_amrlev+1][0];
    const auto cdxinv = cgeom.InvCellSizeArray();
    const auto fdxinv = fgeom.InvCellSizeArray();
    const Box& c_cc_domain = cgeom.Domain();
    Box c_nd_domain = amrex::surroundingNodes(c_cc_domain);

    bool neumann_doubling = m_coarsening_strategy == CoarseningStrategy::Sigma;

#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

    const BoxArray& fba = fine_sol.boxArray();
    const DistributionMapping& fdm = fine_sol.DistributionMap();

    const iMultiFab& fdmsk = *m_dirichlet_mask[crse_amrlev+1][0];
    const auto& stencil    =  m_stencil[crse_amrlev+1][0];

    MultiFab fine_res_for_coarse(amrex::coarsen(fba, 2), fdm, 1, 0);

    applyBC(crse_amrlev+1, 0, fine_res, BCMode::Inhomogeneous, StateMode::Solution);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(fine_res_for_coarse, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const& cfab = fine_res_for_coarse.array(mfi);
        Array4<Real const> const& ffab = fine_res.const_array(mfi);
        Array4<int const> const& mfab = fdmsk.const_array(mfi);
        if (m_coarsening_strategy == CoarseningStrategy::Sigma) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_restriction(i,j,k,cfab,ffab,mfab);
            });
        } else {
            Array4<Real const> const& stfab = stencil->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                mlndlap_restriction_rap(i,j,k,cfab,ffab,stfab,mfab);
            });
        }
    }
    res.ParallelCopy(fine_res_for_coarse, cgeom.periodicity());

    MultiFab fine_contrib(amrex::coarsen(fba, 2), fdm, 1, 0);
    fine_contrib.setVal(0.0);

    const auto& fsigma = *m_sigma[crse_amrlev+1][0][0];

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox sigfab;
        FArrayBox Axfab;
        for (MFIter mfi(fine_contrib,mfi_info); mfi.isValid(); ++mfi)
        {
            const Box& cvbx = mfi.validbox();
            const Box& fvbx = amrex::refine(cvbx,2);
            const Box& cbx = mfi.tilebox();
            const Box& fbx = amrex::refine(cbx,2);

            const Box& cc_fbx = amrex::enclosedCells(fbx);
            const Box& cc_fvbx = amrex::enclosedCells(fvbx);
            const Box& bx_sig = amrex::grow(cc_fbx,2) & amrex::grow(cc_fvbx,1);
            const Box& b = bx_sig & cc_fvbx;

            sigfab.resize(bx_sig, 1);
            Elixir sigeli = sigfab.elixir();
            Array4<Real> const& sigarr = sigfab.array();
            Array4<Real const> const& sigarr_orig = fsigma.const_array(mfi);
            AMREX_HOST_DEVICE_FOR_3D(bx_sig, i, j, k,
            {
                if (b.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                    sigarr(i,j,k) = sigarr_orig(i,j,k);
                } else {
                    sigarr(i,j,k) = 0.0;
                }
            });

            const Box& bx_Ax = amrex::grow(fbx,1);
            const Box& b2 = bx_Ax & amrex::grow(fvbx,-1);
            Axfab.resize(bx_Ax);
            Elixir Axeli = Axfab.elixir();
            Array4<Real> const& Axarr = Axfab.array();
            Array4<Real const> const& rhsarr = fine_rhs.const_array(mfi);
            Array4<Real const> const& resarr = fine_res.const_array(mfi);
            Array4<Real const> const& solarr = fine_sol.const_array(mfi);
#if (AMREX_SPACEDIM == 2)
            AMREX_HOST_DEVICE_FOR_3D(bx_Ax, i, j, k,
            {
                if (b2.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                    Axarr(i,j,k) = rhsarr(i,j,k) - resarr(i,j,k);
                } else {
                    Axarr(i,j,k) = 0.0;
                }
                mlndlap_res_fine_Ax(i,j,k, fvbx, Axarr, solarr, sigarr, is_rz, fdxinv);
            });
#else
            AMREX_HOST_DEVICE_FOR_3D(bx_Ax, i, j, k,
            {
                if (b2.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                    Axarr(i,j,k) = rhsarr(i,j,k) - resarr(i,j,k);
                } else {
                    Axarr(i,j,k) = 0.0;
                }
                mlndlap_res_fine_Ax(i,j,k, fvbx, Axarr, solarr, sigarr, fdxinv);
            });
#endif
            Array4<Real> const& farr = fine_contrib.array(mfi);
            Array4<int const> const& marr = fdmsk.const_array(mfi);
            AMREX_HOST_DEVICE_FOR_3D(cbx, i, j, k,
            {
                mlndlap_res_fine_contrib(i,j,k,farr,Axarr,marr);
            });
        }
    }

    MultiFab fine_contrib_on_crse(crse_sol.boxArray(), crse_sol.DistributionMap(), 1, 0);
    fine_contrib_on_crse.setVal(0.0);
    fine_contrib_on_crse.ParallelAdd(fine_contrib, cgeom.periodicity());

    const iMultiFab& cdmsk = *m_dirichlet_mask[crse_amrlev][0];
    const auto& nd_mask     = m_nd_fine_mask[crse_amrlev];
    const auto& cc_mask     = m_cc_fine_mask[crse_amrlev];
    const auto& has_fine_bndry = m_has_fine_bndry[crse_amrlev];

    const auto lobc = LoBC();
    const auto hibc = HiBC();

    const auto& csigma = *m_sigma[crse_amrlev][0][0];

#ifdef _OPENMP
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
            Array4<Real const> const& csigarr = csigma.const_array(mfi);
            Array4<int const> const& cdmskarr = cdmsk.const_array(mfi);
            Array4<int const> const& ndmskarr = nd_mask->const_array(mfi);
            Array4<int const> const& ccmskarr = cc_mask->const_array(mfi);
            Array4<Real const> const& fcocarr = fine_contrib_on_crse.const_array(mfi);

#if (AMREX_SPACEDIM == 2)
            AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
            {
                mlndlap_res_cf_contrib(i,j,k,resarr,csolarr,crhsarr,csigarr,
                                       cdmskarr,ndmskarr,ccmskarr,fcocarr,
                                       cdxinv,c_nd_domain,
                                       is_rz,
                                       lobc,hibc, neumann_doubling);
            });
#else
            AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
            {
                mlndlap_res_cf_contrib(i,j,k,resarr,csolarr,crhsarr,csigarr,
                                       cdmskarr,ndmskarr,ccmskarr,fcocarr,
                                       cdxinv,c_nd_domain,
                                       lobc,hibc, neumann_doubling);
            });
#endif
        }
    }
#ifdef AMREX_USE_EB
    // Make sure to zero out the residual on any nodes completely surrounded by covered cells
    amrex::EB_set_covered(res,0.0);
#endif
}

#ifdef AMREX_USE_EB
void
MLNodeLaplacian::buildIntegral ()
{
    if (m_integral_built) return;

    BL_PROFILE("MLNodeLaplacian::buildIntegral()");

    m_integral_built = true;

#if (AMREX_SPACEDIM == 2)
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        MultiFab* intg = m_integral[amrlev].get();

        auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        if (factory)
        {
            const int ncomp = intg->nComp();
            const auto& flags = factory->getMultiEBCellFlagFab();
            const auto& vfrac = factory->getVolFrac();
            const auto& area = factory->getAreaFrac();
            const auto& bcent = factory->getBndryCent();

            MFItInfo mfi_info;
            if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*intg,mfi_info); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox();
                Array4<Real> const& garr = intg->array(mfi);
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
                    Array4<Real const> const& vfracarr = vfrac.const_array(mfi);
                    Array4<Real const> const& axarr = area[0]->const_array(mfi);
                    Array4<Real const> const& ayarr = area[1]->const_array(mfi);
                    Array4<Real const> const& bcarr = bcent.const_array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
                    {
                        mlndlap_set_integral_eb(i,j,k,garr,flagarr,vfracarr,axarr,ayarr,bcarr);
                    });
                }
            }
        }
    }
#else
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        amrex::algoim::compute_integrals(*m_integral[amrlev]);
    }
#endif
}
#endif

void
MLNodeLaplacian::checkPoint (std::string const& file_name) const
{
    if (ParallelContext::IOProcessorSub())
    {
        UtilCreateCleanDirectory(file_name, false);
        {
            std::string HeaderFileName(file_name+"/Header");
            std::ofstream HeaderFile;
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) {
                FileOpenFailed(HeaderFileName);
            }
            
            HeaderFile.precision(17);

            // MLLinop stuff
            HeaderFile << "verbose = " << verbose << "\n"
                       << "nlevs = " << NAMRLevels() << "\n"
                       << "do_agglomeration = " << info.do_agglomeration << "\n"
                       << "do_consolidation = " << info.do_consolidation << "\n"
                       << "agg_grid_size = " << info.agg_grid_size << "\n"
                       << "con_grid_size = " << info.con_grid_size << "\n"
                       << "has_metric_term = " << info.has_metric_term << "\n"
                       << "max_coarsening_level = " << info.max_coarsening_level << "\n";
#if (AMREX_SPACEDIM == 1)
            HeaderFile << "lobc = " << static_cast<int>(m_lobc[0][0]) << "\n";
#elif (AMREX_SPACEDIM == 2)
            HeaderFile << "lobc = " << static_cast<int>(m_lobc[0][0])
                       << " "       << static_cast<int>(m_lobc[0][1]) << "\n";
#else
            HeaderFile << "lobc = " << static_cast<int>(m_lobc[0][0])
                       << " "       << static_cast<int>(m_lobc[0][1])
                       << " "       << static_cast<int>(m_lobc[0][2]) << "\n";
#endif
#if (AMREX_SPACEDIM == 1)
            HeaderFile << "hibc = " << static_cast<int>(m_hibc[0][0]) << "\n";
#elif (AMREX_SPACEDIM == 2)
            HeaderFile << "hibc = " << static_cast<int>(m_hibc[0][0])
                       << " "       << static_cast<int>(m_hibc[0][1]) << "\n";
#else
            HeaderFile << "hibc = " << static_cast<int>(m_hibc[0][0])
                       << " "       << static_cast<int>(m_hibc[0][1])
                       << " "       << static_cast<int>(m_hibc[0][2]) << "\n";
#endif
            // m_coarse_data_for_bc: not used
            HeaderFile << "maxorder = " << getMaxOrder() << "\n";

            // MLNodeLaplacian stuff
            HeaderFile << "is_rz = " << m_is_rz << "\n";
            HeaderFile << "use_gauss_seidel = " << m_use_gauss_seidel << "\n";
            HeaderFile << "use_harmonic_average = " << m_use_harmonic_average << "\n";
            HeaderFile << "coarsen_strategy = " << static_cast<int>(m_coarsening_strategy) << "\n";
            // No level bc multifab
        }

        for (int ilev = 0; ilev < NAMRLevels(); ++ilev)
        {
            UtilCreateCleanDirectory(file_name+"/Level_"+std::to_string(ilev), false);
            std::string HeaderFileName(file_name+"/Level_"+std::to_string(ilev)+"/Header");
            std::ofstream HeaderFile;
            HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                    std::ofstream::trunc |
                                                    std::ofstream::binary);
            if( ! HeaderFile.good()) {
                FileOpenFailed(HeaderFileName);
            }
            
            HeaderFile.precision(17);

            HeaderFile << Geom(ilev) << "\n";
            m_grids[ilev][0].writeOn(HeaderFile);  HeaderFile << "\n";
        }
    }

    ParallelContext::BarrierSub();

    for (int ilev = 0; ilev < NAMRLevels(); ++ilev)
    {
        VisMF::Write(*m_sigma[ilev][0][0], file_name+"/Level_"+std::to_string(ilev)+"/sigma");
    }
}

#ifdef AMREX_USE_HYPRE

#if (AMREX_SPACEDIM == 2)
void
MLNodeLaplacian::fillIJMatrix (MFIter const& mfi, Array4<HypreNodeLap::Int const> const& nid,
                               Array4<int const> const& owner,
                               Vector<HypreNodeLap::Int>& ncols, Vector<HypreNodeLap::Int>& rows,
                               Vector<HypreNodeLap::Int>& cols, Vector<Real>& mat) const
{
    const int amrlev = 0;
    const int mglev  = m_num_mg_levels[amrlev]-1;

    const Box& ndbx = mfi.validbox();
    const auto lo = amrex::lbound(ndbx);
    const auto hi = amrex::ubound(ndbx);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE((mglev == 0) || (m_coarsening_strategy == CoarseningStrategy::RAP),
                                     "Coarsening strategy must be RAP to use hypre at mglev > 0");

    const auto& sten = m_stencil[amrlev][mglev]->array(mfi);

    constexpr int k = 0;
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if (nid(i,j,k) >= 0 && owner(i,j,k))
            {
                rows.push_back(nid(i,j,k));
                cols.push_back(nid(i,j,k));
                mat.push_back(sten(i,j,k,0));
                HypreNodeLap::Int nc = 1;

                if                (nid(i-1,j-1,k) >= 0) {
                    cols.push_back(nid(i-1,j-1,k));                  
                    mat.push_back(sten(i-1,j-1,k,3));
                    ++nc;
                }

                if                (nid(i,j-1,k) >= 0) {
                    cols.push_back(nid(i,j-1,k));
                    mat.push_back(sten(i,j-1,k,2));
                    ++nc;
                }

                if                (nid(i+1,j-1,k) >= 0) {
                    cols.push_back(nid(i+1,j-1,k));
                    mat.push_back(sten(i  ,j-1,k,3));
                    ++nc;
                }

                if                (nid(i-1,j,k) >= 0) {
                    cols.push_back(nid(i-1,j,k));
                    mat.push_back(sten(i-1,j,k,1));
                    ++nc;
                }

                if                (nid(i+1,j,k) >= 0) {
                    cols.push_back(nid(i+1,j,k));
                    mat.push_back(sten(i  ,j,k,1));
                    ++nc;
                }

                if                (nid(i-1,j+1,k) >= 0) {
                    cols.push_back(nid(i-1,j+1,k));
                    mat.push_back(sten(i-1,j  ,k,3));
                    ++nc;
                }

                if                (nid(i,j+1,k) >= 0) {
                    cols.push_back(nid(i,j+1,k));
                    mat.push_back(sten(i,j  ,k,2));
                    ++nc;
                }

                if                (nid(i+1,j+1,k) >= 0) {
                    cols.push_back(nid(i+1,j+1,k));
                    mat.push_back(sten(i  ,j  ,k,3));
                    ++nc;
                }

                ncols.push_back(nc);
            }
        }
    }
}
#else
void
MLNodeLaplacian::fillIJMatrix (MFIter const& mfi, Array4<HypreNodeLap::Int const> const& nid,
                               Array4<int const> const& owner,
                               Vector<HypreNodeLap::Int>& ncols, Vector<HypreNodeLap::Int>& rows,
                               Vector<HypreNodeLap::Int>& cols, Vector<Real>& mat) const
{
    const int amrlev = 0;
    const int mglev  = m_num_mg_levels[amrlev]-1;

    const Box& ndbx = mfi.validbox();
    const auto lo = amrex::lbound(ndbx);
    const auto hi = amrex::ubound(ndbx);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE((mglev == 0) || (m_coarsening_strategy == CoarseningStrategy::RAP),
                                     "Coarsening strategy must be RAP to use hypre at mglev > 0");

    const auto& sten = m_stencil[amrlev][mglev]->array(mfi);

    constexpr int ist_000 = 1-1;
    constexpr int ist_p00 = 2-1;
    constexpr int ist_0p0 = 3-1;
    constexpr int ist_00p = 4-1;
    constexpr int ist_pp0 = 5-1;
    constexpr int ist_p0p = 6-1;
    constexpr int ist_0pp = 7-1;
    constexpr int ist_ppp = 8-1;

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (nid(i,j,k) >= 0 && owner(i,j,k))
                {
                    rows.push_back(nid(i,j,k));
                    cols.push_back(nid(i,j,k));
                    mat.push_back(sten(i,j,k,ist_000));
                    HypreNodeLap::Int nc = 1;

                    if                (nid(i-1,j-1,k-1) >= 0) {
                        cols.push_back(nid(i-1,j-1,k-1));                  
                        mat.push_back(sten(i-1,j-1,k-1,ist_ppp));
                        ++nc;
                    }

                    if                (nid(i,j-1,k-1) >= 0) {
                        cols.push_back(nid(i,j-1,k-1));
                        mat.push_back(sten(i,j-1,k-1,ist_0pp));
                        ++nc;
                    }

                    if                (nid(i+1,j-1,k-1) >= 0) {
                        cols.push_back(nid(i+1,j-1,k-1));
                        mat.push_back(sten(i,j-1,k-1,ist_ppp));
                        ++nc;
                    }

                    if                (nid(i-1,j,k-1) >= 0) {
                        cols.push_back(nid(i-1,j,k-1));
                        mat.push_back(sten(i-1,j,k-1,ist_p0p));
                        ++nc;
                    }

                    if                (nid(i,j,k-1) >= 0) {
                        cols.push_back(nid(i,j,k-1));
                        mat.push_back(sten(i,j,k-1,ist_00p));
                        ++nc;
                    }

                    if                (nid(i+1,j,k-1) >= 0) {
                        cols.push_back(nid(i+1,j,k-1));
                        mat.push_back(sten(i,j,k-1,ist_p0p));
                        ++nc;
                    }

                    if                (nid(i-1,j+1,k-1) >= 0) {
                        cols.push_back(nid(i-1,j+1,k-1));
                        mat.push_back(sten(i-1,j,k-1,ist_ppp));
                        ++nc;
                    }

                    if                (nid(i,j+1,k-1) >= 0) {
                        cols.push_back(nid(i,j+1,k-1));
                        mat.push_back(sten(i,j,k-1,ist_0pp));
                        ++nc;
                    }

                    if                (nid(i+1,j+1,k-1) >= 0) {
                        cols.push_back(nid(i+1,j+1,k-1));
                        mat.push_back(sten(i,j,k-1,ist_ppp));
                        ++nc;
                    }

                    if                (nid(i-1,j-1,k) >= 0) {
                        cols.push_back(nid(i-1,j-1,k));
                        mat.push_back(sten(i-1,j-1,k,ist_pp0));
                        ++nc;
                    }

                    if                (nid(i,j-1,k) >= 0) {
                        cols.push_back(nid(i,j-1,k));
                        mat.push_back(sten(i,j-1,k,ist_0p0));
                        ++nc;
                    }

                    if                (nid(i+1,j-1,k) >= 0) {
                        cols.push_back(nid(i+1,j-1,k));
                        mat.push_back(sten(i,j-1,k,ist_pp0));
                        ++nc;
                    }

                    if                (nid(i-1,j,k) >= 0) {
                        cols.push_back(nid(i-1,j,k));
                        mat.push_back(sten(i-1,j,k,ist_p00));
                        ++nc;
                    }

                    if                (nid(i+1,j,k) >= 0) {
                        cols.push_back(nid(i+1,j,k));
                        mat.push_back(sten(i,j,k,ist_p00));
                        ++nc;
                    }

                    if                (nid(i-1,j+1,k) >= 0) {
                        cols.push_back(nid(i-1,j+1,k));
                        mat.push_back(sten(i-1,j,k,ist_pp0));
                        ++nc;
                    }

                    if                (nid(i,j+1,k) >= 0) {
                        cols.push_back(nid(i,j+1,k));
                        mat.push_back(sten(i,j,k,ist_0p0));
                        ++nc;
                    }

                    if                (nid(i+1,j+1,k) >= 0) {
                        cols.push_back(nid(i+1,j+1,k));
                        mat.push_back(sten(i,j,k,ist_pp0));
                        ++nc;
                    }

                    if                (nid(i-1,j-1,k+1) >= 0) {
                        cols.push_back(nid(i-1,j-1,k+1));
                        mat.push_back(sten(i-1,j-1,k,ist_ppp));
                        ++nc;
                    }

                    if                (nid(i,j-1,k+1) >= 0) {
                        cols.push_back(nid(i,j-1,k+1));
                        mat.push_back(sten(i,j-1,k,ist_0pp));
                        ++nc;
                    }

                    if                (nid(i+1,j-1,k+1) >= 0) {
                        cols.push_back(nid(i+1,j-1,k+1));
                        mat.push_back(sten(i,j-1,k,ist_ppp));
                        ++nc;
                    }

                    if                (nid(i-1,j,k+1) >= 0) {
                        cols.push_back(nid(i-1,j,k+1));
                        mat.push_back(sten(i-1,j,k,ist_p0p));
                        ++nc;
                    }

                    if                (nid(i,j,k+1) >= 0) {
                        cols.push_back(nid(i,j,k+1));
                        mat.push_back(sten(i,j,k,ist_00p));
                        ++nc;
                    }

                    if                (nid(i+1,j,k+1) >= 0) {
                        cols.push_back(nid(i+1,j,k+1));
                        mat.push_back(sten(i,j,k,ist_p0p));
                        ++nc;
                    }

                    if                (nid(i-1,j+1,k+1) >= 0) {
                        cols.push_back(nid(i-1,j+1,k+1));
                        mat.push_back(sten(i-1,j,k,ist_ppp));
                        ++nc;
                    }

                    if                (nid(i,j+1,k+1) >= 0) {
                        cols.push_back(nid(i,j+1,k+1));
                        mat.push_back(sten(i,j,k,ist_0pp));
                        ++nc;
                    }

                    if                (nid(i+1,j+1,k+1) >= 0) {
                        cols.push_back(nid(i+1,j+1,k+1));
                        mat.push_back(sten(i,j,k,ist_ppp));
                        ++nc;
                    }

                    ncols.push_back(nc);
                }
            }
        }
    }
}
#endif

#endif

}
