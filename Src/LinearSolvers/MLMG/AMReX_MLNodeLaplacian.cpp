
#include <limits>

#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLNodeLap_F.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BaseFab_f.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

namespace {
    const Real bogus_value = std::numeric_limits<Real>::quiet_NaN();
}

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
    Vector<FabFactory<FArrayBox> const*> _factory;
    for (auto x : a_factory) {
        _factory.push_back(static_cast<FabFactory<FArrayBox> const*>(x));
    }
    define(a_geom, a_grids, a_dmap, a_info, _factory);
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
            (new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1));
        m_sigma[amrlev][mglev][idim]->setVal(0.0);
    }

#ifdef AMREX_USE_EB
    const int ncomp_c = (AMREX_SPACEDIM == 2) ? 6 : 24;
    const int ncomp_i = AMREX_SPACEDIM;
    m_connection.resize(m_num_amr_levels);
    m_integral.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_integral[amrlev].reset
            (new MultiFab(m_grids[amrlev][0], m_dmap[amrlev][0], ncomp_i, 1));
        m_integral[amrlev]->setVal(0.0);

        m_connection[amrlev].reset
            (new MultiFab(m_grids[amrlev][0], m_dmap[amrlev][0], ncomp_c, 1));
        m_connection[amrlev]->setVal(1.0);
    }
#endif

#if (AMREX_SPACEDIM == 2)
    m_is_rz = Geometry::IsRZ();
#endif
}

void
MLNodeLaplacian::setSigma (int amrlev, const MultiFab& a_sigma)
{
    MultiFab::Copy(*m_sigma[amrlev][0][0], a_sigma, 0, 0, 1, 0);
}

void
MLNodeLaplacian::compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel,
                          const Vector<const MultiFab*>& rhnd,
                          const Vector<MultiFab*>& a_rhcc)
{
    BL_PROFILE("MLNodeLaplacian::compRHS()");

    if (!m_masks_built) buildMasks();

#if AMREX_USE_EB
    if (!m_connection_built) buildConnection();
#endif

    Vector<std::unique_ptr<MultiFab> > rhcc(m_num_amr_levels);
    Vector<std::unique_ptr<MultiFab> > rhs_cc(m_num_amr_levels);

    for (int ilev = 0; ilev < m_num_amr_levels; ++ilev)
    {
        const Geometry& geom = m_geom[ilev][0];
        AMREX_ASSERT(vel[ilev]->nComp() >= AMREX_SPACEDIM);
        AMREX_ASSERT(vel[ilev]->nGrow() >= 1);
        vel[ilev]->FillBoundary(0, AMREX_SPACEDIM, geom.periodicity());

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

        const Real* dxinv = geom.InvCellSize();
        const Box& nddom = amrex::surroundingNodes(geom.Domain());

        const iMultiFab& dmsk = *m_dirichlet_mask[ilev][0];

#ifdef AMREX_USE_EB
        auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[ilev][0].get());
        const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
        const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
        const MultiFab* intg = m_integral[ilev].get();
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*rhs[ilev], MFItInfo().EnableTiling().SetDynamic(true)); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

#ifdef AMREX_USE_EB
            bool regular = !factory;
            if (factory)
            {
                const auto& flag = (*flags)[mfi];
                const auto& ccbx = amrex::grow(amrex::enclosedCells(bx),1);
                const auto& typ = flag.getType(ccbx);
                if (typ == FabType::covered)
                {
                    (*rhs[ilev])[mfi].setVal(0.0, bx);
                }
                else if (typ == FabType::singlevalued)
                {
                    amrex_mlndlap_divu_eb(BL_TO_FORTRAN_BOX(bx),
                                          BL_TO_FORTRAN_ANYD((*rhs[ilev])[mfi]),
                                          BL_TO_FORTRAN_ANYD((*vel[ilev])[mfi]),
                                          BL_TO_FORTRAN_ANYD((*vfrac)[mfi]),
                                          BL_TO_FORTRAN_ANYD((*intg)[mfi]),
                                          BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                          dxinv);
                }
                else
                {
                    regular = true;
                }
            }
            if (regular)
#endif
            {
                amrex_mlndlap_divu(BL_TO_FORTRAN_BOX(bx),
                                   BL_TO_FORTRAN_ANYD((*rhs[ilev])[mfi]),
                                   BL_TO_FORTRAN_ANYD((*vel[ilev])[mfi]),
                                   BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                   dxinv);
            }

            if (m_coarsening_strategy == CoarseningStrategy::Sigma) {
                amrex_mlndlap_impose_neumann_bc(BL_TO_FORTRAN_BOX(bx),
                                                BL_TO_FORTRAN_ANYD((*rhs[ilev])[mfi]),
                                                BL_TO_FORTRAN_BOX(nddom),
                                                m_lobc.data(), m_hibc.data());
            }

            if (rhcc[ilev])
            {
                amrex_mlndlap_rhcc(BL_TO_FORTRAN_BOX(bx),
                                   BL_TO_FORTRAN_ANYD((*rhs_cc[ilev])[mfi]),
                                   BL_TO_FORTRAN_ANYD((*rhcc[ilev])[mfi]),
                                   BL_TO_FORTRAN_ANYD(dmsk[mfi]));

                if (m_coarsening_strategy == CoarseningStrategy::Sigma) {
                    amrex_mlndlap_impose_neumann_bc(BL_TO_FORTRAN_BOX(bx),
                                                    BL_TO_FORTRAN_ANYD((*rhs_cc[ilev])[mfi]),
                                                    BL_TO_FORTRAN_BOX(nddom),
                                                    m_lobc.data(), m_hibc.data());
                }
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
        const Real* fdxinv = fgeom.InvCellSize();
        const iMultiFab& fdmsk = *m_dirichlet_mask[ilev+1][0];

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            FArrayBox vfab;
            FArrayBox rfab;
            for (MFIter mfi(*frhs[ilev], MFItInfo().EnableTiling().SetDynamic(true));
                 mfi.isValid(); ++mfi)
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
                    if (m_lobc[idim] == LinOpBCType::inflow)
                    {
                        if (b.smallEnd(idim) == cccdom.smallEnd(idim)) {
                            b.growLo(idim, 1);
                        }
                    }
                    if (m_hibc[idim] == LinOpBCType::inflow)
                    {
                        if (b.bigEnd(idim) == cccdom.bigEnd(idim)) {
                            b.growHi(idim, 1);
                        }
                    }
                }
                vfab.resize(bx_vel, AMREX_SPACEDIM);
                vfab.setVal(0.0);
                vfab.copy((*vel[ilev+1])[mfi], b, 0, b, 0, AMREX_SPACEDIM);

                const Box& bx_rhs = amrex::grow(fbx,1);
                const Box& b2 = bx_rhs & amrex::grow(fvbx,-1);
                rfab.resize(bx_rhs);
                rfab.setVal(0.0);
                rfab.copy((*rhs[ilev+1])[mfi], b2, 0, b2, 0, 1);

                amrex_mlndlap_divu_fine_contrib(BL_TO_FORTRAN_BOX(cbx),
                                                BL_TO_FORTRAN_BOX(cvbx),
                                                BL_TO_FORTRAN_ANYD((*frhs[ilev])[mfi]),
                                                BL_TO_FORTRAN_ANYD(vfab),
                                                BL_TO_FORTRAN_ANYD(rfab),
                                                BL_TO_FORTRAN_ANYD(fdmsk[mfi]),
                                                fdxinv);

                if (rhcc[ilev+1])
                {
                    const Box& bx_rhcc = amrex::grow(cc_fbx,2);
                    const Box& b3 = bx_rhcc & cc_fvbx;
                    FArrayBox* rhcc_fab = &vfab;
                    rhcc_fab->resize(bx_rhcc);
                    rhcc_fab->setVal(0.0);
                    rhcc_fab->copy((*rhcc[ilev+1])[mfi], b3, 0, b3, 0, 1);
                    amrex_mlndlap_rhcc_fine_contrib(BL_TO_FORTRAN_BOX(cbx),
                                                    BL_TO_FORTRAN_BOX(cvbx),
                                                    BL_TO_FORTRAN_ANYD((*frhs[ilev])[mfi]),
                                                    BL_TO_FORTRAN_ANYD(*rhcc_fab),
                                                    BL_TO_FORTRAN_ANYD(fdmsk[mfi]));
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
        const Real* cdxinv = cgeom.InvCellSize();
        const iMultiFab& cdmsk = *m_dirichlet_mask[ilev][0];
        const iMultiFab& c_nd_mask = *m_nd_fine_mask[ilev];
        const iMultiFab& c_cc_mask = *m_cc_fine_mask[ilev];
        const auto& has_fine_bndry = *m_has_fine_bndry[ilev];

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*rhs[ilev], MFItInfo().EnableTiling().SetDynamic(true));
             mfi.isValid(); ++mfi)
        {
            if (has_fine_bndry[mfi])
            {
                const Box& bx = mfi.tilebox();

                if (rhcc[ilev])
                {
                    amrex_mlndlap_rhcc_crse_contrib(BL_TO_FORTRAN_BOX(bx),
                                                    BL_TO_FORTRAN_ANYD(crhs[mfi]),
                                                    BL_TO_FORTRAN_ANYD((*rhcc[ilev])[mfi]),
                                                    BL_TO_FORTRAN_ANYD(cdmsk[mfi]),
                                                    BL_TO_FORTRAN_ANYD(c_nd_mask[mfi]),
                                                    BL_TO_FORTRAN_ANYD(c_cc_mask[mfi]));
                }

                amrex_mlndlap_divu_cf_contrib(BL_TO_FORTRAN_BOX(bx),
                                              BL_TO_FORTRAN_ANYD((*rhs[ilev])[mfi]),
                                              BL_TO_FORTRAN_ANYD((*vel[ilev])[mfi]),
                                              BL_TO_FORTRAN_ANYD(cdmsk[mfi]),
                                              BL_TO_FORTRAN_ANYD(c_nd_mask[mfi]),
                                              BL_TO_FORTRAN_ANYD(c_cc_mask[mfi]),
                                              BL_TO_FORTRAN_ANYD(crhs[mfi]),
                                              cdxinv, BL_TO_FORTRAN_BOX(cnddom),
                                              m_lobc.data(), m_hibc.data());
            }
        }
    }

    for (int ilev = 0; ilev < m_num_amr_levels; ++ilev)
    {
        if (ilev < rhnd.size() && rhnd[ilev]) {
            MultiFab::Add(*rhs[ilev], *rhnd[ilev], 0, 0, 1, 0);
        }
    }
}

void
MLNodeLaplacian::updateVelocity (const Vector<MultiFab*>& vel, const Vector<MultiFab const*>& sol) const
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        const auto& sigma = *m_sigma[amrlev][0][0];
        const Real* dxinv = m_geom[amrlev][0].InvCellSize();
#ifdef AMREX_USE_EB
        auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
        const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
        const MultiFab* intg = m_integral[amrlev].get();
#endif
        for (MFIter mfi(*vel[amrlev], true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto& vfab = (*vel[amrlev])[mfi];
#ifdef AMREX_USE_EB
            bool regular = !factory;
            if (factory)
            {
                auto type = (*flags)[mfi].getType(bx);
                if (type == FabType::covered)
                {
                    vfab.setVal(0.0, bx, 0, AMREX_SPACEDIM);
                }
                else if (type == FabType::singlevalued)
                {
                    amrex_mlndlap_mknewu_eb(BL_TO_FORTRAN_BOX(bx),
                                            BL_TO_FORTRAN_ANYD(vfab),
                                            BL_TO_FORTRAN_ANYD((*sol[amrlev])[mfi]),
                                            BL_TO_FORTRAN_ANYD(sigma[mfi]),
                                            BL_TO_FORTRAN_ANYD((*vfrac)[mfi]),
                                            BL_TO_FORTRAN_ANYD((*intg)[mfi]),
                                            dxinv);
                }
                else
                {
                    regular = true;
                }
            }
            if (regular)
#endif
            {
                amrex_mlndlap_mknewu(BL_TO_FORTRAN_BOX(bx),
                                     BL_TO_FORTRAN_ANYD(vfab),
                                     BL_TO_FORTRAN_ANYD((*sol[amrlev])[mfi]),
                                     BL_TO_FORTRAN_ANYD(sigma[mfi]),
                                     dxinv);
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
    amrex::average_down(*m_sigma[flev][mglev][idim], *m_sigma[flev-1][mglev][idim], 0, 1,
                        m_amr_ref_ratio[flev-1]);
}

void
MLNodeLaplacian::averageDownCoeffsSameAmrLevel (int amrlev)
{
    if (m_coarsening_strategy != CoarseningStrategy::Sigma) return;

    const int nsigma = (m_use_harmonic_average) ? AMREX_SPACEDIM : 1;

    for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
    {
        for (int idim = 0; idim < nsigma; ++idim)
        {
            const MultiFab& fine = *m_sigma[amrlev][mglev-1][idim];
            MultiFab& crse = *m_sigma[amrlev][mglev][idim];
            bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
            MultiFab cfine;
            if (need_parallel_copy) {
                const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
                cfine.define(ba, fine.DistributionMap(), 1, 0);
            }

            MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(*pcrse, true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                amrex_mlndlap_avgdown_coeff(BL_TO_FORTRAN_BOX(bx),
                                            BL_TO_FORTRAN_ANYD((*pcrse)[mfi]),
                                            BL_TO_FORTRAN_ANYD(fine[mfi]),
                                            &idim);
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

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(sigma, MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
        {
            if (!domain.contains(mfi.fabbox()))
            {
                amrex_mlndlap_fillbc_cc(BL_TO_FORTRAN_ANYD(sigma[mfi]),
                                        BL_TO_FORTRAN_BOX(domain),
                                        m_lobc.data(), m_hibc.data());
            }
        }
    }
}

void
MLNodeLaplacian::buildMasks ()
{
    if (m_masks_built) return;

    BL_PROFILE("MLNodeLaplacian::buildMasks()");

    m_masks_built = true;

    m_is_bottom_singular = false;
    auto itlo = std::find(m_lobc.begin(), m_lobc.end(), BCType::Dirichlet);
    auto ithi = std::find(m_hibc.begin(), m_hibc.end(), BCType::Dirichlet);
    if (itlo == m_lobc.end() && ithi == m_hibc.end())
    {  // No Dirichlet
        m_is_bottom_singular = m_domain_covered[0];
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector< std::pair<int,Box> > isects;
        IArrayBox ccfab;

        for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
        {
            for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
            {
                const Geometry& geom = m_geom[amrlev][mglev];
                const auto& period = geom.periodicity();
                const Box& ccdomain = geom.Domain();
                const Box& nddomain = amrex::surroundingNodes(ccdomain);
                const std::vector<IntVect>& pshifts = period.shiftIntVect();

                Box ccdomain_p = ccdomain;
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    if (Geometry::isPeriodic(idim)) {
                        ccdomain_p.grow(idim, 1);
                    }
                }

                {
                    auto& dmask = *m_dirichlet_mask[amrlev][mglev];
                    const BoxArray& ccba = m_grids[amrlev][mglev];

                    for (MFIter mfi(dmask, MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
                    {
                        const Box& ndbx = mfi.validbox();
                        const Box& ccbx = amrex::enclosedCells(ndbx);
                        const Box& ccbxg1 = amrex::grow(ccbx,1);
                        IArrayBox& mskfab = dmask[mfi];
                        
                        ccfab.resize(ccbxg1);
                        ccfab.setVal(1);
                        ccfab.setComplement(2,ccdomain_p,0,1);

                        for (const auto& iv : pshifts)
                        {
                            ccba.intersections(ccbxg1+iv, isects);
                            for (const auto& is : isects)
                            {
                                ccfab.setVal(0, is.second-iv, 0, 1);
                            }
                        }
                        
                        amrex_mlndlap_set_dirichlet_mask(BL_TO_FORTRAN_ANYD(mskfab),
                                                         BL_TO_FORTRAN_ANYD(ccfab),
                                                         BL_TO_FORTRAN_BOX(nddomain),
                                                         m_lobc.data(), m_hibc.data());
                    }
                }
            }
        }
    }

    for (int amrlev = 0; amrlev < m_num_amr_levels-1; ++amrlev)
    {
        iMultiFab& cc_mask = *m_cc_fine_mask[amrlev];
        iMultiFab& nd_mask = *m_nd_fine_mask[amrlev];
        LayoutData<int>& has_cf = *m_has_fine_bndry[amrlev];
        const BoxArray& fba = m_grids[amrlev+1][0];
        const BoxArray& cfba = amrex::coarsen(fba, AMRRefRatio(amrlev));

        const Box& ccdom = m_geom[amrlev][0].Domain();

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMRRefRatio(amrlev) == 2, "ref_ratio != 0 not supported");

        cc_mask.setVal(0);  // coarse by default

        const std::vector<IntVect>& pshifts = m_geom[amrlev][0].periodicity().shiftIntVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            std::vector< std::pair<int,Box> > isects;

            for (MFIter mfi(cc_mask); mfi.isValid(); ++mfi)
            {
                has_cf[mfi] = 0;
                IArrayBox& fab = cc_mask[mfi];
                const Box& bx = fab.box();
                for (const auto& iv : pshifts)
                {
                    cfba.intersections(bx+iv, isects);
                    for (const auto& is : isects)
                    {
                        fab.setVal(1, is.second-iv, 0, 1);
                    }
                    if (!isects.empty()) has_cf[mfi] = 1;
                }

                amrex_mlndlap_fillbc_cc_i(BL_TO_FORTRAN_ANYD(fab),
                                            BL_TO_FORTRAN_BOX(ccdom),
                                            m_lobc.data(), m_hibc.data());
            }
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(nd_mask,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            amrex_mlndlap_set_nodal_mask(BL_TO_FORTRAN_BOX(bx),
                                         BL_TO_FORTRAN_ANYD(nd_mask[mfi]),
                                         BL_TO_FORTRAN_ANYD(cc_mask[mfi]));
        }
    }

    auto& has_cf = *m_has_fine_bndry[m_num_amr_levels-1];
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(has_cf); mfi.isValid(); ++mfi)
    {
        has_cf[mfi] = 0;
    }

    {
        int amrlev = 0;
        int mglev = m_num_mg_levels[amrlev]-1;
        const iMultiFab& omask = *m_owner_mask[amrlev][mglev];
        m_bottom_dot_mask.define(omask.boxArray(), omask.DistributionMap(), 1, 0);

        const Geometry& geom = m_geom[amrlev][mglev];
        Box nddomain = amrex::surroundingNodes(geom.Domain());

        if (m_coarsening_strategy != CoarseningStrategy::Sigma) {
            nddomain.grow(1000); // hack to avoid masks being modified at Neuman boundary
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(m_bottom_dot_mask,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto& dfab = m_bottom_dot_mask[mfi];
            const auto& sfab = omask[mfi];
            amrex_mlndlap_set_dot_mask(BL_TO_FORTRAN_BOX(bx),
                                       BL_TO_FORTRAN_ANYD(dfab),
                                       BL_TO_FORTRAN_ANYD(sfab),
                                       BL_TO_FORTRAN_BOX(nddomain),
                                       m_lobc.data(), m_hibc.data());
        }
    }

    if (m_is_bottom_singular)
    {
        int amrlev = 0;
        int mglev = 0;
        const iMultiFab& omask = *m_owner_mask[amrlev][mglev];
        m_coarse_dot_mask.define(omask.boxArray(), omask.DistributionMap(), 1, 0);

        const Geometry& geom = m_geom[amrlev][mglev];
        Box nddomain = amrex::surroundingNodes(geom.Domain());

        if (m_coarsening_strategy != CoarseningStrategy::Sigma) {
            nddomain.grow(1000); // hack to avoid masks being modified at Neuman boundary
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(m_coarse_dot_mask,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto& dfab = m_coarse_dot_mask[mfi];
            const auto& sfab = omask[mfi];
            amrex_mlndlap_set_dot_mask(BL_TO_FORTRAN_BOX(bx),
                                       BL_TO_FORTRAN_ANYD(dfab),
                                       BL_TO_FORTRAN_ANYD(sfab),
                                       BL_TO_FORTRAN_BOX(nddomain),
                                       m_lobc.data(), m_hibc.data());
        }
    }
}

void
MLNodeLaplacian::buildStencil ()
{
    m_stencil.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_stencil[amrlev].resize(m_num_mg_levels[amrlev]);
    }
    
    if (m_coarsening_strategy != CoarseningStrategy::RAP) return;

    const int ncomp = (AMREX_SPACEDIM == 2) ? 5 : 15;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMREX_SPACEDIM != 1,
                                     "MLNodeLaplacian::buildStencil: 1d not supported");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!Geometry::IsRZ(),
                                     "MLNodeLaplacian::buildStencil: cylindrical not supported for RAP");

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_stencil[amrlev][mglev].reset
                (new MultiFab(amrex::convert(m_grids[amrlev][mglev],
                                             IntVect::TheNodeVector()),
                              m_dmap[amrlev][mglev], ncomp, 4));
            m_stencil[amrlev][mglev]->setVal(0.0);
        }

        {
            const Geometry& geom = m_geom[amrlev][0];
            const Real* dxinv = geom.InvCellSize();

#ifdef AMREX_USE_EB
            auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
            const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
            const MultiFab* conn = (factory) ? m_connection[amrlev].get() : nullptr;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                FArrayBox sgfab;
                FArrayBox cnfab;
                for (MFIter mfi(*m_stencil[amrlev][0], MFItInfo().EnableTiling().SetDynamic(true));
                     mfi.isValid(); ++mfi)
                {
                    Box vbx = mfi.validbox();
                    AMREX_D_TERM(vbx.growLo(0,1);, vbx.growLo(1,1);, vbx.growLo(2,1));
                    Box bx = mfi.growntilebox(1);
                    bx &= vbx;
                    const Box& ccbx = amrex::enclosedCells(bx);
                    FArrayBox& stfab = (*m_stencil[amrlev][0])[mfi];
                    const FArrayBox& sgfab_orig = (*m_sigma[amrlev][0][0])[mfi];
                    
#ifdef AMREX_USE_EB
                    bool regular = !factory;
                    if (factory)
                    {
                        const auto& flag = (*flags)[mfi];
                        const auto& ccbxg1 = amrex::grow(ccbx,1);
                        const auto& typ = flag.getType(ccbxg1);
                        if (typ == FabType::covered)
                        {
                            stfab.setVal(0.0, bx, 0, ncomp);
                        }
                        else if (typ == FabType::singlevalued)
                        {
                            const Box& bx2 = amrex::grow(ccbx,1);
                            sgfab.resize(bx2);
                            const FArrayBox& cnfab_orig = (*conn)[mfi];
                            cnfab.resize(bx2,cnfab_orig.nComp());

                            sgfab.setVal(0.0);
                            cnfab.setVal(0.0);

                            Box btmp = bx2 & sgfab_orig.box();
                            sgfab.copy(sgfab_orig, btmp, 0, btmp, 0, 1);

                            btmp = bx2 & cnfab_orig.box();
                            cnfab.copy(cnfab_orig, btmp, 0, btmp, 0, cnfab.nComp());

                            amrex_mlndlap_set_stencil_eb(BL_TO_FORTRAN_BOX(bx),
                                                         BL_TO_FORTRAN_ANYD(stfab),
                                                         BL_TO_FORTRAN_ANYD(sgfab),
                                                         BL_TO_FORTRAN_ANYD(cnfab),
                                                         dxinv);
                        }
                        else
                        {
                            regular = true;
                        }
                    }
                    if (regular)
#endif
                    {
                        Box bx2 = amrex::grow(ccbx,1);
                        sgfab.resize(bx2);
                        sgfab.setVal(0.0);
                        bx2 &= sgfab_orig.box();
                        sgfab.copy(sgfab_orig, bx2, 0, bx2, 0, 1);
                        amrex_mlndlap_set_stencil(BL_TO_FORTRAN_BOX(bx),
                                                  BL_TO_FORTRAN_ANYD(stfab),
                                                  BL_TO_FORTRAN_ANYD(sgfab),
                                                  dxinv);
                    }
                }

                for (MFIter mfi(*m_stencil[amrlev][0],true); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    FArrayBox& stfab = (*m_stencil[amrlev][0])[mfi];
                    
                    amrex_mlndlap_set_stencil_s0(BL_TO_FORTRAN_BOX(bx),
                                                 BL_TO_FORTRAN_ANYD(stfab));
                }
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
#pragma omp parallel
#endif
            {
                for (MFIter mfi(*pcrse, true); mfi.isValid(); ++mfi)
                {
                    Box vbx = mfi.validbox();
                    AMREX_D_TERM(vbx.growLo(0,1);, vbx.growLo(1,1);, vbx.growLo(2,1));
                    Box bx = mfi.growntilebox(1);
                    bx &= vbx;
                    if (m_use_simple_interp)
                    {
                        amrex_mlndlap_stencil_rap_sp(BL_TO_FORTRAN_BOX(bx),
                                                     BL_TO_FORTRAN_ANYD((*pcrse)[mfi]),
                                                     BL_TO_FORTRAN_ANYD(fine[mfi]));
                    }
                    else
                    {
                        amrex_mlndlap_stencil_rap(BL_TO_FORTRAN_BOX(bx),
                                                  BL_TO_FORTRAN_ANYD((*pcrse)[mfi]),
                                                  BL_TO_FORTRAN_ANYD(fine[mfi]));
                    }
                }

                for (MFIter mfi(*pcrse,true); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    FArrayBox& stfab = (*pcrse)[mfi];
                    
                    amrex_mlndlap_set_stencil_s0(BL_TO_FORTRAN_BOX(bx),
                                                 BL_TO_FORTRAN_ANYD(stfab));
                }
            }

            if (need_parallel_copy) {
                crse.ParallelCopy(cfine);
            }

            m_stencil[amrlev][mglev]->FillBoundary(m_geom[amrlev][mglev].periodicity());
        }
    }
}

void
MLNodeLaplacian::fixUpResidualMask (int amrlev, iMultiFab& resmsk)
{
    if (!m_masks_built) buildMasks();

    const iMultiFab& cfmask = *m_nd_fine_mask[amrlev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(resmsk,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        amrex_mlndlap_fixup_res_mask(BL_TO_FORTRAN_BOX(bx),
                                     BL_TO_FORTRAN_ANYD(resmsk[mfi]),
                                     BL_TO_FORTRAN_ANYD(cfmask[mfi]));
    }
}

void
MLNodeLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLNodeLaplacian::prepareForSolve()");

    MLNodeLinOp::prepareForSolve();

    buildMasks();

    averageDownCoeffs();

#if (AMREX_SPACEDIM == 2)
    amrex_mlndlap_set_rz(&m_is_rz);
#endif

#ifdef AMREX_USE_EB
    buildConnection();
#endif

    buildStencil();
}

void
MLNodeLaplacian::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
    BL_PROFILE("MLNodeLaplacian::restriction()");

    applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous);

    const Box& nd_domain = amrex::surroundingNodes(m_geom[amrlev][cmglev].Domain());

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
        cfine.define(ba, fine.DistributionMap(), 1, 0);
    }

    MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][cmglev-1];

    const auto& stencil = m_stencil[amrlev][cmglev-1];

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*pcrse, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        if (m_coarsening_strategy == CoarseningStrategy::Sigma || m_use_simple_interp)
        {
            amrex_mlndlap_restriction(BL_TO_FORTRAN_BOX(bx),
                                      BL_TO_FORTRAN_ANYD((*pcrse)[mfi]),
                                      BL_TO_FORTRAN_ANYD(fine[mfi]),
                                      BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                      BL_TO_FORTRAN_BOX(nd_domain),
                                      m_lobc.data(), m_hibc.data());
        }
        else
        {
            amrex_mlndlap_restriction_rap(BL_TO_FORTRAN_BOX(bx),
                                          BL_TO_FORTRAN_ANYD((*pcrse)[mfi]),
                                          BL_TO_FORTRAN_ANYD(fine[mfi]),
                                          BL_TO_FORTRAN_ANYD((*stencil)[mfi]),
                                          BL_TO_FORTRAN_ANYD(dmsk[mfi]));
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

    const Box& nd_domain = amrex::surroundingNodes(m_geom[amrlev][fmglev].Domain());

    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][fmglev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox tmpfab;
        for (MFIter mfi(fine, true); mfi.isValid(); ++mfi)
        {
            const Box& fbx = mfi.tilebox();
            const Box& cbx = amrex::coarsen(fbx,2);
            const Box& tmpbx = amrex::refine(cbx,2);
            tmpfab.resize(tmpbx);

            if (m_coarsening_strategy == CoarseningStrategy::RAP)
            {
                if (m_use_simple_interp)
                {
                    amrex_mlndlap_interpolation_rap_sp(BL_TO_FORTRAN_BOX(cbx),
                                                       BL_TO_FORTRAN_ANYD(tmpfab),
                                                       BL_TO_FORTRAN_ANYD((*cmf)[mfi]),
                                                       BL_TO_FORTRAN_ANYD((*stencil)[mfi]),
                                                       BL_TO_FORTRAN_ANYD(dmsk[mfi]));
                }
                else
                {
                    amrex_mlndlap_interpolation_rap(BL_TO_FORTRAN_BOX(cbx),
                                                    BL_TO_FORTRAN_ANYD(tmpfab),
                                                    BL_TO_FORTRAN_ANYD((*cmf)[mfi]),
                                                    BL_TO_FORTRAN_ANYD((*stencil)[mfi]),
                                                    BL_TO_FORTRAN_ANYD(dmsk[mfi]));
                }
            }
            else if (m_use_harmonic_average && fmglev > 0)
            {
                AMREX_D_TERM(const FArrayBox& sxfab = (*sigma[0])[mfi];,
                             const FArrayBox& syfab = (*sigma[1])[mfi];,
                             const FArrayBox& szfab = (*sigma[2])[mfi];);
                amrex_mlndlap_interpolation_ha(BL_TO_FORTRAN_BOX(cbx),
                                               BL_TO_FORTRAN_ANYD(tmpfab),
                                               BL_TO_FORTRAN_ANYD((*cmf)[mfi]),
                                               AMREX_D_DECL(BL_TO_FORTRAN_ANYD(sxfab),
                                                            BL_TO_FORTRAN_ANYD(syfab),
                                                            BL_TO_FORTRAN_ANYD(szfab)),
                                               BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                               BL_TO_FORTRAN_BOX(nd_domain),
                                               m_lobc.data(), m_hibc.data());
            }
            else
            {
                const FArrayBox& sfab = (*sigma[0])[mfi];
                amrex_mlndlap_interpolation_aa(BL_TO_FORTRAN_BOX(cbx),
                                               BL_TO_FORTRAN_ANYD(tmpfab),
                                               BL_TO_FORTRAN_ANYD((*cmf)[mfi]),
                                               BL_TO_FORTRAN_ANYD(sfab),
                                               BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                               BL_TO_FORTRAN_BOX(nd_domain),
                                               m_lobc.data(), m_hibc.data());
            }
            fine[mfi].plus(tmpfab,fbx,fbx,0,0,1);
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
    const Box& c_cc_domain = cgeom.Domain();
    const Box& c_nd_domain = amrex::surroundingNodes(c_cc_domain);

    const iMultiFab& fdmsk = *m_dirichlet_mask[camrlev+1][0];

    MultiFab cfine(amrex::coarsen(fba, 2), fdm, 1, 0);

    applyBC(camrlev+1, 0, *frhs, BCMode::Inhomogeneous);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(cfine, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        amrex_mlndlap_restriction(BL_TO_FORTRAN_BOX(bx),
                                  BL_TO_FORTRAN_ANYD(cfine[mfi]),
                                  BL_TO_FORTRAN_ANYD((*frhs)[mfi]),
                                  BL_TO_FORTRAN_ANYD(fdmsk[mfi]),
                                  BL_TO_FORTRAN_BOX(c_nd_domain),
                                  m_lobc.data(), m_hibc.data());
    }

    MultiFab tmp_crhs(crhs.boxArray(), crhs.DistributionMap(), 1, 0);
    tmp_crhs.setVal(0.0);
    tmp_crhs.ParallelCopy(cfine, cgeom.periodicity());

    const iMultiFab& c_nd_mask = *m_nd_fine_mask[camrlev];
    const auto& has_fine_bndry = *m_has_fine_bndry[camrlev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(crhs, MFItInfo().EnableTiling().SetDynamic(true));
         mfi.isValid(); ++mfi)
    {
        if (has_fine_bndry[mfi])
        {
            const Box& bx = mfi.tilebox();
            const auto& mfab = c_nd_mask[mfi];
            auto& dfab = crhs[mfi];
            const auto& sfab = tmp_crhs[mfi];
            amrex_mlndlap_copy_fine_node(BL_TO_FORTRAN_BOX(bx),
                                         BL_TO_FORTRAN_ANYD(dfab),
                                         BL_TO_FORTRAN_ANYD(sfab),
                                         BL_TO_FORTRAN_ANYD(mfab));
        }
    }
}

void
MLNodeLaplacian::applyBC (int amrlev, int mglev, MultiFab& phi, BCMode/* bc_mode*/,
                          bool skip_fillboundary) const
{
    BL_PROFILE("MLNodeLaplacian::applyBC()");

    const Geometry& geom = m_geom[amrlev][mglev];
    const Box& nd_domain = amrex::surroundingNodes(geom.Domain());

    if (!skip_fillboundary) {
        phi.FillBoundary(geom.periodicity());
    }

//    int inhom = (bc_mode == BCMode::Inhomogeneous);

    if (m_coarsening_strategy == CoarseningStrategy::Sigma)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(phi); mfi.isValid(); ++mfi)
        {
            if (!nd_domain.strictly_contains(mfi.fabbox()))
            {
                amrex_mlndlap_applybc(BL_TO_FORTRAN_ANYD(phi[mfi]),
                                      BL_TO_FORTRAN_BOX(nd_domain),
                                      m_lobc.data(), m_hibc.data());
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
    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

    const Box& domain_box = amrex::surroundingNodes(m_geom[amrlev][mglev].Domain());
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(out,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const FArrayBox& xfab = in[mfi];
        FArrayBox& yfab = out[mfi];

        if (m_coarsening_strategy == CoarseningStrategy::RAP)
        {
            amrex_mlndlap_adotx_sten(BL_TO_FORTRAN_BOX(bx),
                                     BL_TO_FORTRAN_ANYD(yfab),
                                     BL_TO_FORTRAN_ANYD(xfab),
                                     BL_TO_FORTRAN_ANYD((*stencil)[mfi]),
                                     BL_TO_FORTRAN_ANYD(dmsk[mfi]));
        }
        else if (m_use_harmonic_average && mglev > 0)
        {
            AMREX_D_TERM(const FArrayBox& sxfab = (*sigma[0])[mfi];,
                         const FArrayBox& syfab = (*sigma[1])[mfi];,
                         const FArrayBox& szfab = (*sigma[2])[mfi];);

            amrex_mlndlap_adotx_ha(BL_TO_FORTRAN_BOX(bx),
                                   BL_TO_FORTRAN_ANYD(yfab),
                                   BL_TO_FORTRAN_ANYD(xfab),
                                   AMREX_D_DECL(BL_TO_FORTRAN_ANYD(sxfab),
                                                BL_TO_FORTRAN_ANYD(syfab),
                                                BL_TO_FORTRAN_ANYD(szfab)),
                                   BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                   dxinv, BL_TO_FORTRAN_BOX(domain_box),
                                   m_lobc.data(), m_hibc.data());
        }
        else
        {
            const FArrayBox& sfab = (*sigma[0])[mfi];

            amrex_mlndlap_adotx_aa(BL_TO_FORTRAN_BOX(bx),
                                   BL_TO_FORTRAN_ANYD(yfab),
                                   BL_TO_FORTRAN_ANYD(xfab),
                                   BL_TO_FORTRAN_ANYD(sfab),
                                   BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                   dxinv, BL_TO_FORTRAN_BOX(domain_box),
                                   m_lobc.data(), m_hibc.data());
        }
    }
}

void
MLNodeLaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const
{
    BL_PROFILE("MLNodeLaplacian::Fsmooth()");

    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];

    if (m_use_gauss_seidel)
    {
        const auto& sigma = m_sigma[amrlev][mglev];
        const auto& stencil = m_stencil[amrlev][mglev];
        const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

        const Box& domain_box = amrex::surroundingNodes(m_geom[amrlev][mglev].Domain());

        if (m_coarsening_strategy == CoarseningStrategy::RAP)
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(sol); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                amrex_mlndlap_gauss_seidel_sten(BL_TO_FORTRAN_BOX(bx),
                                                BL_TO_FORTRAN_ANYD(sol[mfi]),
                                                BL_TO_FORTRAN_ANYD(rhs[mfi]),
                                                BL_TO_FORTRAN_ANYD((*stencil)[mfi]),
                                                BL_TO_FORTRAN_ANYD(dmsk[mfi]));
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
                AMREX_D_TERM(const FArrayBox& sxfab = (*sigma[0])[mfi];,
                             const FArrayBox& syfab = (*sigma[1])[mfi];,
                             const FArrayBox& szfab = (*sigma[2])[mfi];);

                amrex_mlndlap_gauss_seidel_ha(BL_TO_FORTRAN_BOX(bx),
                                              BL_TO_FORTRAN_ANYD(sol[mfi]),
                                              BL_TO_FORTRAN_ANYD(rhs[mfi]),
                                              AMREX_D_DECL(BL_TO_FORTRAN_ANYD(sxfab),
                                                           BL_TO_FORTRAN_ANYD(syfab),
                                                           BL_TO_FORTRAN_ANYD(szfab)),
                                              BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                              dxinv, BL_TO_FORTRAN_BOX(domain_box),
                                              m_lobc.data(), m_hibc.data());
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
                const FArrayBox& sfab = (*sigma[0])[mfi];

                amrex_mlndlap_gauss_seidel_aa(BL_TO_FORTRAN_BOX(bx),
                                              BL_TO_FORTRAN_ANYD(sol[mfi]),
                                              BL_TO_FORTRAN_ANYD(rhs[mfi]),
                                              BL_TO_FORTRAN_ANYD(sfab),
                                              BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                              dxinv, BL_TO_FORTRAN_BOX(domain_box),
                                              m_lobc.data(), m_hibc.data());
            }
        }

        nodalSync(amrlev, mglev, sol);
    }
    else
    {
        MultiFab Ax(sol.boxArray(), sol.DistributionMap(), 1, 0);
        Fapply(amrlev, mglev, Ax, sol);

        const auto& sigma = m_sigma[amrlev][mglev];
        const auto& stencil = m_stencil[amrlev][mglev];
        const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

        const Box& domain_box = amrex::surroundingNodes(m_geom[amrlev][mglev].Domain());

        if (m_coarsening_strategy == CoarseningStrategy::RAP)
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(sol); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.validbox();
                amrex_mlndlap_jacobi_sten(BL_TO_FORTRAN_BOX(bx),
                                          BL_TO_FORTRAN_ANYD(sol[mfi]),
                                          BL_TO_FORTRAN_ANYD(Ax[mfi]),
                                          BL_TO_FORTRAN_ANYD(rhs[mfi]),
                                          BL_TO_FORTRAN_ANYD((*stencil)[mfi]),
                                          BL_TO_FORTRAN_ANYD(dmsk[mfi]));
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
                AMREX_D_TERM(const FArrayBox& sxfab = (*sigma[0])[mfi];,
                             const FArrayBox& syfab = (*sigma[1])[mfi];,
                             const FArrayBox& szfab = (*sigma[2])[mfi];);

                amrex_mlndlap_jacobi_ha(BL_TO_FORTRAN_BOX(bx),
                                        BL_TO_FORTRAN_ANYD(sol[mfi]),
                                        BL_TO_FORTRAN_ANYD(Ax[mfi]),
                                        BL_TO_FORTRAN_ANYD(rhs[mfi]),
                                        AMREX_D_DECL(BL_TO_FORTRAN_ANYD(sxfab),
                                                     BL_TO_FORTRAN_ANYD(syfab),
                                                     BL_TO_FORTRAN_ANYD(szfab)),
                                        BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                        dxinv, BL_TO_FORTRAN_BOX(domain_box),
                                        m_lobc.data(), m_hibc.data());
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
                const FArrayBox& sfab = (*sigma[0])[mfi];

                amrex_mlndlap_jacobi_aa(BL_TO_FORTRAN_BOX(bx),
                                        BL_TO_FORTRAN_ANYD(sol[mfi]),
                                        BL_TO_FORTRAN_ANYD(Ax[mfi]),
                                        BL_TO_FORTRAN_ANYD(rhs[mfi]),
                                        BL_TO_FORTRAN_ANYD(sfab),
                                        BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                        dxinv, BL_TO_FORTRAN_BOX(domain_box),
                                        m_lobc.data(), m_hibc.data());
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
    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox& fab = mf[mfi];
        if (m_coarsening_strategy == CoarseningStrategy::RAP)
        {
            amrex_mlndlap_normalize_sten(BL_TO_FORTRAN_BOX(bx),
                                         BL_TO_FORTRAN_ANYD(fab),
                                         BL_TO_FORTRAN_ANYD((*stencil)[mfi]),
                                         BL_TO_FORTRAN_ANYD(dmsk[mfi]));
        }
        else if (m_use_harmonic_average && mglev > 0)
        {
            AMREX_D_TERM(const FArrayBox& sxfab = (*sigma[0])[mfi];,
                         const FArrayBox& syfab = (*sigma[1])[mfi];,
                         const FArrayBox& szfab = (*sigma[2])[mfi];);

            amrex_mlndlap_normalize_ha(BL_TO_FORTRAN_BOX(bx),
                                       BL_TO_FORTRAN_ANYD(fab),
                                       AMREX_D_DECL(BL_TO_FORTRAN_ANYD(sxfab),
                                                    BL_TO_FORTRAN_ANYD(syfab),
                                                    BL_TO_FORTRAN_ANYD(szfab)),
                                       BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                       dxinv);
        }
        else
        {
            const FArrayBox& sfab = (*sigma[0])[mfi];

            amrex_mlndlap_normalize_aa(BL_TO_FORTRAN_BOX(bx),
                                       BL_TO_FORTRAN_ANYD(fab),
                                       BL_TO_FORTRAN_ANYD(sfab),
                                       BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                       dxinv);
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
    const BoxArray& cc_fba = amrex::coarsen(fine_grids, ref_ratio);

    iMultiFab crse_cc_mask(ccba, dmap, 1, 1); // cell-center, 1: coarse; 0: covered by fine

    const int owner = 1;
    const int nonowner = 0;

    crse_cc_mask.setVal(owner);

    const Box& ccdom = geom.Domain();

    const std::vector<IntVect>& pshifts = geom.periodicity().shiftIntVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector< std::pair<int,Box> > isects;
        
        for (MFIter mfi(crse_cc_mask); mfi.isValid(); ++mfi)
        {
            IArrayBox& fab = crse_cc_mask[mfi];
            const Box& bx = fab.box();
            for (const auto& iv: pshifts)
            {
                cc_fba.intersections(bx+iv, isects);
                for (const auto& is : isects)
                {
                    fab.setVal(nonowner, is.second-iv, 0, 1);
                }
            }

            amrex_mlndlap_fillbc_cc_i(BL_TO_FORTRAN_ANYD(fab),
                                      BL_TO_FORTRAN_BOX(ccdom),
                                      m_lobc.data(), m_hibc.data());
        }
    }

    MultiFab phi(ndba, dmap, 1, 1);
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(phi, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Box& gbx = mfi.growntilebox();
        FArrayBox& fab = phi[mfi];
        fab.setVal(0.0, gbx);
        fab.copy(a_phi[mfi], bx, 0, bx, 0, 1);
        amrex_mlndlap_zero_fine(BL_TO_FORTRAN_BOX(bx),
                                BL_TO_FORTRAN_ANYD(fab),
                                BL_TO_FORTRAN_ANYD(crse_cc_mask[mfi]),
                                &nonowner);
    }

    const auto& nddom = amrex::surroundingNodes(ccdom);

    if (m_coarsening_strategy == CoarseningStrategy::Sigma)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(phi); mfi.isValid(); ++mfi)
        {
            if (!nddom.strictly_contains(mfi.fabbox()))
            {
                amrex_mlndlap_applybc(BL_TO_FORTRAN_ANYD(phi[mfi]),
                                      BL_TO_FORTRAN_BOX(nddom),
                                      m_lobc.data(), m_hibc.data());
            }
        }
    }

    const Real* dxinv = geom.InvCellSize();

    const MultiFab& sigma_orig = *m_sigma[0][0][0];
    const iMultiFab& dmsk = *m_dirichlet_mask[0][0];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox rhs, rhs2;
        FArrayBox u, sigma;
        for (MFIter mfi(sync_resid, MFItInfo().EnableTiling().SetDynamic(true)); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& bxg1 = amrex::grow(bx,1);
            const Box& ccbxg1 = amrex::enclosedCells(bxg1);
            if (amrex_mlndlap_any_fine_sync_cells(BL_TO_FORTRAN_BOX(ccbxg1),
                                                  BL_TO_FORTRAN_ANYD(crse_cc_mask[mfi]),
                                                  &nonowner))
            {
                const Box& ccvbx = amrex::enclosedCells(mfi.validbox());

                u.resize(ccbxg1, AMREX_SPACEDIM);
                u.setVal(0.0);
                Box b = ccbxg1 & ccvbx;
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                {
                    if (m_lobc[idim] == LinOpBCType::inflow)
                    {
                        if (b.smallEnd(idim) == ccdom.smallEnd(idim)) {
                            b.growLo(idim, 1);
                        }
                    }
                    if (m_hibc[idim] == LinOpBCType::inflow)
                    {
                        if (b.bigEnd(idim) == ccdom.bigEnd(idim)) {
                            b.growHi(idim, 1);
                        }
                    }
                }
                u.copy(vold[mfi], b, 0, b, 0, AMREX_SPACEDIM);

                amrex_fab_setval_ifnot(BL_TO_FORTRAN_BOX(ccbxg1),
                                       BL_TO_FORTRAN_FAB(u),
                                       BL_TO_FORTRAN_ANYD(crse_cc_mask[mfi]),
                                       0.0);

                rhs.resize(bx);
                amrex_mlndlap_divu(BL_TO_FORTRAN_BOX(bx),
                                   BL_TO_FORTRAN_ANYD(rhs),
                                   BL_TO_FORTRAN_ANYD(u),
                                   BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                   dxinv);

                if (rhcc)
                {
                    FArrayBox* rhcc_fab = &u;
                    rhcc_fab->resize(ccbxg1);
                    rhcc_fab->setVal(0.0);
                    const Box& b2 = ccbxg1 & ccvbx;
                    rhcc_fab->copy((*rhcc)[mfi], b2, 0, b2, 0, 1);

                    amrex_fab_setval_ifnot(BL_TO_FORTRAN_BOX(ccbxg1),
                                           BL_TO_FORTRAN_FAB(*rhcc_fab),
                                           BL_TO_FORTRAN_ANYD(crse_cc_mask[mfi]),
                                           0.0);                    

                    rhs2.resize(bx);
                    amrex_mlndlap_rhcc(BL_TO_FORTRAN_BOX(bx),
                                       BL_TO_FORTRAN_ANYD(rhs2),
                                       BL_TO_FORTRAN_ANYD(*rhcc_fab),
                                       BL_TO_FORTRAN_ANYD(dmsk[mfi]));
                    rhs.plus(rhs2);
                }

                sigma.resize(ccbxg1);
                sigma.setVal(0, ccbxg1, 0, 1);
                const Box& ibx = ccbxg1 & amrex::enclosedCells(mfi.validbox());
                sigma.copy(sigma_orig[mfi], ibx, 0, ibx, 0, 1);
                amrex_fab_setval_ifnot(BL_TO_FORTRAN_BOX(ccbxg1),
                                       BL_TO_FORTRAN_FAB(sigma),
                                       BL_TO_FORTRAN_ANYD(crse_cc_mask[mfi]),
                                       0.0);

                amrex_mlndlap_adotx_aa(BL_TO_FORTRAN_BOX(bx),
                                       BL_TO_FORTRAN_ANYD(sync_resid[mfi]),
                                       BL_TO_FORTRAN_ANYD(phi[mfi]),
                                       BL_TO_FORTRAN_ANYD(sigma),
                                       BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                       dxinv, BL_TO_FORTRAN_BOX(nddom),
                                       m_lobc.data(), m_hibc.data());

                amrex_mlndlap_crse_resid(BL_TO_FORTRAN_BOX(bx),
                                         BL_TO_FORTRAN_ANYD(sync_resid[mfi]),
                                         BL_TO_FORTRAN_ANYD(rhs),
                                         BL_TO_FORTRAN_ANYD(crse_cc_mask[mfi]),
                                         BL_TO_FORTRAN_BOX(nddom),
                                         m_lobc.data(), m_hibc.data());
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

    const Geometry& geom = m_geom[0][0];
    const Box& ccdom = geom.Domain();
    const auto& nddom = amrex::surroundingNodes(ccdom);

    const Real* dxinv = geom.InvCellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox rhs, rhs2;
        FArrayBox u;
        IArrayBox tmpmask;

        for (MFIter mfi(sync_resid, MFItInfo().EnableTiling().SetDynamic(true)); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& gbx = mfi.growntilebox();
            const Box& vbx = mfi.validbox();
            const Box& ccvbx = amrex::enclosedCells(vbx);
            const Box& bxg1 = amrex::grow(bx,1);
            const Box& ccbxg1 = amrex::enclosedCells(bxg1);

            u.resize(ccbxg1, AMREX_SPACEDIM);
            u.setVal(0.0, ccbxg1, 0, AMREX_SPACEDIM);
            Box ovlp = ccvbx & ccbxg1;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                if (m_lobc[idim] == LinOpBCType::inflow)
                {
                    if (ovlp.smallEnd(idim) == ccdom.smallEnd(idim)) {
                        ovlp.growLo(idim, 1);
                    }
                }
                if (m_hibc[idim] == LinOpBCType::inflow)
                {
                    if (ovlp.bigEnd(idim) == ccdom.bigEnd(idim)) {
                        ovlp.growHi(idim, 1);
                    }
                }
            }
            u.copy(vold[mfi], ovlp, 0, ovlp, 0, AMREX_SPACEDIM);

            tmpmask.resize(bx);
            tmpmask.copy(dmsk[mfi], bx, 0, bx, 0, 1);
            tmpmask -= 1;
            tmpmask *= -1;  //  0 in dmsk --> 1 in tmpmask, and 1 in dmsk --> 0 in tmpmask

            rhs.resize(bx);
            amrex_mlndlap_divu(BL_TO_FORTRAN_BOX(bx),
                               BL_TO_FORTRAN_ANYD(rhs),
                               BL_TO_FORTRAN_ANYD(u),
                               BL_TO_FORTRAN_ANYD(tmpmask),
                               dxinv);

            if (rhcc)
            {
                FArrayBox* rhcc_fab = &u;
                rhcc_fab->resize(ccbxg1);
                rhcc_fab->setVal(0.0);
                const Box& ovlp3 = ccvbx & ccbxg1;
                rhcc_fab->copy((*rhcc)[mfi], ovlp3, 0, ovlp3, 0, 1);

                rhs2.resize(bx);
                amrex_mlndlap_rhcc(BL_TO_FORTRAN_BOX(bx),
                                   BL_TO_FORTRAN_ANYD(rhs2),
                                   BL_TO_FORTRAN_ANYD(*rhcc_fab),
                                   BL_TO_FORTRAN_ANYD(tmpmask));
                rhs.plus(rhs2);
            }

            FArrayBox* sigma = &u;
            sigma->resize(ccbxg1);
            sigma->setVal(0.0, ccbxg1, 0, 1);
            const Box& ovlp2 = ccvbx & ccbxg1;
            sigma->copy(sigma_orig[mfi], ovlp2, 0, ovlp2, 0, 1);

            sync_resid[mfi].setVal(0.0, gbx, 0, 1);

            amrex_mlndlap_adotx_aa(BL_TO_FORTRAN_BOX(bx),
                                   BL_TO_FORTRAN_ANYD(sync_resid[mfi]),
                                   BL_TO_FORTRAN_ANYD(phi[mfi]),
                                   BL_TO_FORTRAN_ANYD(*sigma),
                                   BL_TO_FORTRAN_ANYD(tmpmask),
                                   dxinv, BL_TO_FORTRAN_BOX(nddom),
                                   m_lobc.data(), m_hibc.data());

            sync_resid[mfi].xpay(-1.0, rhs, bx, bx, 0, 0, 1);

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
    const Real* cdxinv = cgeom.InvCellSize();
    const Real* fdxinv = fgeom.InvCellSize();
    const Box& c_cc_domain = cgeom.Domain();
    const Box& c_nd_domain = amrex::surroundingNodes(c_cc_domain);

    const BoxArray& fba = fine_sol.boxArray();
    const DistributionMapping& fdm = fine_sol.DistributionMap();

    const iMultiFab& fdmsk = *m_dirichlet_mask[crse_amrlev+1][0];

    MultiFab fine_res_for_coarse(amrex::coarsen(fba, 2), fdm, 1, 0);

    applyBC(crse_amrlev+1, 0, fine_res, BCMode::Inhomogeneous);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(fine_res_for_coarse, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        amrex_mlndlap_restriction(BL_TO_FORTRAN_BOX(bx),
                                  BL_TO_FORTRAN_ANYD(fine_res_for_coarse[mfi]),
                                  BL_TO_FORTRAN_ANYD(fine_res[mfi]),
                                  BL_TO_FORTRAN_ANYD(fdmsk[mfi]),
                                  BL_TO_FORTRAN_BOX(c_nd_domain),
                                  m_lobc.data(), m_hibc.data());
    }
    res.ParallelCopy(fine_res_for_coarse, cgeom.periodicity());

    MultiFab fine_contrib(amrex::coarsen(fba, 2), fdm, 1, 0);
    fine_contrib.setVal(0.0);

    const auto& fsigma = *m_sigma[crse_amrlev+1][0][0];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox sigfab;
        FArrayBox Axfab;
        for (MFIter mfi(fine_contrib, MFItInfo().EnableTiling().SetDynamic(true));
             mfi.isValid(); ++mfi)
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
            sigfab.setVal(0.0);
            sigfab.copy(fsigma[mfi], b, 0, b, 0, 1);

            const Box& bx_Ax = amrex::grow(fbx,1);
            const Box& b2 = bx_Ax & amrex::grow(fvbx,-1);
            Axfab.resize(bx_Ax);
            Axfab.setVal(0.0);
            Axfab.copy(fine_rhs[mfi], b2, 0, b2, 0, 1);
            Axfab.minus(fine_res[mfi], b2, 0, 0, 1);

            amrex_mlndlap_res_fine_contrib(BL_TO_FORTRAN_BOX(cbx),
                                           BL_TO_FORTRAN_BOX(cvbx),
                                           BL_TO_FORTRAN_ANYD(fine_contrib[mfi]),
                                           BL_TO_FORTRAN_ANYD(fine_sol[mfi]),
                                           BL_TO_FORTRAN_ANYD(sigfab),
                                           BL_TO_FORTRAN_ANYD(Axfab),
                                           BL_TO_FORTRAN_ANYD(fdmsk[mfi]),
                                           fdxinv);
        }
    }

    MultiFab fine_contrib_on_crse(crse_sol.boxArray(), crse_sol.DistributionMap(), 1, 0);
    fine_contrib_on_crse.setVal(0.0);
    fine_contrib_on_crse.ParallelAdd(fine_contrib, cgeom.periodicity());

    const iMultiFab& cdmsk = *m_dirichlet_mask[crse_amrlev][0];
    const auto& nd_mask     = m_nd_fine_mask[crse_amrlev];
    const auto& cc_mask     = m_cc_fine_mask[crse_amrlev];
    const auto& has_fine_bndry = m_has_fine_bndry[crse_amrlev];

    const auto& csigma = *m_sigma[crse_amrlev][0][0];

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(res, MFItInfo().EnableTiling().SetDynamic(true)); mfi.isValid(); ++mfi)
    {
        if ((*has_fine_bndry)[mfi])
        {
            const Box& bx = mfi.tilebox();
            amrex_mlndlap_res_cf_contrib(BL_TO_FORTRAN_BOX(bx),
                                         BL_TO_FORTRAN_ANYD(res[mfi]),
                                         BL_TO_FORTRAN_ANYD(crse_sol[mfi]),
                                         BL_TO_FORTRAN_ANYD(crse_rhs[mfi]),
                                         BL_TO_FORTRAN_ANYD(csigma[mfi]),
                                         BL_TO_FORTRAN_ANYD(cdmsk[mfi]),
                                         BL_TO_FORTRAN_ANYD((*nd_mask)[mfi]),
                                         BL_TO_FORTRAN_ANYD((*cc_mask)[mfi]),
                                         BL_TO_FORTRAN_ANYD(fine_contrib_on_crse[mfi]),
                                         cdxinv, BL_TO_FORTRAN_BOX(c_nd_domain),
                                         m_lobc.data(), m_hibc.data());
        }
    }
}

#ifdef AMREX_USE_EB

void
MLNodeLaplacian::buildConnection ()
{
    if (m_connection_built) return;

    BL_PROFILE("MLNodeLaplacian::buildConnection()");

    m_connection_built = true;

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        auto& conn = *m_connection[amrlev];
        MultiFab* intg = m_integral[amrlev].get();
        const int ncomp = conn.nComp();

        auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        if (factory)
        {
            const auto& flags = factory->getMultiEBCellFlagFab();
            const auto& vfrac = factory->getVolFrac();
            const auto& area = factory->getAreaFrac();
            const auto& bcent = factory->getBndryCent();
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(conn, MFItInfo().EnableTiling().SetDynamic(true));
                 mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox();
                auto& cfab = conn[mfi];
                auto& gfab = (*intg)[mfi];                    
                const auto& flag = flags[mfi];
                auto typ = flag.getType(bx);
                
                if (typ == FabType::covered) {
                    cfab.setVal(0.0, bx, 0, ncomp);
                    gfab.setVal(0.0, bx, 0, AMREX_SPACEDIM);
                } else if (typ == FabType::singlevalued) {
                    amrex_mlndlap_set_connection(BL_TO_FORTRAN_BOX(bx),
                                                 BL_TO_FORTRAN_ANYD(cfab),
                                                 BL_TO_FORTRAN_ANYD(gfab),
                                                 BL_TO_FORTRAN_ANYD(flag),
                                                 BL_TO_FORTRAN_ANYD(vfrac[mfi]),
                                                 AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*area[0])[mfi]),
                                                              BL_TO_FORTRAN_ANYD((*area[1])[mfi]),
                                                              BL_TO_FORTRAN_ANYD((*area[2])[mfi])),
                                                 BL_TO_FORTRAN_ANYD(bcent[mfi]));
                }
            }
        }
    }
}

#endif

}
