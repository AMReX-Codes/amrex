
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
                                  const LPInfo& a_info)
{
    define(a_geom, a_grids, a_dmap, a_info);
}

MLNodeLaplacian::~MLNodeLaplacian ()
{}

void
MLNodeLaplacian::define (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info)
{
    BL_PROFILE("MLNodeLaplacian::define()");

    // This makes sure grids are cell-centered;
    Vector<BoxArray> cc_grids = a_grids;
    for (auto& ba : cc_grids) {
        ba.enclosedCells();
    }

    MLNodeLinOp::define(a_geom, cc_grids, a_dmap, a_info);

    m_sigma.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_sigma[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            int idim = 0;
            m_sigma[amrlev][mglev][idim].reset
                (new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1));
        }
    }
}

void
MLNodeLaplacian::setSigma (int amrlev, const MultiFab& a_sigma)
{
    MultiFab::Copy(*m_sigma[amrlev][0][0], a_sigma, 0, 0, 1, 0);
}

void
MLNodeLaplacian::compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel,
                          const Vector<const MultiFab*>& rhnd,
                          const Vector<MultiFab*>& rhcc)
{
    if (!m_masks_built) buildMasks();

    for (int ilev = 0; ilev < m_num_amr_levels; ++ilev)
    {
        const Geometry& geom = m_geom[ilev][0];
        vel[ilev]->FillBoundary(geom.periodicity());

        const Real* dxinv = geom.InvCellSize();
        const Box& nddom = amrex::surroundingNodes(geom.Domain());

        const iMultiFab& dmsk = *m_dirichlet_mask[ilev][0];

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*rhs[ilev], true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            amrex_mlndlap_divu(BL_TO_FORTRAN_BOX(bx),
                               BL_TO_FORTRAN_ANYD((*rhs[ilev])[mfi]),
                               BL_TO_FORTRAN_ANYD((*vel[ilev])[mfi]),
                               BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                               dxinv, BL_TO_FORTRAN_BOX(nddom),
                               m_lobc.data(), m_hibc.data());
        }
    }

    for (int ilev = m_num_amr_levels-2; ilev >= 0; --ilev)
    {
        const Geometry& fgeom = m_geom[ilev+1][0];
        const Geometry& cgeom = m_geom[ilev  ][0];

        MultiFab frhs(rhs[ilev+1]->boxArray(), rhs[ilev+1]->DistributionMap(), 1, 1);
        frhs.setVal(0.0);
        MultiFab::Copy(frhs, *rhs[ilev+1], 0, 0, 1, 0);
        frhs.FillBoundary(fgeom.periodicity());

        MultiFab crhs(amrex::coarsen(frhs.boxArray(),2), frhs.DistributionMap(), 1, 0);

        const iMultiFab& fdmsk = *m_dirichlet_mask[ilev+1][0];
        const Box& cnddom = amrex::surroundingNodes(cgeom.Domain());
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(crhs, true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            amrex_mlndlap_restriction(BL_TO_FORTRAN_BOX(bx),
                                      BL_TO_FORTRAN_ANYD(crhs[mfi]),
                                      BL_TO_FORTRAN_ANYD(frhs[mfi]),
                                      BL_TO_FORTRAN_ANYD(fdmsk[mfi]),
                                      BL_TO_FORTRAN_BOX(cnddom),
                                      m_lobc.data(), m_hibc.data());
        }

        rhs[ilev]->ParallelCopy(crhs, cgeom.periodicity());
    }

    for (int ilev = m_num_amr_levels-2; ilev >= 0; --ilev)
    {
        const Geometry& cgeom = m_geom[ilev  ][0];
        const Geometry& fgeom = m_geom[ilev+1][0];

        MultiFab frhs(amrex::coarsen(rhs[ilev+1]->boxArray(),2),
                      rhs[ilev+1]->DistributionMap(), 1, 0);
        frhs.setVal(0.0);

        const Box& fnddom = amrex::surroundingNodes(fgeom.Domain());
        const Real* fdxinv = fgeom.InvCellSize();
        const iMultiFab& fdmsk = *m_dirichlet_mask[ilev+1][0];

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(frhs, MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
        {
            const Box& cbx = mfi.validbox();
            const Box& fbx = amrex::refine(cbx,2);
            amrex_mlndlap_divu_fine_contrib(BL_TO_FORTRAN_BOX(cbx),
                                            BL_TO_FORTRAN_BOX(fbx),
                                            BL_TO_FORTRAN_ANYD(frhs[mfi]),
                                            BL_TO_FORTRAN_ANYD((*vel[ilev+1])[mfi]),
                                            BL_TO_FORTRAN_ANYD((*rhs[ilev+1])[mfi]),
                                            BL_TO_FORTRAN_ANYD(fdmsk[mfi]),
                                            fdxinv, BL_TO_FORTRAN_BOX(fnddom),
                                            m_lobc.data(), m_hibc.data());
        }

        MultiFab crhs(rhs[ilev]->boxArray(), rhs[ilev]->DistributionMap(), 1, 0);
        crhs.setVal(0.0);
        crhs.ParallelAdd(frhs, cgeom.periodicity());

        const Box& cnddom = amrex::surroundingNodes(cgeom.Domain());
        const Real* cdxinv = cgeom.InvCellSize();
        const iMultiFab& cdmsk = *m_dirichlet_mask[ilev][0];
        const iMultiFab& cfmask = *m_crsefine_mask[ilev];
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
                amrex_mlndlap_divu_cf_contrib(BL_TO_FORTRAN_BOX(bx),
                                              BL_TO_FORTRAN_ANYD((*rhs[ilev])[mfi]),
                                              BL_TO_FORTRAN_ANYD((*vel[ilev])[mfi]),
                                              BL_TO_FORTRAN_ANYD(cdmsk[mfi]),
                                              BL_TO_FORTRAN_ANYD(cfmask[mfi]),
                                              BL_TO_FORTRAN_ANYD(crhs[mfi]),
                                              cdxinv, BL_TO_FORTRAN_BOX(cnddom),
                                              m_lobc.data(), m_hibc.data());
            }
        }
    }    

    for (int ilev = 0; ilev < m_num_amr_levels; ++ilev)
    {
        if (rhcc[ilev])
        {
            rhcc[ilev]->FillBoundary(m_geom[ilev][0].periodicity());

            const Box& ccdom = m_geom[ilev][0].Domain();
            const iMultiFab& dmsk = *m_dirichlet_mask[ilev][0];

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(*rhcc[ilev], MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
            {
                if (!ccdom.contains(mfi.fabbox()))
                {
                    amrex_mlndlap_fillbc_cc(BL_TO_FORTRAN_ANYD((*rhcc[ilev])[mfi]),
                                            BL_TO_FORTRAN_BOX(ccdom),
                                            m_lobc.data(), m_hibc.data());
                }
            }

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(*rhs[ilev], true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                amrex_mlndlap_add_rhcc(BL_TO_FORTRAN_BOX(bx),
                                       BL_TO_FORTRAN_ANYD((*rhs[ilev])[mfi]),
                                       BL_TO_FORTRAN_ANYD((*rhcc[ilev])[mfi]),
                                       BL_TO_FORTRAN_ANYD(dmsk[mfi]));
            }
        }

        if (rhnd[ilev]) {
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
        for (MFIter mfi(*vel[amrlev], true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            amrex_mlndlap_mknewu(BL_TO_FORTRAN_BOX(bx),
                                 BL_TO_FORTRAN_ANYD((*vel[amrlev])[mfi]),
                                 BL_TO_FORTRAN_ANYD((*sol[amrlev])[mfi]),
                                 BL_TO_FORTRAN_ANYD(sigma[mfi]),
                                 dxinv);
        }
    }
}

void
MLNodeLaplacian::averageDownCoeffs ()
{
    BL_PROFILE("MLNodeLaplacian::averageDownCoeffs()");

    if (m_use_harmonic_average)
    {
        for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
        {
            for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
            {
                for (int idim = 1; idim < AMREX_SPACEDIM; ++idim)
                {
                    if (m_sigma[amrlev][mglev][idim] == nullptr) {
                        if (mglev == 0) {
                            m_sigma[amrlev][mglev][idim].reset
                                (new MultiFab(*m_sigma[amrlev][mglev][0], amrex::make_alias, 0, 1));
                        } else {
                            m_sigma[amrlev][mglev][idim].reset
                                (new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1));
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
                    FillBoundaryCoeff(*m_sigma[amrlev][mglev][idim], m_geom[amrlev][mglev]);
                }
            }
        } else {
            int idim = 0;
            for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
            {
                FillBoundaryCoeff(*m_sigma[amrlev][mglev][idim], m_geom[amrlev][mglev]);
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
    sigma.FillBoundary(geom.periodicity());

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

void
MLNodeLaplacian::buildMasks ()
{
    if (m_masks_built) return;

    m_masks_built = true;

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
                const Box& nddomain = amrex::surroundingNodes(geom.Domain());
                const std::vector<IntVect>& pshifts = period.shiftIntVect();

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
        iMultiFab& mask = *m_crsefine_mask[amrlev];
        MultiFab& weight = *m_crsefine_weight[amrlev];
        LayoutData<int>& has_cf = *m_has_fine_bndry[amrlev];
        const BoxArray& fba = m_grids[amrlev+1][0];
        const BoxArray& cfba = amrex::coarsen(fba, AMRRefRatio(amrlev));

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMRRefRatio(amrlev) == 2, "ref_ratio != 0 not supported");

        mask.setVal(0);  // coarse by default

        const std::vector<IntVect>& pshifts = m_geom[amrlev][0].periodicity().shiftIntVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            std::vector< std::pair<int,Box> > isects;

            for (MFIter mfi(mask); mfi.isValid(); ++mfi)
            {
                has_cf[mfi] = 0;
                IArrayBox& fab = mask[mfi];
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
            }
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(weight,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            amrex_mlndlap_set_cf_weight(BL_TO_FORTRAN_BOX(bx),
                                        BL_TO_FORTRAN_ANYD(weight[mfi]),
                                        BL_TO_FORTRAN_ANYD(mask[mfi]));
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
}

void
MLNodeLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLNodeLaplacian::prepareForSolve()");

    MLNodeLinOp::prepareForSolve();

    averageDownCoeffs();

    buildMasks();
}

void
MLNodeLaplacian::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
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

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*pcrse, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        amrex_mlndlap_restriction(BL_TO_FORTRAN_BOX(bx),
                                  BL_TO_FORTRAN_ANYD((*pcrse)[mfi]),
                                  BL_TO_FORTRAN_ANYD(fine[mfi]),
                                  BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                  BL_TO_FORTRAN_BOX(nd_domain),
                                  m_lobc.data(), m_hibc.data());
    }

    if (need_parallel_copy) {
        crse.ParallelCopy(cfine);
    }
}

void
MLNodeLaplacian::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
    const auto& sigma = m_sigma[amrlev][fmglev];

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
            if (m_use_harmonic_average && fmglev > 0)
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
    const auto amrrr = AMRRefRatio(camrlev);
    amrex::average_down(fine_sol, crse_sol, 0, 1, amrrr);
}

void
MLNodeLaplacian::applyBC (int amrlev, int mglev, MultiFab& phi, BCMode/* bc_mode*/,
                          bool skip_fillboundary) const
{
    const Geometry& geom = m_geom[amrlev][mglev];
    const Box& nd_domain = amrex::surroundingNodes(geom.Domain());

    if (!skip_fillboundary) {
        phi.FillBoundary(geom.periodicity());
    }

//    int inhom = (bc_mode == BCMode::Inhomogeneous);

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

void
MLNodeLaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    const auto& sigma = m_sigma[amrlev][mglev];
    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

    const Box& domain_box = amrex::surroundingNodes(m_geom[amrlev][mglev].Domain());
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox dg;
        for (MFIter mfi(out,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const FArrayBox& xfab = in[mfi];
            FArrayBox& yfab = out[mfi];
            const Box& bxg1 = amrex::grow(bx,1);
            dg.resize(bxg1,AMREX_SPACEDIM);
            if (m_use_harmonic_average && mglev > 0)
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
                                       BL_TO_FORTRAN_ANYD(dg),
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
                                       BL_TO_FORTRAN_ANYD(dg),
                                       BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                       dxinv, BL_TO_FORTRAN_BOX(domain_box),
                                       m_lobc.data(), m_hibc.data());
            }
        }
    }
}

void
MLNodeLaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const
{
    const iMultiFab& dmsk = *m_dirichlet_mask[amrlev][mglev];

    if (m_use_gauss_seidel)
    {
        const auto& sigma = m_sigma[amrlev][mglev];
        const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

        const Box& domain_box = amrex::surroundingNodes(m_geom[amrlev][mglev].Domain());

        if (m_use_harmonic_average && mglev > 0)
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
        const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

        const Box& domain_box = amrex::surroundingNodes(m_geom[amrlev][mglev].Domain());

        if (m_use_harmonic_average && mglev > 0)
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
MLNodeLaplacian::compSyncResidualCoarse (MultiFab& sync_resid, const MultiFab& a_phi,
                                         const MultiFab& vold, const BoxArray& fine_grids,
                                         const IntVect& ref_ratio)
{
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
                                BL_TO_FORTRAN_ANYD(crse_cc_mask[mfi]));
    }

    const auto& nddom = amrex::surroundingNodes(geom.Domain());

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

    MultiFab u(ccba, dmap, AMREX_SPACEDIM, 1);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(u,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = u[mfi];
        const IArrayBox& ifab = crse_cc_mask[mfi];
        const Box& bx = mfi.tilebox();
        const Box& gbx = mfi.growntilebox();
        fab.copy(vold[mfi], gbx, 0, gbx, 0, AMREX_SPACEDIM);
        amrex_fab_setval_ifnot (BL_TO_FORTRAN_BOX(bx),
                                BL_TO_FORTRAN_FAB(fab),
                                BL_TO_FORTRAN_ANYD(ifab),
                                0.0);
    }
    
    const Real* dxinv = geom.InvCellSize();

    const MultiFab& sigma_orig = *m_sigma[0][0][0];
    const iMultiFab& dmsk = *m_dirichlet_mask[0][0];

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox rhs, sigma, dg;
        for (MFIter mfi(sync_resid, MFItInfo().EnableTiling().SetDynamic(true)); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& bxg1 = amrex::grow(bx,1);
            const Box& ccbxg1 = amrex::enclosedCells(bxg1);
            if (amrex_mlndlap_any_zero(BL_TO_FORTRAN_BOX(bx),
                                       BL_TO_FORTRAN_ANYD(crse_cc_mask[mfi])))
            {
                rhs.resize(bx);
                amrex_mlndlap_divu(BL_TO_FORTRAN_BOX(bx),
                                   BL_TO_FORTRAN_ANYD(rhs),
                                   BL_TO_FORTRAN_ANYD(u[mfi]),
                                   BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                   dxinv, BL_TO_FORTRAN_BOX(nddom),
                                   m_lobc.data(), m_hibc.data());

                sigma.resize(ccbxg1);
                sigma.setVal(0, ccbxg1, 0, 1);
                const Box& ibx = ccbxg1 & amrex::enclosedCells(mfi.validbox());
                sigma.copy(sigma_orig[mfi], ibx, 0, ibx, 0, 1);
                amrex_fab_setval_ifnot(BL_TO_FORTRAN_BOX(ccbxg1),
                                       BL_TO_FORTRAN_FAB(sigma),
                                       BL_TO_FORTRAN_ANYD(crse_cc_mask[mfi]),
                                       0.0);

                dg.resize(bxg1,AMREX_SPACEDIM);
                amrex_mlndlap_adotx_aa(BL_TO_FORTRAN_BOX(bx),
                                       BL_TO_FORTRAN_ANYD(sync_resid[mfi]),
                                       BL_TO_FORTRAN_ANYD(phi[mfi]),
                                       BL_TO_FORTRAN_ANYD(sigma),
                                       BL_TO_FORTRAN_ANYD(dg),
                                       BL_TO_FORTRAN_ANYD(dmsk[mfi]),
                                       dxinv, BL_TO_FORTRAN_BOX(nddom),
                                       m_lobc.data(), m_hibc.data());

                amrex_mlndlap_crse_resid(BL_TO_FORTRAN_BOX(bx),
                                         BL_TO_FORTRAN_ANYD(sync_resid[mfi]),
                                         BL_TO_FORTRAN_ANYD(rhs),
                                         BL_TO_FORTRAN_ANYD(crse_cc_mask[mfi]));
            }
        }
    }
}

void
MLNodeLaplacian::compSyncResidualFine (MultiFab& sync_resid, const MultiFab& phi, const MultiFab& vold)
{
    const MultiFab& sigma_orig = *m_sigma[0][0][0];
    const iMultiFab& dmsk = *m_dirichlet_mask[0][0];

    const Geometry& geom = m_geom[0][0];
    const auto& nddom = amrex::surroundingNodes(geom.Domain());
    const auto& fake_nddom = amrex::grow(nddom, 1);

    const Real* dxinv = geom.InvCellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox rhs, u, sigma, dg;
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
            const Box& ovlp = ccvbx & ccbxg1;
            u.copy(vold[mfi], ovlp, 0, ovlp, 0, AMREX_SPACEDIM);

            tmpmask.resize(bx);
            tmpmask.copy(dmsk[mfi], bx, 0, bx, 0, 1);
            tmpmask -= 1;

            rhs.resize(bx);
            amrex_mlndlap_divu(BL_TO_FORTRAN_BOX(bx),
                               BL_TO_FORTRAN_ANYD(rhs),
                               BL_TO_FORTRAN_ANYD(u),
                               BL_TO_FORTRAN_ANYD(tmpmask),
                               dxinv, BL_TO_FORTRAN_BOX(fake_nddom),
                               m_lobc.data(), m_hibc.data());

            sigma.resize(ccbxg1);
            sigma.setVal(0.0, ccbxg1, 0, 1);
            sigma.copy(sigma_orig[mfi], ovlp, 0, ovlp, 0, 1);

            dg.resize(bxg1, AMREX_SPACEDIM);

            sync_resid[mfi].setVal(0.0, gbx, 0, 1);

            // What do we do at physical boundary?
            amrex_mlndlap_adotx_aa(BL_TO_FORTRAN_BOX(bx),
                                   BL_TO_FORTRAN_ANYD(sync_resid[mfi]),
                                   BL_TO_FORTRAN_ANYD(phi[mfi]),
                                   BL_TO_FORTRAN_ANYD(sigma),
                                   BL_TO_FORTRAN_ANYD(dg),
                                   BL_TO_FORTRAN_ANYD(tmpmask),
                                   dxinv, BL_TO_FORTRAN_BOX(nddom),
                                   m_lobc.data(), m_hibc.data());

            sync_resid[mfi].xpay(-1.0, rhs, bx, bx, 0, 0, 1);
        }
    }
}

void
MLNodeLaplacian::reflux (int crse_amrlev, MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                         MultiFab& fine_res, MultiFab& fine_sol) const
{
    const Geometry& cgeom = m_geom[crse_amrlev  ][0];
    const Geometry& fgeom = m_geom[crse_amrlev+1][0];
    const Real* cdxinv = cgeom.InvCellSize();
    const Real* fdxinv = fgeom.InvCellSize();
    const Box& c_nd_domain = amrex::surroundingNodes(cgeom.Domain());

    const BoxArray& fba = fine_sol.boxArray();
    const DistributionMapping& fdm = fine_sol.DistributionMap();

    const iMultiFab& fdmsk = *m_dirichlet_mask[crse_amrlev+1][0];

    MultiFab fine_res_for_coarse(amrex::coarsen(fba, 2), fdm, 1, 0);
    fine_res.FillBoundary(fgeom.periodicity());

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
    const auto& ovmsk  = *m_overlap_mask[crse_amrlev+1];

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(fine_contrib, MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
    {
        const Box& cbx = mfi.validbox();
        const Box& fbx = amrex::refine(cbx,2);
        amrex_mlndlap_res_fine_contrib(BL_TO_FORTRAN_BOX(cbx),
                                       BL_TO_FORTRAN_BOX(fbx),
                                       BL_TO_FORTRAN_ANYD(fine_contrib[mfi]),
                                       BL_TO_FORTRAN_ANYD(fine_sol[mfi]),
                                       BL_TO_FORTRAN_ANYD(fsigma[mfi]),
                                       BL_TO_FORTRAN_ANYD(fine_res[mfi]),
                                       BL_TO_FORTRAN_ANYD(fdmsk[mfi]),
                                       BL_TO_FORTRAN_ANYD(ovmsk[mfi]),
                                       fdxinv);
    }

    MultiFab fine_contrib_on_crse(crse_sol.boxArray(), crse_sol.DistributionMap(), 1, 0);
    fine_contrib_on_crse.setVal(0.0);
    fine_contrib_on_crse.ParallelAdd(fine_contrib, cgeom.periodicity());

    const iMultiFab& cdmsk = *m_dirichlet_mask[crse_amrlev][0];
    const auto& cfmask     = m_crsefine_mask[crse_amrlev];
    const auto& cfweight   = m_crsefine_weight[crse_amrlev];
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
                                         BL_TO_FORTRAN_ANYD((*cfmask)[mfi]),
                                         BL_TO_FORTRAN_ANYD((*cfweight)[mfi]),
                                         BL_TO_FORTRAN_ANYD(fine_contrib_on_crse[mfi]),
                                         cdxinv);
        }
    }
}

}
