
#include <limits>

#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLNodeLap_F.H>

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
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                auto iv = IntVect::TheUnitVector();
                iv[idim] = 0;
                const BoxArray& ba = amrex::convert(m_grids[amrlev][mglev], iv);
                m_sigma[amrlev][mglev][idim].define(ba, m_dmap[amrlev][mglev], 1, 1);
                m_sigma[amrlev][mglev][idim].setVal(bogus_value);
            }
        }
    }
}

void
MLNodeLaplacian::setSigma (int amrlev, const MultiFab& a_sigma)
{
    auto& nd_sigma = m_sigma[amrlev][0];

    std::array<IntVect, AMREX_SPACEDIM> nodal{AMREX_D_DECL(nd_sigma[0].ixType().ixType(),
                                                           nd_sigma[1].ixType().ixType(),
                                                           nd_sigma[2].ixType().ixType())};
    std::array<IntVect, AMREX_SPACEDIM> ngrow{AMREX_D_DECL(IntVect::TheDimensionVector(0),
                                                           IntVect::TheDimensionVector(1),
                                                           IntVect::TheDimensionVector(2))};

    MultiFab cc_sigma(a_sigma.boxArray(), a_sigma.DistributionMap(), 1, 1);
    cc_sigma.setVal(bogus_value);
    MultiFab::Copy(cc_sigma, a_sigma, 0, 0, 1, 0);
    cc_sigma.FillBoundary(m_geom[amrlev][0].periodicity());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(cc_sigma, true); mfi.isValid(); ++mfi)
    {
        AMREX_D_TERM(const Box& xbx = mfi.tilebox(nodal[0],ngrow[0]);,
                     const Box& ybx = mfi.tilebox(nodal[1],ngrow[1]);,
                     const Box& zbx = mfi.tilebox(nodal[2],ngrow[2]););
        AMREX_D_TERM(FArrayBox& xfab = nd_sigma[0][mfi];,
                     FArrayBox& yfab = nd_sigma[1][mfi];,
                     FArrayBox& zfab = nd_sigma[2][mfi];);
        const FArrayBox& ccfab = cc_sigma[mfi];

        amrex_mlndlap_sigma_cctoedge(AMREX_D_DECL(BL_TO_FORTRAN_BOX(xbx),
                                                  BL_TO_FORTRAN_BOX(ybx),
                                                  BL_TO_FORTRAN_BOX(zbx)),
                                     AMREX_D_DECL(BL_TO_FORTRAN_ANYD(xfab),
                                                  BL_TO_FORTRAN_ANYD(yfab),
                                                  BL_TO_FORTRAN_ANYD(zfab)),
                                     BL_TO_FORTRAN_ANYD(ccfab));
    }
}

void
MLNodeLaplacian::compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel,
                          const Vector<const MultiFab*>& rhnd,
                          const Vector<const MultiFab*>& rhcc) const
{
    AMREX_ALWAYS_ASSERT(rhs.size() == 1);
    AMREX_ALWAYS_ASSERT(rhcc[0] == nullptr);

    for (int ilev = 0; ilev < m_num_amr_levels; ++ilev)
    {
        vel[ilev]->FillBoundary(m_geom[ilev][0].periodicity());

        const Real* dxinv = m_geom[ilev][0].InvCellSize();

        const auto& nddom = amrex::surroundingNodes(m_geom[ilev][0].Domain());

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*rhs[ilev], true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            amrex_mlndlap_divu(BL_TO_FORTRAN_BOX(bx),                               
                               BL_TO_FORTRAN_ANYD((*rhs[ilev])[mfi]),
                               BL_TO_FORTRAN_ANYD((*vel[ilev])[mfi]),
                               dxinv, BL_TO_FORTRAN_BOX(nddom),
                               m_lobc.data(), m_hibc.data());
        }

        if (rhnd[ilev]) {
            MultiFab::Copy(*rhs[ilev], *rhnd[ilev], 0, 0, 1, 0);
        }
    }
}

// average codes


void
MLNodeLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLNodeLaplacian::prepareForSolve()");

    MLNodeLinOp::prepareForSolve();

#if (AMREX_SPACEDIM != 3)
    // applyMetricTermsCoeffs();
#endif

    // averageDownCoeffs();

    m_is_singular.clear();
    m_is_singular.resize(m_num_amr_levels, false);
    auto itlo = std::find(m_lobc.begin(), m_lobc.end(), BCType::Dirichlet);
    auto ithi = std::find(m_hibc.begin(), m_hibc.end(), BCType::Dirichlet);
    if (itlo == m_lobc.end() && ithi == m_hibc.end())
    {  // No Dirichlet
        for (int alev = 0; alev < m_num_amr_levels; ++alev)
        {
            if (m_domain_covered[alev])
            {
                m_is_singular[alev] = true;
            }    
        }
    }
}

void
MLNodeLaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{

}

void
MLNodeLaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const
{

}

void
MLNodeLaplacian::FFlux (int amrlev, const MFIter& mfi,
                        const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only) const
{

}

Real
MLNodeLaplacian::Anorm (int amrlev, int mglev) const
{

}

}
