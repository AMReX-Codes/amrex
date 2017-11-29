
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLNodeLap_F.H>

namespace amrex {

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

    MLLinOp::define(a_geom, cc_grids, a_dmap, a_info);

    m_sigma.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_sigma.resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                auto iv = IntVect::TheUnitVector();
                iv[idim] = 0;
                const BoxArray& ba = amrex::convert(m_grids[amrlev][mglev], iv);
                m_sigma[amrlev][mglev][idim].define(ba, m_dmap[amrlev][mglev], 1, 1);
            }
        }
    }
}

void
MLNodeLaplacian::setSigma (int amrlev, const std::array<MultiFab const*,AMREX_SPACEDIM>& a_sigma)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        MultiFab::Copy(m_sigma[amrlev][0][idim], *a_sigma[idim], 0, 0, 1, 0);
    }
}


// average codes


void
MLNodeLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLNodeLaplacian::prepareForSolve()");

    MLLinOp::prepareForSolve();

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
