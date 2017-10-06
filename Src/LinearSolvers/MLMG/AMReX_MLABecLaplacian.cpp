
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MultiFabUtil.H>

namespace amrex {

MLABecLaplacian::MLABecLaplacian (const Vector<Geometry>& a_geom,
                                  const Vector<BoxArray>& a_grids,
                                  const Vector<DistributionMapping>& a_dmap)
{
    define(a_geom, a_grids, a_dmap);
}

void
MLABecLaplacian::define (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap)
{
    MLLinOp::define(a_geom, a_grids, a_dmap);

    m_a_coeffs.resize(m_num_amr_levels);
    m_b_coeffs.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_a_coeffs[amrlev].resize(m_num_mg_levels[amrlev]);
        m_b_coeffs[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_a_coeffs[amrlev][mglev].define(m_grids[amrlev][mglev],
                                             m_dmap[amrlev][mglev],
                                             1, 0);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& ba = amrex::convert(m_grids[amrlev][mglev], IntVect::TheDimensionVector(idim));
                m_b_coeffs[amrlev][mglev][idim].define(ba,
                                                       m_dmap[amrlev][mglev],
                                                       1, 0);
            }
        }
    }
}

MLABecLaplacian::~MLABecLaplacian ()
{}

void
MLABecLaplacian::setACoeffs (int amrlev, const MultiFab& alpha)
{
    MultiFab::Copy(m_a_coeffs[amrlev][0], alpha, 0, 0, 1, 0);
}

void
MLABecLaplacian::setBCoeffs (int amrlev, int direction, const MultiFab& beta)
{
    MultiFab::Copy(m_b_coeffs[amrlev][0][direction], beta, 0, 0, 1, 0);
}

void
MLABecLaplacian::averageDownCoeffs ()
{
    for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
        auto& fine_a_coeffs = m_a_coeffs[amrlev];
        auto& fine_b_coeffs = m_b_coeffs[amrlev];
        auto& crse_a_coeffs = m_a_coeffs[amrlev-1][0];
        auto& crse_b_coeffs = m_b_coeffs[amrlev-1][0];
        auto& crse_geom     = m_geom    [amrlev-1][0];

        averageDownCoeffsSameAmrLevel(fine_a_coeffs, fine_b_coeffs);
        
        crse_a_coeffs.ParallelCopy(fine_a_coeffs.back());
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            crse_b_coeffs[idim].ParallelCopy(fine_b_coeffs.back()[idim], crse_geom.periodicity());
        }
    }

    averageDownCoeffsSameAmrLevel(m_a_coeffs[0], m_b_coeffs[0]);
}

void
MLABecLaplacian::averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a,
                                                Vector<std::array<MultiFab,AMREX_SPACEDIM> >& b)
{
    int nmglevs = a.size();
    for (int mglev = 1; mglev < nmglevs; ++mglev)
    {
        amrex::average_down(a[mglev-1], a[mglev], 0, 1, mg_coarsen_ratio);
        
        Vector<const MultiFab*> fine {AMREX_D_DECL(&(b[mglev-1][0]),
                                                   &(b[mglev-1][1]),
                                                   &(b[mglev-1][2]))};
        Vector<MultiFab*> crse {AMREX_D_DECL(&(b[mglev][0]),
                                             &(b[mglev][1]),
                                             &(b[mglev][2]))};
        IntVect ratio {AMREX_D_DECL(mg_coarsen_ratio, mg_coarsen_ratio, mg_coarsen_ratio)};
        amrex::average_down_faces(fine, crse, ratio, 0);
    }
}

}
