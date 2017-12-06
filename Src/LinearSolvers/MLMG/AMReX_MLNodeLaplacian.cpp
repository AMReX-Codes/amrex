
#include <limits>

#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLNodeLap_F.H>
#include <AMReX_MultiFabUtil.H>

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
            if (mglev == 0) {
                int idim = 0;
                m_sigma[amrlev][mglev][idim].reset
                    (new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1));
                m_sigma[amrlev][mglev][idim]->setVal(0.0);
                for (idim = 1; idim < AMREX_SPACEDIM; ++idim)
                {
                    m_sigma[amrlev][mglev][idim].reset
                        (new MultiFab(*m_sigma[amrlev][mglev][0], amrex::make_alias, 0, 1));
                }
            } else {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                {
                    m_sigma[amrlev][mglev][idim].reset
                        (new MultiFab(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 1));
                    m_sigma[amrlev][mglev][idim]->setVal(0.0);
                }                
            }
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

void
MLNodeLaplacian::averageDownCoeffs ()
{
    BL_PROFILE("MLNodeLaplacian::averageDownCoeffs()");

    for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
        averageDownCoeffsSameAmrLevel(amrlev);
        averageDownCoeffsToCoarseAmrLevel(amrlev);
    }

    averageDownCoeffsSameAmrLevel(0);

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        int mglev = 0;
        m_sigma[amrlev][mglev][0]->FillBoundary(m_geom[amrlev][mglev].periodicity());
        for (mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_sigma[amrlev][mglev][idim]->FillBoundary(m_geom[amrlev][mglev].periodicity());
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
    for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const MultiFab& fine = *m_sigma[amrlev][mglev-1][idim];
            MultiFab& crse = *m_sigma[amrlev][mglev][idim];
            if (amrex::isMFIterSafe(crse, fine))
            {
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (MFIter mfi(crse, true); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    amrex_mlndlap_avgdown_coeff(BL_TO_FORTRAN_BOX(bx),
                                                BL_TO_FORTRAN_ANYD(crse[mfi]),
                                                BL_TO_FORTRAN_ANYD(fine[mfi]),
                                                &idim);
                }
            }
            else
            {
                BoxArray ba = fine.boxArray();
                ba.coarsen(2);
                MultiFab tmp(ba, fine.DistributionMap(), 1, 0);
#ifdef _OPENMP
#pragma omp parallel
#endif
                for (MFIter mfi(tmp, true); mfi.isValid(); ++mfi)
                {
                    const Box& bx = mfi.tilebox();
                    amrex_mlndlap_avgdown_coeff(BL_TO_FORTRAN_BOX(bx),
                                                BL_TO_FORTRAN_ANYD(tmp[mfi]),
                                                BL_TO_FORTRAN_ANYD(fine[mfi]),
                                                &idim);
                }
                crse.ParallelCopy(tmp);
            }
        }
    }
}

void
MLNodeLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLNodeLaplacian::prepareForSolve()");

    MLNodeLinOp::prepareForSolve();

    averageDownCoeffs();
}

void
MLNodeLaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    const auto& sigma = m_sigma[amrlev][mglev];
    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();
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
            AMREX_D_TERM(const FArrayBox& sxfab = (*sigma[0])[mfi];,
                         const FArrayBox& syfab = (*sigma[1])[mfi];,
                         const FArrayBox& szfab = (*sigma[2])[mfi];);
            amrex_mlndlap_adotx(BL_TO_FORTRAN_BOX(bx),
                                BL_TO_FORTRAN_ANYD(yfab),
                                BL_TO_FORTRAN_ANYD(xfab),
                                AMREX_D_DECL(BL_TO_FORTRAN_ANYD(sxfab),
                                             BL_TO_FORTRAN_ANYD(syfab),
                                             BL_TO_FORTRAN_ANYD(szfab)),
                                BL_TO_FORTRAN_ANYD(dg),
                                dxinv);
        }
    }
}

void
MLNodeLaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const
{
    amrex::Abort("MLNodeLaplacian::Fsmooth to be implemented");
}

void
MLNodeLaplacian::FFlux (int amrlev, const MFIter& mfi,
                        const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only) const
{
    amrex::Abort("MLNodeLaplacian::FFlux to be implemented");
}

Real
MLNodeLaplacian::Anorm (int amrlev, int mglev) const
{
    amrex::Abort("MLNodeLaplacian::Anorm to be implemented");
}

}
