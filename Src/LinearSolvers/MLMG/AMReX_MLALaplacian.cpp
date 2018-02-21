
#include <AMReX_MLALaplacian.H>
#include <AMReX_MultiFabUtil.H>

#include <AMReX_MLALap_F.H>

namespace amrex {

MLALaplacian::MLALaplacian (const Vector<Geometry>& a_geom,
                            const Vector<BoxArray>& a_grids,
                            const Vector<DistributionMapping>& a_dmap,
                            const LPInfo& a_info)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMREX_SPACEDIM == 3, "MLALaplacian: only 3d is supported");
    define(a_geom, a_grids, a_dmap, a_info);
}

void
MLALaplacian::define (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap,
                      const LPInfo& a_info)
{
    BL_PROFILE("MLALaplacian::define()");

    MLCellLinOp::define(a_geom, a_grids, a_dmap, a_info);

    m_a_coeffs.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_a_coeffs[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_a_coeffs[amrlev][mglev].define(m_grids[amrlev][mglev],
                                             m_dmap[amrlev][mglev],
                                             1, 0);
        }
    }
}

MLALaplacian::~MLALaplacian ()
{}

void
MLALaplacian::setScalars (Real a, Real b)
{
    m_a_scalar = a;
    m_b_scalar = b;
    if (a == 0.0)
    {
        for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
        {
            m_a_coeffs[amrlev][0].setVal(0.0);
        }
    }
}

void
MLALaplacian::setACoeffs (int amrlev, const MultiFab& alpha)
{
    MultiFab::Copy(m_a_coeffs[amrlev][0], alpha, 0, 0, 1, 0);
}

void
MLALaplacian::averageDownCoeffs ()
{
    BL_PROFILE("MLALaplacian::averageDownCoeffs()");

    for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
        auto& fine_a_coeffs = m_a_coeffs[amrlev];

        averageDownCoeffsSameAmrLevel(fine_a_coeffs);
        averageDownCoeffsToCoarseAmrLevel(amrlev);
    }

    averageDownCoeffsSameAmrLevel(m_a_coeffs[0]);
}

void
MLALaplacian::averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a)
{
    int nmglevs = a.size();
    for (int mglev = 1; mglev < nmglevs; ++mglev)
    {
        if (m_a_scalar == 0.0)
        {
            a[mglev].setVal(0.0);
        }
        else
        {
            amrex::average_down(a[mglev-1], a[mglev], 0, 1, mg_coarsen_ratio);
        }
    }
}

void
MLALaplacian::averageDownCoeffsToCoarseAmrLevel (int flev)
{
    auto& fine_a_coeffs = m_a_coeffs[flev  ].back();
    auto& crse_a_coeffs = m_a_coeffs[flev-1].front();

    if (m_a_scalar != 0.0) {
        amrex::average_down(fine_a_coeffs, crse_a_coeffs, 0, 1, mg_coarsen_ratio);
    }     
}

void
MLALaplacian::prepareForSolve ()
{
    BL_PROFILE("MLALaplacian::prepareForSolve()");

    MLCellLinOp::prepareForSolve();

    averageDownCoeffs();

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
                if (m_a_scalar == 0.0)
                {
                    m_is_singular[alev] = true;
                }
                else
                {
                    Real asum = m_a_coeffs[alev].back().sum();
                    Real amax = m_a_coeffs[alev].back().norm0();
                    m_is_singular[alev] = (asum <= amax * 1.e-12);
                }
            }
        }
    }
}

void
MLALaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLALaplacian::Fapply()");

    const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(out, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const FArrayBox& xfab = in[mfi];
        FArrayBox& yfab = out[mfi];
        const FArrayBox& afab = acoef[mfi];

#if (AMREX_SPACEDIM != 3)
        const auto& mfac = *m_metric_factor[amrlev][mglev];
        const auto& rc = mfac.cellCenters(mfi);
        const auto& re = mfac.cellEdges(mfi);
        const Box& vbx = mfi.validbox();
#endif

        amrex_mlalap_adotx(BL_TO_FORTRAN_BOX(bx),
                           BL_TO_FORTRAN_ANYD(yfab),
                           BL_TO_FORTRAN_ANYD(xfab),
                           BL_TO_FORTRAN_ANYD(afab),
#if (AMREX_SPACEDIM != 3)
                           rc.data(), re.data(), vbx.loVect(), vbx.hiVect(),
#endif
                           dxinv, m_a_scalar, m_b_scalar);

    }
}

void
MLALaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const
{
    BL_PROFILE("MLALaplacian::Fsmooth()");

    const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
    const auto& undrrelxr = m_undrrelxr[amrlev][mglev];
    const auto& maskvals  = m_maskvals [amrlev][mglev];

    OrientationIter oitr;

    const FabSet& f0 = undrrelxr[oitr()]; ++oitr;
    const FabSet& f1 = undrrelxr[oitr()]; ++oitr;
#if (AMREX_SPACEDIM > 1)
    const FabSet& f2 = undrrelxr[oitr()]; ++oitr;
    const FabSet& f3 = undrrelxr[oitr()]; ++oitr;
#if (AMREX_SPACEDIM > 2)
    const FabSet& f4 = undrrelxr[oitr()]; ++oitr;
    const FabSet& f5 = undrrelxr[oitr()]; ++oitr;
#endif
#endif

    const MultiMask& mm0 = maskvals[0];
    const MultiMask& mm1 = maskvals[1];
#if (AMREX_SPACEDIM > 1)
    const MultiMask& mm2 = maskvals[2];
    const MultiMask& mm3 = maskvals[3];
#if (AMREX_SPACEDIM > 2)
    const MultiMask& mm4 = maskvals[4];
    const MultiMask& mm5 = maskvals[5];
#endif
#endif

    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(sol,MFItInfo().EnableTiling().SetDynamic(true));
         mfi.isValid(); ++mfi)
    {
	const Mask& m0 = mm0[mfi];
        const Mask& m1 = mm1[mfi];
#if (AMREX_SPACEDIM > 1)
        const Mask& m2 = mm2[mfi];
        const Mask& m3 = mm3[mfi];
#if (AMREX_SPACEDIM > 2)
        const Mask& m4 = mm4[mfi];
        const Mask& m5 = mm5[mfi];
#endif
#endif

	const Box&       tbx     = mfi.tilebox();
        const Box&       vbx     = mfi.validbox();
        FArrayBox&       solnfab = sol[mfi];
        const FArrayBox& rhsfab  = rhs[mfi];
        const FArrayBox& afab    = acoef[mfi];

        const FArrayBox& f0fab = f0[mfi];
        const FArrayBox& f1fab = f1[mfi];
#if (AMREX_SPACEDIM > 1)
        const FArrayBox& f2fab = f2[mfi];
        const FArrayBox& f3fab = f3[mfi];
#if (AMREX_SPACEDIM > 2)
        const FArrayBox& f4fab = f4[mfi];
        const FArrayBox& f5fab = f5[mfi];
#endif
#endif

#if (AMREX_SPACEDIM == 1)
#endif

#if (AMREX_SPACEDIM == 2)
#endif

#if (AMREX_SPACEDIM == 3)
        amrex_mlalap_gsrb(BL_TO_FORTRAN_BOX(tbx),
                          BL_TO_FORTRAN_ANYD(solnfab),
                          BL_TO_FORTRAN_ANYD(rhsfab),
                          BL_TO_FORTRAN_ANYD(f0fab),
                          BL_TO_FORTRAN_ANYD(f1fab),
                          BL_TO_FORTRAN_ANYD(f2fab),
                          BL_TO_FORTRAN_ANYD(f3fab),
                          BL_TO_FORTRAN_ANYD(f4fab),
                          BL_TO_FORTRAN_ANYD(f5fab),
                          BL_TO_FORTRAN_ANYD(m0),
                          BL_TO_FORTRAN_ANYD(m1),
                          BL_TO_FORTRAN_ANYD(m2),
                          BL_TO_FORTRAN_ANYD(m3),
                          BL_TO_FORTRAN_ANYD(m4),
                          BL_TO_FORTRAN_ANYD(m5),
                          BL_TO_FORTRAN_ANYD(afab),
                          BL_TO_FORTRAN_BOX(vbx), dxinv,
                          &m_a_scalar, &m_b_scalar,
                          redblack);
#endif
    }
}

void
MLALaplacian::FFlux (int amrlev, const MFIter& mfi,
                        const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only) const
{
    BL_PROFILE("MLALaplacian::FFlux()");

    const int mglev = 0;
    const Box& box = mfi.tilebox();
    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

    amrex_mlalap_flux(BL_TO_FORTRAN_BOX(box),
                      AMREX_D_DECL(BL_TO_FORTRAN_ANYD(*flux[0]),
                                   BL_TO_FORTRAN_ANYD(*flux[1]),
                                   BL_TO_FORTRAN_ANYD(*flux[2])),
                      BL_TO_FORTRAN_ANYD(sol),
                      dxinv, m_b_scalar, face_only);
}

}
