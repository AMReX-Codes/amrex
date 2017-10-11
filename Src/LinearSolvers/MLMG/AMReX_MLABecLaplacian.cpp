
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MultiFabUtil.H>

#include <AMReX_MLABecLap_F.H>
#include <AMReX_ABec_F.H>

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
MLABecLaplacian::setBCoeffs (int amrlev,
                             const std::array<MultiFab const*,AMREX_SPACEDIM>& beta)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        MultiFab::Copy(m_b_coeffs[amrlev][0][idim], *beta[idim], 0, 0, 1, 0);
    }
}

void
MLABecLaplacian::averageDownCoeffs ()
{
    for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
        auto& fine_a_coeffs = m_a_coeffs[amrlev];
        auto& fine_b_coeffs = m_b_coeffs[amrlev];

        averageDownCoeffsSameAmrLevel(fine_a_coeffs, fine_b_coeffs);
        averageDownCoeffsToCoarseAmrLevel(amrlev);
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

void
MLABecLaplacian::averageDownCoeffsToCoarseAmrLevel (int flev)
{
    auto& fine_a_coeffs = m_a_coeffs[flev  ].back();
    auto& fine_b_coeffs = m_b_coeffs[flev  ].back();
    auto& crse_a_coeffs = m_a_coeffs[flev-1].front();
    auto& crse_b_coeffs = m_b_coeffs[flev-1].front();
    auto& crse_geom     = m_geom    [flev-1][0];

    amrex::average_down(fine_a_coeffs, crse_a_coeffs, 0, 1, mg_coarsen_ratio);
     
    std::array<MultiFab,AMREX_SPACEDIM> bb;
    Vector<MultiFab*> crse(AMREX_SPACEDIM);
    Vector<MultiFab const*> fine(AMREX_SPACEDIM);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        BoxArray ba = fine_b_coeffs[idim].boxArray();
        ba.coarsen(mg_coarsen_ratio);
        bb[idim].define(ba, fine_b_coeffs[idim].DistributionMap(), 1, 0);
        crse[idim] = &bb[idim];
        fine[idim] = &fine_b_coeffs[idim];
    }
    IntVect ratio {AMREX_D_DECL(mg_coarsen_ratio, mg_coarsen_ratio, mg_coarsen_ratio)};
    amrex::average_down_faces(fine, crse, ratio, 0);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        crse_b_coeffs[idim].ParallelCopy(bb[idim], crse_geom.periodicity());
    }
}

void
MLABecLaplacian::prepareForSolve ()
{
    averageDownCoeffs();
}

void
MLABecLaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
    AMREX_D_TERM(const MultiFab& bxcoef = m_b_coeffs[amrlev][mglev][0];,
                 const MultiFab& bycoef = m_b_coeffs[amrlev][mglev][1];,
                 const MultiFab& bzcoef = m_b_coeffs[amrlev][mglev][2];);

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
        AMREX_D_TERM(const FArrayBox& bxfab = bxcoef[mfi];,
                     const FArrayBox& byfab = bycoef[mfi];,
                     const FArrayBox& bzfab = bzcoef[mfi];);

        amrex_mlabeclap_adotx(BL_TO_FORTRAN_BOX(bx),
                              BL_TO_FORTRAN_ANYD(yfab),
                              BL_TO_FORTRAN_ANYD(xfab),
                              BL_TO_FORTRAN_ANYD(afab),
                              AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bxfab),
                                           BL_TO_FORTRAN_ANYD(byfab),
                                           BL_TO_FORTRAN_ANYD(bzfab)),
                              dxinv, m_a_scalar, m_b_scalar);

    }
}

void
MLABecLaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const
{
    const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
    AMREX_D_TERM(const MultiFab& bxcoef = m_b_coeffs[amrlev][mglev][0];,
                 const MultiFab& bycoef = m_b_coeffs[amrlev][mglev][1];,
                 const MultiFab& bzcoef = m_b_coeffs[amrlev][mglev][2];);
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

    const int nc = 1;
    const Real* h = m_geom[amrlev][mglev].CellSize();

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

        AMREX_D_TERM(const FArrayBox& bxfab = bxcoef[mfi];,
                     const FArrayBox& byfab = bycoef[mfi];,
                     const FArrayBox& bzfab = bzcoef[mfi];);

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
        amrex::Abort("MLABecLaplacian::Fsmooth: 1d not supported");
#endif

#if (AMREX_SPACEDIM == 2)
        FORT_GSRB(solnfab.dataPtr(), ARLIM(solnfab.loVect()),ARLIM(solnfab.hiVect()),
                  rhsfab.dataPtr(), ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
                  &m_a_scalar, &m_b_scalar,
                  afab.dataPtr(), ARLIM(afab.loVect()),    ARLIM(afab.hiVect()),
                  bxfab.dataPtr(), ARLIM(bxfab.loVect()),   ARLIM(bxfab.hiVect()),
                  byfab.dataPtr(), ARLIM(byfab.loVect()),   ARLIM(byfab.hiVect()),
                  f0fab.dataPtr(), ARLIM(f0fab.loVect()),   ARLIM(f0fab.hiVect()),
                  m0.dataPtr(), ARLIM(m0.loVect()),   ARLIM(m0.hiVect()),
                  f1fab.dataPtr(), ARLIM(f1fab.loVect()),   ARLIM(f1fab.hiVect()),
                  m1.dataPtr(), ARLIM(m1.loVect()),   ARLIM(m1.hiVect()),
                  f2fab.dataPtr(), ARLIM(f2fab.loVect()),   ARLIM(f2fab.hiVect()),
                  m2.dataPtr(), ARLIM(m2.loVect()),   ARLIM(m2.hiVect()),
                  f3fab.dataPtr(), ARLIM(f3fab.loVect()),   ARLIM(f3fab.hiVect()),
                  m3.dataPtr(), ARLIM(m3.loVect()),   ARLIM(m3.hiVect()),
                  tbx.loVect(), tbx.hiVect(), vbx.loVect(), vbx.hiVect(),
                  &nc, h, &redblack);
#endif

#if (AMREX_SPACEDIM == 3)
        FORT_GSRB(solnfab.dataPtr(), ARLIM(solnfab.loVect()),ARLIM(solnfab.hiVect()),
                  rhsfab.dataPtr(), ARLIM(rhsfab.loVect()), ARLIM(rhsfab.hiVect()),
                  &m_a_scalar, &m_b_scalar,
                  afab.dataPtr(), ARLIM(afab.loVect()), ARLIM(afab.hiVect()),
                  bxfab.dataPtr(), ARLIM(bxfab.loVect()), ARLIM(bxfab.hiVect()),
                  byfab.dataPtr(), ARLIM(byfab.loVect()), ARLIM(byfab.hiVect()),
                  bzfab.dataPtr(), ARLIM(bzfab.loVect()), ARLIM(bzfab.hiVect()),
                  f0fab.dataPtr(), ARLIM(f0fab.loVect()), ARLIM(f0fab.hiVect()),
                  m0.dataPtr(), ARLIM(m0.loVect()), ARLIM(m0.hiVect()),
                  f1fab.dataPtr(), ARLIM(f1fab.loVect()), ARLIM(f1fab.hiVect()),
                  m1.dataPtr(), ARLIM(m1.loVect()), ARLIM(m1.hiVect()),
                  f2fab.dataPtr(), ARLIM(f2fab.loVect()), ARLIM(f2fab.hiVect()),
                  m2.dataPtr(), ARLIM(m2.loVect()), ARLIM(m2.hiVect()),
                  f3fab.dataPtr(), ARLIM(f3fab.loVect()), ARLIM(f3fab.hiVect()),
                  m3.dataPtr(), ARLIM(m3.loVect()), ARLIM(m3.hiVect()),
                  f4fab.dataPtr(), ARLIM(f4fab.loVect()), ARLIM(f4fab.hiVect()),
                  m4.dataPtr(), ARLIM(m4.loVect()), ARLIM(m4.hiVect()),
                  f5fab.dataPtr(), ARLIM(f5fab.loVect()), ARLIM(f5fab.hiVect()),
                  m5.dataPtr(), ARLIM(m5.loVect()), ARLIM(m5.hiVect()),
                  tbx.loVect(), tbx.hiVect(), vbx.loVect(), vbx.hiVect(),
                  &nc, h, &redblack);
#endif
    }
}

}
