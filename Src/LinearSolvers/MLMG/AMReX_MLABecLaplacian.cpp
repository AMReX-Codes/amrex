
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MultiFabUtil.H>

#include <AMReX_MLABecLap_K.H>

namespace amrex {

MLABecLaplacian::MLABecLaplacian (const Vector<Geometry>& a_geom,
                                  const Vector<BoxArray>& a_grids,
                                  const Vector<DistributionMapping>& a_dmap,
                                  const LPInfo& a_info,
                                  const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

void
MLABecLaplacian::define (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info,
                         const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("MLABecLaplacian::define()");

    MLCellABecLap::define(a_geom, a_grids, a_dmap, a_info, a_factory);

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
                                             1, 0, MFInfo(), *m_factory[amrlev][mglev]);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& ba = amrex::convert(m_grids[amrlev][mglev],
                                                    IntVect::TheDimensionVector(idim));
                m_b_coeffs[amrlev][mglev][idim].define(ba,
                                                       m_dmap[amrlev][mglev],
                                                       1, 0, MFInfo(), *m_factory[amrlev][mglev]);
            }
        }
    }
}

MLABecLaplacian::~MLABecLaplacian ()
{}

void
MLABecLaplacian::setScalars (Real a, Real b) noexcept
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
MLABecLaplacian::setACoeffs (int amrlev, const MultiFab& alpha)
{
    MultiFab::Copy(m_a_coeffs[amrlev][0], alpha, 0, 0, 1, 0);
    m_needs_update = true;
}

void
MLABecLaplacian::setBCoeffs (int amrlev,
                             const Array<MultiFab const*,AMREX_SPACEDIM>& beta)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        MultiFab::Copy(m_b_coeffs[amrlev][0][idim], *beta[idim], 0, 0, 1, 0);
    }
    m_needs_update = true;
}

void
MLABecLaplacian::averageDownCoeffs ()
{
    BL_PROFILE("MLABecLaplacian::averageDownCoeffs()");

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
                                                Vector<Array<MultiFab,AMREX_SPACEDIM> >& b)
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
        
        Vector<const MultiFab*> fine {AMREX_D_DECL(&(b[mglev-1][0]),
                                                   &(b[mglev-1][1]),
                                                   &(b[mglev-1][2]))};
        Vector<MultiFab*> crse {AMREX_D_DECL(&(b[mglev][0]),
                                             &(b[mglev][1]),
                                             &(b[mglev][2]))};
        IntVect ratio {mg_coarsen_ratio};
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

    if (m_a_scalar != 0.0) {
        amrex::average_down(fine_a_coeffs, crse_a_coeffs, 0, 1, mg_coarsen_ratio);
    }
     
    Array<MultiFab,AMREX_SPACEDIM> bb;
    Vector<MultiFab*> crse(AMREX_SPACEDIM);
    Vector<MultiFab const*> fine(AMREX_SPACEDIM);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        BoxArray ba = fine_b_coeffs[idim].boxArray();
        ba.coarsen(mg_coarsen_ratio);
        bb[idim].define(ba, fine_b_coeffs[idim].DistributionMap(), 1, 0);
        crse[idim] = &bb[idim];
        fine[idim] = &fine_b_coeffs[idim];
    }
    IntVect ratio {mg_coarsen_ratio};
    amrex::average_down_faces(fine, crse, ratio, 0);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        crse_b_coeffs[idim].ParallelCopy(bb[idim], crse_geom.periodicity());
    }
}

void
MLABecLaplacian::applyMetricTermsCoeffs ()
{
#if (AMREX_SPACEDIM != 3)
    for (int alev = 0; alev < m_num_amr_levels; ++alev)
    {
        const int mglev = 0;
        applyMetricTerm(alev, mglev, m_a_coeffs[alev][mglev]);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            applyMetricTerm(alev, mglev, m_b_coeffs[alev][mglev][idim]);
        }
    }
#endif
}

void
MLABecLaplacian::prepareForSolve ()
{
    BL_PROFILE("MLABecLaplacian::prepareForSolve()");

    MLCellABecLap::prepareForSolve();

#if (AMREX_SPACEDIM != 3)
    applyMetricTermsCoeffs();
#endif

    averageDownCoeffs();

    m_is_singular.clear();
    m_is_singular.resize(m_num_amr_levels, false);
    auto itlo = std::find(m_lobc[0].begin(), m_lobc[0].end(), BCType::Dirichlet);
    auto ithi = std::find(m_hibc[0].begin(), m_hibc[0].end(), BCType::Dirichlet);
    if (itlo == m_lobc[0].end() && ithi == m_hibc[0].end())
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

    m_needs_update = false;
}

void
MLABecLaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLABecLaplacian::Fapply()");

    const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
    AMREX_D_TERM(const MultiFab& bxcoef = m_b_coeffs[amrlev][mglev][0];,
                 const MultiFab& bycoef = m_b_coeffs[amrlev][mglev][1];,
                 const MultiFab& bzcoef = m_b_coeffs[amrlev][mglev][2];);

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();

    const Real ascalar = m_a_scalar;
    const Real bscalar = m_b_scalar;

    const int ncomp = getNComp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& xfab = in.array(mfi);
        const auto& yfab = out.array(mfi);
        const auto& afab = acoef.array(mfi);
        AMREX_D_TERM(const auto& bxfab = bxcoef.array(mfi);,
                     const auto& byfab = bycoef.array(mfi);,
                     const auto& bzfab = bzcoef.array(mfi););

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
        {
            mlabeclap_adotx(tbx, yfab, xfab, afab, AMREX_D_DECL(bxfab,byfab,bzfab),
                            dxinv, ascalar, bscalar, ncomp);
        });
    }
}

void
MLABecLaplacian::normalize (int amrlev, int mglev, MultiFab& mf) const
{
    BL_PROFILE("MLABecLaplacian::normalize()");

    const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
    AMREX_D_TERM(const MultiFab& bxcoef = m_b_coeffs[amrlev][mglev][0];,
                 const MultiFab& bycoef = m_b_coeffs[amrlev][mglev][1];,
                 const MultiFab& bzcoef = m_b_coeffs[amrlev][mglev][2];);

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();

    const Real ascalar = m_a_scalar;
    const Real bscalar = m_b_scalar;

    const int ncomp = getNComp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& fab = mf.array(mfi);
        const auto& afab = acoef.array(mfi);
        AMREX_D_TERM(const auto& bxfab = bxcoef.array(mfi);,
                     const auto& byfab = bycoef.array(mfi);,
                     const auto& bzfab = bzcoef.array(mfi););

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
        {
            mlabeclap_normalize(tbx, fab, afab, AMREX_D_DECL(bxfab,byfab,bzfab),
                                dxinv, ascalar, bscalar, ncomp);
        });
    }
}

void
MLABecLaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const
{
    BL_PROFILE("MLABecLaplacian::Fsmooth()");

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

    const int nc = getNComp();
    const Real* h = m_geom[amrlev][mglev].CellSize();
    AMREX_D_TERM(const Real dhx = m_b_scalar/(h[0]*h[0]);,
                 const Real dhy = m_b_scalar/(h[1]*h[1]);,
                 const Real dhz = m_b_scalar/(h[2]*h[2]));
    const Real alpha = m_a_scalar;

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sol,mfi_info); mfi.isValid(); ++mfi)
    {
	const auto& m0 = mm0.array(mfi);
        const auto& m1 = mm1.array(mfi);
#if (AMREX_SPACEDIM > 1)
        const auto& m2 = mm2.array(mfi);
        const auto& m3 = mm3.array(mfi);
#if (AMREX_SPACEDIM > 2)
        const auto& m4 = mm4.array(mfi);
        const auto& m5 = mm5.array(mfi);
#endif
#endif

	const Box& tbx = mfi.tilebox();
        const Box& vbx = mfi.validbox();
        const auto& solnfab = sol.array(mfi);
        const auto& rhsfab  = rhs.array(mfi);
        const auto& afab    = acoef.array(mfi);

        AMREX_D_TERM(const auto& bxfab = bxcoef.array(mfi);,
                     const auto& byfab = bycoef.array(mfi);,
                     const auto& bzfab = bzcoef.array(mfi););

        const auto& f0fab = f0.array(mfi);
        const auto& f1fab = f1.array(mfi);
#if (AMREX_SPACEDIM > 1)
        const auto& f2fab = f2.array(mfi);
        const auto& f3fab = f3.array(mfi);
#if (AMREX_SPACEDIM > 2)
        const auto& f4fab = f4.array(mfi);
        const auto& f5fab = f5.array(mfi);
#endif
#endif

#if (AMREX_SPACEDIM == 1)
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
        {
            abec_gsrb(thread_box, solnfab, rhsfab, alpha, dhx,
                      afab, bxfab,
                      f0fab, m0,
                      f1fab, m1,
                      vbx, nc, redblack);
        });
#endif

#if (AMREX_SPACEDIM == 2)
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
        {
            abec_gsrb(thread_box, solnfab, rhsfab, alpha, dhx, dhy,
                      afab, bxfab, byfab,
                      f0fab, m0,
                      f1fab, m1,
                      f2fab, m2,
                      f3fab, m3,
                      vbx, nc, redblack);
        });
#endif

#if (AMREX_SPACEDIM == 3)
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
        {
            abec_gsrb(thread_box, solnfab, rhsfab, alpha, dhx, dhy, dhz,
                      afab, bxfab, byfab, bzfab,
                      f0fab, m0,
                      f1fab, m1,
                      f2fab, m2,
                      f3fab, m3,
                      f4fab, m4,
                      f5fab, m5,
                      vbx, nc, redblack);
        });
#endif
    }
}

void
MLABecLaplacian::FFlux (int amrlev, const MFIter& mfi,
                        const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, Location, const int face_only) const
{
    BL_PROFILE("MLABecLaplacian::FFlux()");

    const int mglev = 0;
    const Box& box = mfi.tilebox();
    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();
    const int ncomp = getNComp();
    FFlux(box, dxinv, m_b_scalar,
          Array<FArrayBox const*,AMREX_SPACEDIM>{AMREX_D_DECL(&(m_b_coeffs[amrlev][mglev][0][mfi]),
                                                              &(m_b_coeffs[amrlev][mglev][1][mfi]),
                                                              &(m_b_coeffs[amrlev][mglev][2][mfi]))},
          flux, sol, face_only, ncomp);
}

void
MLABecLaplacian::FFlux (Box const& box, Real const* dxinv, Real bscalar,
                        Array<FArrayBox const*, AMREX_SPACEDIM> const& bcoef,
                        Array<FArrayBox*,AMREX_SPACEDIM> const& flux,
                        FArrayBox const& sol, int face_only, int ncomp)
{
    AMREX_D_TERM(const auto bx = bcoef[0]->array();,
                 const auto by = bcoef[1]->array();,
                 const auto bz = bcoef[2]->array(););
    AMREX_D_TERM(const auto& fxarr = flux[0]->array();,
                 const auto& fyarr = flux[1]->array();,
                 const auto& fzarr = flux[2]->array(););
    const auto& solarr = sol.array();

    if (face_only)
    {
        Real fac = bscalar*dxinv[0];
        Box blo = amrex::bdryLo(box, 0);
        int blen = box.length(0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlabeclap_flux_xface(tbox, fxarr, solarr, bx, fac, blen, ncomp);
        });
#if (AMREX_SPACEDIM >= 2)
        fac = bscalar*dxinv[1];
        blo = amrex::bdryLo(box, 1);
        blen = box.length(1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlabeclap_flux_yface(tbox, fyarr, solarr, by, fac, blen, ncomp);
        });
#endif
#if (AMREX_SPACEDIM == 3)
        fac = bscalar*dxinv[2];
        blo = amrex::bdryLo(box, 2);
        blen = box.length(2);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlabeclap_flux_zface(tbox, fzarr, solarr, bz, fac, blen, ncomp);
        });
#endif
    }
    else
    {
        Real fac = bscalar*dxinv[0];
        Box bflux = amrex::surroundingNodes(box, 0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlabeclap_flux_x(tbox, fxarr, solarr, bx, fac, ncomp);
        });
#if (AMREX_SPACEDIM >= 2)
        fac = bscalar*dxinv[1];
        bflux = amrex::surroundingNodes(box, 1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlabeclap_flux_y(tbox, fyarr, solarr, by, fac, ncomp);
        });
#endif
#if (AMREX_SPACEDIM == 3)
        fac = bscalar*dxinv[2];
        bflux = amrex::surroundingNodes(box, 2);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlabeclap_flux_z(tbox, fzarr, solarr, bz, fac, ncomp);
        });
#endif
    }
}

void
MLABecLaplacian::update ()
{
    if (MLCellABecLap::needsUpdate()) MLCellABecLap::update();

#if (AMREX_SPACEDIM != 3)
    applyMetricTermsCoeffs();
#endif

    averageDownCoeffs();

    m_is_singular.clear();
    m_is_singular.resize(m_num_amr_levels, false);
    auto itlo = std::find(m_lobc[0].begin(), m_lobc[0].end(), BCType::Dirichlet);
    auto ithi = std::find(m_hibc[0].begin(), m_hibc[0].end(), BCType::Dirichlet);
    if (itlo == m_lobc[0].end() && ithi == m_hibc[0].end())
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

    m_needs_update = false;
}

}
