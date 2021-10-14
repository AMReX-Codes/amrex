
#include <AMReX_MLALaplacian.H>
#include <AMReX_MLALap_K.H>
#include <AMReX_MultiFabUtil.H>

namespace amrex {

MLALaplacian::MLALaplacian (const Vector<Geometry>& a_geom,
                            const Vector<BoxArray>& a_grids,
                            const Vector<DistributionMapping>& a_dmap,
                            const LPInfo& a_info,
                            const Vector<FabFactory<FArrayBox> const*>& a_factory,
                            const int a_ncomp)
    : m_ncomp(a_ncomp)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

void
MLALaplacian::define (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap,
                      const LPInfo& a_info,
                      const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("MLALaplacian::define()");

    MLCellABecLap::define(a_geom, a_grids, a_dmap, a_info, a_factory);

    const int ncomp = getNComp();

    m_a_coeffs.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_a_coeffs[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_a_coeffs[amrlev][mglev].define(m_grids[amrlev][mglev],
                                             m_dmap[amrlev][mglev], ncomp, 0);
        }
    }
}

MLALaplacian::~MLALaplacian ()
{}

void
MLALaplacian::setScalars (Real a, Real b) noexcept
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
    const int ncomp = getNComp();
    MultiFab::Copy(m_a_coeffs[amrlev][0], alpha, 0, 0, ncomp, 0);
    m_needs_update = true;
}

void
MLALaplacian::averageDownCoeffs ()
{
    BL_PROFILE("MLALaplacian::averageDownCoeffs()");

    for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
        auto& fine_a_coeffs = m_a_coeffs[amrlev];

        averageDownCoeffsSameAmrLevel(amrlev, fine_a_coeffs);
        averageDownCoeffsToCoarseAmrLevel(amrlev);
    }

    averageDownCoeffsSameAmrLevel(0, m_a_coeffs[0]);
}

void
MLALaplacian::averageDownCoeffsSameAmrLevel (int amrlev, Vector<MultiFab>& a)
{
    const int ncomp = getNComp();
    const int nmglevs = a.size();
    for (int mglev = 1; mglev < nmglevs; ++mglev)
    {
        if (m_a_scalar == 0.0)
        {
            a[mglev].setVal(0.0);
        }
        else
        {
            AMREX_ASSERT(amrlev == 0 || !hasHiddenDimension());
            IntVect ratio = (amrlev > 0) ? IntVect(mg_coarsen_ratio) : mg_coarsen_ratio_vec[mglev-1];
            amrex::average_down(a[mglev-1], a[mglev], 0, ncomp, ratio);
        }
    }
}

void
MLALaplacian::averageDownCoeffsToCoarseAmrLevel (int flev)
{
    const int ncomp = getNComp();
    auto& fine_a_coeffs = m_a_coeffs[flev  ].back();
    auto& crse_a_coeffs = m_a_coeffs[flev-1].front();

    if (m_a_scalar != 0.0) {
        // We coarsen from the back of flev to the front of flev-1.
        // So we use mg_coarsen_ratio.
        amrex::average_down(fine_a_coeffs, crse_a_coeffs, 0, ncomp, mg_coarsen_ratio);
    }
}

void
MLALaplacian::updateSingularFlag ()
{
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
                    // We are only testing component 0 here, assuming the others
                    // are similar.
                    Real asum = m_a_coeffs[alev].back().sum();
                    Real amax = m_a_coeffs[alev].back().norm0();
                    m_is_singular[alev] = (asum <= amax * 1.e-12);
                }
            }
        }
    }
}

void
MLALaplacian::prepareForSolve ()
{
    BL_PROFILE("MLALaplacian::prepareForSolve()");
    MLCellABecLap::prepareForSolve();
    averageDownCoeffs();
    updateSingularFlag();
    m_needs_update = false;
}

void
MLALaplacian::update ()
{
    if (MLCellABecLap::needsUpdate()) MLCellABecLap::update();
    averageDownCoeffs();
    updateSingularFlag();
    m_needs_update = false;
}

void
MLALaplacian::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLALaplacian::Fapply()");

    const int ncomp = getNComp();

    const MultiFab& acoef = m_a_coeffs[amrlev][mglev];

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
#if (AMREX_SPACEDIM < 3)
    const Real dx = m_geom[amrlev][mglev].CellSize(0);
    const Real probxlo = m_geom[amrlev][mglev].ProbLo(0);
#endif

#if (AMREX_SPACEDIM == 3)
    GpuArray<Real,2> dhinv {get_d0(dxinv[0], dxinv[1], dxinv[2]),
                            get_d1(dxinv[0], dxinv[1], dxinv[2])};
#endif

    const Real ascalar = m_a_scalar;
    const Real bscalar = m_b_scalar;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& xfab = in.array(mfi);
        const auto& yfab = out.array(mfi);
        const auto& afab = acoef.array(mfi);

#if (AMREX_SPACEDIM != 3)
        if (m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlalap_adotx_m(tbx, yfab, xfab, afab, dxinv, ascalar, bscalar, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlalap_adotx(tbx, yfab, xfab, afab, dxinv, ascalar, bscalar, ncomp);
            });
        }
#else
        if (hasHiddenDimension()) {
            Box const& bx2d = compactify(bx);
            const auto& xfab2d = compactify(xfab);
            const auto& yfab2d = compactify(yfab);
            const auto& afab2d = compactify(afab);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx2d, tbx2d,
            {
                TwoD::mlalap_adotx(tbx2d, yfab2d, xfab2d, afab2d, dhinv, ascalar, bscalar, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlalap_adotx(tbx, yfab, xfab, afab, dxinv, ascalar, bscalar, ncomp);
            });
        }
#endif
    }
}

void
MLALaplacian::normalize (int amrlev, int mglev, MultiFab& mf) const
{
    BL_PROFILE("MLALaplacian::normalize()");

    const int ncomp = getNComp();

    const MultiFab& acoef = m_a_coeffs[amrlev][mglev];

    const auto dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
#if (AMREX_SPACEDIM < 3)
    const Real dx = m_geom[amrlev][mglev].CellSize(0);
    const Real probxlo = m_geom[amrlev][mglev].ProbLo(0);
#endif

#if (AMREX_SPACEDIM == 3)
    GpuArray<Real,2> dhinv {get_d0(dxinv[0], dxinv[1], dxinv[2]),
                            get_d1(dxinv[0], dxinv[1], dxinv[2])};
#endif

    const Real ascalar = m_a_scalar;
    const Real bscalar = m_b_scalar;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& fab = mf.array(mfi);
        const auto& afab = acoef.array(mfi);

#if (AMREX_SPACEDIM != 3)
        if (m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlalap_normalize_m(tbx, fab, afab, dxinv, ascalar, bscalar, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlalap_normalize(tbx, fab, afab, dxinv, ascalar, bscalar, ncomp);
            });
        }
#else
        if (hasHiddenDimension()) {
            Box const& bx2d = compactify(bx);
            const auto&  fab2d = compactify(fab);
            const auto& afab2d = compactify(afab);
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx2d, tbx2d,
            {
                TwoD::mlalap_normalize(tbx2d, fab2d, afab2d, dhinv, ascalar, bscalar, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
            {
                mlalap_normalize(tbx, fab, afab, dxinv, ascalar, bscalar, ncomp);
            });
        }
#endif
    }
}

void
MLALaplacian::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const
{
    BL_PROFILE("MLALaplacian::Fsmooth()");

    const int ncomp = getNComp();

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
    AMREX_D_TERM(const Real dhx = m_b_scalar*dxinv[0]*dxinv[0];,
                 const Real dhy = m_b_scalar*dxinv[1]*dxinv[1];,
                 const Real dhz = m_b_scalar*dxinv[2]*dxinv[2];);

#if (AMREX_SPACEDIM == 3)
    Real dh0 = get_d0(dhx, dhy, dhz);
    Real dh1 = get_d1(dhx, dhy, dhz);
#endif

#if (AMREX_SPACEDIM < 3)
    const Real dx = m_geom[amrlev][mglev].CellSize(0);
    const Real probxlo = m_geom[amrlev][mglev].ProbLo(0);
#endif

    const Real alpha = m_a_scalar;

    MFItInfo mfi_info;
    if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);

#ifdef AMREX_USE_OMP
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
        if (m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
            {
                mlalap_gsrb_m(thread_box, solnfab, rhsfab, alpha, dhx,
                              afab,
                              f0fab, m0,
                              f1fab, m1,
                              vbx, redblack,
                              dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
            {
                mlalap_gsrb(thread_box, solnfab, rhsfab, alpha, dhx,
                            afab,
                            f0fab, m0,
                            f1fab, m1,
                            vbx, redblack, ncomp);
            });
        }

#endif

#if (AMREX_SPACEDIM == 2)
        if (m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
            {
                mlalap_gsrb_m(thread_box, solnfab, rhsfab, alpha, dhx, dhy,
                              afab,
                              f0fab, m0,
                              f1fab, m1,
                              f2fab, m2,
                              f3fab, m3,
                              vbx, redblack,
                              dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
            {
                mlalap_gsrb(thread_box, solnfab, rhsfab, alpha, dhx, dhy,
                            afab,
                            f0fab, m0,
                            f1fab, m1,
                            f2fab, m2,
                            f3fab, m3,
                            vbx, redblack, ncomp);
            });
        }
#endif

#if (AMREX_SPACEDIM == 3)
        if (hasHiddenDimension()) {
            Box const& tbx_2d = compactify(tbx);
            Box const& vbx_2d = compactify(vbx);
            const auto& solnfab_2d = compactify(solnfab);
            const auto& rhsfab_2d = compactify(rhsfab);
            const auto& afab_2d = compactify(afab);
            const auto& f0fab_2d = compactify(get_d0(f0fab,f1fab,f2fab));
            const auto& f1fab_2d = compactify(get_d1(f0fab,f1fab,f2fab));
            const auto& f2fab_2d = compactify(get_d0(f3fab,f4fab,f5fab));
            const auto& f3fab_2d = compactify(get_d1(f3fab,f4fab,f5fab));
            const auto& m0_2d = compactify(get_d0(m0,m1,m2));
            const auto& m1_2d = compactify(get_d1(m0,m1,m2));
            const auto& m2_2d = compactify(get_d0(m3,m4,m5));
            const auto& m3_2d = compactify(get_d1(m3,m4,m5));
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx_2d, thread_box,
            {
                TwoD::mlalap_gsrb(thread_box, solnfab_2d, rhsfab_2d, alpha, dh0, dh1,
                                  afab_2d,
                                  f0fab_2d, m0_2d,
                                  f1fab_2d, m1_2d,
                                  f2fab_2d, m2_2d,
                                  f3fab_2d, m3_2d,
                                  vbx_2d, redblack, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
            {
                mlalap_gsrb(thread_box, solnfab, rhsfab, alpha, dhx, dhy, dhz,
                            afab,
                            f0fab, m0,
                            f1fab, m1,
                            f2fab, m2,
                            f3fab, m3,
                            f4fab, m4,
                            f5fab, m5,
                            vbx, redblack, ncomp);
            });
        }
#endif
    }
}

void
MLALaplacian::FFlux (int amrlev, const MFIter& mfi,
                     const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                     const FArrayBox& sol, Location, const int face_only) const
{
    AMREX_ASSERT(!hasHiddenDimension());

    BL_PROFILE("MLALaplacian::FFlux()");

    const int ncomp = getNComp();
    const int mglev = 0;
    const Box& box = mfi.tilebox();
    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

    AMREX_D_TERM(const auto& fxarr = flux[0]->array();,
                 const auto& fyarr = flux[1]->array();,
                 const auto& fzarr = flux[2]->array(););
    const auto& solarr = sol.array();

#if (AMREX_SPACEDIM != 3)
    const Real dx = m_geom[amrlev][mglev].CellSize(0);
    const Real probxlo = m_geom[amrlev][mglev].ProbLo(0);
#endif

#if (AMREX_SPACEDIM == 3)
    if (face_only) {
        Real fac = m_b_scalar * dxinv[0];
        Box blo = amrex::bdryLo(box, 0);
        int blen = box.length(0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlalap_flux_xface(tbox, fxarr, solarr, fac, blen, ncomp);
        });
        fac = m_b_scalar * dxinv[1];
        blo = amrex::bdryLo(box, 1);
        blen = box.length(1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlalap_flux_yface(tbox, fyarr, solarr, fac, blen, ncomp);
        });
        fac = m_b_scalar * dxinv[2];
        blo = amrex::bdryLo(box, 2);
        blen = box.length(2);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlalap_flux_zface(tbox, fzarr, solarr, fac, blen, ncomp);
        });
    } else {
        Real fac = m_b_scalar * dxinv[0];
        Box bflux = amrex::surroundingNodes(box, 0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlalap_flux_x(tbox, fxarr, solarr, fac, ncomp);
        });
        fac = m_b_scalar * dxinv[1];
        bflux = amrex::surroundingNodes(box, 1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlalap_flux_y(tbox, fyarr, solarr, fac, ncomp);
        });
        fac = m_b_scalar * dxinv[2];
        bflux = amrex::surroundingNodes(box, 2);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlalap_flux_z(tbox, fzarr, solarr, fac, ncomp);
        });
    }
#elif (AMREX_SPACEDIM == 2)
    if (face_only) {
        Real fac = m_b_scalar * dxinv[0];
        Box blo = amrex::bdryLo(box, 0);
        int blen = box.length(0);
        if (m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlalap_flux_xface_m(tbox, fxarr, solarr, fac, blen, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlalap_flux_xface(tbox, fxarr, solarr, fac, blen, ncomp);
            });
        }
        fac = m_b_scalar * dxinv[1];
        blo = amrex::bdryLo(box, 1);
        blen = box.length(1);
        if (m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlalap_flux_yface_m(tbox, fyarr, solarr, fac, blen, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlalap_flux_yface(tbox, fyarr, solarr, fac, blen, ncomp);
            });
        }
    } else {
        Real fac = m_b_scalar * dxinv[0];
        Box bflux = amrex::surroundingNodes(box, 0);
        if (m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlalap_flux_x_m(tbox, fxarr, solarr, fac, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlalap_flux_x(tbox, fxarr, solarr, fac, ncomp);
            });
        }
        fac = m_b_scalar * dxinv[1];
        bflux = amrex::surroundingNodes(box, 1);
        if (m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlalap_flux_y_m(tbox, fyarr, solarr, fac, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlalap_flux_y(tbox, fyarr, solarr, fac, ncomp);
            });
        }
    }
#else
    if (face_only) {
        Real fac = m_b_scalar * dxinv[0];
        Box blo = amrex::bdryLo(box, 0);
        int blen = box.length(0);
        if (m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlalap_flux_xface_m(tbox, fxarr, solarr, fac, blen, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
            {
                mlalap_flux_xface(tbox, fxarr, solarr, fac, blen, ncomp);
            });
        }
    } else {
        Real fac = m_b_scalar * dxinv[0];
        Box bflux = amrex::surroundingNodes(box, 0);
        if (m_has_metric_term) {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlalap_flux_x_m(tbox, fxarr, solarr, fac, dx, probxlo, ncomp);
            });
        } else {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
            {
                mlalap_flux_x(tbox, fxarr, solarr, fac, ncomp);
            });
        }
    }
#endif
}

}
