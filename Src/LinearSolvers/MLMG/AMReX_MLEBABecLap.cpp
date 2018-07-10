
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFArrayBox.H>

#include <AMReX_MLEBABecLap_F.H>
#include <AMReX_MLABecLap_F.H>
#include <AMReX_ABec_F.H>
#include <AMReX_MG_F.H>

namespace amrex {

MLEBABecLap::MLEBABecLap (const Vector<Geometry>& a_geom,
                          const Vector<BoxArray>& a_grids,
                          const Vector<DistributionMapping>& a_dmap,
                          const LPInfo& a_info,
                          const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

std::unique_ptr<FabFactory<FArrayBox> >
MLEBABecLap::makeFactory (int amrlev, int mglev) const
{
    return makeEBFabFactory(m_geom[amrlev][mglev],
                            m_grids[amrlev][mglev],
                            m_dmap[amrlev][mglev],
                            {1,1,1}, EBSupport::full);
}

void
MLEBABecLap::define (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info,
                     const Vector<EBFArrayBoxFactory const*>& a_factory)
{
    BL_PROFILE("MLEBABecLap::define()");

    Vector<FabFactory<FArrayBox> const*> _factory;
    for (auto x : a_factory) {
        _factory.push_back(static_cast<FabFactory<FArrayBox> const*>(x));
    }

    MLCellLinOp::define(a_geom, a_grids, a_dmap, a_info, _factory);

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
                const int ng = 1;
                m_b_coeffs[amrlev][mglev][idim].define(ba,
                                                       m_dmap[amrlev][mglev],
                                                       1, ng, MFInfo(), *m_factory[amrlev][mglev]);
                m_b_coeffs[amrlev][mglev][idim].setVal(0.0);
            }
        }
    }
}

MLEBABecLap::~MLEBABecLap ()
{}

void
MLEBABecLap::setScalars (Real a, Real b)
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
MLEBABecLap::setACoeffs (int amrlev, const MultiFab& alpha)
{
    MultiFab::Copy(m_a_coeffs[amrlev][0], alpha, 0, 0, 1, 0);
}

void
MLEBABecLap::setBCoeffs (int amrlev, const std::array<MultiFab const*,AMREX_SPACEDIM>& beta)
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        MultiFab::Copy(m_b_coeffs[amrlev][0][idim], *beta[idim], 0, 0, 1, 0);
        m_b_coeffs[amrlev][0][idim].FillBoundary(m_geom[amrlev][0].periodicity());
    }
}

void
MLEBABecLap::averageDownCoeffs ()
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
MLEBABecLap::averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a,
                                            Vector<std::array<MultiFab,AMREX_SPACEDIM> >& b)
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
            amrex::EB_average_down(a[mglev-1], a[mglev], 0, 1, mg_coarsen_ratio);
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
MLEBABecLap::averageDownCoeffsToCoarseAmrLevel (int flev)
{
    amrex::Abort("averageDownCoeffsToCoarseAmrLevel: todo");
}

void
MLEBABecLap::prepareForSolve ()
{
    BL_PROFILE("MLABecLaplacian::prepareForSolve()");

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
MLEBABecLap::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLEBABecLap::Fapply()");

    const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
    AMREX_D_TERM(const MultiFab& bxcoef = m_b_coeffs[amrlev][mglev][0];,
                 const MultiFab& bycoef = m_b_coeffs[amrlev][mglev][1];,
                 const MultiFab& bzcoef = m_b_coeffs[amrlev][mglev][2];);

    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
    auto area = (factory) ? factory->getAreaFrac()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    auto fcent = (factory) ? factory->getFaceCent()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    const MultiCutFab* centroid = (factory) ? &(factory->getCentroid()) : nullptr; 
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(out, MFItInfo().EnableTiling().SetDynamic(true)); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const FArrayBox& xfab = in[mfi];
        FArrayBox& yfab = out[mfi];
        const FArrayBox& afab = acoef[mfi];
        AMREX_D_TERM(const FArrayBox& bxfab = bxcoef[mfi];,
                     const FArrayBox& byfab = bycoef[mfi];,
                     const FArrayBox& bzfab = bzcoef[mfi];);

        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;

        if (fabtyp == FabType::covered) {
            yfab.setVal(0.0, bx, 0, 1);
        } else if (fabtyp == FabType::regular) {
            amrex_mlabeclap_adotx(BL_TO_FORTRAN_BOX(bx),
                                  BL_TO_FORTRAN_ANYD(yfab),
                                  BL_TO_FORTRAN_ANYD(xfab),
                                  BL_TO_FORTRAN_ANYD(afab),
                                  AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bxfab),
                                               BL_TO_FORTRAN_ANYD(byfab),
                                               BL_TO_FORTRAN_ANYD(bzfab)),
                                  dxinv, m_a_scalar, m_b_scalar);
        } else {
            amrex_mlebabeclap_adotx(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_ANYD(yfab),
                                    BL_TO_FORTRAN_ANYD(xfab),
                                    BL_TO_FORTRAN_ANYD(afab),
                                    AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bxfab),
                                                 BL_TO_FORTRAN_ANYD(byfab),
                                                 BL_TO_FORTRAN_ANYD(bzfab)),
                                    BL_TO_FORTRAN_ANYD((*flags)[mfi]),
                                    BL_TO_FORTRAN_ANYD((*vfrac)[mfi]),
                                    AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*area[0])[mfi]),
                                                 BL_TO_FORTRAN_ANYD((*area[1])[mfi]),
                                                 BL_TO_FORTRAN_ANYD((*area[2])[mfi])),
                                    AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*fcent[0])[mfi]),
                                                 BL_TO_FORTRAN_ANYD((*fcent[1])[mfi]),
                                                 BL_TO_FORTRAN_ANYD((*fcent[2])[mfi])),
                                    BL_TO_FORTRAN_ANYD((*centroid)[mfi]),
                                    dxinv, m_a_scalar, m_b_scalar);
        }
    }
}

void
MLEBABecLap::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const
{
    BL_PROFILE("MLEBABecLap::Fsmooth()");

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
    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
    auto area = (factory) ? factory->getAreaFrac()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    auto fcent = (factory) ? factory->getFaceCent()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    const MultiCutFab* centroid = (factory) ? &(factory->getCentroid()) : nullptr; 
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(sol,MFItInfo().SetDynamic(true));
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

        const Box&       vbx     = mfi.validbox();
        const Box&       tbx     = vbx;
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

        auto fabtyp = (flags) ? (*flags)[mfi].getType(tbx) : FabType::regular;

        if (fabtyp == FabType::covered)
        {
            solnfab.setVal(0.0, tbx, 0, nc);
        }
        else if (fabtyp == FabType::regular)
        {
#if (AMREX_SPACEDIM == 1)
            amrex_abec_linesolve (solnfab.dataPtr(), AMREX_ARLIM(solnfab.loVect()),AMREX_ARLIM(solnfab.hiVect()),
                                  rhsfab.dataPtr(), AMREX_ARLIM(rhsfab.loVect()), AMREX_ARLIM(rhsfab.hiVect()),
                                  &m_a_scalar, &m_b_scalar,
                                  afab.dataPtr(), AMREX_ARLIM(afab.loVect()),    AMREX_ARLIM(afab.hiVect()),
                                  bxfab.dataPtr(), AMREX_ARLIM(bxfab.loVect()),   AMREX_ARLIM(bxfab.hiVect()),
                                  f0fab.dataPtr(), AMREX_ARLIM(f0fab.loVect()),   AMREX_ARLIM(f0fab.hiVect()),
                                  m0.dataPtr(), AMREX_ARLIM(m0.loVect()),   AMREX_ARLIM(m0.hiVect()),
                                  f1fab.dataPtr(), AMREX_ARLIM(f1fab.loVect()),   AMREX_ARLIM(f1fab.hiVect()),
                                  m1.dataPtr(), AMREX_ARLIM(m1.loVect()),   AMREX_ARLIM(m1.hiVect()),
                                  tbx.loVect(), tbx.hiVect(), &nc, h);
#endif

#if (AMREX_SPACEDIM == 2)
            amrex_abec_gsrb(solnfab.dataPtr(), AMREX_ARLIM(solnfab.loVect()),AMREX_ARLIM(solnfab.hiVect()),
                            rhsfab.dataPtr(), AMREX_ARLIM(rhsfab.loVect()), AMREX_ARLIM(rhsfab.hiVect()),
                            &m_a_scalar, &m_b_scalar,
                            afab.dataPtr(), AMREX_ARLIM(afab.loVect()),    AMREX_ARLIM(afab.hiVect()),
                            bxfab.dataPtr(), AMREX_ARLIM(bxfab.loVect()),   AMREX_ARLIM(bxfab.hiVect()),
                            byfab.dataPtr(), AMREX_ARLIM(byfab.loVect()),   AMREX_ARLIM(byfab.hiVect()),
                            f0fab.dataPtr(), AMREX_ARLIM(f0fab.loVect()),   AMREX_ARLIM(f0fab.hiVect()),
                            m0.dataPtr(), AMREX_ARLIM(m0.loVect()),   AMREX_ARLIM(m0.hiVect()),
                            f1fab.dataPtr(), AMREX_ARLIM(f1fab.loVect()),   AMREX_ARLIM(f1fab.hiVect()),
                            m1.dataPtr(), AMREX_ARLIM(m1.loVect()),   AMREX_ARLIM(m1.hiVect()),
                            f2fab.dataPtr(), AMREX_ARLIM(f2fab.loVect()),   AMREX_ARLIM(f2fab.hiVect()),
                            m2.dataPtr(), AMREX_ARLIM(m2.loVect()),   AMREX_ARLIM(m2.hiVect()),
                            f3fab.dataPtr(), AMREX_ARLIM(f3fab.loVect()),   AMREX_ARLIM(f3fab.hiVect()),
                            m3.dataPtr(), AMREX_ARLIM(m3.loVect()),   AMREX_ARLIM(m3.hiVect()),
                            tbx.loVect(), tbx.hiVect(), vbx.loVect(), vbx.hiVect(),
                            &nc, h, &redblack);
#endif

#if (AMREX_SPACEDIM == 3)
            amrex_abec_gsrb(solnfab.dataPtr(), AMREX_ARLIM(solnfab.loVect()),AMREX_ARLIM(solnfab.hiVect()),
                            rhsfab.dataPtr(), AMREX_ARLIM(rhsfab.loVect()), AMREX_ARLIM(rhsfab.hiVect()),
                            &m_a_scalar, &m_b_scalar,
                            afab.dataPtr(), AMREX_ARLIM(afab.loVect()), AMREX_ARLIM(afab.hiVect()),
                            bxfab.dataPtr(), AMREX_ARLIM(bxfab.loVect()), AMREX_ARLIM(bxfab.hiVect()),
                            byfab.dataPtr(), AMREX_ARLIM(byfab.loVect()), AMREX_ARLIM(byfab.hiVect()),
                            bzfab.dataPtr(), AMREX_ARLIM(bzfab.loVect()), AMREX_ARLIM(bzfab.hiVect()),
                            f0fab.dataPtr(), AMREX_ARLIM(f0fab.loVect()), AMREX_ARLIM(f0fab.hiVect()),
                            m0.dataPtr(), AMREX_ARLIM(m0.loVect()), AMREX_ARLIM(m0.hiVect()),
                            f1fab.dataPtr(), AMREX_ARLIM(f1fab.loVect()), AMREX_ARLIM(f1fab.hiVect()),
                            m1.dataPtr(), AMREX_ARLIM(m1.loVect()), AMREX_ARLIM(m1.hiVect()),
                            f2fab.dataPtr(), AMREX_ARLIM(f2fab.loVect()), AMREX_ARLIM(f2fab.hiVect()),
                            m2.dataPtr(), AMREX_ARLIM(m2.loVect()), AMREX_ARLIM(m2.hiVect()),
                            f3fab.dataPtr(), AMREX_ARLIM(f3fab.loVect()), AMREX_ARLIM(f3fab.hiVect()),
                            m3.dataPtr(), AMREX_ARLIM(m3.loVect()), AMREX_ARLIM(m3.hiVect()),
                            f4fab.dataPtr(), AMREX_ARLIM(f4fab.loVect()), AMREX_ARLIM(f4fab.hiVect()),
                            m4.dataPtr(), AMREX_ARLIM(m4.loVect()), AMREX_ARLIM(m4.hiVect()),
                            f5fab.dataPtr(), AMREX_ARLIM(f5fab.loVect()), AMREX_ARLIM(f5fab.hiVect()),
                            m5.dataPtr(), AMREX_ARLIM(m5.loVect()), AMREX_ARLIM(m5.hiVect()),
                            tbx.loVect(), tbx.hiVect(), vbx.loVect(), vbx.hiVect(),
                            &nc, h, &redblack);
#endif
        }
        else
        {
            amrex_mlebabeclap_gsrb(BL_TO_FORTRAN_BOX(tbx),
                                   BL_TO_FORTRAN_ANYD(solnfab),
                                   BL_TO_FORTRAN_ANYD(rhsfab),
                                   BL_TO_FORTRAN_ANYD(afab),
                                   AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bxfab),
                                                BL_TO_FORTRAN_ANYD(byfab),
                                                BL_TO_FORTRAN_ANYD(bzfab)),
                                   AMREX_D_DECL(BL_TO_FORTRAN_ANYD(m0),
                                                BL_TO_FORTRAN_ANYD(m2),
                                                BL_TO_FORTRAN_ANYD(m4)),
                                   AMREX_D_DECL(BL_TO_FORTRAN_ANYD(m1),
                                                BL_TO_FORTRAN_ANYD(m3),
                                                BL_TO_FORTRAN_ANYD(m5)),
                                   AMREX_D_DECL(BL_TO_FORTRAN_ANYD(f0fab),
                                                BL_TO_FORTRAN_ANYD(f2fab),
                                                BL_TO_FORTRAN_ANYD(f4fab)),
                                   AMREX_D_DECL(BL_TO_FORTRAN_ANYD(f1fab),
                                                BL_TO_FORTRAN_ANYD(f3fab),
                                                BL_TO_FORTRAN_ANYD(f5fab)),
                                   BL_TO_FORTRAN_ANYD((*flags)[mfi]),
                                   BL_TO_FORTRAN_ANYD((*vfrac)[mfi]),
                                   AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*area[0])[mfi]),
                                                BL_TO_FORTRAN_ANYD((*area[1])[mfi]),
                                                BL_TO_FORTRAN_ANYD((*area[2])[mfi])),
                                   AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*fcent[0])[mfi]),
                                                BL_TO_FORTRAN_ANYD((*fcent[1])[mfi]),
                                                BL_TO_FORTRAN_ANYD((*fcent[2])[mfi])),
                                   BL_TO_FORTRAN_ANYD((*centroid)[mfi]),
                                   dxinv, m_a_scalar, m_b_scalar, redblack);
        }
    }
}

void
MLEBABecLap::FFlux (int amrlev, const MFIter& mfi, const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
                    const FArrayBox& sol, const int face_only) const
{
    amrex::Abort("FFlux: todo");
}

void
MLEBABecLap::normalize (int amrlev, int mglev, MultiFab& mf) const
{
    const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
    AMREX_D_TERM(const MultiFab& bxcoef = m_b_coeffs[amrlev][mglev][0];,
                 const MultiFab& bycoef = m_b_coeffs[amrlev][mglev][1];,
                 const MultiFab& bzcoef = m_b_coeffs[amrlev][mglev][2];);

    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][mglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;
    auto area = (factory) ? factory->getAreaFrac()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    auto fcent = (factory) ? factory->getFaceCent()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    const MultiCutFab* centroid = (factory) ? &(factory->getCentroid()) : nullptr; 

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox& fab = mf[mfi];
        const FArrayBox& afab = acoef[mfi];
        AMREX_D_TERM(const FArrayBox& bxfab = bxcoef[mfi];,
                     const FArrayBox& byfab = bycoef[mfi];,
                     const FArrayBox& bzfab = bzcoef[mfi];);

        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;

        if (fabtyp == FabType::regular)
        {
            amrex_mlabeclap_normalize(BL_TO_FORTRAN_BOX(bx),
                                      BL_TO_FORTRAN_ANYD(fab),
                                      BL_TO_FORTRAN_ANYD(afab),
                                      AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bxfab),
                                                   BL_TO_FORTRAN_ANYD(byfab),
                                                   BL_TO_FORTRAN_ANYD(bzfab)),
                                      dxinv, m_a_scalar, m_b_scalar);
        }
        else if (fabtyp == FabType::singlevalued)
        {
            amrex_mlebabeclap_normalize(BL_TO_FORTRAN_BOX(bx),
                                        BL_TO_FORTRAN_ANYD(fab),
                                        BL_TO_FORTRAN_ANYD(afab),
                                        AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bxfab),
                                                     BL_TO_FORTRAN_ANYD(byfab),
                                                     BL_TO_FORTRAN_ANYD(bzfab)),
                                        BL_TO_FORTRAN_ANYD((*flags)[mfi]),
                                        BL_TO_FORTRAN_ANYD((*vfrac)[mfi]),
                                        AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*area[0])[mfi]),
                                                     BL_TO_FORTRAN_ANYD((*area[1])[mfi]),
                                                     BL_TO_FORTRAN_ANYD((*area[2])[mfi])),
                                        AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*fcent[0])[mfi]),
                                                     BL_TO_FORTRAN_ANYD((*fcent[1])[mfi]),
                                                     BL_TO_FORTRAN_ANYD((*fcent[2])[mfi])),
                                        BL_TO_FORTRAN_ANYD((*centroid)[mfi]), 
                                        dxinv, m_a_scalar, m_b_scalar);
        }
    }
}

void
MLEBABecLap::restriction (int, int, MultiFab& crse, MultiFab& fine) const
{
    const int ncomp = getNComp();
    amrex::EB_average_down(fine, crse, 0, ncomp, 2);
}

void
MLEBABecLap::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][fmglev].get());
    const FabArray<EBCellFlagFab>* flags = (factory) ? &(factory->getMultiEBCellFlagFab()) : nullptr;
    const MultiFab* vfrac = (factory) ? &(factory->getVolFrac()) : nullptr;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(crse,true); mfi.isValid(); ++mfi)
    {
        const Box&         bx    = mfi.tilebox();
        const int          ncomp = getNComp();
        const FArrayBox& cfab    = crse[mfi];
        FArrayBox&       ffab    = fine[mfi];

        auto fabtyp = (flags) ? (*flags)[mfi].getType(amrex::refine(bx,2)) : FabType::regular;

        if (fabtyp == FabType::regular)
        {
            amrex_mg_interp(ffab.dataPtr(),
                            AMREX_ARLIM(ffab.loVect()), AMREX_ARLIM(ffab.hiVect()),
                            cfab.dataPtr(),
                            AMREX_ARLIM(cfab.loVect()), AMREX_ARLIM(cfab.hiVect()),
                            bx.loVect(), bx.hiVect(), &ncomp);
        }
        else if (fabtyp == FabType::singlevalued)
        {
            amrex_eb_mg_interp(BL_TO_FORTRAN_BOX(bx),
                               BL_TO_FORTRAN_ANYD(ffab),
                               BL_TO_FORTRAN_ANYD(cfab),
                               BL_TO_FORTRAN_ANYD((*flags)[mfi]),
                               &ncomp);
        }
    }
}

void
MLEBABecLap::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                     const MultiFab& fine_sol, const MultiFab& fine_rhs)
{
    const auto amrrr = AMRRefRatio(camrlev);
    const int ncomp = getNComp();
    amrex::average_down(fine_sol, crse_sol, 0, ncomp, amrrr);
    amrex::EB_average_down(fine_rhs, crse_rhs, 0, ncomp, amrrr);
}

}
