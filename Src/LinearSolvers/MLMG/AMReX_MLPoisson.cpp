#include <AMReX_MLPoisson.H>
#include <AMReX_MLPoisson_F.H>

namespace amrex {

MLPoisson::MLPoisson (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap)
{
    define(a_geom, a_grids, a_dmap);
}

void
MLPoisson::define (const Vector<Geometry>& a_geom,
                   const Vector<BoxArray>& a_grids,
                   const Vector<DistributionMapping>& a_dmap)
{
    BL_PROFILE("MLPoisson::define()");
    MLLinOp::define(a_geom, a_grids, a_dmap);
}

MLPoisson::~MLPoisson ()
{}

void
MLPoisson::prepareForSolve ()
{
    BL_PROFILE("MLPoisson::prepareForSolve()");

    m_Anorm.resize(m_num_amr_levels);
    for (int alev = 0; alev < m_num_amr_levels; ++alev)
    {
        m_Anorm[alev].assign(m_num_mg_levels[alev], -1.0);
    }

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
MLPoisson::Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const
{
    BL_PROFILE("MLPoisson::Fapply()");

    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(out, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const FArrayBox& xfab = in[mfi];
        FArrayBox& yfab = out[mfi];

        amrex_mlpoisson_adotx(BL_TO_FORTRAN_BOX(bx),
                              BL_TO_FORTRAN_ANYD(yfab),
                              BL_TO_FORTRAN_ANYD(xfab),
                              dxinv);
    }
}

void
MLPoisson::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const
{
#if 0
    BL_PROFILE("MLPoisson::Fsmooth()");

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
        amrex::Abort("MLPoisson::Fsmooth: 1d not supported");
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
#endif
}

void
MLPoisson::FFlux (int amrlev, const MFIter& mfi,
                        std::array<FArrayBox,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only) const
{
    BL_PROFILE("MLPoisson::FFlux()");

    const int mglev = 0;
    const Box& box = mfi.tilebox();
    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

    amrex_mlpoisson_flux(BL_TO_FORTRAN_BOX(box),
                         AMREX_D_DECL(BL_TO_FORTRAN_ANYD(flux[0]),
                                      BL_TO_FORTRAN_ANYD(flux[1]),
                                      BL_TO_FORTRAN_ANYD(flux[2])),
                         BL_TO_FORTRAN_ANYD(sol),
                         dxinv, face_only);
}

Real
MLPoisson::Anorm (int amrlev, int mglev) const
{
#if 0
    BL_PROFILE("MLPoisson::Anorm()");

    if (m_Anorm[amrlev][mglev] < 0.0)
    {
        const MultiFab& acoef = m_a_coeffs[amrlev][mglev];
        AMREX_D_TERM(const MultiFab& bxcoef = m_b_coeffs[amrlev][mglev][0];,
                     const MultiFab& bycoef = m_b_coeffs[amrlev][mglev][1];,
                     const MultiFab& bzcoef = m_b_coeffs[amrlev][mglev][2];);
        const Real* dx = m_geom[amrlev][mglev].CellSize();

        const int nc = 1;
        Real res = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(max:res)
#endif
        {
            for (MFIter mfi(acoef,true); mfi.isValid(); ++mfi)
            {
                Real tres;
	    
                const Box&       tbx  = mfi.tilebox();
                const FArrayBox& afab = acoef[mfi];
                AMREX_D_TERM(const FArrayBox& bxfab = bxcoef[mfi];,
                             const FArrayBox& byfab = bycoef[mfi];,
                             const FArrayBox& bzfab = bzcoef[mfi];);

#if (BL_SPACEDIM==2)
                FORT_NORMA(&tres,
                           &m_a_scalar, &m_b_scalar,
                           afab.dataPtr(),  ARLIM(afab.loVect()), ARLIM(afab.hiVect()),
                           bxfab.dataPtr(), ARLIM(bxfab.loVect()), ARLIM(bxfab.hiVect()),
                           byfab.dataPtr(), ARLIM(byfab.loVect()), ARLIM(byfab.hiVect()),
                           tbx.loVect(), tbx.hiVect(), &nc, dx);
#elif (BL_SPACEDIM==3)
                
                FORT_NORMA(&tres,
                           &m_a_scalar, &m_b_scalar,
                           afab.dataPtr(),  ARLIM(afab.loVect()), ARLIM(afab.hiVect()),
                           bxfab.dataPtr(), ARLIM(bxfab.loVect()), ARLIM(bxfab.hiVect()),
                           byfab.dataPtr(), ARLIM(byfab.loVect()), ARLIM(byfab.hiVect()),
                           bzfab.dataPtr(), ARLIM(bzfab.loVect()), ARLIM(bzfab.hiVect()),
                           tbx.loVect(), tbx.hiVect(), &nc, dx);
#endif
                
                res = std::max(res, tres);
            }
        }
        
        ParallelDescriptor::ReduceRealMax(res, acoef.color());
        m_Anorm[amrlev][mglev] = res;
    }
#endif

    return m_Anorm[amrlev][mglev];
}

}
