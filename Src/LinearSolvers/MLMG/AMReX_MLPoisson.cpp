
#include <AMReX_MLPoisson.H>
#include <AMReX_MLPoisson_F.H>
#include <AMReX_MLALaplacian.H>

namespace amrex {

MLPoisson::MLPoisson (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap,
                      const LPInfo& a_info,
                      const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

void
MLPoisson::define (const Vector<Geometry>& a_geom,
                   const Vector<BoxArray>& a_grids,
                   const Vector<DistributionMapping>& a_dmap,
                   const LPInfo& a_info,
                   const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("MLPoisson::define()");
    MLCellLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

MLPoisson::~MLPoisson ()
{}

void
MLPoisson::prepareForSolve ()
{
    BL_PROFILE("MLPoisson::prepareForSolve()");

    MLCellLinOp::prepareForSolve();

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

#if (AMREX_SPACEDIM != 3)
        const auto& mfac = *m_metric_factor[amrlev][mglev];
        const auto& rc = mfac.cellCenters(mfi);
        const auto& re = mfac.cellEdges(mfi);
        const Box& vbx = mfi.validbox();
#endif
        amrex_mlpoisson_adotx(BL_TO_FORTRAN_BOX(bx),
                              BL_TO_FORTRAN_ANYD(yfab),
                              BL_TO_FORTRAN_ANYD(xfab),
#if (AMREX_SPACEDIM != 3)
                              rc.data(), re.data(), vbx.loVect(), vbx.hiVect(),
#endif
                              dxinv);                                         
    }
}

void
MLPoisson::normalize (int amrlev, int mglev, MultiFab& mf) const
{
#if (AMREX_SPACEDIM != 3)
    BL_PROFILE("MLPoisson::normalize()");

    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox& fab = mf[mfi];

        const auto& mfac = *m_metric_factor[amrlev][mglev];
        const auto& rc = mfac.cellCenters(mfi);
        const auto& re = mfac.cellEdges(mfi);
        const Box& vbx = mfi.validbox();
        amrex_mlpoisson_normalize(BL_TO_FORTRAN_BOX(bx),
                                  BL_TO_FORTRAN_ANYD(fab),
                                  rc.data(), re.data(), vbx.loVect(), vbx.hiVect(),
                                  dxinv);                                     
    }
#endif
}

void
MLPoisson::Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const
{
    BL_PROFILE("MLPoisson::Fsmooth()");

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
        const auto& mfac = *m_metric_factor[amrlev][mglev];
        const auto& rc = mfac.cellCenters(mfi);
        const auto& re = mfac.cellEdges(mfi);
        
        amrex_mlpoisson_gsrb(BL_TO_FORTRAN_BOX(tbx),
                             BL_TO_FORTRAN_ANYD(solnfab),
                             BL_TO_FORTRAN_ANYD(rhsfab),
                             BL_TO_FORTRAN_ANYD(f0fab),
                             BL_TO_FORTRAN_ANYD(f1fab),
                             BL_TO_FORTRAN_ANYD(m0),
                             BL_TO_FORTRAN_ANYD(m1),
                             rc.data(), re.data(),
                             BL_TO_FORTRAN_BOX(vbx), dxinv, redblack);            
#endif

#if (AMREX_SPACEDIM == 2)
        const auto& mfac = *m_metric_factor[amrlev][mglev];
        const auto& rc = mfac.cellCenters(mfi);
        const auto& re = mfac.cellEdges(mfi);
        
        amrex_mlpoisson_gsrb(BL_TO_FORTRAN_BOX(tbx),
                             BL_TO_FORTRAN_ANYD(solnfab),
                             BL_TO_FORTRAN_ANYD(rhsfab),
                             BL_TO_FORTRAN_ANYD(f0fab),
                             BL_TO_FORTRAN_ANYD(f1fab),
                             BL_TO_FORTRAN_ANYD(f2fab),
                             BL_TO_FORTRAN_ANYD(f3fab),
                             BL_TO_FORTRAN_ANYD(m0),
                             BL_TO_FORTRAN_ANYD(m1),
                             BL_TO_FORTRAN_ANYD(m2),
                             BL_TO_FORTRAN_ANYD(m3),
                             rc.data(), re.data(),
                             BL_TO_FORTRAN_BOX(vbx), dxinv, redblack);            
#endif

#if (AMREX_SPACEDIM == 3)
        amrex_mlpoisson_gsrb(BL_TO_FORTRAN_BOX(tbx),
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
                             BL_TO_FORTRAN_BOX(vbx), dxinv, redblack);
#endif
    }
}

void
MLPoisson::FFlux (int amrlev, const MFIter& mfi,
                  const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
                  const FArrayBox& sol, const int face_only) const
{
    BL_PROFILE("MLPoisson::FFlux()");

    const int mglev = 0;
    const Box& box = mfi.tilebox();
    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

#if (AMREX_SPACEDIM != 3)
    const auto& mfac = *m_metric_factor[amrlev][mglev];
    const auto& rc = mfac.cellCenters(mfi);
    const auto& re = mfac.cellEdges(mfi);
    const Box& vbx = m_grids[amrlev][mglev][mfi];
#endif
    amrex_mlpoisson_flux(BL_TO_FORTRAN_BOX(box),
                         AMREX_D_DECL(BL_TO_FORTRAN_ANYD(*flux[0]),
                                      BL_TO_FORTRAN_ANYD(*flux[1]),
                                      BL_TO_FORTRAN_ANYD(*flux[2])),
                         BL_TO_FORTRAN_ANYD(sol),
#if (AMREX_SPACEDIM != 3)
                         rc.data(), re.data(), vbx.loVect(), vbx.hiVect(),
#endif
                         dxinv, face_only);
}

std::unique_ptr<MLLinOp>
MLPoisson::makeNLinOp (int grid_size) const
{
    const Geometry& geom = m_geom[0].back();
    const BoxArray& ba = makeNGrids(grid_size);

    DistributionMapping dm;
    {
        const std::vector<std::vector<int> >& sfc = DistributionMapping::makeSFC(ba);
        Vector<int> pmap(ba.size());
        AMREX_ALWAYS_ASSERT(ParallelContext::CommunicatorSub() == ParallelDescriptor::Communicator());
        const int nprocs = ParallelDescriptor::NProcs();
        for (int iproc = 0; iproc < nprocs; ++iproc) {
            for (int ibox : sfc[iproc]) {
                pmap[ibox] = iproc;
            }
        }
        dm.define(pmap);
    }

    LPInfo minfo{};
    minfo.has_metric_term = info.has_metric_term;

    std::unique_ptr<MLLinOp> r{new MLALaplacian({geom}, {ba}, {dm}, minfo)};

    MLALaplacian* nop = dynamic_cast<MLALaplacian*>(r.get());

    nop->m_parent = this;

    nop->setMaxOrder(maxorder);
    nop->setVerbose(verbose);

    nop->setDomainBC(m_lobc, m_hibc);

    if (needsCoarseDataForBC())
    {
        const Real* dx0 = m_geom[0][0].CellSize();
        const Real fac = 0.5*m_coarse_data_crse_ratio;
        RealVect cbloc {AMREX_D_DECL(dx0[0]*fac, dx0[1]*fac, dx0[2]*fac)};
        nop->setCoarseFineBCLocation(cbloc);
    }

    nop->setScalars(1.0, -1.0);

    const BoxArray& myba = m_grids[0].back();

    const Real* dxinv = geom.InvCellSize();
    Real dxscale = dxinv[0];
#if (AMREX_SPACEDIM >= 2)
    dxscale = std::max(dxscale,dxinv[1]);
#endif
#if (AMREX_SPACEDIM == 3)
    dxscale = std::max(dxscale,dxinv[2]);
#endif

    MultiFab alpha(ba, dm, 1, 0);
    alpha.setVal(1.e30*dxscale*dxscale);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::vector< std::pair<int,Box> > isects;
        
        for (MFIter mfi(alpha, MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = alpha[mfi];
            myba.intersections(fab.box(), isects);
            for (const auto& is : isects)
            {
                fab.setVal(0.0, is.second, 0, 1);
            }
        }
    }

    nop->setACoeffs(0, alpha);

    return r;    
}

}
