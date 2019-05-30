
#include <AMReX_MLPoisson.H>
#include <AMReX_MLPoisson_K.H>
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
    MLCellABecLap::define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

MLPoisson::~MLPoisson ()
{}

void
MLPoisson::prepareForSolve ()
{
    BL_PROFILE("MLPoisson::prepareForSolve()");

    MLCellABecLap::prepareForSolve();

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

    AMREX_D_TERM(const Real dhx = dxinv[0]*dxinv[0];,
                 const Real dhy = dxinv[1]*dxinv[1];,
                 const Real dhz = dxinv[2]*dxinv[2];);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& xfab = in.array(mfi);
        const auto& yfab = out.array(mfi);

#if (AMREX_SPACEDIM == 3)
        AMREX_HOST_DEVICE_FOR_3D (bx, i, j, k,
        {
            mlpoisson_adotx(i, j, k, yfab, xfab, dhx, dhy, dhz);
        });
#else
        const auto& mfac = *m_metric_factor[amrlev][mglev];
        const auto& re = mfac.cellEdges(mfi);
        AsyncArray<Real> aa_re(re.data(), re.size());
        Real const* rep = aa_re.data();
        const int rlo = mfi.validbox().smallEnd(0);

#if (AMREX_SPACEDIM == 2)
        const auto& rc = mfac.cellCenters(mfi);
        AsyncArray<Real> aa_rc(rc.data(), rc.size());
        Real const* rcp = aa_rc.data();
        AMREX_HOST_DEVICE_FOR_3D (bx, i, j, k,
        {
            mlpoisson_adotx(i, j, k, yfab, xfab, dhx, dhy, rcp, rep, rlo);
        });
#elif (AMREX_SPACEDIM == 1)
        AMREX_HOST_DEVICE_FOR_3D (bx, i, j, k,
        {
            mlpoisson_adotx(i, j, k, yfab, xfab, dhx, rep, rlo);
        });
#endif
#endif
    }
}

void
MLPoisson::normalize (int amrlev, int mglev, MultiFab& mf) const
{
#if (AMREX_SPACEDIM != 3)
    BL_PROFILE("MLPoisson::normalize()");

    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();
    AMREX_D_TERM(const Real dhx = dxinv[0]*dxinv[0];,
                 const Real dhy = dxinv[1]*dxinv[1];,
                 const Real dhz = dxinv[2]*dxinv[2];);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const auto& fab = mf.array(mfi);

        const auto& mfac = *m_metric_factor[amrlev][mglev];
        const auto& re = mfac.cellEdges(mfi);
        AsyncArray<Real> aa_re(re.data(), re.size());
        Real const* rep = aa_re.data();
        const int rlo = mfi.validbox().smallEnd(0);

#if (AMREX_SPACEDIM == 2)
        const auto& rc = mfac.cellCenters(mfi);
        AsyncArray<Real> aa_rc(rc.data(), rc.size());
        Real const* rcp = aa_rc.data();
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, thread_box,
        {
            mlpoisson_normalize(thread_box, fab, rcp, rep, rlo, dhx, dhy);
        });
#else
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, thread_box,
        {
            mlpoisson_normalize(thread_box, fab, rep, rlo, dhx);
        });
#endif
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
    AMREX_D_TERM(const Real dhx = dxinv[0]*dxinv[0];,
                 const Real dhy = dxinv[1]*dxinv[1];,
                 const Real dhz = dxinv[2]*dxinv[2];);

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

#if (AMREX_SPACEDIM != 3)
        const auto& mfac = *m_metric_factor[amrlev][mglev];
        const auto& re = mfac.cellEdges(mfi);
        AsyncArray<Real> aa_re(re.data(), re.size());
        Real const* rep = aa_re.data();
        const int rlo = mfi.validbox().smallEnd(0);
#endif

#if (AMREX_SPACEDIM == 1)
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
        {
            mlpoisson_gsrb(thread_box, solnfab, rhsfab, dhx,
                           f0fab, m0,
                           f1fab, m1,
                           vbx, redblack,
                           rep, rlo);
        });
#endif

#if (AMREX_SPACEDIM == 2)
        const auto& rc = mfac.cellCenters(mfi);
        AsyncArray<Real> aa_rc(rc.data(), rc.size());
        Real const* rcp = aa_rc.data();
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
        {
            mlpoisson_gsrb(thread_box, solnfab, rhsfab, dhx, dhy,
                           f0fab, m0,
                           f1fab, m1,
                           f2fab, m2,
                           f3fab, m3,
                           vbx, redblack,
                           rcp, rep, rlo);
        });
#endif

#if (AMREX_SPACEDIM == 3)
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( tbx, thread_box,
        {
            mlpoisson_gsrb(thread_box, solnfab, rhsfab, dhx, dhy, dhz,
                           f0fab, m0,
                           f1fab, m1,
                           f2fab, m2,
                           f3fab, m3,
                           f4fab, m4,
                           f5fab, m5,
                           vbx, redblack);
        });
#endif
    }
}

void
MLPoisson::FFlux (int amrlev, const MFIter& mfi,
                  const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                  const FArrayBox& sol, Location, const int face_only) const
{
    BL_PROFILE("MLPoisson::FFlux()");

    const int mglev = 0;
    const Box& box = mfi.tilebox();
    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

    AMREX_D_TERM(const auto& fxarr = flux[0]->array();,
                 const auto& fyarr = flux[1]->array();,
                 const auto& fzarr = flux[2]->array(););
    const auto& solarr = sol.array();

#if (AMREX_SPACEDIM != 3)
    const auto& mfac = *m_metric_factor[amrlev][mglev];
    const auto& re = mfac.cellEdges(mfi);
    AsyncArray<Real> aa_re(re.data(), re.size());
    Real const* rep = aa_re.data();
    const int rlo = m_grids[amrlev][mglev][mfi].smallEnd(0);
#endif
#if (AMREX_SPACEDIM == 2)
    const auto& rc = mfac.cellCenters(mfi);
    AsyncArray<Real> aa_rc(rc.data(), rc.size());
    Real const* rcp = aa_rc.data();
#endif

#if (AMREX_SPACEDIM == 3)
    if (face_only) {
        Real fac = dxinv[0];
        Box blo = amrex::bdryLo(box, 0);
        int blen = box.length(0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlpoisson_flux_xface(tbox, fxarr, solarr, fac, blen);
        });
        fac = dxinv[1];
        blo = amrex::bdryLo(box, 1);
        blen = box.length(1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlpoisson_flux_yface(tbox, fyarr, solarr, fac, blen);
        });
        fac = dxinv[2];
        blo = amrex::bdryLo(box, 2);
        blen = box.length(2);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlpoisson_flux_zface(tbox, fzarr, solarr, fac, blen);
        });
    } else {
        Real fac = dxinv[0];
        Box bflux = amrex::surroundingNodes(box, 0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlpoisson_flux_x(tbox, fxarr, solarr, fac);
        });
        fac = dxinv[1];
        bflux = amrex::surroundingNodes(box, 1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlpoisson_flux_y(tbox, fyarr, solarr, fac);
        });
        fac = dxinv[2];
        bflux = amrex::surroundingNodes(box, 2);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlpoisson_flux_z(tbox, fzarr, solarr, fac);
        });
    }
#elif (AMREX_SPACEDIM == 2)
    if (face_only) {
        Real fac = dxinv[0];
        Box blo = amrex::bdryLo(box, 0);
        int blen = box.length(0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlpoisson_flux_xface(tbox, fxarr, solarr, fac, blen, rep, rlo);
        });
        fac = dxinv[1];
        blo = amrex::bdryLo(box, 1);
        blen = box.length(1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlpoisson_flux_yface(tbox, fyarr, solarr, fac, blen, rcp, rlo);
        });
    } else {
        Real fac = dxinv[0];
        Box bflux = amrex::surroundingNodes(box, 0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlpoisson_flux_x(tbox, fxarr, solarr, fac, rep, rlo);
        });
        fac = dxinv[1];
        bflux = amrex::surroundingNodes(box, 1);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlpoisson_flux_y(tbox, fyarr, solarr, fac, rcp, rlo);
        });
    }
#else
    if (face_only) {
        Real fac = dxinv[0];
        Box blo = amrex::bdryLo(box, 0);
        int blen = box.length(0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( blo, tbox,
        {
            mlpoisson_flux_xface(tbox, fxarr, solarr, fac, blen, rep, rlo);
        });
    } else {
        Real fac = dxinv[0];
        Box bflux = amrex::surroundingNodes(box, 0);
        AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bflux, tbox,
        {
            mlpoisson_flux_x(tbox, fxarr, solarr, fac, rep, rlo);
        });
    }
#endif
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
        dm.define(std::move(pmap));
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
