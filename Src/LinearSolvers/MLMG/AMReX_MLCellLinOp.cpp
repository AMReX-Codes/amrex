
#include <AMReX_MLCellLinOp.H>
#include <AMReX_MLLinOp_F.H>

namespace amrex {

MLCellLinOp::MLCellLinOp () {}

MLCellLinOp::~MLCellLinOp () {}

void
MLCellLinOp::define (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info)
{
    MLLinOp::define(a_geom, a_grids, a_dmap, a_info);
    defineAuxData();
    defineBC();
}

void
MLCellLinOp::defineAuxData ()
{
    m_undrrelxr.resize(m_num_amr_levels);
    m_maskvals.resize(m_num_amr_levels);
    m_fluxreg.resize(m_num_amr_levels-1);

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_undrrelxr[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_undrrelxr[amrlev][mglev].define(m_grids[amrlev][mglev],
                                              m_dmap[amrlev][mglev],
                                              1, 0, 0, 1);
        }
    }
    
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_maskvals[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            for (OrientationIter oitr; oitr; ++oitr)
            {
                const Orientation face = oitr();
                const int ngrow = 1;
                m_maskvals[amrlev][mglev][face].define(m_grids[amrlev][mglev],
                                                       m_dmap[amrlev][mglev],
                                                       m_geom[amrlev][mglev],
                                                       face, 0, ngrow, 0, 1, true);
            }
        }
    }

    for (int amrlev = 0; amrlev < m_num_amr_levels-1; ++amrlev)
    {
        const IntVect ratio{m_amr_ref_ratio[amrlev]};
        m_fluxreg[amrlev].define(m_grids[amrlev+1][0], m_grids[amrlev][0],
                                 m_dmap[amrlev+1][0], m_dmap[amrlev][0],
                                 m_geom[amrlev+1][0], m_geom[amrlev][0],
                                 ratio, amrlev+1, 1);
    }

#if (AMREX_SPACEDIM != 3)
    bool no_metric_term = Geometry::IsCartesian() || !info.has_metric_term;
    m_metric_factor.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_metric_factor[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_metric_factor[amrlev][mglev].reset(new MetricFactor(m_grids[amrlev][mglev],
                                                                  m_dmap[amrlev][mglev],
                                                                  m_geom[amrlev][mglev],
                                                                  no_metric_term));
        }
    }
#endif
}

void
MLCellLinOp::defineBC ()
{
    m_bndry_sol.resize(m_num_amr_levels);
    m_crse_sol_br.resize(m_num_amr_levels);

    m_bndry_cor.resize(m_num_amr_levels);
    m_crse_cor_br.resize(m_num_amr_levels);

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_bndry_sol[amrlev].reset(new MLMGBndry(m_grids[amrlev][0], m_dmap[amrlev][0],
                                                1, m_geom[amrlev][0]));
    }

    for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
    {
        const int ncomp = 1;
        const int in_rad = 0;
        const int out_rad = 1;
        const int extent_rad = 2;
        const int crse_ratio = m_amr_ref_ratio[amrlev-1];
        BoxArray cba = m_grids[amrlev][0];
        cba.coarsen(crse_ratio);
        m_crse_sol_br[amrlev].reset(new BndryRegister(cba, m_dmap[amrlev][0],
                                                      in_rad, out_rad, extent_rad, ncomp));
    }

    for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
    {
        const int ncomp = 1;
        const int in_rad = 0;
        const int out_rad = 1;
        const int extent_rad = 2;
        const int crse_ratio = m_amr_ref_ratio[amrlev-1];
        BoxArray cba = m_grids[amrlev][0];
        cba.coarsen(crse_ratio);
        m_crse_cor_br[amrlev].reset(new BndryRegister(cba, m_dmap[amrlev][0],
                                                      in_rad, out_rad, extent_rad, ncomp));
        m_crse_cor_br[amrlev]->setVal(0.0);
    }

    // This has be to done after m_crse_cor_br is defined.
    for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_bndry_cor[amrlev].reset(new MLMGBndry(m_grids[amrlev][0], m_dmap[amrlev][0],
                                                1, m_geom[amrlev][0]));
        MultiFab bc_data(m_grids[amrlev][0], m_dmap[amrlev][0], 1, 1);
        bc_data.setVal(0.0);
        m_bndry_cor[amrlev]->setBndryValues(*m_crse_cor_br[amrlev], 0, bc_data, 0, 0, 1,
                                            m_amr_ref_ratio[amrlev-1], BCRec());
        m_bndry_cor[amrlev]->setLOBndryConds({AMREX_D_DECL(BCType::Dirichlet,
                                                           BCType::Dirichlet,
                                                           BCType::Dirichlet)},
                                             {AMREX_D_DECL(BCType::Dirichlet,
                                                           BCType::Dirichlet,
                                                           BCType::Dirichlet)},
                                              m_amr_ref_ratio[amrlev-1]);
    }

    m_bcondloc.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_bcondloc[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_bcondloc[amrlev][mglev].reset(new BndryCondLoc(m_grids[amrlev][mglev],
                                                             m_dmap[amrlev][mglev]));
        } 
    }
}

void
MLCellLinOp::setLevelBC (int amrlev, const MultiFab* a_levelbcdata)
{
    BL_PROFILE("MLCellLinOp::setLevelBC()");

    AMREX_ALWAYS_ASSERT(amrlev >= 0 && amrlev < m_num_amr_levels);

    MultiFab zero;
    if (a_levelbcdata == nullptr) {
        zero.define(m_grids[amrlev][0], m_dmap[amrlev][0], 1, 1);
        zero.setVal(0.0);
    } else {
        AMREX_ALWAYS_ASSERT(a_levelbcdata->nGrow() >= 1);
    }
    const MultiFab& bcdata = (a_levelbcdata == nullptr) ? zero : *a_levelbcdata;

    int br_ref_ratio;

    if (amrlev == 0)
    {
        if (needsCoarseDataForBC())
        {
            if (m_crse_sol_br[amrlev] == nullptr)
            {
                const int ncomp = 1;
                const int in_rad = 0;
                const int out_rad = 1;
                const int extent_rad = 2;
                const int crse_ratio = m_coarse_data_crse_ratio;
                BoxArray cba = m_grids[amrlev][0];
                cba.coarsen(crse_ratio);
                m_crse_sol_br[amrlev].reset(new BndryRegister(cba, m_dmap[amrlev][0],
                                                              in_rad, out_rad,
                                                              extent_rad, ncomp));
            }
            if (m_coarse_data_for_bc != nullptr) {
                const Box& cbx = amrex::coarsen(m_geom[0][0].Domain(), m_coarse_data_crse_ratio);
                m_crse_sol_br[amrlev]->copyFrom(*m_coarse_data_for_bc, 0, 0, 0, 1,
                                                Geometry::periodicity(cbx));
            } else {
                m_crse_sol_br[amrlev]->setVal(0.0);
            }
            m_bndry_sol[amrlev]->setBndryValues(*m_crse_sol_br[amrlev], 0,
                                                bcdata, 0, 0, 1,
                                                m_coarse_data_crse_ratio, BCRec());
            br_ref_ratio = m_coarse_data_crse_ratio;
        }
        else
        {
            m_bndry_sol[amrlev]->setBndryValues(bcdata,0,0,1,BCRec());
            br_ref_ratio = 1;
        }
    }
    else
    {
        m_bndry_sol[amrlev]->setBndryValues(bcdata,0,0,1, m_amr_ref_ratio[amrlev-1], BCRec());
        br_ref_ratio = m_amr_ref_ratio[amrlev-1];
    }

    m_bndry_sol[amrlev]->setLOBndryConds(m_lobc, m_hibc, br_ref_ratio);

    const Real* dx = m_geom[amrlev][0].CellSize();
    for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
    {
        m_bcondloc[amrlev][mglev]->setLOBndryConds(m_geom[amrlev][mglev], dx,
                                                   m_lobc, m_hibc, br_ref_ratio);
    }
}

void
MLCellLinOp::apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                    const MLMGBndry* bndry) const
{
    BL_PROFILE("MLCellLinOp::apply()");
    applyBC(amrlev, mglev, in, bc_mode, bndry);
    Fapply(amrlev, mglev, out, in);
}

void
MLCellLinOp::smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                     bool skip_fillboundary) const
{
    BL_PROFILE("MLCellLinOp::smooth()");
    for (int redblack = 0; redblack < 2; ++redblack)
    {
        applyBC(amrlev, mglev, sol, BCMode::Homogeneous, nullptr, skip_fillboundary);
        Fsmooth(amrlev, mglev, sol, rhs, redblack);
        skip_fillboundary = false;
    }
}

void
MLCellLinOp::updateSolBC (int amrlev, const MultiFab& crse_bcdata) const
{
    BL_PROFILE("MLCellLinOp::updateSolBC()");

    AMREX_ALWAYS_ASSERT(amrlev > 0);
    m_crse_sol_br[amrlev]->copyFrom(crse_bcdata, 0, 0, 0, 1, m_geom[amrlev-1][0].periodicity());
    m_bndry_sol[amrlev]->updateBndryValues(*m_crse_sol_br[amrlev], 0, 0, 1, m_amr_ref_ratio[amrlev-1]);
}

void
MLCellLinOp::updateCorBC (int amrlev, const MultiFab& crse_bcdata) const
{
    BL_PROFILE("MLCellLinOp::updateCorBC()");
    AMREX_ALWAYS_ASSERT(amrlev > 0);
    m_crse_cor_br[amrlev]->copyFrom(crse_bcdata, 0, 0, 0, 1, m_geom[amrlev-1][0].periodicity());
    m_bndry_cor[amrlev]->updateBndryValues(*m_crse_cor_br[amrlev], 0, 0, 1, m_amr_ref_ratio[amrlev-1]);
}

void
MLCellLinOp::solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                           const MultiFab* crse_bcdata)
{
    BL_PROFILE("MLCellLinOp::solutionResidual()");
    if (crse_bcdata != nullptr) {
        updateSolBC(amrlev, *crse_bcdata);
    }
    const int mglev = 0;
    apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, m_bndry_sol[amrlev].get());
    MultiFab::Xpay(resid, -1.0, b, 0, 0, resid.nComp(), 0);
}

void
MLCellLinOp::fillSolutionBC (int amrlev, MultiFab& sol, const MultiFab* crse_bcdata)
{
    BL_PROFILE("MLCellLinOp::fillSolutionBC()");
    if (crse_bcdata != nullptr) {
        updateSolBC(amrlev, *crse_bcdata);
    }
    const int mglev = 0;
    applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, m_bndry_sol[amrlev].get());    
}

void
MLCellLinOp::correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                             BCMode bc_mode, const MultiFab* crse_bcdata)
{
    BL_PROFILE("MLCellLinOp::correctionResidual()");
    if (bc_mode == BCMode::Inhomogeneous)
    {
        if (crse_bcdata)
        {
            AMREX_ALWAYS_ASSERT(mglev == 0);
            AMREX_ALWAYS_ASSERT(amrlev > 0);
            updateCorBC(amrlev, *crse_bcdata);
        }
        apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, m_bndry_cor[amrlev].get());
    }
    else
    {
        AMREX_ALWAYS_ASSERT(crse_bcdata == nullptr);
        apply(amrlev, mglev, resid, x, BCMode::Homogeneous, nullptr);
    }

    MultiFab::Xpay(resid, -1.0, b, 0, 0, resid.nComp(), 0);
}

void
MLCellLinOp::applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode,
                      const MLMGBndry* bndry, bool skip_fillboundary) const
{
    BL_PROFILE("MLCellLinOp::applyBC()");
    // No coarsened boundary values, cannot apply inhomog at mglev>0.
    BL_ASSERT(mglev == 0 || bc_mode == BCMode::Homogeneous);
    BL_ASSERT(bndry != nullptr || bc_mode == BCMode::Homogeneous);

    const bool cross = true;
    if (!skip_fillboundary) {
        in.FillBoundary(0, 1, m_geom[amrlev][mglev].periodicity(), cross);
    }

    int flagbc = (bc_mode == BCMode::Homogeneous) ? 0 : 1;

    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

    const auto& maskvals = m_maskvals[amrlev][mglev];

    const auto& bcondloc = *m_bcondloc[amrlev][mglev];

    FArrayBox foo(Box::TheUnitBox());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(in, MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
    {
        const Box& vbx   = mfi.validbox();
        FArrayBox& iofab = in[mfi];

        const RealTuple & bdl = bcondloc.bndryLocs(mfi);
        const BCTuple   & bdc = bcondloc.bndryConds(mfi);

        for (OrientationIter oitr; oitr; ++oitr)
        {
            const Orientation ori = oitr();

            int  cdr = ori;
            Real bcl = bdl[ori];
            int  bct = bdc[ori];

            const FArrayBox& fsfab = (bndry != nullptr) ? bndry->bndryValues(ori)[mfi] : foo;

            const Mask& m = maskvals[ori][mfi];

            amrex_mllinop_apply_bc(BL_TO_FORTRAN_BOX(vbx),
                                   BL_TO_FORTRAN_ANYD(iofab),
                                   BL_TO_FORTRAN_ANYD(m),
                                   cdr, bct, bcl,
                                   BL_TO_FORTRAN_ANYD(fsfab),
                                   maxorder, dxinv, flagbc);
        }
    }
}

void
MLCellLinOp::reflux (int crse_amrlev, MultiFab& res,
                     const MultiFab& crse_sol, MultiFab& fine_sol) const
{
    BL_PROFILE("MLCellLinOp::reflux()");
    YAFluxRegister& fluxreg = m_fluxreg[crse_amrlev];
    fluxreg.reset();

    const int fine_amrlev = crse_amrlev+1;

    Real dt = 1.0;
    const Real* crse_dx = m_geom[crse_amrlev][0].CellSize();
    const Real* fine_dx = m_geom[fine_amrlev][0].CellSize();

    const int mglev = 0;
    applyBC(fine_amrlev, mglev, fine_sol, BCMode::Inhomogeneous, m_bndry_sol[fine_amrlev].get());

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::array<FArrayBox,AMREX_SPACEDIM> flux;
        std::array<FArrayBox*,AMREX_SPACEDIM> pflux { AMREX_D_DECL(&flux[0], &flux[1], &flux[2]) };
        std::array<FArrayBox const*,AMREX_SPACEDIM> cpflux { AMREX_D_DECL(&flux[0], &flux[1], &flux[2]) };

        for (MFIter mfi(crse_sol, MFItInfo().EnableTiling().SetDynamic(true));  mfi.isValid(); ++mfi)
        {
            if (fluxreg.CrseHasWork(mfi))
            {
                const Box& tbx = mfi.tilebox();
                AMREX_D_TERM(flux[0].resize(amrex::surroundingNodes(tbx,0));,
                             flux[1].resize(amrex::surroundingNodes(tbx,1));,
                             flux[2].resize(amrex::surroundingNodes(tbx,2)););
                FFlux(crse_amrlev, mfi, pflux, crse_sol[mfi]);
                fluxreg.CrseAdd(mfi, cpflux, crse_dx, dt);
            }
        }

#ifdef _OPENMP
#pragma omp barrier
#endif

        for (MFIter mfi(fine_sol, MFItInfo().EnableTiling().SetDynamic(true));  mfi.isValid(); ++mfi)
        {
            if (fluxreg.FineHasWork(mfi))
            {
                const Box& tbx = mfi.tilebox();
                AMREX_D_TERM(flux[0].resize(amrex::surroundingNodes(tbx,0));,
                             flux[1].resize(amrex::surroundingNodes(tbx,1));,
                             flux[2].resize(amrex::surroundingNodes(tbx,2)););
                const int face_only = true;
                FFlux(fine_amrlev, mfi, pflux, fine_sol[mfi], face_only);
                fluxreg.FineAdd(mfi, cpflux, fine_dx, dt);            
            }
        }
    }

    fluxreg.Reflux(res);
}

void
MLCellLinOp::compFlux (int amrlev, const std::array<MultiFab*,AMREX_SPACEDIM>& fluxes, MultiFab& sol) const
{
    BL_PROFILE("MLCellLinOp::compFlux()");

    const int mglev = 0;
    applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, m_bndry_sol[amrlev].get());

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::array<FArrayBox,AMREX_SPACEDIM> flux;
        std::array<FArrayBox*,AMREX_SPACEDIM> pflux { AMREX_D_DECL(&flux[0], &flux[1], &flux[2]) };
        for (MFIter mfi(sol, MFItInfo().EnableTiling().SetDynamic(true));  mfi.isValid(); ++mfi)
        {
            const Box& tbx = mfi.tilebox();
            AMREX_D_TERM(flux[0].resize(amrex::surroundingNodes(tbx,0));,
                         flux[1].resize(amrex::surroundingNodes(tbx,1));,
                         flux[2].resize(amrex::surroundingNodes(tbx,2)););
            FFlux(amrlev, mfi, pflux, sol[mfi]);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const Box& nbx = mfi.nodaltilebox(idim);
                (*fluxes[idim])[mfi].copy(flux[idim], nbx, 0, nbx, 0, 1);
            }
        }
    }
}

void
MLCellLinOp::compGrad (int amrlev, const std::array<MultiFab*,AMREX_SPACEDIM>& grad, MultiFab& sol) const
{
    BL_PROFILE("MLCellLinOp::compGrad()");

    const int mglev = 0;
    applyBC(amrlev, mglev, sol, BCMode::Inhomogeneous, m_bndry_sol[amrlev].get());

    const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(sol, MFItInfo().EnableTiling().SetDynamic(true));  mfi.isValid(); ++mfi)
    {
        AMREX_D_TERM(const Box& xbx = mfi.nodaltilebox(0);,
                     const Box& ybx = mfi.nodaltilebox(1);,
                     const Box& zbx = mfi.nodaltilebox(2););
        amrex_mllinop_grad(AMREX_D_DECL(BL_TO_FORTRAN_BOX(xbx),
                                        BL_TO_FORTRAN_BOX(ybx),
                                        BL_TO_FORTRAN_BOX(zbx)),
                           BL_TO_FORTRAN_ANYD(sol[mfi]),
                           AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*grad[0])[mfi]),
                                        BL_TO_FORTRAN_ANYD((*grad[1])[mfi]),
                                        BL_TO_FORTRAN_ANYD((*grad[2])[mfi])),
                           dxinv);
    }
}

void
MLCellLinOp::prepareForSolve ()
{
    BL_PROFILE("MLCellLinOp::prepareForSolve()");

    for (int amrlev = 0;  amrlev < m_num_amr_levels; ++amrlev)
    {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            const auto& bcondloc = *m_bcondloc[amrlev][mglev];
            const auto& maskvals = m_maskvals[amrlev][mglev];
            const Real* dxinv = m_geom[amrlev][mglev].InvCellSize();

            BndryRegister& undrrelxr = m_undrrelxr[amrlev][mglev];
            MultiFab foo(m_grids[amrlev][mglev], m_dmap[amrlev][mglev], 1, 0, MFInfo().SetAlloc(false));
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(foo, MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
            {
                const Box& vbx = mfi.validbox();

                const RealTuple & bdl = bcondloc.bndryLocs(mfi);
                const BCTuple   & bdc = bcondloc.bndryConds(mfi);

                for (OrientationIter oitr; oitr; ++oitr)
                {
                    const Orientation ori = oitr();
                    
                    int  cdr = ori;
                    Real bcl = bdl[ori];
                    int  bct = bdc[ori];
                    
                    FArrayBox& ffab = undrrelxr[ori][mfi];
                    const Mask& m   =  maskvals[ori][mfi];

                    amrex_mllinop_comp_interp_coef0(BL_TO_FORTRAN_BOX(vbx),
                                                    BL_TO_FORTRAN_ANYD(ffab),
                                                    BL_TO_FORTRAN_ANYD(m),
                                                    cdr, bct, bcl, maxorder, dxinv);
                }
            }
        }
    }
}

MLCellLinOp::BndryCondLoc::BndryCondLoc (const BoxArray& ba, const DistributionMapping& dm)
    : bcond(ba, dm),
      bcloc(ba, dm)
{
}

void
MLCellLinOp::BndryCondLoc::setLOBndryConds (const Geometry& geom, const Real* dx,
                                            const std::array<BCType,AMREX_SPACEDIM>& lobc,
                                            const std::array<BCType,AMREX_SPACEDIM>& hibc,
                                            int ratio)
{
    const Box&  domain = geom.Domain();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(bcloc); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        RealTuple & bloc  = bcloc[mfi];
        BCTuple   & bctag = bcond[mfi];

        MLMGBndry::setBoxBC(bloc, bctag, bx, domain, lobc, hibc, dx, ratio);
    }
}

void
MLCellLinOp::applyMetricTerm (int amrlev, int mglev, MultiFab& rhs) const
{
#if (AMREX_SPACEDIM != 3)
    
    if (Geometry::IsCartesian() || !info.has_metric_term) return;

    const auto& mfac = *m_metric_factor[amrlev][mglev];

    int nextra = rhs.ixType().cellCentered(0) ? 0 : 1;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rhs,true); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        const Vector<Real>& r = (nextra==0) ? mfac.cellCenters(mfi) : mfac.cellEdges(mfi);
        const Box& vbx = mfi.validbox();
        const int rlo = vbx.loVect()[0];
        const int rhi = vbx.hiVect()[0] + nextra;
        amrex_mllinop_apply_metric(BL_TO_FORTRAN_BOX(tbx),
                                   BL_TO_FORTRAN_ANYD(rhs[mfi]),
                                   r.data(), &rlo, &rhi);
    }
#endif
}

void
MLCellLinOp::unapplyMetricTerm (int amrlev, int mglev, MultiFab& rhs) const
{
#if (AMREX_SPACEDIM != 3)

    if (Geometry::IsCartesian() || !info.has_metric_term) return;

    const auto& mfac = *m_metric_factor[amrlev][mglev];

    int nextra = rhs.ixType().cellCentered(0) ? 0 : 1;
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rhs,true); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        const Vector<Real>& r = (nextra==0) ? mfac.invCellCenters(mfi) : mfac.invCellEdges(mfi);
        const Box& vbx = mfi.validbox();
        const int rlo = vbx.loVect()[0];
        const int rhi = vbx.hiVect()[0] + nextra;
        amrex_mllinop_apply_metric(BL_TO_FORTRAN_BOX(tbx),
                                   BL_TO_FORTRAN_ANYD(rhs[mfi]),
                                   r.data(), &rlo, &rhi);
    }
#endif
}

#if (AMREX_SPACEDIM != 3)

MLCellLinOp::MetricFactor::MetricFactor (const BoxArray& ba, const DistributionMapping& dm,
                                         const Geometry& geom, bool null_metric)
    : r_cellcenter(ba,dm),
      r_celledge(ba,dm),
      inv_r_cellcenter(ba,dm),
      inv_r_celledge(ba,dm)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(r_cellcenter); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        if (null_metric)
        {
            const int N = bx.length(0);
            auto& rcc = r_cellcenter[mfi];
            rcc.resize(N, 1.0);
            auto& rce = r_celledge[mfi];
            rce.resize(N+1, 1.0);
            auto& ircc = inv_r_cellcenter[mfi];
            ircc.resize(N, 1.0);
            auto& irce = inv_r_celledge[mfi];
            irce.resize(N+1, 1.0);
        }
        else
        {
            auto& rcc = r_cellcenter[mfi];
            geom.GetCellLoc(rcc, bx, 0);
            
            auto& rce = r_celledge[mfi];
            geom.GetEdgeLoc(rce, bx, 0);
            
            auto& ircc = inv_r_cellcenter[mfi];
            const int N = rcc.size();
            ircc.resize(N);
            for (int i = 0; i < N; ++i) {
                ircc[i] = 1.0/rcc[i];
            }
            
            auto& irce = inv_r_celledge[mfi];
            irce.resize(N+1);
            if (rce[0] == 0.0) {
                irce[0] = 0.0;
            } else {
                irce[0] = 1.0/rce[0];
            }
            for (int i = 1; i < N+1; ++i) {
                irce[i] = 1.0/rce[i];
            }

            if (Geometry::IsSPHERICAL()) {
                for (auto& x : rcc) {
                    x = x*x;
                }
                for (auto& x : rce) {
                    x = x*x;
                }
                for (auto& x : ircc) {
                    x = x*x;
                }
                for (auto& x : irce) {
                    x = x*x;
                }                    
            }
        }
    }
}

#endif

}
