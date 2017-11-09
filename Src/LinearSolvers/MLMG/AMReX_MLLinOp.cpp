
#include <AMReX_VisMF.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_MLLinOp_F.H>
#include <AMReX_LO_BCTYPES.H>

namespace amrex {

constexpr int MLLinOp::mg_coarsen_ratio;
constexpr int MLLinOp::mg_box_min_width;
int MLLinOp::do_agglomeration = 1;
int MLLinOp::do_consolidation = 1;
#if (AMREX_SPACEDIM == 1)
int MLLinOp::agg_grid_size = 32;
#elif (AMREX_SPACEDIM == 2)
int MLLinOp::agg_grid_size = 16;
#elif (AMREX_SPACEDIM == 3)
int MLLinOp::agg_grid_size = 8;
#endif

MLLinOp::MLLinOp (const Vector<Geometry>& a_geom,
                  const Vector<BoxArray>& a_grids,
                  const Vector<DistributionMapping>& a_dmap)
{
    define(a_geom, a_grids, a_dmap);
}

void
MLLinOp::define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap)
{
    BL_PROFILE("MLLinOp::define()");
    defineGrids(a_geom, a_grids, a_dmap);
    defineAuxData();
    defineBC();
}

void
MLLinOp::defineGrids (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap)
{
    m_num_amr_levels = a_geom.size();

    m_amr_ref_ratio.resize(m_num_amr_levels);
    m_num_mg_levels.resize(m_num_amr_levels);

    m_geom.resize(m_num_amr_levels);
    m_grids.resize(m_num_amr_levels);
    m_dmap.resize(m_num_amr_levels);

    m_default_comm = ParallelDescriptor::Communicator();

    // fine amr levels
    for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
        m_num_mg_levels[amrlev] = 1;
        m_geom[amrlev].push_back(a_geom[amrlev]);
        m_grids[amrlev].push_back(a_grids[amrlev]);
        m_dmap[amrlev].push_back(a_dmap[amrlev]);

        int rr = mg_coarsen_ratio;
        const Box& dom = a_geom[amrlev].Domain();
        for (int i = 0; i < 30; ++i)
        {
            if (!dom.coarsenable(rr)) amrex::Abort("MLLinOp: Uncoarsenable domain");

            const Box& cdom = amrex::coarsen(dom,rr);
            if (cdom == a_geom[amrlev-1].Domain()) break;

            ++(m_num_mg_levels[amrlev]);

            m_geom[amrlev].emplace_back(cdom);

            m_grids[amrlev].push_back(a_grids[amrlev]);
            AMREX_ASSERT(m_grids[amrlev].back().coarsenable(rr));
            m_grids[amrlev].back().coarsen(rr);

            m_dmap[amrlev].push_back(a_dmap[amrlev]);

            rr *= mg_coarsen_ratio;
        }

        m_amr_ref_ratio[amrlev-1] = rr;
    }

    // coarsest amr level
    m_num_mg_levels[0] = 1;
    m_geom[0].push_back(a_geom[0]);
    m_grids[0].push_back(a_grids[0]);
    m_dmap[0].push_back(a_dmap[0]);

    m_domain_covered.resize(m_num_amr_levels, false);
    m_domain_covered[0] = (m_grids[0][0].numPts() == m_geom[0][0].Domain().numPts());
    for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
    {
        if (!m_domain_covered[amrlev-1]) break;
        m_domain_covered[amrlev] = (m_grids[amrlev][0].numPts() == m_geom[amrlev][0].Domain().numPts());
    }

    if (do_agglomeration && m_domain_covered[0])
    {
        Vector<Box> domainboxes;
        Box dbx = m_geom[0][0].Domain();
        Real nbxs = static_cast<Real>(m_grids[0][0].size());
        Real threshold_npts = static_cast<Real>(AMREX_D_TERM(agg_grid_size,*agg_grid_size,*agg_grid_size));
        Vector<int> agg_flag;
        domainboxes.push_back(dbx);
        agg_flag.push_back(false);
        while (dbx.coarsenable(mg_coarsen_ratio,mg_box_min_width))
        {
            dbx.coarsen(mg_coarsen_ratio);
            domainboxes.push_back(dbx);
            bool to_agg = (dbx.d_numPts() / nbxs) < 0.999*threshold_npts;
            agg_flag.push_back(to_agg);
        }

        int first_agglev = std::distance(agg_flag.begin(),
                                         std::find(agg_flag.begin(),agg_flag.end(),1));
        int nmaxlev = domainboxes.size();
        int rr = mg_coarsen_ratio;
        for (int lev = 1; lev < nmaxlev; ++lev)
        {
            if (lev >= first_agglev or !a_grids[0].coarsenable(rr,mg_box_min_width))
            {
                m_geom[0].emplace_back(domainboxes[lev]);
            
                m_grids[0].emplace_back(domainboxes[lev]);
                m_grids[0].back().maxSize(agg_grid_size);

                const auto& dm = makeConsolidatedDMap(m_grids[0].back());
                m_dmap[0].push_back(dm);
            }
            else
            {
                m_geom[0].emplace_back(amrex::coarsen(a_geom[0].Domain(),rr));
                
                m_grids[0].push_back(a_grids[0]);
                m_grids[0].back().coarsen(rr);
            
                m_dmap[0].push_back(a_dmap[0]);
            }

            ++(m_num_mg_levels[0]);
            rr *= mg_coarsen_ratio;
        }
    }
    else
    {
        int rr = mg_coarsen_ratio;
        Real avg_npts, threshold_npts;
        if (do_consolidation) {
            avg_npts = static_cast<Real>(a_grids[0].d_numPts()) / static_cast<Real>(ParallelDescriptor::NProcs());
            threshold_npts = static_cast<Real>(AMREX_D_TERM(agg_grid_size,*agg_grid_size,*agg_grid_size));
        }
        while (a_geom[0].Domain().coarsenable(rr)
               and a_grids[0].coarsenable(rr, mg_box_min_width))
        {
            m_geom[0].emplace_back(amrex::coarsen(a_geom[0].Domain(),rr));
            
            m_grids[0].push_back(a_grids[0]);
            m_grids[0].back().coarsen(rr);

            if (do_consolidation)
            {
                if (avg_npts/(AMREX_D_TERM(rr,*rr,*rr)) < 0.999*threshold_npts)
                {
                    const auto& dm = makeConsolidatedDMap(m_dmap[0].back());
                    m_dmap[0].push_back(dm);
                }
                else
                {
                    auto dm = m_dmap[0].back();
                    m_dmap[0].push_back(dm);
                }
            }
            else
            {
                m_dmap[0].push_back(a_dmap[0]);
            }
            
            ++(m_num_mg_levels[0]);
            rr *= mg_coarsen_ratio;
        }
    }

    if (do_agglomeration || do_consolidation)
    {
        m_bottom_comm = makeSubCommunicator(m_dmap[0].back());
    }
    else
    {
        m_bottom_comm = m_default_comm;
    }
}

void
MLLinOp::defineAuxData ()
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
    m_metric_factor.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_metric_factor[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_metric_factor[amrlev][mglev].reset(new MetricFactor(m_grids[amrlev][mglev],
                                                                  m_dmap[amrlev][mglev],
                                                                  m_geom[amrlev][mglev]));
        }
    }
#endif
}

void
MLLinOp::defineBC ()
{
    m_bndry_sol.resize(m_num_amr_levels);
    m_crse_sol_br.resize(m_num_amr_levels);

    m_bndry_cor.resize(m_num_amr_levels);
    m_crse_cor_br.resize(m_num_amr_levels);

    m_needs_coarse_data_for_bc = !m_domain_covered[0];
        
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_bndry_sol[amrlev].reset(new MacBndry(m_grids[amrlev][0], m_dmap[amrlev][0],
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
        m_bndry_cor[amrlev].reset(new MacBndry(m_grids[amrlev][0], m_dmap[amrlev][0],
                                              1, m_geom[amrlev][0]));
        // this will make it Dirichlet, what we want for correction
        BCRec phys_bc{AMREX_D_DECL(Outflow,Outflow,Outflow),
                      AMREX_D_DECL(Outflow,Outflow,Outflow)};
        
        MultiFab bc_data(m_grids[amrlev][0], m_dmap[amrlev][0], 1, 1);
        bc_data.setVal(0.0);
        m_bndry_cor[amrlev]->setBndryValues(*m_crse_cor_br[amrlev], 0, bc_data, 0, 0, 1,
                                            m_amr_ref_ratio[amrlev-1], phys_bc);
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

MLLinOp::~MLLinOp ()
{}

void
MLLinOp::make (Vector<Vector<MultiFab> >& mf, int nc, int ng) const
{
    mf.clear();
    mf.resize(m_num_amr_levels);
    for (int alev = 0; alev < m_num_amr_levels; ++alev)
    {
        mf[alev].resize(m_num_mg_levels[alev]);
        for (int mlev = 0; mlev < m_num_mg_levels[alev]; ++mlev)
        {
            mf[alev][mlev].define(m_grids[alev][mlev], m_dmap[alev][mlev], nc, ng);
        }
    }
}

void
MLLinOp::setDomainBC (const std::array<BCType,AMREX_SPACEDIM>& a_lobc,
                      const std::array<BCType,AMREX_SPACEDIM>& a_hibc)
{
    m_lobc = a_lobc;
    m_hibc = a_hibc;
}

void
MLLinOp::setBCWithCoarseData (const MultiFab* crse, int crse_ratio)
{
    m_coarse_data_for_bc = crse;
    m_coarse_data_crse_ratio = crse_ratio;
}

void
MLLinOp::setLevelBC (int amrlev, const MultiFab* a_levelbcdata)
{
    BL_PROFILE("MLLinOp::setLevelBC()");

    BCRec phys_bc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        AMREX_ALWAYS_ASSERT(m_lobc[idim] != BCType::Bogus);
        if (m_lobc[idim] == BCType::Dirichlet) {
            phys_bc.setLo(idim, Outflow);  // MacBndry treats it as Dirichlet
        } else {
            phys_bc.setLo(idim, Inflow);   // MacBndry treats it as Neumann
        }
        AMREX_ALWAYS_ASSERT(m_hibc[idim] != BCType::Bogus);
        if (m_hibc[idim] == BCType::Dirichlet) {
            phys_bc.setHi(idim, Outflow);  // MacBndry treats it as Dirichlet
        } else {
            phys_bc.setHi(idim, Inflow);   // MacBndry treats it as Neumann
        }        
    }

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
                m_crse_sol_br[amrlev]->copyFrom(*m_coarse_data_for_bc, 0, 0, 0, 1);
            } else {
                m_crse_sol_br[amrlev]->setVal(0.0);
            }
            m_bndry_sol[amrlev]->setBndryValues(*m_crse_sol_br[amrlev], 0,
                                                bcdata, 0, 0, 1,
                                                m_coarse_data_crse_ratio, phys_bc);
            br_ref_ratio = m_coarse_data_crse_ratio;
        }
        else
        {
            m_bndry_sol[amrlev]->setBndryValues(bcdata,0,0,1,phys_bc);
            br_ref_ratio = 1;
        }
    }
    else
    {
        m_bndry_sol[amrlev]->setBndryValues(bcdata,0,0,1, m_amr_ref_ratio[amrlev-1], phys_bc);
        br_ref_ratio = m_amr_ref_ratio[amrlev-1];
    }

    const Real* dx = m_geom[amrlev][0].CellSize();
    for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
    {
        m_bcondloc[amrlev][mglev]->setBndryConds(m_geom[amrlev][mglev], dx,
                                                 phys_bc, br_ref_ratio);
    }
}

void
MLLinOp::updateSolBC (int amrlev, const MultiFab& crse_bcdata) const
{
    BL_PROFILE("MLLinOp::updateSolBC()");

    AMREX_ALWAYS_ASSERT(amrlev > 0);
    m_crse_sol_br[amrlev]->copyFrom(crse_bcdata, 0, 0, 0, 1);
    m_bndry_sol[amrlev]->updateBndryValues(*m_crse_sol_br[amrlev], 0, 0, 1, m_amr_ref_ratio[amrlev-1]);
}

void
MLLinOp::updateCorBC (int amrlev, const MultiFab& crse_bcdata) const
{
    BL_PROFILE("MLLinOp::updateCorBC()");
    AMREX_ALWAYS_ASSERT(amrlev > 0);
    m_crse_cor_br[amrlev]->copyFrom(crse_bcdata, 0, 0, 0, 1);
    m_bndry_cor[amrlev]->updateBndryValues(*m_crse_cor_br[amrlev], 0, 0, 1, m_amr_ref_ratio[amrlev-1]);
}

void
MLLinOp::solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                           const MultiFab* crse_bcdata)
{
    BL_PROFILE("MLLinOp::solutionResidual()");
    if (crse_bcdata != nullptr) {
        updateSolBC(amrlev, *crse_bcdata);
    }
    const int mglev = 0;
    apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, m_bndry_sol[amrlev].get());
    MultiFab::Xpay(resid, -1.0, b, 0, 0, resid.nComp(), 0);
}

void
MLLinOp::correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                             BCMode bc_mode, const MultiFab* crse_bcdata)
{
    BL_PROFILE("MLLinOp::correctionResidual()");
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
MLLinOp::apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                const MacBndry* bndry) const
{
    BL_PROFILE("MLLinOp::apply()");
    applyBC(amrlev, mglev, in, bc_mode, bndry);
    Fapply(amrlev, mglev, out, in);
}

void
MLLinOp::applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode,
                  const MacBndry* bndry, bool skip_fillboundary) const
{
    BL_PROFILE("MLLinOp::applyBC()");
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
MLLinOp::smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                 bool skip_fillboundary) const
{
    BL_PROFILE("MLLinOp::smooth()");
    for (int redblack = 0; redblack < 2; ++redblack)
    {
        applyBC(amrlev, mglev, sol, BCMode::Homogeneous, nullptr, skip_fillboundary);
        Fsmooth(amrlev, mglev, sol, rhs, redblack);
        skip_fillboundary = false;
    }
}

void
MLLinOp::reflux (int crse_amrlev, MultiFab& res,
                 const MultiFab& crse_sol, MultiFab& fine_sol) const
{
    BL_PROFILE("MLLinOp::reflux()");
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
MLLinOp::compFlux (int amrlev, const std::array<MultiFab*,AMREX_SPACEDIM>& fluxes, MultiFab& sol) const
{
    BL_PROFILE("MLLinOp::compFlux()");

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
MLLinOp::compGrad (int amrlev, const std::array<MultiFab*,AMREX_SPACEDIM>& grad, MultiFab& sol) const
{
    BL_PROFILE("MLLinOp::compGrad()");

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
MLLinOp::prepareForSolve ()
{
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

MLLinOp::BndryCondLoc::BndryCondLoc (const BoxArray& ba, const DistributionMapping& dm)
    : bcond(ba, dm),
      bcloc(ba, dm)
{
}

void
MLLinOp::BndryCondLoc::setBndryConds (const Geometry& geom, const Real* dx,
                                      const BCRec& phys_bc, int ratio)
{
    // Same as MacBndry::setBndryConds

    const Box&  domain = geom.Domain();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(bcloc); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        RealTuple & bloc  = bcloc[mfi];
        BCTuple   & bctag = bcond[mfi];
        for (OrientationIter fi; fi; ++fi)
        {
            const Orientation face = fi();
            const int         dir  = face.coordDir();

            if (domain[face] == bx[face] && !geom.isPeriodic(dir))
            {
                //
                // All physical bc values are located on face.
                //
                const int p_bc  = (face.isLow() ? phys_bc.lo(dir) : phys_bc.hi(dir));

                bctag[face] = (p_bc == Outflow) ? LO_DIRICHLET : LO_NEUMANN;
                bloc[face]  = 0;
            }
            else
            {
                //
                // Internal bndry.
                //
                const Real delta = dx[dir]*ratio;

                bctag[face] = LO_DIRICHLET;
		bloc[face]  = 0.5*delta;
            }

        }
    }
}

DistributionMapping
MLLinOp::makeConsolidatedDMap (const BoxArray& ba)
{
    const std::vector< std::vector<int> >& sfc = DistributionMapping::makeSFC(ba);

    const int nprocs = ParallelDescriptor::NProcs();
    AMREX_ASSERT(static_cast<int>(sfc.size()) == nprocs);
    const int nboxes = ba.size();

    Vector<int> pmap(ba.size());
    if (nboxes >= nprocs)
    {
        for (int iproc = 0; iproc < nprocs; ++iproc) {
            for (int ibox : sfc[iproc]) {
                pmap[ibox] = iproc;
            }
        }
    }
    else
    {
        for (int i = 0; i < nboxes; ++i) { // after nboxes sfc[i] is empty
            for (int ibox : sfc[i]) {
                const int iproc = i;
                pmap[ibox] = iproc;
            }
        }
    }

    return DistributionMapping{pmap};
}

DistributionMapping
MLLinOp::makeConsolidatedDMap (const DistributionMapping& fdm)
{
    Vector<int> pmap = fdm.ProcessorMap();
    for (auto& x: pmap) {
        x /= 2;
    }
    return DistributionMapping{pmap};
}

MPI_Comm
MLLinOp::makeSubCommunicator (const DistributionMapping& dm)
{
#ifdef BL_USE_MPI
    MPI_Comm newcomm;
    MPI_Group defgrp, newgrp;

    MPI_Comm_group(m_default_comm, &defgrp);

    Array<int> newgrp_ranks = dm.ProcessorMap();
    std::sort(newgrp_ranks.begin(), newgrp_ranks.end());
    auto last = std::unique(newgrp_ranks.begin(), newgrp_ranks.end());
    newgrp_ranks.erase(last, newgrp_ranks.end());
    
    MPI_Group_incl(defgrp, newgrp_ranks.size(), newgrp_ranks.data(), &newgrp);

    MPI_Comm_create(m_default_comm, newgrp, &newcomm);   
    m_raii_comm.reset(new CommContainer(newcomm));

    MPI_Group_free(&defgrp);
    MPI_Group_free(&newgrp);

    return newcomm;
#else
    return m_default_comm;
#endif
}

void
MLLinOp::applyMetricTerm (int amrlev, int mglev, MultiFab& rhs) const
{
#if (AMREX_SPACEDIM != 3)
    
    if (Geometry::IsCartesian()) return;

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
MLLinOp::unapplyMetricTerm (int amrlev, int mglev, MultiFab& rhs) const
{
#if (AMREX_SPACEDIM != 3)

    if (Geometry::IsCartesian()) return;

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

MLLinOp::MetricFactor::MetricFactor (const BoxArray& ba, const DistributionMapping& dm,
                                     const Geometry& geom)
    : r_cellcenter(ba,dm),
      r_celledge(ba,dm),
      inv_r_cellcenter(ba,dm),
      inv_r_celledge(ba,dm)
{
    bool is_cart = Geometry::IsCartesian();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(r_cellcenter); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        if (is_cart)
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
        }
    }
}

#endif

}
