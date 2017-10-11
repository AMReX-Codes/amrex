
#include <AMReX_MLLinOp.H>
#include <AMReX_LO_F.H>

namespace amrex {

constexpr int MLLinOp::mg_coarsen_ratio;
constexpr int MLLinOp::mg_box_min_width;
constexpr int MLLinOp::maxorder;

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
    m_num_amr_levels = a_geom.size();

    m_amr_ref_ratio.resize(m_num_amr_levels);
    m_num_mg_levels.resize(m_num_amr_levels);

    m_geom.resize(m_num_amr_levels);
    m_grids.resize(m_num_amr_levels);
    m_dmap.resize(m_num_amr_levels);

    m_undrrelxr.resize(m_num_amr_levels);

    m_maskvals.resize(m_num_amr_levels);

    m_macbndry.resize(m_num_amr_levels);

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

    int rr = mg_coarsen_ratio;
    while (a_geom[0].Domain().coarsenable(rr)
           and a_grids[0].coarsenable(rr, mg_box_min_width))
    {
        ++(m_num_mg_levels[0]);

        m_geom[0].emplace_back(amrex::coarsen(a_geom[0].Domain(),rr));

        m_grids[0].push_back(a_grids[0]);
        m_grids[0].back().coarsen(rr);

        m_dmap[0].push_back(a_dmap[0]);
        
        rr *= mg_coarsen_ratio;
    }

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

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_macbndry[amrlev].reset(new MacBndry(m_grids[amrlev][0], m_dmap[amrlev][0],
                                              1, m_geom[amrlev][0]));
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
MLLinOp::setBC (int amrlev, const MultiFab* bc_data, const MultiFab* crse_bcdata)
{
    // xxxxx this will make it Dirichlet
    BCRec phys_bc{AMREX_D_DECL(Outflow,Outflow,Outflow),
                  AMREX_D_DECL(Outflow,Outflow,Outflow)};

    if (bc_data == nullptr && crse_bcdata == nullptr)
    {   // xxxxx ?????
        m_macbndry[amrlev]->setHomogValues(phys_bc, IntVect::TheZeroVector());
    }
    else if (bc_data && crse_bcdata == nullptr)
    {
        AMREX_ALWAYS_ASSERT(amrlev == 0);
        m_macbndry[amrlev]->setBndryValues(*bc_data,0,0,1,phys_bc);
    }
    else if (bc_data && crse_bcdata)
    {
        const int ncomp = 1;
        const int in_rad = 0;
        const int out_rad = 1;
        const int extent_rad = 2;
        const int crse_ratio = m_amr_ref_ratio[amrlev-1];

        BoxArray cba = m_grids[amrlev][0];
        cba.coarsen(crse_ratio);
        BndryRegister crse_br(cba, m_dmap[amrlev][0], in_rad, out_rad, extent_rad, ncomp);
        crse_br.copyFrom(*crse_bcdata, crse_bcdata->nGrow(), 0, 0, ncomp);

        m_macbndry[amrlev]->setBndryValues(crse_br, 0, *bc_data, 0, 0, ncomp, crse_ratio, phys_bc); 
    }
    else
    {
        amrex::Abort("MLLinOp::setBC(): don't know what to do");
    }
}

void
MLLinOp::updateBC (int amrlev, const MultiFab& crse_bcdata)
{
    AMREX_ALWAYS_ASSERT(amrlev > 0);

    // xxxxx this will make it Dirichlet
    BCRec phys_bc{AMREX_D_DECL(Outflow,Outflow,Outflow),
                  AMREX_D_DECL(Outflow,Outflow,Outflow)};

    static bool first = true;
    if (first) {
        if (ParallelDescriptor::IOProcessor()) {
            amrex::Warning("MLLinOp::updateBC() not implemented");
        }
        first = false;
    }
}

void
MLLinOp::residual (int amrlev, int mglev,
                   MultiFab& resid, MultiFab& sol, const MultiFab& rhs,
                   BCMode bc_mode) const
{
    apply(amrlev, mglev, resid, sol, bc_mode);
    MultiFab::Xpay(resid, -1.0, rhs, 0, 0, resid.nComp(), 0);
}

void
MLLinOp::apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode) const
{
    applyBC(amrlev, mglev, in, bc_mode);
    Fapply(amrlev, mglev, out, in);
}

void
MLLinOp::applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode) const
{
    // No coarsened boundary values, cannot apply inhomog at mglev>0.
    BL_ASSERT(mglev == 0 || bc_mode == BCMode::Homogeneous);

    const bool cross = true;
    in.FillBoundary(0, 1, m_geom[amrlev][mglev].periodicity(), cross);

    int flagbc = (bc_mode == BCMode::Homogeneous) ? 0 : 1;

    int flagden = 0;

    const int num_comp = 1;
    const Real* h = m_geom[amrlev][mglev].CellSize();

    const MacBndry& macbndry = *m_macbndry[amrlev];
    BndryRegister& undrrelxr = m_undrrelxr[amrlev][mglev];
    const auto& maskvals = m_maskvals[amrlev][mglev];

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(in, MFItInfo().SetDynamic(true));
         mfi.isValid(); ++mfi)
    {
        const BndryData::RealTuple      & bdl = macbndry.bndryLocs(mfi);
        const Vector<Vector<BoundCond> >& bdc = macbndry.bndryConds(mfi);

        for (OrientationIter oitr; oitr; ++oitr)
        {
            const Orientation ori = oitr();

            int           cdr = ori;
            Real          bcl = bdl[ori];
            int           bct = bdc[ori][0];

            const Box&       vbx   = mfi.validbox();
            FArrayBox&       iofab = in[mfi];

            FArrayBox&       ffab  = undrrelxr[ori][mfi];
            const FArrayBox& fsfab = macbndry.bndryValues(ori)[mfi];

            const Mask&   m   = maskvals[ori][mfi];

            FORT_APPLYBC(&flagden, &flagbc, &maxorder,
                         iofab.dataPtr(),
                         ARLIM(iofab.loVect()), ARLIM(iofab.hiVect()),
                         &cdr, &bct, &bcl,
                         fsfab.dataPtr(), 
                         ARLIM(fsfab.loVect()), ARLIM(fsfab.hiVect()),
                         m.dataPtr(),
                         ARLIM(m.loVect()), ARLIM(m.hiVect()),
                         ffab.dataPtr(),
                         ARLIM(ffab.loVect()), ARLIM(ffab.hiVect()),
                         vbx.loVect(),
                         vbx.hiVect(), &num_comp, h);
        }
    }

}

void
MLLinOp::smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, BCMode bc_mode) const
{
    for (int redblack = 0; redblack < 2; ++redblack)
    {
        applyBC(amrlev, mglev, sol, bc_mode);
        Fsmooth(amrlev, mglev, sol, rhs, redblack);
    }
}

void
MLLinOp::reflux (int crse_amrlev, MultiFab& res,
                 const MultiFab& crse_sol, const MultiFab& fine_sol) const
{

}

}
