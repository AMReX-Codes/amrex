
#include <AMReX_MLLinOp.H>

namespace amrex {

constexpr int MLLinOp::mg_coarsen_ratio;
constexpr int MLLinOp::mg_box_min_width;

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
MLLinOp::residual (int amrlev, int mglev,
                   MultiFab& resid, MultiFab& sol, const MultiFab& rhs,
                   BCMode bc_mode)
{
    apply(amrlev, mglev, resid, sol, bc_mode);
    MultiFab::Xpay(resid, -1.0, rhs, 0, 0, resid.nComp(), 0);
}

void
MLLinOp::apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode)
{
    applyBC(amrlev, mglev, in, bc_mode);
    Fapply(amrlev, mglev, out, in);
}

void
MLLinOp::applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode)
{
    
}

}
