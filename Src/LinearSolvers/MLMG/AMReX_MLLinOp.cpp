
#include <AMReX_MLLinOp.H>
#include <AMReX_ParmParse.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBTower.H>
#endif

namespace amrex {

constexpr int MLLinOp::mg_coarsen_ratio;
constexpr int MLLinOp::mg_box_min_width;

namespace {
    // experimental features
    bool initialized = false;
    int consolidation_ratio = 2;
    int consolidation_strategy = 3;
}

MLLinOp::MLLinOp () {}

MLLinOp::~MLLinOp () {}

void
MLLinOp::define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("MLLinOp::define()");

    if (!initialized) {
	ParmParse pp("mg");
	pp.query("consolidation_ratio", consolidation_ratio);
	pp.query("consolidation_strategy", consolidation_strategy);
	initialized = true;
    }

    info = a_info;
    defineGrids(a_geom, a_grids, a_dmap, a_factory);
    defineAuxData();
    defineBC();
}

void
MLLinOp::defineGrids (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap,
                      const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("MLLinOp::defineGrids()");

    m_num_amr_levels = a_geom.size();

    m_amr_ref_ratio.resize(m_num_amr_levels);
    m_num_mg_levels.resize(m_num_amr_levels);

    m_geom.resize(m_num_amr_levels);
    m_grids.resize(m_num_amr_levels);
    m_dmap.resize(m_num_amr_levels);
    m_factory.resize(m_num_amr_levels);

    m_default_comm = ParallelDescriptor::Communicator();

    // fine amr levels
    for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
        m_num_mg_levels[amrlev] = 1;
        m_geom[amrlev].push_back(a_geom[amrlev]);
        m_grids[amrlev].push_back(a_grids[amrlev]);
        m_dmap[amrlev].push_back(a_dmap[amrlev]);
        if (amrlev < a_factory.size()) {
            m_factory[amrlev].emplace_back(a_factory[amrlev]->clone());
        } else {
            m_factory[amrlev].emplace_back(new FArrayBoxFactory());
        }

        int rr = mg_coarsen_ratio;
        const Box& dom = a_geom[amrlev].Domain();
        for (int i = 0; i < info.max_coarsening_level; ++i)
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
    if (a_factory.size() > 0) {
        m_factory[0].emplace_back(a_factory[0]->clone());
    } else {
        m_factory[0].emplace_back(new FArrayBoxFactory());
    }

    m_domain_covered.resize(m_num_amr_levels, false);
    auto npts0 = m_grids[0][0].numPts();
    m_domain_covered[0] = (npts0 == m_geom[0][0].Domain().numPts());
    for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
    {
        if (!m_domain_covered[amrlev-1]) break;
        m_domain_covered[amrlev] = (m_grids[amrlev][0].numPts() == m_geom[amrlev][0].Domain().numPts());
    }

    Box aggbox;
    bool aggable = false;

    if (info.do_agglomeration)
    {
        if (m_domain_covered[0])
        {
            aggbox = m_geom[0][0].Domain();
            aggable = true;
        }
        else
        {
            aggbox = m_grids[0][0].minimalBox();
            aggable = (aggbox.numPts() == npts0);
        }
    }

    bool agged = false;
    bool coned = false;

    if (info.do_agglomeration && aggable)
    {
        Vector<Box> domainboxes;
        Vector<Box> boundboxes;
        Box dbx = m_geom[0][0].Domain();
        Box bbx = aggbox;
        Real nbxs = static_cast<Real>(m_grids[0][0].size());
        Real threshold_npts = static_cast<Real>(AMREX_D_TERM(info.agg_grid_size,
                                                             *info.agg_grid_size,
                                                             *info.agg_grid_size));
        Vector<int> agg_flag;
        domainboxes.push_back(dbx);
        boundboxes.push_back(bbx);
        agg_flag.push_back(false);
        while (    dbx.coarsenable(mg_coarsen_ratio,mg_box_min_width)
               and bbx.coarsenable(mg_coarsen_ratio,mg_box_min_width))
        {
            dbx.coarsen(mg_coarsen_ratio);
            domainboxes.push_back(dbx);
            bbx.coarsen(mg_coarsen_ratio);
            boundboxes.push_back(bbx);
            bool to_agg = (bbx.d_numPts() / nbxs) < 0.999*threshold_npts;
            agg_flag.push_back(to_agg);
        }

        int first_agglev = std::distance(agg_flag.begin(),
                                         std::find(agg_flag.begin(),agg_flag.end(),1));
        int nmaxlev = std::min(static_cast<int>(domainboxes.size()),
                               info.max_coarsening_level + 1);
        int rr = mg_coarsen_ratio;
        for (int lev = 1; lev < nmaxlev; ++lev)
        {
            if (lev >= first_agglev or !a_grids[0].coarsenable(rr,mg_box_min_width))
            {
                m_geom[0].emplace_back(domainboxes[lev]);
            
                m_grids[0].emplace_back(boundboxes[lev]);
                m_grids[0].back().maxSize(info.agg_grid_size);

                m_dmap[0].push_back(DistributionMapping());
                agged = true;
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
        if (info.do_consolidation) {
            avg_npts = static_cast<Real>(a_grids[0].d_numPts()) / static_cast<Real>(ParallelDescriptor::NProcs());
            threshold_npts = static_cast<Real>(AMREX_D_TERM(info.con_grid_size,
                                                            *info.con_grid_size,
                                                            *info.con_grid_size));
        }
        while (m_num_mg_levels[0] < info.max_coarsening_level + 1
               and a_geom[0].Domain().coarsenable(rr)
               and a_grids[0].coarsenable(rr, mg_box_min_width))
        {
            m_geom[0].emplace_back(amrex::coarsen(a_geom[0].Domain(),rr));
            
            m_grids[0].push_back(a_grids[0]);
            m_grids[0].back().coarsen(rr);

            if (info.do_consolidation)
            {
                if (avg_npts/(AMREX_D_TERM(rr,*rr,*rr)) < 0.999*threshold_npts)
                {
                    coned = true;
                    m_dmap[0].push_back(DistributionMapping());
                }
                else
                {
                    m_dmap[0].push_back(m_dmap[0].back());
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

    if (agged)
    {
        makeAgglomeratedDMap(m_grids[0], m_dmap[0]);
    }
    else if (coned)
    {
        makeConsolidatedDMap(m_grids[0], m_dmap[0], consolidation_ratio, consolidation_strategy);
    }

    if (info.do_agglomeration || info.do_consolidation)
    {
        m_bottom_comm = makeSubCommunicator(m_dmap[0].back());
    }
    else
    {
        m_bottom_comm = m_default_comm;
    }

    m_do_agglomeration = agged;
    m_do_consolidation = coned;

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_factory[amrlev].emplace_back(new FArrayBoxFactory());
        }
    }
}

void
MLLinOp::defineAuxData ()
{
}

void
MLLinOp::defineBC ()
{
    m_needs_coarse_data_for_bc = !m_domain_covered[0];
}

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
            const auto& ba = amrex::convert(m_grids[alev][mlev], m_ixtype);
            mf[alev][mlev].define(ba, m_dmap[alev][mlev], nc, ng);
        }
    }
}

void
MLLinOp::setDomainBC (const std::array<BCType,AMREX_SPACEDIM>& a_lobc,
                      const std::array<BCType,AMREX_SPACEDIM>& a_hibc)
{
    m_lobc = a_lobc;
    m_hibc = a_hibc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (Geometry::isPeriodic(idim)) {
            AMREX_ALWAYS_ASSERT(m_lobc[idim] == BCType::Periodic);
            AMREX_ALWAYS_ASSERT(m_hibc[idim] == BCType::Periodic);
        }
        if (m_lobc[idim] == BCType::Periodic or
            m_hibc[idim] == BCType::Periodic) {
            AMREX_ALWAYS_ASSERT(Geometry::isPeriodic(idim));
        }
    }
}

void
MLLinOp::setCoarseFineBC (const MultiFab* crse, int crse_ratio)
{
    m_coarse_data_for_bc = crse;
    m_coarse_data_crse_ratio = crse_ratio;
}

MPI_Comm
MLLinOp::makeSubCommunicator (const DistributionMapping& dm)
{
    BL_PROFILE("MLLinOp::makeSubCommunicator()");

#ifdef BL_USE_MPI
    MPI_Comm newcomm;
    MPI_Group defgrp, newgrp;

    MPI_Comm_group(m_default_comm, &defgrp);

    Vector<int> newgrp_ranks = dm.ProcessorMap();
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
MLLinOp::makeAgglomeratedDMap (const Vector<BoxArray>& ba, Vector<DistributionMapping>& dm)
{
    BL_PROFILE("MLLinOp::makeAgglomeratedDMap");

    BL_ASSERT(!dm[0].empty());
    for (int i = 1, N=ba.size(); i < N; ++i)
    {
        if (dm[i].empty())
        {
            const std::vector< std::vector<int> >& sfc = DistributionMapping::makeSFC(ba[i]);
            
            const int nprocs = ParallelDescriptor::NProcs();
            AMREX_ASSERT(static_cast<int>(sfc.size()) == nprocs);
            
            Vector<int> pmap(ba[i].size());
            for (int iproc = 0; iproc < nprocs; ++iproc) {
                for (int ibox : sfc[iproc]) {
                    pmap[ibox] = iproc;
                }
            }
            dm[i].define(pmap);
        }
    }
}


void
MLLinOp::makeConsolidatedDMap (const Vector<BoxArray>& ba, Vector<DistributionMapping>& dm,
                               int ratio, int strategy)
{
    BL_PROFILE("MLLinOp::makeConsolidatedDMap()");

    int factor = 1;
    BL_ASSERT(!dm[0].empty());
    for (int i = 1, N=ba.size(); i < N; ++i)
    {
        if (dm[i].empty())
        {
            factor *= ratio;

            const int nprocs = ParallelDescriptor::NProcs();
            Vector<int> pmap = dm[i-1].ProcessorMap();
            if (strategy == 1) {
                for (auto& x: pmap) {
                    x /= ratio;
                }
            } else if (strategy == 2) {
                int nprocs_con = static_cast<int>(std::ceil(static_cast<Real>(nprocs)
                                                            / static_cast<Real>(factor)));
                for (auto& x: pmap) {
                    auto d = std::div(x,nprocs_con);
                    x = d.rem;
                }
            } else if (strategy == 3) {
                if (factor == ratio) {
                    const std::vector< std::vector<int> >& sfc = DistributionMapping::makeSFC(ba[i]);
                    for (int iproc = 0; iproc < nprocs; ++iproc) {
                        for (int ibox : sfc[iproc]) {
                            pmap[ibox] = iproc;
                        }
                    }
                }
                for (auto& x: pmap) {
                    x /= ratio;
                }
            }
            dm[i].define(pmap);
        }
    }
}

void
MLLinOp::apply (const Vector<MultiFab*>& out, const Vector<MultiFab*>& in)
{
    Vector<MultiFab> rhs(m_num_amr_levels);
    for (int alev = m_num_amr_levels-1; alev >= 0; --alev)
    {
        rhs[alev].define(in[alev]->boxArray(), in[alev]->DistributionMap(),
                         in[alev]->nComp(), 0);
        rhs[alev].setVal(0.0);

        const MultiFab* crse_bcdata = (alev > 0) ? in[alev-1] : nullptr;

        solutionResidual(alev, *out[alev], *in[alev], rhs[alev], crse_bcdata);

        if (alev < m_num_amr_levels-1) {
            reflux(alev, *out[alev], *in[alev], rhs[alev],
                   *out[alev+1], *in[alev+1], rhs[alev+1]);
        }
    }

    for (int alev = 0; alev < m_num_amr_levels; ++alev) {
        out[alev]->negate();
    }
}

}
