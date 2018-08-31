
#include <cmath>
#include <algorithm>
#include <AMReX_MLLinOp.H>
#include <AMReX_ParmParse.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>
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
#if AMREX_USE_EB
    if (!a_factory.empty()){
        auto f = dynamic_cast<EBFArrayBoxFactory const*>(a_factory[0]);
        if (f) {
            info.max_coarsening_level = std::min(info.max_coarsening_level,
                                                 f->maxCoarseningLevel());
        }
    }
#endif
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

    m_default_comm = ParallelContext::CommunicatorSub();

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
        for (int i = 0; i < 2; ++i)
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

        // We may have to agglomerate earlier because the original
        // BoxArray has to be coarsenable to the first agglomerated
        // level or the bottom level and in amrex::average_down the
        // fine BoxArray needs to be coarsenable (unless we make
        // average_down more general).
        int last_coarsenableto_lev = 0;
        for (int lev = std::min(nmaxlev,first_agglev); lev >= 1; --lev) {
            int ratio = static_cast<int>(std::pow(rr,lev));
            if (a_grids[0].coarsenable(ratio, mg_box_min_width)) {
                last_coarsenableto_lev = lev;
                break;
            }
        }

        // We now know we could coarsen the original BoxArray to at
        // least last_coarsenableto_lev.  last_coarsenableto_lev == 0
        // means the original BoxArray is not coarsenable even once.

        if (last_coarsenableto_lev > 0)
        {
            // last_coarsenableto_lev will be the first agglomeration level, except
            if (last_coarsenableto_lev == nmaxlev-1 && first_agglev > nmaxlev-1) {
                // then there is no reason to agglomerate
                last_coarsenableto_lev = nmaxlev;
            }

            for (int lev = 1; lev < last_coarsenableto_lev; ++lev)
            {
                m_geom[0].emplace_back(amrex::coarsen(a_geom[0].Domain(),rr));
                
                m_grids[0].push_back(a_grids[0]);
                m_grids[0].back().coarsen(rr);
            
                m_dmap[0].push_back(a_dmap[0]);
                
                rr *= mg_coarsen_ratio;
            }

            for (int lev = last_coarsenableto_lev; lev < nmaxlev; ++lev)
            {
                m_geom[0].emplace_back(domainboxes[lev]);
            
                m_grids[0].emplace_back(boundboxes[lev]);
                m_grids[0].back().maxSize(info.agg_grid_size);

                m_dmap[0].push_back(DistributionMapping());
                agged = true;
            }

            m_num_mg_levels[0] = m_grids[0].size();
        }
    }
    else
    {
        int rr = mg_coarsen_ratio;
        Real avg_npts, threshold_npts;
        if (info.do_consolidation) {
            avg_npts = static_cast<Real>(a_grids[0].d_numPts()) / static_cast<Real>(ParallelContext::NProcsSub());
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
            m_factory[amrlev].emplace_back(makeFactory(amrlev,mglev));
        }
    }

    for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
    {
        AMREX_ASSERT_WITH_MESSAGE(m_grids[amrlev][0].coarsenable(m_amr_ref_ratio[amrlev-1]),
                                  "MLLinOp: grids not coarsenable between AMR levels");
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
            mf[alev][mlev].define(ba, m_dmap[alev][mlev], nc, ng, MFInfo(), *m_factory[alev][mlev]);
        }
    }
}

void
MLLinOp::setDomainBC (const Array<BCType,AMREX_SPACEDIM>& a_lobc,
                      const Array<BCType,AMREX_SPACEDIM>& a_hibc)
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
    
    if (ParallelContext::CommunicatorSub() == ParallelDescriptor::Communicator()) {
        MPI_Group_incl(defgrp, newgrp_ranks.size(), newgrp_ranks.data(), &newgrp);
    } else {
        Vector<int> local_newgrp_ranks(newgrp_ranks.size());
        ParallelContext::global_to_local_rank(local_newgrp_ranks.data(),
                                              newgrp_ranks.data(), newgrp_ranks.size());
        MPI_Group_incl(defgrp, local_newgrp_ranks.size(), local_newgrp_ranks.data(), &newgrp);
    }

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
            
            const int nprocs = ParallelContext::NProcsSub();
            AMREX_ASSERT(static_cast<int>(sfc.size()) == nprocs);
            
            Vector<int> pmap(ba[i].size());
            for (int iproc = 0; iproc < nprocs; ++iproc) {
                int grank = ParallelContext::local_to_global_rank(iproc);
                for (int ibox : sfc[iproc]) {
                    pmap[ibox] = grank;
                }
            }
            dm[i].define(std::move(pmap));
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

            const int nprocs = ParallelContext::NProcsSub();
            const auto& pmap_fine = dm[i-1].ProcessorMap();
            Vector<int> pmap(pmap_fine.size());
            ParallelContext::global_to_local_rank(pmap.data(), pmap_fine.data(), pmap.size()); 
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

            if (ParallelContext::CommunicatorSub() == ParallelDescriptor::Communicator()) {
                dm[i].define(std::move(pmap));
            } else {
                Vector<int> pmap_g(pmap.size());
                ParallelContext::local_to_global_rank(pmap_g.data(), pmap.data(), pmap.size());
                dm[i].define(std::move(pmap_g));
            }
        }
    }
}

}
