
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <AMReX_Utility.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_MLCellLinOp.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Machine.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>
#endif

#ifdef AMREX_USE_PETSC
#include <petscksp.h>
#include <AMReX_PETSc.H>
#endif

namespace amrex {

constexpr int MLLinOp::mg_coarsen_ratio;
constexpr int MLLinOp::mg_box_min_width;

namespace {
    // experimental features
    bool initialized = false; // track initialization of static state
    int consolidation_threshold = -1;
    int consolidation_ratio = 2;
    int consolidation_strategy = 3;

    int flag_verbose_linop = 0;
    int flag_comm_cache = 0;
    int flag_use_mota = 0;
    int remap_nbh_lb = 1;

#ifdef BL_USE_MPI
    class CommCache
    {
      public:
        CommCache () = default;
        CommCache (const CommCache&) = delete;
        CommCache (CommCache&&) = delete;
        void operator= (const CommCache&) = delete;
        void operator= (CommCache&&) = delete;

        ~CommCache () {
            for (auto & p : cache) {
                if (p.second != MPI_COMM_NULL) {
                    MPI_Comm_free(&p.second);
                }
            }
        }
        void add (size_t key, MPI_Comm comm) {
            AMREX_ASSERT(cache.count(key) == 0);
            cache[key] = comm;
        }
        bool get (size_t key, MPI_Comm &comm) {
            bool result = cache.count(key) > 0;
            if (result) {
                comm = cache.at(key);
            }
            return result;
        }

      private:
        std::unordered_map<size_t, MPI_Comm> cache;
    };

    std::unique_ptr<CommCache> comm_cache;
#endif

    Vector<int> get_subgroup_ranks ()
    {
        int rank_n = ParallelContext::NProcsSub();
        Vector<int> lranks(rank_n);
        for (int i = 0; i < rank_n; ++i) {
            lranks[i] = i;
        }

        Vector<int> granks(rank_n);
        ParallelContext::local_to_global_rank(granks.data(), lranks.data(), rank_n);
        return granks;
    }
}

// static member function
void MLLinOp::Initialize ()
{
    ParmParse pp("mg");
    pp.query("consolidation_threshold", consolidation_threshold);
    pp.query("consolidation_ratio", consolidation_ratio);
    pp.query("consolidation_strategy", consolidation_strategy);
    pp.query("verbose_linop", flag_verbose_linop);
    pp.query("comm_cache", flag_comm_cache);
    pp.query("mota", flag_use_mota);
    pp.query("remap_nbh_lb", remap_nbh_lb);

#ifdef BL_USE_MPI
    comm_cache.reset(new CommCache());
#endif
    amrex::ExecOnFinalize(MLLinOp::Finalize);
    initialized = true;
}

// static member function
void MLLinOp::Finalize ()
{
    initialized = false;
#ifdef BL_USE_MPI
    comm_cache.reset();
#endif
#ifdef AMREX_SOFT_PERF_COUNTERS
    MLCellLinOp::perf_counters.reset();
#endif
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
        Initialize();
    }

    info = a_info;
#ifdef AMREX_USE_EB
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
    int agg_lev, con_lev;

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
                agg_lev = last_coarsenableto_lev;
            }

            m_num_mg_levels[0] = m_grids[0].size();
        }
    }
    else
    {
        int rr = mg_coarsen_ratio;
        Real avg_npts;
        if (info.do_consolidation) {
            avg_npts = static_cast<Real>(a_grids[0].d_numPts()) / static_cast<Real>(ParallelContext::NProcsSub());
            if (consolidation_threshold == -1) {
                consolidation_threshold = static_cast<Real>(AMREX_D_TERM(info.con_grid_size,
                                                                         *info.con_grid_size,
                                                                         *info.con_grid_size));
            }
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
                if (avg_npts/(AMREX_D_TERM(rr,*rr,*rr)) < 0.999*consolidation_threshold)
                {
                    coned = true;
                    con_lev = m_dmap[0].size();
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

    if (flag_use_mota && (agged || coned))
    {
        remapNeighborhoods(m_dmap[0]);
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

    if (flag_verbose_linop) {
        if (agged) {
            Print() << "MLLinOp::defineGrids(): agglomerated AMR level 0 starting at MG level "
                    << agg_lev << " of " << m_num_mg_levels[0] << std::endl;
        } else if (coned) {
            Print() << "MLLinOp::defineGrids(): consolidated AMR level 0 starting at MG level "
                    << con_lev << " of " << m_num_mg_levels[0]
                    << " (ratio = " << consolidation_ratio << ")" << std::endl;
        } else {
            Print() << "MLLinOp::defineGrids(): no agglomeration or consolidation of AMR level 0" << std::endl;
        }
    }

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
                      const Array<BCType,AMREX_SPACEDIM>& a_hibc) noexcept
{
    const int ncomp = getNComp();
    m_lobc.clear();
    m_hibc.clear();
    m_lobc.resize(ncomp,a_lobc);
    m_hibc.resize(ncomp,a_hibc);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (Geometry::isPeriodic(idim)) {
            AMREX_ALWAYS_ASSERT(a_lobc[idim] == BCType::Periodic);
            AMREX_ALWAYS_ASSERT(a_hibc[idim] == BCType::Periodic);
        }
        if (a_lobc[idim] == BCType::Periodic or
            a_hibc[idim] == BCType::Periodic) {
            AMREX_ALWAYS_ASSERT(Geometry::isPeriodic(idim));
        }
    }
}

void
MLLinOp::setDomainBC (const Vector<Array<BCType,AMREX_SPACEDIM> >& a_lobc,
                      const Vector<Array<BCType,AMREX_SPACEDIM> >& a_hibc) noexcept
{
    const int ncomp = getNComp();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncomp == a_lobc.size() && ncomp == a_hibc.size(),
                                     "MLLinOp::setDomainBC: wrong size");
    m_lobc = a_lobc;
    m_hibc = a_hibc;
    for (int icomp = 0; icomp < ncomp; ++icomp) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (Geometry::isPeriodic(idim)) {
                AMREX_ALWAYS_ASSERT(m_lobc[icomp][idim] == BCType::Periodic);
                AMREX_ALWAYS_ASSERT(m_hibc[icomp][idim] == BCType::Periodic);
            }
            if (m_lobc[icomp][idim] == BCType::Periodic or
                m_hibc[icomp][idim] == BCType::Periodic) {
                AMREX_ALWAYS_ASSERT(Geometry::isPeriodic(idim));
            }
        }
    }
}

void
MLLinOp::setCoarseFineBC (const MultiFab* crse, int crse_ratio) noexcept
{
    m_coarse_data_for_bc = crse;
    m_coarse_data_crse_ratio = crse_ratio;
}

MPI_Comm
MLLinOp::makeSubCommunicator (const DistributionMapping& dm)
{
    BL_PROFILE("MLLinOp::makeSubCommunicator()");

#ifdef BL_USE_MPI

    Vector<int> newgrp_ranks = dm.ProcessorMap();
    std::sort(newgrp_ranks.begin(), newgrp_ranks.end());
    auto last = std::unique(newgrp_ranks.begin(), newgrp_ranks.end());
    newgrp_ranks.erase(last, newgrp_ranks.end());
    
    if (flag_verbose_linop) {
        Print() << "MLLinOp::makeSubCommunicator(): called for " << newgrp_ranks.size() << " ranks" << std::endl;
    }

    MPI_Comm newcomm;
    bool cache_hit = false;
    uint64_t key = 0;
    if (flag_comm_cache) {
        AMREX_ASSERT(comm_cache);
        key = hash_vector(newgrp_ranks, hash_vector(get_subgroup_ranks()));
        cache_hit = comm_cache->get(key, newcomm);
        if (cache_hit && flag_verbose_linop) {
            Print() << "MLLinOp::makeSubCommunicator(): found subcomm in cache" << std::endl;
        }
    }

    if (!flag_comm_cache || !cache_hit) {
        MPI_Group defgrp, newgrp;
        MPI_Comm_group(m_default_comm, &defgrp);
        if (ParallelContext::CommunicatorSub() == ParallelDescriptor::Communicator()) {
            MPI_Group_incl(defgrp, newgrp_ranks.size(), newgrp_ranks.data(), &newgrp);
        } else {
            Vector<int> local_newgrp_ranks(newgrp_ranks.size());
            ParallelContext::global_to_local_rank(local_newgrp_ranks.data(),
                                                  newgrp_ranks.data(), newgrp_ranks.size());
            MPI_Group_incl(defgrp, local_newgrp_ranks.size(), local_newgrp_ranks.data(), &newgrp);
        }

        if (flag_verbose_linop) {
            Print() << "MLLinOp::makeSubCommunicator(): MPI_Comm_create: (";
            for (auto rank : newgrp_ranks) {
                Print() << rank << ",";
            }
            Print() << "\b)" << std::endl;
        }

        MPI_Comm_create(m_default_comm, newgrp, &newcomm);

        if (flag_comm_cache) {
            comm_cache->add(key, newcomm);
        } else {
            m_raii_comm.reset(new CommContainer(newcomm));
        }

        MPI_Group_free(&defgrp);
        MPI_Group_free(&newgrp);
    }

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

void
MLLinOp::remapNeighborhoods (Vector<DistributionMapping> & dms)
{
    BL_PROFILE("MLLinOp::remapNeighborhoods()");

    if (flag_verbose_linop) {
        Print() << "Remapping ranks to neighborhoods ..." << std::endl;
    }

    for (int j = 1; j < dms.size(); ++j)
    {
        const Vector<int> & pmap = dms[j].ProcessorMap();
        std::set<int> g_ranks_set(pmap.begin(), pmap.end());
        auto lev_rank_n = g_ranks_set.size();
        if (lev_rank_n >= remap_nbh_lb && lev_rank_n < ParallelContext::NProcsSub())
        {
            // find best neighborhood with lev_rank_n ranks
            auto nbh_g_ranks = machine::find_best_nbh(lev_rank_n);
            AMREX_ASSERT(nbh_g_ranks.size() == lev_rank_n);

            // construct mapping from original global rank to neighborhood global rank
            int idx = 0;
            std::unordered_map<int, int> rank_mapping;
            for (auto orig_g_rank : g_ranks_set) {
                AMREX_ASSERT(idx < nbh_g_ranks.size());
                rank_mapping[orig_g_rank] = nbh_g_ranks[idx++];
            }

            // remap and create new DM
            Vector<int> nbh_pmap(pmap.size());
            for (int i = 0; i < pmap.size(); ++i) {
                nbh_pmap[i] = rank_mapping.at(pmap[i]);
            }
            dms[j] = DistributionMapping(std::move(nbh_pmap));
        }
    }
}

#ifdef AMREX_USE_PETSC
std::unique_ptr<PETScABecLap>
MLLinOp::makePETSc () const
{
    amrex::Abort("MLLinOp::makePETSc: How did we get here?");
    return {nullptr};
}
#endif

}
