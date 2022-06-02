
#include <AMReX_Utility.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_MLCellLinOp.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Machine.H>
#include <AMReX_MultiFabUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBMultiFabUtil.H>
#endif

#ifdef AMREX_USE_PETSC
#include <petscksp.h>
#include <AMReX_PETSc.H>
#endif

#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_map>


namespace amrex {

#if __cplusplus < 201703L
constexpr int MLLinOp::mg_coarsen_ratio;
constexpr int MLLinOp::mg_box_min_width;
constexpr int MLLinOp::mg_domain_min_width;
#endif

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
#endif
}

// static member function
void MLLinOp::Initialize ()
{
    ParmParse pp("mg");
    pp.queryAdd("consolidation_threshold", consolidation_threshold);
    pp.queryAdd("consolidation_ratio", consolidation_ratio);
    pp.queryAdd("consolidation_strategy", consolidation_strategy);
    pp.queryAdd("verbose_linop", flag_verbose_linop);
    pp.queryAdd("comm_cache", flag_comm_cache);
    pp.queryAdd("mota", flag_use_mota);
    pp.queryAdd("remap_nbh_lb", remap_nbh_lb);

#ifdef BL_USE_MPI
    comm_cache = std::make_unique<CommCache>();
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
                 const Vector<FabFactory<FArrayBox> const*>& a_factory,
                 bool eb_limit_coarsening)
{
    amrex::ignore_unused(eb_limit_coarsening);

    BL_PROFILE("MLLinOp::define()");

    if (!initialized) {
        Initialize();
    }

    info = a_info;
#ifdef AMREX_USE_GPU
    if (Gpu::notInLaunchRegion())
    {
        if (info.agg_grid_size <= 0) info.agg_grid_size = AMREX_D_PICK(32, 16, 8);
        if (info.con_grid_size <= 0) info.con_grid_size = AMREX_D_PICK(32, 16, 8);
    }
    else
#endif
    {
        if (info.agg_grid_size <= 0) info.agg_grid_size = LPInfo::getDefaultAgglomerationGridSize();
        if (info.con_grid_size <= 0) info.con_grid_size = LPInfo::getDefaultConsolidationGridSize();
    }

#ifdef AMREX_USE_EB
    if (!a_factory.empty() && eb_limit_coarsening) {
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

    const RealBox& rb = a_geom[0].ProbDomain();
    const int coord = a_geom[0].Coord();
    const Array<int,AMREX_SPACEDIM>& is_per = a_geom[0].isPeriodic();

    IntVect mg_coarsen_ratio_v(mg_coarsen_ratio);
    IntVect mg_box_min_width_v(mg_box_min_width);
    IntVect mg_domain_min_width_v(mg_domain_min_width);
    if (hasHiddenDimension()) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMREX_SPACEDIM == 3 && m_num_amr_levels == 1,
                                         "Hidden direction only supported for 3d level solve");
        mg_coarsen_ratio_v[info.hidden_direction] = 1;
        mg_box_min_width_v[info.hidden_direction] = 0;
        mg_domain_min_width_v[info.hidden_direction] = 0;
    }

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
            m_factory[amrlev].push_back(std::make_unique<FArrayBoxFactory>());
        }

        IntVect rr = mg_coarsen_ratio_v;
        const Box& dom = a_geom[amrlev].Domain();
        for (int i = 0; i < 2; ++i)
        {
            if (!dom.coarsenable(rr)) amrex::Abort("MLLinOp: Uncoarsenable domain");

            const Box& cdom = amrex::coarsen(dom,rr);
            if (cdom == a_geom[amrlev-1].Domain()) break;

            ++(m_num_mg_levels[amrlev]);

            m_geom[amrlev].emplace_back(cdom, rb, coord, is_per);

            m_grids[amrlev].push_back(a_grids[amrlev]);
            AMREX_ASSERT(m_grids[amrlev].back().coarsenable(rr));
            m_grids[amrlev].back().coarsen(rr);

            m_dmap[amrlev].push_back(a_dmap[amrlev]);

            rr *= mg_coarsen_ratio_v;
        }

        if (hasHiddenDimension()) {
            m_amr_ref_ratio[amrlev-1] = rr[AMREX_SPACEDIM-info.hidden_direction];
        } else {
            m_amr_ref_ratio[amrlev-1] = rr[0];
        }
    }

    // coarsest amr level
    m_num_mg_levels[0] = 1;
    m_geom[0].push_back(a_geom[0]);
    m_grids[0].push_back(a_grids[0]);
    m_dmap[0].push_back(a_dmap[0]);
    if (a_factory.size() > 0) {
        m_factory[0].emplace_back(a_factory[0]->clone());
    } else {
        m_factory[0].push_back(std::make_unique<FArrayBoxFactory>());
    }

    m_domain_covered.resize(m_num_amr_levels, false);
    auto npts0 = m_grids[0][0].numPts();
    m_domain_covered[0] = (npts0 == compactify(m_geom[0][0].Domain()).numPts());
    for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
    {
        if (!m_domain_covered[amrlev-1]) break;
        m_domain_covered[amrlev] = (m_grids[amrlev][0].numPts() ==
                                    compactify(m_geom[amrlev][0].Domain()).numPts());
    }

    Box aggbox;
    bool aggable = false;

    if (info.do_agglomeration)
    {
        if (m_domain_covered[0])
        {
            aggbox = m_geom[0][0].Domain();
            if (hasHiddenDimension()) {
                aggbox.makeSlab(hiddenDirection(), m_grids[0][0][0].smallEnd(hiddenDirection()));
            }
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
    int agg_lev = 0, con_lev = 0;

    AMREX_ALWAYS_ASSERT( ! (info.do_semicoarsening && info.hasHiddenDimension())
                         && info.semicoarsening_direction >= -1
                         && info.semicoarsening_direction < AMREX_SPACEDIM );

    if (info.do_agglomeration && aggable)
    {
        Box dbx = m_geom[0][0].Domain();
        Box bbx = aggbox;
        Real const nbxs = static_cast<Real>(m_grids[0][0].size());
        Real const threshold_npts = static_cast<Real>(AMREX_D_TERM(info.agg_grid_size,
                                                                  *info.agg_grid_size,
                                                                  *info.agg_grid_size));
        Vector<Box> domainboxes{dbx};
        Vector<Box> boundboxes{bbx};
        Vector<int> agg_flag{false};
        Vector<IntVect> accum_coarsen_ratio{IntVect(1)};
        int numsclevs = 0;

        for (int lev = 0; lev < info.max_coarsening_level; ++lev)
        {
            IntVect rr_level = mg_coarsen_ratio_v;
            bool const do_semicoarsening_level = info.do_semicoarsening
                && numsclevs < info.max_semicoarsening_level;
            if (do_semicoarsening_level
                && info.semicoarsening_direction != -1)
            {
                rr_level[info.semicoarsening_direction] = 1;
            }
            IntVect is_coarsenable;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                IntVect rr_dir(1);
                rr_dir[idim] = rr_level[idim];
                is_coarsenable[idim] = dbx.coarsenable(rr_dir, mg_domain_min_width_v)
                    && bbx.coarsenable(rr_dir, mg_box_min_width_v);
                if (!is_coarsenable[idim] && do_semicoarsening_level
                    && info.semicoarsening_direction == -1)
                {
                    is_coarsenable[idim] = true;
                    rr_level[idim] = 1;
                }
            }
            if (is_coarsenable != IntVect(1) || rr_level == IntVect(1)) {
                break;
            }
            if (do_semicoarsening_level && info.semicoarsening_direction == -1) {
                // make sure there is at most one direction that is not coarsened
                int n_ones = AMREX_D_TERM(  static_cast<int>(rr_level[0] == 1),
                                          + static_cast<int>(rr_level[1] == 1),
                                          + static_cast<int>(rr_level[2] == 1));
                if (n_ones > 1) { break; }
            }
            if (rr_level != mg_coarsen_ratio_v) {
                ++numsclevs;
            }

            accum_coarsen_ratio.push_back(accum_coarsen_ratio.back()*rr_level);
            domainboxes.push_back(dbx.coarsen(rr_level));
            boundboxes.push_back(bbx.coarsen(rr_level));
            bool to_agg = (bbx.d_numPts() / nbxs) < 0.999*threshold_npts;
            agg_flag.push_back(to_agg);
        }

        for (int lev = 1, nlevs = domainboxes.size(); lev < nlevs; ++lev) {
            if (!agged && !agg_flag[lev] &&
                a_grids[0].coarsenable(accum_coarsen_ratio[lev], mg_box_min_width_v))
            {
                m_grids[0].push_back(amrex::coarsen(a_grids[0], accum_coarsen_ratio[lev]));
                m_dmap[0].push_back(a_dmap[0]);
            } else {
                IntVect cr = domainboxes[lev-1].length() / domainboxes[lev].length();
                if (!m_grids[0].back().coarsenable(cr)) {
                    break; // average_down would fail if fine boxarray is not coarsenable.
                }
                m_grids[0].emplace_back(boundboxes[lev]);
                IntVect max_grid_size(info.agg_grid_size);
                if (info.do_semicoarsening && info.max_semicoarsening_level >= lev
                    && info.semicoarsening_direction != -1)
                {
                    IntVect blen = amrex::enclosedCells(boundboxes[lev]).size();
                    AMREX_D_TERM(int mgs_0 = (max_grid_size[0]+blen[0]-1) / blen[0];,
                                 int mgs_1 = (max_grid_size[1]+blen[1]-1) / blen[1];,
                                 int mgs_2 = (max_grid_size[2]+blen[2]-1) / blen[2]);
                    max_grid_size[info.semicoarsening_direction]
                        *= AMREX_D_TERM(mgs_0, *mgs_1, *mgs_2);
                }
                m_grids[0].back().maxSize(max_grid_size);
                m_dmap[0].push_back(DistributionMapping());
                if (!agged) {
                    agged = true;
                    agg_lev = lev;
                }
            }
            m_geom[0].emplace_back(domainboxes[lev],rb,coord,is_per);
        }
    }
    else
    {
        Real avg_npts = 0.0;
        if (info.do_consolidation) {
            avg_npts = static_cast<Real>(a_grids[0].d_numPts()) / static_cast<Real>(ParallelContext::NProcsSub());
            if (consolidation_threshold == -1) {
                consolidation_threshold = AMREX_D_TERM(info.con_grid_size,
                                                      *info.con_grid_size,
                                                      *info.con_grid_size);
            }
        }

        Box const& dom0 = a_geom[0].Domain();
        IntVect rr_vec(1);
        int numsclevs = 0;
        for (int lev = 0; lev < info.max_coarsening_level; ++lev)
        {
            IntVect rr_level = mg_coarsen_ratio_v;
            bool do_semicoarsening_level = info.do_semicoarsening
                && numsclevs < info.max_semicoarsening_level;
            if (do_semicoarsening_level
                && info.semicoarsening_direction != -1)
            {
                rr_level[info.semicoarsening_direction] = 1;
            }
            IntVect is_coarsenable;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                IntVect rr_dir(1);
                rr_dir[idim] = rr_vec[idim] * rr_level[idim];
                is_coarsenable[idim] = dom0.coarsenable(rr_dir, mg_domain_min_width_v)
                    && a_grids[0].coarsenable(rr_dir, mg_box_min_width_v);
                if (!is_coarsenable[idim] && do_semicoarsening_level
                    && info.semicoarsening_direction == -1)
                {
                    is_coarsenable[idim] = true;
                    rr_level[idim] = 1;
                }
            }
            if (is_coarsenable != IntVect(1) || rr_level == IntVect(1)) {
                break;
            }
            if (do_semicoarsening_level && info.semicoarsening_direction == -1) {
                // make sure there is at most one direction that is not coarsened
                int n_ones = AMREX_D_TERM(  static_cast<int>(rr_level[0] == 1),
                                          + static_cast<int>(rr_level[1] == 1),
                                          + static_cast<int>(rr_level[2] == 1));
                if (n_ones > 1) { break; }
            }
            if (rr_level != mg_coarsen_ratio_v) {
                ++numsclevs;
            }
            rr_vec *= rr_level;

            m_geom[0].emplace_back(amrex::coarsen(dom0, rr_vec), rb, coord, is_per);
            m_grids[0].push_back(amrex::coarsen(a_grids[0], rr_vec));

            if (info.do_consolidation)
            {
                if (avg_npts/(AMREX_D_TERM(rr_vec[0], *rr_vec[1], *rr_vec[2]))
                    < Real(0.999)*consolidation_threshold)
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
        }
    }

    m_num_mg_levels[0] = m_grids[0].size();

    for (int mglev = 0; mglev < m_num_mg_levels[0] - 1; mglev++){
        const Box& fine_domain = m_geom[0][mglev].Domain();
        const Box& crse_domain = m_geom[0][mglev+1].Domain();
        mg_coarsen_ratio_vec.push_back(fine_domain.length()/crse_domain.length());
    }

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        if (AMRRefRatio(amrlev) == 4 && mg_coarsen_ratio_vec.size() == 0) {
            mg_coarsen_ratio_vec.push_back(IntVect(2));
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

    if (agged || coned)
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
MLLinOp::make (Vector<Vector<Any> >& mf, IntVect const& ng) const
{
    mf.clear();
    mf.resize(m_num_amr_levels);
    for (int alev = 0; alev < m_num_amr_levels; ++alev)
    {
        mf[alev].resize(m_num_mg_levels[alev]);
        for (int mlev = 0; mlev < m_num_mg_levels[alev]; ++mlev)
        {
            mf[alev][mlev] = AnyMake(alev, mlev, ng);
        }
    }
}

void
MLLinOp::setDomainBC (const Array<BCType,AMREX_SPACEDIM>& a_lobc,
                      const Array<BCType,AMREX_SPACEDIM>& a_hibc) noexcept
{
    const int ncomp = getNComp();
    setDomainBC(Vector<Array<BCType,AMREX_SPACEDIM> >(ncomp,a_lobc),
                Vector<Array<BCType,AMREX_SPACEDIM> >(ncomp,a_hibc));
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
    m_lobc_orig = m_lobc;
    m_hibc_orig = m_hibc;
    for (int icomp = 0; icomp < ncomp; ++icomp) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (m_geom[0][0].isPeriodic(idim)) {
                AMREX_ALWAYS_ASSERT(m_lobc[icomp][idim] == BCType::Periodic &&
                                    m_hibc[icomp][idim] == BCType::Periodic);
            } else {
                AMREX_ALWAYS_ASSERT(m_lobc[icomp][idim] != BCType::Periodic &&
                                    m_hibc[icomp][idim] != BCType::Periodic);
            }

            if (m_lobc[icomp][idim] == LinOpBCType::inhomogNeumann ||
                m_lobc[icomp][idim] == LinOpBCType::Robin)
            {
                m_lobc[icomp][idim] = LinOpBCType::Neumann;
            }

            if (m_hibc[icomp][idim] == LinOpBCType::inhomogNeumann ||
                m_hibc[icomp][idim] == LinOpBCType::Robin)
            {
                m_hibc[icomp][idim] = LinOpBCType::Neumann;
            }
        }
    }

    if (hasHiddenDimension()) {
        const int hd = hiddenDirection();
        for (int n = 0; n < ncomp; ++n) {
            m_lobc[n][hd] = LinOpBCType::Neumann;
            m_hibc[n][hd] = LinOpBCType::Neumann;
        }
    }

    if (hasInhomogNeumannBC() && !supportInhomogNeumannBC()) {
        amrex::Abort("Inhomogeneous Neumann BC not supported");
    }
    if (hasRobinBC() && !supportRobinBC()) {
        amrex::Abort("Robin BC not supported");
    }
}

bool
MLLinOp::hasInhomogNeumannBC () const noexcept
{
    int ncomp = m_lobc_orig.size();
    for (int n = 0; n < ncomp; ++n) {
        for (int idim = 0; idim <AMREX_SPACEDIM; ++idim) {
            if (m_lobc_orig[n][idim] == BCType::inhomogNeumann ||
                m_hibc_orig[n][idim] == BCType::inhomogNeumann)
            {
                return true;
            }
        }
    }
    return false;
}

bool
MLLinOp::hasRobinBC () const noexcept
{
    int ncomp = m_lobc_orig.size();
    for (int n = 0; n < ncomp; ++n) {
        for (int idim = 0; idim <AMREX_SPACEDIM; ++idim) {
            if (m_lobc_orig[n][idim] == BCType::Robin ||
                m_hibc_orig[n][idim] == BCType::Robin)
            {
                return true;
            }
        }
    }
    return false;
}

void
MLLinOp::setDomainBCLoc (const Array<Real,AMREX_SPACEDIM>& lo_bcloc,
                         const Array<Real,AMREX_SPACEDIM>& hi_bcloc) noexcept
{
    m_domain_bloc_lo = lo_bcloc;
    m_domain_bloc_hi = hi_bcloc;
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
            m_raii_comm = std::make_unique<CommContainer>(newcomm);
        }

        MPI_Group_free(&defgrp);
        MPI_Group_free(&newgrp);
    }

    return newcomm;
#else
    amrex::ignore_unused(dm);
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
#ifdef AMREX_USE_MPI
    BL_PROFILE("MLLinOp::remapNeighborhoods()");

    if (flag_verbose_linop) {
        Print() << "Remapping ranks to neighborhoods ..." << std::endl;
    }

    for (int j = 1; j < dms.size(); ++j)
    {
        const Vector<int> & pmap = dms[j].ProcessorMap();
        std::set<int> g_ranks_set(pmap.begin(), pmap.end());
        int lev_rank_n = g_ranks_set.size();
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
#else
    amrex::ignore_unused(dms);
#endif
}

Box
MLLinOp::compactify (Box const& b) const noexcept
{
#if (AMREX_SPACEDIM == 3)
    if (info.hasHiddenDimension()) {
        const auto& lo = b.smallEnd();
        const auto& hi = b.bigEnd();
        if (info.hidden_direction == 0) {
            return Box(IntVect(lo[1],lo[2],0), IntVect(hi[1],hi[2],0), b.ixType());
        } else if (info.hidden_direction == 1) {
            return Box(IntVect(lo[0],lo[2],0), IntVect(hi[0],hi[2],0), b.ixType());
        } else {
            return Box(IntVect(lo[0],lo[1],0), IntVect(hi[0],hi[1],0), b.ixType());
        }
    } else
#endif
    {
        return b;
    }
}

void
MLLinOp::resizeMultiGrid (int new_size)
{
    if (new_size <= 0 || new_size >= m_num_mg_levels[0]) { return; }

    m_num_mg_levels[0] = new_size;

    m_geom[0].resize(new_size);
    m_grids[0].resize(new_size);
    m_dmap[0].resize(new_size);
    m_factory[0].resize(new_size);

    if (m_bottom_comm != m_default_comm) {
        m_bottom_comm = makeSubCommunicator(m_dmap[0].back());
    }
}

Any
MLLinOp::AnyMake (int amrlev, int mglev, IntVect const& ng) const
{
    return Any(MultiFab(amrex::convert(m_grids[amrlev][mglev], m_ixtype),
                        m_dmap[amrlev][mglev], getNComp(), ng, MFInfo(),
                        *m_factory[amrlev][mglev]));
}

Any
MLLinOp::AnyMakeCoarseMG (int amrlev, int mglev, IntVect const& ng) const
{
    BoxArray cba = m_grids[amrlev][mglev];
    IntVect ratio = (amrlev > 0) ? IntVect(2) : mg_coarsen_ratio_vec[mglev];
    cba.coarsen(ratio);
    cba.convert(m_ixtype);
    return Any(MultiFab(cba, m_dmap[amrlev][mglev], getNComp(), ng));
}

Any
MLLinOp::AnyMakeCoarseAmr (int famrlev, IntVect const& ng) const
{
    BoxArray cba = m_grids[famrlev][0];
    IntVect ratio(AMRRefRatio(famrlev-1));
    cba.coarsen(ratio);
    cba.convert(m_ixtype);
    return Any(MultiFab(cba, m_dmap[famrlev][0], getNComp(), ng));
}

Any
MLLinOp::AnyMakeAlias (Any const& a) const
{
    AMREX_ASSERT(a.is<MultiFab>());
    MultiFab const& mf = a.get<MultiFab>();
    return Any(MultiFab(mf, amrex::make_alias, 0, mf.nComp()));
}

IntVect
MLLinOp::AnyGrowVect (Any const& a) const
{
    AMREX_ASSERT(a.is<MultiFab>());
    MultiFab const& mf = a.get<MultiFab>();
    return mf.nGrowVect();
}

void
MLLinOp::AnySetToZero (Any& a) const
{
    AMREX_ASSERT(a.is<MultiFab>());
    MultiFab& mf = a.get<MultiFab>();
    mf.setVal(0._rt);
}

void
MLLinOp::AnySetBndryToZero (Any& a) const
{
    AMREX_ASSERT(a.is<MultiFab>());
    MultiFab& mf = a.get<MultiFab>();
    mf.setBndry(0._rt, 0, getNComp());
}

#ifdef AMREX_USE_EB
void
MLLinOp::AnySetCoveredToZero (Any& a) const
{
    AMREX_ASSERT(a.is<MultiFab>());
    auto& mf = a.get<MultiFab>();
    EB_set_covered(mf, 0, getNComp(), 0, 0._rt);
}
#endif

void
MLLinOp::AnyCopy (Any& dst, Any const& src, IntVect const& ng) const
{
    AMREX_ASSERT(dst.is<MultiFab>() && src.is<MultiFab>());
    MultiFab& dmf = dst.get<MultiFab>();
    MultiFab const& smf = src.get<MultiFab>();
    MultiFab::Copy(dmf, smf, 0, 0, getNComp(), ng);
}

void
MLLinOp::AnyAdd (Any& dst, Any const& src, IntVect const& ng) const
{
    AMREX_ASSERT(dst.is<MultiFab>() && src.is<MultiFab>());
    MultiFab& dmf = dst.get<MultiFab>();
    MultiFab const& smf = src.get<MultiFab>();
    MultiFab::Add(dmf, smf, 0, 0, getNComp(), ng);
}

void
MLLinOp::AnyAverageDownSolutionRHS (int camrlev, Any& a_crse_sol, Any& a_crse_rhs,
                                    const Any& a_fine_sol, const Any& a_fine_rhs)
{
    AMREX_ASSERT(a_crse_sol.is<MultiFab>() &&
                 a_crse_rhs.is<MultiFab>() &&
                 a_fine_sol.is<MultiFab>() &&
                 a_fine_rhs.is<MultiFab>());
    auto& crse_sol = a_crse_sol.get<MultiFab>();
    auto& crse_rhs = a_crse_rhs.get<MultiFab>();
    auto& fine_sol = a_fine_sol.get<MultiFab>();
    auto& fine_rhs = a_fine_rhs.get<MultiFab>();
    averageDownSolutionRHS(camrlev, crse_sol, crse_rhs, fine_sol, fine_rhs);
}

void
MLLinOp::AnyParallelCopy (Any& dst, Any const& src,
                          IntVect const& src_nghost, IntVect const& dst_nghost,
                          Periodicity const& period) const
{
    AMREX_ASSERT(dst.is<MultiFab>());
    MultiFab& dmf = dst.get<MultiFab>();
    MultiFab const& smf = src.get<MultiFab>();
    dmf.ParallelCopy(smf, 0, 0, getNComp(), src_nghost, dst_nghost, period);
}

Real
MLLinOp::AnyNormInf (Any& a) const
{
    AMREX_ASSERT(a.is<MultiFab>());
    return a.get<MultiFab>().norminf();
}

void
MLLinOp::AnySolutionResidual (int amrlev, Any& resid, Any& x, Any const& b,
                              Any const* crse_bcdata)
{
    AMREX_ASSERT(x.is<MultiFab>());
    solutionResidual(amrlev, resid.get<MultiFab>(), x.get<MultiFab>(), b.get<MultiFab>(),
                     (crse_bcdata) ? &(crse_bcdata->get<MultiFab>()) : nullptr);
}

void
MLLinOp::AnyCorrectionResidual (int amrlev, int mglev, Any& resid, Any& x, const Any& b,
                                BCMode bc_mode, const Any* crse_bcdata)
{
    AMREX_ASSERT(x.is<MultiFab>());
    correctionResidual(amrlev, mglev, resid.get<MultiFab>(), x.get<MultiFab>(),
                       b.get<MultiFab>(), bc_mode,
                       (crse_bcdata) ? &(crse_bcdata->get<MultiFab>()) : nullptr);
}

void
MLLinOp::AnyReflux (int clev, Any& res, const Any& crse_sol, const Any& crse_rhs,
                    Any& fine_res, Any& fine_sol, const Any& fine_rhs)
{
    AMREX_ASSERT(res.is<MultiFab>());
    reflux(clev,res.get<MultiFab>(), crse_sol.get<MultiFab>(), crse_rhs.get<MultiFab>(),
           fine_res.get<MultiFab>(), fine_sol.get<MultiFab>(), fine_rhs.get<MultiFab>());
}

Real
MLLinOp::MFNormInf (MultiFab const& mf, iMultiFab const* fine_mask, bool local) const
{
    const int ncomp = getNComp();
    Real norm = 0._rt;

    if (fine_mask == nullptr) {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            auto const& ma = mf.const_arrays();
            norm = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{},
                             mf, IntVect(0), ncomp,
                             [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n)
                                 -> GpuTuple<Real>
                             {
                                 return amrex::Math::abs(ma[box_no](i,j,k,n));
                             });
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:norm)
#endif
            for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                auto const& fab = mf.const_array(mfi);
                AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
                {
                    norm = std::max(norm, amrex::Math::abs(fab(i,j,k,n)));
                });
            }
        }
    } else {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            auto const& ma = mf.const_arrays();
            auto const& mask_ma = fine_mask->const_arrays();
            norm = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{},
                             mf, IntVect(0), ncomp,
                             [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n)
                                 -> GpuTuple<Real>
                             {
                                 if (mask_ma[box_no](i,j,k)) {
                                     return amrex::Math::abs(ma[box_no](i,j,k,n));
                                 } else {
                                     return Real(0.0);
                                 }
                             });
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:norm)
#endif
            for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                auto const& fab = mf.const_array(mfi);
                auto const& mask = fine_mask->const_array(mfi);
                AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
                {
                    if (mask(i,j,k)) {
                        norm = std::max(norm, amrex::Math::abs(fab(i,j,k,n)));
                    }
                });
            }
        }
    }

    if (!local) ParallelAllReduce::Max(norm, ParallelContext::CommunicatorSub());
    return norm;
}

void
MLLinOp::AnyAvgDownResMG (int clev, Any& cres, Any const& fres) const
{
    AMREX_ASSERT(cres.is<MultiFab>());
#ifdef AMREX_USE_EB
    amrex::EB_average_down
#else
    amrex::average_down
#endif
        (fres.get<MultiFab>(), cres.get<MultiFab>(), 0, getNComp(),
         mg_coarsen_ratio_vec[clev-1]);
}

void
MLLinOp::AnySmooth (int amrlev, int mglev, Any& sol, const Any& rhs,
                    bool skip_fillboundary) const
{
    AMREX_ASSERT(sol.is<MultiFab>() && rhs.is<MultiFab>());
    smooth(amrlev, mglev, sol.get<MultiFab>(), rhs.get<MultiFab>(), skip_fillboundary);
}

void
MLLinOp::AnyRestriction (int amrlev, int cmglev, Any& crse, Any& fine) const
{
    AMREX_ASSERT(crse.is<MultiFab>() && fine.is<MultiFab>());
    restriction(amrlev, cmglev, crse.get<MultiFab>(), fine.get<MultiFab>());
}

void
MLLinOp::AnyInterpolationMG (int amrlev, int fmglev, Any& fine, const Any& crse) const
{
    AMREX_ASSERT(crse.is<MultiFab>() && fine.is<MultiFab>());
    interpolation(amrlev, fmglev, fine.get<MultiFab>(), crse.get<MultiFab>());
}

void
MLLinOp::AnyInterpAssignMG (int amrlev, int fmglev, Any& fine, Any& crse) const
{
    AMREX_ASSERT(crse.is<MultiFab>() && fine.is<MultiFab>());
    interpAssign(amrlev, fmglev, fine.get<MultiFab>(), crse.get<MultiFab>());
}

bool
MLLinOp::isMFIterSafe (int amrlev, int mglev1, int mglev2) const
{
    return m_dmap[amrlev][mglev1] == m_dmap[amrlev][mglev2]
        && BoxArray::SameRefs(m_grids[amrlev][mglev1], m_grids[amrlev][mglev2]);
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
