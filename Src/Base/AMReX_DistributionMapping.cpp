
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_Utility.H>
#include <AMReX_Morton.H>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <map>
#include <vector>
#include <queue>
#include <algorithm>
#include <numeric>
#include <string>
#include <cstring>
#include <iomanip>

namespace {
int flag_verbose_mapper;
}

namespace amrex {

    bool initialized = false;
    //
    // Set default values for these in Initialize()!!!
    //
    int    verbose;
    int    sfc_threshold;
    Real   max_efficiency;
    int    node_size;

// We default to SFC.
DistributionMapping::Strategy DistributionMapping::m_Strategy = DistributionMapping::SFC;

DistributionMapping::PVMF DistributionMapping::m_BuildMap = nullptr;

const Vector<int>&
DistributionMapping::ProcessorMap () const noexcept
{
    return m_ref->m_pmap;
}

DistributionMapping::Strategy
DistributionMapping::strategy ()
{
    return DistributionMapping::m_Strategy;
}

void
DistributionMapping::strategy (DistributionMapping::Strategy how)
{
    DistributionMapping::m_Strategy = how;

    switch (how)
    {
    case ROUNDROBIN:
        m_BuildMap = &DistributionMapping::RoundRobinProcessorMap;
        break;
    case KNAPSACK:
        m_BuildMap = &DistributionMapping::KnapSackProcessorMap;
        break;
    case SFC:
        m_BuildMap = &DistributionMapping::SFCProcessorMap;
        break;
    case RRSFC:
        m_BuildMap = &DistributionMapping::RRSFCProcessorMap;
        break;
    default:
        amrex::Error("Bad DistributionMapping::Strategy");
    }
}

void
DistributionMapping::SFC_Threshold (int n)
{
    sfc_threshold = std::min(n,1);
}

int
DistributionMapping::SFC_Threshold ()
{
    return sfc_threshold;
}

bool
DistributionMapping::operator== (const DistributionMapping& rhs) const noexcept
{
    return m_ref == rhs.m_ref || m_ref->m_pmap == rhs.m_ref->m_pmap;
}

bool
DistributionMapping::operator!= (const DistributionMapping& rhs) const noexcept
{
    return !operator==(rhs);
}

void
DistributionMapping::Initialize ()
{
    if (initialized) return;
    //
    // Set defaults here!!!
    //
    verbose          = 0;
    sfc_threshold    = 0;
    max_efficiency   = 0.9_rt;
    node_size        = 0;
    flag_verbose_mapper = 0;

    ParmParse pp("DistributionMapping");

    pp.queryAdd("v"      ,             verbose);
    pp.queryAdd("verbose",             verbose);
    pp.queryAdd("efficiency",          max_efficiency);
    pp.queryAdd("sfc_threshold",       sfc_threshold);
    pp.queryAdd("node_size",           node_size);
    pp.queryAdd("verbose_mapper",      flag_verbose_mapper);

    std::string theStrategy;

    if (pp.query("strategy", theStrategy))
    {
        if (theStrategy == "ROUNDROBIN")
        {
            strategy(ROUNDROBIN);
        }
        else if (theStrategy == "KNAPSACK")
        {
            strategy(KNAPSACK);
        }
        else if (theStrategy == "SFC")
        {
            strategy(SFC);
        }
        else if (theStrategy == "RRSFC")
        {
            strategy(RRSFC);
        }
        else
        {
            std::string msg("Unknown strategy: ");
            msg += theStrategy;
            amrex::Warning(msg.c_str());
        }
    }
    else
    {
        strategy(m_Strategy);  // default
    }

    amrex::ExecOnFinalize(DistributionMapping::Finalize);

    initialized = true;
}

void
DistributionMapping::Finalize ()
{
    initialized = false;

    m_Strategy = SFC;

    DistributionMapping::m_BuildMap = nullptr;
}

void
DistributionMapping::Sort (std::vector<LIpair>& vec,
                           bool                 reverse)
{
    if (vec.size() > 1)
    {
        if (reverse) {
            std::stable_sort(vec.begin(), vec.end(), LIpairGT());
        }
        else {
            std::stable_sort(vec.begin(), vec.end(), LIpairLT());
        }
    }
}

void
DistributionMapping::LeastUsedCPUs (int         nprocs,
                                    Vector<int>& result)
{
    result.resize(nprocs);

#ifdef BL_USE_MPI
    BL_PROFILE("DistributionMapping::LeastUsedCPUs()");

    AMREX_ASSERT(nprocs <= ParallelContext::NProcsSub());

    Vector<Long> bytes(ParallelContext::NProcsSub());
    Long thisbyte = amrex::TotalBytesAllocatedInFabs()/1024;
    ParallelAllGather::AllGather(thisbyte, bytes.dataPtr(), ParallelContext::CommunicatorSub());

    std::vector<LIpair> LIpairV;

    LIpairV.reserve(nprocs);

    for (int i(0); i < nprocs; ++i)
    {
        LIpairV.emplace_back(bytes[i],i);
    }

    bytes.clear();

    Sort(LIpairV, false);

    for (int i(0); i < nprocs; ++i)
    {
        result[i] = LIpairV[i].second;
    }

    if (flag_verbose_mapper) {
        Print() << "LeastUsedCPUs:" << std::endl;
        for (const auto &p : LIpairV) {
            Print() << "  Rank " << p.second << " contains " << p.first << std::endl;
        }
    }
#else
    for (int i(0); i < nprocs; ++i)
    {
        result[i] = i;
    }
#endif
}

void
DistributionMapping::LeastUsedTeams (Vector<int>        & rteam,
                                     Vector<Vector<int> >& rworker,
                                     int                 nteams,
                                     int                 nworkers)
{
#ifdef BL_USE_MPI
    BL_PROFILE("DistributionMapping::LeastUsedTeams()");

    AMREX_ALWAYS_ASSERT(ParallelContext::CommunicatorSub() == ParallelDescriptor::Communicator());

    Vector<Long> bytes(ParallelContext::NProcsSub());
    Long thisbyte = amrex::TotalBytesAllocatedInFabs()/1024;
    ParallelAllGather::AllGather(thisbyte, bytes.dataPtr(), ParallelContext::CommunicatorSub());

    std::vector<LIpair> LIpairV;
    std::vector<LIpair> LIworker;

    LIpairV.reserve(nteams);
    LIworker.resize(nworkers);

    rteam.resize(nteams);
    rworker.resize(nteams);

    for (int i(0); i < nteams; ++i)
    {
        rworker[i].resize(nworkers);

        Long teambytes = 0;
        int offset = i*nworkers;
        for (int j = 0; j < nworkers; ++j)
        {
            int globalrank = offset+j;
            Long b = bytes[globalrank];
            teambytes += b;
            LIworker[j] = LIpair(b,j);
        }

        Sort(LIworker, false);

        for (int j = 0; j < nworkers; ++j)
        {
            rworker[i][j] = LIworker[j].second;
        }

        LIpairV.emplace_back(teambytes,i);
    }

    bytes.clear();

    Sort(LIpairV, false);

    for (int i(0); i < nteams; ++i)
    {
        rteam[i] = LIpairV[i].second;
    }
#else
    rteam.clear();
    rteam.push_back(0);
    rworker.clear();
    rworker.push_back(Vector<int>(1,0));
    amrex::ignore_unused(nteams,nworkers);
#endif
}

DistributionMapping::DistributionMapping () noexcept
    :
    m_ref(std::make_shared<Ref>())
{
}

DistributionMapping::DistributionMapping (const Vector<int>& pmap)
    :
    m_ref(std::make_shared<Ref>(pmap))
{
}

DistributionMapping::DistributionMapping (Vector<int>&& pmap) noexcept
    :
    m_ref(std::make_shared<Ref>(std::move(pmap)))
{
}

DistributionMapping::DistributionMapping (const BoxArray& boxes,
                                          int nprocs)
    :
    m_ref(std::make_shared<Ref>(boxes.size()))
{
    define(boxes,nprocs);
}

DistributionMapping::DistributionMapping (const DistributionMapping& d1,
                                          const DistributionMapping& d2)
    :
    m_ref(std::make_shared<Ref>())
{
    m_ref->m_pmap = d1.ProcessorMap();
    const auto& p2 = d2.ProcessorMap();
    m_ref->m_pmap.insert(m_ref->m_pmap.end(), p2.begin(), p2.end());
}

void
DistributionMapping::define (const BoxArray& boxes,
                             int nprocs)
{
    m_ref->clear();
    m_ref->m_pmap.resize(boxes.size());

    BL_ASSERT(m_BuildMap != nullptr);

    (this->*m_BuildMap)(boxes,nprocs);
}

void
DistributionMapping::define (const Vector<int>& pmap)
{
    m_ref->clear();
    m_ref->m_pmap = pmap;
}

void
DistributionMapping::define (Vector<int>&& pmap) noexcept
{
    m_ref->clear();
    m_ref->m_pmap = std::move(pmap);
}

void
DistributionMapping::RoundRobinDoIt (int                  nboxes,
                                     int                 /* nprocs */,
                                     std::vector<LIpair>* LIpairV,
                                     bool                 sort)
{
    if (flag_verbose_mapper) {
        Print() << "DM: RoundRobinDoIt called..." << std::endl;
    }

    int nprocs = ParallelContext::NProcsSub();

    // If team is not use, we are going to treat it as a special case in which
    // the number of teams is nprocs and the number of workers is 1.

    int nteams = nprocs;
    int nworkers = 1;
#if defined(BL_USE_TEAM)
    nteams = ParallelDescriptor::NTeams();
    nworkers = ParallelDescriptor::TeamSize();
#endif

    Vector<int> ord;
    Vector<Vector<int> > wrkerord;

    if (nteams == nprocs)  {
        if (sort) {
            LeastUsedCPUs(nprocs,ord);
        } else {
            ord.resize(nprocs);
            std::iota(ord.begin(), ord.end(), 0);
        }
        wrkerord.resize(nprocs);
        for (int i = 0; i < nprocs; ++i) {
            wrkerord[i].resize(1);
            wrkerord[i][0] = 0;
        }
    } else {
        if (sort) {
            LeastUsedTeams(ord,wrkerord,nteams,nworkers);
        } else {
            ord.resize(nteams);
            std::iota(ord.begin(), ord.end(), 0);
            wrkerord.resize(nteams);
            for (auto& v : wrkerord) {
                v.resize(nworkers);
                std::iota(v.begin(), v.end(), 0);
            }
        }
    }

    Vector<int> w(nteams,0);

    if (LIpairV)
    {
        BL_ASSERT(static_cast<int>(LIpairV->size()) == nboxes);

        for (int i = 0; i < nboxes; ++i)
        {
            int tid = ord[i%nteams];
            int wid = (w[tid]++) % nworkers;
            int rank = tid*nworkers + wrkerord[tid][wid];
            m_ref->m_pmap[(*LIpairV)[i].second] = ParallelContext::local_to_global_rank(rank);
            if (flag_verbose_mapper) {
                Print() << "  Mapping box " << (*LIpairV)[i].second << " of size "
                        << (*LIpairV)[i].first << " to rank " << rank << std::endl;
            }
        }
    }
    else
    {
        for (int i = 0; i < nboxes; ++i)
        {
            int tid = ord[i%nteams];
            int wid = (w[tid]++) % nworkers;
            int rank = tid*nworkers + wrkerord[tid][wid];
            m_ref->m_pmap[i] = ParallelContext::local_to_global_rank(rank);
            if (flag_verbose_mapper) {
                Print() << "  Mapping box " << i << " to rank " << rank << std::endl;
            }
        }
    }
}

void
DistributionMapping::RoundRobinProcessorMap (int nboxes, int nprocs, bool sort)
{
    BL_ASSERT(nboxes > 0);
    m_ref->clear();
    m_ref->m_pmap.resize(nboxes);

    RoundRobinDoIt(nboxes, nprocs, nullptr, sort);
}

void
DistributionMapping::RoundRobinProcessorMap (const BoxArray& boxes, int nprocs)
{
    BL_ASSERT( ! boxes.empty());
    BL_ASSERT(m_ref->m_pmap.size() == boxes.size());
    //
    // Create ordering of boxes from largest to smallest.
    // When we round-robin the boxes we want to go from largest
    // to smallest box, starting from the CPU having the least
    // amount of FAB data to the one having the most.  This "should"
    // help even out the FAB data distribution when running on large
    // numbers of CPUs, where the lower levels of the calculation are
    // using RoundRobin to lay out fewer than NProc boxes across
    // the CPUs.
    //
    std::vector<LIpair> LIpairV;

    const int N = static_cast<int>(boxes.size());

    LIpairV.reserve(N);

    for (int i = 0; i < N; ++i)
    {
        LIpairV.emplace_back(boxes[i].numPts(),i);
    }

    Sort(LIpairV, true);

    RoundRobinDoIt(static_cast<int>(boxes.size()), nprocs, &LIpairV);
}


void
DistributionMapping::RoundRobinProcessorMap (const std::vector<Long>& wgts,
                                             int nprocs, bool sort)
{
    BL_ASSERT( ! wgts.empty());

    m_ref->clear();
    m_ref->m_pmap.resize(wgts.size());

    //
    // Create ordering of boxes from "heaviest" to "lightest".
    // When we round-robin the boxes we want to go from heaviest
    // to lightest box, starting from the CPU having the least
    // amount of FAB data to the one having the most.  This "should"
    // help even out the FAB data distribution when running on large
    // numbers of CPUs, where the lower levels of the calculation are
    // using RoundRobin to lay out fewer than NProc boxes across
    // the CPUs.
    //
    std::vector<LIpair> LIpairV;

    const int N = static_cast<int>(wgts.size());

    LIpairV.reserve(N);

    for (int i = 0; i < N; ++i)
    {
        LIpairV.emplace_back(wgts[i],i);
    }

    Sort(LIpairV, true);

    RoundRobinDoIt(static_cast<int>(wgts.size()), nprocs, &LIpairV, sort);
}

class WeightedBox
{
    int  m_boxid;
    Long m_weight;
public:
    WeightedBox (int b, Long w) : m_boxid(b), m_weight(w) {}
    [[nodiscard]] Long weight () const { return m_weight; }
    [[nodiscard]] int boxid () const { return m_boxid;  }

    bool operator< (const WeightedBox& rhs) const
    {
        return weight() > rhs.weight();
    }
};

struct WeightedBoxList
{
    Vector<WeightedBox>* m_lb     = nullptr;
    Long                 m_weight = 0L;
    int                  m_rank   = -1;
    [[nodiscard]] Long weight () const
    {
        return m_weight;
    }
    [[nodiscard]] int rank () const { return m_rank; }
    void addWeight (Long dw) { m_weight += dw; }
    void erase (Vector<WeightedBox>::iterator& it)
    {
        m_weight -= it->weight();
        m_lb->erase(it);
    }
    void push_back (const WeightedBox& bx)
    {
        m_weight += bx.weight();
        m_lb->push_back(bx);
    }
    [[nodiscard]] int size () const { return static_cast<int>(m_lb->size()); }
    [[nodiscard]] Vector<WeightedBox>::const_iterator begin () const { return m_lb->begin(); }
    [[nodiscard]] Vector<WeightedBox>::iterator       begin ()       { return m_lb->begin(); } // NOLINT(readability-make-member-function-const)
    [[nodiscard]] Vector<WeightedBox>::const_iterator end   () const { return m_lb->end();   }
    [[nodiscard]] Vector<WeightedBox>::iterator       end   ()       { return m_lb->end();   } // NOLINT(readability-make-member-function-const)

    bool operator< (const WeightedBoxList& rhs) const
    {
        return weight() > rhs.weight();
    }
};

static
void
knapsack (const std::vector<Long>&         wgts,
          int                              nprocs,
          std::vector< std::vector<int> >& result,
          Real&                            efficiency,
          bool                             do_full_knapsack,
          int                              nmax)
{
    BL_PROFILE("knapsack()");

    //
    // Sort balls by size largest first.
    //
    result.resize(nprocs);

    Vector<WeightedBox> lb;
    lb.reserve(wgts.size());
    for (int i = 0, N = static_cast<int>(wgts.size()); i < N; ++i)
    {
        lb.emplace_back(i, wgts[i]);
    }
    std::sort(lb.begin(), lb.end());
    //
    // For each ball, starting with heaviest, assign ball to the lightest bin.
    //
    std::priority_queue<WeightedBoxList> wblq;
    Vector<std::unique_ptr<Vector<WeightedBox> > > raii_vwb(nprocs);
    for (int i  = 0; i < nprocs; ++i)
    {
        raii_vwb[i] = std::make_unique<Vector<WeightedBox> >();
        wblq.push(WeightedBoxList({raii_vwb[i].get(),Long(0),-1}));
    }
    Vector<WeightedBoxList> wblv;
    wblv.reserve(nprocs);
    for (unsigned int i = 0, N = wgts.size(); i < N; ++i)
    {
        if (!wblq.empty()) {
            WeightedBoxList wbl = wblq.top();
            wblq.pop();
            wbl.push_back(lb[i]);
            if (wbl.size() < nmax) {
                wblq.push(wbl);
            } else {
                wblv.push_back(wbl);
            }
        } else {
            int ip = static_cast<int>(i) % nprocs;
            wblv[ip].push_back(lb[i]);
        }
    }

    Real max_weight = 0;
    Real sum_weight = 0;
    for (auto const& wbl : wblv)
    {
        Real wgt = static_cast<Real>(wbl.weight());
        sum_weight += wgt;
        max_weight = std::max(wgt, max_weight);
    }

    while (!wblq.empty())
    {
        WeightedBoxList wbl = wblq.top();
        wblq.pop();
        if (wbl.size() > 0) {
            Real wgt = static_cast<Real>(wbl.weight());
            sum_weight += wgt;
            max_weight = std::max(wgt, max_weight);
            wblv.push_back(wbl);
        }
    }

    efficiency = sum_weight/(static_cast<Real>(nprocs)*max_weight);

    std::sort(wblv.begin(), wblv.end());

    if (efficiency < max_efficiency && do_full_knapsack
        && wblv.size() > 1 && wblv.begin()->size() > 1)
    {
        BL_PROFILE_VAR("knapsack()swap", swap);
top: ;

        if (efficiency < max_efficiency && wblv.begin()->size() > 1)
        {
            auto bl_top = wblv.begin();
            auto bl_bottom = wblv.end()-1;
            Long w_top = bl_top->weight();
            Long w_bottom = bl_bottom->weight();
            for (auto ball_1 = bl_top->begin(); ball_1 != bl_top->end(); ++ball_1)
            {
                for (auto ball_2 = bl_bottom->begin(); ball_2 != bl_bottom->end(); ++ball_2)
                {
                    // should we swap ball 1 and ball 2?
                    Long dw = ball_1->weight() - ball_2->weight();
                    Long w_top_new    = w_top    - dw;
                    Long w_bottom_new = w_bottom + dw;
                    if (w_top_new < w_top && w_bottom_new < w_top)
                    {
                        std::swap(*ball_1, *ball_2);
                        bl_top->addWeight(-dw);
                        bl_bottom->addWeight(dw);

                        if (bl_top+1 == bl_bottom)  // they are next to each other
                        {
                            if (*bl_bottom < *bl_top) {
                                std::swap(*bl_top, *bl_bottom);
                            }
                        }
                        else
                        {
                            // bubble up
                            auto it = std::lower_bound(bl_top+1, bl_bottom, *bl_bottom);
                            std::rotate(it, bl_bottom, bl_bottom+1);

                            // sink down
                            it = std::lower_bound(bl_top+1, bl_bottom+1, *bl_top);
                            std::rotate(bl_top, bl_top+1, it);
                        }

                        max_weight = static_cast<Real>(bl_top->weight());
                        efficiency = sum_weight / (static_cast<Real>(nprocs)*max_weight);
                        goto top;
                    }
                }
            }
        }

        BL_ASSERT(std::is_sorted(wblv.begin(), wblv.end()));
    }

    for (int i = 0, N = static_cast<int>(wblv.size()); i < N; ++i)
    {
        const WeightedBoxList& wbl = wblv[i];

        result[i].reserve(wbl.size());
        for (auto const& wb : wbl)
        {
            result[i].push_back(wb.boxid());
        }
    }
}

void
DistributionMapping::KnapSackDoIt (const std::vector<Long>& wgts,
                                   int                    /*  nprocs */,
                                   Real&                    efficiency,
                                   bool                     do_full_knapsack,
                                   int                      nmax,
                                   bool                     sort)
{
    if (flag_verbose_mapper) {
        Print() << "DM: KnapSackDoIt called..." << std::endl;
    }

    BL_PROFILE("DistributionMapping::KnapSackDoIt()");

    int nprocs = ParallelContext::NProcsSub();

    // If team is not use, we are going to treat it as a special case in which
    // the number of teams is nprocs and the number of workers is 1.

    int nteams = nprocs;
    int nworkers = 1;
#if defined(BL_USE_TEAM)
    nteams = ParallelDescriptor::NTeams();
    nworkers = ParallelDescriptor::TeamSize();
#endif

    std::vector< std::vector<int> > vec;

    efficiency = 0;

    knapsack(wgts,nteams,vec,efficiency,do_full_knapsack,nmax);

    if (flag_verbose_mapper) {
        for (int i = 0, ni = static_cast<int>(vec.size()); i < ni; ++i) {
            Print() << "  Bucket " << i << " contains boxes:" << std::endl;
            for (int x : vec[i]) {
                Print() << "    " << x << std::endl;
            }
        }
    }

    BL_ASSERT(static_cast<int>(vec.size()) == nteams);

    std::vector<LIpair> LIpairV;

    LIpairV.reserve(nteams);

    for (int i = 0; i < nteams; ++i)
    {
        Long wgt = 0;
        for (int j : vec[i]) {
            wgt += wgts[j];
        }

        LIpairV.emplace_back(wgt,i);
    }

    if (sort) {Sort(LIpairV, true);}

    if (flag_verbose_mapper) {
        for (const auto &p : LIpairV) {
            Print() << "  Bucket " << p.second << " total weight: " << p.first << std::endl;
        }
    }

    Vector<int> ord;
    Vector<Vector<int> > wrkerord;

    if (nteams == nprocs) {
        if (sort) {
            LeastUsedCPUs(nprocs,ord);
        } else {
            ord.resize(nprocs);
            std::iota(ord.begin(), ord.end(), 0);
        }
    } else {
        if (sort) {
            LeastUsedTeams(ord,wrkerord,nteams,nworkers);
        } else {
            ord.resize(nteams);
            std::iota(ord.begin(), ord.end(), 0);
            wrkerord.resize(nteams);
            for (auto& v : wrkerord) {
                v.resize(nworkers);
                std::iota(v.begin(), v.end(), 0);
            }
        }
    }

    for (int i = 0; i < nteams; ++i)
    {
        const int idx = LIpairV[i].second;
        const int tid = ord[i];

        const std::vector<int>& vi = vec[idx];
        const int N = static_cast<int>(vi.size());

        if (flag_verbose_mapper) {
            Print() << "  Mapping bucket " << idx << " to rank " << tid << std::endl;
        }

        if (nteams == nprocs) {
            for (int j = 0; j < N; ++j)
            {
                m_ref->m_pmap[vi[j]] = ParallelContext::local_to_global_rank(tid);
            }
        } else {
#ifdef BL_USE_TEAM
            int leadrank = tid * nworkers;
            for (int w = 0; w < nworkers; ++w)
            {
                ParallelDescriptor::team_for(0, N, w, [&] (int j) {
                        m_ref->m_pmap[vi[j]] = leadrank + wrkerord[i][w];
                    });
            }
#endif
        }
    }

    if (verbose)
    {
        amrex::Print() << "KNAPSACK efficiency: " << efficiency << '\n';
    }

}

void
DistributionMapping::KnapSackProcessorMap (const std::vector<Long>& wgts,
                                           int                      nprocs,
                                           Real*                    efficiency,
                                           bool                     do_full_knapsack,
                                           int                      nmax,
                                           bool                     sort)
{
    BL_ASSERT( ! wgts.empty());

    m_ref->clear();
    m_ref->m_pmap.resize(wgts.size());

    if (static_cast<int>(wgts.size()) <= nprocs || nprocs < 2)
    {
        RoundRobinProcessorMap(static_cast<int>(wgts.size()),nprocs, sort);

        if (efficiency) *efficiency = 1;
    }
    else
    {
        Real eff = 0;
        KnapSackDoIt(wgts, nprocs, eff, do_full_knapsack, nmax, sort);
        if (efficiency) *efficiency = eff;
    }
}

void
DistributionMapping::KnapSackProcessorMap (const DistributionMapping& olddm,
                                           const std::vector<Long>& wgts,
                                           Real keep_ratio, Real& old_efficiency,
                                           Real& new_efficiency, int nmax)
{
    BL_PROFILE("KnapSack(keep)");

    const int nprocs = ParallelDescriptor::NProcs();
    BL_ASSERT( ! wgts.empty());

    m_ref->clear();
    m_ref->m_pmap.resize(wgts.size(), -1);

    ComputeDistributionMappingEfficiency(olddm, wgts, &old_efficiency);

    if (static_cast<int>(wgts.size()) <= nprocs || nprocs < 2)
    {
        RoundRobinProcessorMap(static_cast<int>(wgts.size()),nprocs, false);
        new_efficiency = Real(1);
        return;
    }
    else
    {
        Vector<Vector<WeightedBox>> keep_balls(nprocs);
        Long sum_weight = 0;
        for (int i = 0, N = static_cast<int>(wgts.size()); i < N; ++i) {
            keep_balls[olddm[i]].emplace_back(i, wgts[i]);
            sum_weight += wgts[i];
        }

        Real avg_weight = static_cast<Real>(sum_weight) / static_cast<Real>(nprocs);
        Real keep_weight = avg_weight * static_cast<Real>(keep_ratio);

        Vector<WeightedBox> lb;
        Vector<Long> base_weight(nprocs);
        for (int iproc = 0; iproc < nprocs; ++iproc) {
            std::sort(keep_balls[iproc].begin(), keep_balls[iproc].end());
            Long w = 0;
            auto& kb = keep_balls[iproc];
            int i = 0;
            for (int N = static_cast<int>(kb.size()); i < N; ++i) {
                auto wi = kb[i].weight();
                if (static_cast<Real>(w+wi) > keep_weight) {
                    break;
                } else {
                    w += wi;
                }
            }
            if (i < static_cast<int>(kb.size())) {
                lb.insert(lb.end(), kb.begin()+i, kb.end());
                kb.erase (          kb.begin()+i, kb.end());
            }
            base_weight[iproc] = w;
        }

        if (lb.empty()) {
            *this = olddm;
            new_efficiency = old_efficiency;
            return;
        } else {
            std::sort(lb.begin(), lb.end());

            // Vector<Vector<WeightedBox>> keep_balls : we keep
            // Vector<Long>                base_weight: weight of balls kept
            // Vector<WeightedBox>         lb         : we still need to assign them

            std::priority_queue<WeightedBoxList> wblq;
            Vector<std::unique_ptr<Vector<WeightedBox>>> raii_vwb(nprocs);
            for (int iproc = 0; iproc < nprocs; ++iproc) {
                raii_vwb[iproc] = std::make_unique<Vector<WeightedBox>>();
                wblq.push(WeightedBoxList({raii_vwb[iproc].get(), base_weight[iproc],
                                           iproc}));
            }
            Vector<WeightedBoxList> wblv;
            wblv.reserve(nprocs);
            for (int i = 0, N = static_cast<int>(lb.size()); i < N; ++i) {
                if (!wblq.empty()) {
                    WeightedBoxList wbl = wblq.top();
                    wblq.pop();
                    wbl.push_back(lb[i]);
                    if (wbl.size() + static_cast<int>(keep_balls[wbl.rank()].size())
                        < nmax) {
                        wblq.push(wbl);
                    } else {
                        wblv.push_back(wbl);
                    }
                } else {
                    int ip = static_cast<int>(i) % nprocs;
                    wblv[ip].push_back(lb[i]);
                }
            }

            Real max_weight = Real(0);
            for (auto const& wbl : wblv) {
                auto wgt = static_cast<Real>(wbl.weight());
                max_weight = std::max(wgt, max_weight);
            }

            while (!wblq.empty()) {
                WeightedBoxList wbl = wblq.top();
                wblq.pop();
                if (wbl.size() > 0) {
                    auto wgt = static_cast<Real>(wbl.weight());
                    max_weight = std::max(wgt, max_weight);
                    wblv.push_back(wbl);
                }
            }

            new_efficiency = avg_weight / max_weight;

            if (new_efficiency < max_efficiency && wblv.size() > 1) {
                BL_PROFILE_VAR("knapsack()swap", swap);

                std::sort(wblv.begin(), wblv.end());

            top: ;
                if (new_efficiency < max_efficiency && wblv.begin()->size() > 1) {
                    auto bl_top = wblv.begin();
                    auto bl_bottom = wblv.end()-1;
                    Long w_top = bl_top->weight();
                    Long w_bottom = bl_bottom->weight();
                    for (auto& ball_1 : *bl_top) {
                        for (auto& ball_2 : *bl_bottom) {
                            // should we swap ball 1 and ball 2?
                            Long dw = ball_1.weight() - ball_2.weight();
                            Long w_top_new    = w_top    - dw;
                            Long w_bottom_new = w_bottom + dw;
                            if (w_top_new < w_top && w_bottom_new < w_top) {
                                std::swap(ball_1, ball_2);
                                bl_top->addWeight(-dw);
                                bl_bottom->addWeight(dw);

                                if (bl_top+1 == bl_bottom) {
                                    // they are next to each other
                                    if (*bl_bottom < *bl_top) {
                                        std::swap(*bl_top, *bl_bottom);
                                    }
                                } else {
                                    // bubble up
                                    auto it = std::lower_bound(bl_top+1, bl_bottom,
                                                               *bl_bottom);
                                    std::rotate(it, bl_bottom, bl_bottom+1);
                                    // sink down
                                    it = std::lower_bound(bl_top+1, bl_bottom+1,
                                                          *bl_top);
                                    std::rotate(bl_top, bl_top+1, it);
                                }

                                max_weight = static_cast<Real>(bl_top->weight());
                                new_efficiency = avg_weight / max_weight;
                                goto top;
                            }
                        }
                    }
                }
            }

            for (auto const& wbl : wblv) {
                for (auto const& wb : wbl) {
                    m_ref->m_pmap[wb.boxid()] = wbl.rank();
                }
            }

            for (int iproc = 0; iproc < nprocs; ++iproc) {
                for (auto const& wb : keep_balls[iproc]) {
                    m_ref->m_pmap[wb.boxid()] = iproc;
                }
            }

            AMREX_ASSERT(std::none_of(m_ref->m_pmap.cbegin(), m_ref->m_pmap.cend(),
                                      [] (int i) { return i < 0; }));
        }
    }
}

void
DistributionMapping::KnapSackProcessorMap (const BoxArray& boxes,
                                           int             nprocs)
{
    BL_ASSERT( ! boxes.empty());

    m_ref->m_pmap.resize(boxes.size());

    if (boxes.size() <= nprocs || nprocs < 2)
    {
        RoundRobinProcessorMap(boxes,nprocs);
    }
    else
    {
        std::vector<Long> wgts(boxes.size());

        for (int i = 0, N = static_cast<int>(boxes.size()); i < N; ++i) {
            wgts[i] = boxes[i].numPts();
        }

        Real effi = 0;
        bool do_full_knapsack = true;
        KnapSackDoIt(wgts, nprocs, effi, do_full_knapsack);
    }
}

namespace
{
    struct SFCToken
    {
        class Compare
        {
        public:
            AMREX_FORCE_INLINE
            bool operator () (const SFCToken& lhs,
                              const SFCToken& rhs) const;
        };
        int m_box;
        Array<uint32_t,AMREX_SPACEDIM> m_morton;
    };

    AMREX_FORCE_INLINE
    bool
    SFCToken::Compare::operator () (const SFCToken& lhs,
                                    const SFCToken& rhs) const
    {
    #if (AMREX_SPACEDIM == 1)
            return lhs.m_morton[0] < rhs.m_morton[0];
    #elif (AMREX_SPACEDIM == 2)
            return (lhs.m_morton[1] <  rhs.m_morton[1]) ||
                  ((lhs.m_morton[1] == rhs.m_morton[1]) &&
                   (lhs.m_morton[0] <  rhs.m_morton[0]));
    #else
            return (lhs.m_morton[2] <  rhs.m_morton[2]) ||
                  ((lhs.m_morton[2] == rhs.m_morton[2]) &&
                  ((lhs.m_morton[1] <  rhs.m_morton[1]) ||
                  ((lhs.m_morton[1] == rhs.m_morton[1]) &&
                   (lhs.m_morton[0] <  rhs.m_morton[0]))));
    #endif
    }
}

namespace {

    AMREX_FORCE_INLINE
    SFCToken makeSFCToken (int box_index, IntVect const& iv)
    {
        SFCToken token;
        token.m_box = box_index;

#if (AMREX_SPACEDIM == 3)

        constexpr int imin = -(1 << 29);
        AMREX_ASSERT_WITH_MESSAGE(AMREX_D_TERM(iv[0] >= imin && iv[0] < -imin,
                                            && iv[1] >= imin && iv[1] < -imin,
                                            && iv[2] >= imin && iv[2] < -imin),
                                  "SFCToken: index out of range");
        uint32_t x = iv[0] - imin;
        uint32_t y = iv[1] - imin;
        uint32_t z = iv[2] - imin;
        // extract lowest 10 bits and make space for interleaving
        token.m_morton[0] = Morton::makeSpace(x & 0x3FF)
                         | (Morton::makeSpace(y & 0x3FF) << 1)
                         | (Morton::makeSpace(z & 0x3FF) << 2);
        x = x >> 10;
        y = y >> 10;
        z = z >> 10;
        token.m_morton[1] = Morton::makeSpace(x & 0x3FF)
                         | (Morton::makeSpace(y & 0x3FF) << 1)
                         | (Morton::makeSpace(z & 0x3FF) << 2);
        x = x >> 10;
        y = y >> 10;
        z = z >> 10;
        token.m_morton[2] = Morton::makeSpace(x & 0x3FF)
                         | (Morton::makeSpace(y & 0x3FF) << 1)
                         | (Morton::makeSpace(z & 0x3FF) << 2);

#elif (AMREX_SPACEDIM == 2)

        constexpr uint32_t offset = 1U << 31;
        static_assert(static_cast<uint32_t>(std::numeric_limits<int>::max())+1 == offset,
                      "INT_MAX != (1<<31)-1");
        uint32_t x = (iv[0] >= 0) ? static_cast<uint32_t>(iv[0]) + offset
            : static_cast<uint32_t>(iv[0]-std::numeric_limits<int>::lowest());
        uint32_t y = (iv[1] >= 0) ? static_cast<uint32_t>(iv[1]) + offset
            : static_cast<uint32_t>(iv[1]-std::numeric_limits<int>::lowest());
        // extract lowest 16 bits and make sapce for interleaving
        token.m_morton[0] = Morton::makeSpace(x & 0xFFFF)
                         | (Morton::makeSpace(y & 0xFFFF) << 1);
        x = x >> 16;
        y = y >> 16;
        token.m_morton[1] = Morton::makeSpace(x) | (Morton::makeSpace(y) << 1);

#elif (AMREX_SPACEDIM == 1)

        constexpr uint32_t offset = 1U << 31;
        static_assert(static_cast<uint32_t>(std::numeric_limits<int>::max())+1 == offset,
                      "INT_MAX != (1<<31)-1");
        token.m_morton[0] = (iv[0] >= 0) ? static_cast<uint32_t>(iv[0]) + offset
            : static_cast<uint32_t>(iv[0]-std::numeric_limits<int>::lowest());

#else
        static_assert(false,"AMREX_SPACEDIM != 1, 2 or 3");
#endif

        return token;
    }
}

static
void
Distribute (const std::vector<SFCToken>&     tokens,
            const std::vector<Long>&         wgts,
            int                              nprocs,
            Real                             volpercpu,
            std::vector< std::vector<int> >& v)

{
    BL_PROFILE("DistributionMapping::Distribute()");

    if (flag_verbose_mapper) {
        Print() << "Distribute:" << std::endl;
        Print() << "  volpercpu: " << volpercpu << std::endl;
        Print() << "  Sorted SFC Tokens:" << std::endl;
        int idx = 0;
        for (const auto &t : tokens) {
            Print() << "    " << idx++ << ": "
                    << t.m_box << ": "
                    << t.m_morton << std::endl;
        }
    }

    BL_ASSERT(static_cast<int>(v.size()) == nprocs);

    int  K        = 0;
    Real totalvol = 0;

    for (int i = 0; i < nprocs; ++i)
    {
        int  cnt = 0;
        Real vol = 0;

        for ( int TSZ = static_cast<int>(tokens.size());
              K < TSZ && (i == (nprocs-1) || (vol < volpercpu));
              ++K)
        {
            vol += static_cast<Real>(wgts[tokens[K].m_box]);
            ++cnt;

            v[i].push_back(tokens[K].m_box);
        }

        totalvol += vol;

        if ((totalvol/static_cast<Real>(i+1)) > volpercpu &&  // Too much for this bin.
            cnt > 1                      &&  // More than one box in this bin.
            i < nprocs-1)                    // Not the last bin, which has to take all.
        {
            --K;
            v[i].pop_back();
            totalvol -= static_cast<Real>(wgts[tokens[K].m_box]);
        }
    }

    if (flag_verbose_mapper) {
        Print() << "Distributed SFC Tokens:" << std::endl;
        int idx = 0;
        for (int i = 0; i < nprocs; ++i) {
            Print() << "  Rank/Team " << i << ":" << std::endl;
            Real rank_vol = 0;
            for (const auto &box : v[i]) {
                amrex::ignore_unused(box);
                const auto &t = tokens[idx];
                BL_ASSERT(box == t.m_box);
                Print() << "    " << idx << ": "
                        << t.m_box << ": "
                        << t.m_morton << std::endl;
                rank_vol += static_cast<Real>(wgts[t.m_box]);
                idx++;
            }
            Print() << "    Total Rank Vol: " << rank_vol << std::endl;
        }
    }

#ifdef AMREX_DEBUG
    std::size_t cnt = 0;
    for (int i = 0; i < nprocs; ++i) {
        cnt += v[i].size();
    }
    BL_ASSERT(cnt == tokens.size());
#endif
}

void
DistributionMapping::SFCProcessorMapDoIt (const BoxArray&          boxes,
                                          const std::vector<Long>& wgts,
                                          int                   /*   nprocs */,
                                          bool                     sort,
                                          Real*                    eff)
{
    if (flag_verbose_mapper) {
        Print() << "DM: SFCProcessorMapDoIt called..." << std::endl;
    }

    BL_PROFILE("DistributionMapping::SFCProcessorMapDoIt()");

    int nprocs = ParallelContext::NProcsSub();

    int nteams = nprocs;
    int nworkers = 1;
#if defined(BL_USE_TEAM)
    nteams = ParallelDescriptor::NTeams();
    nworkers = ParallelDescriptor::TeamSize();
#else
    if (node_size > 0) {
        nteams = nprocs/node_size;
        nworkers = node_size;
        if (nworkers*nteams != nprocs) {
            nteams = nprocs;
            nworkers = 1;
        }
    }
#endif

    if (flag_verbose_mapper) {
        Print() << "  (nprocs, nteams, nworkers) = ("
                << nprocs << ", " << nteams << ", " << nworkers << ")\n";
    }

    const int N = static_cast<int>(boxes.size());
    std::vector<SFCToken> tokens;
    tokens.reserve(N);
    for (int i = 0; i < N; ++i)
    {
        const Box& bx = boxes[i];
        tokens.push_back(makeSFCToken(i, bx.smallEnd()));
    }
    //
    // Put'm in Morton space filling curve order.
    //
    std::sort(tokens.begin(), tokens.end(), SFCToken::Compare());
    //
    // Split'm up as equitably as possible per team.
    //
    Real volperteam = 0;
    for (Long wt : wgts) {
        volperteam += static_cast<Real>(wt);
    }
    volperteam /= static_cast<Real>(nteams);

    std::vector< std::vector<int> > vec(nteams);

    Distribute(tokens,wgts,nteams,volperteam,vec);

    // vec has a size of nteams and vec[] holds a vector of box ids.

    tokens.clear();

    std::vector<LIpair> LIpairV;

    LIpairV.reserve(nteams);

    for (int i = 0; i < nteams; ++i)
    {
        Long wgt = 0;
        const std::vector<int>& vi = vec[i];
        for (int j : vi) {
            wgt += wgts[j];
        }

        LIpairV.emplace_back(wgt,i);
    }

    if (sort) Sort(LIpairV, true);

    if (flag_verbose_mapper) {
        for (const auto &p : LIpairV) {
            Print() << "  Bucket " << p.second << " contains " << p.first << std::endl;
        }
    }

    // LIpairV has a size of nteams and LIpairV[] is pair whose first is weight
    // and second is an index into vec.  LIpairV is sorted by weight such that
    // LIpairV is the heaviest.

    Vector<int> ord;
    Vector<Vector<int> > wrkerord;

    if (nteams == nprocs) {
        if (sort) {
            LeastUsedCPUs(nprocs,ord);
        } else {
            ord.resize(nprocs);
            std::iota(ord.begin(), ord.end(), 0);
        }
    } else {
        if (sort) {
            LeastUsedTeams(ord,wrkerord,nteams,nworkers);
        } else {
            ord.resize(nteams);
            std::iota(ord.begin(), ord.end(), 0);
            wrkerord.resize(nteams);
            for (auto& v : wrkerord) {
                v.resize(nworkers);
                std::iota(v.begin(), v.end(), 0);
            }
        }
    }

    // ord is a vector of process (or team) ids, sorted from least used to more heavily used.
    // wrkerord is a vector of sorted worker ids.

    for (int i = 0; i < nteams; ++i)
    {
        const int tid  = ord[i];                  // tid is team id
        const int ivec = LIpairV[i].second;       // index into vec
        const std::vector<int>& vi = vec[ivec];   // this vector contains boxes assigned to this team
        const int Nbx = static_cast<int>(vi.size());// # of boxes assigned to this team

        if (flag_verbose_mapper) {
            Print() << "Mapping bucket " << LIpairV[i].second << " to rank " << ord[i] << std::endl;
        }

        if (nteams == nprocs) { // In this case, team id is process id.
            for (int j = 0; j < Nbx; ++j)
            {
                m_ref->m_pmap[vi[j]] = ParallelContext::local_to_global_rank(tid);
            }
        }
        else   // We would like to do knapsack within the team workers
        {
            std::vector<Long> local_wgts;
            local_wgts.reserve(Nbx);
            for (int j = 0; j < Nbx; ++j) {
                local_wgts.push_back(wgts[vi[j]]);
            }

            std::vector<std::vector<int> > kpres;
            Real kpeff;
            knapsack(local_wgts, nworkers, kpres, kpeff, true, N);

            // kpres has a size of nworkers. kpres[] contains a vector of indices into vi.

            // sort the knapsacked chunks
            std::vector<LIpair> ww;
            for (int w = 0; w < nworkers; ++w) {
                Long wgt = 0;
                for (int it : kpres[w]) {
                    wgt += local_wgts[it];
                }
                ww.emplace_back(wgt,w);
            }
            Sort(ww,true);

            // ww is a sorted vector of pair whose first is the weight and second is a index
            // into kpres.

            const Vector<int>& sorted_workers = wrkerord[i];

            const int leadrank = tid * nworkers;

            for (int w = 0; w < nworkers; ++w)
            {
                const int cpu = leadrank + sorted_workers[w];
                int ikp = ww[w].second;
                const std::vector<int>& js = kpres[ikp];
                for (int it : js) {
                    m_ref->m_pmap[vi[it]] = cpu;
                }
            }
        }
    }

    if (eff || verbose)
    {
        Long sum_wgt = 0, max_wgt = 0;
        for (int i = 0; i < nteams; ++i)
        {
            const Long W = LIpairV[i].first;
            if (W > max_wgt) max_wgt = W;
            sum_wgt += W;
        }
        Real efficiency = static_cast<Real>(sum_wgt)/static_cast<Real>(nteams*max_wgt);
        if (eff) *eff = efficiency;

        if (verbose)
        {
            amrex::Print() << "SFC efficiency: " << efficiency << '\n';
        }
    }
}

void
DistributionMapping::SFCProcessorMap (const BoxArray& boxes, int nprocs)
{
    BL_ASSERT( ! boxes.empty());

    m_ref->clear();
    m_ref->m_pmap.resize(boxes.size());

    if (boxes.size() < sfc_threshold*nprocs)
    {
        KnapSackProcessorMap(boxes,nprocs);
    }
    else
    {
        std::vector<Long> wgts;

        wgts.reserve(boxes.size());

        for (int i = 0, N = static_cast<int>(boxes.size()); i < N; ++i)
        {
            wgts.push_back(boxes[i].numPts());
        }

        SFCProcessorMapDoIt(boxes,wgts,nprocs);
    }
}

void
DistributionMapping::SFCProcessorMap (const BoxArray&          boxes,
                                      const std::vector<Long>& wgts,
                                      int                      nprocs,
                                      bool                     sort)
{
    BL_ASSERT( ! boxes.empty());
    BL_ASSERT(boxes.size() == static_cast<int>(wgts.size()));

    m_ref->clear();
    m_ref->m_pmap.resize(wgts.size());

    if (boxes.size() < sfc_threshold*nprocs)
    {
        KnapSackProcessorMap(wgts,nprocs);
    }
    else
    {
        SFCProcessorMapDoIt(boxes,wgts,nprocs,sort);
    }
}

void
DistributionMapping::SFCProcessorMap (const BoxArray&          boxes,
                                      const std::vector<Long>& wgts,
                                      int                      nprocs,
                                      Real&                    eff,
                                      bool                     sort)
{
    BL_ASSERT( ! boxes.empty());
    BL_ASSERT(boxes.size() == static_cast<int>(wgts.size()));

    m_ref->clear();
    m_ref->m_pmap.resize(wgts.size());

    if (boxes.size() < sfc_threshold*nprocs)
    {
        KnapSackProcessorMap(wgts,nprocs,&eff);
    }
    else
    {
        SFCProcessorMapDoIt(boxes,wgts,nprocs,sort,&eff);
    }
}

void
DistributionMapping::RRSFCDoIt (const BoxArray&          boxes,
                                int                      nprocs)
{
    BL_PROFILE("DistributionMapping::RRSFCDoIt()");

#if defined (BL_USE_TEAM)
    amrex::Abort("Team support is not implemented yet in RRSFC");
#endif

    const int nboxes = static_cast<int>(boxes.size());
    std::vector<SFCToken> tokens;
    tokens.reserve(nboxes);
    for (int i = 0; i < nboxes; ++i)
    {
        const Box& bx = boxes[i];
        tokens.push_back(makeSFCToken(i, bx.smallEnd()));
    }
    //
    // Put'm in Morton space filling curve order.
    //
    std::sort(tokens.begin(), tokens.end(), SFCToken::Compare());

    Vector<int> ord;

    LeastUsedCPUs(nprocs,ord);

    // Distribute boxes using roundrobin
    for (int i = 0; i < nboxes; ++i) {
        m_ref->m_pmap[i] = ParallelContext::local_to_global_rank(ord[i%nprocs]);
    }
}

void
DistributionMapping::RRSFCProcessorMap (const BoxArray&          boxes,
                                        int                      nprocs)
{
    BL_ASSERT( ! boxes.empty());

    m_ref->clear();
    m_ref->m_pmap.resize(boxes.size());

    RRSFCDoIt(boxes,nprocs);
}

DistributionMapping
DistributionMapping::makeKnapSack (const Vector<Real>& rcost, int nmax)
{
    BL_PROFILE("makeKnapSack");

    DistributionMapping r;

    Vector<Long> cost(rcost.size());

    Real wmax = *std::max_element(rcost.begin(), rcost.end());
    Real scale = (wmax == 0) ? 1.e9_rt : 1.e9_rt/wmax;

    for (int i = 0; i < rcost.size(); ++i) {
        cost[i] = Long(rcost[i]*scale) + 1L;
    }

    int nprocs = ParallelContext::NProcsSub();
    Real eff;

    r.KnapSackProcessorMap(cost, nprocs, &eff, true, nmax);

    return r;
}

DistributionMapping
DistributionMapping::makeKnapSack (const Vector<Real>& rcost, Real& eff, int nmax, bool sort)
{
    BL_PROFILE("makeKnapSack");

    DistributionMapping r;

    Vector<Long> cost(rcost.size());

    Real wmax = *std::max_element(rcost.begin(), rcost.end());
    Real scale = (wmax == 0) ? 1.e9_rt : 1.e9_rt/wmax;

    for (int i = 0; i < rcost.size(); ++i) {
        cost[i] = Long(rcost[i]*scale) + 1L;
    }

    int nprocs = ParallelContext::NProcsSub();

    r.KnapSackProcessorMap(cost, nprocs, &eff, true, nmax, sort);

    return r;
}

DistributionMapping
DistributionMapping::makeKnapSack (const LayoutData<Real>& rcost_local,
                                   Real& currentEfficiency, Real& proposedEfficiency,
                                   int nmax, bool broadcastToAll, int root,
                                   Real keep_ratio)
{
    BL_PROFILE("makeKnapSack");

    // Proposed distribution mapping is computed from global vector of costs on root;
    // required information is gathered on root from the layoutData information
    //
    // Two main steps:
    // 1. collect from rcost_local into the global cost vector rcost; then rcost is
    //    complete (only) on root
    // 2. (optional; default true) Broadcast processor map of the new dm to others

    Vector<Real> rcost(rcost_local.size());
    ParallelDescriptor::GatherLayoutDataToVector<Real>(rcost_local, rcost, root);
    // rcost is now filled out on root

    DistributionMapping r;
    if (ParallelDescriptor::MyProc() == root)
    {
        Vector<Long> cost(rcost.size());

        Real wmax = *std::max_element(rcost.begin(), rcost.end());
        Real scale = (wmax == 0) ? 1.e9_rt : 1.e9_rt/wmax;

        for (int i = 0; i < rcost.size(); ++i) {
            cost[i] = Long(rcost[i]*scale) + 1L;
        }

        if (keep_ratio > Real(0.0)) {
            r.KnapSackProcessorMap(rcost_local.DistributionMap(), cost, keep_ratio,
                                   currentEfficiency, proposedEfficiency, nmax);
        } else {
            int nprocs = ParallelDescriptor::NProcs();
            // `sort` needs to be false here since there's a parallel reduce function
            // in the processor map function, but we are executing only on root
            r.KnapSackProcessorMap(cost,nprocs,&proposedEfficiency,true,nmax,false);
            ComputeDistributionMappingEfficiency(rcost_local.DistributionMap(),
                                                 rcost, &currentEfficiency);
        }
    }

#ifdef BL_USE_MPI
    // Load-balanced distribution mapping is computed on root; broadcast the cost
    // to all proc (optional)
    if (broadcastToAll)
    {
        if (ParallelDescriptor::MyProc() == root)
        {
            auto const& pmap = r.ProcessorMap();
            ParallelDescriptor::Bcast(const_cast<int*>(pmap.data()), pmap.size(), root);
        }
        else
        {
            Vector<int> pmap(rcost_local.DistributionMap().size());
            ParallelDescriptor::Bcast(pmap.data(), pmap.size(), root);
            r = DistributionMapping(std::move(pmap));
        }
    }
#else
    amrex::ignore_unused(broadcastToAll);
#endif

    return r;
}

namespace {
Vector<Long>
gather_weights (const MultiFab& weight)
{
#ifdef AMREX_USE_MPI
    LayoutData<Real> costld(weight.boxArray(),weight.DistributionMap());
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(weight); mfi.isValid(); ++mfi) {
        costld[mfi] = weight[mfi].sum<RunOn::Device>(mfi.validbox(),0);
    }
    Vector<Real> rcost(weight.size());
    ParallelDescriptor::GatherLayoutDataToVector(costld, rcost,
                                                 ParallelContext::IOProcessorNumberSub());
    ParallelDescriptor::Bcast(rcost.data(), rcost.size(), ParallelContext::IOProcessorNumberSub());
    Real wmax = *std::max_element(rcost.begin(), rcost.end());
    Real scale = (wmax == 0) ? 1.e9_rt : 1.e9_rt/wmax;
    Vector<Long> lcost(rcost.size());
    for (int i = 0; i < rcost.size(); ++i) {
        lcost[i] = static_cast<Long>(rcost[i]*scale) + 1L;
    }
    return lcost;
#else
    return Vector<Long>(weight.size(), 1L);
#endif
}
}

DistributionMapping
DistributionMapping::makeKnapSack (const MultiFab& weight, int nmax)
{
    BL_PROFILE("makeKnapSack");
    Vector<Long> cost = gather_weights(weight);
    int nprocs = ParallelContext::NProcsSub();
    Real eff;
    DistributionMapping r;
    r.KnapSackProcessorMap(cost, nprocs, &eff, true, nmax);
    return r;
}

DistributionMapping
DistributionMapping::makeKnapSack (const MultiFab& weight, Real& eff, int nmax)
{
    BL_PROFILE("makeKnapSack");
    Vector<Long> cost = gather_weights(weight);
    int nprocs = ParallelContext::NProcsSub();
    DistributionMapping r;
    r.KnapSackProcessorMap(cost, nprocs, &eff, true, nmax);
    return r;
}

DistributionMapping
DistributionMapping::makeRoundRobin (const MultiFab& weight)
{
    BL_PROFILE("makeRoundRobin");
    Vector<Long> cost = gather_weights(weight);
    int nprocs = ParallelContext::NProcsSub();
    DistributionMapping r;
    r.RoundRobinProcessorMap(cost, nprocs);
    return r;
}

DistributionMapping
DistributionMapping::makeSFC (const MultiFab& weight, bool sort)
{
    BL_PROFILE("makeSFC");
    Vector<Long> cost = gather_weights(weight);
    int nprocs = ParallelContext::NProcsSub();
    DistributionMapping r;
    r.SFCProcessorMap(weight.boxArray(), cost, nprocs, sort);
    return r;
}

DistributionMapping
DistributionMapping::makeSFC (const MultiFab& weight, Real& eff, bool sort)
{
    BL_PROFILE("makeSFC");
    Vector<Long> cost = gather_weights(weight);
    int nprocs = ParallelContext::NProcsSub();
    DistributionMapping r;
    r.SFCProcessorMap(weight.boxArray(), cost, nprocs, eff, sort);
    return r;
}

DistributionMapping
DistributionMapping::makeSFC (const Vector<Real>& rcost, const BoxArray& ba, bool sort)
{
    BL_PROFILE("makeSFC");

    DistributionMapping r;

    Vector<Long> cost(rcost.size());

    Real wmax = *std::max_element(rcost.begin(), rcost.end());
    Real scale = (wmax == 0) ? 1.e9_rt : 1.e9_rt/wmax;

    for (int i = 0; i < rcost.size(); ++i) {
        cost[i] = Long(rcost[i]*scale) + 1L;
    }

    int nprocs = ParallelContext::NProcsSub();

    r.SFCProcessorMap(ba, cost, nprocs, sort);

    return r;
}

DistributionMapping
DistributionMapping::makeSFC (const Vector<Real>& rcost, const BoxArray& ba, Real& eff, bool sort)
{
    BL_PROFILE("makeSFC");

    DistributionMapping r;

    Vector<Long> cost(rcost.size());

    Real wmax = *std::max_element(rcost.begin(), rcost.end());
    Real scale = (wmax == 0) ? 1.e9_rt : 1.e9_rt/wmax;

    for (int i = 0; i < rcost.size(); ++i) {
        cost[i] = Long(rcost[i]*scale) + 1L;
    }

    int nprocs = ParallelContext::NProcsSub();

    r.SFCProcessorMap(ba, cost, nprocs, eff, sort);

    return r;
}

DistributionMapping
DistributionMapping::makeSFC (const LayoutData<Real>& rcost_local,
                              Real& currentEfficiency, Real& proposedEfficiency,
                              bool broadcastToAll, int root)
{
    BL_PROFILE("makeSFC");

    // Proposed distribution mapping is computed from global vector of costs on root;
    // required information is gathered on root from the layoutData information
    //
    // Two main steps:
    // 1. collect from rcost_local into the global cost vector rcost; then rcost is
    //    complete (only) on root
    // 2. (optional; default true) Broadcast processor map of the new dm to others

    Vector<Real> rcost(rcost_local.size());
    ParallelDescriptor::GatherLayoutDataToVector<Real>(rcost_local, rcost, root);
    // rcost is now filled out on root;

    DistributionMapping r;
    if (ParallelDescriptor::MyProc() == root)
    {
        Vector<Long> cost(rcost.size());

        Real wmax = *std::max_element(rcost.begin(), rcost.end());
        Real scale = (wmax == 0) ? 1.e9_rt : 1.e9_rt/wmax;

        for (int i = 0; i < rcost.size(); ++i) {
            cost[i] = Long(rcost[i]*scale) + 1L;
        }

        // `sort` needs to be false here since there's a parallel reduce function
        // in the processor map function, but we are executing only on root
        int nprocs = ParallelDescriptor::NProcs();
        r.SFCProcessorMap(rcost_local.boxArray(), cost, nprocs, proposedEfficiency, false);

        ComputeDistributionMappingEfficiency(rcost_local.DistributionMap(),
                                             rcost,
                                             &currentEfficiency);
    }

#ifdef BL_USE_MPI
    // Load-balanced distribution mapping is computed on root; broadcast the cost
    // to all proc (optional)
    if (broadcastToAll)
    {
        Vector<int> pmap(rcost_local.DistributionMap().size());
        if (ParallelDescriptor::MyProc() == root)
        {
            pmap = r.ProcessorMap();
        }

        // Broadcast vector from which to construct new distribution mapping
        ParallelDescriptor::Bcast(&pmap[0], pmap.size(), root);
        if (ParallelDescriptor::MyProc() != root)
        {
            r = DistributionMapping(pmap);
        }
    }
#else
    amrex::ignore_unused(broadcastToAll);
#endif

    return r;
}

std::vector<std::vector<int> >
DistributionMapping::makeSFC (const BoxArray& ba, bool use_box_vol, const int nprocs)
{
    BL_PROFILE("makeSFC");

    const int N = static_cast<int>(ba.size());
    std::vector<SFCToken> tokens;
    std::vector<Long> wgts;
    tokens.reserve(N);
    wgts.reserve(N);
    Long vol_sum = 0;
    for (int i = 0; i < N; ++i)
    {
        const Box& bx = ba[i];
        tokens.push_back(makeSFCToken(i, bx.smallEnd()));
        const Long v = use_box_vol ? bx.numPts() : Long(1);
        vol_sum += v;
        wgts.push_back(v);
    }
    //
    // Put'm in Morton space filling curve order.
    //
    std::sort(tokens.begin(), tokens.end(), SFCToken::Compare());

    Real volper = static_cast<Real>(vol_sum) / static_cast<Real>(nprocs);

    std::vector< std::vector<int> > r(nprocs);

    Distribute(tokens, wgts, nprocs, volper, r);

    return r;
}

const Vector<int>&
DistributionMapping::getIndexArray ()
{
    if (m_ref->m_index_array.empty())
    {
        int myProc = ParallelDescriptor::MyProc();

        for(int i = 0, N = static_cast<int>(m_ref->m_pmap.size()); i < N; ++i) {
            int rank = m_ref->m_pmap[i];
            if (ParallelDescriptor::sameTeam(rank)) {
                // If Team is not used (i.e., team size == 1), distributionMap[i] == myProc
                m_ref->m_index_array.push_back(i);
                m_ref->m_ownership.push_back(myProc == rank);
            }
        }
    }
    return m_ref->m_index_array;
}

const std::vector<bool>&
DistributionMapping::getOwnerShip ()
{
    if (m_ref->m_ownership.empty())
    {
        int myProc = ParallelDescriptor::MyProc();

        for(int i = 0, N = static_cast<int>(m_ref->m_pmap.size()); i < N; ++i) {
            int rank = m_ref->m_pmap[i];
            if (ParallelDescriptor::sameTeam(rank)) {
                // If Team is not used (i.e., team size == 1), distributionMap[i] == myProc
                m_ref->m_index_array.push_back(i);
                m_ref->m_ownership.push_back(myProc == rank);
            }
        }
    }
    return m_ref->m_ownership;
}

std::ostream&
operator<< (std::ostream&              os,
            const DistributionMapping& pmap)
{
    os << "(DistributionMapping" << '\n';

    for (int i = 0; i < pmap.ProcessorMap().size(); ++i)
    {
        os << "m_pmap[" << i << "] = " << pmap.ProcessorMap()[i] << '\n';
    }

    os << ')' << '\n';

    if (os.fail())
        amrex::Error("operator<<(ostream &, DistributionMapping &) failed");

    return os;
}

std::ostream&
operator<< (std::ostream& os, const DistributionMapping::RefID& id)
{
    os << id.data;
    return os;
}

std::istream&
DistributionMapping::readFrom (std::istream& is)
{
    AMREX_ASSERT(size() == 0);
    m_ref->clear();
    auto& pmap = m_ref->m_pmap;

    int n;
    is.ignore(100000, '(') >> n;
    pmap.resize(n);
    for (auto& x : pmap) {
        is >> x;
    }
    is.ignore(100000, ')');
    if (is.fail()) {
        amrex::Error("DistributionMapping::readFrom(istream&) failed");
    }
    return is;
}

std::ostream&
DistributionMapping::writeOn (std::ostream& os) const
{
    os << '(' << size() << '\n';
    for (int i = 0; i < size(); ++i) {
        os << (*this)[i] << '\n';
    }
    os << ')';
    if (os.fail()) {
        amrex::Error("DistributionMapping::writeOn(ostream&) failed");
    }
    return os;
}

DistributionMapping MakeSimilarDM (const BoxArray& ba, const MultiFab& mf, const IntVect& ng)
{
    const DistributionMapping& mf_dm = mf.DistributionMap();
    const BoxArray& mf_ba = convert(mf.boxArray(),ba.ixType());
    return MakeSimilarDM(ba, mf_ba, mf_dm, ng);
}

DistributionMapping MakeSimilarDM (const BoxArray& ba, const BoxArray& src_ba,
                                   const DistributionMapping& src_dm, const IntVect& ng)
{
    AMREX_ASSERT_WITH_MESSAGE(ba.ixType() == src_ba.ixType(),
                              "input BoxArrays must have the same centering.";);

    Vector<int> pmap(ba.size());
    for (int i = 0; i < static_cast<int>(ba.size()); ++i) {
        Box box = ba[i];
        box.grow(ng);
        bool first_only = false;
        auto isects = src_ba.intersections(box, first_only, ng);
        if (isects.empty()) {
            // no intersection found, revert to round-robin
            int nprocs = ParallelContext::NProcsSub();
            pmap[i] = i % nprocs;
        } else {
            Long max_overlap = 0;
            int max_overlap_index = -1;
            for (const auto& isect : isects) {
                int gid = isect.first;
                const Box& isect_box = isect.second;
                if (isect_box.numPts() > max_overlap) {
                    max_overlap = isect_box.numPts();
                    max_overlap_index = gid;
                }
            }
            AMREX_ASSERT(max_overlap > 0);
            pmap[i] = src_dm[max_overlap_index];
        }
    }
    return DistributionMapping(std::move(pmap));
}

}
