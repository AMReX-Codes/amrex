
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_Utility.H>

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

DistributionMapping::PVMF DistributionMapping::m_BuildMap = 0;

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
    max_efficiency   = 0.9;
    node_size        = 0;
    flag_verbose_mapper = 0;

    ParmParse pp("DistributionMapping");

    pp.query("v"      ,             verbose);
    pp.query("verbose",             verbose);
    pp.query("efficiency",          max_efficiency);
    pp.query("sfc_threshold",       sfc_threshold);
    pp.query("node_size",           node_size);
    pp.query("verbose_mapper",      flag_verbose_mapper);

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

    DistributionMapping::m_BuildMap = 0;
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
        LIpairV.push_back(LIpair(bytes[i],i));
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

        LIpairV.push_back(LIpair(teambytes,i));
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

DistributionMapping::DistributionMapping ()
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

    BL_ASSERT(m_BuildMap != 0);

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
                                     std::vector<LIpair>* LIpairV)
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
	LeastUsedCPUs(nprocs,ord);
	wrkerord.resize(nprocs);
	for (int i = 0; i < nprocs; ++i) {
	    wrkerord[i].resize(1);
	    wrkerord[i][0] = 0;
	}
    } else {
	LeastUsedTeams(ord,wrkerord,nteams,nworkers);
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
DistributionMapping::RoundRobinProcessorMap (int nboxes, int nprocs)
{
    BL_ASSERT(nboxes > 0);
    m_ref->clear();
    m_ref->m_pmap.resize(nboxes);

    RoundRobinDoIt(nboxes, nprocs);
}

void
DistributionMapping::RoundRobinProcessorMap (const BoxArray& boxes, int nprocs)
{
    BL_ASSERT(boxes.size() > 0);
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

    const int N = boxes.size();

    LIpairV.reserve(N);

    for (int i = 0; i < N; ++i)
    {
        LIpairV.push_back(LIpair(boxes[i].numPts(),i));
    }

    Sort(LIpairV, true);

    RoundRobinDoIt(boxes.size(), nprocs, &LIpairV);
}


void
DistributionMapping::RoundRobinProcessorMap (const std::vector<Long>& wgts,
                                             int nprocs)
{
    BL_ASSERT(wgts.size() > 0);

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

    const int N = wgts.size();

    LIpairV.reserve(N);

    for (int i = 0; i < N; ++i)
    {
        LIpairV.push_back(LIpair(wgts[i],i));
    }

    Sort(LIpairV, true);

    RoundRobinDoIt(wgts.size(), nprocs, &LIpairV);
}

class WeightedBox
{
    int  m_boxid;
    Long m_weight;
public:
    WeightedBox (int b, int w) : m_boxid(b), m_weight(w) {}
    Long weight () const { return m_weight; }
    int  boxid ()  const { return m_boxid;  }

    bool operator< (const WeightedBox& rhs) const
    {
        return weight() > rhs.weight();
    }
};

class WeightedBoxList
{
    Vector<WeightedBox>* m_lb;
    Long                 m_weight;
public:
    WeightedBoxList (Long w) : m_lb(nullptr), m_weight(w) {}
    WeightedBoxList (Vector<WeightedBox>* lb) : m_lb(lb), m_weight(0) {}
    Long weight () const
    {
        return m_weight;
    }
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
    int size () const { return m_lb->size(); }
    Vector<WeightedBox>::const_iterator begin () const { return m_lb->begin(); }
    Vector<WeightedBox>::iterator begin ()             { return m_lb->begin(); }
    Vector<WeightedBox>::const_iterator end () const   { return m_lb->end();   }
    Vector<WeightedBox>::iterator end ()               { return m_lb->end();   }

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
    for (unsigned int i = 0, N = wgts.size(); i < N; ++i)
    {
        lb.push_back(WeightedBox(i, wgts[i]));
    }
    std::sort(lb.begin(), lb.end());
    //
    // For each ball, starting with heaviest, assign ball to the lightest bin.
    //
    std::priority_queue<WeightedBoxList> wblq;
    Vector<std::unique_ptr<Vector<WeightedBox> > > raii_vwb(nprocs);
    for (int i  = 0; i < nprocs; ++i)
    {
        raii_vwb[i].reset(new Vector<WeightedBox>);
        wblq.push(WeightedBoxList(raii_vwb[i].get()));
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
        Real wgt = wbl.weight();
        sum_weight += wgt;
        max_weight = std::max(wgt, max_weight);
    }

    while (!wblq.empty())
    {
	WeightedBoxList wbl = wblq.top();
        wblq.pop();
	if (wbl.size() > 0) {
	    Real wgt = wbl.weight();
	    sum_weight += wgt;
	    max_weight = std::max(wgt, max_weight);
	    wblv.push_back(wbl);
	}
    }

    efficiency = sum_weight/(nprocs*max_weight);

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

                        max_weight = bl_top->weight();
                        efficiency = sum_weight / (nprocs*max_weight);
                        goto top;
                    }
                }
            }
        }

        BL_ASSERT(std::is_sorted(wblv.begin(), wblv.end()));
    }

    for (int i = 0, N = wblv.size(); i < N; ++i)
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
        for (int i = 0; i < vec.size(); ++i) {
            Print() << "  Bucket " << i << " contains boxes:" << std::endl;
            for (int j = 0; j < vec[i].size(); ++j) {
                Print() << "    " << vec[i][j] << std::endl;
            }
        }
    }

    BL_ASSERT(static_cast<int>(vec.size()) == nteams);

    std::vector<LIpair> LIpairV;

    LIpairV.reserve(nteams);

    for (int i = 0; i < nteams; ++i)
    {
        Long wgt = 0;
        for (std::vector<int>::const_iterator lit = vec[i].begin(), End = vec[i].end();
             lit != End; ++lit)
        {
            wgt += wgts[*lit];
        }

        LIpairV.push_back(LIpair(wgt,i));
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
	const int N = vi.size();

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
    BL_ASSERT(wgts.size() > 0);

    m_ref->clear();
    m_ref->m_pmap.resize(wgts.size());

    if (static_cast<int>(wgts.size()) <= nprocs || nprocs < 2)
    {
        RoundRobinProcessorMap(wgts.size(),nprocs);

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
DistributionMapping::KnapSackProcessorMap (const BoxArray& boxes,
					   int             nprocs)
{
    BL_ASSERT(boxes.size() > 0);

    m_ref->m_pmap.resize(boxes.size());

    if (boxes.size() <= nprocs || nprocs < 2)
    {
        RoundRobinProcessorMap(boxes,nprocs);
    }
    else
    {
        std::vector<Long> wgts(boxes.size());

        for (unsigned int i = 0, N = boxes.size(); i < N; ++i)
            wgts[i] = boxes[i].numPts();

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
}

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

namespace {
#if (AMREX_SPACEDIM == 3)
    AMREX_FORCE_INLINE
    uint32_t make_space (uint32_t x)
    {
        // x            : 0000,0000,0000,0000,0000,00a9,8765,4321
        x = (x | (x << 16)) & 0x030000FF;
        // x << 16      : 0000,00a9,8765,4321,0000,0000,0000,0000
        // x | (x << 16): 0000,00a9,8765,4321,0000,00a9,8765,4321
        // 0x030000FF   : 0000,0011,0000,0000,0000,0000,1111,1111
        // x            : 0000,00a9,0000,0000,0000,0000,8765,4321
        x = (x | (x <<  8)) & 0x0300F00F;
        // x << 8       : 0000,0000,0000,0000,8765,4321,0000,0000
        // x | (x << 8) : 0000,00a9,0000,0000,8765,4321,8765,4321
        // 0x0300F00F   : 0000,0011,0000,0000,1111,0000,0000,1111
        // x            : 0000,00a9,0000,0000,8765,0000,0000,4321
        x = (x | (x <<  4)) & 0x030C30C3;
        // x << 4       : 00a9,0000,0000,8765,0000,0000,4321,0000
        // x | (x << 4) : 00a9,00a9,0000,8765,8765,0000,4321,4321
        // 0x030C30C3   : 0000,0011,0000,1100,0011,0000,1100,0011
        // x            : 0000,00a9,0000,8700,0065,0000,4300,0021
        x = (x | (x <<  2)) & 0x09249249;
        // x << 2       : 0000,a900,0087,0000,6500,0043,0000,2100
        // x | (x << 2) : 0000,a9a9,0087,8700,6565,0043,4300,2121
        // 0x09249249   : 0000,1001,0010,0100,1001,0010,0100,1001
        // x            : 0000,a009,0080,0700,6005,0040,0300,2001
        return x;
    }
#elif (AMREX_SPACEDIM == 2)
    AMREX_FORCE_INLINE
    uint32_t make_space (uint32_t x)
    {
        // x           : 0000,0000,0000,0000,gfed,cba9,8765,4321
        x = (x | (x << 8)) & 0x00FF00FF;
        // x << 8      : 0000,0000,gfed,cba9,8765,4321,0000,0000
        // x | (x << 8): 0000,0000,gfed,cba9,????,????,8765,4321
        // 0x00FF00FF  : 0000,0000,1111,1111,0000,0000,1111,1111
        // x           : 0000,0000,gfed,cba9,0000,0000,8765,4321
        x = (x | (x << 4)) & 0x0F0F0F0F;
        // x << 4      : 0000,gfed,cba9,0000,0000,8765,4321,0000
        // x | (x << 4): 0000,gfed,????,cba9,0000,8765,????,4321
        // 0x0F0F0F0F  : 0000,1111,0000,1111,0000,1111,0000,1111
        // x           : 0000,gfed,0000,cba9,0000,8765,0000,4321
        x = (x | (x << 2)) & 0x33333333;
        // x << 2      : 00gf,ed00,00cb,a900,0087,6500,0043,2100
        // x | (x << 2): 00gf,??ed,00cb,??a9,0087,??65,0043,??21
        // 0x33333333  : 0011,0011,0011,0011,0011,0011,0011,0011
        // x           : 00gf,00ed,00cb,00a9,0087,0065,0043,0021
        x = (x | (x << 1)) & 0x55555555;
        // x << 1      : 0gf0,0ed0,0cb0,0a90,0870,0650,0430,0210
        // x | (x << 1): 0g?f,0e?d,0c?b,0a?9,08?7,06?5,04?3,02?1
        // 0x55555555  : 0101,0101,0101,0101,0101,0101,0101,0101
        // x           : 0g0f,0e0d,0c0b,0a09,0807,0605,0403,0201
        return x;
    }
#endif

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
        token.m_morton[0] = make_space(x & 0x3FF)
                         | (make_space(y & 0x3FF) << 1)
                         | (make_space(z & 0x3FF) << 2);
        x = x >> 10;
        y = y >> 10;
        z = z >> 10;
        token.m_morton[1] = make_space(x & 0x3FF)
                         | (make_space(y & 0x3FF) << 1)
                         | (make_space(z & 0x3FF) << 2);
        x = x >> 10;
        y = y >> 10;
        z = z >> 10;
        token.m_morton[2] = make_space(x & 0x3FF)
                         | (make_space(y & 0x3FF) << 1)
                         | (make_space(z & 0x3FF) << 2);

#elif (AMREX_SPACEDIM == 2)

        constexpr uint32_t offset = 1u << 31;
        static_assert(static_cast<uint32_t>(std::numeric_limits<int>::max())+1 == offset,
                      "INT_MAX != (1<<31)-1");
        uint32_t x = (iv[0] >= 0) ? static_cast<uint32_t>(iv[0]) + offset
            : static_cast<uint32_t>(iv[0]-std::numeric_limits<int>::lowest());
        uint32_t y = (iv[1] >= 0) ? static_cast<uint32_t>(iv[1]) + offset
            : static_cast<uint32_t>(iv[1]-std::numeric_limits<int>::lowest());
        // extract lowest 16 bits and make sapce for interleaving
        token.m_morton[0] = make_space(x & 0xFFFF)
                         | (make_space(y & 0xFFFF) << 1);
        x = x >> 16;
        y = y >> 16;
        token.m_morton[1] = make_space(x) | (make_space(y) << 1);

#elif (AMREX_SPACEDIM == 1)

        constexpr uint32_t offset = 1u << 31;
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
            vol += wgts[tokens[K].m_box];
            ++cnt;

            v[i].push_back(tokens[K].m_box);
        }

        totalvol += vol;

        if ((totalvol/(i+1)) > volpercpu &&  // Too much for this bin.
            cnt > 1                      &&  // More than one box in this bin.
            i < nprocs-1)                    // Not the last bin, which has to take all.
        {
            --K;
            v[i].pop_back();
            totalvol -= wgts[tokens[K].m_box];
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
                rank_vol += wgts[t.m_box];
                idx++;
            }
            Print() << "    Total Rank Vol: " << rank_vol << std::endl;
        }
    }

#ifdef AMREX_DEBUG
    int cnt = 0;
    for (int i = 0; i < nprocs; ++i) {
        cnt += v[i].size();
    }
    BL_ASSERT(cnt == static_cast<int>(tokens.size()));
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

    const int N = boxes.size();
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
        volperteam += wt;
    }
    volperteam /= nteams;

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
        for (int j = 0, M = vi.size(); j < M; ++j)
            wgt += wgts[vi[j]];

        LIpairV.push_back(LIpair(wgt,i));
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
	const int Nbx = vi.size();                // # of boxes assigned to this team

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
		for (std::vector<int>::const_iterator it = kpres[w].begin();
		     it != kpres[w].end(); ++it)
		{
		    wgt += local_wgts[*it];
		}
		ww.push_back(LIpair(wgt,w));
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
		for (std::vector<int>::const_iterator it = js.begin(); it!=js.end(); ++it)
		    m_ref->m_pmap[vi[*it]] = cpu;
	    }
	}
    }

    if (eff || verbose)
    {
        Real sum_wgt = 0, max_wgt = 0;
        for (int i = 0; i < nteams; ++i)
        {
            const Long W = LIpairV[i].first;
            if (W > max_wgt) max_wgt = W;
            sum_wgt += W;
        }
        Real efficiency = (sum_wgt/(nteams*max_wgt));
        if (eff) *eff = efficiency;

        if (verbose)
        {
            amrex::Print() << "SFC efficiency: " << efficiency << '\n';
        }
    }
}

void
DistributionMapping::SFCProcessorMap (const BoxArray& boxes,
                                      int             nprocs)
{
    BL_ASSERT(boxes.size() > 0);

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

        for (int i = 0, N = boxes.size(); i < N; ++i)
        {
            wgts.push_back(boxes[i].volume());
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
    BL_ASSERT(boxes.size() > 0);
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
    BL_ASSERT(boxes.size() > 0);
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

    const int nboxes = boxes.size();
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
    BL_ASSERT(boxes.size() > 0);

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
    Real scale = (wmax == 0) ? 1.e9 : 1.e9/wmax;

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
    Real scale = (wmax == 0) ? 1.e9 : 1.e9/wmax;

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
                                   int nmax, bool broadcastToAll, int root)
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
        Real scale = (wmax == 0) ? 1.e9 : 1.e9/wmax;

        for (int i = 0; i < rcost.size(); ++i) {
            cost[i] = Long(rcost[i]*scale) + 1L;
        }

        // `sort` needs to be false here since there's a parallel reduce function
        // in the processor map function, but we are executing only on root
        int nprocs = ParallelDescriptor::NProcs();
        r.KnapSackProcessorMap(cost, nprocs, &proposedEfficiency, true, nmax, false);

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

void
DistributionMapping::ComputeDistributionMappingEfficiency (const DistributionMapping& dm,
                                                           const Vector<Real>& cost,
                                                           Real* efficiency)
{
    const int nprocs = ParallelDescriptor::NProcs();
        
    // This will store mapping from processor to the costs of FABs it controls,
    // (proc) --> ([cost_FAB_1, cost_FAB_2, ... ]),
    // for each proc
    Vector<Vector<Real>> rankToCosts(nprocs);

    // Count the number of costs belonging to each rank
    Vector<int> cnt(nprocs);
    for (int i=0; i<dm.size(); ++i)
    {
        ++cnt[dm[i]];
    }
    
    for (int i=0; i<rankToCosts.size(); ++i)
    {
        rankToCosts[i].reserve(cnt[i]);
    }
    
    for (int i=0; i<cost.size(); ++i)
    {
        rankToCosts[dm[i]].push_back(cost[i]);
    }

    Real maxCost = -1.0;

    // This will store mapping from (proc) --> (sum of cost) for each proc
    Vector<Real> rankToCost(nprocs);
    for (int i=0; i<nprocs; ++i)
    {
        const Real rwSum = std::accumulate(rankToCosts[i].begin(),
                                           rankToCosts[i].end(), 0.0);
        rankToCost[i] = rwSum;
        maxCost = std::max(maxCost, rwSum);
    }

    // Write `efficiency` (number between 0 and 1), the mean cost per processor
    // (normalized to the max cost)
    *efficiency = (std::accumulate(rankToCost.begin(),
                                   rankToCost.end(), 0.0) / (nprocs*maxCost));
}

namespace {
Vector<Long>
gather_weights (const MultiFab& weight)
{
#ifdef AMREX_USE_MPI
    LayoutData<Real> costld(weight.boxArray(),weight.DistributionMap());
#ifdef _OPENMP
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
    Real scale = (wmax == 0) ? 1.e9 : 1.e9/wmax;
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
    Real scale = (wmax == 0) ? 1.e9 : 1.e9/wmax;

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
    Real scale = (wmax == 0) ? 1.e9 : 1.e9/wmax;

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
        Real scale = (wmax == 0) ? 1.e9 : 1.e9/wmax;

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

    const int N = ba.size();
    std::vector<SFCToken> tokens;
    std::vector<Long> wgts;
    tokens.reserve(N);
    wgts.reserve(N);
    Long vol_sum = 0;
    for (int i = 0; i < N; ++i)
    {
        const Box& bx = ba[i];
        tokens.push_back(makeSFCToken(i, bx.smallEnd()));
        const Long v = use_box_vol ? bx.volume() : Long(1);
        vol_sum += v;
        wgts.push_back(v);
    }
    //
    // Put'm in Morton space filling curve order.
    //
    std::sort(tokens.begin(), tokens.end(), SFCToken::Compare());

    Real volper;
    volper = vol_sum / nprocs;

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

        for(int i = 0, N = m_ref->m_pmap.size(); i < N; ++i) {
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

        for(int i = 0, N = m_ref->m_pmap.size(); i < N; ++i) {
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

}
