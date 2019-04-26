
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_FArrayBox.H>
#if !defined(BL_NO_FORT)
#include <AMReX_Geometry.H>
#endif
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

    Vector<long> bytes(ParallelContext::NProcsSub());
    long thisbyte = amrex::TotalBytesAllocatedInFabs()/1024;
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

    Vector<long> bytes(ParallelContext::NProcsSub());
    long thisbyte = amrex::TotalBytesAllocatedInFabs()/1024;
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

	long teambytes = 0;
	int offset = i*nworkers;
	for (int j = 0; j < nworkers; ++j)
	{
	    int globalrank = offset+j;
	    long b = bytes[globalrank];
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
DistributionMapping::RoundRobinProcessorMap (const std::vector<long>& wgts,
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
    long m_weight;
public:
    WeightedBox (int b, int w) : m_boxid(b), m_weight(w) {}
    long weight () const { return m_weight; }
    int  boxid ()  const { return m_boxid;  }

    bool operator< (const WeightedBox& rhs) const
    {
        return weight() > rhs.weight();
    }
};

class WeightedBoxList
{
    Vector<WeightedBox>* m_lb;
    long                 m_weight;
public:
    WeightedBoxList (long w) : m_lb(nullptr), m_weight(w) {}
    WeightedBoxList (Vector<WeightedBox>* lb) : m_lb(lb), m_weight(0) {}
    long weight () const
    {
        return m_weight;
    }
    void addWeight (long dw) { m_weight += dw; }
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
knapsack (const std::vector<long>&         wgts,
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
        WeightedBoxList wbl = wblq.top();
        wblq.pop();
        wbl.push_back(lb[i]);
	if (wbl.size() < nmax) {
	    wblq.push(wbl);
	} else {
	    wblv.push_back(wbl);
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
            long w_top = bl_top->weight();
            long w_bottom = bl_bottom->weight();
            for (auto ball_1 = bl_top->begin(); ball_1 != bl_top->end(); ++ball_1)
            {
                for (auto ball_2 = bl_bottom->begin(); ball_2 != bl_bottom->end(); ++ball_2)
                {
                    // should we swap ball 1 and ball 2?
                    long dw = ball_1->weight() - ball_2->weight();
                    long w_top_new    = w_top    - dw;
                    long w_bottom_new = w_bottom + dw;
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
DistributionMapping::KnapSackDoIt (const std::vector<long>& wgts,
                                   int                    /*  nprocs */,
                                   Real&                    efficiency,
                                   bool                     do_full_knapsack,
				   int                      nmax)
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
	long wgt = 0;
        for (std::vector<int>::const_iterator lit = vec[i].begin(), End = vec[i].end();
             lit != End; ++lit)
        {
            wgt += wgts[*lit];
        }

        LIpairV.push_back(LIpair(wgt,i));
    }

    Sort(LIpairV, true);

    if (flag_verbose_mapper) {
        for (const auto &p : LIpairV) {
            Print() << "  Bucket " << p.second << " total weight: " << p.first << std::endl;
        }
    }

    Vector<int> ord;
    Vector<Vector<int> > wrkerord;

    if (nteams == nprocs) {
	LeastUsedCPUs(nprocs,ord);
	wrkerord.resize(nprocs);
	for (int i = 0; i < nprocs; ++i) {
	    wrkerord[i].resize(1);
	    wrkerord[i][0] = 0;
	}
    } else {
	LeastUsedTeams(ord,wrkerord,nteams,nworkers);
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
DistributionMapping::KnapSackProcessorMap (const std::vector<long>& wgts,
                                           int                      nprocs,
                                           Real*                    efficiency,
                                           bool                     do_full_knapsack,
					   int                      nmax)
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
        KnapSackDoIt(wgts, nprocs, eff, do_full_knapsack, nmax);
        if (efficiency) *efficiency = eff;
    }
}

void
DistributionMapping::KnapSackProcessorMap (const BoxArray& boxes,
					   int             nprocs)
{
    BL_ASSERT(boxes.size() > 0);
    BL_ASSERT(m_ref->m_pmap.size() == boxes.size());

    if (boxes.size() <= nprocs || nprocs < 2)
    {
        RoundRobinProcessorMap(boxes,nprocs);
    }
    else
    {
        std::vector<long> wgts(boxes.size());

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
            bool operator () (const SFCToken& lhs,
                              const SFCToken& rhs) const;
        };

        SFCToken (int box, const IntVect& idx, Real vol)
            :
            m_box(box), m_idx(idx), m_vol(vol) {}

        int     m_box;
        IntVect m_idx;
        Real    m_vol;

        static int MaxPower;
    };
}

int SFCToken::MaxPower = 64;

bool
SFCToken::Compare::operator () (const SFCToken& lhs,
                                const SFCToken& rhs) const
{
    for (int i = SFCToken::MaxPower - 1; i >= 0; --i)
    {
        const int N = (1<<i);

        for (int j = AMREX_SPACEDIM-1; j >= 0; --j)
        {
            const int il = lhs.m_idx[j]/N;
            const int ir = rhs.m_idx[j]/N;

            if (il < ir)
            {
                return true;
            }
            else if (il > ir)
            {
                return false;
            }
        }
    }
    return false;
}

static
void
Distribute (const std::vector<SFCToken>&     tokens,
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
                    << t.m_idx << ": "
                    << t.m_vol << std::endl;
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
            vol += tokens[K].m_vol;
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
            totalvol -= tokens[K].m_vol;
        }
    }

    if (flag_verbose_mapper) {
        Print() << "Distributed SFC Tokens:" << std::endl;
        int idx = 0;
        for (int i = 0; i < nprocs; ++i) {
            Print() << "  Rank/Team " << i << ":" << std::endl;
            Real rank_vol = 0;
            for (const auto &box : v[i]) {
                const auto &t = tokens[idx];
                BL_ASSERT(box == t.m_box);
                Print() << "    " << idx << ": "
                        << t.m_box << ": "
                        << t.m_idx << ": "
                        << t.m_vol << std::endl;
                rank_vol += t.m_vol;
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
                                          const std::vector<long>& wgts,
                                          int                   /*   nprocs */,
                                          bool                     sort)
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

    std::vector<SFCToken> tokens;

    const int N = boxes.size();

    tokens.reserve(N);

    int maxijk = 0;

    for (int i = 0; i < N; ++i)
    {
	const Box& bx = boxes[i];
        tokens.push_back(SFCToken(i,bx.smallEnd(),wgts[i]));

        const SFCToken& token = tokens.back();

        AMREX_D_TERM(maxijk = std::max(maxijk, token.m_idx[0]);,
                     maxijk = std::max(maxijk, token.m_idx[1]);,
                     maxijk = std::max(maxijk, token.m_idx[2]););
    }
    //
    // Set SFCToken::MaxPower for BoxArray.
    //
    int m = 0;
    for ( ; (1 << m) <= maxijk; ++m) {
        ;  // do nothing
    }
    SFCToken::MaxPower = m;
    //
    // Put'm in Morton space filling curve order.
    //
    std::sort(tokens.begin(), tokens.end(), SFCToken::Compare());
    //
    // Split'm up as equitably as possible per team.
    //
    Real volperteam = 0;
    for (const SFCToken& tok : tokens) {
        volperteam += tok.m_vol;
    }
    volperteam /= nteams;

    std::vector< std::vector<int> > vec(nteams);

    Distribute(tokens,nteams,volperteam,vec);

    // vec has a size of nteams and vec[] holds a vector of box ids.

    tokens.clear();

    std::vector<LIpair> LIpairV;

    LIpairV.reserve(nteams);

    for (int i = 0; i < nteams; ++i)
    {
	long wgt = 0;
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
	    std::vector<long> local_wgts;
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
		long wgt = 0;
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

    if (verbose)
    {
        Real sum_wgt = 0, max_wgt = 0;
        for (int i = 0; i < nteams; ++i)
        {
            const long W = LIpairV[i].first;
            if (W > max_wgt)
                max_wgt = W;
            sum_wgt += W;
        }

        amrex::Print() << "SFC efficiency: " << (sum_wgt/(nteams*max_wgt)) << '\n';
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
        std::vector<long> wgts;

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
                                      const std::vector<long>& wgts,
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
DistributionMapping::RRSFCDoIt (const BoxArray&          boxes,
				int                      nprocs)
{
    BL_PROFILE("DistributionMapping::RRSFCDoIt()");

#if defined (BL_USE_TEAM)
    amrex::Abort("Team support is not implemented yet in RRSFC");
#endif

    std::vector<SFCToken> tokens;

    const int nboxes = boxes.size();

    tokens.reserve(nboxes);

    int maxijk = 0;

    for (int i = 0; i < nboxes; ++i)
    {
	const Box& bx = boxes[i];
        tokens.push_back(SFCToken(i,bx.smallEnd(),0.0));

        const SFCToken& token = tokens.back();

        AMREX_D_TERM(maxijk = std::max(maxijk, token.m_idx[0]);,
               maxijk = std::max(maxijk, token.m_idx[1]);,
               maxijk = std::max(maxijk, token.m_idx[2]););
    }
    //
    // Set SFCToken::MaxPower for BoxArray.
    //
    int m = 0;
    for ( ; (1 << m) <= maxijk; ++m) {
        ;  // do nothing
    }
    SFCToken::MaxPower = m;
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
DistributionMapping::makeKnapSack (const Vector<Real>& rcost)
{
    BL_PROFILE("makeKnapSack");

    DistributionMapping r;

    Vector<long> cost(rcost.size());

    Real wmax = *std::max_element(rcost.begin(), rcost.end());
    Real scale = (wmax == 0) ? 1.e9 : 1.e9/wmax;

    for (int i = 0; i < rcost.size(); ++i) {
        cost[i] = long(rcost[i]*scale) + 1L;
    }

    int nprocs = ParallelContext::NProcsSub();
    Real eff;

    r.KnapSackProcessorMap(cost, nprocs, &eff, true);

    return r;
}

DistributionMapping
DistributionMapping::makeKnapSack (const MultiFab& weight, int nmax)
{
    BL_PROFILE("makeKnapSack");

    DistributionMapping r;

    Vector<long> cost(weight.size());
#ifdef BL_USE_MPI
    {
	Vector<Real> rcost(cost.size(), 0.0);
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(weight); mfi.isValid(); ++mfi) {
	    int i = mfi.index();
	    rcost[i] = weight[mfi].sum(mfi.validbox(),0);
	}

	ParallelAllReduce::Sum(&rcost[0], rcost.size(), ParallelContext::CommunicatorSub());

	Real wmax = *std::max_element(rcost.begin(), rcost.end());
	Real scale = (wmax == 0) ? 1.e9 : 1.e9/wmax;

	for (int i = 0; i < rcost.size(); ++i) {
	    cost[i] = long(rcost[i]*scale) + 1L;
	}
    }
#endif

    int nprocs = ParallelContext::NProcsSub();
    Real eff;

    r.KnapSackProcessorMap(cost, nprocs, &eff, true, nmax);

    return r;
}

DistributionMapping
DistributionMapping::makeRoundRobin (const MultiFab& weight)
{
    DistributionMapping r;

    Vector<long> cost(weight.size());
#ifdef BL_USE_MPI
    {
	Vector<Real> rcost(cost.size(), 0.0);
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(weight); mfi.isValid(); ++mfi) {
	    int i = mfi.index();
	    rcost[i] = weight[mfi].sum(mfi.validbox(),0);
	}

	ParallelAllReduce::Sum(&rcost[0], rcost.size(), ParallelContext::CommunicatorSub());

	Real wmax = *std::max_element(rcost.begin(), rcost.end());
        Real scale = (wmax == 0) ? 1.e9 : 1.e9/wmax;

	for (int i = 0; i < rcost.size(); ++i) {
	    cost[i] = long(rcost[i]*scale) + 1L;
	}
    }
#endif

    int nprocs = ParallelContext::NProcsSub();

    r.RoundRobinProcessorMap(cost, nprocs);

    return r;
}

DistributionMapping
DistributionMapping::makeSFC (const MultiFab& weight, bool sort)
{
    DistributionMapping r;

    Vector<long> cost(weight.size());
#ifdef BL_USE_MPI
    {
	Vector<Real> rcost(cost.size(), 0.0);
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(weight); mfi.isValid(); ++mfi) {
	    int i = mfi.index();
	    rcost[i] = weight[mfi].sum(mfi.validbox(),0);
	}

	ParallelAllReduce::Sum(&rcost[0], rcost.size(), ParallelContext::CommunicatorSub());

	Real wmax = *std::max_element(rcost.begin(), rcost.end());
        Real scale = (wmax == 0) ? 1.e9 : 1.e9/wmax;

	for (int i = 0; i < rcost.size(); ++i) {
	    cost[i] = long(rcost[i]*scale) + 1L;
	}
    }
#endif

    int nprocs = ParallelContext::NProcsSub();

    r.SFCProcessorMap(weight.boxArray(), cost, nprocs, sort);

    return r;
}

std::vector<std::vector<int> >
DistributionMapping::makeSFC (const BoxArray& ba, bool use_box_vol)
{
    std::vector<SFCToken> tokens;

    const int N = ba.size();

    tokens.reserve(N);

    int maxijk = 0;

    Real vol_sum = 0;
    for (int i = 0; i < N; ++i)
    {
	const Box& bx = ba[i];
        const auto & bx_vol = (use_box_vol ? bx.volume() : 1);
        tokens.push_back(SFCToken(i,bx.smallEnd(),bx_vol));
        vol_sum += bx_vol;

        const SFCToken& token = tokens.back();

        AMREX_D_TERM(maxijk = std::max(maxijk, token.m_idx[0]);,
                     maxijk = std::max(maxijk, token.m_idx[1]);,
                     maxijk = std::max(maxijk, token.m_idx[2]););
    }
    //
    // Set SFCToken::MaxPower for BoxArray.
    //
    int m = 0;
    for ( ; (1 << m) <= maxijk; ++m) {
        ;  // do nothing
    }
    SFCToken::MaxPower = m;
    //
    // Put'm in Morton space filling curve order.
    //
    std::sort(tokens.begin(), tokens.end(), SFCToken::Compare());

    const int nprocs = ParallelContext::NProcsSub();
    Real volper;
    volper = vol_sum / nprocs;

    std::vector< std::vector<int> > r(nprocs);
    Distribute(tokens, nprocs, volper, r);

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

}
