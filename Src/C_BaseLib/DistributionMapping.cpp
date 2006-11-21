
#include <winstd.H>

#include <Profiler.H>
#include <BoxArray.H>
#include <DistributionMapping.H>
#include <ParallelDescriptor.H>
#include <ParmParse.H>

#include <iostream>
#include <cstdlib>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <queue>
#include <algorithm>
#include <numeric>

static int    swap_n_test_count          = 1;
static int    verbose                    = 0;
static int    sfc_threshold              = 8;
static double max_efficiency             = 0.95;
static bool   do_not_minimize_comm_costs = true;
//
// Everyone uses the same Strategy -- defaults to KNAPSACK.
//
DistributionMapping::Strategy
DistributionMapping::m_Strategy = DistributionMapping::KNAPSACK;

DistributionMapping::PVMF
DistributionMapping::m_BuildMap = &DistributionMapping::KnapSackProcessorMap;

const Array<int>&
DistributionMapping::ProcessorMap () const
{
    return m_ref->m_pmap;
}

DistributionMapping::Strategy
DistributionMapping::strategy ()
{
    return DistributionMapping::m_Strategy;
}

int
DistributionMapping::CacheSize ()
{
    return m_Cache.size();
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
    default:
        BoxLib::Error("Bad DistributionMapping::Strategy");
    }
}

//
// We start out uninitialized.
//
bool DistributionMapping::m_Initialized = false;

bool
DistributionMapping::operator== (const DistributionMapping& rhs) const
{
    return m_ref->m_pmap == rhs.m_ref->m_pmap;
}

bool
DistributionMapping::operator!= (const DistributionMapping& rhs) const
{
    return !operator==(rhs);
}

void
DistributionMapping::Initialize ()
{
    DistributionMapping::m_Initialized = true;
        
    ParmParse pp("DistributionMapping");

    pp.query("verbose", verbose);

    pp.query("efficiency", max_efficiency);

    pp.query("do_not_minimize_comm_costs", do_not_minimize_comm_costs);

    pp.query("swap_n_test_count", swap_n_test_count);

    if (swap_n_test_count <= 0)
        BoxLib::Abort("swap_n_test must be integer >= 1");

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

            pp.query("sfc_threshold", sfc_threshold);
        }
        else
        {
            std::string msg("Unknown strategy: ");
            msg += theStrategy;
            BoxLib::Warning(msg.c_str());
        }
    }
}

void
DistributionMapping::Finalize ()
{}

//
// Our cache of processor maps.
//
std::vector< LnClassPtr<DistributionMapping::Ref> > DistributionMapping::m_Cache;

int
DistributionMapping::WhereToStart (int nprocs)
{
    if (m_Cache.size() == 0) return 0;
    //
    // What's the largest proc_map in our cache?
    //
    int N = m_Cache[0]->m_pmap.size();
    for (int i = 1; i < m_Cache.size(); i++)
        if (N < m_Cache[i]->m_pmap.size())
            N = m_Cache[i]->m_pmap.size();

    N--; // Subtract off the extra space reserved for the sentinel value.

    BL_ASSERT(N > 0);

    if (N < nprocs)
        //
        // We've got idle CPU(s); start distribution from there.
        //
        return N;
    //
    // We don't have any idle CPUs; find the one with the fewest boxes.
    //
    std::vector<int> count(nprocs,0);

    for (int i = 0; i < m_Cache.size(); i++)
        for (int j = 0; j < m_Cache[i]->m_pmap.size() - 1; j++)
            count[m_Cache[i]->m_pmap[j]]++;

    N = 0;
    for (int i = 1; i < count.size(); i++)
        if (count[N] > count[i])
            N = i;

    BL_ASSERT(N >= 0 && N < nprocs);

    return N;
}

bool
DistributionMapping::GetMap (const BoxArray& boxes)
{
    const int N = boxes.size();

    BL_ASSERT(m_ref->m_pmap.size() == N + 1);
    //
    // Search from back to front ...
    //
    for (int i = m_Cache.size() - 1; i >= 0; i--)
    {
        if (m_Cache[i]->m_pmap.size() == N + 1)
        {
            m_ref = m_Cache[i];

            BL_ASSERT(m_ref->m_pmap[N] == ParallelDescriptor::MyProc());

            return true;
        }
    }

    return false;
}

DistributionMapping::Ref::Ref () {}

DistributionMapping::DistributionMapping ()
    :
    m_ref(new DistributionMapping::Ref)
{}

DistributionMapping::Ref::Ref (const Array<int>& pmap)
    :
    m_pmap(pmap)
{}

DistributionMapping::DistributionMapping (const Array<int>& pmap)
    :
    m_ref(new DistributionMapping::Ref(pmap))
{}

DistributionMapping::Ref::Ref (int len)
    :
    m_pmap(len)
{}

DistributionMapping::DistributionMapping (const BoxArray& boxes, int nprocs)
    :
    m_ref(new DistributionMapping::Ref(boxes.size() + 1))
{
    define(boxes,nprocs);
}

DistributionMapping::Ref::Ref (const Ref& rhs)
    :
    m_pmap(rhs.m_pmap)
{}

DistributionMapping::DistributionMapping (const DistributionMapping& d1,
                                          const DistributionMapping& d2)
    :
    m_ref(new DistributionMapping::Ref(d1.size() + d2.size() - 1))

{
    const Array<int>& pmap_1 = d1.ProcessorMap();
    const Array<int>& pmap_2 = d2.ProcessorMap();

    const int L1 = pmap_1.size() - 1; // Length not including sentinel.
    const int L2 = pmap_2.size() - 1; // Length not including sentinel.

    for (int i = 0; i < L1; i++)
        m_ref->m_pmap[i] = pmap_1[i];

    for (int i = L1, j = 0; j < L2; i++, j++)
        m_ref->m_pmap[i] = pmap_2[j];
    //
    // Set sentinel equal to our processor number.
    //
    m_ref->m_pmap[m_ref->m_pmap.size()-1] = ParallelDescriptor::MyProc();
}

void
DistributionMapping::define (const BoxArray& boxes, int nprocs)
{
    if (m_ref->m_pmap.size() != boxes.size() + 1)
        m_ref->m_pmap.resize(boxes.size() + 1);

    if (!GetMap(boxes))
    {
        (this->*m_BuildMap)(boxes,nprocs);

        if (ParallelDescriptor::NProcs() > 1)
            //
            // We always append new processor maps.
            //
            DistributionMapping::m_Cache.push_back(m_ref);
    }
}

DistributionMapping::~DistributionMapping () {}

void
DistributionMapping::FlushCache ()
{
    DistributionMapping::CacheStats(std::cout);
    DistributionMapping::m_Cache.clear();
}

void
DistributionMapping::AddToCache (const DistributionMapping& dm)
{
    //
    // No need to maintain a cache when running in serial.
    //
    if (ParallelDescriptor::NProcs() < 2) return;

    bool              doit = true;
    const Array<int>& pmap = dm.ProcessorMap();

    if (pmap.size() > 0)
    {
        BL_ASSERT(pmap[pmap.size()-1] == ParallelDescriptor::MyProc());

        for (unsigned int i = 0; i < m_Cache.size() && doit; i++)
        {
            if (pmap.size() == m_Cache[i]->m_pmap.size())
            {
                BL_ASSERT(pmap == m_Cache[i]->m_pmap);

                doit = false;
            }
        }

        if (doit)
            m_Cache.push_back(dm.m_ref);
    }
}

void
DistributionMapping::RoundRobinProcessorMap (int nboxes, int nprocs)
{
    BL_ASSERT(nboxes > 0);

    m_ref->m_pmap.resize(nboxes+1);

    const int N = WhereToStart(nprocs);

    for (int i = 0; i < nboxes; i++)
        //
        // Start the round-robin at processor N.
        //
        m_ref->m_pmap[i] = (i + N) % nprocs;
    //
    // Set sentinel equal to our processor number.
    //
    m_ref->m_pmap[nboxes] = ParallelDescriptor::MyProc();
}

void
DistributionMapping::RoundRobinProcessorMap (const BoxArray& boxes, int nprocs)
{
    BL_ASSERT(boxes.size() > 0);
    BL_ASSERT(m_ref->m_pmap.size() == boxes.size() + 1);

    const int N = WhereToStart(nprocs);

    for (int i = 0; i < boxes.size(); i++)
        //
        // Start the round-robin at processor N.
        //
        m_ref->m_pmap[i] = (i + N) % nprocs;
    //
    // Set sentinel equal to our processor number.
    //
    m_ref->m_pmap[boxes.size()] = ParallelDescriptor::MyProc();
}

class WeightedBox
{
    int  m_boxid;
    long m_weight;
public:
    WeightedBox () {}
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
    std::list<WeightedBox>* m_lb;
    long                    m_weight;
public:
    WeightedBoxList (std::list<WeightedBox>* lb) : m_lb(lb), m_weight(0) {}
    long weight () const
    {
        return m_weight;
    }
    void erase (std::list<WeightedBox>::iterator& it)
    {
        m_weight -= it->weight();
        m_lb->erase(it);
    }
    void push_back (const WeightedBox& bx)
    {
        m_weight += bx.weight();
        m_lb->push_back(bx);
    }
    std::list<WeightedBox>::const_iterator begin () const { return m_lb->begin(); }
    std::list<WeightedBox>::iterator begin ()             { return m_lb->begin(); }
    std::list<WeightedBox>::const_iterator end () const   { return m_lb->end();   }
    std::list<WeightedBox>::iterator end ()               { return m_lb->end();   }

    bool operator< (const WeightedBoxList& rhs) const
    {
        return weight() > rhs.weight();
    }
};

static
std::vector< std::list<int> >
knapsack (const std::vector<long>& wgts, int nprocs)
{
    BL_PROFILE("knapsack()");

    const Real strttime = ParallelDescriptor::second();
    //
    // Sort balls by size largest first.
    //
    static std::list<int> empty_list;  // Work-around MSVC++ bug :-(

    std::vector< std::list<int> > result(nprocs, empty_list);

    std::vector<WeightedBox> lb;
    lb.reserve(wgts.size());
    for (unsigned int i = 0; i < wgts.size(); ++i)
    {
        lb.push_back(WeightedBox(i, wgts[i]));
    }
    BL_ASSERT(lb.size() == wgts.size());
    std::sort(lb.begin(), lb.end());
    BL_ASSERT(lb.size() == wgts.size());
    //
    // For each ball, starting with heaviest, assign ball to the lightest box.
    //
    std::priority_queue<WeightedBoxList>   wblq;
    std::vector< std::list<WeightedBox>* > vbbs(nprocs);
    for (int i  = 0; i < nprocs; ++i)
    {
        vbbs[i] = new std::list<WeightedBox>;
        wblq.push(WeightedBoxList(vbbs[i]));
    }
    BL_ASSERT(int(wblq.size()) == nprocs);
    for (unsigned int i = 0; i < wgts.size(); ++i)
    {
        WeightedBoxList wbl = wblq.top();
        wblq.pop();
        wbl.push_back(lb[i]);
        wblq.push(wbl);
    }
    BL_ASSERT(int(wblq.size()) == nprocs);
    std::list<WeightedBoxList> wblqg;
    while (!wblq.empty())
    {
        wblqg.push_back(wblq.top());
        wblq.pop();
    }
    BL_ASSERT(int(wblqg.size()) == nprocs);
    wblqg.sort();
    //
    // Compute the max weight and the sum of the weights.
    //
    double max_weight = 0;
    double sum_weight = 0;
    std::list<WeightedBoxList>::iterator it = wblqg.begin();
    for ( ; it != wblqg.end(); ++it)
    {
        long wgt = (*it).weight();
        sum_weight += wgt;
        max_weight = (wgt > max_weight) ? wgt : max_weight;
    }

    double efficiency = sum_weight/(nprocs*max_weight);

top:

    std::list<WeightedBoxList>::iterator it_top = wblqg.begin();

    WeightedBoxList wbl_top = *it_top;
    //
    // For each ball in the heaviest box.
    //
    std::list<WeightedBox>::iterator it_wb = wbl_top.begin();

    if (efficiency > max_efficiency) goto bottom;

    for ( ; it_wb != wbl_top.end(); ++it_wb )
    {
        //
        // For each ball not in the heaviest box.
        //
        std::list<WeightedBoxList>::iterator it_chk = it_top;
        it_chk++;
        for ( ; it_chk != wblqg.end(); ++it_chk)
        {
            WeightedBoxList wbl_chk = *it_chk;
            std::list<WeightedBox>::iterator it_owb = wbl_chk.begin();
            for ( ; it_owb != wbl_chk.end(); ++it_owb)
            {
                //
                // If exchanging these two balls reduces the load balance,
                // then exchange them and go to top.  The way we are doing
                // things, sum_weight cannot change.  So the efficiency will
                // increase if after we switch the two balls *it_wb and
                // *it_owb the max weight is reduced.
                //
                double w_tb = (*it_top).weight() + (*it_owb).weight() - (*it_wb).weight();
                double w_ob = (*it_chk).weight() + (*it_wb).weight() - (*it_owb).weight();
                //
                // If the other ball reduces the weight of the top box when
                // swapped, then it will change the efficiency.
                //
                if (w_tb < (*it_top).weight() && w_ob < (*it_top).weight())
                {
                    //
                    // Adjust the sum weight and the max weight.
                    //
                    WeightedBox wb = *it_wb;
                    WeightedBox owb = *it_owb;
                    wblqg.erase(it_top);
                    wblqg.erase(it_chk);
                    wbl_top.erase(it_wb);
                    wbl_chk.erase(it_owb);
                    wbl_top.push_back(owb);
                    wbl_chk.push_back(wb);
                    std::list<WeightedBoxList> tmp;
                    tmp.push_back(wbl_top);
                    tmp.push_back(wbl_chk);
                    tmp.sort();
                    wblqg.merge(tmp);
                    max_weight = (*wblqg.begin()).weight();
                    efficiency = sum_weight/(nprocs*max_weight);
                    goto top;
                }
            }
        }
    }

 bottom:
    //
    // Here I am "load-balanced".
    //
    std::list<WeightedBoxList>::const_iterator cit = wblqg.begin();
    for (int i = 0; i < nprocs; ++i)
    {
        const WeightedBoxList& wbl = *cit;
        std::list<WeightedBox>::const_iterator it1 = wbl.begin();
        for ( ; it1 != wbl.end(); ++it1)
        {
            result[i].push_back((*it1).boxid());
        }
        ++cit;
    }

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        const Real stoptime = ParallelDescriptor::second() - strttime;

        std::cout << "KNAPSACK efficiency: " << efficiency << '\n';
        std::cout << "KNAPSACK time: " << stoptime << '\n';
    }

    for (int i  = 0; i < nprocs; i++) delete vbbs[i];

    return result;
}

static
void
SwapAndTest (const std::map< int,std::vector<int>,std::greater<int> >& samesize,
             const std::vector< std::vector<int> >&                    nbrs,
             std::vector<int>&                                         procmap,
             std::vector<long>&                                        percpu)
{
    for (std::map< int,std::vector<int>,std::greater<int> >::const_iterator it = samesize.begin();
         it != samesize.end();
         ++it)
    {
        for (std::vector<int>::const_iterator lit1 = it->second.begin();
             lit1 != it->second.end();
             ++lit1)
        {
            std::vector<int>::const_iterator lit2 = lit1;

            const int ilit1 = *lit1;

            lit2++;

            for ( ; lit2 != it->second.end(); ++lit2)
            {
                const int ilit2 = *lit2;

                BL_ASSERT(ilit1 != ilit2);
                //
                // Don't consider Boxes on the same CPU.
                //
                if (procmap[ilit1] == procmap[ilit2]) continue;
                //
                // Will swapping these boxes decrease latency?
                //
                const long percpu_lit1 = percpu[procmap[ilit1]];
                const long percpu_lit2 = percpu[procmap[ilit2]];
                //
                // Now change procmap & redo necessary calculations ...
                //
                std::swap(procmap[ilit1],procmap[ilit2]);

                const int pmap1 = procmap[ilit1];
                const int pmap2 = procmap[ilit2];
                //
                // Update percpu[] in place.
                //
                std::vector<int>::const_iterator end1 = nbrs[ilit1].end();

                for (std::vector<int>::const_iterator it = nbrs[ilit1].begin(); it != end1; ++it)
                {
                    const int pmapstar = procmap[*it];

                    if (pmapstar == pmap2)
                    {
                        percpu[pmap1]++;
                        percpu[pmap2]++;
                    }
                    else if (pmapstar == pmap1)
                    {
                        percpu[pmap1]--;
                        percpu[pmap2]--;
                    }
                    else
                    {
                        percpu[pmap2]--;
                        percpu[pmap1]++;
                    }
                }

                std::vector<int>::const_iterator end2 = nbrs[ilit2].end();

                for (std::vector<int>::const_iterator it = nbrs[ilit2].begin(); it != end2; ++it)
                {
                    const int pmapstar = procmap[*it];

                    if (pmapstar == pmap1)
                    {
                        percpu[pmap1]++;
                        percpu[pmap2]++;
                    }
                    else if (pmapstar == pmap2)
                    {
                        percpu[pmap1]--;
                        percpu[pmap2]--;
                    }
                    else
                    {
                        percpu[pmap1]--;
                        percpu[pmap2]++;
                    }
                }

                const long cost_old = percpu_lit1  + percpu_lit2;
                const long cost_new = percpu[pmap1]+ percpu[pmap2];

                if (cost_new >= cost_old)
                {
                    //
                    // Undo our changes ...
                    //
                    std::swap(procmap[ilit1],procmap[ilit2]);

                    percpu[procmap[ilit1]] = percpu_lit1;
                    percpu[procmap[ilit2]] = percpu_lit2;
                }
            }
        }
    }
}

static
void
Output_CPU_Comm_Counts (const BoxArray&                        ba,
                        const std::vector< std::vector<int> >& nbrs,
                        std::vector<int>&                      procmap)
{
    const int NProcs = ParallelDescriptor::NProcs();
    const int MyProc = ParallelDescriptor::MyProc();
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    //
    // Each CPU must count # of CPUs it communicates with.
    //
    std::set<int> others;

    for (int i = 0; i < ba.size(); i++)
    {
        if (procmap[i] == MyProc)
        {
            for (int j = 0; j < nbrs[i].size(); j++)
            {
                if (procmap[nbrs[i][j]] != MyProc)
                    others.insert(procmap[nbrs[i][j]]);
            }
        }
    }

    int count = others.size();

    Array<int> counts(NProcs, 0);

#ifdef BL_USE_MPI
    MPI_Gather(&count,
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               counts.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               IOProc,
               ParallelDescriptor::Communicator());
#else
    counts[0] = count;
#endif

    typedef std::multimap<int, int, std::greater<int> > IIMap;

    IIMap iimap;

    for (int i = 0; i < NProcs; i++)
    {
        iimap.insert(std::make_pair(counts[i],i));
    }

    if (ParallelDescriptor::IOProcessor())
    {
        const int N = 5;
        int       i = 0;
        std::cout << "High communicators:\n";
        for (IIMap::const_iterator it = iimap.begin(); it != iimap.end() && i < N; ++it, ++i)
        {
            std::cout << "  cpu # " << it->second
                      << " communicates with " << it->first
                      << " other CPUs\n";
        }
    }
}

static
std::vector< std::vector<int> >
CalculateNeighbors (const BoxArray& ba)
{
    std::vector< std::vector<int> > nbrs(ba.size());
    //
    // Our "grow" factor; i.e. how far our tentacles grope for our neighbors.
    //
    const int Ngrow = 1;

    BoxArray grown(ba.size());

    for (int i = 0; i < ba.size(); i++)
    {
        grown.set(i, BoxLib::grow(ba[i],Ngrow));
    }

    for (int i = 0; i < grown.size(); i++)
    {
        std::vector< std::pair<int,Box> > isects = ba.intersections(grown[i]);

        for (int j = 0; j < isects.size(); j++)
            if (isects[j].first != i)
                nbrs[i].push_back(isects[j].first);
    }

    if (verbose > 1 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "The neighbors list:\n";

        for (int i = 0; i < nbrs.size(); i++)
        {
            std::cout << i << "\t:";

            for (std::vector<int>::const_iterator it = nbrs[i].begin(); it != nbrs[i].end(); ++it)
            {
                std::cout << *it << ' ';
            }

            std::cout << "\n";
        }
    }

    return nbrs;
}

//
// Try to "improve" the knapsack()d procmap ...
//

static
void
MinimizeCommCosts (std::vector<int>&        procmap,
                   const BoxArray&          ba,
                   const std::vector<long>& wgts,
                   int                      nprocs)
{
    BL_PROFILE("MinimizeCommCosts()");

    BL_ASSERT(ba.size() == wgts.size());
    BL_ASSERT(procmap.size() >= ba.size());

    if (nprocs < 2 || do_not_minimize_comm_costs) return;

    const Real strttime = ParallelDescriptor::second();

    std::vector< std::vector<int> > nbrs = CalculateNeighbors(ba);
    //
    // Want lists of box IDs having the same size.
    //
    std::map< int,std::vector<int>,std::greater<int> > samesize;

    for (int i = 0; i < wgts.size(); i++)
    {
        samesize[wgts[i]].push_back(i);
    }

    if (verbose > 1 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "Boxes sorted via numPts():\n";

        for (std::map< int,std::vector<int>,std::greater<int> >::const_iterator it = samesize.begin();
             it != samesize.end();
             ++it)
        {
            std::cout << it->first << "\t:";

            for (std::vector<int>::const_iterator lit = it->second.begin();
                 lit != it->second.end();
                 ++lit)
            {
                std::cout << *lit << ' ';
            }

            std::cout << '\n';
        }
    }
    //
    // Build a data structure to maintain the latency count on a per-CPU basis.
    //
    std::vector<long> percpu(nprocs,0L);

    for (int i = 0; i < nbrs.size(); i++)
    {
        for (std::vector<int>::const_iterator it = nbrs[i].begin(); it != nbrs[i].end(); ++it)
        {
            if (procmap[i] != procmap[*it])
                percpu[procmap[*it]]++;
        }
    }

    long initial_conn_count = 0;

    if (verbose > 1 && ParallelDescriptor::IOProcessor())
        for (int i = 0; i < percpu.size(); i++)
            initial_conn_count += percpu[i];

    for (int i = 0; i < swap_n_test_count; i++)
        SwapAndTest(samesize,nbrs,procmap,percpu);

    if (verbose > 1)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            const Real stoptime = ParallelDescriptor::second() - strttime;
            long final_conn_count = 0;
            for (int i = 0; i < percpu.size(); i++) final_conn_count += percpu[i];
            std::cout << "MinimizeCommCosts() time: " << stoptime << '\n';
            std::cout << "Initial off-CPU connection count: " << initial_conn_count << '\n';
            std::cout << "Final   off-CPU connection count: " << final_conn_count   << '\n';
        }

        Output_CPU_Comm_Counts(ba, nbrs, procmap);
    }
}

//
// This version does NOT call the MinimizeCommCosts() stuff.
//
void
DistributionMapping::KnapSackProcessorMap (const std::vector<long>& wgts,
                                           int                      nprocs)
{
    BL_ASSERT(wgts.size() > 0);

    m_ref->m_pmap.resize(wgts.size() + 1);

    if (wgts.size() <= nprocs || nprocs < 2)
    {
        RoundRobinProcessorMap(wgts.size(),nprocs);
    }
    else
    {
        const int N = WhereToStart(nprocs);

        std::vector< std::list<int> > vec = knapsack(wgts,nprocs);

        BL_ASSERT(vec.size() == nprocs);

        for (unsigned int i = 0; i < vec.size(); i++)
        {
            const int where = (i + N) % nprocs;

            for (std::list<int>::iterator lit = vec[i].begin(); lit != vec[i].end(); ++lit)
                m_ref->m_pmap[*lit] = where;
        }
        //
        // Set sentinel equal to our processor number.
        //
        m_ref->m_pmap[wgts.size()] = ParallelDescriptor::MyProc();
    }
}

//
// This version does call the MinimizeCommCosts() stuff.
//
void
DistributionMapping::KnapSackProcessorMap (const BoxArray& boxes,
					   int             nprocs)
{
    BL_ASSERT(boxes.size() > 0);
    BL_ASSERT(m_ref->m_pmap.size() == boxes.size()+1);

    if (boxes.size() <= nprocs || nprocs < 2)
    {
        RoundRobinProcessorMap(boxes,nprocs);
    }
    else
    {
        const int N = WhereToStart(nprocs);

        std::vector<long> wgts(boxes.size());

        for (unsigned int i = 0; i < wgts.size(); i++)
            wgts[i] = boxes[i].numPts();

        std::vector< std::list<int> > vec = knapsack(wgts,nprocs);

        BL_ASSERT(int(vec.size()) == nprocs);

        for (unsigned int i = 0; i < vec.size(); i++)
        {
            const int where = (i + N) % nprocs;

            for (std::list<int>::iterator lit = vec[i].begin(); lit != vec[i].end(); ++lit)
                m_ref->m_pmap[*lit] = where;
        }
        //
        // Set sentinel equal to our processor number.
        //
        m_ref->m_pmap[boxes.size()] = ParallelDescriptor::MyProc();

        if (verbose > 1)
        {
            std::vector< std::vector<int> > nbrs = CalculateNeighbors(boxes);

            Output_CPU_Comm_Counts(boxes, nbrs, m_ref->m_pmap);
        }

	MinimizeCommCosts(m_ref->m_pmap,boxes,wgts,nprocs);
    }
}

struct SFCToken
{
    class Compare
    {
    public:
        bool operator () (const SFCToken& lhs,
                          const SFCToken& rhs) const;
    };

    int     m_box;
    IntVect m_idx;
    Real    m_vol;

    static int MaxPower;
};

int SFCToken::MaxPower = 64;

bool
SFCToken::Compare::operator () (const SFCToken& lhs,
                                const SFCToken& rhs) const
{
    for (int i = SFCToken::MaxPower - 1; i >= 0; --i)
    {
        const int N = (1<<i);

        for (int j = BL_SPACEDIM-1; j >= 0; --j)
        {
            if (lhs.m_idx[j]/N < rhs.m_idx[j]/N)
            {
                return true;
            }
            else if (lhs.m_idx[j]/N > rhs.m_idx[j]/N)
            {
                return false;
            }
        }
    }

    return false;
}

void
DistributionMapping::SFCProcessorMap (const BoxArray& boxes,
                                      int             nprocs)
{
    std::vector<long> wgts(boxes.size());

    for (int i = 0; i < wgts.size(); i++)
        wgts[i] = boxes[i].volume();

    SFCProcessorMap(boxes, wgts, nprocs);
}

void
DistributionMapping::SFCProcessorMap (const BoxArray&          boxes,
                                      const std::vector<long>& wgts,
                                      int                      nprocs)
{
    BL_ASSERT(boxes.size() > 0);
    BL_ASSERT(boxes.size() == wgts.size());

    m_ref->m_pmap.resize(wgts.size()+1);

    if (boxes.size() <= nprocs || nprocs < 2)
    {
        RoundRobinProcessorMap(boxes,nprocs);
    }
    else if (boxes.size() <= sfc_threshold*nprocs)
    {
        KnapSackProcessorMap(wgts,nprocs);
    }
    else
    {
        const Real strttime = ParallelDescriptor::second();

        std::vector<SFCToken> token(boxes.size());

        int maxijk = 0;

        for (int i = 0; i < boxes.size(); i++)
        {
            token[i].m_box = i;
            token[i].m_idx = boxes[i].smallEnd();
            token[i].m_vol = wgts[i];

            for (int j = 0; j < BL_SPACEDIM; j++)
                maxijk = std::max(maxijk, token[i].m_idx[j]);
        }
        //
        // Set SFCToken::MaxPower for BoxArray.
        //
        int m = 0;
        for ( ; (1<<m) <= maxijk; m++)
            ;
        SFCToken::MaxPower = m;
        //
        // Put'm in Morton space filling curve order.
        //
        std::sort(token.begin(), token.end(), SFCToken::Compare());
        //
        // Split'm up as equitably as possible per CPU.
        //
        Real volpercpu = 0;
        for (int i = 0; i < token.size(); i++)
            volpercpu += token[i].m_vol;
        volpercpu /= nprocs;

        const int N = WhereToStart(nprocs);
        int       K = 0;

        for (int i = 0; i < nprocs; i++)
        {
            Real      vol = 0;
            int       cnt = 0;
            const int cpu = (i + N) % nprocs;

            for ( ;
                  K < token.size() && (i == (nprocs-1) || vol < volpercpu);
                  cnt++, K++)
            {
                vol += token[K].m_vol;
                m_ref->m_pmap[token[K].m_box] = cpu;
            }

            if (vol > volpercpu && cnt > 1)
            {
                const Real rdiff = vol - volpercpu;
                const Real ldiff = volpercpu - (vol - token[K-1].m_vol);
                if (rdiff > ldiff)
                    K--;
            }
        }
        //
        // Set sentinel equal to our processor number.
        //
        m_ref->m_pmap[boxes.size()] = ParallelDescriptor::MyProc();

        if (verbose)
        {
            const Real stoptime = ParallelDescriptor::second() - strttime;

            std::vector<Real> wgt(nprocs,0);

            for (int i = 0; i < token.size(); i++)
                wgt[m_ref->m_pmap[token[i].m_box]] += token[i].m_vol;

            Real sum_wgt = 0, max_wgt = 0;
            for (int i = 0; i < wgt.size(); i++)
            {
                if (wgt[i] > max_wgt) max_wgt = wgt[i];
                sum_wgt += wgt[i];
            }

            if (ParallelDescriptor::IOProcessor())
            {
                std::cout << "SFC efficiency: " << (sum_wgt/(nprocs*max_wgt)) << '\n';
                std::cout << "SFC time: " << stoptime << '\n';
            }

            if (verbose > 1)
            {
                std::vector< std::vector<int> > nbrs = CalculateNeighbors(boxes);

                Output_CPU_Comm_Counts(boxes, nbrs, m_ref->m_pmap);
            }
        }
    }
}

void
DistributionMapping::CacheStats (std::ostream& os)
{
    if (false && verbose && ParallelDescriptor::IOProcessor())
    {
        os << "The DistributionMapping cache contains "
           << DistributionMapping::m_Cache.size()
           << " Processor Map(s):\n";

        if (!DistributionMapping::m_Cache.empty())
        {
            for (unsigned int i = 0; i < m_Cache.size(); i++)
            {
                os << "\tMap #"
                   << i
                   << ", size:  "
                   << m_Cache[i]->m_pmap.size()
                   << ", references: "
                   << m_Cache[i].linkCount()
                   << '\n';
            }
            os << '\n';
        }
    }
}

std::ostream&
operator<< (std::ostream&              os,
            const DistributionMapping& pmap)
{
    os << "(DistributionMapping" << '\n';
    //
    // Do not print the sentinel value.
    //
    for (int i = 0; i < pmap.ProcessorMap().size() - 1; i++)
    {
        os << "m_pmap[" << i << "] = " << pmap.ProcessorMap()[i] << '\n';
    }

    os << ')' << '\n';

    if (os.fail())
        BoxLib::Error("operator<<(ostream &, DistributionMapping &) failed");

    return os;
}
