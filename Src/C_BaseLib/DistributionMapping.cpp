//BL_COPYRIGHT_NOTICE

//
// $Id: DistributionMapping.cpp,v 1.28 1998-07-24 22:41:14 lijewski Exp $
//

#include <DistributionMapping.H>
#include <List.H>
#include <ParallelDescriptor.H>
#include <ParmParse.H>
#include <RunStats.H>

#ifdef BL_USE_NEW_HFILES
#include <iostream>
#include <cstdlib>
#include <list>
#include <vector>
#include <queue>
#include <algorithm>
using namespace std;
#else
#include <iostream.h>
#include <stdlib.h>
#include <vector.h>
#if (defined(BL_T3E) && !defined(__KCC)) || defined(__GNUC__) || defined(BL_CRAY)
#include <stack.h>
#else
#include <queue.h>
#endif
#include <list.h>
#if (defined(BL_T3E) && !defined(__KCC)) || defined(__GNUC__) || defined(BL_CRAY)
#include <algo.h>
#else
#include <algorithm.h>
#endif
#endif /*BL_USE_NEW_HFILES*/

//
// Everyone uses the same Strategy -- defaults to ROUNDROBIN.
//
DistributionMapping::Strategy
DistributionMapping::m_Strategy = DistributionMapping::ROUNDROBIN;

DistributionMapping::PVMF
DistributionMapping::m_BuildMap = &DistributionMapping::RoundRobinProcessorMap;

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
    return m_procmap == rhs.m_procmap;
}

bool
DistributionMapping::operator!= (const DistributionMapping& rhs) const
{
    return !operator==(rhs);
}

void
DistributionMapping::init ()
{
    DistributionMapping::m_Initialized = true;
        
    ParmParse pp("DistributionMapping");

    aString theStrategy;

    if (pp.query("strategy", theStrategy))
    {
        DistributionMapping::Strategy how = ROUNDROBIN;

        if (theStrategy == "ROUNDROBIN")
        {
            how = ROUNDROBIN;
        }
        else if (theStrategy == "KNAPSACK")
        {
            how = KNAPSACK;
        }
        else
        {
            aString msg("Unknown strategy: ");
            msg += theStrategy;
            BoxLib::Warning(msg.c_str());
        }

        strategy(how);
    }
}

//
// Our cache of processor maps.
//
vector< Array<int> > DistributionMapping::m_Cache;

bool
DistributionMapping::GetMap (int             nprocs,
                             const BoxArray& boxes)
{
    const int N = boxes.length();

    assert(m_procmap.length() == N+1);
    //
    // Search from back to front ...
    //
    for (int i = m_Cache.size() - 1; i >= 0; i--)
    {
        if (m_Cache[i].length() == N+1)
        {
            const Array<int>& cached_procmap = m_Cache[i];

            for (int i = 0; i <= N; i++)
            {
                m_procmap[i] = cached_procmap[i];
            }

            assert(m_procmap[N] == ParallelDescriptor::MyProc());

            return true;
        }
    }

    return false;
}

DistributionMapping::DistributionMapping ()
{
    if (!m_Initialized)
    {
        DistributionMapping::init();
    }
}

DistributionMapping::DistributionMapping (int             nprocs,
                                          const BoxArray& boxes)
    :
    m_procmap(boxes.length()+1)
{
    if (!m_Initialized)
    {
        DistributionMapping::init();
    }

    define(nprocs, boxes);
}

void
DistributionMapping::define (int             nprocs,
                             const BoxArray& boxes)
{
    static RunStats stats("processor_map");

    stats.start();

    if (!(m_procmap.length() == boxes.length()+1))
    {
        m_procmap.resize(boxes.length()+1);
    }

    if (DistributionMapping::m_Strategy == DistributionMapping::ROUNDROBIN)
    {
        //
        // Don't bother with the cache -- just do it :-)
        //
        (this->*m_BuildMap)(nprocs, boxes);
    }
    else
    {
        if (!GetMap(nprocs, boxes))
        {
            (this->*m_BuildMap)(nprocs, boxes);
            //
            // We always append new processor maps.
            //
            DistributionMapping::m_Cache.push_back(m_procmap);
        }
    }

    stats.end();
}

DistributionMapping::~DistributionMapping () {}

void
DistributionMapping::FlushCache ()
{
    DistributionMapping::m_Cache.clear();
}

void
DistributionMapping::RoundRobinProcessorMap (int             nprocs,
                                             const BoxArray& boxes)
{
    assert(nprocs > 0);
    assert(boxes.length() > 0);
    assert(m_procmap.length() == boxes.length()+1);

    for (int i = 0; i < boxes.length(); i++)
    {
        m_procmap[i] = i % nprocs;
    }
    //
    // Set sentinel equal to our processor number.
    //
    m_procmap[boxes.length()] = ParallelDescriptor::MyProc();
}

#if defined(BL_USE_BSP) || defined(BL_USE_MPI)

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
    list<WeightedBox> m_lb;
    long              m_weight;
public:
    WeightedBoxList() : m_weight(0) {}
    long weight () const
    {
        return m_weight;
    }
    void erase (list<WeightedBox>::iterator& it)
    {
        m_weight -= (*it).weight();
        m_lb.erase(it);
    }
    void push_back (const WeightedBox& bx)
    {
        m_weight += bx.weight();
        m_lb.push_back(bx);
    }
    list<WeightedBox>::const_iterator begin () const { return m_lb.begin(); }
    list<WeightedBox>::iterator begin ()             { return m_lb.begin(); }
    list<WeightedBox>::const_iterator end () const   { return m_lb.end();   }
    list<WeightedBox>::iterator end ()               { return m_lb.end();   }

    bool operator< (const WeightedBoxList& rhs) const
    {
        return weight() > rhs.weight();
    }
};

static
vector< list<int> >
knapsack (const vector<long>& pts,
          int                 nprocs)
{
    assert(nprocs > 0);
    //
    // Sort balls by size largest first.
    //
    static list<int> empty_list;  // Work-around MSVC++ bug :-(

    vector< list<int> > result(nprocs, empty_list);

    vector<WeightedBox> lb;
    lb.reserve(pts.size());
    for (int i = 0; i < pts.size(); ++i)
    {
        lb.push_back(WeightedBox(i, pts[i]));
    }
    assert(lb.size() == pts.size());
    sort(lb.begin(), lb.end());
    assert(lb.size() == pts.size());
    //
    // For each ball, starting with heaviest, assign ball to the lightest box.
    //
    priority_queue<WeightedBoxList> wblq;
    for (int i  = 0; i < nprocs; ++i)
    {
        wblq.push(WeightedBoxList());
    }
    assert(wblq.size() == nprocs);
    for (int i = 0; i < pts.size(); ++i)
    {
        WeightedBoxList wbl = wblq.top();
        wblq.pop();
        wbl.push_back(lb[i]);
        wblq.push(wbl);
    }
    assert(wblq.size() == nprocs);
    list<WeightedBoxList> wblqg;
    while (!wblq.empty())
    {
        wblqg.push_back(wblq.top());
        wblq.pop();
    }
    assert(wblqg.size() == nprocs);
    wblqg.sort();
    //
    // Compute the max weight and the sum of the weights.
    //
    double max_weight = 0;
    double sum_weight = 0;
    list<WeightedBoxList>::iterator it = wblqg.begin();
    for ( ; it != wblqg.end(); ++it)
    {
        long wgt = (*it).weight();
        sum_weight += wgt;
        max_weight = (wgt > max_weight) ? wgt : max_weight;
    }
top:
    list<WeightedBoxList>::iterator it_top = wblqg.begin();
    list<WeightedBoxList>::iterator it_chk = it_top;
    it_chk++;
    WeightedBoxList wbl_top = *it_top;
    //
    // For each ball in the heaviest box.
    //
    list<WeightedBox>::iterator it_wb = wbl_top.begin();
    for ( ; it_wb != wbl_top.end(); ++it_wb )
    {
        //
        // For each ball not in the heaviest box.
        //
        for ( ; it_chk != wblqg.end(); ++it_chk)
        {
            WeightedBoxList wbl_chk = *it_chk;
            list<WeightedBox>::iterator it_owb = wbl_chk.begin();
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
                    list<WeightedBoxList> tmp;
                    tmp.push_back(wbl_top);
                    tmp.push_back(wbl_chk);
                    tmp.sort();
                    wblqg.merge(tmp);
                    max_weight = (*wblqg.begin()).weight();
                    //
                    //efficiency = sum_weight/(nprocs*max_weight);
                    //cout << "Efficiency = " << efficiency << '\n';
                    //cout << "max_weight = " << max_weight << '\n';
                    //
                    goto top;
                }
            }
        }
    }
    //
    // Here I am "load-balanced".
    //
    list<WeightedBoxList>::const_iterator cit = wblqg.begin();
    for (int i = 0; i < nprocs; ++i)
    {
        const WeightedBoxList& wbl = *cit;
        list<WeightedBox>::const_iterator it1 = wbl.begin();
        for ( ; it1 != wbl.end(); ++it1)
        {
            result[i].push_back((*it1).boxid());
        }
        ++cit;
    }

    return result;
}
#endif /*defined(BL_USE_BSP) || defined(BL_USE_MPI)*/

void
DistributionMapping::KnapSackProcessorMap (int             nprocs,
                                           const BoxArray& boxes)
{
    assert(nprocs > 0);
    assert(boxes.length() > 0);
    assert(m_procmap.length() == boxes.length()+1);

    if (boxes.length() <= nprocs || nprocs < 2)
    {
        //
        // Might as well just use ROUNDROBIN.
        //
        RoundRobinProcessorMap(nprocs, boxes);
    }
    else
    {
#if defined(BL_USE_BSP) || defined(BL_USE_MPI)
        vector<long> pts(boxes.length());

        for (int i = 0; i < pts.size(); i++)
        {
            pts[i] = boxes[i].numPts();
        }

        vector< list<int> > vec = knapsack(pts, nprocs);

        assert(vec.size() == nprocs);

        list<int>::iterator lit;

        for (int i = 0; i < vec.size(); i++)
        {
            for (lit = vec[i].begin(); lit != vec[i].end(); ++lit)
            {
                m_procmap[*lit] = i;
            }
        }
        //
        // Set sentinel equal to our processor number.
        //
        m_procmap[boxes.length()] = ParallelDescriptor::MyProc();
#else
        BoxLib::Error("How did this happen?");
#endif
    }
}

void
DistributionMapping::CacheStats (ostream& os)
{
    os << "The DistributionMapping cache contains "
       << DistributionMapping::m_Cache.size()
       << " Processor Map(s):\n";

    if (!DistributionMapping::m_Cache.empty())
    {
        for (int i = 0; i < m_Cache.size(); i++)
        {
            os << "\tMap #"
               << i
               << " is of length "
               << m_Cache[i].length()
               << '\n';
        }
        os << '\n';
    }
}

ostream&
operator<< (ostream&                   os,
            const DistributionMapping& pmap)
{
    os << "(DistributionMapping" << '\n';

    for (int i = 0; i < pmap.ProcessorMap().length(); i++)
    {
        os << "m_procmap[" << i << "] = " << pmap.ProcessorMap()[i] << '\n';
    }

    os << ')' << '\n';

    if (os.fail())
        BoxLib::Error("operator<<(ostream &, DistributionMapping &) failed");

    return os;
}
