//BL_COPYRIGHT_NOTICE

//
// $Id: DistributionMapping.cpp,v 1.10 1997-11-25 19:55:15 car Exp $
//

#include <DistributionMapping.H>
#include <ParallelDescriptor.H>
#include <ParmParse.H>
#include <RunStats.H>

#ifdef BL_USE_NEW_HFILES
#include <iostream>
#include <cstdlib>
#include <vector>
#include <queue>
#include <list>
#include <algorithm>
using namespace std;
#else
#include <iostream.h>
#include <stdlib.h>
#include <vector.h>
#include <queue.h>
#include <list.h>
#include <algorithm.h>
#endif

//
// Everyone uses the same Strategy -- defaults to ROUNDROBIN.
//
DistributionMapping::Strategy
DistributionMapping::distributionStrategy = DistributionMapping::ROUNDROBIN;

void
DistributionMapping::strategy (DistributionMapping::Strategy how)
{
    DistributionMapping::distributionStrategy = how;
}

DistributionMapping::Strategy
DistributionMapping::strategy ()
{
    return DistributionMapping::distributionStrategy;
}

bool DistributionMapping::initialized = false;

//
// Forward declaration.
//
static vector< list<int> > knapsack (const vector<long>&, int);

void
DistributionMapping::init ()
{
    DistributionMapping::initialized = true;
        
    ParmParse pp("DistributionMapping");

    aString theStrategy;

    if (pp.query("strategy", theStrategy))
    {
        if (theStrategy == "ROUNDROBIN")
        {
            DistributionMapping::distributionStrategy = ROUNDROBIN;
        }
        else if (theStrategy == "KNAPSACK")
        {
            DistributionMapping::distributionStrategy = KNAPSACK;
        }
        else if (theStrategy == "RANDOM")
        {
            DistributionMapping::distributionStrategy = RANDOM;
        }
        else if (theStrategy == "SIZEBALANCED")
        {
            DistributionMapping::distributionStrategy = SIZEBALANCED;
        }
        else
        {
            aString msg("Unknown strategy: ");
            msg += theStrategy;
            BoxLib::Warning(msg.c_str());
        }
    }
}

DistributionMapping::DistributionMapping ()
    :
    nProcessors(0)
{
    if (!initialized)
        DistributionMapping::init();

    CreateProcessorMap();
}

DistributionMapping::DistributionMapping (int             nprocessors,
                                          const BoxArray& boxes)
    :
    nProcessors(nprocessors),
    boxarray(boxes),
    processorMap(boxes.length())
{
    if (!initialized)
        DistributionMapping::init();

    CreateProcessorMap();
}

DistributionMapping::~DistributionMapping () {}

void
DistributionMapping::define (int             nprocessors,
                             const BoxArray& boxes)
{
    nProcessors = nprocessors;
    boxarray = boxes;
    processorMap.resize(boxes.length());
    CreateProcessorMap();
}

void
DistributionMapping::CreateProcessorMap ()
{
    RunStats stats("processor_map");

    stats.start();

    switch (DistributionMapping::distributionStrategy)
    {
    case ROUNDROBIN:
    {
        for (int i = 0; i < processorMap.length(); i++)
        {
            processorMap[i] = i % nProcessors;
        }
    }
    break;
    case RANDOM:
        BoxLib::Error("RANDOM not implemented");
        break;
    case KNAPSACK:
    {
        vector<long> pts(boxarray.length());

        for (int i = 0; i < pts.size(); i++)
        {
            pts[i] = boxarray[i].numPts();
        }

        vector< list<int> > vec = knapsack(pts, nProcessors);

        list<int>::iterator lit;

        for (int i = 0; i < vec.size(); i++)
        {
            for (lit = vec[i].begin(); lit != vec[i].end(); ++lit)
            {
                processorMap[*lit] = i;
            }
        }
    }
    break;
    case SIZEBALANCED:
        BoxLib::Error("SIZEBALANCED not implemented");
        break;
    default:
        BoxLib::Error("Bad Strategy");
    }

    stats.end();
}

ostream&
operator<< (ostream&                   os,
            const DistributionMapping& pmap)
{
    os << "(DistributionMapping" << '\n';

    for (int i = 0; i < pmap.ProcessorMap().length(); i++)
    {
        os << "processorMap[" << i << "] = " << pmap.ProcessorMap()[i] << '\n';
    }

    os << ')' << '\n';

    if (os.fail())
        BoxLib::Error("operator<<(ostream &, DistributionMapping &) failed");

    return os;
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
};

inline
bool
operator< (const WeightedBox& lhs,
           const WeightedBox& rhs)
{
    return lhs.weight()  > rhs.weight();
}

class WeightedBoxList
{
    list<WeightedBox> m_lb;
    long m_weight;
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
};

inline
bool
operator< (const WeightedBoxList& lhs,
           const WeightedBoxList& rhs)
{
    return lhs.weight() > rhs.weight();
}

static
vector< list<int> >
knapsack (const vector<long>& pts,
          int                 nProcessors)
{
    //
    // Sort balls by size largest first.
    //
    assert(nProcessors > 0);
    vector< list<int> > result(nProcessors);
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
    for (int i  = 0; i < nProcessors; ++i)
    {
	wblq.push(WeightedBoxList());
    }
    assert(wblq.size() == nProcessors);
    for (int i = 0; i < pts.size(); ++i)
    {
	WeightedBoxList wbl = wblq.top();
	wblq.pop();
	wbl.push_back(lb[i]);
	wblq.push(wbl);
    }
    assert(wblq.size() == nProcessors);
    list<WeightedBoxList> wblqg;
    while (!wblq.empty())
    {
	wblqg.push_back(wblq.top());
        wblq.pop();
    }
    assert(wblqg.size() == nProcessors);
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
    cout << "sum_weight = " << sum_weight << '\n';
    cout << "max_weight = " << max_weight << '\n';
    double efficiency = 0;
    if (max_weight)
        efficiency = sum_weight/nProcessors/max_weight;
    cout << "Efficiency = " << efficiency << '\n';
top:
    list<WeightedBoxList>::iterator it_top = wblqg.begin();
    list<WeightedBoxList>::iterator it_chk = it_top; it_chk++;
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
		    // for(list<WeightedBoxList>::const_iterator ita = wblqg.begin(); ita != wblqg.end(); ++ita)
		    //	cout << ita->weight() << '\n';
		    wblqg.erase(it_top);
		    wblqg.erase(it_chk);
		    // cout << "weight_top = " << wbl_top.weight() << '\n';
		    wbl_top.erase(it_wb);
		    wbl_chk.erase(it_owb);
		    // cout << "weight_top = " << wbl_top.weight() << '\n';
		    wbl_top.push_back(owb);
		    // cout << "weight_top = " << wbl_top.weight() << '\n';
		    wbl_chk.push_back(wb);
		    list<WeightedBoxList> tmp;
		    tmp.push_back(wbl_top);
		    tmp.push_back(wbl_chk);
		    tmp.sort();
		    wblqg.merge(tmp);
		    // for(list<WeightedBoxList>::const_iterator ita = wblqg.begin(); ita != wblqg.end(); ++ita)
		    //	cout << ita->weight() << '\n';
		    max_weight = (*wblqg.begin()).weight();
                    if (max_weight)
                        efficiency = sum_weight/(nProcessors*max_weight);
                    //cout << "Efficiency = " << efficiency << '\n';
                    //cout << "max_weight = " << max_weight << '\n';
		    goto top;
		}
	    }
	}
    }
    //
    // Here I am "load-balanced".
    //
    list<WeightedBoxList>::const_iterator cit = wblqg.begin();
    for (int i = 0; i < nProcessors; ++i)
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
