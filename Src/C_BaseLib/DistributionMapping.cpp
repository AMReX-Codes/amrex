//BL_COPYRIGHT_NOTICE

//
// $Id: DistributionMapping.cpp,v 1.7 1997-11-25 17:08:25 lijewski Exp $
//

#include <DistributionMapping.H>
#include <ParallelDescriptor.H>

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

DistributionMapping::DistributionMapping ()
    :
    nProcessors(0),
    boxarray(),
    distributionStrategy(ROUNDROBIN),
    processorMap(),
    objectsPerProcessor(),
    nPtsPerProcessor()
{
    CreateProcessorMap();
}

DistributionMapping::DistributionMapping (int                  nprocessors,
                                          const BoxArray&      boxes,
                                          DistributionStrategy strategy)
    :
    nProcessors(nprocessors),
    boxarray(boxes),
    distributionStrategy(strategy),
    processorMap(boxes.length()),
    objectsPerProcessor(nprocessors),
    nPtsPerProcessor(nprocessors)
{
    CreateProcessorMap();
}

DistributionMapping::~DistributionMapping () {}

void
DistributionMapping::define (int                  nprocessors,
                             const BoxArray&      boxes,
                             DistributionStrategy strategy)
{
    nProcessors = nprocessors;
    boxarray = boxes;
    distributionStrategy = strategy;
    processorMap.resize(boxes.length());
    objectsPerProcessor.resize(nprocessors, 0);
    nPtsPerProcessor.resize(nprocessors);
    CreateProcessorMap();
}

void
DistributionMapping::CreateProcessorMap ()
{
    int i;
    switch (distributionStrategy)
    {
    case ROUNDROBIN:
        for (i = 0; i < processorMap.length(); i++)
        {
            processorMap[i] = i % nProcessors;
            ++objectsPerProcessor[processorMap[i]];
        }
        break;
    case RANDOM:
        BoxLib::Error("RANDOM not implemented");
        break;
    case KNAPSACK:
	// Bill's Algorithm:
	BoxLib::Error("KNAPSACK not implement");
	break;
    case SIZEBALANCED:
        BoxLib::Error("SIZEBALANCED not implemented");
        break;
    default:
        BoxLib::Error("Bad DistributionStrategy");
    }
}

int
DistributionMapping::operator () (int        level,
                                  const Box& box) const
{
    return -1;
}

ostream&
operator<< (ostream&                   os,
            const DistributionMapping& pmap)
{
    int i;

    os << "(DistributionMapping" << '\n';
    for (i = 0; i < pmap.processorMap.length(); i++)
    {
        os << "processorMap[" << i << "] = " << pmap.processorMap[i] << '\n';
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
public:
    long weight () const
    {
	long wgt = 0;
        list<WeightedBox>::const_iterator it = m_lb.begin();
	for( ; it != m_lb.end(); ++it)
	{
	    wgt += (*it).weight();
	}
	return wgt;
    }
    void erase (list<WeightedBox>::iterator& it)
    {
	m_lb.erase(it);
    }
    void push_back (const WeightedBox& bx)           { m_lb.push_back(bx);  }
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

vector< list<int> >
knapsack (const vector<long>& pts,
          int                 nProcessors)
{
    //
    // Sort balls by size largest first.
    //
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
    double efficiency = sum_weight/nProcessors/max_weight;
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
		if ( w_tb < (*it_top).weight() && w_ob < (*it_top).weight() )
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
		    efficiency = sum_weight/(nProcessors*max_weight);
		    // cout << "Efficiency = " << efficiency << '\n';
		    // cout << "max_weight = " << max_weight << '\n';
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
