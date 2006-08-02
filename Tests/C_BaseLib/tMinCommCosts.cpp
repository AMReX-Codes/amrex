#include <iostream>
#include <fstream>
#include <BoxArray.H>
#include <BoxDomain.H>
#include <ParallelDescriptor.H>

#include <iostream>
#include <cstdlib>
#include <list>
#include <map>
#include <vector>
#include <queue>
#include <algorithm>
#include <numeric>

static const int verbose = 1;
static double max_efficiency             = 0.95;
static bool   do_not_minimize_comm_costs = false;

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
knapsack (const std::vector<long>& pts, int nprocs)
{
    const Real strttime = ParallelDescriptor::second();
    //
    // Sort balls by size largest first.
    //
    static std::list<int> empty_list;  // Work-around MSVC++ bug :-(

    std::vector< std::list<int> > result(nprocs, empty_list);

    std::vector<WeightedBox> lb;
    lb.reserve(pts.size());
    for (unsigned int i = 0; i < pts.size(); ++i)
    {
        lb.push_back(WeightedBox(i, pts[i]));
    }
    BL_ASSERT(lb.size() == pts.size());
    std::sort(lb.begin(), lb.end());
    BL_ASSERT(lb.size() == pts.size());
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
    for (unsigned int i = 0; i < pts.size(); ++i)
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

    std::cout << "knapsack initial efficiency: " << efficiency << std::endl;
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

    const Real stoptime = ParallelDescriptor::second() - strttime;

    std::cout << "knapsack final efficiency: " << efficiency << std::endl;
    std::cout << "knapsack time: " << stoptime << std::endl;

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

//
// Try to "improve" the knapsack()d procmap ...
//

static
void
MinimizeCommCosts (std::vector<int>&        procmap,
                   const BoxArray&          ba,
                   const std::vector<long>& pts,
                   int                      nprocs)
{
    BL_PROFILE("MinimizeCommCosts()");

    BL_ASSERT(ba.size() == pts.size());
    BL_ASSERT(procmap.size() >= ba.size());

    if (nprocs < 2 || do_not_minimize_comm_costs) return;

    const Real strttime = ParallelDescriptor::second();
    //
    // Build a data structure that'll tell us who are our neighbors.
    //
    std::vector< std::vector<int> > nbrs(ba.size());
    //
    // Our "grow" factor; i.e. how far our tentacles grope for our neighbors.
    //
    const int Ngrow = 1;

    BoxArray grown(ba.size());

    for (int i = 0; i < ba.size(); i++)
        grown.set(i,BoxLib::grow(ba[i],Ngrow));

    for (int i = 0; i < grown.size(); i++)
    {
        std::list<int> li;

        std::vector< std::pair<int,Box> > isects = ba.intersections(grown[i]);

        for (int j = 0; j < isects.size(); j++)
            if (isects[j].first != i)
                li.push_back(isects[j].first);

        nbrs[i].resize(li.size());

        int k = 0;
        for (std::list<int>::const_iterator it = li.begin();
             it != li.end();
             ++it, ++k)
        {
            nbrs[i][k] = *it;
        }
    }

    if (verbose > 1 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "The neighbors list:\n";

        for (int i = 0; i < nbrs.size(); i++)
        {
            std::cout << i << "\t:";

            for (std::vector<int>::const_iterator it = nbrs[i].begin();
                 it != nbrs[i].end();
                 ++it)
            {
                std::cout << *it << ' ';
            }

            std::cout << "\n";
        }
    }
    //
    // Want lists of box IDs having the same size.
    //
    std::map< int,std::vector<int>,std::greater<int> > samesize;

    for (int i = 0; i < pts.size(); i++)
        samesize[pts[i]].push_back(i);

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

            std::cout << "\n";
        }
    }
    //
    // Build a data structure to maintain the latency count on a per-CPU basis.
    //
    std::vector<long> percpu(nprocs,0L);

    for (int i = 0; i < nbrs.size(); i++)
    {
        for (std::vector<int>::const_iterator it = nbrs[i].begin();
             it != nbrs[i].end();
             ++it)
        {
            if (procmap[i] != procmap[*it]) percpu[procmap[*it]]++;
        }
    }

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        long cnt = 0;
        for (int i = 0; i < percpu.size(); i++) cnt += percpu[i];
        std::cout << "Initial off-CPU connection count: " << cnt << '\n';
    }
    //
    // Originally I called SwapAndTest() until no links were changed.
    // This turned out to be very costly.  Next I tried calling it no
    // more than three times, or until no links were changed.  But after
    // testing a bunch of quite large meshes, it appears that the first
    // call gets "most" of the benefit of multiple calls.
    //
    SwapAndTest(samesize,nbrs,procmap,percpu);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        long cnt = 0;
        for (int i = 0; i < percpu.size(); i++) cnt += percpu[i];
        std::cout << "Final   off-CPU connection count: " << cnt << '\n';
    }

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        const Real stoptime = ParallelDescriptor::second() - strttime;

        std::cout << "MinimizeCommCosts() time: " << stoptime << '\n';
    }
}

int
main (int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

//    std::ifstream ifs("ba.60", std::ios::in);
//    std::ifstream ifs("ba.213", std::ios::in);
//    std::ifstream ifs("ba.5034", std::ios::in);
//    std::ifstream ifs("ba.15784", std::ios::in);
    std::ifstream ifs("ba.25600", std::ios::in);

    int nprocs = 512;

//    if (argc > 1) nprocs = atoi(argv[1]);

    std::cout << "nprocs = " << nprocs << std::endl;

    BoxArray ba;

    ba.readFrom(ifs);

//    ba.maxSize(32);

//    ba.writeOn(std::cout); std::cout << std::endl;

    std::cout << "ba.size() = " << ba.size() << std::endl;

    Array<int> procmap(ba.size());

    int N = 0;

    std::vector<long> pts(ba.size());

    for (unsigned int i = 0; i < pts.size(); i++)
        pts[i] = ba[i].numPts();

    std::vector< std::list<int> > vec = knapsack(pts,nprocs);

    std::list<int>::iterator lit;

    for (unsigned int i = 0; i < vec.size(); i++)
    {
        int where = (i + N) % nprocs;

        for (lit = vec[i].begin(); lit != vec[i].end(); ++lit)
            procmap[*lit] = where;
    }

    MinimizeCommCosts(procmap,ba,pts,nprocs);

    BoxLib::Finalize();
}
