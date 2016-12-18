#include <iostream>
#include <fstream>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxDomain.H>
#include <AMReX_ParallelDescriptor.H>
#include <map>

using namespace amrex;

static
BoxArray
GetBndryCells_old (const BoxArray& ba,
                   int             ngrow,
                   const Box&      domain)
{
    std::cout << "GetBndryCells_old():" << std::endl;

    Real beg = ParallelDescriptor::second();

    const BoxList blgrids = BoxList(ba);

    BoxDomain bd;
    for (int i = 0; i < ba.size(); ++i)
    {
	BoxList gCells = amrex::boxDiff(amrex::grow(ba[i],ngrow),ba[i]);

	for (BoxList::iterator bli = gCells.begin(); bli != gCells.end(); ++bli)
	    bd.add(amrex::complementIn(*bli,blgrids));
    }

    BoxList bl;
    for (BoxDomain::const_iterator bdi = bd.begin(); bdi != bd.end(); ++bdi)
    {
        bl.push_back(*bdi);
    }

    std::cout << "    size before simplify() = " << bl.size() << std::endl;
    bl.simplify();
    Real end = ParallelDescriptor::second() - beg;
    std::cout << "    size after simplify() = " << bl.size() << std::endl;
    std::cout << "    time = " << end << std::endl;

    return BoxArray(bl);
}

static
BoxArray
GetBndryCells_new (const BoxArray& ba,
                   int             ngrow,
                   const Box&      domain)
{
    std::cout << "GetBndryCells_new():" << std::endl;

    Real beg = ParallelDescriptor::second();
    //
    // First get list of all ghost cells.
    //
    BoxList gcells;
    for (int i = 0; i < ba.size(); ++i)
    {
	gcells.join(amrex::boxDiff(amrex::grow(ba[i],ngrow),ba[i]));
    }
    std::cout << "    size of ghostcell list: " << gcells.size() << std::endl;
    //
    // Now strip out intersections with original BoxArray.
    //
    BoxList bcells;
    for (BoxList::const_iterator it = gcells.begin(); it != gcells.end(); ++it)
    {
        std::vector< std::pair<int,Box> > isects = ba.intersections(*it);

        if (isects.empty())
        {
            bcells.push_back(*it);
        }
        else
        {
            //
            // Collect all the intersection pieces.
            //
            BoxList pieces;
            for (int i = 0; i < isects.size(); i++)
                pieces.push_back(isects[i].second);
            BoxList leftover = amrex::complementIn(*it,pieces);
            bcells.catenate(leftover);
        }
    }

    std::cout << "    size of bndrycell list: " << bcells.size() << std::endl;
    //
    // Now strip out overlaps.
    //
    gcells.clear();

    gcells = amrex::removeOverlap(bcells);

    std::cout << "    size before simplify(): " << gcells.size() << std::endl;

    int min[BL_SPACEDIM] = {10000};
    int max[BL_SPACEDIM] = {0};

    std::map<int,int> dist;

    for (BoxList::const_iterator it = gcells.begin(); it != gcells.end(); ++it)
    {
        dist[it->numPts()]++;

        for (int i = 0; i < BL_SPACEDIM; i++)
        {
            int l = it->length(i);
            if (l > max[i]) max[i] = l;
            if (l < min[i]) min[i] = l;
        }
    }

    std::cout << "min,max length:\n";
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        std::cout << "i: " << i << ' ' << min[i] << ' ' << max[i] << '\n';
    }
    std::cout << std::endl;

    std::cout << "numPts() distribution before simplify():\n";
    for (std::map<int,int>::const_iterator it = dist.begin();
         it != dist.end();
         ++it)
    {
        std::cout << it->first << ' ' << it->second << '\n';
    }

    gcells.simplify();
    std::cout << "    size after simplify() = " << gcells.size() << std::endl;

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        min[i] = 10000; max[i] = 0;
    }
    dist.clear();

    for (BoxList::const_iterator it = gcells.begin(); it != gcells.end(); ++it)
    {
        dist[it->numPts()]++;

        for (int i = 0; i < BL_SPACEDIM; i++)
        {
            int l = it->length(i);
            if (l > max[i]) max[i] = l;
            if (l < min[i]) min[i] = l;
        }
    }

    std::cout << "min,max length:\n";
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        std::cout << "i: " << i << ' ' << min[i] << ' ' << max[i] << '\n';
    }
    std::cout << std::endl;

    std::cout << "numPts() distribution after simplify():\n";
    for (std::map<int,int>::const_iterator it = dist.begin();
         it != dist.end();
         ++it)
    {
        std::cout << it->first << ' ' << it->second << '\n';
    }

    gcells.maxSize(64);
    std::cout << "    size after maxSize() = " << gcells.size() << std::endl;
    Real end = ParallelDescriptor::second() - beg;
    std::cout << "    time = " << end << std::endl;

    dist.clear();
    for (BoxList::const_iterator it = gcells.begin(); it != gcells.end(); ++it)
    {
        dist[it->numPts()]++;
    }
    std::cout << "numPts() distribution after maxSize():\n";
    for (std::map<int,int>::const_iterator it = dist.begin();
         it != dist.end();
         ++it)
    {
        std::cout << it->first << ' ' << it->second << '\n';
    }

    
    return BoxArray(gcells);
}

static
void
intersections_old (const BoxArray& ba)
{
    const Real beg = ParallelDescriptor::second();
 
    int cnt = 0;

    const int ngrow = 1;

    Box isect;

    for (int j = 0; j < ba.size(); j++)
    {
        const Box& bx = amrex::grow(ba[j], ngrow);

        for (int i = 0; i < ba.size(); i++)
        {
            isect = bx & ba[i];

            if (isect.ok())
                cnt++;
        }
    }

    Real end = ParallelDescriptor::second() - beg;

    std::cout << "old cnt = " << cnt << ", time = " << end << std::endl;
}

static
void
intersections_new (const BoxArray& ba)
{
    const Real beg = ParallelDescriptor::second();

    int cnt = 0;

    const int ngrow = 1;
    
    for (int j = 0; j < ba.size(); j++)
    {
        std::vector< std::pair<int,Box> > v = ba.intersections(amrex::grow(ba[j], ngrow));

        cnt += v.size();
    }

    Real end = ParallelDescriptor::second() - beg;

    std::cout << "new cnt = " << cnt << ", time = " << end << std::endl;
}

static
BoxList
newComplementIn_old (const Box&     b,
                     const BoxList& bl)
{
    BoxList newb(b.ixType());
    newb.push_back(b);
    for (BoxList::const_iterator bli = bl.begin(); bli != bl.end() && newb.isNotEmpty(); ++bli)
    {
        for (BoxList::iterator newbli = newb.begin(); newbli != newb.end(); )
        {
            if (newbli->intersects(*bli))
            {
                BoxList tm = amrex::boxDiff(*newbli, *bli);
                newb.catenate(tm);
                newb.remove(newbli++);
            }
            else
            {
                ++newbli;
            }
        }
    }
    return newb;
}

static
void
Print (const BoxList& bl, const char* str)
{
    std::cout << str << ", size = " << bl.size() << " :\n";

    for (BoxList::const_iterator bli = bl.begin(); bli != bl.end(); ++bli)
    {
        std::cout << *bli << '\n';
    }
}

int
main ()
{
//    std::ifstream ifs("ba.60", std::ios::in);
//    std::ifstream ifs("ba.213", std::ios::in);
//    std::ifstream ifs("ba.1000", std::ios::in);
//    std::ifstream ifs("ba.5034", std::ios::in);
    std::ifstream ifs("ba.15456", std::ios::in);
//    std::ifstream ifs("ba.mac.294", std::ios::in);
//    std::ifstream ifs("ba.3865", std::ios::in);

    std::cout << "Got Here" << std::endl;

    BoxArray ba;

    ba.readFrom(ifs);

    std::cout << "Got Here 2" << std::endl;

//    ba.writeOn(std::cout); std::cout << std::endl;

    Box bb = ba.minimalBox();
    std::cout << "First Minimal box: " << bb << std::endl;

//    ba.refine(2);
//    ba.maxSize(32);
    
//    std::cout << "ba.size() = " << ba.size() << std::endl;
//    bb = ba.minimalBox();
//    std::cout << "Second Minimal box: " << bb << std::endl;

//    for (int i = 0; i < ba.size(); i++)
//        std::cout << ba[i] << '\n';

    if (ba.isDisjoint())
        std::cout << "The new BoxArray is disjoint" << std::endl;
    else
        std::cout << "The new BoxArray is NOT disjoint" << std::endl;

//    exit(0);

    bb.grow(4);

    BoxList bl;
    for (int i = 0; i < ba.size(); i++)
        bl.push_back(ba[i]);

//    bl.push_back(Box(IntVect(43,0,0),IntVect(52,31,31)));
//    bl.push_back(Box(IntVect(53,3,3),IntVect(70,28,28)));
//    bb = Box(IntVect(46,-2,-2),IntVect(51,4,33));

//    intersections_old(ba);
//    intersections_new(ba);

    BoxList bl1, bl2;

    {
    const Real beg = ParallelDescriptor::second();
    bl1 = amrex::complementIn(bb, bl);
    const Real end = ParallelDescriptor::second() - beg;
    std::cout << "complementIn(), size = " << bl1.size() << " time = " << end << std::endl;
    bl1.simplify();
    std::cout << "complementIn(), size after simplify() = " << bl1.size() << std::endl;
    bl1.minimize();
    std::cout << "complementIn(), size after minimize() = " << bl1.size() << std::endl;
    }

    {
    const Real beg = ParallelDescriptor::second();
    bl2 = newComplementIn_old(bb, bl);
    const Real end = ParallelDescriptor::second() - beg;
    std::cout << "newComplementIn_old(), size = " << bl2.size() << " time = " << end << std::endl;
    }

//    bl1.simplify();
//    bl2.simplify();

    BoxArray nba1(bl1);
    BoxArray nba2(bl2);

    std::cout << "nba1.numPts = " << nba1.numPts() << std::endl;
    std::cout << "nba2.numPts = " << nba2.numPts() << std::endl;

    if (nba1.contains(nba2) && nba2.contains(nba1))
    {
        std::cout << "nba1 & nba2 cover the same area" << std::endl;
    }
    else
    {
        std::cout << "nba1 & nba2 do NOT cover the same area" << std::endl;
        Print(bl1, "bl1");
        Print(bl2, "bl2");
    }

    exit(0);

    nba2 = GetBndryCells_new(ba, 1, bb);
    nba1 = GetBndryCells_old(ba, 1, bb);

    if (nba2.isDisjoint())
        std::cout << "The new BoxArray is disjoint" << std::endl;
    else
        std::cout << "The new BoxArray is NOT disjoint" << std::endl;
    
    if (nba1.contains(nba2) && nba2.contains(nba1))
        std::cout << "nba1 & nba2 cover the same area" << std::endl;
    else
        std::cout << "nba1 & nba2 do NOT cover the same area" << std::endl;
}
