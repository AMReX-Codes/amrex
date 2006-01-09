#include <iostream>
#include <fstream>
#include <BoxArray.H>
#include <BoxDomain.H>
#include <ParallelDescriptor.H>

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
	BoxList gCells = BoxLib::boxDiff(BoxLib::grow(ba[i],ngrow),ba[i]);

	for (BoxList::iterator bli = gCells.begin(); bli != gCells.end(); ++bli)
	    bd.add(BoxLib::complementIn(*bli,blgrids));
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
	gcells.join(BoxLib::boxDiff(BoxLib::grow(ba[i],ngrow),ba[i]));
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
            BoxList leftover = BoxLib::complementIn(*it,pieces);
            bcells.catenate(leftover);
        }
    }

    std::cout << "    size of bndrycell list: " << bcells.size() << std::endl;
    //
    // Now strip out overlaps.
    //
    gcells.clear();

    gcells = BoxLib::removeOverlap(bcells);

//    std::cout << "    size before simplify(): " << gcells.size() << std::endl;
//    gcells.simplify();
    std::cout << "    size after simplify() = " << gcells.size() << std::endl;
    Real end = ParallelDescriptor::second() - beg;
    std::cout << "    time = " << end << std::endl;
    
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
        const Box bx = BoxLib::grow(ba[j], ngrow);

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
        std::vector< std::pair<int,Box> > v = ba.intersections(BoxLib::grow(ba[j], ngrow));

        cnt += v.size();
    }

    Real end = ParallelDescriptor::second() - beg;

    std::cout << "new cnt = " << cnt << ", time = " << end << std::endl;
}

static
BoxList
newComplementIn (const Box&     b,
		 const BoxList& bl)
{
    BoxList newb(b.ixType());

    Box minbox = bl.minimalBox();
    BoxList tmpbl = BoxLib::boxDiff(b,minbox);
    newb.catenate(tmpbl);

    BoxList mesh;
    mesh.push_back(minbox);
    mesh.maxSize(64);

    BoxArray ba(bl);

    for (BoxList::const_iterator bli = mesh.begin(); bli != mesh.end(); ++bli)
    {
        std::vector< std::pair<int,Box> > isects = ba.intersections(*bli);

        if (!isects.empty())
        {
            BoxList tmpbl;

            for (int i = 0; i < isects.size(); i++)
                tmpbl.push_back(isects[i].second);

            BoxList tm = BoxLib::complementIn(*bli,tmpbl);

            newb.catenate(tm);
        }
        else
        {
            newb.push_back(*bli);
        }
    }

    return newb;
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
                BoxList tm = BoxLib::boxDiff(*newbli, *bli);
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
    std::ifstream ifs("ba.15784", std::ios::in);

    BoxArray ba;

    ba.readFrom(ifs);

//    ba.writeOn(std::cout); std::cout << std::endl;

    Box bb = ba.minimalBox();

    bb.grow(2);

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
    bl1 = BoxLib::complementIn(bb, bl);
    const Real end = ParallelDescriptor::second() - beg;
    std::cout << "complementIn(), size = " << bl1.size() << " time = " << end << std::endl;
    bl1.simplify();
    std::cout << "complementIn(), size after simplify() = " << bl1.size() << std::endl;
    bl1.minimize();
    std::cout << "complementIn(), size after minimize() = " << bl1.size() << std::endl;
    }

#if 0
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
        std::cout << "nba1 & nba2 cover the same area" << std::endl;
    else
    {
        std::cout << "nba1 & nba2 do NOT cover the same area" << std::endl;

        Print(bl1, "bl1");
        Print(bl2, "bl2");
    }

    BoxArray nba2 = GetBndryCells_new(ba, 1, bb);
    BoxArray nba1 = GetBndryCells_old(ba, 1, bb);

    if (nba2.isDisjoint())
        std::cout << "The new BoxArray is disjoint" << std::endl;
    else
        std::cout << "The new BoxArray is NOT disjoint" << std::endl;
    
    if (nba1.contains(nba2) && nba2.contains(nba1))
        std::cout << "nba1 & nba2 cover the same area" << std::endl;
    else
        std::cout << "nba1 & nba2 do NOT cover the same area" << std::endl;
#endif
}
