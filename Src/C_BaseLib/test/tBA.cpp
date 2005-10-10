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

int
main ()
{
//    std::ifstream ifs("ba.60", std::ios::in);
//    std::ifstream ifs("ba.213", std::ios::in);
    std::ifstream ifs("ba.1000", std::ios::in);
//    std::ifstream ifs("ba.5034", std::ios::in);

    BoxArray ba;

    ba.readFrom(ifs);

//    ba.writeOn(std::cout); std::cout << std::endl;

    Box bb = ba.minimalBox();

//    intersections_old(ba);
//    intersections_new(ba);

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
}
