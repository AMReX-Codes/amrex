#include <iostream>
#include <fstream>
#include <BoxArray.H>
#include <BoxDomain.H>
#include <ParallelDescriptor.H>

static
BoxArray
GetBndryCells(const BoxArray& ba,
              int             ngrow,
              const Box&      domain)
{
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

    Real end = ParallelDescriptor::second() - beg;

    std::cout << "BoxArray(bl).size() = " << bl.size() << " for ngrow = " << ngrow << std::endl;

    std::cout << "GetBndryCells() time = " << end << std::endl;

    beg = ParallelDescriptor::second();

    bl.simplify();

    end = ParallelDescriptor::second() - beg;

    std::cout << "GetBndryCells() simplify() time = " << end << std::endl;

    return BoxArray(bl);
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
    std::ifstream ifs("ba.5034", std::ios::in);

    BoxArray ba;

    ba.readFrom(ifs);

//    ba.writeOn(std::cout); std::cout << std::endl;

    Box bb = ba.minimalBox();

    intersections_old(ba);
    intersections_new(ba);

    BoxArray nba;
    nba = GetBndryCells(ba, 1, bb);
    nba = GetBndryCells(ba, 2, bb);
    nba = GetBndryCells(ba, 3, bb);
}
