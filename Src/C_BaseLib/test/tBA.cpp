#include <iostream>
#include <fstream>
#include <BoxArray.H>
#include <ParallelDescriptor.H>

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

    intersections_old(ba);
    intersections_new(ba);
}
