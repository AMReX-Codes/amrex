#include <iostream>
#include <fstream>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxDomain.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DistributionMapping.H>

using namespace amrex;

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
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

//    std::ifstream ifs("ba.60", std::ios::in);
//    std::ifstream ifs("ba.213", std::ios::in);
//    std::ifstream ifs("ba.1000", std::ios::in);
//    std::ifstream ifs("ba.5034", std::ios::in);
//    std::ifstream ifs("ba.15456", std::ios::in);
//    std::ifstream ifs("ba.mac.294", std::ios::in);
//    std::ifstream ifs("ba.3865", std::ios::in);
    std::ifstream ifs("ba.23925", std::ios::in);

    BoxArray ba;

    ba.readFrom(ifs);

    std::cout << "# of grids: " << ba.size() << '\n';

    for (int nprocs = 2; nprocs < 5000; nprocs *= 2)
    {
        std::cout << "\nnprocs = " << nprocs << '\n';

        DistributionMapping::strategy(DistributionMapping::SFC);
        DistributionMapping dm1(ba,nprocs);
        DistributionMapping::FlushCache();

        std::cout << '\n';

        DistributionMapping::strategy(DistributionMapping::KNAPSACK);
        DistributionMapping dm2(ba,nprocs);
        DistributionMapping::FlushCache();
    }

    amrex::Finalize();
}
