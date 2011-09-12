#include <iostream>
#include <fstream>
#include <BoxArray.H>
#include <MultiFab.H>
#include <ParallelDescriptor.H>

static
void
copy_test_one (const BoxArray& fba)
{
    //
    // In this test we'll copy from 
    //
}

int
main (int argc, char* argv[])
{
    BoxLib::Initialize(argc, argv);
    //
    // Use Space Filling Curve algorithm for distributing grids.
    //
    DistributionMapping::strategy(DistributionMapping::SFC);
    //
    // If there are >= "2" grids per CPU on average use Space Filling
    // Curse algorithm for distributing grids, otherwise use KnapSack
    // algorithm.
    //
    DistributionMapping::SFC_Threshold(2);

    const char* file = "ba.23925";

    std::ifstream ifs(file, std::ios::in);

    if (!ifs.good())
        BoxLib::Error("Unable to open file");

    BoxArray fba;

    fba.readFrom(ifs);

    if (!ifs.good())
        BoxLib::Error("Read of BoxArray failed");

    if (!fba.ok())
        BoxLib::Error("BoxArray is not OK");

    if (!fba.isDisjoint())
        BoxLib::Error("BoxArray is not disjoint");

    if (ParallelDescriptor::IOProcessor())
        std::cout << "number of grids in fba: " << fba.size() << '\n';
    //
    // Let do a simple example of copy()ing from one MultiFab to another
    // that cover the same area, but that have a different processor distribution.
    //
    BoxList cbl = fba.boxList();

    cbl.coarsen(2);

    cbl.simplify();

    cbl.maxSize(64);

    BoxArray cba(cbl);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "number of grids in cba: " << cba.size() << '\n';
    //
    // If you want to make the copy do more work increase NComp.
    //
    const int NComp = 1;

    MultiFab fmf(fba, NComp, 0);
    MultiFab cmf(cba, NComp, 0);

    fmf.setVal(1.23e45);

    cmf.copy(fmf);

    BoxLib::Finalize();

    return 0;
}
