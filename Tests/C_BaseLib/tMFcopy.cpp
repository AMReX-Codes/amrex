#include <iostream>
#include <fstream>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    //
    // Use Space Filling Curve algorithm for distributing grids.
    //
    DistributionMapping::strategy(DistributionMapping::SFC);
    //
    // If there are >= "2" grids per CPU on average use Space Filling
    // Curve algorithm for distributing grids, otherwise use KnapSack
    // algorithm.
    //
    DistributionMapping::SFC_Threshold(2);

    const char* file = "ba.23925";

    std::ifstream ifs(file, std::ios::in);

    if (!ifs.good())
        amrex::Error("Unable to open file");

    BoxArray fba;

    fba.readFrom(ifs);

    if (!ifs.good())
        amrex::Error("Read of BoxArray failed");

    if (!fba.ok())
        amrex::Error("BoxArray is not OK");

    if (!fba.isDisjoint())
        amrex::Error("BoxArray is not disjoint");

    fba.maxSize(32);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "number of grids in fba: " << fba.size() << '\n';
    //
    // Let do a simple example of copy()ing from one MultiFab to another
    // that it covers. We'll build a BoxArray that contains every other Box
    // in fba.
    //
    BoxList cbl;

    for (int i = 0; i < fba.size(); i++)
        if (i%2 == 0)
            cbl.push_back(fba[i]);

    cbl.simplify();

    cbl.maxSize(32);

    BoxArray cba(cbl);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "number of grids in cba: " << cba.size() << std::endl;
    //
    // If you want to make the copy do more work increase NComp.
    //
    const int NComp = 1;

    DistributionMapping fdm{fba};
    DistributionMapping cdm{cba};

    MultiFab fmf(fba, fdm, NComp, 0);
    MultiFab cmf(cba, cdm, NComp, 0);

    fmf.setVal(1.23e45);

    cmf.copy(fmf);

    if (cdm[0] == ParallelDescriptor::MyProc())
        std::cout << cmf[0] << std::endl;

    amrex::Finalize();

    return 0;
}
