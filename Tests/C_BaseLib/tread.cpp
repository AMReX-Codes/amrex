//
// Read in and sum up FABs on command line.
//
// This assumes they're all the same size and shape.
//

#include <fstream>
#include <iostream>

#include <AMReX_Utility.H>
#include <AMReX_FArrayBox.H>

using namespace amrex;

int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv);

    if (argc < 2)
    {
        std::cerr << "tread: fab1 fab2 fab3 ..." << std::endl;
        exit(1);
    }

    int idx = 1;

    FArrayBox sum;

    std::ifstream ifs;

    ifs.open(argv[idx], std::ios::in);

    if (!ifs.good())
        amrex::FileOpenFailed(argv[idx]);

    std::cout << "Reading " << argv[idx] << " ..." << std::endl;

    sum.readFrom(ifs);

    idx++;

    for ( ; idx < argc; idx++)
    {


        FArrayBox tmp;

        std::ifstream ifs;

        ifs.open(argv[idx], std::ios::in);

        if (!ifs.good())
            amrex::FileOpenFailed(argv[idx]);

        std::cout << "Reading " << argv[idx] << " ..." << std::endl;

        tmp.readFrom(ifs);

        sum += tmp;
    }

    std::ofstream ofs;

    ofs.open("SUM",std::ios::out|std::ios::trunc);

    if (!ofs.good()) amrex::FileOpenFailed("SUM");

    sum.writeOn(ofs);

    amrex::Finalize();

    return 0;
}
