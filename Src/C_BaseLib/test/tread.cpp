
//
// $Id: tread.cpp,v 1.1 2001-10-11 16:33:38 lijewski Exp $
//
// Read in and sum up FABs on command line.
//
// This assumes they're all the same size and shape.
//

#include <fstream>
#include <iostream>

#include <Utility.H>
#include <FArrayBox.H>

int
main (int argc, char** argv)
{
    BoxLib::Initialize(argc, argv);

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
        BoxLib::FileOpenFailed(argv[idx]);

    std::cout << "Reading " << argv[idx] << " ..." << std::endl;

    sum.readFrom(ifs);

    idx++;

    for ( ; idx < argc; idx++)
    {


        FArrayBox tmp;

        std::ifstream ifs;

        ifs.open(argv[idx], std::ios::in);

        if (!ifs.good())
            BoxLib::FileOpenFailed(argv[idx]);

        std::cout << "Reading " << argv[idx] << " ..." << std::endl;

        tmp.readFrom(ifs);

        sum += tmp;
    }

    std::ofstream ofs;

    ofs.open("SUM",std::ios::out|std::ios::trunc);

    if (!ofs.good()) BoxLib::FileOpenFailed("SUM");

    sum.writeOn(ofs);

    BoxLib::Finalize();

    return 0;
}
