
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#include <AMReX_Box.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_VisMF.H>

static
void
PrintUsage (const char* progName)
{
    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "This utility reads a fixed grids file and coarsens the\n"
             << "grids by a specified amount.  If there is more than one\n"
             << "level of refinement, then different coarsening ratios\n"
             << "are input for each level.  This is useful in the process\n"
             << "of converting from output grid format to the fixed_grids\n"
             << "format.\n";
        std::cout << '\n';
        std::cout << "Usage:" << '\n';
        std::cout << progName << '\n';
        std::cout << "    infile  = inputFileName" << '\n';
        std::cout << "  crsratio  = Coarsening ratio for each level" << '\n';
        std::cout << "   [-help]" << '\n';
        std::cout << '\n';
    }
    exit(1);
}

int
main (int   argc,
      char* argv[])
{
    if (argc == 1)
        PrintUsage(argv[0]);

    amrex::Initialize(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    std::string iFile;

    pp.query("infile", iFile);                             // Input File
    if (iFile.empty() && ParallelDescriptor::IOProcessor())
        amrex::Abort("You must specify `infile'");

    int nCrsRatio = pp.countval("crsratio");
    if (nCrsRatio == 0)
        amrex::Abort("You must specify `crsratio'");

    Vector<int> crsRatio(nCrsRatio);
    for (int n = 0; n < nCrsRatio; n++)
        pp.get("crsratio", crsRatio[n], n);


    //------------------------------------//
    // Open the Grid File and Read Header //
    //------------------------------------//

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    std::ifstream is;
    is.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    is.open(iFile.c_str(), std::ios::in);
    if (is.fail())
        amrex::FileOpenFailed(iFile);

    int nRefLevels;
    is >> nRefLevels;

    if (nCrsRatio != nRefLevels)
        amrex::Abort("nCrsRatio != nRefLevels");

    std::cout << nRefLevels << std::endl;
        

    //----------------------------------------------------//
    // Loop Through Refined Levels and Generate MultiFabs //
    //----------------------------------------------------//

    int nGrds;
    for (int lev = 0; lev < nRefLevels; lev++)
    {
        is >> nGrds;
        std::cout << nGrds << std::endl;

        Box inBox;
        for (int bx=0; bx < nGrds; bx++)
        {
            is >> inBox;
            inBox.coarsen(crsRatio[lev]);
            std::cout << inBox << std::endl;
        }
    }

    amrex::Finalize();

    return 0;
}
