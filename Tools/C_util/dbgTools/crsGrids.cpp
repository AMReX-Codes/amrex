//BL_COPYRIGHT_NOTICE

//
// $ $
//

#ifdef BL_USE_NEW_HFILES
#include <new>
#include <iostream>
#include <cstdlib>
#include <cstring>
using std::ios;
using std::set_new_handler;
#else
#include <new.h>
#include <iostream.h>
#include <stdlib.h>
#include <string.h>
#endif

#include <unistd.h>

#include "Box.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "Utility.H"
#include "ArrayLim.H"
#include "VisMF.H"

//
// This MUST be defined if don't have pubsetbuf() in I/O Streams Library.
//
#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif



// ============================================================
// ===== Routine:  void PrintUsage (const char* progName) =====
// ============================================================

static
void
PrintUsage (const char* progName)
{
    if (ParallelDescriptor::IOProcessor())
    {
        cout << "This utility reads a fixed grids file and coarsens the\n"
             << "grids by a specified amount.  If there is more than one\n"
             << "level of refinement, then different coarsening ratios\n"
             << "are input for each level.  This is useful in the process\n"
             << "of converting from output grid format to the fixed_grids\n"
             << "format.\n";
        cout << '\n';
        cout << "Usage:" << '\n';
        cout << progName << '\n';
        cout << "    infile  = inputFileName" << '\n';
        cout << "  crsratio  = Coarsening ratio for each level" << '\n';
        cout << "   [-help]" << '\n';
        cout << '\n';
    }
    exit(1);
}



// =======================================================
// ===== Routine:  int main (int argc, char* argv[]) =====
// =======================================================

int
main (int   argc,
      char* argv[])
{
    //-------------//
    // Print Usage //
    //-------------//

    if (argc == 1)
        PrintUsage(argv[0]);


    //----------------------------------//
    // Make sure to catch new failures. //
    //----------------------------------//

    set_new_handler(Utility::OutOfMemory);


    //-----------------------------------//
    // Initialize and Parse Command Line //
    //-----------------------------------//

    // Parallel Initialization

    ParallelDescriptor::StartParallel(&argc, &argv);


    // Parse Command Line

    ParmParse pp(argc-1,argv+1);

    if (pp.contains("help"))                          // Print usage
        PrintUsage(argv[0]);


    // Scan the arguments.

    aString iFile;

    pp.query("infile", iFile);                             // Input File
    if (iFile.isNull() && ParallelDescriptor::IOProcessor())
        BoxLib::Abort("You must specify `infile'");

    int nCrsRatio = pp.countval("crsratio");
    if (nCrsRatio == 0)
        BoxLib::Abort("You must specify `crsratio'");

    Array<int> crsRatio(nCrsRatio);
    for (int n = 0; n < nCrsRatio; n++)
        pp.get("crsratio", crsRatio[n], n);


    //------------------------------------//
    // Open the Grid File and Read Header //
    //------------------------------------//

    ifstream is;
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    is.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.length());
    is.open(iFile.c_str(), ios::in);
    if (is.fail())
        Utility::FileOpenFailed(iFile);

    int nRefLevels;
    is >> nRefLevels;

    if (nCrsRatio != nRefLevels)
        BoxLib::Abort("nCrsRatio != nRefLevels");

    cout << nRefLevels << endl;
        

    //----------------------------------------------------//
    // Loop Through Refined Levels and Generate MultiFabs //
    //----------------------------------------------------//

    int nGrds;
    for (int lev = 0; lev < nRefLevels; lev++)
    {
        is >> nGrds;
        cout << nGrds << endl;

        Box inBox;
        for (int bx=0; bx < nGrds; bx++)
        {
            is >> inBox;
            inBox.coarsen(crsRatio[lev]);
            cout << inBox << endl;
        }
    }

    //
    // Close the Input File
    //
    if (ParallelDescriptor::IOProcessor())
        is.close();

}
