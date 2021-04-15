
#include <new>
#include <iostream>
#include <cstdlib>
#include <cstring>
using std::ios;
using std::set_new_handler;

#include <unistd.h>

#include <AMReX_Box.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_VisMF.H>


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
        cout << "    lodims  = lower bound of the intersecting box" << '\n';
        cout << "    hidims  = upper bound of the intersecting box" << '\n';
        cout << "  refratio  = Coarsening ratio for each level" << '\n';
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
        amrex::Abort("You must specify `infile'");

    Vector<int> loDims(BL_SPACEDIM);
    if (pp.countval("lodims") != BL_SPACEDIM)
        amrex::Abort("You must specify BL_SPACEDIM `lodims'");
    for (int n = 0; n < BL_SPACEDIM; n++)
        pp.get("lodims", loDims[n], n);

    Vector<int> hiDims(BL_SPACEDIM);
    if (pp.countval("hidims") != BL_SPACEDIM)
        amrex::Abort("You must specify BL_SPACEDIM `hidims'");
    for (int n = 0; n < BL_SPACEDIM; n++)
        pp.get("hidims", hiDims[n], n);

    IntVect loIV(AMREX_D_DECL(loDims[0],loDims[1],loDims[2]));
    IntVect hiIV(AMREX_D_DECL(hiDims[0],hiDims[1],hiDims[2]));
    Box intersectBox(loIV, hiIV, IndexType::TheCellType());

    int nCrsRatio = pp.countval("refratio");
    if (nCrsRatio == 0)
        amrex::Abort("You must specify `refratio'");

    Vector<int> refRatio(nCrsRatio);
    for (int n = 0; n < nCrsRatio; n++)
        pp.get("refratio", refRatio[n], n);


    //------------------------------------//
    // Open the Grid File and Read Header //
    //------------------------------------//

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    ifstream is;
    is.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.length());
    is.open(iFile.c_str(), ios::in);
    if (is.fail())
        Utility::FileOpenFailed(iFile);

    int nRefLevels;
    is >> nRefLevels;

    if (nCrsRatio != nRefLevels)
        amrex::Abort("nCrsRatio != nRefLevels");

    cout << nRefLevels << endl;


    //----------------------------------------------------//
    // Loop Through Refined Levels and Generate MultiFabs //
    //----------------------------------------------------//

    int nGrds;
    for (int lev = 0; lev < nRefLevels; lev++)
    {
        is >> nGrds;

        Vector<Box> oBoxes(nGrds);

        int nOBoxes = 0;
        Box inBox;
        for (int bx=0; bx < nGrds; bx++)
        {
            is >> inBox;

            if (intersectBox.intersects(inBox))
            {
                oBoxes.set(nOBoxes, inBox);
                nOBoxes++;
            }
        }

        cout << nOBoxes << endl;
        for (int n = 0; n < nOBoxes; ++n)
            cout << oBoxes[n] << endl;

        intersectBox = intersectBox.refine(refRatio[lev]);
    }

    //
    // Close the Input File
    //
    if (ParallelDescriptor::IOProcessor())
        is.close();

}
