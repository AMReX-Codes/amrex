#ifdef BL_ARCH_CRAY
#ifdef BL_USE_DOUBLE
#error "DOUBLE PRECISION NOT ALLOWED ON CRAY"
#endif
#endif

#ifndef        WIN32
#include <unistd.h>
#endif

#include "MultiFab.H"
#include "ArrayView.H"
#include "ParmParse.H"
#include "Utility.H"
#include "ParallelDescriptor.H"
#include "TV_TempWrite.H"

#ifdef BL_USE_NEW_HFILES
#include <new>
using std::setprecision;
#include <iostream>
#ifndef WIN32
using std::set_new_handler;
#endif
#else
#include <new.h>
#include <iostream.h>
#endif


static
void 
PrintUsage(int argc, char *argv[])
{
    cout << "Usage: " << endl;
    cout << argv[0] << endl;
    cout << "     mfab = MultiFab" << endl;
    exit(0);
}


int
main (int   argc,
      char* argv[])
{
    //
    // Make sure to catch new failures.
    //
#ifndef WIN32
    set_new_handler(Utility::OutOfMemory);
#endif

    ParallelDescriptor::StartParallel(&argc, &argv);

//
//  Parse the command line
//
    if (argc == 1)
        PrintUsage(argc, argv);

    ParmParse pp(argc-1,argv+1);

    if (pp.contains("help"))
        PrintUsage(argc, argv);

    aString iFile;

    pp.query("mfab", iFile);
    if (iFile.isNull())
        BoxLib::Abort("You must specify `mfab'");


//
//  Read multifab
//
    MultiFab mf;
    readMF(mf,iFile.c_str());
    
    return ArrayViewMultiFab(&mf);
}
