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

#include <new>
using std::setprecision;
#include <iostream>
#ifndef WIN32
using std::set_new_handler;
#endif

static
void 
PrintUsage(int argc, char *argv[])
{
    cout << "Usage: " << endl;
    cout << argv[0] << " infile [options] \n\tOptions:" << endl;
    cout << "\t   -ascii   \t[if set, dump ascii mfab to stdout]" << endl;
    cout << "\t   ngrow=<#>\t[number of grow cells to include in output]" << endl;
    cout << endl;
    exit(0);
}


int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);
//
//  Parse the command line
//
    if (argc < 2)
        PrintUsage(argc,argv);

    if (argv[1][0] == '-')
    {
        cerr << "input file must be first argument\n";
        PrintUsage(argc, argv);
    }

    ParmParse pp;
    
    if (pp.contains("help"))
        PrintUsage(argc, argv);
    
    std::string iFile = argv[1];
    
    bool ascii = false;
    if (pp.contains("ascii"))
        ascii = true;
//
//  Read multifab
//
    MultiFab mf;
    readMF(mf,iFile.c_str());
    
    int ngrow = mf.nGrow();
    pp.query("ngrow",ngrow);
    ngrow = std::min(ngrow,mf.nGrow());
    
    MultiFab tmp(mf.boxArray(),mf.nComp(),ngrow,Fab_allocate);
    MultiFab::Copy(tmp,mf,0,0,mf.nComp(),ngrow);
    if (ascii)
    {
        for (MFIter mfi(tmp); mfi.isValid(); ++mfi)
        {
            cout << "FAB: " << mfi.index() << endl;
            cout << tmp[mfi] << endl;
        }
        return true;
    }

    BoxLib::Finalize();
    
    return ArrayViewMultiFab(&tmp);
}
