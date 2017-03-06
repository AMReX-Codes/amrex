
#include <unistd.h>

#include <AMReX_TagBox.H>
#include <AMReX_MultiFab.H>
#include <ArrayView.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>
#include <TV_TempWrite.H>

#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::cerr;

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
    amrex::Initialize(argc,argv);
//
//  Parse the command line
//
    if (argc < 2)
        PrintUsage(argc,argv);

    ParmParse pp;
    
    std::string iFile; pp.get("iFile", iFile);
    int ascii = 0; pp.query("ascii",ascii);

//
//  Read multifab
//
    MultiFab mf;
    VisMF::Read(mf,iFile);
    
    int ngrow = mf.nGrow();
    pp.query("ngrow",ngrow);
    ngrow = std::min(ngrow,mf.nGrow());
    
    MultiFab tmp(mf.boxArray(),mf.nComp(),ngrow,Fab_allocate);
    MultiFab::Copy(tmp,mf,0,0,mf.nComp(),ngrow);
    if (ascii == 1)
    {
        for (MFIter mfi(tmp); mfi.isValid(); ++mfi)
        {
            cout << "FAB: " << mfi.index() << " on " << tmp[mfi].box() << endl;
            cout << tmp[mfi] << endl;
        }
        return true;
    }

    amrex::Finalize();
    
    return ArrayViewMultiFab(&tmp);
}
