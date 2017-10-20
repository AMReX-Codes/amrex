#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>

#include <new>
#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::cerr;

using std::set_new_handler;

static
void 
PrintUsage(int argc, char *argv[])
{
    cout << "Usage: " << endl;
    cout << argv[0] << " infile [options] \n\tOptions:" << endl;
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

//
//  Read multifab
//
    MultiFab mf;
    VisMF::Read(mf,iFile);

    Vector<int> comps;
    if (int nc = pp.countval("comps"))
    {
        comps.resize(nc);
        pp.getarr("comps",comps,0,nc);
    }
    else
    {
        int sComp = 0;
        pp.query("sComp",sComp);
        int nComp = mf.nComp();
        pp.query("nComp",nComp);
        BL_ASSERT(sComp+nComp < mf.nComp());
        comps.resize(nComp);
        for (int i=0; i<nComp; ++i)
            comps[i] = sComp + i;
    }

    
    int ngrow = mf.nGrow();
    pp.query("ngrow",ngrow);
    ngrow = std::min(ngrow,mf.nGrow());

    if (ParallelDescriptor::IOProcessor()) {
        for (int n=0; n<comps.size(); ++n)
            cout << "Comp: " << comps[n] << ", (min,max): " << mf.min(comps[n]) << ", " << mf.max(comps[n])  << endl;
    }
    
    amrex::Finalize();
}
