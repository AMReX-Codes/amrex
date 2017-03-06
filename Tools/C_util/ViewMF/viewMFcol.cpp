
#include <unistd.h>

#include <AMReX_MultiFab.H>
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
//
//  Read multifab
//
    MultiFab mf;
    VisMF::Read(mf,iFile);

    
    int sComp = 0; pp.query("sComp",sComp); sComp=std::min(sComp,mf.nComp()-1);
    int nComp = mf.nComp(); pp.query("nComp",nComp);
    if (nComp==-1) nComp = mf.nComp();
    nComp = std::min(mf.nComp()-sComp,nComp);
    
    Box domain = mf.boxArray().minimalBox();
    int dir = 1; pp.query("dir",dir);

    int idx = (int) (0.5*(domain.smallEnd()[dir] + domain.bigEnd()[dir])); pp.query("idx",idx);
    
    for (int d=0; d<BL_SPACEDIM; ++d) {
        if (d!=dir) {            
            domain.setSmall(dir,idx);
            domain.setBig(dir,idx);
        }
    }

    FArrayBox fab(domain,nComp);
    mf.copy(fab,sComp,0,nComp);

    cout << "Components: " << sComp << " : " << sComp + nComp - 1 << endl;
    cout << fab << endl;

    amrex::Finalize();
    
    return 0;
}
