
#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
using std::ios;
using std::set_new_handler;

#include <unistd.h>
#include "WritePlotFile.H"
#include "ComputeAmrDataStat.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"

#ifndef NDEBUG
#include "TV_TempWrite.H"
#endif

static
void
PrintUsage (const char* progName)
{
    std::cout << '\n';
    std::cout << "This routine reads a pltfile and determines its statistics" << std::endl
	      << std::endl;
    std::cout << "Usage:" << '\n';
    std::cout << progName << '\n';
    std::cout << "   infile=inputFileName" << '\n';
    std::cout << "   [outfile=outputFileName]" << '\n';
    std::cout << "   [-help]" << '\n';
    std::cout << "   [-verbose]" << '\n';
    std::cout << '\n';
    std::cout << " Note: outfile required if verbose used" << '\n';
    exit(1);
}

int
main (int   argc,
      char* argv[])
{
    if (argc == 1)
        PrintUsage(argv[0]);
    
    BoxLib::Initialize(argc,argv);
    ParmParse pp;

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    std::string iFile;
    std::string outfile;
    std::string hdrFile = "rdm_";//("rdm_0_plt00");
    std::string File;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }

    pp.query("infile", iFile);
    if (iFile.empty())
        BoxLib::Abort("You must specify `infile'");

    pp.query("outfile", outfile);
    if (outfile.empty())
        BoxLib::Abort("You must specify `outfile'");

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    MultiFab mean;
    MultiFab variance;

    int finestLevel;
    int nmax = 100;
    Real nmaxinv = 1.0/double(nmax);
    Real nmaxinvm1 = 1.0/double(nmax-1);
    std::cout << nmaxinv << std::endl;
    int nComp;

    for (int i = 0; i < nmax; i++)
    {

      char idx[4];
      sprintf(idx,"%i",i);
          
      std::string idxs(idx);
      File = hdrFile + idxs + iFile;
      //File = hdrFile + idxs;
      //if (i == 0) File +=idxs;
      std::cout << File << std::endl;
      
      DataServices dataServices(File, fileType);

      if (!dataServices.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);

      AmrData& amrData = dataServices.AmrDataRef();
      
      nComp = amrData.NComp();
      Array<string> names = amrData.PlotVarNames();
      Array<int> destcomp(names.size());

      for (int j=0; j<names.size();j++) 
	destcomp[j] = j;

      if (i == 0) {

	finestLevel = amrData.FinestLevel();
	Box tmpbox(amrData.ProbDomain()[finestLevel]);
	BoxArray ba(tmpbox);
	mean.define(ba,nComp,0,Fab_allocate);
	variance.define(ba,nComp,0,Fab_allocate);
	mean.setVal(0);
	variance.setVal(0);

      }
    	
      MultiFab tmpmean(mean.boxArray(),nComp,0);
      MultiFab tmpvar(mean.boxArray(),nComp,0);
      amrData.FillVar(tmpmean,finestLevel,names,destcomp);
      amrData.FillVar(tmpvar,finestLevel,names,destcomp);

      MultiFab::Add(mean,tmpmean,0,0,nComp,0);

      tmpvar[0].mult(tmpvar[0]);
      MultiFab::Add(variance,tmpvar,0,0,nComp,0);
    }


    mean.mult(nmaxinv);
    MultiFab tmpmean(mean.boxArray(),nComp,0);
    tmpmean.copy(mean);
    tmpmean[0].mult(tmpmean[0]);
    tmpmean.mult(double(nmax));
    variance.minus(tmpmean,0,nComp,0);
    variance.mult(nmaxinvm1);

    VisMF::Write(mean,"mean");
    VisMF::Write(variance,"variance");

    BoxLib::Finalize();
}


