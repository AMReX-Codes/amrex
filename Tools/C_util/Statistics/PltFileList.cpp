
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
    std::string hdrFile = ("rdm_0_plt00");
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
    Array<MultiFab*> mean;
    Array<MultiFab*> variance;

    int finestLevel;
    int nmax = 10;
    Real nmaxinv = 1.0/double(nmax);
    Real nmaxinvm1 = 1.0/double(nmax-1);
    std::cout << nmaxinv << std::endl;
    for (int i = 0; i < nmax; i++)
    {

      char idx[4];
      sprintf(idx,"%i",i*10);
          
      std::string idxs(idx);
      //File = hdrFile + idxs+ iFile;
      File = hdrFile + idxs;
      if (i == 0) File +=idxs;
      std::cout << File << std::endl;
      
      DataServices dataServices(File, fileType);

      if (!dataServices.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);

      AmrData& amrData = dataServices.AmrDataRef();
      int nComp = amrData.NComp();

      if (i == 0) {

	finestLevel = amrData.FinestLevel();

	mean.resize(finestLevel+1);
	variance.resize(finestLevel+1);

	for (int iLevel=0; iLevel<=finestLevel; ++iLevel) {
	  const BoxArray& crseBA = amrData.boxArray(iLevel);
	  mean[iLevel] = new MultiFab(crseBA,nComp,0);
	  variance[iLevel] = new MultiFab(crseBA,nComp,0);
	  mean[iLevel]->setVal(0.0);
	  variance[iLevel]->setVal(0.0);
	}
      }
	
      ComputeAmrDataList(amrData, mean, variance, 0, 1);
      
    }

    File = hdrFile + "00";
    DataServices dataServices(File, fileType);
    AmrData& amrData = dataServices.AmrDataRef();
    int nComp = amrData.NComp();

    for (int iLevel=0; iLevel<=finestLevel; ++iLevel) {

      const BoxArray& crseBA = amrData.boxArray(iLevel);
      MultiFab tmp(crseBA,nComp,0);
      
      mean[iLevel]->mult(nmaxinv,0);

      for (MFIter mfi(*mean[iLevel]);mfi.isValid();++mfi) {
	FArrayBox& fab = (*mean[iLevel])[mfi];
	const Box& fabbox = mfi.validbox();

	FArrayBox musquare(fabbox,nComp);
	
	musquare.copy(fab,0,0,nComp);
	musquare.mult(fab,0,0,nComp);
	musquare.mult(double(nmax));
	(*variance[iLevel])[mfi].minus(musquare,0,0,nComp);
	(*variance[iLevel])[mfi].mult(nmaxinvm1);
      }
      
    }

    // output solns to plot files
    File = outfile + "_mean";
    WritePlotFile(mean,amrData,File,verbose);

    File = outfile + "_var";
    WritePlotFile(variance,amrData,File,verbose);
    
    //
    // This calls ParallelDescriptor::EndParallel() and exit()
    //
    //DataServices::Dispatch(DataServices::ExitRequest, NULL);
    BoxLib::Finalize();
}


