
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <unistd.h>

#include <WritePlotFile.H>
#include <ComputeAmrDataStat.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <AMReX_Utility.H>

#ifdef AMREX_DEBUG
#include <TV_TempWrite.H>
#endif

using std::ios;

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
    
    amrex::Initialize(argc,argv);
    ParmParse pp;

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    std::string iFile;
    std::string outfile;
    std::string hdrFile = "sol";//("rdm_0_plt00");
    std::string File;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }

    pp.query("infile", iFile);
    if (iFile.empty())
        amrex::Abort("You must specify `infile'");

    pp.query("outfile", outfile);
    if (outfile.empty())
        amrex::Abort("You must specify `outfile'");

    int iFile_type = 0;
    pp.query("infile_type", iFile_type);
  
    int nmax = 100;
    pp.query("nmax", nmax);



    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    MultiFab mean;
    MultiFab variance;

    int finestLevel;
    Real nmaxinv = 1.0/double(nmax);
    Real nmaxinvm1 = 1.0/double(nmax-1);
    int nComp;
    bool first = true;


    if (iFile_type == 0) 
    {
      for (int i = 1; i <= nmax; i++)
      {
        File  = amrex::Concatenate(hdrFile, i, 1);
        File += '/';
        File += iFile;
      
	DataServices dataServices(File, fileType);
	
	if (!dataServices.AmrDataOk())
	  //
	  // This calls ParallelDescriptor::EndParallel() and exit()
	  //
	  DataServices::Dispatch(DataServices::ExitRequest, NULL);

	AmrData& amrData = dataServices.AmrDataRef();
      
	nComp = amrData.NComp();
	Vector<string> names = amrData.PlotVarNames();
	Vector<int> destcomp(names.size());

	for (int j=0; j<names.size();j++) 
	  destcomp[j] = j;

	if (first) {
	  
	  finestLevel = amrData.FinestLevel();
	  Box tmpbox(amrData.ProbDomain()[finestLevel]);
	  BoxArray ba(tmpbox);
	  ba.maxSize(128);
	  mean.define(ba,nComp,0,Fab_allocate);
	  variance.define(ba,nComp,0,Fab_allocate);
	  mean.setVal(0);
	  variance.setVal(0);
	  first = false;
	}
    	
	MultiFab tmpmean(mean.boxArray(),nComp,0);
	MultiFab tmpvar(mean.boxArray(),nComp,0);
	amrData.FillVar(tmpmean,finestLevel,names,destcomp);
	amrData.FillVar(tmpvar,finestLevel,names,destcomp);
      
	MultiFab::Add(mean,tmpmean,0,0,nComp,0);
	for (MFIter mfi(tmpvar); mfi.isValid(); ++mfi)
	{ 
	  const int idx = mfi.index();
	  tmpvar[idx].mult(tmpvar[idx]);
	}
	MultiFab::Add(variance,tmpvar,0,0,nComp,0);
      }
    }
    else if (iFile_type == 1) 
    {
      for (int i = 1; i <= nmax; i++)
      {
        File  = amrex::Concatenate(hdrFile, i, 1);
        File += '/';
        File += iFile;

	MultiFab tmpmean;
	VisMF::Read(tmpmean,File);
	MultiFab tmpvar(tmpmean.boxArray(),tmpmean.nComp(),0);
	tmpvar.copy(tmpmean);

	if (first) {
	  nComp = tmpmean.nComp();
	  mean.define(tmpmean.boxArray(),nComp,0,Fab_allocate);
	  variance.define(tmpmean.boxArray(),nComp,0,Fab_allocate);
	  mean.setVal(0.);
	  variance.setVal(0.);
	  first = false;
	}
      
	MultiFab::Add(mean,tmpmean,0,0,nComp,0);
	for (MFIter mfi(tmpvar); mfi.isValid(); ++mfi)
	{ 
	  const int idx = mfi.index();
	  tmpvar[idx].mult(tmpvar[idx]);
	}
	MultiFab::Add(variance,tmpvar,0,0,nComp,0);
      }

    }


    mean.mult(nmaxinv);
    MultiFab tmpmean(mean.boxArray(),nComp,0);
    tmpmean.copy(mean);

    for (MFIter mfi(tmpmean); mfi.isValid(); ++mfi)
    { 
      const int idx = mfi.index();
      tmpmean[idx].mult(tmpmean[idx]);
    }
    tmpmean.mult(double(nmax));
    variance.minus(tmpmean,0,nComp,0);
    variance.mult(nmaxinvm1);

    // The I/O processor makes the directory if it doesn't already exist.
    if (ParallelDescriptor::IOProcessor())
      if (!amrex::UtilCreateDirectory(outfile, 0755))
	amrex::CreateDirectoryFailed(outfile);
    ParallelDescriptor::Barrier();
    std::string mfile = outfile + "/" + iFile + "_mean";
    std::string vfile = outfile + "/" + iFile + "_var";

    VisMF::Write(mean,mfile);
    VisMF::Write(variance,vfile);

    amrex::Finalize();
}


