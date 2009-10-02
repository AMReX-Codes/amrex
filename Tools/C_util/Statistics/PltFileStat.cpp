
#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
using std::ios;
using std::set_new_handler;

#include <unistd.h>

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

    bool verbose = false;
    if (pp.contains("verbose"))
    {
      verbose = true;
      AmrData::SetVerbose(true);
    }
    pp.query("infile", iFile);
    if (iFile.empty())
      BoxLib::Abort("You must specify `infile'");

    int analysis_type = 0;
    pp.query("analysis", analysis_type);
    if (analysis_type == 0)
      BoxLib::Abort("Analysis type is not specified.");

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    DataServices dataServices(iFile, fileType);
    if (!dataServices.AmrDataOk())
      //
      // This calls ParallelDescriptor::EndParallel() and exit()
      //
      DataServices::Dispatch(DataServices::ExitRequest, NULL);

    AmrData& amrData = dataServices.AmrDataRef();
    int finestLevel = amrData.FinestLevel();
    int Nlev = finestLevel + 1;

    // Limit to selected components
    Array<std::string> cNames;
    Array<std::string> varNames = amrData.PlotVarNames();
    if (int nx=pp.countval("cNames"))
    {
      Array<int> var_match(nx,0);
      pp.getarr("cNames",cNames,0,nx);
      for (int ic=0; ic<nx; ic++) {
	for (int iv=0; iv<varNames.size(); iv++){
	  if (cNames[ic].compare(varNames[iv]) != 0) {
	    var_match[ic] = 1;
	  }
	}
	if (var_match[ic] == 0) {
	  std::cout << "Component " << cNames[ic] << " is not defined!\n";
	  BoxLib::Abort("Please check that the specified component is defined.");
	}
      }
    }
    else {
      cNames.resize(varNames.size());
      for (int iv=0;iv<varNames.size();iv++)
	cNames[iv] = varNames[iv];
    }
    int nComp = cNames.size();
    
    // Limit to a smaller region    
    Box domain = amrData.ProbDomain()[0];
    vector<Real> bbll,bbur;
    bool do_bounds = 0;
    if (int nx=pp.countval("bounds"))
    {
      do_bounds = 1;
      Array<Real> barr;
      pp.getarr("bounds",barr,0,nx);
      int d=BL_SPACEDIM;
      BL_ASSERT(barr.size()==2*d);
      bbll.resize(d);
      bbur.resize(d);
      for (int i=0; i<d; ++i)
      {
	bbll[i] = barr[i];
	bbur[i] = barr[d+i];
      }

      // Find coarse-grid coordinates of bounding box, round outwardly
      for (int i=0; i<BL_SPACEDIM; ++i) {
	const Real dx = amrData.ProbSize()[i] / amrData.ProbDomain()[0].length(i);            
	domain.setSmall(i,std::max(domain.smallEnd()[i], 
				   (int)((bbll[i]-amrData.ProbLo()[i]+.0001*dx)/dx)));
	domain.setBig(i,std::min(domain.bigEnd()[i], 
				 (int)((bbur[i]-amrData.ProbLo()[i]-.0001*dx)/dx)));
      }
    }

    // Build boxarrays for fillvar call
    Box levelDomain = domain;
    Array<BoxArray> bas(Nlev);
    for (int iLevel=0; (iLevel<=finestLevel)&&(bas.size()==Nlev); ++iLevel)
    {
        BoxArray baThisLev = BoxLib::intersect(amrData.boxArray(iLevel),levelDomain);

        if (baThisLev.size() > 0) {
            bas.set(iLevel,baThisLev);
            if (iLevel < finestLevel) {
                levelDomain.refine(amrData.RefRatio()[iLevel]);
            }
        }
        else
        {
            bas.resize(iLevel);
        }
    }
    if (analysis_type == 1) { // mean and variance of a plot file
      Array<Real> mean, variance;
      ComputeAmrDataMeanVar(amrData, mean, variance, 0, nComp, verbose);

      for (int i=0;i<nComp;i++) {
	std::cout << " comp= " << i 
	          << " mean = " << mean[i] 
	          << " variance = " << variance[i] << std::endl;
      }
    }

    else if (analysis_type == 2) {// determine the pdf of a plot file

      int nBin = 0;
      pp.query("nBin", nBin);
      if (nBin == 0) {
	std::cout << "nBin was not specified.  Setting nBin to 100.\n";
	nBin = 100;
      }

      Real* icount[nComp];
      for (int ic=0;ic<nComp;ic++)
	icount[ic]=new Real[nBin];

      if (do_bounds)
	ComputeAmrDataPDF(amrData,icount,nBin,cNames,bas);
      else
	ComputeAmrDataPDF(amrData,icount,nBin,cNames);

      std::string oFile(iFile + "_PDF");
      pp.query("outfile",oFile);

      char outputFileName[oFile.size()];
      for (int ic=0; ic<oFile.size(); ic++)
	outputFileName[ic] = oFile[ic];

      std::ofstream outputFile;
      outputFile.open(outputFileName,std::ios::out);
      if (outputFile.fail()) return (1);
      for (int ic=0; ic<nComp; ic++) {
	for (int ib=0; ib<nBin; ib++) {		
	  outputFile << icount[ic][ib] << " ";
	}
	outputFile << "\n";
      }
      outputFile.close();

      for (int ic=0;ic<nComp;ic++)
	delete [] icount[ic];
    }
     else if (analysis_type == 3) {// determine the correlation in space
       
       std::cout << "In progress\n";
       
     }

    else 
      std::cout << "Analysis Type is undefined.\n";
      
    
    BoxLib::Finalize();
}


