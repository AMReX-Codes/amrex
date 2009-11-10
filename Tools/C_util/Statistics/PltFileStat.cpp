
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
    std::string tmpFile;
    pp.query("infile", iFile);
    if (iFile.empty())
      BoxLib::Abort("You must specify `infile'");
    else
      tmpFile = iFile;

    int analysis_type = 0;
    pp.query("analysis", analysis_type);
    if (analysis_type == 0)
      BoxLib::Abort("Analysis type is not specified.");

    if (analysis_type == 4) {

      int nstart = 10;
      pp.query("nstart",nstart);

      char buftmp[64];
      sprintf(buftmp,"%05d",nstart);
          
      std::string idxs(buftmp);
      tmpFile = iFile + buftmp;
    }

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    DataServices dataServices(tmpFile, fileType);
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
    Array<Real> barr;
    vector<Real> bbll,bbur;
    bool do_bounds = 0;
    if (int nx=pp.countval("bounds"))
    {
      do_bounds = 1;
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

    int nBin = 0;
    pp.query("nBin", nBin);

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

      std::string oFile(tmpFile + "_PDF");
      pp.query("outfile",oFile);

      std::ofstream outputFile;
      outputFile.open(oFile.c_str(),std::ios::out);
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
    else if (analysis_type == 3) {// determine the simple variogram
       
      if (nBin == 0) {
	std::cout << "nBin was not specified.  Setting nBin to 100.\n";
	nBin = 100;
      }
      
      std::string oFile(tmpFile + "_VAR");
      pp.query("outfile",oFile);

      ComputeAmrDataVAR(amrData,nBin,cNames,barr,oFile);
      
    }

    else if (analysis_type == 4) {// same as 3 but for a list of files
       
      if (nBin == 0) {
	std::cout << "nBin was not specified.  Setting nBin to 100.\n";
	nBin = 100;
      }

      int nstart = 1;
      pp.query("nstart",nstart);
      int nmax = 100;
      pp.query("nfile",nmax);
      int nfac = 10;
      pp.query("nfac",nfac);

      for (int i = nstart; i<nmax;i++) {
	char buf[64];
	sprintf(buf,"%05d",i*nfac);
	std::string idxs(buf);
	tmpFile = iFile + idxs;


	DataServices tmpdataServices(tmpFile, fileType);
	if (!tmpdataServices.AmrDataOk())
	  DataServices::Dispatch(DataServices::ExitRequest, NULL);

	AmrData& tmpamrData = tmpdataServices.AmrDataRef();

	std::string oFile = tmpFile + "_VAR";

	ComputeAmrDataVAR(tmpamrData,nBin,cNames,barr,oFile);
      }      
    }


    else if (analysis_type == 5) 
    {
      std::cout << "GSIB variogram calculations.\n";
      if (nBin == 0) {
	std::cout << "nlag was not specified.  Setting nlag to 100.\n";
	nBin = 100;
      }
      
      std::string oFile(tmpFile + "_VAR");
      pp.query("outfile",oFile);

      int nvarg = pp.countname("varg");
      if (nvarg == 0)
	BoxLib::Abort("No variogram is specified");

      Array< Array<int> > ivoption(nvarg);
      for (int i=0; i<nvarg; i++) {
	int nopt = pp.countval("varg");
	ivoption[i].resize(nopt);
	pp.queryktharr("varg",i,ivoption[i],0,nopt);
      }

      int isill = 0;
      Array<Real> mean(nComp), variance(nComp);
      pp.query("isill",isill);
      if (isill == 1) 
	ComputeAmrDataMeanVar (amrData,cNames,bas,mean,variance);
      
      VariogramUniform(amrData,cNames,barr,ivoption,nBin,isill,variance,oFile);
      
    }
    else 
      std::cout << "Analysis Type is undefined.\n";
      
    
    BoxLib::Finalize();
}


