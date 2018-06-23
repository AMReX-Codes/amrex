
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <unistd.h>

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

    bool verbose = false;
    if (pp.contains("verbose"))
    {
      verbose = true;
      AmrData::SetVerbose(true);
    }
    std::string tmpFile;
    pp.query("infile", iFile);
    if (iFile.empty())
      amrex::Abort("You must specify `infile'");
    else
      tmpFile = iFile;

    int analysis_type = 0;
    pp.query("analysis", analysis_type);
    if (analysis_type == 0)
      amrex::Abort("Analysis type is not specified.");

    if (analysis_type == 4) {

      int nstart = 10;
      pp.query("nstart",nstart);

      tmpFile = amrex::Concatenate(iFile, nstart, 5);
    }

    int iFile_type = 0;
    pp.query("infile_type", iFile_type);

    int nBin = 0;
    pp.query("nBin", nBin);
    if (nBin == 0) 
    {
      std::cout << "nBin was not specified.  Setting nBin to 100.\n";
      nBin = 100;
    }

    // file is a plotfile
    if (iFile_type == 0) 
    {
      Amrvis::FileType fileType(Amrvis::NEWPLT);
      DataServices::SetBatchMode();   
      DataServices dataServices(tmpFile,fileType);

      if (!dataServices.AmrDataOk())
	DataServices::Dispatch(DataServices::ExitRequest, NULL);

      AmrData& amrData = dataServices.AmrDataRef();
      int finestLevel = amrData.FinestLevel();
      int Nlev = finestLevel + 1;

      // limit analysis to components specified
      Vector<std::string> cNames;
      Vector<std::string> varNames = amrData.PlotVarNames();
      if (int nx=pp.countval("cNames"))
      {
	Vector<int> var_match(nx,0);
	pp.getarr("cNames",cNames,0,nx);
	for (int ic=0; ic<nx; ic++) {
	  for (int iv=0; iv<varNames.size(); iv++){
	    if (cNames[ic].compare(varNames[iv]) != 0) {
	      var_match[ic] = 1;
	    }
	  }
	  if (var_match[ic] == 0) {
	    std::cout << "Component " << cNames[ic] << " is not defined!\n";
	    amrex::Abort("Please check that the specified component is defined.");
	  }
	}
      }
      else 
      {
	cNames.resize(varNames.size());
	for (int iv=0;iv<varNames.size();iv++)
	  cNames[iv] = varNames[iv];
      }
      int nComp = cNames.size();
    
      // limit to a smaller region    
      Box domain = amrData.ProbDomain()[0];
      Vector<Real> barr;
      Vector<BoxArray> bas;
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
	  const Real dx = amrData.ProbSize()[i] / 
	                  amrData.ProbDomain()[0].length(i);            
	  domain.setSmall(i,std::max(domain.smallEnd()[i], 
		     (int)((bbll[i]-amrData.ProbLo()[i]+.0001*dx)/dx)));
	  domain.setBig(i,std::min(domain.bigEnd()[i], 
		     (int)((bbur[i]-amrData.ProbLo()[i]-.0001*dx)/dx)));
	}
      }
      else
      {
	barr.resize(2*BL_SPACEDIM);
	for (int i=0;i<BL_SPACEDIM;++i)
	{
	  barr[i] = domain.smallEnd()[i];
	  barr[BL_SPACEDIM+i] = domain.bigEnd()[i];
	}
      }
      // Build boxarrays for fillvar call
      Box levelDomain = domain;
      bas.resize(Nlev);
      for (int iLevel=0; (iLevel<=finestLevel)&&(bas.size()==Nlev); ++iLevel)
      {
        BoxArray baThisLev = amrex::intersect(amrData.boxArray(iLevel),levelDomain);

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
    
      // analysis options
      // 1: mean and variance of a plot file
      // 2: determine the pdf of a plot file
      // 3: 
      
      if (analysis_type == 1) { // mean and variance of a plot file
	Vector<Real> mean, variance;
	AmrData& amrData =  dataServices.AmrDataRef();
	ComputeAmrDataMeanVar(amrData, mean, variance, 0, nComp, verbose);

	for (int i=0;i<nComp;i++) 
	{
	  std::cout << " comp= " << i 
		    << " mean = " << mean[i] 
		    << " variance = " << variance[i] << std::endl;
	}
      }

      else if (analysis_type == 2) {// determine the pdf of a plot file
	
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

      else if (analysis_type == 3) // determine the simple variogram
      {
	std::string oFile(tmpFile + "_VAR");
	pp.query("outfile",oFile);
	ComputeAmrDataVAR(amrData,nBin,cNames,barr,oFile);
      }

      else if (analysis_type == 4) // same as 3 but for a list of files
      {	
	int nstart = 1;
	pp.query("nstart",nstart);
	int nmax = 100;
	pp.query("nfile",nmax);
	int nfac = 10;
	pp.query("nfac",nfac);

	for (int i = nstart; i<nmax;i++) {

          tmpFile = amrex::Concatenate(iFile, i*nfac, 5);

	  DataServices tmpdataServices(tmpFile, fileType);
	  if (!tmpdataServices.AmrDataOk())
	    DataServices::Dispatch(DataServices::ExitRequest, NULL);
	  AmrData& tmpamrData = tmpdataServices.AmrDataRef();

	  std::string oFile = tmpFile + "_VAR";
	  ComputeAmrDataVAR(tmpamrData,nBin,cNames,barr,oFile);
	}      
      }


      else if (analysis_type == 5) // determine the variogram based on GSLIB
      {
	std::cout << "GSLIB variogram calculations.\n";
      
	std::string oFile(tmpFile + "_VAR");
	pp.query("outfile",oFile);
	
	int nvarg = pp.countname("varg");
	if (nvarg == 0)
	  amrex::Abort("No variogram is specified");

	Vector< Vector<int> > ivoption(nvarg);
	for (int i=0; i<nvarg; i++) {
	  int nopt = pp.countval("varg");
	  ivoption[i].resize(nopt);
	  pp.queryktharr("varg",i,ivoption[i],0,nopt);
	}
	
	int isill = 0;
	pp.query("isill",isill);
	
	Vector<Real> mean(nComp), variance(nComp);
	if (isill == 1) 
	  ComputeAmrDataMeanVar (amrData,cNames,bas,mean,variance);

	VariogramUniform(amrData,cNames,barr,
			 ivoption,nBin,isill,variance,oFile);
      }
    

      else if (analysis_type == 6) // compare coarse and fine solution
      {
	std::cout << "Analysis 6: take difference on fine grid.\n";

	std::string crsefile;
	pp.query("crsefile",crsefile);
	if (crsefile.empty()) 
	  amrex::Abort("You must specify `crsefile'");

	std::string oFile = iFile + "_diffine";
	pp.query("outfile",oFile);
      
	DataServices crseDataServices(crsefile,fileType);
	if (!crseDataServices.AmrDataOk())
	  DataServices::Dispatch(DataServices::ExitRequest, NULL);
	AmrData& amrDataCrse = crseDataServices.AmrDataRef();
	const Real dtcrse = amrDataCrse.Time();
	
	// The I/O processor makes the directory if it doesn't already exist.
	if (ParallelDescriptor::IOProcessor())
	  if (!amrex::UtilCreateDirectory(oFile, 0755))
	  amrex::CreateDirectoryFailed(oFile);
	ParallelDescriptor::Barrier();
	std::string mFile = iFile + "_dif";
	TakeDifferenceFine(amrData,amrDataCrse,cNames,barr,mFile);
      }

      else if (analysis_type == 7) // compare coarse and fine solution
      {
	std::cout << "Analysis 6: coarse-fine comparison.\n";

	std::string crsefile;
	pp.query("crsefile",crsefile);
	if (crsefile.empty()) 
	  amrex::Abort("You must specify `crsefile'");
	
	std::string oFile;
	pp.query("outfile",oFile);
	
	DataServices crseDataServices(crsefile,fileType);
	if (!crseDataServices.AmrDataOk())
	  DataServices::Dispatch(DataServices::ExitRequest, NULL);
	AmrData& amrDataCrse = crseDataServices.AmrDataRef();
	const Real dtcrse = amrDataCrse.Time();
	
	// The I/O processor makes the directory if it doesn't already exist.
	if (ParallelDescriptor::IOProcessor())
	  if (!amrex::UtilCreateDirectory(oFile, 0755))
	    amrex::CreateDirectoryFailed(oFile);
	ParallelDescriptor::Barrier();
	std::string mFile = iFile + "_dif";
	TakeDifferenceCrse(amrData,amrDataCrse,cNames,barr,mFile);
	
      }

      else if (analysis_type == 8) // compare coarse and fine solution
      {
	std::string crsefile;
	pp.query("crsefile",crsefile);
	if (crsefile.empty()) 
	  amrex::Abort("You must specify `crsefile'");

	std::string oFile;
	pp.query("outfile",oFile);
	
	DataServices crseDataServices(crsefile,fileType);
	if (!crseDataServices.AmrDataOk())
	  DataServices::Dispatch(DataServices::ExitRequest, NULL);
	AmrData& amrDataCrse = crseDataServices.AmrDataRef();
	const Real dtcrse = amrDataCrse.Time();

	// The I/O processor makes the directory if it doesn't already exist.
	if (ParallelDescriptor::IOProcessor())
	  if (!amrex::UtilCreateDirectory(oFile, 0755))
	    amrex::CreateDirectoryFailed(oFile);
	ParallelDescriptor::Barrier();
	std::string mFile = iFile + "_dif";
	TakeDifferenceSum(amrData,amrDataCrse,cNames,barr,mFile);
      
      }

      else if (analysis_type == 9) // determine the variogram based on GSLIB
      {
	std::cout << "GSLIB cross-variogram calculations.\n";
      
	std::string oFile(tmpFile + "_VAR");
	pp.query("outfile",oFile);
	
	std::string sFile;
	pp.query("secfile",sFile);
	
	int nvarg = pp.countname("varg");
	if (nvarg == 0)
	  amrex::Abort("No variogram is specified");

	Vector< Vector<int> > ivoption(nvarg);
	for (int i=0; i<nvarg; i++) {
	  int nopt = pp.countval("varg");
	  ivoption[i].resize(nopt);
	  pp.queryktharr("varg",i,ivoption[i],0,nopt);
	}
	
	int isill = 0;
	pp.query("isill",isill);

	MultiFab secmf;
	VisMF::Read(secmf,sFile);
	
	Vector<Real> mean, variance;
	Vector<Real> secmean, secvariance;
	if (isill == 1) {

	  mean.resize(nComp+secmf.nComp());
	  variance.resize(nComp+secmf.nComp());
	  secmean.resize(secmf.nComp());
	  secvariance.resize(secmf.nComp());

	  ComputeAmrDataMeanVar(amrData,cNames,bas,mean,variance);
	  ComputeMeanVarMF(secmf,secmean,secvariance);

	  for (int i=0;i < secmf.nComp(); i++) {
	    mean[nComp+i] = secmean[i];
	    variance[nComp+i] = secvariance[i];
	  }

	}

	VariogramCross(amrData,cNames,secmf,barr,
		       ivoption,nBin,isill,variance,oFile);
      }

      else if (analysis_type == 10) // compare coarse and fine solution
      {
	std::cout << "Analysis 10: difference from mean.\n";
	
	std::string oFile;
	pp.query("outfile",oFile);

	Vector<int> rratio(BL_SPACEDIM,0);
	if (int nx=pp.countval("rratio"))
	  pp.getarr("rratio",rratio,0,BL_SPACEDIM);
	for (int i=0;i<BL_SPACEDIM; i++)
	  if (rratio[i] == 0) amrex::Abort("rratio must be nonzero.");

	// The I/O processor makes the directory if it doesn't already exist.
	if (ParallelDescriptor::IOProcessor())
	  if (!amrex::UtilCreateDirectory(oFile, 0755))
	    amrex::CreateDirectoryFailed(oFile);
	ParallelDescriptor::Barrier();
	std::string mFile = iFile + "_dif";
	TakeDifferenceMean(amrData,cNames,barr,rratio,mFile);
	
      }
      else 
	std::cout << "Analysis Type is undefined.\n";

    }

    if (iFile_type == 1)
    {
      if (analysis_type == 5) 
      {
	std::cout << "GSLIB variogram calculations.\n";
      
	std::string oFile(tmpFile + "_VAR");
	pp.query("outfile",oFile);
	
	int nvarg = pp.countname("varg");
	if (nvarg == 0)
	  amrex::Abort("No variogram is specified");

	Vector< Vector<int> > ivoption(nvarg);
	for (int i=0; i<nvarg; i++) {
	  int nopt = pp.countval("varg");
	  ivoption[i].resize(nopt);
	  pp.queryktharr("varg",i,ivoption[i],0,nopt);
	}
	
	int isill = 0;
	pp.query("isill",isill);
	

	MultiFab mfIn;
	VisMF::Read(mfIn,iFile);
	  
	Vector<Real> mean(mfIn.nComp(),0.0), variance(mfIn.nComp(),0.0);
	Vector<Real> dx(BL_SPACEDIM);
	dx[0] = 2;
	dx[1] = 2;

	if (isill == 1)
	  ComputeMeanVarMF (mfIn,mean,variance);
	VariogramUniformMF(mfIn,dx,ivoption,nBin,isill,variance,oFile); 
      }     
    }

    amrex::Finalize();
    std::cout << "Done ...\n";

}


