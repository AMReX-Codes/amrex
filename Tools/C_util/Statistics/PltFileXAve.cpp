
#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
using std::ios;
using std::set_new_handler;

#include <unistd.h>
#include "PltFileXAve.H"
#include "WritePlotFile.H"
#include "ComputeAmrDataStat.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
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

static
Real
SumThisComp(AmrData &amrData, int iComp)
{
  Real sum(0.0);
  Real phi(1.0);
  int finestLevel = amrData.FinestLevel();      
  int nComp = amrData.NComp();

  Array<int> refMult(finestLevel + 1, 1);
  for (int iLevel=finestLevel-1; iLevel>=0; --iLevel)
  {
    int ref_ratio = amrData.RefRatio()[iLevel];
    int vol = 1;
    for (int i=0; i<BL_SPACEDIM; ++i)
      vol *= ref_ratio;
    for (int jLevel=0; jLevel<=iLevel; ++jLevel)
      refMult[jLevel] *= vol;
  }

  //
  // Compute the sum and sum-squares
  //
  long total_volume = 0;
  Array<MultiFab*> error(finestLevel+1);
  
  for (int iLevel = finestLevel; iLevel>=0; --iLevel)
  {
    const BoxArray& ba = amrData.boxArray(iLevel);
    const Real *dx     = amrData.DxLevel()[iLevel].dataPtr();
    error[iLevel]      = new MultiFab(ba, nComp, 0);
    for (int i=0; i<nComp; ++i)
    {
      MultiFab& data = amrData.GetGrids(iLevel,i);
      error[iLevel]->copy(data,0,i,1);
    }

    // Zero out the error covered by fine grid
    long covered_volume = 0;
    if (iLevel != finestLevel)
    {
      int ref_ratio = amrData.RefRatio()[iLevel];	    
      BoxArray baF = ::BoxArray(amrData.boxArray(iLevel+1)).coarsen(ref_ratio);
      for (MFIter mfi(*error[iLevel]); mfi.isValid(); ++mfi)
      {
	for (int iGrid=0; iGrid<baF.size(); ++iGrid)
	{
	  Box ovlp = baF[iGrid] & mfi.validbox();
	  if (ovlp.ok())
	  {
	    (*error[iLevel])[mfi].setVal(0.0, ovlp, 0, nComp);
	    covered_volume += ovlp.numPts()*refMult[iLevel];
	  }
	}
      }
      ParallelDescriptor::ReduceLongSum(covered_volume);
    }

    // Compute volume at this level
    Real level_volume = 0.0;
    for (int iGrid=0; iGrid<ba.size(); ++iGrid)
      level_volume += ba[iGrid].numPts()*refMult[iLevel];
    level_volume -= covered_volume;

    if (level_volume > 0.0)
    {
      // Convert volume in numPts to volume in number of fine cells
      total_volume += long(level_volume);
      
      // Get norms at this level
      Real n1 = 0.0;
      
      for (MFIter mfi(*error[iLevel]); mfi.isValid(); ++mfi)
      {
	FArrayBox& fab = (*error[iLevel])[mfi];
	const Box& fabbox = mfi.validbox();
	const int* lo = fab.loVect();
	const int* hi = fab.hiVect();

#if (BL_SPACEDIM == 2)	
	for (int iy=lo[1]; iy<hi[1]+1; iy++) {
	  for (int ix=lo[0]; ix<hi[0]+1; ix++) {
	    Real sum = fab(IntVect(ix,iy),0) + fab(IntVect(ix,iy),1);
	    
	    if (sum > 0.0)
	      n1 += fab(IntVect(ix,iy),iComp)*dx[0]*dx[1];
	  }
	}
#else
	for (int iz=lo[2]; iz<hi[2]+1; iz++) {
	  for (int iy=lo[1]; iy<hi[1]+1; iy++) {
	    for (int ix=lo[0]; ix<hi[0]+1; ix++) {
	      Real sum = fab(IntVect(ix,iy,iz),0) + fab(IntVect(ix,iy,iz),1);
	    
	      if (sum > 0.0)
		n1 += fab(IntVect(ix,iy,iz),iComp)*dx[0]*dx[1]*dx[2];
	    }
	  }
	}
#endif
      }

      // Do necessary communication, then blend this level's norms
      //  in with the running global values
      ParallelDescriptor::ReduceRealSum(n1);
      
      sum += n1;
    }
  }
  

  // Clean up memory
  for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    delete error[iLevel];

  return sum;
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

    int nstart = 1;
    pp.query("nstart",nstart);
    int nmax = 100;
    pp.query("nfile", nmax);
    int nfac = 10;
    pp.query("nfac", nfac);
    //pp.query("outfile", outfile);
    //if (outfile.empty())
    //    BoxLib::Abort("You must specify `outfile'");

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);

    Real phi = 0.3;

    int analysis;
    pp.query("analysis",analysis);

    if (analysis == 0) 
      compute_flux_all(nstart, nmax, nfac, iFile);
    else if (analysis == 1) {

      Real flux  = 0;
      Real dtnew = 0;
      Real dtold = 0;
      Real dt;
      for (int i = nstart; i < nmax; i++) {
	char buf[64];
	sprintf(buf,"%05d",i*nfac);
          
	std::string idxs(buf);
	File = iFile + idxs;
      
	DataServices dataServices(File, fileType);

	if (!dataServices.AmrDataOk())
	  //
	  // This calls ParallelDescriptor::EndParallel() and exit()
	  //
	  DataServices::Dispatch(DataServices::ExitRequest, NULL);

	AmrData& amrData = dataServices.AmrDataRef();

	Array<string> cNames(2);
	cNames[0] = amrData.PlotVarNames()[0];
	cNames[1] = amrData.PlotVarNames()[1];
	
	dtnew = amrData.Time();
	dt    = dtnew - dtold;
	dtold = dtnew;
	
	compute_flux(amrData,0,cNames,dt,flux,phi);
	
	std::cout << dtnew <<  " " << flux << std::endl;
      }
    }
    else 
      BoxLib::Abort("Analysis type not defined");

    BoxLib::Finalize();

}

void 
compute_flux_all(int nstart,
		 int nmax,
		 int nfac,
		 std::string iFile)
{

  DataServices::SetBatchMode();
  FileType fileType(NEWPLT);

  int finestLevel;
  int nComp;
  Box dmn;
  Array<Real> xnew;
  Array<Real> xold;
  Array<Real> FLs;
  Real sumold,sumnew;
  MultiFab tmpmean;
  Real dtnew, dtold;
  dtold = 0.;
  sumold = 0.;
  Real phi = 0.3;
    
  for (int i = nstart; i < nmax; i++)
  {
    char buf[64];
    sprintf(buf,"%05d",i*nfac);
          
    std::string idxs(buf);
    std::string File = iFile + idxs;
    //File = hdrFile + idxs + iFile;
    //File = hdrFile + idxs;
    //if (i == 0) File +=idxs;
    //std::cout << File << std::endl;
      
    DataServices dataServices(File, fileType);

    if (!dataServices.AmrDataOk())
      //
      // This calls ParallelDescriptor::EndParallel() and exit()
      //
      DataServices::Dispatch(DataServices::ExitRequest, NULL);

    AmrData& amrData = dataServices.AmrDataRef();
    dtnew = amrData.Time();

    nComp = 2;
    Array<string> names(2);
    Array<int>    destcomp(2);

    names[0] = amrData.PlotVarNames()[0];
    names[1] = amrData.PlotVarNames()[1];
    destcomp[0] = 0;
    destcomp[1] = 1;
	

    if (i == nstart) {
      //finestLevel = amrData.FinestLevel();
      finestLevel = 0;
      BoxArray ba = amrData.boxArray(finestLevel);
      
      tmpmean.define(ba,nComp,0,Fab_allocate);

      dmn = amrData.ProbDomain()[finestLevel];
      //
      // Currently we assume dmn starts at zero.
      //
      xnew.resize(dmn.bigEnd(BL_SPACEDIM-1)+1);
      xold.resize(dmn.bigEnd(BL_SPACEDIM-1)+1);
      FLs.resize(dmn.bigEnd(BL_SPACEDIM-1)+1);

      for (int iy=0;iy<xold.size();iy++)
	xold[iy] = 0.0;
    }

    const Array<Real>& dx = amrData.DxLevel()[finestLevel];
    
    //      if (ParallelDescriptor::IOProcessor())
    //std::cout << "Filling tmpmean ... " << std::flush;
    
    amrData.FillVar(tmpmean,finestLevel,names,destcomp);
    
    for (int n = 0; n < nComp; n++)
      amrData.FlushGrids(destcomp[n]);
    
    Real dt = dtnew - dtold;
    dtold = dtnew;
    
    for (int ix = 0; ix < xnew.size(); ix++)
      xnew[ix] = 0;
    
    for (MFIter mfi(tmpmean); mfi.isValid(); ++mfi)
    {
      const int* lo = tmpmean[mfi].loVect();
      const int* hi = tmpmean[mfi].hiVect();
      
#if (BL_SPACEDIM == 2)
      for (int iy=lo[1]; iy<=hi[1]; iy++)
	for (int ix=lo[0]; ix<=hi[0]; ix++)
	  xnew[iy] += tmpmean[mfi](IntVect(ix,iy),0)*phi;
#else
      for (int iz=lo[2]; iz<=hi[2]; iz++)
	for (int iy=lo[1]; iy<=hi[1]; iy++)
	  for (int ix=lo[0]; ix<=hi[0]; ix++)
	    xnew[iz] += tmpmean[mfi](IntVect(ix,iy,iz),0)*phi;
#endif
    }

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    
    ParallelDescriptor::ReduceRealSum(xnew.dataPtr(), xnew.size(), IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
      Real FL = 0;

#if (BL_SPACEDIM == 2)
      for (int iy = 0; iy < xnew.size(); iy++)
      {
	xnew[iy] = xnew[iy]/dmn.length(0);
	Real ct = (xnew[iy]-xold[iy])/dt;
	FL += ct*dx[1];
	FLs[iy] = FL;
	xold[iy] = xnew[iy];
      }
#else
      for (int iz = 0; iz < xnew.size(); iz++)
      {
	xnew[iz] = xnew[iz]/(dmn.length(0)*dmn.length(1));
	Real ct = (xnew[iz]-xold[iz])/dt;
	FL += ct*dx[2];
	FLs[iz] = FL;
	xold[iz] = xnew[iz];
      }
#endif

      //sumnew = SumThisComp(amrData, 0);
      //Real dF = phi*(sumnew - sumold)/dt;
      //sumold = sumnew;
      //std::cout << dtnew << " " << FL << " " << FLs[896] << " " << FLs[768] << " " << FLs[512] << std::endl;
      std::cout << dtnew <<  " " << FL << std::endl;
    }

    //sumnew = SumThisComp(amrData, 0);
    //Real dF = phi*(sumnew - sumold)/dt;
    //sumold = sumnew;
    
  }    

}

void
compute_flux(AmrData&           amrData, 
	     int                dir, 
	     Array<std::string> cNames,
	     Real               dt,
	     Real               flux,
	     Real               phi,
 	     Real*              barr)
{ 
  Array<Real> xnew;
  Array<Real> xold;
  Array<Real> FLs;
  MultiFab tmpmean;
  Real sumold,sumnew;
  int finestLevel = 0;
  int nComp = cNames.size();

  Array<int> destFillComps(nComp);
  for (int i=0; i<nComp; ++i)
    destFillComps[i] = i;

  BoxArray ba;
  Box domain = amrData.ProbDomain()[0];
  if (barr == 0) {
    ba = amrData.boxArray(finestLevel);
  }
  else {
    vector<Real> bbll,bbur;
    BL_ASSERT(barr.size()==2*BL_SPACEDIM);
    bbll.resize(BL_SPACEDIM);
    bbur.resize(BL_SPACEDIM);
    for (int i=0; i<BL_SPACEDIM; ++i)
    {
      bbll[i] = barr[i];
      bbur[i] = barr[BL_SPACEDIM+i];
    }

    // Find coarse-grid coordinates of bounding box, round outwardly
    for (int i=0; i<BL_SPACEDIM; ++i) {
      const Real dx = amrData.ProbSize()[i]/ 
	amrData.ProbDomain()[0].length(i);
      domain.setSmall(i,std::max(domain.smallEnd()[i], 
	(int)((bbll[i]-amrData.ProbLo()[i]+.0001*dx)/dx)));
      domain.setBig(i,std::min(domain.bigEnd()[i], 
	(int)((bbur[i]-amrData.ProbLo()[i]-.0001*dx)/dx)));
    }

   //for (int i=1; i<=finestLevel; i++)
   //    domain.refine(amrData.RefRatio()[i]);

    ba = BoxLib::intersect(amrData.boxArray(0),domain);
  }
  tmpmean.define(ba,nComp,0,Fab_allocate);
  xnew.resize(domain.bigEnd(BL_SPACEDIM-1)-domain.smallEnd(BL_SPACEDIM-1)+1);
  xold.resize(domain.bigEnd(BL_SPACEDIM-1)-domain.smallEnd(BL_SPACEDIM-1)+1);
  FLs.resize (domain.bigEnd(BL_SPACEDIM-1)-domain.smallEnd(BL_SPACEDIM-1)+1);

  for (int iy=0;iy<xold.size();iy++)
    xold[iy] = 0.0;


  const Array<Real>& dx = amrData.DxLevel()[finestLevel];
      
  amrData.FillVar(tmpmean,finestLevel,cNames,destFillComps);

  for (int n = 0; n < nComp; n++)
    amrData.FlushGrids(destFillComps[n]);

  for (int ix = 0; ix < xnew.size(); ix++)
    xnew[ix] = 0;

  for (MFIter mfi(tmpmean); mfi.isValid(); ++mfi)
  {
    const int* lo = tmpmean[mfi].loVect();
    const int* hi = tmpmean[mfi].hiVect();
      
#if (BL_SPACEDIM == 2)
    for (int iy=lo[1]; iy<=hi[1]; iy++)
      for (int ix=lo[0]; ix<=hi[0]; ix++)
	xnew[iy] += tmpmean[mfi](IntVect(ix,iy),0)*phi;
#else
    for (int iz=lo[2]; iz<=hi[2]; iz++)
      for (int iy=lo[1]; iy<=hi[1]; iy++)
	for (int ix=lo[0]; ix<=hi[0]; ix++)
	  xnew[iz] += tmpmean[mfi](IntVect(ix,iy,iz),0)*phi;
#endif
  }

  const int IOProc = ParallelDescriptor::IOProcessorNumber();

  ParallelDescriptor::ReduceRealSum(xnew.dataPtr(), xnew.size(), IOProc);

  if (ParallelDescriptor::IOProcessor())
  {
    Real FL = 0;

#if (BL_SPACEDIM == 2)
    for (int iy = 0; iy < xnew.size(); iy++)
    {
      xnew[iy] = xnew[iy]/domain.length(0);
      Real ct = (xnew[iy]-xold[iy])/dt;
      FL += ct*dx[1];
      FLs[iy] = FL;
      xold[iy] = xnew[iy];
    }
#else
    for (int iz = 0; iz < xnew.size(); iz++)
    {
      xnew[iz] = xnew[iz]/(domain.length(0)*domain.length(1));
      Real ct = (xnew[iz]-xold[iz])/dt;
      FL += ct*dx[2];
      FLs[iz] = FL;
      xold[iz] = xnew[iz];
    }
#endif
    flux = FL;  
  }
}


