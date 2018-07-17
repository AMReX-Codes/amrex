
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <unistd.h>

#include <PltFileFluxAve.H>
#include <WritePlotFile.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
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
    std::cout << "This routine reads a set of pltfiles "
              << "and determines x/y fluxes. " << std::endl
	      << std::endl;
    std::cout << "Usage:" << '\n';
    std::cout << progName << '\n';
    std::cout << "   infile=inputFileName" << '\n';
    std::cout << "   [outfile=outputFileName]" << '\n';
    std::cout << "   [-help]" << '\n';
    std::cout << "   [-verbose]" << '\n';
    std::cout << '\n';
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

  Vector<int> refMult(finestLevel + 1, 1);
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
  Vector<MultiFab*> error(finestLevel+1);
  
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
      amrex::Abort("You must specify `infile'");

    int nstart = 1;
    pp.query("nstart",nstart);
    int nmax = 100;
    pp.query("nfile", nmax);
    int nfac = 10;
    pp.query("nfac", nfac);
    int dir = 1;
    pp.query("dir", dir);


    // limit to a smaller region   
    Vector<Real> barr;
    bool do_bounds = 0;
    if (int nx=pp.countval("bounds"))
    {
      do_bounds = 1;
      pp.getarr("bounds",barr,0,nx);
    }

    std::string pfile;
    MultiFab phidata;
    pp.query("phifile",pfile);
    if (!pfile.empty()) {
      pfile += "/pp";
      VisMF::Read(phidata,pfile);
    }
    Real phi = 0.;
    pp.query("phi",phi);

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    

    int analysis;
    pp.query("analysis",analysis);

    if (analysis == 0) { 
      if (pfile.empty()) {
	compute_flux_all(nstart, nmax, nfac, dir, iFile, phi);
      }
      else
	compute_flux_all(nstart, nmax, nfac, dir, iFile, phidata);
    }
    else if (analysis == 1) {

      bool do_init = false;
      Real flux  = 0;
      Real dtnew = 0;
      Real dtold = 0;
      Real dt;
      Vector<string> cNames(2);
      Vector<Real> xold;
      for (int i = nstart; i < nmax; i++) {

        File = amrex::Concatenate(iFile, i*nfac, 5);
      
	DataServices dataServices(File, fileType);

	if (!dataServices.AmrDataOk())
	  //
	  // This calls ParallelDescriptor::EndParallel() and exit()
	  //
	  DataServices::Dispatch(DataServices::ExitRequest, NULL);

	AmrData& amrData = dataServices.AmrDataRef();
	dtnew = amrData.Time();
	dt    = dtnew - dtold;
	dtold = dtnew;

	if (i == nstart) {
	  cNames[0] = amrData.PlotVarNames()[0];
	  cNames[1] = amrData.PlotVarNames()[1];
	  do_init = true;
	}
	else
	  do_init = false;
	
	compute_flux(amrData,dir,cNames,dt,xold,flux,phi,do_init,barr);
		
	std::cout << dtnew <<  " " << flux << std::endl;
      }
    }
    else 
      amrex::Abort("Analysis type not defined");

    amrex::Finalize();

}

void 
compute_flux_all(int nstart,
		 int nmax,
		 int nfac,
		 int dir,
		 std::string iFile,
		 Real phi)
{
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);

  int finestLevel;
  int nComp;
  Box dmn;
  Vector<Real> xnew;
  Vector<Real> xold;
  Vector<Real> FLs;
  Real sumnew;
  MultiFab tmpmean;
  Real dtnew, dtold;
  dtold = 0.;
    
  for (int i = nstart; i < nmax; i++)
  {
    std::string File = amrex::Concatenate(iFile, i*nfac, 5);
      
    DataServices dataServices(File, fileType);

    if (!dataServices.AmrDataOk())
      //
      // This calls ParallelDescriptor::EndParallel() and exit()
      //
      DataServices::Dispatch(DataServices::ExitRequest, NULL);

    AmrData& amrData = dataServices.AmrDataRef();
    dtnew = amrData.Time();

    nComp = 2;
    Vector<string> names(2);
    Vector<int>    destcomp(2);

    names[0] = amrData.PlotVarNames()[0];
    names[1] = amrData.PlotVarNames()[1];
    destcomp[0] = 0;
    destcomp[1] = 1;
	

    if (i == nstart) {
      finestLevel = 0;
      BoxArray ba = amrData.boxArray(finestLevel);
      
      tmpmean.define(ba,nComp,0,Fab_allocate);

      dmn = amrData.ProbDomain()[finestLevel];
      //
      // Currently we assume dmn starts at zero.
      //
      xnew.resize(dmn.bigEnd(dir)+1);
      xold.resize(dmn.bigEnd(dir)+1);
      FLs.resize(dmn.bigEnd(dir)+1);

      for (int iy=0;iy<xold.size();iy++)
	xold[iy] = 0.0;
    }

    const Vector<Real>& dx = amrData.DxLevel()[finestLevel];
    
    amrData.FillVar(tmpmean,finestLevel,names,destcomp);
    
    for (int n = 0; n < nComp; n++)
      amrData.FlushGrids(destcomp[n]);
    
    Real dt = dtnew - dtold;
    dtold = dtnew;
    
    for (int ix = 0; ix < xnew.size(); ix++)
      xnew[ix] = 0.0;
    
    Real dmn_length;
#if (BL_SPACEDIM == 2)
    if (dir == 0)dmn_length = dmn.length(1);
    if (dir == 1)dmn_length = dmn.length(0);
#else
    if (dir == 0)dmn_length = dmn.length(1)*dmn.length(2);
    if (dir == 1)dmn_length = dmn.length(0)*dmn.length(2);
    if (dir == 2)dmn_length = dmn.length(0)*dmn.length(1);
#endif

    for (MFIter mfi(tmpmean); mfi.isValid(); ++mfi)
    {
      const int* lo = tmpmean[mfi].loVect();
      const int* hi = tmpmean[mfi].hiVect();
      
#if (BL_SPACEDIM == 2)
      if (dir == 0) {
	for (int ix=lo[0]; ix<=hi[0]; ix++){
	  for (int iy=lo[1]; iy<=hi[1]; iy++) {
	    xnew[ix] += tmpmean[mfi](IntVect(ix,iy),0)*phi;
	  }
	}
      }
      if (dir == 1) {
	for (int iy=lo[1]; iy<=hi[1]; iy++){
	  for (int ix=lo[0]; ix<=hi[0]; ix++) {
	    xnew[iy] += tmpmean[mfi](IntVect(ix,iy),0)*phi;
	  }
	}
      }
#else
      if (dir == 2) {
	for (int iz=lo[2]; iz<=hi[2]; iz++)
	  for (int iy=lo[1]; iy<=hi[1]; iy++)
	    for (int ix=lo[0]; ix<=hi[0]; ix++)
	      xnew[iz] += tmpmean[mfi](IntVect(ix,iy,iz),0)*phi;
      }
#endif
    }

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
     
    ParallelDescriptor::ReduceRealSum(xnew.dataPtr(), xnew.size(), IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
      Real FL = 0;

      for (int iy = xnew.size()-1; iy >=0; iy--)
      {
	xnew[iy] = xnew[iy]/dmn_length;
	Real ct = (xnew[iy]-xold[iy])/dt;
	FL += ct*dx[dir];
	FLs[iy] = FL;
	xold[iy] = xnew[iy];
      }

      std::cout << dtnew <<  " " << FL << std::endl;
    }

    //sumnew = SumThisComp(amrData, 0);
    //Real dF = phi*(sumnew - sumold)/dt;
    //sumold = sumnew;
    
  }    

}
void 
compute_flux_all(int nstart,
		 int nmax,
		 int nfac,
		 int dir,
		 std::string iFile,
		 MultiFab& phidata)
{
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);

  int finestLevel;
  int nComp=2;
  Box dmn;
  Vector<Real> xnew;
  Vector<Real> xold;
  Vector<Real> FLs;
  Vector<string> names(2);
  Vector<int>    destcomp(2);
  Real sumnew;
  MultiFab tmpmean;
  MultiFab tmpphi; 
  Real dtnew, dtold;
  dtold = 0.;
    
  for (int i = nstart; i < nmax; i++)
  {
    std::string File = amrex::Concatenate(iFile, i*nfac, 5);
      
    DataServices dataServices(File, fileType);

    if (!dataServices.AmrDataOk())
      //
      // This calls ParallelDescriptor::EndParallel() and exit()
      //
      DataServices::Dispatch(DataServices::ExitRequest, NULL);

    AmrData& amrData = dataServices.AmrDataRef();
    dtnew = amrData.Time();
    
    if (i == nstart) {

      names[0] = amrData.PlotVarNames()[0];
      names[1] = amrData.PlotVarNames()[1];
      destcomp[0] = 0;
      destcomp[1] = 1;

      finestLevel = amrData.FinestLevel();
      dmn = amrData.ProbDomain()[finestLevel];
      BoxArray ba(dmn);
      ba.maxSize(128);

      //int baseLevel = 0;
      //BoxArray ba = amrData.boxArray(baseLevel);
      
      int ng_twoexp = 1;
      for (int ii = 0; ii<finestLevel; ii++) {
	ng_twoexp *= amrData.RefRatio()[ii];
      }	
      ng_twoexp = ng_twoexp*3;
      BoxArray ba2(ba.size());
      for (int ii=0; ii<ba.size(); ii++) {
	Box bx=ba[ii];
	bx.grow(ng_twoexp);
	ba2.set(ii,bx);
      }
      
      MultiFab mftmp(ba2,1,0);
      mftmp.copy(phidata);
      
      tmpmean.define(ba,nComp,0,Fab_allocate);
      tmpphi.define(ba,nComp,0,Fab_allocate);
      for (MFIter mfi(mftmp);mfi.isValid();++mfi)
	tmpphi[mfi].copy(mftmp[mfi]);
      mftmp.clear();      

      //
      // Currently we assume dmn starts at zero.
      //
      xnew.resize(dmn.bigEnd(dir)+1);
      xold.resize(dmn.bigEnd(dir)+1);
      FLs.resize(dmn.bigEnd(dir)+1);

      for (int iy=0;iy<xold.size();iy++)
	xold[iy] = 0.0;
    }

    const Vector<Real>& dx = amrData.DxLevel()[finestLevel];
    
    amrData.FillVar(tmpmean,finestLevel,names,destcomp);
    
    for (int n = 0; n < nComp; n++)
      amrData.FlushGrids(destcomp[n]);
    
    Real dt = dtnew - dtold;
    dtold = dtnew;
    
    for (int ix = 0; ix < xnew.size(); ix++)
      xnew[ix] = 0.0;

    Real dmn_length;
#if (BL_SPACEDIM == 2)
    if (dir == 0)dmn_length = dmn.length(1);
    if (dir == 1)dmn_length = dmn.length(0);
#else
    if (dir == 0)dmn_length = dmn.length(1)*dmn.length(2);
    if (dir == 1)dmn_length = dmn.length(0)*dmn.length(2);
    if (dir == 2)dmn_length = dmn.length(0)*dmn.length(1);
#endif
    
    for (MFIter mfi(tmpmean); mfi.isValid(); ++mfi)
    {
      const int* lo = tmpmean[mfi].loVect();
      const int* hi = tmpmean[mfi].hiVect();
      
#if (BL_SPACEDIM == 2)
      if (dir == 0) {
	for (int ix=lo[0]; ix<=hi[0]; ix++) {
	  for (int iy=lo[1]; iy<=hi[1]; iy++) {
	    xnew[ix] += tmpmean[mfi](IntVect(ix,iy),0)*
	      tmpphi[mfi](IntVect(ix,iy),0);
	  }
	}
      }
      if (dir == 1) {
	for (int iy=lo[1]; iy<=hi[1]; iy++) {
	  for (int ix=lo[0]; ix<=hi[0]; ix++) {
	    xnew[iy] += tmpmean[mfi](IntVect(ix,iy),0)*
	      tmpphi[mfi](IntVect(ix,iy),0);
	  }
	}
      }
#else
      if (dir == 2) {
	for (int iz=lo[2]; iz<=hi[2]; iz++)
	  for (int iy=lo[1]; iy<=hi[1]; iy++)
	    for (int ix=lo[0]; ix<=hi[0]; ix++)
	      xnew[iz] += tmpmean[mfi](IntVect(ix,iy,iz),0)*
		tmpphi[mfi](IntVect(ix,iy),0);
      }
#endif
    }

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    
    ParallelDescriptor::ReduceRealSum(xnew.dataPtr(), xnew.size(), IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
      Real FL = 0;

      for (int iy = 0; iy < xnew.size(); iy++)
      {
	xnew[iy] = xnew[iy]/dmn_length;
	Real ct = (xnew[iy]-xold[iy])/dt;
	FL += ct*dx[dir];
	FLs[iy] = FL;
	xold[iy] = xnew[iy];
	std::cout << ct << " " <<dx[dir] << " " << FL << std::endl;
      }

      std::cout << dtnew <<  " " << FL << std::endl;

    }
  }    
}

void
compute_flux(AmrData&           amrData, 
	     int                dir, 
	     Vector<std::string> cNames,
	     Real               dt,
	     Vector<Real>&       xold,
	     Real&              flux,
	     Real               phi,
	     bool               do_init,
 	     Vector<Real>&       barr)
{ 
  Vector<Real> xnew;
  Vector<Real> FLs;
  MultiFab tmpmean;
  Real sumnew;
  int finestLevel = 0;
  int nComp = cNames.size();

  Vector<int> destFillComps(nComp);
  for (int i=0; i<nComp; ++i)
    destFillComps[i] = i;
  
  Box domain = amrData.ProbDomain()[finestLevel];
  xnew.resize(domain.bigEnd(dir)-domain.smallEnd(dir)+1);
  for (int ix = 0; ix < xnew.size(); ix++)
    xnew[ix] = 0.0;

  BoxArray ba;
  if (barr.empty()) {
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

    ba = amrex::intersect(amrData.boxArray(0),domain);
  }

  tmpmean.define(ba,nComp,0,Fab_allocate);
  FLs.resize (domain.bigEnd(dir)-domain.smallEnd(dir)+1);

  if (do_init && ParallelDescriptor::IOProcessor()) {
    xold.resize(domain.bigEnd(dir)-domain.smallEnd(dir)+1);
    for (int iy=0;iy<xold.size();iy++)
      xold[iy] = 0.0;
  }
  const Vector<Real>& dx = amrData.DxLevel()[finestLevel];
      
  amrData.FillVar(tmpmean,finestLevel,cNames,destFillComps);

  for (int n = 0; n < nComp; n++)
    amrData.FlushGrids(destFillComps[n]);

    Real dmn_length;
#if (BL_SPACEDIM == 2)
    if (dir == 0)dmn_length = domain.length(1);
    if (dir == 1)dmn_length = domain.length(0);
#else
    if (dir == 0)dmn_length = domain.length(1)*domain.length(2);
    if (dir == 1)dmn_length = domain.length(0)*domain.length(2);
    if (dir == 2)dmn_length = domain.length(0)*domain.length(1);
#endif

  for (MFIter mfi(tmpmean); mfi.isValid(); ++mfi)
  {
    const int* lo = tmpmean[mfi].loVect();
    const int* hi = tmpmean[mfi].hiVect();
      
#if (BL_SPACEDIM == 2)
    if (dir == 0) {
      for (int ix=lo[0]; ix<=hi[0]; ix++) {
	for (int iy=lo[1]; iy<=hi[1]; iy++) {
	  xnew[ix] += tmpmean[mfi](IntVect(ix,iy),0)*phi;
	}
      }
    }
    if (dir == 1) {
      for (int iy=lo[1]; iy<=hi[1]; iy++)
	for (int ix=lo[0]; ix<=hi[0]; ix++)
	  xnew[iy] += tmpmean[mfi](IntVect(ix,iy),0)*phi;
    }
#else
    if (dir == 2) {
      for (int iz=lo[2]; iz<=hi[2]; iz++)
	for (int iy=lo[1]; iy<=hi[1]; iy++)
	  for (int ix=lo[0]; ix<=hi[0]; ix++)
	    xnew[iz] += tmpmean[mfi](IntVect(ix,iy,iz),0)*phi;
    }
#endif
  }

  const int IOProc = ParallelDescriptor::IOProcessorNumber();

  ParallelDescriptor::ReduceRealSum(xnew.dataPtr(),xnew.size(),IOProc);

  if (ParallelDescriptor::IOProcessor())
  {
    Real FL = 0;
    for (int iy = xold.size()-1; iy >=0; iy--)
    //for (int iy = 0; iy < xnew.size(); iy++)
    {
      Real xtmp = xnew[iy+domain.smallEnd()[dir]]/dmn_length;
      Real ct   = (xtmp-xold[iy])/dt;
      FL += ct*dx[dir];
      FLs[iy]  = FL;
      xold[iy] = xtmp;
    }
    flux = FL;  
  }
}

