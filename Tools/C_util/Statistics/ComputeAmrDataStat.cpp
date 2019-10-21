
#include <string>
#include <iostream>
#include <cmath>

#include <ComputeAmrDataStat.H>
#include <WritePlotFile.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

// Determine the volume-averaged mean for an AMR data.
void
ComputeAmrDataMeanVar  (AmrData&     amrData,
			Vector<Real>& mean,
			Vector<Real>& variance,
			int          sComp,
			int          nComp,
			bool         verbose)
{
    std::string oFile, iFileDir, oFileDir;
    
    if (verbose)
    {
      ParmParse pp;
      pp.query("outfile", oFile);
      if (oFile.empty())
	amrex::Abort("You must specify `outfile' if run in verbose mode");
    }
    
    int finestLevel = amrData.FinestLevel();

    mean.clear(); mean.resize(nComp,0.0);
    variance.clear(); variance.resize(nComp,0.0);
    
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
      
      error[iLevel] = new MultiFab(ba, nComp, 0);
      for (int iComp=0; iComp<nComp; ++iComp) {
	MultiFab& data = amrData.GetGrids(iLevel,iComp);
	error[iLevel]->copy(data,0,iComp+sComp,1);
      }

      // Zero out the error covered by fine grid
      long covered_volume = 0;
      if (iLevel != finestLevel) 
      {
	int ref_ratio = amrData.RefRatio()[iLevel];	    
	BoxArray baF  = ::BoxArray(amrData.boxArray(iLevel+1)).coarsen(ref_ratio);
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
	Vector<Real> n1(nComp,0.0), n2(nComp,0.0);

	for (MFIter mfi(*error[iLevel]); mfi.isValid(); ++mfi)
	{
	  FArrayBox& fab = (*error[iLevel])[mfi];
	  const Box& fabbox = mfi.validbox();

	  FArrayBox vwFab(fabbox,nComp);
	  FArrayBox vwFabSqrd(fabbox,nComp);
	  
	  // sum
	  vwFab.copy(fab,0,0,nComp);
	  vwFab.mult(refMult[iLevel]);
	  
	  //sum-squared
	  vwFabSqrd.copy(fab,0,0,nComp);
	  vwFabSqrd.mult(fab,0,0,nComp);
	  vwFabSqrd.mult(refMult[iLevel]);
	  
	  for (int iComp=0; iComp<nComp+sComp; ++iComp)
	  {
	    n1[iComp] += vwFab.norm(fabbox, 1, iComp, 1);
	    n2[iComp] += vwFabSqrd.norm(fabbox, 1, iComp, 1);
	  }
	}
	    
	// Do necessary communication, then blend this level's norms
	//  in with the running global values
	for (int iComp=0; iComp<nComp+sComp; ++iComp)
	{
	  ParallelDescriptor::ReduceRealSum(n1[iComp]);
	  ParallelDescriptor::ReduceRealSum(n2[iComp]);
		
	  mean[iComp] += n1[iComp];
	  variance[iComp] += n2[iComp];
	}
      }
    }
    if (ParallelDescriptor::IOProcessor()) {
      for (int iComp=0; iComp<nComp+sComp; ++iComp)
      {
	mean[iComp] /= total_volume;
	variance[iComp] = variance[iComp]/total_volume - 
	                  mean[iComp]*mean[iComp];
      }
    }

    if (ParallelDescriptor::IOProcessor() && verbose)
    {
      std::cout << "Writing zeroed state to " << oFile << '\n';
      WritePlotFile(error, amrData, oFile, verbose);
    }

    // Clean up memory
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
      delete error[iLevel];
}



// take a list of files and find their mean and variance.
// this is not really correct.
void
ComputeAmrDataList  (AmrData&         amrData,
		     Vector<MultiFab*> mean,
		     Vector<MultiFab*> variance,
		     int              sComp,
		     int              nComp)
{
    std::string oFile, iFileDir, oFileDir;
    
    int finestLevel = amrData.FinestLevel();

    //
    // Compute the sum and sum-squares
    //
    for (int iLevel = finestLevel; iLevel>=0; --iLevel)
    {
        const BoxArray& ba = amrData.boxArray(iLevel);

	MultiFab error(ba, nComp, 0);
	for (int iComp=0; iComp<nComp; ++iComp)
	{
	    MultiFab& data = amrData.GetGrids(iLevel,iComp);
	    error.copy(data,0,iComp+sComp,1);
	}

	for (MFIter mfi(error); mfi.isValid(); ++mfi)
	{
	  FArrayBox& fab = error[mfi];
	  const Box& fabbox = mfi.validbox();
	  
	  FArrayBox vwFab(fabbox,nComp);
	  FArrayBox vwFabSqrd(fabbox,nComp);
	  
	  // sum
	  vwFab.copy(fab,0,0,nComp);
	  
	  //sum-squared
	  vwFabSqrd.copy(fab,0,0,nComp);
	  vwFabSqrd.mult(fab,0,0,nComp);
		
	  (*mean[iLevel])[mfi].plus(vwFab,0,0,1);
	  (*variance[iLevel])[mfi].plus(vwFabSqrd,0,0,1);
	}	    
    }
}

// Determine the mean and variance for an AMR data.
void
ComputeAmrDataMeanVar (AmrData&           amrData,
		       Vector<std::string> cNames,
		       Vector<BoxArray>    bas,
		       Vector<Real>&       mean,
		       Vector<Real>&       variance)
{
    std::string oFile, iFileDir, oFileDir;

    int finestLevel = amrData.FinestLevel();
    int nComp = cNames.size();

    Vector<int> refMult(finestLevel + 1, 1);
    for (int iLevel=finestLevel-1; iLevel>=0; --iLevel)
    {
	int ref_ratio = amrData.RefRatio()[iLevel];
	int area = ref_ratio;
	for (int jLevel=0; jLevel<=iLevel; ++jLevel)
	    refMult[jLevel] *= area;
    }

    //
    // Compute the sum and sum-squares
    //
    long total_volume = 0;
    Vector<int> destFillComps(nComp);
    for (int i=0; i<nComp; ++i)
      destFillComps[i] = i;
    for (int iLevel = finestLevel; iLevel>=0; --iLevel)
    { 
      MultiFab mf(bas[iLevel],nComp,0,Fab_allocate);
      amrData.FillVar(mf,iLevel,cNames,destFillComps);

      // Zero out the error covered by fine grid
      long covered_volume = 0;
      if (iLevel != finestLevel)
      {
	int ref_ratio = amrData.RefRatio()[iLevel];	    
	BoxArray baF  = BoxArray(bas[iLevel+1]).coarsen(ref_ratio);
	for (MFIter mfi(mf); mfi.isValid(); ++mfi)
	{
	  for (int iGrid=0; iGrid<baF.size(); ++iGrid)
	  {
	    Box ovlp = baF[iGrid] & mfi.validbox();
	    if (ovlp.ok())
	    {
	      mf[mfi].setVal(0.0, ovlp, 0, nComp);
	      covered_volume += ovlp.numPts()*refMult[iLevel];
	    }
	  }
	}
	ParallelDescriptor::ReduceLongSum(covered_volume);
      }


      // Compute volume at this level
      Real level_volume = 0.0;
      for (int iGrid=0; iGrid<bas.size(); ++iGrid)
	level_volume += bas[iGrid].numPts()*refMult[iLevel];
      level_volume -= covered_volume;

      if (level_volume > 0.0)
      {
	// Convert volume in numPts to volume in number of fine cells
	total_volume += long(level_volume);
	    
	// Get norms at this level
	Vector<Real> n1(nComp,0.0), n2(nComp,0.0);

	for (MFIter mfi(mf); mfi.isValid(); ++mfi)
	{
	  FArrayBox& fab = mf[mfi];
	  const Box& fabbox = mfi.validbox();

	  FArrayBox vwFab(fabbox,nComp);
	  FArrayBox vwFabSqrd(fabbox,nComp);
		
	  // sum
	  vwFab.copy(fab,0,0,nComp);
	  vwFab.mult(refMult[iLevel]);

	  //sum-squared
	  vwFabSqrd.copy(fab,0,0,nComp);
	  vwFabSqrd.mult(fab,0,0,nComp);
	  vwFabSqrd.mult(refMult[iLevel]);
		
	  for (int iComp=0; iComp<nComp; ++iComp)
	  {
	    n1[iComp] += vwFab.norm(fabbox, 1, iComp, 1);
	    n2[iComp] += vwFabSqrd.norm(fabbox, 1, iComp, 1);
	  }
	}

	// Do necessary communication, then blend this level's norms
	//  in with the running global values
	const int IOProc = ParallelDescriptor::IOProcessorNumber();
	ParallelDescriptor::ReduceRealSum(n1.dataPtr(),nComp,IOProc);
	ParallelDescriptor::ReduceRealSum(n2.dataPtr(),nComp,IOProc);

	if (ParallelDescriptor::IOProcessor()) {
	  for (int iComp=0; iComp<nComp; iComp++) {
	    mean[iComp] += n1[iComp];
	    variance[iComp] += n2[iComp];
	  }
	}
      }
    }

    if (ParallelDescriptor::IOProcessor()) {
      for (int iComp=0; iComp<nComp; ++iComp)
      {
	mean[iComp] /= total_volume;
	variance[iComp] = variance[iComp]/total_volume - 
	  mean[iComp]*mean[iComp];

	std::cout << " comp= " << iComp 
	          << " mean = " << mean[iComp] 
	          << " variance = " << variance[iComp] << std::endl;

      }
    }
}


// Determine the mean and variance of a multifab.
void
ComputeMeanVarMF (MultiFab&          mf,
		  Vector<Real>&       mean,
		  Vector<Real>&       variance)
{
    
    int nComp = mf.nComp();
    BoxArray  ba = mf.boxArray();
    //
    // Compute the sum and sum-squares
    //
    long total_volume = 0;
    for (int iGrid=0; iGrid<ba.size(); ++iGrid)
	total_volume += ba[iGrid].numPts();

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
      FArrayBox& fab = mf[mfi];
      const Box& fabbox = mfi.validbox();

      FArrayBox vwFab(fabbox,nComp);
      FArrayBox vwFabSqrd(fabbox,nComp);
		
      // sum
      vwFab.copy(fab,0,0,nComp);
      
      //sum-squared
      vwFabSqrd.copy(fab,0,0,nComp);
      vwFabSqrd.mult(fab,0,0,nComp);
		
      for (int iComp=0; iComp<nComp; ++iComp)
      {
	mean[iComp] += vwFab.norm(fabbox, 1, iComp, 1);
	variance[iComp] += vwFabSqrd.norm(fabbox, 1, iComp, 1);
      }
    }

    // Do necessary communication, then blend this level's norms
    //  in with the running global values
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealSum(mean.dataPtr(),nComp,IOProc);
    ParallelDescriptor::ReduceRealSum(variance.dataPtr(),nComp,IOProc);

    if (ParallelDescriptor::IOProcessor()) {
      for (int iComp=0; iComp<nComp; ++iComp)
      {
	mean[iComp] /= total_volume;
	variance[iComp] = variance[iComp]/total_volume - 
	  mean[iComp]*mean[iComp];

	std::cout << " comp= " << iComp 
	          << " mean = " << mean[iComp] 
	          << " variance = " << variance[iComp] << std::endl;

      }
    }
}

// determine the pdf of an AMR data.
void
ComputeAmrDataPDF (AmrData&           amrData,
		   Real**             icount,
		   int                nBin,
		   Vector<std::string> cNames,
		   Vector<BoxArray>    bas)
{
    std::string oFile, iFileDir, oFileDir;

    int finestLevel = amrData.FinestLevel();
    int nComp = cNames.size();
    for (int iComp=0; iComp<nComp; iComp++) {
      for (int iBin=0; iBin<nBin; iBin++){
	icount[iComp][iBin] = 0;
      }
    }

    // determine min and max
    Vector<Real> smin(nComp,1.e20), smax(nComp,-1.e20);
    for (int iComp = 0; iComp < nComp; iComp++) 
      amrData.MinMax(amrData.ProbDomain()[finestLevel], 
		     cNames[iComp], finestLevel, smin[iComp], smax[iComp]);
   
    Real scount[nComp][nBin+1];
    for (int iComp=0; iComp < nComp; iComp++) {
      Real ds = (smax[iComp]-smin[iComp])/nBin;
      scount[iComp][0] = smin[iComp];
      for (int iBin=1;iBin <=nBin; iBin++) {
	scount[iComp][iBin] = scount[iComp][iBin-1]+ds;
      }
    }

    Vector<int> refMult(finestLevel + 1, 1);
    for (int iLevel=finestLevel-1; iLevel>=0; --iLevel)
    {
	int ref_ratio = amrData.RefRatio()[iLevel];
	int area = ref_ratio;
	for (int jLevel=0; jLevel<=iLevel; ++jLevel)
	    refMult[jLevel] *= area;
    }

    //
    // Compute the sum and sum-squares
    //
    long total_volume = 0;
    Vector<int> destFillComps(nComp);
    for (int i=0; i<nComp; ++i)
      destFillComps[i] = i;
    for (int iLevel = finestLevel; iLevel>=0; --iLevel)
    { 
      MultiFab mf(bas[iLevel],nComp,0,Fab_allocate);
      amrData.FillVar(mf,iLevel,cNames,destFillComps);

      // Zero out the error covered by fine grid
      long covered_volume = 0;
      if (iLevel != finestLevel)
      {
	int ref_ratio = amrData.RefRatio()[iLevel];	    
	BoxArray baF  = BoxArray(bas[iLevel+1]).coarsen(ref_ratio);
	for (MFIter mfi(mf); mfi.isValid(); ++mfi)
	{
	  for (int iGrid=0; iGrid<baF.size(); ++iGrid)
	  {
	    Box ovlp = baF[iGrid] & mfi.validbox();
	    if (ovlp.ok())
	    {
	      mf[mfi].setVal(0.0, ovlp, 0, nComp);
	      covered_volume += ovlp.numPts()*refMult[iLevel];
	    }
	  }
	}
	ParallelDescriptor::ReduceLongSum(covered_volume);
      }


      // Compute volume at this level
      Real level_volume = 0.0;
      for (int iGrid=0; iGrid<bas.size(); ++iGrid)
	level_volume += bas[iGrid].numPts()*refMult[iLevel];
      level_volume -= covered_volume;

      if (level_volume > 0.0)
      {
	// Convert volume in numPts to volume in number of fine cells
	total_volume += long(level_volume);
	    
	// Get counts at this level
	int n1[nComp][nBin];
	for (int iComp=0; iComp<nComp; iComp++)
	  for (int iBin=0; iBin<nBin; iBin++)
	    n1[iComp][iBin] = 0;

	for (MFIter mfi(mf); mfi.isValid(); ++mfi)
	{
	  const int* lo = mf[mfi].loVect();
	  const int* hi = mf[mfi].hiVect();

	  Real tmp;
#if (BL_SPACEDIM == 2)
	  for (int iy=lo[1]; iy<=hi[1]; iy++) {
	    for (int ix=lo[0]; ix<=hi[0]; ix++) {
	      for (int ic=0; ic<nComp; ic++) {
		tmp = mf[mfi](IntVect(ix,iy),ic);
		if (tmp != 0) {
		  for (int ib=0; ib<nBin; ib++) {
		    if (tmp >= scount[ic][ib] && tmp <  scount[ic][ib+1])
		      n1[ic][ib] += refMult[iLevel];
		  }
		}
	      }
	    }
	  }
#else 
	  for (int iz=lo[1]; iz<=hi[1]; iz++) {
	    for (int iy=lo[1]; iy<=hi[1]; iy++) {
	      for (int ix=lo[0]; ix<=hi[0]; ix++) {
		for (int ic=0; ic<nComp; ic++) {
		  tmp = mf[mfi](IntVect(ix,iy,iz),ic);
		  if (tmp != 0) {
		    for (int ib=0; ib<nBin; ib++) 
		      if (tmp >= scount[ic][ib] && tmp < scount[ic][ib+1])
			n1[ic][ib] += refMult[iLevel];
		  }
		}
	      }
	    }
	  }
#endif
	}

	// Do necessary communication, then blend this level's norms
	//  in with the running global values
	for (int iComp=0; iComp<nComp; iComp++) {
	  for (int iBin=0; iBin<nBin; iBin++) {
	    ParallelDescriptor::ReduceIntSum(n1[iComp][iBin]);
	    icount[iComp][iBin] += n1[iComp][iBin];
	  }
	}
      }
    }
    if (0) {
    for (int iComp=0; iComp<nComp; ++iComp) {
      for (int iBin=0; iBin<nBin; iBin++) {
	icount[iComp][iBin] /= total_volume;    
      }
    }
    }
}


// determine the pdf of a portion of a plotfile
void
ComputeAmrDataPDF (AmrData&     amrData,
		   Real**       icount,
		   int          nBin,
		   Vector<std::string> cNames)
{
    std::string oFile, iFileDir, oFileDir;
        
    int finestLevel = amrData.FinestLevel();
    int nComp = cNames.size();

    for (int iComp=0; iComp<nComp; iComp++) {
      for (int iBin=0; iBin<nBin; iBin++){
	icount[iComp][iBin] = 0;
      }
    }

    // determine min and max
    Vector<Real> smin(nComp,1.e20), smax(nComp,-1.e20);
    Vector<std::string>  VarNames = amrData.PlotVarNames();
    for (int iComp = 0; iComp < nComp; iComp++) 
      amrData.MinMax(amrData.ProbDomain()[finestLevel], 
		     VarNames[iComp], finestLevel, smin[iComp], smax[iComp]);
   
    Real scount[nComp][nBin+1];
    for (int iComp=0; iComp < nComp; iComp++) {
      Real ds = (smax[iComp]-smin[iComp])/nBin;
      scount[iComp][0] = smin[iComp];
      for (int iBin=1;iBin <=nBin; iBin++) {
	scount[iComp][iBin] = scount[iComp][iBin-1]+ds;
      }
    }

    Vector<int> refMult(finestLevel + 1, 1);
    for (int iLevel=finestLevel-1; iLevel>=0; --iLevel)
    {
      int ref_ratio = amrData.RefRatio()[iLevel];
      int area = ref_ratio;
      for (int jLevel=0; jLevel<=iLevel; ++jLevel)
	refMult[jLevel] *= area;
    }

    //
    // Compute the sum and sum-squares
    //
    long total_volume = 0;
    Vector<int> destFillComps(nComp);
    for (int i=0; i<nComp; ++i)
      destFillComps[i] = i;
    
    for (int iLevel = finestLevel; iLevel>=0; --iLevel)
    {
      const BoxArray& ba = amrData.boxArray(iLevel);
      MultiFab mf(ba,nComp,0,Fab_allocate);
      amrData.FillVar(mf,iLevel,cNames,destFillComps);

      // Zero out the error covered by fine grid
      long covered_volume = 0;
      if (iLevel != finestLevel)
      {
	int ref_ratio = amrData.RefRatio()[iLevel];	    
	BoxArray baF = BoxArray(ba[iLevel+1]).coarsen(ref_ratio);
	for (MFIter mfi(mf); mfi.isValid(); ++mfi)
	{
	  for (int iGrid=0; iGrid<baF.size(); ++iGrid)
	  {
	    Box ovlp = baF[iGrid] & mfi.validbox();
	    if (ovlp.ok())
	    {
	      mf[mfi].setVal(0.0, ovlp, 0, nComp);
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
	    
	// Get counts at this level
	int n1[nComp][nBin];
	for (int iComp=0; iComp<nComp; iComp++)
	  for (int iBin=0; iBin<nBin; iBin++)
	    n1[iComp][iBin] = 0;

	for (MFIter mfi(mf); mfi.isValid(); ++mfi)
	{
	  const int* lo = mf[mfi].loVect();
	  const int* hi = mf[mfi].hiVect();
	  
	  Real tmp;
#if (BL_SPACEDIM == 2)
	  for (int iy=lo[1]; iy<=hi[1]; iy++) {
	    for (int ix=lo[0]; ix<=hi[0]; ix++) {
	      for (int ic=0; ic<nComp; ic++) {
		tmp = mf[mfi](IntVect(ix,iy),ic);
		if (tmp != 0) {
		  for (int ib=0; ib<nBin; ib++) {
		    if (tmp >= scount[ic][ib] && 
			tmp <  scount[ic][ib+1])
		      n1[ic][ib] += refMult[iLevel];
		  }
		}
	      }
	    }
	  }
#else 
	  for (int iz=lo[1]; iz<=hi[1]; iz++) {
	    for (int iy=lo[1]; iy<=hi[1]; iy++) {
	      for (int ix=lo[0]; ix<=hi[0]; ix++) {
		for (int ic=0; ic<nComp; ic++) {
		  tmp = mf[mfi](IntVect(ix,iy,iz),ic);
		  if (tmp != 0) {
		    for (int ib=0; ib<nBin; ib++) 
		      if (tmp >= scount[ic][ib] && 
			  tmp < scount[ic][ib+1])
			n1[ic][ib] += refMult[iLevel];
		  }
		}
	      }
	    }
	  }
#endif
	}

	// Do necessary communication, then blend this level's norms
	//  in with the running global values
	for (int iComp=0; iComp<nComp; iComp++) {
	  for (int iBin=0; iBin<nBin; iBin++) {
	    ParallelDescriptor::ReduceIntSum(n1[iComp][iBin]);
	    icount[iComp][iBin] += n1[iComp][iBin];
	  }
	}
      }
    }
    for (int iComp=0; iComp<nComp; ++iComp) {
      for (int iBin=0; iBin<nBin; iBin++) {
	icount[iComp][iBin] /= total_volume;    
      }
    }

}

// Compute the correlation in space
void
ComputeAmrDataVAR (AmrData&           amrData,
		   int                nBin,
		   Vector<std::string> cNames,
		   Vector<Real>        barr,
		   std::string        oFile)
{

    bool do_normalize = false;

    int finestLevel = amrData.FinestLevel();
    int nComp = cNames.size();

    Box domain = amrData.ProbDomain()[0];
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
      const Real dx = amrData.ProbSize()[i] / amrData.ProbDomain()[0].length(i);            
      domain.setSmall(i,std::max(domain.smallEnd()[i], 
				 (int)((bbll[i]-amrData.ProbLo()[i]+.0001*dx)/dx)));
      domain.setBig(i,std::min(domain.bigEnd()[i], 
			       (int)((bbur[i]-amrData.ProbLo()[i]-.0001*dx)/dx)));
    }
    
    for (int i=1; i<=finestLevel; i++)
      domain.refine(amrData.RefRatio()[i]);

    IntVect sm = domain.smallEnd();
    IntVect bg = domain.bigEnd();

    // Initialize covariance, correlation and variogram matrices.
    Real* covx[nComp], * covy[nComp];
    Real* corx[nComp], * cory[nComp];
    Real* varx[nComp], * vary[nComp]; 
    for (int ic=0; ic<nComp; ic++) {
      covx[ic]=new Real[nBin+1];
      corx[ic]=new Real[nBin+1];
      varx[ic]=new Real[nBin+1];
      covy[ic]=new Real[nBin+1];
      cory[ic]=new Real[nBin+1];
      vary[ic]=new Real[nBin+1];
      
      for (int iBin=0; iBin<=nBin; iBin++) { 
	covx[ic][iBin] = 0.0;
	corx[ic][iBin] = 0.0;
	varx[ic][iBin] = 0.0;
	covy[ic][iBin] = 0.0;
	cory[ic][iBin] = 0.0;
	vary[ic][iBin] = 0.0;
      }
    }

    // for 2D only for now
    // In x-direction
    int sz_in_y = bg[1]-sm[1]+1;
    BoxArray bax(sz_in_y);
    IntVect smx = sm;
    IntVect bgx = bg;
    for (int i=0; i<sz_in_y; i++) {
      smx[1] = sm[1]+i;
      bgx[1] = sm[1]+i;
      Box bx(sm,bg);
      bax.set(i,bx);
    }

    Vector<int> destFillComps(nComp);
    for (int i=0; i<nComp; ++i)
      destFillComps[i] = i;
    MultiFab tmpx(bax,nComp,0);
    amrData.FillVar(tmpx,amrData.FinestLevel(),cNames,destFillComps);

    int  Nh[nBin+1];
    Real m0[nComp][nBin+1], mh[nComp][nBin+1];
    Real s0[nComp][nBin+1], sh[nComp][nBin+1];

    for (int iBin=0; iBin<=nBin; iBin++) { 
      Nh[iBin] = 0;
      for (int iComp=0; iComp<nComp; iComp++) {
	m0[iComp][iBin] = 0.0;
	mh[iComp][iBin] = 0.0;
	s0[iComp][iBin] = 0.0;
	sh[iComp][iBin] = 0.0;
      }
    }
    Vector<Real> mean(nComp,0.0), variance(nComp,0.0);
    Real volume = 0.0;
    for (MFIter mfi(tmpx); mfi.isValid(); ++mfi) {

      const int* lo = tmpx[mfi].loVect();
      const int* hi = tmpx[mfi].hiVect();

      FArrayBox& fab = (tmpx)[mfi];
      const Box& fabbox = mfi.validbox();

      FArrayBox vwFab(fabbox,nComp);
      FArrayBox vwFabSqrd(fabbox,nComp);

      // sum and sum-squared
      vwFab.copy(fab,0,0,nComp);
      vwFabSqrd.copy(fab,0,0,nComp);
      
      for (int iComp=0; iComp<nComp; ++iComp)
      {
	mean[iComp] += vwFab.norm(fabbox, 1, iComp, 1);
	variance[iComp] += vwFabSqrd.norm(fabbox, 1, iComp, 1);
      }

      volume += fabbox.numPts();
      
      for (int i=0; i<=nBin; i++) {

	Real tmp, tmpp1;

	for (int iy=lo[1]; iy<=hi[1]; iy++) {
	  for (int ix=lo[0]; ix<hi[0]-i; ix++) {

	    Nh[i] = Nh[i] + 1;
	    for (int ic=0; ic<nComp; ic++) {
	      tmp = tmpx[mfi](IntVect(ix,iy),ic);
	      tmpp1 = tmpx[mfi](IntVect(ix+i,iy),ic);
	      
	      covx[ic][i] = covx[ic][i]+ tmp*tmpp1;
	      varx[ic][i] = varx[ic][i]+ (tmpp1 - tmp)*(tmpp1 - tmp);
	      m0[ic][i]   = m0[ic][i] + tmp;
	      mh[ic][i]   = mh[ic][i] + tmpp1;
	      s0[ic][i]   = s0[ic][i] + tmp*tmp;
	      sh[ic][i]   = sh[ic][i] + tmpp1*tmpp1;

	    }
	  }
	}
      }
    }

    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    ParallelDescriptor::ReduceIntSum(Nh,nBin+1,IOProc);
    ParallelDescriptor::ReduceRealSum(volume);
    for (int iComp=0; iComp<nComp; iComp++) {
      ParallelDescriptor::ReduceRealSum(covx[iComp],nBin+1,IOProc);
      ParallelDescriptor::ReduceRealSum(m0[iComp],nBin+1,IOProc);
      ParallelDescriptor::ReduceRealSum(mh[iComp],nBin+1,IOProc);
      ParallelDescriptor::ReduceRealSum(s0[iComp],nBin+1,IOProc);
      ParallelDescriptor::ReduceRealSum(sh[iComp],nBin+1,IOProc);

      ParallelDescriptor::ReduceRealSum(mean[iComp]);
      ParallelDescriptor::ReduceRealSum(variance[iComp]);

    }



    if (ParallelDescriptor::IOProcessor()) {
      for (int iComp=0; iComp<nComp; iComp++) {
	
	mean[iComp] /= volume;
	variance[iComp] = variance[iComp]/volume - mean[iComp]*mean[iComp];

	for (int iBin=0; iBin<=nBin; iBin++) { 
	  m0[iComp][iBin] /= Nh[iBin];
	  mh[iComp][iBin] /= Nh[iBin];
	  s0[iComp][iBin] = s0[iComp][iBin]/Nh[iBin] - m0[iComp][iBin];
	  sh[iComp][iBin] = sh[iComp][iBin]/Nh[iBin] - mh[iComp][iBin];

	  covx[iComp][iBin] = covx[iComp][iBin]/Nh[iBin] -
	    m0[iComp][iBin]*mh[iComp][iBin];
	  corx[iComp][iBin] = covx[iComp][iBin]/sqrt(s0[iComp][iBin]*sh[iComp][iBin]);
	  varx[iComp][iBin] = 0.5*varx[iComp][iBin]/Nh[iBin];

	  // normalize with mean
	  if (do_normalize) {
	    covx[iComp][iBin] /= mean[iComp];
	    corx[iComp][iBin] /= mean[iComp];
	    varx[iComp][iBin] /= mean[iComp];
	  }
	}
      }
    }

    bax.clear();
    tmpx.clear();

    // in y direction
    int sz_in_x = bg[0]-sm[0]+1;
    BoxArray bay(sz_in_x);
    smx = sm;
    bgx = bg;
    for (int i=0; i<sz_in_x; i++) {
      smx[0] = sm[0]+i;
      bgx[0] = sm[0]+i;
      Box bx(sm,bg);
      bay.set(i,bx);
    }

    MultiFab tmpy(bay,nComp,0);
    amrData.FillVar(tmpy,amrData.FinestLevel(),cNames,destFillComps);

    for (int iBin=0; iBin<=nBin; iBin++) { 
      Nh[iBin] = 0;
      for (int iComp=0; iComp<nComp; iComp++) {
	m0[iComp][iBin] = 0.0;
	mh[iComp][iBin] = 0.0;
	s0[iComp][iBin] = 0.0;
	sh[iComp][iBin] = 0.0;
      }
    }

    for (MFIter mfi(tmpy); mfi.isValid(); ++mfi) {

      const int* lo = tmpy[mfi].loVect();
      const int* hi = tmpy[mfi].hiVect();

      for (int i=0; i<=nBin; i++) {

	Real tmp, tmpp1;

	for (int ix=lo[0]; ix<=hi[0]; ix++) {
	  for (int iy=lo[1]; iy<hi[1]-i; iy++) {

	    Nh[i] = Nh[i] + 1;
	    for (int ic=0; ic<nComp; ic++) {
	      tmp = tmpy[mfi](IntVect(ix,iy),ic);
	      tmpp1 = tmpy[mfi](IntVect(ix+i,iy),ic);
	      
	      covy[ic][i] = covy[ic][i]+ tmp*tmpp1;
	      vary[ic][i] = vary[ic][i]+ (tmpp1 - tmp)*(tmpp1 - tmp);
	      m0[ic][i]   = m0[ic][i] + tmp;
	      mh[ic][i]   = mh[ic][i] + tmpp1;
	      s0[ic][i]   = s0[ic][i] + tmp*tmp;
	      sh[ic][i]   = sh[ic][i] + tmpp1*tmpp1;

	    }
	  }
	}
      }
    }

    bay.clear();
    tmpy.clear();

    ParallelDescriptor::ReduceIntSum(Nh,nBin+1,IOProc);
    for (int iComp=0; iComp<nComp; iComp++) {
      ParallelDescriptor::ReduceRealSum(covy[iComp],nBin+1,IOProc);
      ParallelDescriptor::ReduceRealSum(m0[iComp],nBin+1,IOProc);
      ParallelDescriptor::ReduceRealSum(mh[iComp],nBin+1,IOProc);
      ParallelDescriptor::ReduceRealSum(s0[iComp],nBin+1,IOProc);
      ParallelDescriptor::ReduceRealSum(sh[iComp],nBin+1,IOProc);
    }

    if (ParallelDescriptor::IOProcessor()) {
      for (int iComp=0; iComp<nComp; iComp++) {
	for (int iBin=0; iBin<=nBin; iBin++) { 
	  m0[iComp][iBin] /= Nh[iBin];
	  mh[iComp][iBin] /= Nh[iBin];
	  s0[iComp][iBin] = s0[iComp][iBin]/Nh[iBin] - m0[iComp][iBin];
	  sh[iComp][iBin] = sh[iComp][iBin]/Nh[iBin] - mh[iComp][iBin];

	  covy[iComp][iBin] = covy[iComp][iBin]/Nh[iBin] -
	    m0[iComp][iBin]*mh[iComp][iBin];
	  cory[iComp][iBin] = covy[iComp][iBin]/sqrt(s0[iComp][iBin]*sh[iComp][iBin]);
	  vary[iComp][iBin] = 0.5*vary[iComp][iBin]/Nh[iBin];
	  
	  // normalize with mean
	  if (do_normalize) {
	    covy[iComp][iBin] /= mean[iComp];
	    cory[iComp][iBin] /= mean[iComp];
	    vary[iComp][iBin] /= mean[iComp];
	  }
	}
      }
    }


    if (ParallelDescriptor::IOProcessor()) {
      
      std::ofstream outputFile;
      outputFile.open(oFile.c_str(),std::ios::out);
      if (outputFile.fail()) amrex::Abort("Output file cannot be opened");
      for (int ic=0; ic<nComp; ic++) {
	for (int ib=0; ib<=nBin; ib++) 		
	  outputFile << covx[ic][ib] << " ";
	outputFile << "\n";
	for (int ib=0; ib<=nBin; ib++) 		
	  outputFile << corx[ic][ib] << " ";
	outputFile << "\n";
	for (int ib=0; ib<=nBin; ib++) 		
	  outputFile << varx[ic][ib] << " ";
	outputFile << "\n";
	for (int ib=0; ib<=nBin; ib++) 		
	  outputFile << covy[ic][ib] << " ";
	outputFile << "\n";
	for (int ib=0; ib<=nBin; ib++) 		
	  outputFile << cory[ic][ib] << " ";
	outputFile << "\n";
	for (int ib=0; ib<=nBin; ib++) 		
	  outputFile << vary[ic][ib] << " ";
	outputFile << "\n";


      }
      outputFile.close();
    }

    for (int ic=0;ic<nComp;ic++) {
      delete [] covx[ic];
      delete [] corx[ic];
      delete [] varx[ic];
      delete [] covy[ic];
      delete [] cory[ic];
      delete [] vary[ic];
    }
}

// Compute variograms based on GSLIB for a uniform grid.
// If grid is not uniform, values will be interpolated to
// the grid at the finest level.
void
VariogramUniform (AmrData&             amrData,
		  Vector<std::string>   cNames,
		  Vector<Real>          barr,
		  Vector< Vector<int> >  ivoption,
		  int                  nlag,
		  int                  isill,
		  Vector<Real>          sills,
		  std::string          oFile)
{

    bool firsttime = true;

    int finestLevel = amrData.FinestLevel();
    int nComp = cNames.size();

    Box domain = amrData.ProbDomain()[0];
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
    Vector<Real> dx(BL_SPACEDIM);
    for (int i=0; i<BL_SPACEDIM; ++i) {
      dx[i] = amrData.ProbSize()[i]/
	amrData.ProbDomain()[0].length(i);            
      domain.setSmall(i,std::max(domain.smallEnd()[i], 
	(int)((bbll[i]-amrData.ProbLo()[i]+.0001*dx[i])/dx[i])));
      domain.setBig(i,std::min(domain.bigEnd()[i], 
	(int)((bbur[i]-amrData.ProbLo()[i]-.0001*dx[i])/dx[i])));
    }
    
    for (int i=1; i<=finestLevel; i++)
      domain.refine(amrData.RefRatio()[i]);

    IntVect sm = domain.smallEnd();
    IntVect bg = domain.bigEnd();

    BoxArray ba(domain);
    Vector<int> destFillComps(nComp);
    for (int i=0; i<nComp; ++i) destFillComps[i] = i;
    MultiFab mf(ba,nComp,0);
    amrData.FillVar(mf,amrData.FinestLevel(),cNames,destFillComps);

    VariogramUniformMFG(mf,dx,sm,bg,ivoption,nlag,isill,sills,oFile);

}

// Compute variograms based on GSLIB for a uniform grid.
// Cross correlation is done with respect to a multifab.
void
VariogramCross(AmrData&             amrData,
	       Vector<std::string>   cNames,
	       MultiFab&            mf,
	       Vector<Real>          barr,
	       Vector< Vector<int> >  ivoption,
	       int                  nlag,
	       int                  isill,
	       Vector<Real>          sills,
	       std::string          oFile)
{
    int finestLevel = amrData.FinestLevel();
    int nComp = cNames.size();

    Box domain = amrData.ProbDomain()[0];
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
    Vector<Real> dx(BL_SPACEDIM);
    for (int i=0; i<BL_SPACEDIM; ++i) {
      dx[i] = amrData.ProbSize()[i]/
	amrData.ProbDomain()[0].length(i);            
      domain.setSmall(i,std::max(domain.smallEnd()[i], 
	(int)((bbll[i]-amrData.ProbLo()[i]+.0001*dx[i])/dx[i])));
      domain.setBig(i,std::min(domain.bigEnd()[i], 
	(int)((bbur[i]-amrData.ProbLo()[i]-.0001*dx[i])/dx[i])));
    }
    
    for (int i=1; i<=finestLevel; i++)
      domain.refine(amrData.RefRatio()[i]);

    IntVect sm = domain.smallEnd();
    IntVect bg = domain.bigEnd();

    BoxArray ba(domain);
    Vector<int> destFillComps(nComp);
    for (int i=0; i<nComp; ++i) destFillComps[i] = i;
    MultiFab mf_amr(ba,nComp+1,0);
    mf_amr.setVal(0.);
    amrData.FillVar(mf_amr,amrData.FinestLevel(),cNames,destFillComps);
    mf_amr.copy(mf,0,nComp,1);

    VariogramUniformMFG(mf_amr,dx,sm,bg,ivoption,nlag,isill,sills,oFile);

}

void
VariogramUniformMF (const MultiFab&      mf,
		    Vector<Real>          dx,
		    Vector< Vector<int> >  ivoption,
		    int                  nlag,
		    int                  isill,
		    Vector<Real>          sills,
		    std::string          oFile)
{

  Vector<int> domloc(BL_SPACEDIM), domhic(BL_SPACEDIM);
  for (MFIter mfi(mf); mfi.isValid(); ++mfi)
  {
    const int* k_lo  = mf[mfi].loVect();
    const int* k_hi  = mf[mfi].hiVect();
    
    for (int i=0;i<BL_SPACEDIM; i++) {
      domloc[i] = k_lo[i];
      domhic[i] = k_hi[i];
    }
  }
  ParallelDescriptor::ReduceIntMin(domloc.dataPtr(),BL_SPACEDIM);
  ParallelDescriptor::ReduceIntMax(domhic.dataPtr(),BL_SPACEDIM);
  
  IntVect sm(domloc.dataPtr());
  IntVect bg(domhic.dataPtr());

  VariogramUniformMFG (mf,dx,sm,bg,ivoption,nlag,isill,sills,oFile);
  
}

void
VariogramUniformMFG (const MultiFab&      mf,
		     Vector<Real>          dx,
		     IntVect              sm,
		     IntVect              bg,
		     Vector< Vector<int> >  ivoption,
		     int                  nlag,
		     int                  isill,
		     Vector<Real>          sills,
		     std::string          oFile)
{

    bool firsttime = true;

    for (int dir=0; dir<BL_SPACEDIM; dir++) 
    {
      // Set up BoxArray
      BoxArray ba;
      if (dir == 0) {
	int sz_in = bg[1] - sm[1] + 1;
	ba.resize(sz_in);
	IntVect smx = sm;
	IntVect bgx = bg;
	for (int i=0; i<sz_in; i++) {
	  smx[1] = sm[1] + i;
	  bgx[1] = sm[1] + i;
	  Box bx(smx,bgx);
	  ba.set(i,bx);
	}
      }
      
      else if (dir == 1) {
	int sz_in = bg[0] - sm[0] + 1;
	ba.resize(sz_in);
	IntVect smx = sm;
	IntVect bgx = bg;
	for (int i=0; i<sz_in; i++) {
	  smx[0] = sm[0] + i;
	  bgx[0] = sm[0] + i;
	  Box bx(smx,bgx);
	  ba.set(i,bx);
	}
      }

      // Fill tmpx with the data
      MultiFab tmpx(ba,mf.nComp(),0);
      tmpx.copy(mf);

      int nvarg = ivoption.size();
      for (int iv=0; iv<nvarg; iv++) 
      {
	int ivtail = ivoption[iv][0];
	int ivhead = ivoption[iv][1];
	int ivtype = ivoption[iv][2];

      	Vector<Real> np(nlag,0.0);
	Vector<Real> gam(nlag,0.0);
	Vector<Real> hm(nlag,0.0);
	Vector<Real> tm(nlag,0.0);
	Vector<Real> hv(nlag,0.0);
	Vector<Real> tv(nlag,0.0);

	for (MFIter mfi(tmpx); mfi.isValid(); ++mfi) {
	  
	  const int* lo = tmpx[mfi].loVect();
	  const int* hi = tmpx[mfi].hiVect();

	  Vector<int> lod(BL_SPACEDIM),hid(BL_SPACEDIM);
	  if (dir == 0) {
	    lod[0] = lo[0]; 
	    lod[1] = lo[1];
	    hid[0] = hi[0];
	    hid[1] = hi[1];
	  }
	  else if (dir == 1){
	    lod[0] = lo[1]; 
	    lod[1] = lo[0];
	    hid[0] = hi[1];
	    hid[1] = hi[0];
	  }

	  for (int i=0; i<nlag; i++) {
	    for (int iy=lod[1]; iy<=hid[1]; iy++) {
	      for (int ix=lod[0]; ix<hid[0]-i-1; ix++) {

		Real vrt, vrh;
		Real vrtpr, vrhpr;

		if (dir == 0) {
		  vrt = tmpx[mfi](IntVect(ix,iy),ivhead);
		  vrh = tmpx[mfi](IntVect(ix+i,iy),ivtail);
		
		  if (ivtype == 2) {
		    vrhpr = tmpx[mfi](IntVect(ix,iy),ivtail);
		    vrtpr = tmpx[mfi](IntVect(ix+i,iy),ivhead);
		  }
		}
		else if (dir == 1) {
		  vrt = tmpx[mfi](IntVect(iy,ix),ivhead);
		  vrh = tmpx[mfi](IntVect(iy,ix+i),ivtail);
		
		  if (ivtype == 2) {
		    vrhpr = tmpx[mfi](IntVect(iy,ix),ivtail);
		    vrtpr = tmpx[mfi](IntVect(iy,ix+i),ivhead);
		  }
		}

		np[i] = np[i] + 1.0;
		tm[i] = tm[i] + vrt;
		hm[i] = hm[i] + vrh;
		
		if (ivtype == 1 | ivtype == 9)
		  gam[i] = gam[i] + (vrh-vrt)*(vrh-vrt);
		else if (ivtype == 2)
		  gam[i] = gam[i] + (vrhpr-vrh)*(vrt-vrtpr);
		else if (ivtype == 3)
		  gam[i] = gam[i] + vrh*vrt;
		else if (ivtype == 4) {
		  gam[i] = gam[i] + vrh*vrt;
		  hv[i]  = hv[i]  + vrh*vrh;
		  tv[i]  = tv[i]  + vrt*vrt;
		}
		else if (ivtype == 5)
		  gam[i] = gam[i] + (vrh-vrt)*(vrh-vrt);
		else if (ivtype == 6) {
		  if ((vrt+vrh) < 1.e-20) {
		    np[i] = np[i] - 1;
		    tm[i] = tm[i] - vrt;
		    hm[i] = hm[i] - vrh;
		  }
		  else {
		    Real tempvar = 2.0*(vrt-vrh)/(vrt+vrh);
		    gam[i] = gam[i] + tempvar*tempvar;
		  }
		}
		else if (ivtype == 7) {
		  if (vrt < 1.e-20 || vrh < 1.e-20) {
		    np[i] = np[i] - 1;
		    tm[i] = tm[i] - vrt;
		    hm[i] = hm[i] - vrh;
		  }
		  else {
		    Real tempvar = log(vrt) - log(vrh);
		    gam[i] = gam[i] + tempvar*tempvar;
		  }
		}
		else if (ivtype == 8)
		  gam[i] = gam[i] + std::abs(vrt-vrh);
 
	      }
	    }
	  }
	}
	const int IOProc = ParallelDescriptor::IOProcessorNumber();
	ParallelDescriptor::ReduceRealSum(np.dataPtr(),nlag,IOProc);
	ParallelDescriptor::ReduceRealSum(gam.dataPtr(),nlag,IOProc);
	ParallelDescriptor::ReduceRealSum(hm.dataPtr(),nlag,IOProc);
	ParallelDescriptor::ReduceRealSum(tm.dataPtr(),nlag,IOProc);
	ParallelDescriptor::ReduceRealSum(hv.dataPtr(),nlag,IOProc);
	ParallelDescriptor::ReduceRealSum(tv.dataPtr(),nlag,IOProc);

	if (ParallelDescriptor::IOProcessor()) {
	  for (int i=0; i<nlag; i++) {
	    gam[i] = gam[i]/np[i];
	    hm[i]  = hm[i]/np[i];
	    tm[i]  = hm[i]/np[i];
	    hv[i]  = hm[i]/np[i];
	    tv[i]  = hm[i]/np[i];
	    
	    if (isill == 1 && ivtail == ivhead) {
	      if ((ivtype == 1 || ivtype == 9) && sills[ivtail] > 0.0) {
		gam[i] = gam[i]/sills[ivtail];
	      }
	    }
	    
	    if (ivtype == 1 || ivtype == 2 || ivtype == 6) {
	      gam[i] = 0.5*gam[i];
	    }
	    else if (ivtype == 3) 
	      gam[i] = gam[i] - hm[i]*hm[i];
	    else if (ivtype == 4) {
	      hv[i] = hv[i] - hm[i]*hm[i];
	      if (hv[i] <= 0.0) hv[i] = 0.0;
	      hv[i] = sqrt(hv[i]);
	      tv[i] = tv[i] - tm[i]*tm[i];
	      if (tv[i] <= 0.0) tv[i] = 0.0;
	      tv[i] = sqrt(tv[i]);
	      if (hv[i]*tv[i] < 1.e-20) 
		gam[i] = 0.0;
	      else
		gam[i] = (gam[i] - hm[i]*tm[i])/(hv[i]*tv[i]);
	      hv[i] = hv[i]*hv[i];
	      tv[i] = tv[i]*tv[i];
	    }
	    else if (ivtype == 5) {
	      Real htave = 0.5*hm[i]*tm[i];
	      htave *= htave;
	      if (htave < 1.e-20) 
		gam[i] = 0.0;
	      else
		gam[i] = gam[i]/htave;
	    }
	  }
	  std::cout << ivtail  << " " << ivhead << " " 
		    << ivtype  << " " << dir << " " 
		    << dx[dir] << " " << nlag << std::endl;
	  // write out to datafile
	  std::ofstream outputFile;
	  if (firsttime) {
	    outputFile.open(oFile.c_str(),std::ios::out);
	    if (outputFile.fail()) 
	      amrex::Abort("Output file cannot be opened");
	    firsttime = false;
	  }
	  else 
	    outputFile.open(oFile.c_str(),std::ios::app);

	  // ivtail,ivhead,ivtype,dir,nlag
	  outputFile << ivtail  << " " << ivhead << " " 
		     << ivtype  << " " << dir << " " 
		     << dx[dir] << " " << nlag << std::endl;
	  
	  for (int i=0; i<nlag; i++)
	    outputFile << gam[i] << " ";
	  outputFile << "\n";

	  outputFile.close();

	}
      }
    }

    std::cout << "DONE ...\n";
}



// fine solution - coarse solution on the grid finest level specified.
void
TakeDifferenceFine(AmrData&             amrDataf,
		   AmrData&             amrDatac,
		   Vector<std::string>   cNames,
		   Vector<Real>          barr,
		   std::string          oFile)
{
    int nComp = cNames.size();

    Box domainf = amrDataf.ProbDomain()[0];
    Box domainc = amrDatac.ProbDomain()[0];
    for (int i=0;i<BL_SPACEDIM;i++) {
      if (domainf.smallEnd()[i] != domainc.smallEnd()[i])
	amrex::Abort("Domains of coarse and fine do not match!");
    }

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
    Vector <int> rratio(BL_SPACEDIM);
    Vector<Real> dxc(BL_SPACEDIM),dxf(BL_SPACEDIM);
    for (int i=0; i<BL_SPACEDIM; ++i) {
      dxf[i] = amrDataf.ProbSize()[i]/
	amrDataf.ProbDomain()[0].length(i);            
      domainf.setSmall(i,std::max(domainf.smallEnd()[i], 
	(int)((bbll[i]-amrDataf.ProbLo()[i]+.0001*dxf[i])/dxf[i])));
      domainf.setBig(i,std::min(domainf.bigEnd()[i], 
	(int)((bbur[i]-amrDataf.ProbLo()[i]-.0001*dxf[i])/dxf[i])));

      dxc[i] = amrDatac.ProbSize()[i]/
	amrDatac.ProbDomain()[0].length(i);            
      domainc.setSmall(i,std::max(domainc.smallEnd()[i], 
	(int)((bbll[i]-amrDatac.ProbLo()[i]+.0001*dxc[i])/dxc[i])));
      domainc.setBig(i,std::min(domainc.bigEnd()[i], 
	(int)((bbur[i]-amrDatac.ProbLo()[i]-.0001*dxc[i])/dxc[i])));

      rratio[i] = (int)((dxc[i]/(amrDatac.FinestLevel()+1))/
		        (dxf[i]/(amrDataf.FinestLevel()+1)));
    }

    for (int i=0;i<BL_SPACEDIM;i++) {
      if (domainf.smallEnd()[i] != domainc.smallEnd()[i])
	amrex::Abort("Domains of coarse and fine do not match!");
    }

    for (int i=1; i<=amrDataf.FinestLevel(); i++)
      domainf.refine(amrDataf.RefRatio()[i]);
    for (int i=1; i<=amrDatac.FinestLevel(); i++)
      domainc.refine(amrDatac.RefRatio()[i]);

    BoxArray baf(domainf), bac(domainc);
    
    // Fill tmpx with the data
    Vector<int> destFillComps(nComp);
    Vector<std::string> destNames(nComp);
    for (int i=0; i<nComp; ++i) 
      destFillComps[i] = i;
    MultiFab tmpc(bac,nComp,0), tmpf(baf,nComp,0);

    amrDatac.FillVar(tmpc,amrDatac.FinestLevel(),cNames,destFillComps);
    amrDataf.FillVar(tmpf,amrDataf.FinestLevel(),cNames,destFillComps);

    for (MFIter mfi(tmpc); mfi.isValid(); ++mfi) {
	  
      const int* lo = tmpc[mfi].loVect();
      const int* hi = tmpc[mfi].hiVect();

      for (int joff=0; joff<rratio[1]; joff++) {
	for (int jc=lo[1]; jc<=hi[1]; jc++) {
	  int j = jc*rratio[1] +joff;
	  for (int ioff=0; ioff<rratio[0]; ioff++) {
	    for (int ic=lo[0]; ic<=hi[0]; ic++) {  
	      int i = ic*rratio[0] + ioff;
	      for (int n=0; n<nComp; n++) {
		tmpf[mfi](IntVect(i,j),n) -= 
		  tmpc[mfi](IntVect(ic,jc),n);
	      }
	    }
	  }
	}
      }
    }
    VisMF::Write(tmpf,oFile);
}

// fine solution - coarse solution on the coarse finest level specified.
void
TakeDifferenceCrse(AmrData&             amrDataf,
		   AmrData&             amrDatac,
		   Vector<std::string>   cNames,
		   Vector<Real>          barr,
		   std::string          oFile)
{
    int nComp = cNames.size();

    Box domainf = amrDataf.ProbDomain()[0];
    Box domainc = amrDatac.ProbDomain()[0];
    for (int i=0;i<BL_SPACEDIM;i++) {
      if (domainf.smallEnd()[i] != domainc.smallEnd()[i])
	amrex::Abort("Domains of coarse and fine do not match!");
    }

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
    Vector <int> rratio(BL_SPACEDIM);
    Vector<Real> dxc(BL_SPACEDIM),dxf(BL_SPACEDIM);
    for (int i=0; i<BL_SPACEDIM; ++i) {
      dxf[i] = amrDataf.ProbSize()[i]/
	amrDataf.ProbDomain()[0].length(i);            
      domainf.setSmall(i,std::max(domainf.smallEnd()[i], 
	(int)((bbll[i]-amrDataf.ProbLo()[i]+.0001*dxf[i])/dxf[i])));
      domainf.setBig(i,std::min(domainf.bigEnd()[i], 
	(int)((bbur[i]-amrDataf.ProbLo()[i]-.0001*dxf[i])/dxf[i])));

      dxc[i] = amrDatac.ProbSize()[i]/
	amrDatac.ProbDomain()[0].length(i);            
      domainc.setSmall(i,std::max(domainc.smallEnd()[i], 
	(int)((bbll[i]-amrDatac.ProbLo()[i]+.0001*dxc[i])/dxc[i])));
      domainc.setBig(i,std::min(domainc.bigEnd()[i], 
	(int)((bbur[i]-amrDatac.ProbLo()[i]-.0001*dxc[i])/dxc[i])));

      rratio[i] = (int)((dxc[i]/(amrDatac.FinestLevel()+1))/
		        (dxf[i]/(amrDataf.FinestLevel()+1)));
    }

    for (int i=0;i<BL_SPACEDIM;i++) {
      if (domainf.smallEnd()[i] != domainc.smallEnd()[i])
	amrex::Abort("Domains of coarse and fine do not match!");
    }

    for (int i=1; i<=amrDataf.FinestLevel(); i++)
      domainf.refine(amrDataf.RefRatio()[i]);
    for (int i=1; i<=amrDatac.FinestLevel(); i++)
      domainc.refine(amrDatac.RefRatio()[i]);

    BoxArray baf(domainf), bac(domainc);
    
    // Fill tmpx with the data
    Vector<int> destFillComps(nComp);
    Vector<std::string> destNames(nComp);
    for (int i=0; i<nComp; ++i) 
      destFillComps[i] = i;
    MultiFab tmpc(bac,nComp,0), tmpf(baf,nComp,0);

    amrDatac.FillVar(tmpc,amrDatac.FinestLevel(),cNames,destFillComps);
    amrDataf.FillVar(tmpf,amrDataf.FinestLevel(),cNames,destFillComps);

    Real fv = 1.0/(rratio[0]*rratio[1]);

    for (MFIter mfi(tmpc); mfi.isValid(); ++mfi) {
	  
      const int* lo = tmpc[mfi].loVect();
      const int* hi = tmpc[mfi].hiVect();

      for (int joff=0; joff<rratio[1]; joff++) {
	for (int jc=lo[1]; jc<=hi[1]; jc++) {
	  int j = jc*rratio[1] +joff;
	  for (int ioff=0; ioff<rratio[0]; ioff++) {
	    for (int ic=lo[0]; ic<=hi[0]; ic++) {  
	      int i = ic*rratio[0] + ioff;
	      for (int n=0; n<nComp; n++) {
		tmpc[mfi](IntVect(ic,jc),n) -= fv* 
		  tmpf[mfi](IntVect(i,j),n);
	      }
	    }
	  }
	}
      }

      tmpc[mfi].mult(-1.0);
    }

    VisMF::Write(tmpc,oFile);
}

// fine solution - coarse solution on the grid finest level specified.
void
TakeDifferenceSum(AmrData&             amrDataf,
		  AmrData&             amrDatac,
		  Vector<std::string>   cNames,
		  Vector<Real>          barr,
		  std::string          oFile)
{
    int nComp = cNames.size();

    Box domainf = amrDataf.ProbDomain()[0];
    Box domainc = amrDatac.ProbDomain()[0];
    for (int i=0;i<BL_SPACEDIM;i++) {
      if (domainf.smallEnd()[i] != domainc.smallEnd()[i])
	amrex::Abort("Domains of coarse and fine do not match!");
    }

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
    Vector <int> rratio(BL_SPACEDIM);
    Vector<Real> dxc(BL_SPACEDIM),dxf(BL_SPACEDIM);
    for (int i=0; i<BL_SPACEDIM; ++i) {
      dxf[i] = amrDataf.ProbSize()[i]/
	amrDataf.ProbDomain()[0].length(i);            
      domainf.setSmall(i,std::max(domainf.smallEnd()[i], 
	(int)((bbll[i]-amrDataf.ProbLo()[i]+.0001*dxf[i])/dxf[i])));
      domainf.setBig(i,std::min(domainf.bigEnd()[i], 
	(int)((bbur[i]-amrDataf.ProbLo()[i]-.0001*dxf[i])/dxf[i])));

      dxc[i] = amrDatac.ProbSize()[i]/
	amrDatac.ProbDomain()[0].length(i);            
      domainc.setSmall(i,std::max(domainc.smallEnd()[i], 
	(int)((bbll[i]-amrDatac.ProbLo()[i]+.0001*dxc[i])/dxc[i])));
      domainc.setBig(i,std::min(domainc.bigEnd()[i], 
	(int)((bbur[i]-amrDatac.ProbLo()[i]-.0001*dxc[i])/dxc[i])));

      rratio[i] = (int)((dxc[i]/(amrDatac.FinestLevel()+1))/
		        (dxf[i]/(amrDataf.FinestLevel()+1)));
    }

    for (int i=0;i<BL_SPACEDIM;i++) {
      if (domainf.smallEnd()[i] != domainc.smallEnd()[i])
	amrex::Abort("Domains of coarse and fine do not match!");
    }

    for (int i=1; i<=amrDataf.FinestLevel(); i++)
      domainf.refine(amrDataf.RefRatio()[i]);
    for (int i=1; i<=amrDatac.FinestLevel(); i++)
      domainc.refine(amrDatac.RefRatio()[i]);

    BoxArray baf(domainf), bac(domainc);
    
    // Fill tmpx with the data
    Vector<int> destFillComps(nComp);
    Vector<std::string> destNames(nComp);
    for (int i=0; i<nComp; ++i) 
      destFillComps[i] = i;
    MultiFab tmpc(bac,nComp,0), tmpf(baf,nComp,0);
    MultiFab tmpfl(bac,2,0);
    tmpfl.setVal(0.);
    amrDatac.FillVar(tmpc,amrDatac.FinestLevel(),cNames,destFillComps);
    amrDataf.FillVar(tmpf,amrDataf.FinestLevel(),cNames,destFillComps);

    for (MFIter mfi(tmpc); mfi.isValid(); ++mfi) {
	  
      const int* lo = tmpc[mfi].loVect();
      const int* hi = tmpc[mfi].hiVect();

      for (int joff=0; joff<rratio[1]; joff++) {
	for (int jc=lo[1]; jc<=hi[1]; jc++) {
	  int j = jc*rratio[1] +joff;
	  for (int ic=lo[0]; ic<=hi[0]; ic++) {  
	    int i = ic*rratio[0];
	    tmpfl[mfi](IntVect(ic,jc),0) += 
	      tmpf[mfi](IntVect(i,j),0);
	  }
	}
      }

      
      for (int jc=lo[1]; jc<=hi[1]; jc++) {
	int j = jc*rratio[1];
	for (int ioff=0; ioff<rratio[0]; ioff++) {
	  for (int ic=lo[0]; ic<=hi[0]; ic++) {  
	    int i = ic*rratio[0] + ioff;
	    for (int n=0; n<nComp; n++) {
	      tmpfl[mfi](IntVect(ic,jc),1) += 
		tmpf[mfi](IntVect(i,j),1);
	    }
	  }
	}
      }
    }

    tmpfl.minus(tmpc,0,2,0);  
    VisMF::Write(tmpfl,oFile);
}
// fine solution - coarse solution on the coarse finest level specified.
void
TakeDifferenceMean(AmrData&             amrDataf,
		   Vector<std::string>   cNames,
		   Vector<Real>          barr,
		   Vector<int>           rratio,
		   std::string          oFile)
{
    int nComp = cNames.size();

    Box domainf = amrDataf.ProbDomain()[0];
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
    Vector<Real> dxc(BL_SPACEDIM),dxf(BL_SPACEDIM);
    for (int i=0; i<BL_SPACEDIM; ++i) {
      dxf[i] = amrDataf.ProbSize()[i]/
	amrDataf.ProbDomain()[0].length(i);            
      domainf.setSmall(i,std::max(domainf.smallEnd()[i], 
	(int)((bbll[i]-amrDataf.ProbLo()[i]+.0001*dxf[i])/dxf[i])));
      domainf.setBig(i,std::min(domainf.bigEnd()[i], 
	(int)((bbur[i]-amrDataf.ProbLo()[i]-.0001*dxf[i])/dxf[i])));
    }

    for (int i=1; i<=amrDataf.FinestLevel(); i++)
      domainf.refine(amrDataf.RefRatio()[i]);

    Box domainc = domainf;
    for (int i=0; i<BL_SPACEDIM; i++) {
      domainc.setSmall(i,domainc.smallEnd()[i]/rratio[i]);
      domainc.setBig(i,domainf.bigEnd()[i]/rratio[i]);
    }

    BoxArray baf(domainf), bac(domainc);
    
    // Fill tmpx with the data
    Vector<int> destFillComps(nComp);
    Vector<std::string> destNames(nComp);
    for (int i=0; i<nComp; ++i) 
      destFillComps[i] = i;
    MultiFab tmpc(bac,nComp,0), tmpf(baf,nComp,0);
    amrDataf.FillVar(tmpf,amrDataf.FinestLevel(),cNames,destFillComps);
    tmpc.setVal(0.);

    Real fv = 1.0/(rratio[0]*rratio[1]);

    for (MFIter mfi(tmpc); mfi.isValid(); ++mfi) {
	  
      const int* lo = tmpc[mfi].loVect();
      const int* hi = tmpc[mfi].hiVect();

      for (int joff=0; joff<rratio[1]; joff++) {
	for (int jc=lo[1]; jc<=hi[1]; jc++) {
	  int j = jc*rratio[1] +joff;
	  for (int ioff=0; ioff<rratio[0]; ioff++) {
	    for (int ic=lo[0]; ic<=hi[0]; ic++) {  
	      int i = ic*rratio[0] + ioff;
	      for (int n=0; n<nComp; n++) {
		tmpc[mfi](IntVect(ic,jc),n) += fv*
		  tmpf[mfi](IntVect(i,j),n);
	      }
	    }
	  }
	}
      }

      for (int joff=0; joff<rratio[1]; joff++) {
	for (int jc=lo[1]; jc<=hi[1]; jc++) {
	  int j = jc*rratio[1] +joff;
	  for (int ioff=0; ioff<rratio[0]; ioff++) {
	    for (int ic=lo[0]; ic<=hi[0]; ic++) {  
	      int i = ic*rratio[0] + ioff;
	      for (int n=0; n<nComp; n++) {
		tmpf[mfi](IntVect(i,j),n) -= 
		  tmpc[mfi](IntVect(ic,jc),n);
	      }
	    }
	  }
	}
      }
    }

    VisMF::Write(tmpf,oFile);
}

