
#include <string>

#include "ComputeAmrDataStat.H"

#include "WritePlotFile.H"
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "Utility.H"
#include "VisMF.H"

void
ComputeAmrDataMeanVar  (AmrData&     amrData,
			Array<Real>& mean,
			Array<Real>& variance,
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
	    BoxLib::Abort("You must specify `outfile' if run in verbose mode");
    }
    
    int finestLevel = amrData.FinestLevel();

    mean.clear(); mean.resize(nComp,0.0);
    variance.clear(); variance.resize(nComp,0.0);
    
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

	error[iLevel] = new MultiFab(ba, nComp, 0);
	for (int iComp=0; iComp<nComp; ++iComp)
	{
	    MultiFab& data = amrData.GetGrids(iLevel,iComp);
	    error[iLevel]->copy(data,0,iComp+sComp,1);
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
	    Array<Real> n1(nComp,0.0), n2(nComp,0.0);

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
	variance[iComp] = variance[iComp]/total_volume - mean[iComp]*mean[iComp];
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
		     Array<MultiFab*> mean,
		     Array<MultiFab*> variance,
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

// determine the pdf of a plotfile
void
ComputeAmrDataPDF (AmrData&           amrData,
		   Real**             icount,
		   int                nBin,
		   Array<std::string> cNames,
		   Array<BoxArray>    bas)
{
    std::string oFile, iFileDir, oFileDir;
    
    // if (verbose)
//     {
// 	ParmParse pp;
// 	pp.query("outfile", oFile);
// 	if (oFile.empty())
// 	    BoxLib::Abort("You must specify `outfile' if run in verbose mode");
//     }

    int finestLevel = amrData.FinestLevel();
    int nComp = cNames.size();
    for (int iComp=0; iComp<nComp; iComp++) {
      for (int iBin=0; iBin<nBin; iBin++){
	icount[iComp][iBin] = 0;
      }
    }

    // determine min and max
    Array<Real> smin(nComp,1.e20), smax(nComp,-1.e20);
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

    Array<int> refMult(finestLevel + 1, 1);
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
    Array<int> destFillComps(nComp);
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

//     if (ParallelDescriptor::IOProcessor() && verbose)
//     {
// 	std::cout << "Writing zeroed state to " << oFile << '\n';
// 	WritePlotFile(error, amrData, oFile, verbose);
//     }

}


// determine the pdf of a portion of a plotfile
void
ComputeAmrDataPDF (AmrData&     amrData,
		   Real**       icount,
		   int          nBin,
		   Array<std::string> cNames)
{
    std::string oFile, iFileDir, oFileDir;
    
    // if (verbose)
//     {
// 	ParmParse pp;
// 	pp.query("outfile", oFile);
// 	if (oFile.empty())
// 	    BoxLib::Abort("You must specify `outfile' if run in verbose mode");
//     }

    
    
    int finestLevel = amrData.FinestLevel();
    int nComp = cNames.size();

    for (int iComp=0; iComp<nComp; iComp++) {
      for (int iBin=0; iBin<nBin; iBin++){
	icount[iComp][iBin] = 0;
      }
    }

    // determine min and max
    Array<Real> smin(nComp,1.e20), smax(nComp,-1.e20);
    Array<std::string>  VarNames = amrData.PlotVarNames();
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

    Array<int> refMult(finestLevel + 1, 1);
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
    Array<int> destFillComps(nComp);
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
    for (int iComp=0; iComp<nComp; ++iComp) {
      for (int iBin=0; iBin<nBin; iBin++) {
	icount[iComp][iBin] /= total_volume;    
      }
    }

//     if (ParallelDescriptor::IOProcessor() && verbose)
//     {
// 	std::cout << "Writing zeroed state to " << oFile << '\n';
// 	WritePlotFile(error, amrData, oFile, verbose);
//     }
}

// Compute the correlation in space
void
ComputeAmrDataVAR (AmrData&           amrData,
		   int                nBin,
		   Array<std::string> cNames,
		   Array<Real>        barr,
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

    Array<int> destFillComps(nComp);
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
    Array<Real> mean(nComp,0.0), variance(nComp,0.0);
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
      if (outputFile.fail()) BoxLib::Abort("Output file cannot be opened");
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
