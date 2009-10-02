
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
		
		Real* tmp = 0;
		int   tmplen = 0;
		

		for (int iComp=0; iComp<nComp+sComp; ++iComp)
		{
		  n1[iComp] += vwFab.norm(fabbox, 1, iComp, 1);
		  n2[iComp] += vwFabSqrd.norm(fabbox, 1, iComp, 1);

		  delete [] tmp;

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
    for (int iComp=0; iComp<nComp+sComp; ++iComp)
    {
	mean[iComp] /= total_volume;
	variance[iComp] = variance[iComp]/total_volume - mean[iComp]*mean[iComp];
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
