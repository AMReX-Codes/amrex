#include "ComputeAmrDataNorms.H"

#include "WritePlotFile.H"
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "Boolean.H"
#include "ParallelDescriptor.H"
#include "Utility.H"
#include "VisMF.H"

void
ComputeAmrDataNorms (AmrData&     amrData,
		     Array<Real>& norm0,
		     Array<Real>& norm1,
		     Array<Real>& norm2,
		     bool         verbose)
{
    aString oFile, iFileDir, oFileDir;
    
    if (verbose)
    {
	ParmParse pp;
	pp.query("outfile", oFile);
	if (oFile.isNull())
	    BoxLib::Abort("You must specify `outfile' if run in verbose mode");
    }
    
    int finestLevel = amrData.FinestLevel();
    int nComp       = amrData.NComp();

    norm0.clear(); norm0.resize(nComp,0.0);
    norm1.clear(); norm1.resize(nComp,0.0);
    norm2.clear(); norm2.resize(nComp,0.0);
    
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
    // Compute the norms
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
	    error[iLevel]->copy(data,0,iComp,1);
	}

	// Zero out the error covered by fine grid
	long covered_volume = 0;
	if (iLevel != finestLevel)
	{
	    int ref_ratio = amrData.RefRatio()[iLevel];	    
	    BoxArray baF = ::BoxArray(amrData.boxArray(iLevel+1)).coarsen(ref_ratio);
	    for (MultiFabIterator mfi(*error[iLevel]); mfi.isValid(); ++mfi)
	    {
		for (int iGrid=0; iGrid<baF.length(); ++iGrid)
		{
		    Box ovlp = baF[iGrid] & mfi.validbox();
		    if (ovlp.ok())
		    {
			mfi().setVal(0.0, ovlp, 0, nComp);
			covered_volume += ovlp.numPts()*refMult[iLevel];
		    }
		}
	    }
	    ParallelDescriptor::ReduceLongSum(covered_volume);
	}

	// Compute volume at this level
	Real level_volume = 0.0;
	for (int iGrid=0; iGrid<ba.length(); ++iGrid)
	    level_volume += ba[iGrid].numPts()*refMult[iLevel];
	level_volume -= covered_volume;

	if (level_volume > 0.0)
	{
	    // Convert volume in numPts to volume in number of fine cells
	    total_volume += level_volume;
	    
	    // Get norms at this level
	    Array<Real> n0(nComp,0.0), n1(nComp,0.0), n2(nComp,0.0);
	    for (MultiFabIterator mfi(*error[iLevel]); mfi.isValid(); ++mfi)
	    {
		FArrayBox& fab = mfi();
		const Box& fabbox = mfi.validbox();
		FArrayBox vwFab(fabbox,nComp);
		FArrayBox vwFabSqrd(fabbox,nComp);
		
		// compute volume-weighted error
		vwFab.copy(fab,0,0,nComp);
		vwFab.mult(refMult[iLevel]);

		vwFabSqrd.copy(fab,0,0,nComp);
		vwFabSqrd.mult(fab,0,0,nComp);
		vwFabSqrd.mult(refMult[iLevel]);
		
		for (int iComp=0; iComp<nComp; ++iComp)
		{
		    n0[iComp] = Max(n0[iComp], fab.norm(fabbox, 0, iComp, 1));
		    n1[iComp] += vwFab.norm(fabbox, 1, iComp, 1);
		    n2[iComp] += vwFabSqrd.norm(fabbox, 1, iComp, 1);
		}
	    }

	    // Do necessary communication, then blend this level's norms
	    //  in with the running global values
	    for (int iComp=0; iComp<nComp; ++iComp)
	    {
		ParallelDescriptor::ReduceRealMax(n0[iComp]);
		ParallelDescriptor::ReduceRealSum(n1[iComp]);
		ParallelDescriptor::ReduceRealSum(n2[iComp]);
		
		norm0[iComp] = Max(norm0[iComp],n0[iComp]);
		norm1[iComp] += n1[iComp];
		norm2[iComp] += n2[iComp];
	    }
	}
    }
    for (int iComp=0; iComp<nComp; ++iComp)
    {
	norm1[iComp] /= total_volume;
	norm2[iComp] = sqrt(norm2[iComp] / total_volume);
    }
    
    if (ParallelDescriptor::IOProcessor() && verbose)
    {
	cout << "Writing zeroed state to " << oFile << '\n';
	WritePlotFile(error, amrData, oFile, verbose);
    }

    // Clean up memory
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
	delete error[iLevel];
}
