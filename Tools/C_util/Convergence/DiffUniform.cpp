//BL_COPYRIGHT_NOTICE

//
// $Id: DiffUniform.cpp,v 1.4 1999-05-10 18:54:27 car Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
using std::ios;
using std::set_new_handler;
#else
#include <new.h>
#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#endif

#include <unistd.h>

#include "WritePlotFile.H"
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "VisMF.H"

#ifndef NDEBUG
#include "TV_TempWrite.H"
#endif

#include "AVGDOWN_F.H"

//
// This MUST be defined if don't have pubsetbuf() in I/O Streams Library.
//
#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

#define GARBAGE 666.e+40

static
void
PrintUsage (const char* progName)
{
    cout << '\n';
    cout << "Usage:" << '\n';
    cout << progName << '\n';
    cout << "    infile  = inputFileName" << '\n';
    cout << "    exact   = exactFileName" << '\n';
    cout << "    outfile = outputFileName" << '\n';
    cout << "   [-help]" << '\n';
    cout << "   [-verbose]" << '\n';
    cout << '\n';
    exit(1);
}

int
finestLevelCoveringDomain(const AmrData& amrData);

IntVect
getRefRatio(const Box& crse,
	    const Box& fine);

bool
amrDatasHaveSameDerives(const AmrData& amrd1,
			const AmrData& amrd2);
int
main (int   argc,
      char* argv[])
{
    if (argc == 1)
        PrintUsage(argv[0]);
    //
    // Make sure to catch new failures.
    //
    set_new_handler(Utility::OutOfMemory);

    ParallelDescriptor::StartParallel(&argc, &argv);

    ParmParse pp(argc-1,argv+1);

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    aString iFileDir, iFile, eFile, oFile, oFileDir;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }
    pp.query("infile", iFile);
    if (iFile.isNull())
        BoxLib::Abort("You must specify `infile'");

    pp.query("exact", eFile);
    if (eFile.isNull())
        BoxLib::Abort("You must specify `exact' file");

    pp.query("outfile", oFile);
    if (oFile.isNull())
        BoxLib::Abort("You must specify `outfile'");

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    
    DataServices dataServicesC(iFile, fileType);
    DataServices dataServicesF(eFile, fileType);

    if (!dataServicesC.AmrDataOk() || !dataServicesF.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);

    AmrData& amrDataC = dataServicesC.AmrDataRef();
    AmrData& amrDataF = dataServicesF.AmrDataRef();

    BL_ASSERT(amrDatasHaveSameDerives(amrDataC,amrDataF));
    int exact_level = finestLevelCoveringDomain(amrDataF);
    if (exact_level < 0)
    {
	cout << "Exact data does not contain a level covering the domain" << '\n';
	BoxLib::Abort();
    }
    if (verbose)
	cout << "Using level = " << exact_level << " in 'exact' file" << '\n';
	
    int finestLevel = amrDataC.FinestLevel();
    
    //
    // Compute the error
    //
    Array<MultiFab*> error(finestLevel+1);
    
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        const BoxArray& ba = amrDataC.boxArray(iLevel);
	int nComp          = amrDataC.NComp();
	const Box& domainC = amrDataC.ProbDomain()[iLevel];
	const Box& domainF = amrDataF.ProbDomain()[exact_level];
	IntVect refine_ratio = getRefRatio(domainC, domainF);
	if (refine_ratio == IntVect())
	    BoxLib::Error("Cannot find refinement ratio from data to exact");

	error[iLevel] = new MultiFab(ba, nComp, 0);
	error[iLevel]->setVal(GARBAGE);

	for (int iComp=0; iComp<nComp; ++iComp)
	{
	    for (int iGrid=0; iGrid<ba.length(); ++iGrid)
	    {
		const Box& dataGrid = ba[iGrid];
		const Box fineGrid = ::Box(dataGrid).refine(refine_ratio);
		
		MultiFab& exact = amrDataF.GetGrids(finestLevel,iComp,fineGrid);
		int nc = exact.nComp();
		BoxArray cba = ::BoxArray(exact.boxArray()).coarsen(refine_ratio);
		MultiFab exactAvg(cba,nc,0);
		for (MultiFabIterator mfi(exact); mfi.isValid(); ++mfi)
		{
		    DependentMultiFabIterator amfi(mfi, exactAvg);
		    const FArrayBox& srcFab = mfi();
		    FArrayBox& dstFab = amfi();
		    const Box& box = amfi.validbox();
		    FORT_CV_AVGDOWN(dstFab.dataPtr(),
				    ARLIM(dstFab.loVect()), ARLIM(dstFab.hiVect()), &nc,
				    srcFab.dataPtr(),
				    ARLIM(srcFab.loVect()), ARLIM(srcFab.hiVect()),
				    box.loVect(), box.hiVect(), refine_ratio.getVect());
		}

		// Copy result of coarsening into error as temporary storage
		error[iLevel]->copy(exactAvg,0,iComp,exact.nComp());
	    }

	    // Subtract coarse data from coarsened exact data
	    MultiFab& data = amrDataC.GetGrids(iLevel,iComp);
	    for (MultiFabIterator dmfi(data); dmfi.isValid(); ++dmfi)
	    {
		DependentMultiFabIterator emfi(dmfi, *error[iLevel]);
		emfi().minus(dmfi(),0,iComp,1);
	    }
	}
    }

    WritePlotFile(error, amrDataC, oFile, verbose);
    
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
	delete error[iLevel];

    //
    // This calls ParallelDescriptor::EndParallel() and exit()
    //
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
}


int
finestLevelCoveringDomain(const AmrData& amr_data)
{
    // Find the finest level covering the entire domain.  Return
    // -1 if there isn't one suitable
    int finest_level = amr_data.FinestLevel();
    const Array<Box>& domain_array = amr_data.ProbDomain();

    for (int iLevel=finest_level; iLevel>=0; --iLevel)
    {
	const BoxArray& ba = amr_data.boxArray(iLevel);
	BoxDomain bd;
	bd.add(BoxList(ba));
	BoxDomain complement = ::complementIn(domain_array[iLevel],bd);
	if (complement.isEmpty())
	    return iLevel;
    }
    return -1;
}

IntVect
getRefRatio(const Box& crse,
	    const Box& fine)
{
    // Compute refinement ratio between crse and fine boxes, return invalid
    // IntVect if there is none suitable
    IntVect ref_ratio;
    for (int i=0; i<BL_SPACEDIM; ++i)
	ref_ratio[i] = fine.length(i) / crse.length(i);

    // Check results
    Box test1 = ::Box(fine).coarsen(ref_ratio);
    Box test2 = ::Box(test1).refine(ref_ratio);
    if (test1 != crse  ||  test2 != fine)
	ref_ratio = IntVect();
    return ref_ratio;
}

bool
amrDatasHaveSameDerives(const AmrData& amrd1,
			const AmrData& amrd2)
{
    const Array<aString>& derives1 = amrd1.PlotVarNames();
    const Array<aString>& derives2 = amrd2.PlotVarNames();
    int length = derives1.length();
    if (length != derives2.length())
	return false;
    for (int i=0; i<length; ++i)
	if (derives1[i] != derives2[i])
	    return false;
    return true;
}
    
