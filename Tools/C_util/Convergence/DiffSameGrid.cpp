//BL_COPYRIGHT_NOTICE

//
// $Id: DiffSameGrid.cpp,v 1.3 1999-05-10 17:18:54 car Exp $
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
    cout << "    infile1 = inputFileName1" << '\n';
    cout << "    infile2 = inputFileName2" << '\n';
    cout << "    diffile = differenceFileName" << '\n';
    cout << "              (If not specified no file is written)" << '\n';
    cout << "       norm = integer norm (Ie. default is 2 for L2 norm)" << '\n';
    cout << "   [-help]" << '\n';
    cout << "   [-verbose]" << '\n';
    cout << '\n';
    exit(1);
}

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
    aString iFile1, iFile2, difFile;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }
    pp.query("infile1", iFile1);
    if (iFile1.isNull())
        BoxLib::Abort("You must specify `infile1'");

    pp.query("infile2", iFile2);
    if (iFile2.isNull())
        BoxLib::Abort("You must specify `infile2'");

    pp.query("diffile", difFile);

    int norm = 2;
    pp.query("norm", norm);

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    
    DataServices dataServicesC(iFile1, fileType);
    DataServices dataServicesF(iFile2, fileType);

    if (!dataServicesC.AmrDataOk() || !dataServicesF.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);


    //
    // Generate AmrData Objects 
    //
    AmrData& amrDataI = dataServicesC.AmrDataRef();
    AmrData& amrDataE = dataServicesF.AmrDataRef();

    //
    // Initial Tests 
    //
    BLassert(amrDatasHaveSameDerives(amrDataI,amrDataE));
    BLassert(amrDataI.FinestLevel() == amrDataE.FinestLevel());

    int nComp          = amrDataI.NComp();
    int finestLevel = amrDataI.FinestLevel();
    const Array<aString>& derives = amrDataI.PlotVarNames();
    

    //
    // Compute the error
    //
    Array<MultiFab*> error(finestLevel+1);
    
    cout << "Level  L"<< norm << " norm of Error in Each Component" << endl
         << "-----------------------------------------------" << endl;

    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        const BoxArray& baI = amrDataI.boxArray(iLevel);
        const BoxArray& baE = amrDataE.boxArray(iLevel);

        BLassert(baI.length() == baE.length());

	error[iLevel] = new MultiFab(baI, nComp, 0);
	error[iLevel]->setVal(GARBAGE);

        MultiFab dataI(baI, nComp, 0);
        MultiFab dataE(baE, nComp, 0);

	for (int iGrid=0; iGrid<baI.length(); ++iGrid)
	{
            BLassert (baI[iGrid] == baE[iGrid]);

	    for (int iComp=0; iComp<nComp; ++iComp)
	    {
		const Box& dataGrid = baI[iGrid];
                FArrayBox tmpFab(dataGrid,1);

                amrDataI.FillVar(&tmpFab, dataGrid,
                                 iLevel, derives[iComp], 0);
                dataI[iGrid].copy(tmpFab,0,iComp,1);

                amrDataE.FillVar(&tmpFab, dataGrid,
                                 iLevel, derives[iComp], 0);
                dataE[iGrid].copy(tmpFab,0,iComp,1);
	    }

	}

        (*error[iLevel]).copy(dataI);
        (*error[iLevel]).minus(dataE, 0, nComp, 0);

   
        //
        // Output Statistics
        //
        cout << "  " << iLevel << "    ";
        for (int iComp=0; iComp<nComp; ++iComp)
        {
            Real L2 = 0.0;
            for (int iGrid=0; iGrid<baI.length(); ++iGrid)
            {
                Real grdL2 = (*error[iLevel])[iGrid].norm(norm,iComp,1);
                L2 = L2 + pow(grdL2, norm);
            }
            L2 = pow(L2, (1.0/norm));

            cout << L2 << "  ";
        }
        cout << endl;
    }

    if (!difFile.isNull())
        WritePlotFile(error, amrDataI, difFile, verbose);
    
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
	delete error[iLevel];

    //
    // This calls ParallelDescriptor::EndParallel() and exit()
    //
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
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
    
