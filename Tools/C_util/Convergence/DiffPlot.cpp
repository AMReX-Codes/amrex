//BL_COPYRIGHT_NOTICE

//
// $Id: DiffPlot.cpp,v 1.4 1999-05-10 18:54:27 car Exp $
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
    cout << "    infile  = inputFileName" << '\n';
    cout << "    exact   = exactFileName" << '\n';
    cout << "    outfile = outputFileName" << '\n';
    cout << "       norm = integer norm (Ie. default is 2 for Ln norm)" << '\n';
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

    int norm = 2;
    pp.query("norm", norm);

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    
    DataServices dataServicesC(iFile, fileType);
    DataServices dataServicesF(eFile, fileType);

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
    BL_ASSERT(amrDatasHaveSameDerives(amrDataI,amrDataE));
    if (amrDataI.ProbDomain()[0] != amrDataE.ProbDomain()[0])
    {
        cerr << "ERROR: ProbDomain[0] not the same:\n"
             << "                  amrDataI.ProbDomain()[0] = " 
             <<                    amrDataI.ProbDomain()[0] << endl
             << "                  amrDataE.ProbDomain()[0] = " 
             <<                    amrDataE.ProbDomain()[0] << endl;
        BL_ASSERT(amrDataI.ProbDomain()[0] == amrDataE.ProbDomain()[0]);
        
    }

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
        Array<Real> delI = amrDataI.DxLevel()[iLevel];

	error[iLevel] = new MultiFab(baI, nComp, 0);
	error[iLevel]->setVal(GARBAGE);

        MultiFab dataI(baI, nComp, 0);
        MultiFab dataE(baI, nComp, 0);

	for (int iGrid=0; iGrid<baI.length(); ++iGrid)
	{
	    for (int iComp=0; iComp<nComp; ++iComp)
	    {
		const Box& dataGrid = baI[iGrid];
                FArrayBox tmpFab(dataGrid,1);

                amrDataI.FillVar(&tmpFab, dataGrid,
                                 finestLevel, derives[iComp], 0);
                dataI[iGrid].copy(tmpFab,0,iComp,1);

                amrDataE.FillVar(&tmpFab, dataGrid,
                                 finestLevel, derives[iComp], 0);
                dataE[iGrid].copy(tmpFab,0,iComp,1);
	    }

	}

        (*error[iLevel]).copy(dataI);
        (*error[iLevel]).minus(dataE, 0, nComp, 0);

   
        //
        // Output Statistics
        //
        Real cellvol = 1.0;
        for (int i=0; i<BL_SPACEDIM; i++)
           cellvol = cellvol * delI[i];

        Real delAvg = pow(cellvol, (1.0/BL_SPACEDIM));

        cout << "  " << iLevel << " " << delAvg << "    ";
        for (int iComp=0; iComp<nComp; ++iComp)
        {
            Real Ln = 0.0;
            for (int iGrid=0; iGrid<baI.length(); ++iGrid)
            {
                Real grdLn = (*error[iLevel])[iGrid].norm(norm,iComp,1);
                Ln = Ln + pow(grdLn, norm) * cellvol;
            }
            Ln = pow(Ln, (1.0/norm));

            cout << Ln << "  ";
        }
        cout << endl;
    }

    WritePlotFile(error, amrDataI, oFile, verbose);
    
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
    
