//BL_COPYRIGHT_NOTICE

//
// $Id: DiffSameGrid.cpp,v 1.6 1999-09-24 21:26:34 sstanley Exp $
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
        BoxLib::Abort("ERROR: Dataservices not OK");


    //
    // Generate AmrData Objects 
    //
    AmrData& amrDataI = dataServicesC.AmrDataRef();
    AmrData& amrDataE = dataServicesF.AmrDataRef();

    //
    // Initial Tests 
    //
    if (!amrDatasHaveSameDerives(amrDataI,amrDataE))
        BoxLib::Abort("ERROR: Plotfiles do not have the same state variables");

    if (amrDataI.FinestLevel() != amrDataE.FinestLevel())
        BoxLib::Abort("ERROR: Finest level is not the same in the two plotfiles");

    int nComp       = amrDataI.NComp();
    int finestLevel = amrDataI.FinestLevel();
    const Array<aString>& derives = amrDataI.PlotVarNames();
    Array<int> destComps(nComp);
    for (int i = 0; i < nComp; i++) 
        destComps[i] = i;
    

    //
    // Compute the error
    //
    Array<MultiFab*> error(finestLevel+1);
    
    if (ParallelDescriptor::IOProcessor())
        cout << "Level  L"<< norm << " norm of Error in Each Component" << endl
             << "-----------------------------------------------" << endl;

    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        const BoxArray& baI = amrDataI.boxArray(iLevel);
        const BoxArray& baE = amrDataE.boxArray(iLevel);

        if (baI.length() != baE.length())
        {
            cout << "ERROR: BoxArray lengths are not the same at level " 
                 << iLevel << endl;
            ParallelDescriptor::Abort();
        }

	error[iLevel] = new MultiFab(baI, nComp, 0);
	error[iLevel]->setVal(GARBAGE);

        MultiFab dataI(baI, nComp, 0);
        MultiFab dataE(baE, nComp, 0);

        amrDataI.FillVar(dataI, iLevel, derives, destComps);
        amrDataE.FillVar(dataE, iLevel, derives, destComps);

        (*error[iLevel]).copy(dataI);
        (*error[iLevel]).minus(dataE, 0, nComp, 0);

   
        //
        // Output Statistics
        //
        if (ParallelDescriptor::IOProcessor())
            cout << "  " << iLevel << "    ";

        Array<Real> norms(nComp);
        for (int iComp = 0; iComp < nComp; iComp++)
            norms[iComp] = 0.0;

        for (MultiFabIterator mfi(*error[iLevel]); mfi.isValid(); ++mfi)
        {
            for (int iComp = 0; iComp < nComp; iComp++)
            {
                Real grdL2 = mfi().norm(norm, iComp, 1);

                if (norm != 0)
                {
                    norms[iComp] = norms[iComp] + pow(grdL2, norm);
                }
                else
                {
                    norms[iComp] = Max(norms[iComp], grdL2);
                }
                
            }
        }


#ifdef BL_USE_MPI
        MPI_Datatype datatype = mpi_data_type(norms.dataPtr());
        if (ParallelDescriptor::IOProcessor())
        {
            Array<Real> tmp(nComp);
            for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
                if (proc != ParallelDescriptor::IOProcessorNumber())
                {
                    MPI_Status stat;
                    int rc = MPI_Recv(tmp.dataPtr(), nComp, datatype, 
                                      MPI_ANY_SOURCE, proc, MPI_COMM_WORLD, 
                                      &stat);

                    if (rc != MPI_SUCCESS)
                        ParallelDescriptor::Abort(rc);

                    for (int iComp = 0; iComp < nComp; iComp++)
                        if (norm != 0)
                        {
                            norms[iComp] = norms[iComp] + tmp[iComp];
                        }
                        else
                        {
                            norms[iComp] = Max(norms[iComp], tmp[iComp]);
                        }
                }
        }
        else
        {
            int rc = MPI_Send(norms.dataPtr(), nComp, datatype, 
                              ParallelDescriptor::IOProcessorNumber(),
                              ParallelDescriptor::MyProc(),
                              MPI_COMM_WORLD);

            if (rc != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
        }
#endif


        if (ParallelDescriptor::IOProcessor())
        {
            for (int iComp = 0; iComp < nComp; iComp++)
            {
                if (norm != 0)
                    norms[iComp] = pow(norms[iComp], (1.0/norm));

                cout << norms[iComp] << " ";
            }
            cout << endl;
        }
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

