
//
// $Id: PltFileNormB.cpp,v 1.4 2001-10-17 17:53:33 lijewski Exp $
//

#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
using std::ios;
using std::set_new_handler;

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

#define GARBAGE 666.e+40

static
void
PrintUsage (const char* progName)
{
    cout << '\n';
    cout << "Usage:" << '\n';
    cout << progName << '\n';
    cout << "    infile = inputFileName" << '\n';
    cout << "      norm = integer norm (Ie. default is 2 for L2 norm)" << '\n';
    cout << "   [-help]" << '\n';
    cout << "   [-verbose]" << '\n';
    cout << '\n';
    exit(1);
}

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
    aString iFile;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }
    pp.query("infile", iFile);
    if (iFile.isNull())
        BoxLib::Abort("You must specify `infile'");

    int norm = 2;
    pp.query("norm", norm);

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    
    DataServices dataServicesC(iFile, fileType);

    if (!dataServicesC.AmrDataOk())
        BoxLib::Abort("ERROR: Dataservices not OK");


    //
    // Generate AmrData Objects 
    //
    AmrData& amrDataI = dataServicesC.AmrDataRef();

    //
    // Initial Tests 
    //
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

        MultiFab dataI(baI, nComp, 0);

        amrDataI.FillVar(dataI, iLevel, derives, destComps);

        //
        // Output Statistics
        //
        if (ParallelDescriptor::IOProcessor())
            cout << "  " << iLevel << "    ";

        Array<Real> norms(nComp);
        for (int iComp = 0; iComp < nComp; iComp++)
            norms[iComp] = 0.0;

        for (MultiFabIterator mfi(dataI); mfi.isValid(); ++mfi)
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
                                      MPI_ANY_SOURCE, proc, ParallelDescriptor::Communicator(), 
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
                              ParallelDescriptor::Communicator());

            if (rc != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
        }
#endif


        Real vol = 1.0;
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
            vol *= amrDataI.DxLevel()[iLevel][dir];

        if (ParallelDescriptor::IOProcessor())
        {
            for (int iComp = 0; iComp < nComp; iComp++)
            {
                if (norm != 0)
                {
                    norms[iComp] = norms[iComp] * vol;
                    norms[iComp] = pow(norms[iComp], (1.0/norm));
                }

                cout << norms[iComp] << " ";
            }
            cout << endl;
        }
    }


    //
    // This calls ParallelDescriptor::EndParallel() and exit()
    //
    DataServices::Dispatch(DataServices::ExitRequest, NULL);
}

