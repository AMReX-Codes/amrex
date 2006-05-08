
//
// $Id: DiffSameGrid.cpp,v 1.14 2006-05-08 21:41:33 lijewski Exp $
//

#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
using std::ios;

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
    cout << "This utility performs a diff operation between two"     << endl
         << "plotfiles which have the exact same grids."             << endl;
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

static
bool
amrDatasHaveSameDerives(const AmrData& amrd1,
			const AmrData& amrd2)
{
    const Array<std::string>& derives1 = amrd1.PlotVarNames();
    const Array<std::string>& derives2 = amrd2.PlotVarNames();
    int length = derives1.size();
    if (length != derives2.size())
	return false;
    for (int i=0; i<length; ++i)
	if (derives1[i] != derives2[i])
	    return false;
    return true;
}

int
main (int   argc,
      char* argv[])
{
    if (argc == 1)
        PrintUsage(argv[0]);

    BoxLib::Initialize(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    std::string iFile1, iFile2, difFile;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }
    pp.query("infile1", iFile1);
    if (iFile1.empty())
        BoxLib::Abort("You must specify `infile1'");

    pp.query("infile2", iFile2);
    if (iFile2.empty())
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
    const Array<std::string>& derives = amrDataI.PlotVarNames();
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

        if (baI.size() != baE.size())
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

        for (MFIter mfi(*error[iLevel]); mfi.isValid(); ++mfi)
        {
            for (int iComp = 0; iComp < nComp; iComp++)
            {
                Real grdL2 = (*error[iLevel])[mfi].norm(norm, iComp, 1);

                if (norm != 0)
                {
                    norms[iComp] = norms[iComp] + pow(grdL2, norm);
                }
                else
                {
                    norms[iComp] = std::max(norms[iComp], grdL2);
                }
                
            }
        }


#ifdef BL_USE_MPI
        MPI_Datatype datatype = ParallelDescriptor::Mpi_typemap<Real>::type(),
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
                            norms[iComp] = std::max(norms[iComp], tmp[iComp]);
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


    if (!difFile.empty())
        WritePlotFile(error, amrDataI, difFile, verbose);
    
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
	delete error[iLevel];

    //
    // This calls ParallelDescriptor::EndParallel() and exit()
    //
    DataServices::Dispatch(DataServices::ExitRequest, NULL);

    BoxLib::Finalize();

    return 0;
}

