
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>

#include <WritePlotFile.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <AMReX_AmrData.H>

#ifdef AMREX_DEBUG
#include <TV_TempWrite.H>
#endif

#define GARBAGE 666.e+40
using namespace amrex;

static
void
PrintUsage (const char* progName)
{
    std::cout << "\nThis utility performs a diff operation between two"     << std::endl
         << "plotfiles which have the exact same grids."             << std::endl;
    std::cout << '\n';
    std::cout << "Usage:" << '\n';
    std::cout << progName << '\n';
    std::cout << "    infile1 = inputFileName1" << '\n';
    std::cout << "    infile2 = inputFileName2" << '\n';
    std::cout << "    diffile = differenceFileName" << '\n';
    std::cout << "              (If not specified no file is written)" << '\n';
    std::cout << "       norm = integer norm (Ie. default is 2 for L2 norm)" << '\n';
    std::cout << "   [help=t|f]" << '\n';
    std::cout << "   [verbose=t|f]" << '\n';
    std::cout << '\n';
    exit(1);
}

static
bool
amrDatasHaveSameDerives(const AmrData& amrd1,
			const AmrData& amrd2)
{
    const Vector<std::string>& derives1 = amrd1.PlotVarNames();
    const Vector<std::string>& derives2 = amrd2.PlotVarNames();
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
    amrex::Initialize(argc,argv);

    if (argc == 1)
        PrintUsage(argv[0]);

    ParmParse pp;

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    std::string iFile1, iFile2, difFile;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }

    pp.query("infile1", iFile1);
    if (iFile1.empty())
        amrex::Abort("You must specify `infile1'");

    pp.query("infile2", iFile2);
    if (iFile2.empty())
        amrex::Abort("You must specify `infile2'");

    pp.query("diffile", difFile);

    int norm = 2;
    pp.query("norm", norm);

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    
    DataServices dataServicesC(iFile1, fileType);
    DataServices dataServicesF(iFile2, fileType);

    if (!dataServicesC.AmrDataOk() || !dataServicesF.AmrDataOk())
        amrex::Abort("ERROR: Dataservices not OK");
    //
    // Generate AmrData Objects 
    //
    AmrData& amrDataI = dataServicesC.AmrDataRef();
    AmrData& amrDataE = dataServicesF.AmrDataRef();
    //
    // Initial Tests 
    //
    if (!amrDatasHaveSameDerives(amrDataI,amrDataE))
        amrex::Abort("ERROR: Plotfiles do not have the same state variables");

    if (amrDataI.FinestLevel() != amrDataE.FinestLevel())
        amrex::Abort("ERROR: Finest level is not the same in the two plotfiles");

    int nComp       = amrDataI.NComp();
    int finestLevel = amrDataI.FinestLevel();
    const Vector<std::string>& derives = amrDataI.PlotVarNames();
    Vector<int> destComps(nComp);
    for (int i = 0; i < nComp; i++) 
        destComps[i] = i;
    //
    // Compute the error
    //
    Vector<MultiFab*> error(finestLevel+1);
    
    if (ParallelDescriptor::IOProcessor())
        std::cout << "L"<< norm << " norm of Error in Each Component" << std::endl
             << "-----------------------------------------" << std::endl;

    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        const BoxArray& baI = amrDataI.boxArray(iLevel);
        const BoxArray& baE = amrDataE.boxArray(iLevel);

        if (baI != baE)
        {
            std::cout << "ERROR: BoxArrays are not the same at level " << iLevel << std::endl;
            ParallelDescriptor::Abort();
        }

        DistributionMapping dmI(baI);
	error[iLevel] = new MultiFab(baI, dmI, nComp, 0);
	error[iLevel]->setVal(GARBAGE);

        DistributionMapping dmE(baE);
        MultiFab dataI(baI, dmI, nComp, 0);
        MultiFab dataE(baE, dmE, nComp, 0);

        amrDataI.FillVar(dataI, iLevel, derives, destComps);
        amrDataE.FillVar(dataE, iLevel, derives, destComps);

        for (int i = 0; i < destComps.size(); i++)
        {
            amrDataI.FlushGrids(destComps[i]);
            amrDataE.FlushGrids(destComps[i]);
        }

        (*error[iLevel]).copy(dataI);
        (*error[iLevel]).minus(dataE, 0, nComp, 0);
        //
        // Output Statistics
        //
        if (ParallelDescriptor::IOProcessor())
	  std::cout << "Level:  " << iLevel << std::endl;

        Vector<Real> norms(nComp,0);

        for (MFIter mfi(*error[iLevel]); mfi.isValid(); ++mfi)
        {
            for (int iComp = 0; iComp < nComp; iComp++)
            {
                const Real grdL2 = (*error[iLevel])[mfi].norm(norm, iComp, 1);

                if (norm != 0)
                    norms[iComp] = norms[iComp] + pow(grdL2, norm);
                else
                    norms[iComp] = std::max(norms[iComp], grdL2);
            }
        }

#ifdef BL_USE_MPI
        MPI_Datatype datatype = ParallelDescriptor::Mpi_typemap<Real>::type();

        if (ParallelDescriptor::IOProcessor())
        {
            Vector<Real> tmp(nComp);
            for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
            {
                if (proc != ParallelDescriptor::IOProcessorNumber())
                {
                    MPI_Status stat;
                    int rc = MPI_Recv(tmp.dataPtr(), nComp, datatype, 
                                      MPI_ANY_SOURCE, proc, ParallelDescriptor::Communicator(), 
                                      &stat);

                    if (rc != MPI_SUCCESS)
                        ParallelDescriptor::Abort(rc);

                    for (int iComp = 0; iComp < nComp; iComp++)
                    {
                        if (norm != 0)
                            norms[iComp] = norms[iComp] + tmp[iComp];
                        else
                            norms[iComp] = std::max(norms[iComp], tmp[iComp]);
                    }
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

                std::cout << "  " << derives[iComp] << ": " << norms[iComp] << std::endl;
            }
            std::cout << std::endl;
        }
    }

    if (!difFile.empty())
        WritePlotFile(error, amrDataI, difFile, verbose);
    
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
	delete error[iLevel];

    amrex::Finalize();

    return 0;
}

