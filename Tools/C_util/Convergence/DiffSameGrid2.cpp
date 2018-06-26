
#include <iostream>
#include <iomanip>
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
#include <AMReX_DistributionMapping.H>

#ifdef AMREX_DEBUG
#include <TV_TempWrite.H>
#endif

#define GARBAGE 666.e+40


using namespace amrex;
static
void
PrintUsage (const char* progName)
{
    std::cout << "\nThis utility performs a diff operation between two"      << std::endl
	      << "plotfiles which have the exact same grids."                << std::endl
	      << "Both absolute and relative differences are reported."      << std::endl
	      << "Note: no attempt has been made to prevent a divide by 0, " << std::endl
	      << "so some relative errors may not be well defined."          << std::endl;
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
    // Compute the absolute and relative errors
    //
    Vector<MultiFab*> aerror(finestLevel+1);
    Vector<MultiFab*> rerror(finestLevel+1);

    // keep track of the absolute errors -- if any of them are different than 0, then
    // then we failed the comparison
    bool failed = false;

    if (ParallelDescriptor::IOProcessor())
        std::cout << "L"<< norm << " norm of Absolute and Relative Error in Each Component" << std::endl
		  << std::setfill('-') << std::setw(80) << "-" << std::setfill(' ') << std::endl;

    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        const BoxArray& baI = amrDataI.boxArray(iLevel);
        const BoxArray& baE = amrDataE.boxArray(iLevel);

        DistributionMapping dm {baI}; 
        if (baI != baE)
        {
            std::cout << "ERROR: BoxArrays are not the same at level " << iLevel << std::endl;
            ParallelDescriptor::Abort();
        }

	aerror[iLevel] = new MultiFab(baI, dm, nComp, 0);
	aerror[iLevel]->setVal(GARBAGE);

	rerror[iLevel] = new MultiFab(baI, dm, nComp, 0);
	rerror[iLevel]->setVal(GARBAGE);

        MultiFab dataI(baI, dm, nComp, 0);
        MultiFab dataE(baE, dm, nComp, 0);

        amrDataI.FillVar(dataI, iLevel, derives, destComps);
        amrDataE.FillVar(dataE, iLevel, derives, destComps);

	// this part needs to be checked -- it appears that
	// FlushGrids works on all levels at once -- but we are 
	// looping over levels
        for (int i = 0; i < destComps.size(); i++)
        {
            amrDataI.FlushGrids(destComps[i]);
            amrDataE.FlushGrids(destComps[i]);
        }

	// aerror will contain the absolute errors
        (*aerror[iLevel]).copy(dataI);
        (*aerror[iLevel]).minus(dataE, 0, nComp, 0);

	// rerror will contain the relative errors
        (*rerror[iLevel]).copy(dataI);
        (*rerror[iLevel]).minus(dataE, 0, nComp, 0);
	(*rerror[iLevel]).divide(dataI, 0, nComp, 0);


        //
        // Output Statistics
        //
        if (ParallelDescriptor::IOProcessor())
	  std::cout << "Level:  " << iLevel << std::endl;

        Vector<Real> anorms(nComp,0);
        Vector<Real> rnorms(nComp,0);

        for (MFIter mfi(*aerror[iLevel]); mfi.isValid(); ++mfi)
        {
            for (int iComp = 0; iComp < nComp; iComp++)
            {
	      
	      // compute the norm of the absolute and relative errors for the
	      // current FAB
	      const Real agrdL2 = (*aerror[iLevel])[mfi].norm(norm, iComp, 1);
	      const Real rgrdL2 = (*rerror[iLevel])[mfi].norm(norm, iComp, 1);

	      // Adding the current norm to that from the previous FABs.
	      // note: right now, we are storing anorm**norm, so that the 
	      // summation makes sense.
	      if (norm != 0) {
		anorms[iComp] = anorms[iComp] + pow(agrdL2, norm);
		rnorms[iComp] = rnorms[iComp] + pow(rgrdL2, norm);
	      }
	      else {
		anorms[iComp] = std::max(anorms[iComp], agrdL2);
		rnorms[iComp] = std::max(rnorms[iComp], rgrdL2);
	      }

            }
        }

#ifdef BL_USE_MPI
        MPI_Datatype datatype = ParallelDescriptor::Mpi_typemap<Real>::type();

        if (ParallelDescriptor::IOProcessor())
        {
            Vector<Real> atmp(nComp);
            Vector<Real> rtmp(nComp);
            for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
            {
                if (proc != ParallelDescriptor::IOProcessorNumber())
                {
                    MPI_Status stat;
                    int rc = MPI_Recv(atmp.dataPtr(), nComp, datatype, 
                                      MPI_ANY_SOURCE, 
				      proc, 
				      ParallelDescriptor::Communicator(), 
                                      &stat);

                    if (rc != MPI_SUCCESS)
                        ParallelDescriptor::Abort(rc);

                    rc = MPI_Recv(rtmp.dataPtr(), nComp, datatype, 
                                      MPI_ANY_SOURCE, 
				      ParallelDescriptor::NProcs() + proc, 
				      ParallelDescriptor::Communicator(), 
                                      &stat);

                    if (rc != MPI_SUCCESS)
                        ParallelDescriptor::Abort(rc);


                    for (int iComp = 0; iComp < nComp; iComp++)
                    {
		      if (norm != 0) {
			anorms[iComp] = anorms[iComp] + atmp[iComp];
			rnorms[iComp] = rnorms[iComp] + rtmp[iComp];
		      }
		      else {
			anorms[iComp] = std::max(anorms[iComp], atmp[iComp]);
			rnorms[iComp] = std::max(rnorms[iComp], rtmp[iComp]);
		      }

                    }
                }
            }
        }
        else
        {
            int rc = MPI_Send(anorms.dataPtr(), nComp, datatype, 
                              ParallelDescriptor::IOProcessorNumber(),
                              ParallelDescriptor::MyProc(),
                              ParallelDescriptor::Communicator());

            if (rc != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);


            rc = MPI_Send(rnorms.dataPtr(), nComp, datatype, 
                              ParallelDescriptor::IOProcessorNumber(),
                              ParallelDescriptor::NProcs() + ParallelDescriptor::MyProc(),
                              ParallelDescriptor::Communicator());

            if (rc != MPI_SUCCESS)
                ParallelDescriptor::Abort(rc);
        }
#endif

	// normalize the norms and print out
        Real vol = 1.0;
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
            vol *= amrDataI.DxLevel()[iLevel][dir];

        if (ParallelDescriptor::IOProcessor())
        {
            for (int iComp = 0; iComp < nComp; iComp++)
            {
                if (norm != 0)
                {
		  anorms[iComp] = anorms[iComp] * vol;
		  anorms[iComp] = pow(anorms[iComp], (1.0/norm));

		  rnorms[iComp] = rnorms[iComp] * vol;
		  rnorms[iComp] = pow(rnorms[iComp], (1.0/norm));
                }

                if (anorms[iComp] != 0.0) failed = true;

                std::cout << "  " << std::setw(32) << derives[iComp] 
			  << ": " << std::setw(20) << anorms[iComp] 
			  << std::setw(20) << rnorms[iComp] << std::endl;
            }
            std::cout << std::endl;
        }
    }

    if (!failed) {
      std::cout << "PLOTFILES AGREE" << std::endl;
    }

    // optionally dump out a plotfile containing the absolute errors
    if (!difFile.empty())
        WritePlotFile(aerror, amrDataI, difFile, verbose);
    
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel) {
      delete aerror[iLevel];
      delete rerror[iLevel];
    }

    amrex::Finalize();

    return 0;
}

