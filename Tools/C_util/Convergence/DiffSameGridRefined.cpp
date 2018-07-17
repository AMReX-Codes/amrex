
#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
using std::ios;

#include <unistd.h>

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
#include <AVGDOWN_F.H>

#ifdef AMREX_DEBUG
#include <TV_TempWrite.H>
#endif

#define GARBAGE 666.e+40
using namespace amrex;
static
void
PrintUsage (const char* progName)
{
    std::cout << "This utility performs a diff operation between two"   << std::endl
         << "plotfiles which have the same geometrical domain and grid" << std::endl
         << "configuration but a factor of refinement between"          << std::endl
         << "the cells from each plotfile at the same level."           << std::endl
         << "For instance, it works to diff two plotfiles having"       << std::endl
         << "  Plotfile 1: 25x25 base grid, Ref_Ratio = 2"              << std::endl
         << "  Plotfile 2: 50x50 base grid, Ref_Ratio = 2"              << std::endl
         << "Should also work for,"                                     << std::endl
         << "  Plotfile 1: 25x25 base grid, Ref_Ratio = 2"              << std::endl
         << "  Plotfile 2: 25x25 base grid, Ref_Ratio = 4"              << std::endl
         << "In both cases, the geometrical region which is refined"    << std::endl
         << "must be the same.  So, this is generally good for"         << std::endl
         << "comparing cases ran using a fixed grid file."              << std::endl;
    std::cout << '\n';
    std::cout << "Usage:" << '\n';
    std::cout << progName << '\n';
    std::cout << "    infile1 = inputFileName1" << '\n';
    std::cout << "    reffile = refinedPlotFile" << '\n';
    std::cout << "    diffile = differenceFileName" << '\n';
    std::cout << "              (If not specified no file is written)" << '\n';
    std::cout << "       norm = integer norm (Ie. default is 2 for L2 norm)" << '\n';
    std::cout << "   [-help]" << '\n';
    std::cout << "   [-verbose]" << '\n';
    std::cout << '\n';
    exit(1);
}

bool
amrDatasHaveSameDerives(const AmrData& amrd1,
			const AmrData& amrd2);

IntVect
getRefRatio(const Box& crse,
            const Box& fine);


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

    FArrayBox::setFormat(FABio::FAB_IEEE_32);
    //
    // Scan the arguments.
    //
    std::string iFile1, iFile2, difFile;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
    }
    pp.query("infile1", iFile1);
    if (iFile1.empty())
        amrex::Abort("You must specify `infile1'");

    pp.query("reffile", iFile2);
    if (iFile2.empty())
        amrex::Abort("You must specify `reffile'");

    pp.query("diffile", difFile);

    int norm = 2;
    pp.query("norm", norm);

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    
    DataServices dataServices1(iFile1, fileType);
    DataServices dataServices2(iFile2, fileType);

    if (!dataServices1.AmrDataOk() || !dataServices2.AmrDataOk())
        amrex::Abort("ERROR: Dataservices not OK");


    //
    // Generate AmrData Objects 
    //
    AmrData& amrData1 = dataServices1.AmrDataRef();
    AmrData& amrData2 = dataServices2.AmrDataRef();

    //
    // Initial Tests 
    //
    if (!amrDatasHaveSameDerives(amrData1,amrData2))
        amrex::Abort("ERROR: Plotfiles do not have the same state variables");

    if (amrData1.FinestLevel() != amrData2.FinestLevel())
        amrex::Abort("ERROR: Finest level is not the same in the two plotfiles");

    int nComp       = amrData1.NComp();
    int finestLevel = amrData1.FinestLevel();
    const Vector<std::string>& derives = amrData1.PlotVarNames();
    Vector<int> destComps(nComp);
    for (int i = 0; i < nComp; i++) 
        destComps[i] = i;
    

    //
    // Compute the error
    //
    Vector<MultiFab*> error(finestLevel+1);
    
    if (ParallelDescriptor::IOProcessor())
        std::cout << "Level  L"<< norm << " norm of Error in Each Component" << std::endl
             << "-----------------------------------------------" << std::endl;

    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        //
        // Construct level box array and check that the lengths are the same
        //
        const BoxArray& ba1 = amrData1.boxArray(iLevel);
        const BoxArray& ba2 = amrData2.boxArray(iLevel);

        if (ba1.size() != ba2.size())
        {
            std::cout << "ERROR: BoxArray lengths are not the same at level " 
                 << iLevel << std::endl;
            ParallelDescriptor::Abort();
        }

        //
        // Construct MultiFab for errors
        //
        DistributionMapping dm(ba1);
	error[iLevel] = new MultiFab(ba1, dm, nComp, 0);
	error[iLevel]->setVal(GARBAGE);

        //
        // Construct refinement ratio, build the coarsened boxarray
        // (amrData2 is the refined plotfile)
        //
        const Box& domain1     = amrData1.ProbDomain()[iLevel];
        const Box& domain2     = amrData2.ProbDomain()[iLevel];
        IntVect refine_ratio   = getRefRatio(domain1, domain2);
        if (refine_ratio == IntVect())
            amrex::Error("Cannot find refinement ratio from data to exact");

        if (verbose)
            std::cerr << "level = " << iLevel << "  Ref_Ratio = " << refine_ratio
                                                                  << std::endl;

        BoxArray ba2Coarse(ba2);
        ba2Coarse.coarsen(refine_ratio);

        //
        // For each component, average the fine fields down and calculate
        // the errors
        //
        Vector<Real> norms(nComp);
        for (int iComp = 0; iComp < nComp; iComp++)
            norms[iComp] = 0.0;

        for (int iComp = 0; iComp < nComp; ++iComp)
        {
            MultiFab& data1 = amrData1.GetGrids(iLevel, iComp);
            MultiFab& data2Fine = amrData2.GetGrids(iLevel, iComp);


            //
            // Calculate the errors  for each FAB in the MultiFab
            //
            for (MFIter mfi(data1); mfi.isValid(); ++mfi)
            {
                //
                // Create the Coarsened version of data2
                //
                int index = mfi.index();

		const Box& bx = ba2Coarse[index];
                FArrayBox data2Coarse(bx, 1);
                int ncCoarse = data2Coarse.nComp();

                FORT_CV_AVGDOWN(data2Coarse.dataPtr(),
                                  ARLIM(bx.loVect()), ARLIM(bx.hiVect()),
                                &ncCoarse,
                                data2Fine[mfi].dataPtr(),
                                  ARLIM(data2Fine[mfi].loVect()),
                                  ARLIM(data2Fine[mfi].hiVect()), 
                                bx.loVect(), bx.hiVect(),
                                refine_ratio.getVect());


                //
                // Calculate the errors on this FAB for this component
                //
                (*error[iLevel])[mfi].copy(data1[mfi], 0, iComp, 1);
                (*error[iLevel])[mfi].minus(data2Coarse, 0, iComp, 1);

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


        //
        // Output Statistics
        //
        if (ParallelDescriptor::IOProcessor())
            std::cout << "  " << iLevel << "    ";


#ifdef BL_USE_MPI

        MPI_Datatype datatype = ParallelDescriptor::Mpi_typemap<Real>::type();

        if (ParallelDescriptor::IOProcessor())
        {
            Vector<Real> tmp(nComp);
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
            vol *= amrData1.DxLevel()[iLevel][dir];

        if (ParallelDescriptor::IOProcessor())
        {
            for (int iComp = 0; iComp < nComp; iComp++)
            {
                if (norm != 0)
                {
                    norms[iComp] = norms[iComp] * vol;
                    norms[iComp] = pow(norms[iComp], (1.0/norm));
                }

                std::cout << norms[iComp] << " ";
            }
            std::cout << std::endl;
        }
    }


    if (!difFile.empty())
        WritePlotFile(error, amrData1, difFile, verbose);
    
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
	delete error[iLevel];

    amrex::Finalize();
}


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


IntVect
getRefRatio(const Box& crse,
            const Box& fine)
{
    // Compute refinement ratio between crse and fine boxes, return invalid
    // IntVect if there is none suitable
    ParmParse pp("");
    Vector<int> rr_in(BL_SPACEDIM,-1);
    int Nrr = pp.countval("ref_ratio");
    BL_ASSERT(Nrr==0 || Nrr==BL_SPACEDIM || Nrr==1);
    if (Nrr>0)
    {
        pp.queryarr("ref_ratio",rr_in,0,Nrr);
        return IntVect(rr_in);
    }

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

