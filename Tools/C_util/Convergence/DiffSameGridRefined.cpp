
//
// $Id: DiffSameGridRefined.cpp,v 1.6 2001-12-03 22:39:22 lijewski Exp $
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
#include "AVGDOWN_F.H"

#ifndef NDEBUG
#include "TV_TempWrite.H"
#endif

#define GARBAGE 666.e+40

static
void
PrintUsage (const char* progName)
{
    cout << "This utility performs a diff operation between two"     << endl
         << "plotfiles which have the same geometrical grid"         << endl
         << "configurations but a factor of refinement between"      << endl
         << "the grids from each plotfile at the same level."        << endl
         << "For instance, it works to diff two plotfiles having"    << endl
         << "  Plotfile 1: 25x25 base grid, Ref_Ratio = 2"           << endl
         << "  Plotfile 2: 50x50 base grid, Ref_Ratio = 2"           << endl
         << "Should also work for,"                                  << endl
         << "  Plotfile 1: 25x25 base grid, Ref_Ratio = 2"           << endl
         << "  Plotfile 2: 25x25 base grid, Ref_Ratio = 4"           << endl
         << "In both cases, the geometrical region which is refined" << endl
         << "must be the same.  So, this is generally good for"      << endl
         << "comparing cases ran using a fixed grid file."           << endl;
    cout << '\n';
    cout << "Usage:" << '\n';
    cout << progName << '\n';
    cout << "    infile1 = inputFileName1" << '\n';
    cout << "    reffile = refinedPlotFile" << '\n';
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

IntVect
getRefRatio(const Box& crse,
            const Box& fine);


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
//      AmrData::SetVerbose(true);
    }
    pp.query("infile1", iFile1);
    if (iFile1.isNull())
        BoxLib::Abort("You must specify `infile1'");

    pp.query("reffile", iFile2);
    if (iFile2.isNull())
        BoxLib::Abort("You must specify `reffile'");

    pp.query("diffile", difFile);

    int norm = 2;
    pp.query("norm", norm);

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    
    DataServices dataServices1(iFile1, fileType);
    DataServices dataServices2(iFile2, fileType);

    if (!dataServices1.AmrDataOk() || !dataServices2.AmrDataOk())
        BoxLib::Abort("ERROR: Dataservices not OK");


    //
    // Generate AmrData Objects 
    //
    AmrData& amrData1 = dataServices1.AmrDataRef();
    AmrData& amrData2 = dataServices2.AmrDataRef();

    //
    // Initial Tests 
    //
    if (!amrDatasHaveSameDerives(amrData1,amrData2))
        BoxLib::Abort("ERROR: Plotfiles do not have the same state variables");

    if (amrData1.FinestLevel() != amrData2.FinestLevel())
        BoxLib::Abort("ERROR: Finest level is not the same in the two plotfiles");

    int nComp       = amrData1.NComp();
    int finestLevel = amrData1.FinestLevel();
    const Array<aString>& derives = amrData1.PlotVarNames();
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
        //
        // Construct level box array and check that the lengths are the same
        //
        const BoxArray& ba1 = amrData1.boxArray(iLevel);
        const BoxArray& ba2 = amrData2.boxArray(iLevel);

        if (ba1.length() != ba2.length())
        {
            cout << "ERROR: BoxArray lengths are not the same at level " 
                 << iLevel << endl;
            ParallelDescriptor::Abort();
        }

        //
        // Construct MultiFab for errors
        //
	error[iLevel] = new MultiFab(ba1, nComp, 0);
	error[iLevel]->setVal(GARBAGE);

        //
        // Construct refinement ratio, build the coarsened boxarray
        // (amrData2 is the refined plotfile)
        //
        const Box& domain1     = amrData1.ProbDomain()[iLevel];
        const Box& domain2     = amrData2.ProbDomain()[iLevel];
        IntVect refine_ratio   = getRefRatio(domain1, domain2);
        if (refine_ratio == IntVect())
            BoxLib::Error("Cannot find refinement ratio from data to exact");

        if (verbose)
            cerr << "level = " << iLevel << "  Ref_Ratio = " << refine_ratio
                                                             << endl;

        BoxArray ba2Coarse(ba2);
        ba2Coarse.coarsen(refine_ratio);

        //
        // For each component, average the fine fields down and calculate
        // the errors
        //
        Array<Real> norms(nComp);
        for (int iComp = 0; iComp < nComp; iComp++)
            norms[iComp] = 0.0;

        for (int iComp = 0; iComp < nComp; ++iComp)
        {
            MultiFab& data1 = amrData1.GetGrids(iLevel, iComp);
            MultiFab& data2Fine = amrData2.GetGrids(iLevel, iComp);


            //
            // Calculate the errors  for each FAB in the MultiFab
            //
            for (MultiFabIterator mfi1(data1); mfi1.isValid(); ++mfi1)
            {
                DependentMultiFabIterator mfi2Fine(mfi1, data2Fine);
                DependentMultiFabIterator errMfi(mfi1, *error[iLevel]);
    
                //
                // Create the Coarsened version of data2
                //
                FArrayBox data2Coarse(ba2Coarse[mfi1.index()], 1);
                int ncCoarse = data2Coarse.nComp();

                FORT_CV_AVGDOWN(data2Coarse.dataPtr(),
                                  ARLIM(data2Coarse.loVect()),
                                  ARLIM(data2Coarse.hiVect()),
                                &ncCoarse,
                                mfi2Fine().dataPtr(),
                                  ARLIM(mfi2Fine().loVect()),
                                  ARLIM(mfi2Fine().hiVect()), 
                                ba2Coarse[mfi1.index()].loVect(), 
                                ba2Coarse[mfi1.index()].hiVect(),
                                refine_ratio.getVect());


                //
                // Calculate the errors on this FAB for this component
                //
                errMfi().copy(mfi1(), 0, iComp, 1);
                errMfi().minus(data2Coarse, 0, iComp, 1);

   
                Real grdL2 = errMfi().norm(norm, iComp, 1);

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


        //
        // Output Statistics
        //
        if (ParallelDescriptor::IOProcessor())
            cout << "  " << iLevel << "    ";


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

                cout << norms[iComp] << " ";
            }
            cout << endl;
        }
    }


    if (!difFile.isNull())
        WritePlotFile(error, amrData1, difFile, verbose);
    
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


IntVect
getRefRatio(const Box& crse,
            const Box& fine)
{
    // Compute refinement ratio between crse and fine boxes, return invalid
    // IntVect if there is none suitable
    ParmParse pp("");
    Array<int> rr_in(BL_SPACEDIM,-1);
    int Nrr = pp.countval("ref_ratio",Nrr);
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

