#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <string>

#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_PlotFileUtil.H>

#include "AugmentPlotfile_F.H"

#ifndef NDEBUG
#include <TV_TempWrite.H>
#endif

#define GARBAGE 666.e+40


using namespace amrex;

static void PrintUsage (const char* progName)
{
    std::cout << "\nThis utility adds the requested components to the plotfile " << std::endl;
    std::cout << "specified as an input and writes out the original data and   " << std::endl;
    std::cout << "addtional components to the plotfile specified as an output. " << std::endl;
    std::cout << '\n';
    std::cout << "Usage:" << '\n';
    std::cout << progName << '\n';
    std::cout << "    infile  = inputFileName" << '\n';
    std::cout << "    outfile = outputFileName" << '\n';
    std::cout << "    components = ..." << '\n';
    std::cout << "    [help=t|f]" << '\n';
    std::cout << "    [verbose=t|f]" << '\n';
    std::cout << '\n';
    exit(1);
}

static Vector<int> findVelocityComponents(const AmrData& amrd)
{
    Vector<int> velComps(3);
    velComps[0] = amrd.StateNumber("x_velocity");
    velComps[1] = amrd.StateNumber("y_velocity");
    velComps[2] = amrd.StateNumber("z_velocity");
    return velComps;
}

template <class T>
static Vector<T> concatVectors(const Vector<T>& a, const Vector<T>& b)
{
    Vector<T> ab(a.size() + b.size());
    for (int i = 0; i < a.size(); i++)
        ab[i] = a[i];
    for (int i = 0; i < b.size(); i++)
        ab[i + a.size()] = b[i];
    return ab;
}

template <class T>
static Vector<T> concatVectorsDestructive(Vector<T>& a, const Vector<T>& b)
{
    std::cout << "begin dest concat" << '\n';
    int oldSize = a.size();
    a.resize(oldSize + b.size());
    for (int i = 0; i < b.size(); i++)
        a[i + oldSize] = b[i];
    std::cout << "end dest concat" << '\n';
    return a;
}


int
main (int argc, char* argv[])
{
    //// Initialize AMRex
    amrex::Initialize(argc, argv);

    //// Parse Input Arguments
    if (argc == 1)
        PrintUsage(argv[0]);

    ParmParse pp;

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    std::string infile, outfile;

    bool verbose = false;
    if (pp.contains("verbose"))
    {
        verbose = true;
        AmrData::SetVerbose(true);
    }

    pp.query("infile", infile);
    if (infile.empty())
        amrex::Abort("You must specify `infile'");

    pp.query("outfile", outfile);
    if (outfile.empty())
        amrex::Abort("You must specify `outfile'");

    //// Set Up Data Services
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServicesC(infile, fileType);
    if (!dataServicesC.AmrDataOk())
        amrex::Abort("ERROR: Dataservices not OK");

    //// Get information from old plotfile
    AmrData& amrDataIn = dataServicesC.AmrDataRef();
    Real time = amrDataIn.Time();
    int nOldComps   = amrDataIn.NComp();
    int finestLevel = amrDataIn.FinestLevel();
    auto level_steps = amrDataIn.LevelSteps();
    // amrData gives one ref ratio per level; the write utility wants three
    // dimensional IntVects at each level.
    Vector<IntVect> ref_ratio(finestLevel + 1);
    {
        Vector<int> trr = amrDataIn.RefRatio(); // trr = temp_ref_ratio
        for (int i = 0; i <= finestLevel; ++i)
            ref_ratio[i] = IntVect(AMREX_D_DECL(trr[i], trr[i], trr[i]));
            //ref_ratio[i] = AMREX_D_DECL(trr[i], trr[i], trr[i]);
    }
    const Vector<std::string>& oldComps = amrDataIn.PlotVarNames();
    Vector<Geometry> geoms(finestLevel + 1);
    {
        Vector<Box> domains = amrDataIn.ProbDomain();
        Vector<Real> probSizes = amrDataIn.ProbSize();
        Vector<Real> ProbLo = amrDataIn.ProbLo();
        Vector<Real> ProbHi = amrDataIn.ProbHi();
        auto rbox = RealBox(
                {AMREX_D_DECL(ProbLo[0], ProbLo[1], ProbLo[2])},
                {AMREX_D_DECL(ProbHi[0], ProbHi[1], ProbHi[2])}
                );
        std::array<int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(1, 1, 1)};
        // TODO fix above; we may not not actually want periodicity
        for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
            geoms[iLevel] = Geometry(
                    domains[iLevel],
                    &rbox,
                    0,  // cartesian coords
                    is_periodic.data()
                    );
    }







    int nGhost = 1;


    // Construct information for new plotfile
    Vector<std::string> divuCompNames = {"testdivu"};
    Vector<std::string> vortCompNames = {"testvortx", "testvorty", "testvortz"};
    Vector<std::string> newCompNames = concatVectors(
            divuCompNames,
            vortCompNames
            );



    int nComps = nOldComps + newCompNames.size();
    const Vector<std::string>& allCompNames = concatVectors(
            amrDataIn.PlotVarNames(),
            newCompNames
            );
    auto velComps = findVelocityComponents(amrDataIn);
    Vector<int> divuComps = {nOldComps};
    Vector<int> vortComps = {nOldComps + 1, nOldComps + 2, nOldComps + 3};

    // Utility arrays with numbers of components
    Vector<int> oldCompNums(nOldComps);
    for (int i = 0; i < nOldComps; i++) 
        oldCompNums[i] = i;
    Vector<int> allCompNums(nComps);
    for (int i = 0; i < nComps; i++) 
        allCompNums[i] = i;





    // Make list of output multifabs
    Vector<MultiFab*> dataOut(finestLevel + 1);



    //// Print IO Info
    //if (ParallelDescriptor::IOProcessor())
    //    std::cout << "L"<< norm << " norm of Absolute and Relative Error in Each Component" << std::endl
    //        << std::setfill('-') << std::setw(80) << "-" << std::setfill(' ') << std::endl;




    //// Iterate through levels
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        // Copy structure from read file
        const BoxArray& ba = amrDataIn.boxArray(iLevel);
        DistributionMapping dm {ba};

        // Get cell size from read file
        Vector<Real> delta = amrDataIn.CellSize(iLevel);

        // Copy data from read file with interpolated ghosts
        MultiFab dataIn(ba, dm, nComps, nGhost);
        amrDataIn.FillVar(dataIn, iLevel, oldComps, oldCompNums);

        // Allocate space for new arrays
        dataOut[iLevel] = new MultiFab(ba, dm, nComps, 1);
        dataOut[iLevel]->setVal(GARBAGE); // TODO do we need this?








        // TODO
        // this part needs to be checked -- it appears that
        // FlushGrids works on all levels at once -- but we are 
        // looping over levels
        for (int i = 0; i < oldCompNums.size(); i++)
        {
            amrDataIn.FlushGrids(oldCompNums[i]);
        }







//        for (auto comp = oldCompNums.begin(); comp != oldCompNums.end(); comp++)
//            MultiFab::Copy(*dataOut[iLevel],
//                           dataIn,
//                           nOldComps,
//                           nComps,
//                           *comp,
//                           nGhost);








        // Copy Old Data
        {
            MFIter omfi(dataIn);                // old
            MFIter nmfi(*dataOut[iLevel]);      // new
            while (omfi.isValid() && nmfi.isValid())
            {
                const Box& nbox = nmfi.validbox();
                FArrayBox& nfab = (*dataOut[iLevel])[nmfi];
                const Box& nabox = nfab.box();

                const Box& obox = omfi.validbox();
                FArrayBox& ofab = dataIn[omfi];
                const Box& oabox = ofab.box();

                for (auto comp = oldCompNums.begin(); comp != oldCompNums.end(); comp++)
                    FORT_COPY_3D(
                            obox.loVect(), obox.hiVect(),
                            ofab.dataPtr(), &nOldComps,
                            oabox.loVect(), oabox.hiVect(),
                            nbox.loVect(), nbox.hiVect(),
                            nfab.dataPtr(), &nComps,
                            nabox.loVect(), nabox.hiVect(),
                            &(*comp), &(*comp)
                            );

                ++omfi;
                ++nmfi;
            }
        }
        auto period = Periodicity(amrDataIn.ProbDomain()[0].size());
        dataOut[iLevel]->FillBoundary(period);

        // Generate new components
        for (MFIter mfi(*dataOut[iLevel]); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox();
            FArrayBox& fab = (*dataOut[iLevel])[mfi];
            const Box& abox = fab.box();

            FORT_DIVU_3D (
                    box.loVect(), box.hiVect(),
                    fab.dataPtr(), &nComps,
                    abox.loVect(), abox.hiVect(),
                    &velComps[0],
                    &divuComps[0],
                    &delta[0]);

            FORT_VORT_3D (
                    box.loVect(), box.hiVect(),
                    fab.dataPtr(), &nComps,
                    abox.loVect(), abox.hiVect(),
                    &velComps[0],
                    &vortComps[0],
                    &delta[0]);
        }






        // Output Progress
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Level:  " << iLevel << std::endl;



    }





    // Write out a plotfile
    WriteMultiLevelPlotfile(
            outfile,
            finestLevel + 1,
            GetVecOfConstPtrs(dataOut),
            allCompNames,
            geoms,
            time,
            level_steps,
            ref_ratio
            );









    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel) {
        delete dataOut[iLevel];
    }

    amrex::Finalize();

    return 0;
}



