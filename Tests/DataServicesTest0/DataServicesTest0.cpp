// --------------------------------------------------------------------
// DataServicesTest0.cpp
// --------------------------------------------------------------------
//   this file does the following:
//     .................
//     tests DataServices and AmrData.
//     .................
// --------------------------------------------------------------------
#include <iostream>
#include <cstdlib>
#include <unistd.h>

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_AmrData.H>
#include <AMReX_DataServices.H>
#include <AMReX_Utility.H>

using namespace amrex;

const unsigned int msps(1000000);

const int XDIR(0);
const int YDIR(1);
const int ZDIR(2);

// --------------------------------------------------------------------
static void PrintUsage(const char *progName) {
    std::cout << '\n';
    std::cout << "Usage:" << '\n';
    std::cout << progName << "    plotFileName" << '\n';
    std::cout << '\n';
    std::cout << '\n';
    exit(1);
}


// --------------------------------------------------------------------
int main(int argc, char *argv[]) {
    if(argc != 2) {
      PrintUsage(argv[0]);
    }

    bool bInitParmParse(false);
    amrex::Initialize(argc, argv, bInitParmParse);

    int myProc(ParallelDescriptor::MyProc());

    std::string infile(argv[1]);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "argv[1] = " << argv[1] << std::endl;
      std::cout << "infile  = " << infile << std::endl;
    }

    // ---- open pltfile and get amrData reference
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData &amrData = dataServices.AmrDataRef();

    int finestLevel(amrData.FinestLevel());
    int numberOfLevels(finestLevel + 1);
    const Vector<std::string> &plotVarNames = amrData.PlotVarNames();

    // ---- print some information about the plot file
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "finestLevel = " << finestLevel << std::endl;
      std::cout << "nuberOfLevels = " << numberOfLevels << std::endl;
      for(int i(0); i < amrData.ProbDomain().size(); ++i) {
        std::cout << "ProbDomain[" << i << "] = " << amrData.ProbDomain()[i]  << std::endl;
      }
      std::cout << "NComp       = " << amrData.NComp()       << std::endl;
      for(int i(0); i < plotVarNames.size(); ++i) {
        std::cout << "plotVarNames[" << i << "] = " << plotVarNames[i] << std::endl;
      }
      for(int i(0); i < numberOfLevels; ++i) {
        //std::cout << "BoxArray[" << i << "] = " << amrData.boxArray(i) << std::endl;
      }
    }


    BoxArray fillBoxes(amrData.ProbDomain()[0]);
    fillBoxes.maxSize(128);  // ---- break the boxarray into smaller boxes
    int nVar(1), nGrow(0);
    DistributionMapping dmap(fillBoxes);
    MultiFab fillMF(fillBoxes, dmap, nVar, nGrow);  // ---- one component
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "fillBoxes = " << fillBoxes << std::endl;
      std::cout << "filling multifab for " << plotVarNames[0] << std::endl;
    }
    int fillLevelZero(0);
    int compZero(0);
    amrData.FillVar(fillMF, fillLevelZero, plotVarNames[0], compZero);
    VisMF::Write(fillMF, "fillMF");

    for(MFIter mfi(fillMF); mfi.isValid(); ++mfi) {
      FArrayBox &currentFAB = fillMF[mfi];
      const Box &fabBox = currentFAB.box();
      Real *dataPtr = currentFAB.dataPtr();  // ---- this points to the data
      usleep(myProc * msps / 10.0);  // ---- to make the output readable
      std::cout << myProc << ":: fabBox = " << fabBox << std::endl;
      std::cout << myProc << ":: dataPtr[0] = " << dataPtr[0] << std::endl;
    }
    usleep(msps);  // ---- to make the output readable

    amrex::Finalize();

    return 0;
}
// --------------------------------------------------------------------
