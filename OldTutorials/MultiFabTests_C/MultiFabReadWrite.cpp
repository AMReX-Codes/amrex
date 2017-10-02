// --------------------------------------------------------------------------
// MultiFabReadWrite.cpp
// --------------------------------------------------------------------------
//  this file writes and reads multifabs.
// --------------------------------------------------------------------------

#include <new>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <unistd.h>

#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

const int maxGrid(8);
const int pdHi(63);
const int nComp(32);
const int nGhost(0);
const int nFiles(64);

// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {

    amrex::Initialize(argc,argv);    

    BL_PROFILE_VAR("main()", pmain);
    BL_PROFILE_REGION_START("main");

    VisMF::SetNOutFiles(nFiles);  // ---- this will enforce the range [1, nprocs]

    // ---- make a box, then a boxarray with maxSize
    Box bDomain(IntVect(0,0,0), IntVect(pdHi,pdHi,pdHi));
    BoxArray ba(bDomain);
    ba.maxSize(maxGrid);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "ba = " << ba << std::endl;
    }

    DistributionMapping dmap{ba};

    // ---- make a multifab, set interior to the index
    // ---- set the ghost regions to -index
    std::string outfile = "MF_Out";
    MultiFab mf(ba, dmap, nComp, nGhost);
    for(int i(0); i < mf.nComp(); ++i) {
      mf.setVal(static_cast<Real> (i), i, 1);
      mf.setBndry(static_cast<Real> (-i), i, 1);
    }
    VisMF::Write(mf, outfile);  // ---- write the multifab to MF_Out

    // uncomment the following to write out each fab
    /*
    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const int index(mfi.index());
      FArrayBox &fab = mf[mfi];
      std::string fname = amrex::Concatenate("FAB_", index, 4);
      std::ofstream fabs(fname.c_str());
      fab.writeOn(fabs);
      fabs.close();
    }
    */

    // ---- make a new multifab and read in the one we just wrote
    std::string infile(outfile);
    MultiFab mfInOut;
    VisMF::Read(mfInOut, infile);

    std::string inoutfile = "MF_InOut";
    BL_PROFILE_REGION_START("VisMF::Write()");
    VisMF::Write(mfInOut, inoutfile);
    BL_PROFILE_REGION_STOP("VisMF::Write()");

    // ---- make a new distribution map
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "dmap = " << dmap << std::endl;
      std::cout << "dmap.size() = " << dmap.size() << std::endl;
    }

    const Array<int> procMap = dmap.ProcessorMap();
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "procMap.size() = " << procMap.size() << std::endl;
    }

    // ---- initialize it to (oldmap + 1) % nprocs
    Array<int> newMap(procMap.size());
    for(int i(0); i < procMap.size(); ++i) {
      newMap[i] = (procMap[i] + 1) % ParallelDescriptor::NProcs();
    }
    DistributionMapping newDMap(newMap);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "newDMap = " << newDMap << std::endl;
    }

    // ---- make a new multifab with the new map and copy from mf
    MultiFab mfNewMap(ba, newDMap, nComp, nGhost);
    mfNewMap.setVal(-42.0);

    // ---- now copy from mf
    BL_PROFILE_REGION_START("MFCopy");
    mfNewMap.copy(mf);
    BL_PROFILE_REGION_STOP("MFCopy");

    std::string mfnmoutfile = "MF_NewMap";
    VisMF::Write(mfNewMap, mfnmoutfile);

    // ---- all three multifabs should be the same

    // ---- this will fill the ghost regions from intersecting fabs
    mf.FillBoundary();

    BL_PROFILE_REGION_STOP("main");
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
