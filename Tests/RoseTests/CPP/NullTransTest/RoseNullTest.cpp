// --------------------------------------------------------------------------
// MultiFabReadWrite.cpp
// --------------------------------------------------------------------------
//  this file writes and reads multifabs.
// --------------------------------------------------------------------------
#include <winstd.H>

#include <new>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#ifndef WIN32
#include <unistd.h>
#endif

#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <VisMF.H>
#include <LevelBld.H>

#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

const int maxGrid(64);
const int pdHi(127);
const int nComp(5);
const int nGhost(2);
const int nFiles(64);

LevelBld *getLevelBld() {
  return 0;
}


// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {

    BoxLib::Initialize(argc,argv);    

    VisMF::SetNOutFiles(nFiles);  // ---- this will enforce the range [1, nprocs]

    // ---- make a box, then a boxarray with maxSize
    Box bDomain(IntVect(0,0,0), IntVect(pdHi,pdHi,pdHi));
    BoxArray ba(bDomain);
    ba.maxSize(maxGrid);
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "ba = " << ba << std::endl;
    }

    // ---- make a multifab, set interior to the index
    // ---- set the ghost regions to -index
    std::string outfile = "MF_Out";
    MultiFab mf(ba, nComp, nGhost);
    for(int i(0); i < mf.nComp(); ++i) {
      mf.setVal(static_cast<Real> (i), i, 1);
      mf.setBndry(static_cast<Real> (-i), i, 1);
    }
    VisMF::Write(mf, outfile);  // ---- write the multifab to MF_Out

    // uncomment the following to write out each fab
    /*
    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const int index(mfi.index());
      FArrayBox &fab = mf[index];
      std::string fname = BoxLib::Concatenate("FAB_", index, 4);
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
    VisMF::Write(mfInOut, inoutfile);

    // ---- make a new distribution map
    DistributionMapping dmap(mf.DistributionMap());
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "dmap = " << dmap << std::endl;
      std::cout << "dmap.size() = " << dmap.size() << std::endl;
    }

    // ------------------------------------------------------------------
    // ----- very important:  here we are copying a procmap,
    // -----                  but if you just make your own Array<int>
    // -----                  it must have an extra value at the end
    // -----                  set to ParallelDescriptor::MyProc()
    // -----                  see DistributionMapping.H
    // ------------------------------------------------------------------
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
    MultiFab mfNewMap;
    mfNewMap.define(ba, nComp, nGhost, newDMap, Fab_allocate);
    mfNewMap.setVal(-42.0);

    // ---- now copy from mf
    mfNewMap.copy(mf);

    std::string mfnmoutfile = "MF_NewMap";
    VisMF::Write(mfNewMap, mfnmoutfile);

    // ---- all three multifabs should be the same

    // ---- this will fill the ghost regions from intersecting fabs
    mf.FillBoundary();

    BoxLib::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
