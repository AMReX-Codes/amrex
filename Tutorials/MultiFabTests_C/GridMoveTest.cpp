// --------------------------------------------------------------------------
// GridMoveTest.cpp
// --------------------------------------------------------------------------
//  this file tests the performance for moving grids
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

#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

const int maxGrid(64);
const int nComp(128);
const int nGhost(0);
const int XDIR(0);
const int nTimes(16);

// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {

    BoxLib::Initialize(argc,argv);    

    BL_PROFILE_VAR("main()", pmain);
    BL_PROFILE_REGION_START("main");

    int nProcs(ParallelDescriptor::NProcs());

    // ---- make a box, then a boxarray with maxSize
    Box baseBox(IntVect(0,0,0), IntVect(maxGrid - 1, maxGrid - 1, maxGrid - 1));
    Array<Box> boxes(nProcs);

    for(int p(0); p < boxes.size(); ++p) {
      boxes[p] = baseBox.shift(XDIR, maxGrid);
    }
    BoxArray ba(boxes.dataPtr(), boxes.size());
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "ba.size() = " << ba.size() << std::endl;
      std::cout << "ba[0] = " << ba[0] << std::endl;
      std::cout << "ba[1] = " << ba[1] << std::endl;
      std::cout << "ba[last] = " << ba[ba.size() - 1] << std::endl;
    }

    // ---- make a multifab, setval to the index
    MultiFab mf(ba, nComp, nGhost);
    for(int i(0); i < mf.nComp(); ++i) {
      mf.setVal(static_cast<Real> (i), i, 1);
    }

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

    // ------------------ copy entire multifab to nProcs/2) % nProcs
    {
      // ---- initialize it to (oldmap + (nProcs/2)) % nProcs
      Array<int> newMap(procMap.size());
      for(int i(0); i < procMap.size(); ++i) {
        newMap[i] = (procMap[i] + (nProcs/2)) % nProcs;
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
      BL_PROFILE_REGION_START("MFCopy_ProcPlusNPd2");
      BL_PROFILE_VAR("MFCopy_ProcPlusNPd2", mfcppnpd2);
      for(int n(0); n < nTimes; ++n) {
        mfNewMap.copy(mf);
      }
      BL_PROFILE_VAR_STOP(mfcppnpd2);
      BL_PROFILE_REGION_STOP("MFCopy_ProcPlusNPd2");
    }

    // ------------------ copy a random set of grids elsewhere
    {
      int nCopies(nProcs/10);
      nCopies = std::min(nCopies, nProcs/2);
      int nRanks(nCopies * 2);

      Array<int> copyArray;
      if(ParallelDescriptor::IOProcessor()) {
        std::set<int> copySet;  // ---- a unique set of random numbers
        while(copySet.size() < nRanks) {
          int r(BoxLib::Random_int(nProcs));
	  if(copySet.find(r) == copySet.end()) {
	    copySet.insert(r);
	    copyArray.push_back(r);
	  }
        }
        for(int i(0); i < copyArray.size(); ++i) {
          std::cout << "copyArray[" << i << "]  = " << copyArray[i] << std::endl;
	}
      } else {
        copyArray.resize(nRanks);
      }

      ParallelDescriptor::Barrier();
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "Starting random copy." << std::endl;
      }
      ParallelDescriptor::Bcast(copyArray.dataPtr(), copyArray.size());
        
      Array<int> newMap(nCopies + 1);
      BoxArray baCopy(nCopies);
      for(int i(0); i < nCopies; ++i) {
        newMap[i] = copyArray[i];
	baCopy.set(i, ba[copyArray[i + nCopies]]);
      }
      newMap[nCopies] = ParallelDescriptor::MyProc();

      DistributionMapping newDMap(newMap);
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "newDMap = " << newDMap << std::endl;
      }

      // ---- make a new multifab with the new map and copy from mf
      MultiFab mfNewMap;
      mfNewMap.define(baCopy, nComp, nGhost, newDMap, Fab_allocate);
      mfNewMap.setVal(-42.0);

      // ---- now copy from mf
      BL_PROFILE_REGION_START("MFCopy_RandomSet");
      BL_PROFILE_VAR("MFCopy_RandomSet", mfcprs);
      for(int n(0); n < nTimes; ++n) {
        mfNewMap.copy(mf);
      }
      BL_PROFILE_VAR_STOP(mfcprs);
      BL_PROFILE_REGION_STOP("MFCopy_RandomSet");
    }


    BL_PROFILE_REGION_STOP("main");
    BL_PROFILE_VAR_STOP(pmain);

    BoxLib::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
