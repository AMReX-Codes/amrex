// --------------------------------------------------------------------------
// MoveAllFabsTest.cpp
// --------------------------------------------------------------------------
//  this file tests the performance for moving grids
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

using std::cout;
using std::endl;

using namespace amrex;

const int maxGrid(32);
const int nComp(8);
const int nGhost(1);
const int XDIR(0);
const int nTimes(4);

#define SHOWVAL(val) { amrex::Print() << #val << " = " << val << std::endl; }

// --------------------------------------------------------------------------
void SetFabValsToPMap(MultiFab &mf) {
  const Array<int> &newDMA = mf.DistributionMap().ProcessorMap();
  for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const int index(mfi.index());
    FArrayBox &fab = mf[mfi];
    for(int i(0); i < fab.nComp(); ++i) {  // ---- setVal to distmap index
      Real val(newDMA[index] + ((static_cast<Real> (i)) / fab.nComp()));
      fab.setVal(val, i);
    }
  }
}


// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {

    amrex::Initialize(argc,argv);    

    BL_PROFILE_VAR("main()", pmain);
    BL_PROFILE_REGION_START("main");

    //int nProcs(ParallelDescriptor::NProcs());
    int nBoxes(8);

    // ---- make a box, then a boxarray with maxSize
    Box baseBox(IntVect(0,0,0), IntVect(maxGrid - 1, maxGrid - 1, maxGrid - 1));
    Array<Box> boxes(nBoxes);

    for(int p(0); p < boxes.size(); ++p) {
      boxes[p] = baseBox.shift(XDIR, maxGrid);
    }
    BoxArray ba(boxes.dataPtr(), boxes.size());

    Array<int> ai0(ba.size()), ai1(ba.size()), aiMove(ba.size());
    for(int i(0); i < ai0.size(); ++i) {
      ai0[i] = i;
      ai1[i] = ai0.size() - 1 - i;
      aiMove[i] = i % (ai0.size() / 2);
    }

    DistributionMapping dm0(ai0);
    DistributionMapping dm1(ai1);
    DistributionMapping dmMove(aiMove);

    SHOWVAL(DistributionMapping::SameRefs(dm0, dm1));
    SHOWVAL(DistributionMapping::SameRefs(dm0, dmMove));

    MultiFab mf0(ba, dm0, nComp, nGhost);
    MultiFab mf1(ba, dm1, nComp, nGhost);
    //MultiFab mf0(ba, dmMove, nComp, nGhost);
    //MultiFab mf1(ba, dmMove, nComp, nGhost);
    for(int i(0); i < mf0.nComp(); ++i) {
      mf0.setVal(static_cast<Real> (i), i, 1);
    }
    for(int i(0); i < mf1.nComp(); ++i) {
      mf1.setVal(static_cast<Real> (i + 100), i, 1);
    }

    MultiFab::BDKey bdkey0 = mf0.getBDKey();
    MultiFab::BDKey bdkey1 = mf1.getBDKey();
    MultiFab::BDKey bdkey_ba_dm0(ba.getRefID(), dm0.getRefID());
    MultiFab::BDKey bdkey_ba_dmMove(ba.getRefID(), dmMove.getRefID());
    SHOWVAL((bdkey0 == bdkey1));
    SHOWVAL((bdkey0 == bdkey_ba_dm0));
    SHOWVAL((bdkey0 == bdkey_ba_dmMove));


    amrex::Print() << "**** mf0.dm().linkCount = " << mf0.DistributionMap().linkCount() << std::endl;
    amrex::Print() << "**** mf1.dm().linkCount = " << mf1.DistributionMap().linkCount() << std::endl;

    //amrex::Print() << "mf0.dmap = " <<  mf0.DistributionMap() << std::endl;
    //amrex::Print() << "mf1.dmap = " <<  mf1.DistributionMap() << std::endl;
    //amrex::Print() << "dmMove   = " <<  dmMove << std::endl;

    //MultiFab::MoveAllFabs(dmMove.ProcessorMap());
    MultiFab::MoveAllFabs(dmMove);

    mf0.updateBDKey();
    mf1.updateBDKey();

    MultiFab::BDKey bdkey0new = mf0.getBDKey();
    MultiFab::BDKey bdkey1new = mf1.getBDKey();
    SHOWVAL((bdkey0new == bdkey1new));
    SHOWVAL((bdkey0 == bdkey0new));
    SHOWVAL((bdkey1 == bdkey1new));
    SHOWVAL((bdkey0new == bdkey_ba_dmMove));
    MultiFab::BDKey bdkey_ba_dmMove_A(ba.getRefID(), dmMove.getRefID());
    SHOWVAL((bdkey0new == bdkey_ba_dmMove_A));

    amrex::Print() << "-------- after MoveAllFabs." << std::endl;
    amrex::Print() << "mf0.dmap = " <<  mf0.DistributionMap() << std::endl;
    amrex::Print() << "mf1.dmap = " <<  mf1.DistributionMap() << std::endl;


    BL_PROFILE_REGION_STOP("main");
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
