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

// Define show value and show value including the expected result (for ease of debugging & testing).
#define SHOWVAL(val) { amrex::Print() << #val << " = " << val << std::endl; }
#define SHOWEXPCTVAL(val,expct) { amrex::Print() << #val << " = (" << val << " ? " << expct << ")" << std::endl; }

// --------------------------------------------------------------------------
void SetFabValsToPMap(MultiFab &mf) {
  const Vector<int> &newDMA = mf.DistributionMap().ProcessorMap();
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

    amrex::system::verbose = 0;

    BL_PROFILE_VAR("main()", pmain);
    BL_PROFILE_REGION_START("main");

    amrex::Print() << std::endl << std::endl;
    amrex::Print() << " ++++++++++++++++++++++++++++++++++++++++++++++++++++ " << std::endl;
    amrex::Print() << " ++++++++++++ Beginning MoveAllFabsTest +++++++++++++ " << std::endl;

    int nProcs(ParallelDescriptor::NProcs());
    int nBoxes(8);

    // ---- make a box, then a boxarray with maxSize
    Box baseBox(IntVect(0,0,0), IntVect(maxGrid - 1, maxGrid - 1, maxGrid - 1));
    Vector<Box> boxes(nBoxes);

    for(int p(0); p < boxes.size(); ++p) {
      boxes[p] = baseBox.shift(XDIR, maxGrid);
    }
    BoxArray ba(boxes.dataPtr(), boxes.size());

    // Build 3 distribution maps, increasing, decreasing and paired.
    // =================
    Vector<int> ai0(ba.size()), ai1(ba.size()), aiMove(ba.size());
    for(int i(0); i < ai0.size(); ++i) {
      ai0[i] = i%nProcs;
      ai1[i] = (ai0.size() - 1 - i)%nProcs;
      aiMove[i] = (i % (ai0.size() / 2))%nProcs;
    }

    DistributionMapping dm0(ai0);
    DistributionMapping dm1(ai1);
    DistributionMapping dmMove(aiMove);

    amrex::Print() << "dm0: = " << dm0 << std::endl;
    amrex::Print() << "dm1: = " << dm1 << std::endl;
    amrex::Print() << "dmMove: = " << dmMove << std::endl;

    // Confirm all 3 distribution maps are different
    // =================
    SHOWEXPCTVAL(DistributionMapping::SameRefs(dm0, dm0), 1);
    SHOWEXPCTVAL(DistributionMapping::SameRefs(dm0, dm1), 0);
    SHOWEXPCTVAL(DistributionMapping::SameRefs(dm0, dmMove), 0);

    amrex::Print() << "**** dm0().linkCount = " << dm0.linkCount() << std::endl;
    amrex::Print() << "**** dm1().linkCount = " << dm1.linkCount() << std::endl;

    // Make 3 MultiFabs, with 0 & 2 being identical
    // =================
    MultiFab mf0(ba, dm0, nComp, nGhost);
    MultiFab mf1(ba, dm1, nComp, nGhost);
    MultiFab mf2(ba, dm0, nComp, nGhost);
    for(int i(0); i < mf0.nComp(); ++i) {
      mf0.setVal(static_cast<Real> (i), i, 1);
    }
    for(int i(0); i < mf1.nComp(); ++i) {
      mf1.setVal(static_cast<Real> (i + 100), i, 1);
    }
    for(int i(0); i < mf2.nComp(); ++i) {
      mf2.setVal(static_cast<Real> (i + 50), i, 1);
    }

    // Pull and test BDKeys of Multifabs before the move
    // =================
    MultiFab::BDKey bdkey0 = mf0.getBDKey();
    MultiFab::BDKey bdkey1 = mf1.getBDKey();
    MultiFab::BDKey bdkey2 = mf2.getBDKey();
    MultiFab::BDKey bdkey_ba_dm0(ba.getRefID(), dm0.getRefID());
    MultiFab::BDKey bdkey_ba_dmMove(ba.getRefID(), dmMove.getRefID());
    SHOWEXPCTVAL((bdkey0 == bdkey1), 0);
    SHOWEXPCTVAL((bdkey0 == bdkey2), 1);
    SHOWEXPCTVAL((bdkey0 == bdkey_ba_dm0), 1);
    SHOWEXPCTVAL((bdkey0 == bdkey_ba_dmMove), 0);

    amrex::Print() << "**** mf0.dm().linkCount = " << mf0.DistributionMap().linkCount() << std::endl;
    amrex::Print() << "**** mf1.dm().linkCount = " << mf1.DistributionMap().linkCount() << std::endl;
    amrex::Print() << "**** mf2.dm().linkCount = " << mf2.DistributionMap().linkCount() << std::endl;

    //amrex::Print() << "mf0.dmap = " <<  mf0.DistributionMap() << std::endl;
    //amrex::Print() << "mf1.dmap = " <<  mf1.DistributionMap() << std::endl;
    //amrex::Print() << "dmMove   = " <<  dmMove << std::endl;

    // Exchange any fab with map 0 to map move.
    // (Should change Multifabs 0 & 2, leaving 1 unchanged).
    // =================
    amrex::Print() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;
    amrex::Print() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;
    MultiFab::MoveAllFabs(dm0, dmMove);
    amrex::Print() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;
    amrex::Print() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;
 
    // Confirm new Keys are correct and print distribution maps 
    // =================
    MultiFab::BDKey bdkey0new = mf0.getBDKey();
    MultiFab::BDKey bdkey1new = mf1.getBDKey();
    SHOWEXPCTVAL((bdkey0new == bdkey1new), 0);
    SHOWEXPCTVAL((bdkey0 == bdkey0new), 0);
    SHOWEXPCTVAL((bdkey1 == bdkey1new), 1);
    SHOWEXPCTVAL((bdkey0new == bdkey_ba_dmMove), 1);
    MultiFab::BDKey bdkey_ba_dmMove_A(ba.getRefID(), dmMove.getRefID());
    SHOWEXPCTVAL((bdkey0new == bdkey_ba_dmMove_A), 1);

    amrex::Print() << "-------- after MoveAllFabs." << std::endl;
    amrex::Print() << "mf0.dmap = " <<  mf0.DistributionMap() << std::endl;
    amrex::Print() << "mf1.dmap = " <<  mf1.DistributionMap() << std::endl;
    amrex::Print() << "mf2.dmap = " <<  mf2.DistributionMap() << std::endl;

    // Now, exchange ALL fabs of appropriate size (all of them) to 
    // map 0 (which none of them should currently be). 
    // =================
    amrex::Print() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;
    amrex::Print() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;
    MultiFab::MoveAllFabs(dm0);
    amrex::Print() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;
    amrex::Print() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;

    // Confirm result
    // =================
    amrex::Print() << " =============================== " << std::endl;

    bdkey0new = mf0.getBDKey();
    bdkey1new = mf1.getBDKey();
    SHOWEXPCTVAL((bdkey0new == bdkey1new), 1);
    SHOWEXPCTVAL((bdkey0 == bdkey0new), 1);
    SHOWEXPCTVAL((bdkey1 == bdkey1new), 0);
    SHOWEXPCTVAL((bdkey0new == bdkey_ba_dmMove), 0);
    bdkey_ba_dmMove_A= MultiFab::BDKey(ba.getRefID(), dmMove.getRefID());
    SHOWEXPCTVAL((bdkey0new == bdkey_ba_dmMove_A), 0);

    amrex::Print() << "-------- after MoveAllFabs." << std::endl;
    amrex::Print() << "mf0.dmap = " <<  mf0.DistributionMap() << std::endl;
    amrex::Print() << "mf1.dmap = " <<  mf1.DistributionMap() << std::endl;
    amrex::Print() << "mf2.dmap = " <<  mf2.DistributionMap() << std::endl;

    BL_PROFILE_REGION_STOP("main");
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
