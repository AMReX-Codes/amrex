// --------------------------------------------------------------------------
// BcastClasses.cpp
// --------------------------------------------------------------------------
//  this file tests functions to broadcast and serialize classes.
// --------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
using std::cout;
using std::endl;


// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    amrex::Initialize(argc, argv);
    
    bool bIOP(amrex::ParallelDescriptor::IOProcessor());
    int  myProc(amrex::ParallelDescriptor::MyProc());
    int  ioProcNum(amrex::ParallelDescriptor::IOProcessorNumber());

    amrex::Box bControl(amrex::IntVect(0,0,0), amrex::IntVect(63,63,63));
    amrex::Box bBcast;

    // BroadcastBox Test ---------------------------------------------

    amrex::Vector<int> aSB;
    if(bIOP) {
      aSB = amrex::SerializeBox(bControl);
    }
    
    amrex::BroadcastArray(aSB, myProc, ioProcNum, amrex::ParallelDescriptor::CommunicatorAll());

    amrex::ParallelDescriptor::Barrier();
    amrex::USleep(myProc/10.0);
    for(int i(0); i < aSB.size(); ++i) {
      cout << myProc << "::  aSB[" << i << "] = " << aSB[i] << endl;
    }

    amrex::Box uSB = amrex::UnSerializeBox(aSB);
    amrex::ParallelDescriptor::Barrier();
    amrex::USleep(myProc/10.0);
    cout << myProc << "::  uSB = " << uSB << endl;

    // BroadcastBoxArray Test -----------------------------------------

    amrex::ParallelDescriptor::Barrier();
    amrex::BoxArray baControl(bControl);
    baControl.maxSize(32);
    amrex::BoxArray baBcast;
    if(bIOP) {
      baBcast = baControl;
      cout << "baBcast = " << baBcast << endl;
    }

    amrex::Vector<int> aBASerial;
    if(bIOP) {
      aBASerial = amrex::SerializeBoxArray(baControl);
      for(int i(0); i < aSB.size(); ++i) {
        cout << myProc << "::  aBASerial[" << i << "] = " << aBASerial[i] << endl;
      }
    }
    amrex::ParallelDescriptor::Barrier();

    amrex::BroadcastBoxArray(baBcast, myProc, ioProcNum,
                              amrex::ParallelDescriptor::CommunicatorAll());

    amrex::ParallelDescriptor::Barrier();
    amrex::USleep(myProc/10.0);
    cout << myProc << "::  baBcast = " << baBcast << endl;
    if(baBcast != baControl) {
      cout << myProc << "::  **** Error:  bad BoxArrayBroadcast:  baBcast baControl = "
           << baBcast << "  " << baControl << endl;
    }

    amrex::ParallelDescriptor::Barrier();
    if(bIOP) {
      cout << endl;
    }

   // BroadcastString Test -----------------------------------------

    amrex::ParallelDescriptor::Barrier();
    std::string strTest("");
    std::string strControl("This is a string test.");
    if(bIOP) {
      strTest = strControl;
      cout << "strTest = " << strTest << endl;
    }

    amrex::ParallelDescriptor::Barrier();

    amrex::BroadcastString(strTest, myProc, ioProcNum,
                              amrex::ParallelDescriptor::CommunicatorAll());

    amrex::ParallelDescriptor::Barrier();
    amrex::USleep(myProc/10.0);
    cout << myProc << "::  strTest = " << strTest << endl;
    if(strTest != strControl) {
      cout << myProc << "::  **** Error:  bad BroadcastString:  strTest strControl = "
           << strTest << "  " << strControl << endl;
    }

    amrex::ParallelDescriptor::Barrier();
    if(bIOP) {
      cout << endl;
    }

    // BroadcastStringArray Test -----------------------------------------

    amrex::ParallelDescriptor::Barrier();
    amrex::Vector<std::string> strArrayTest;
    amrex::Vector<std::string> strArrayControl;
    strArrayControl.resize(5);
    strArrayControl[0] = "This is a string array test.";
    strArrayControl[1] = "Do not adjust your Linux settings.";
    strArrayControl[2] = "AMReX is in control.";
    strArrayControl[3] = "We are now in control of the Regridding,";
    strArrayControl[4] = "and the LoadBalancing.";

    if(bIOP) {
      strArrayTest = strArrayControl;
      for(int i(0); i<strArrayTest.size(); ++i)
      {
        cout << strArrayTest[i] << endl;
      }
    }

    amrex::ParallelDescriptor::Barrier();

    amrex::BroadcastStringArray(strArrayTest, myProc, ioProcNum,
                                amrex::ParallelDescriptor::CommunicatorAll());

    amrex::ParallelDescriptor::Barrier();
    amrex::USleep(myProc/10.0);
    for(int i(0); i<strArrayTest.size(); ++i)
    {
      cout << myProc << ":: strArrayTest[i] =" << strArrayTest[i] << endl;
    }
    if(strArrayTest != strArrayControl) {
      cout << myProc << "::  **** Error:  bad strArrayBroadcast:  strArrayTest strArrayControl";
    }

/*
    if( ! fba.ok()) {
        amrex::Error("BoxArray is not OK");
    }

    if( ! fba.isDisjoint()) {
      amrex::Error("BoxArray is not disjoint");
    }

    fba.maxSize(32);

    if(bIOP) {
      std::cout << "number of grids in fba: " << fba.size() << '\n';
    }
*/

    amrex::Finalize();

    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
