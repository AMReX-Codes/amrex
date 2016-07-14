// --------------------------------------------------------------------------
// BcastClasses.cpp
// --------------------------------------------------------------------------
//  this file tests functions to broadcast and serialize classes.
// --------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <BoxArray.H>
#include <MultiFab.H>
#include <ParallelDescriptor.H>
using std::cout;
using std::endl;


// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    BoxLib::Initialize(argc, argv);
    
    bool bIOP(ParallelDescriptor::IOProcessor());
    int  myProc(ParallelDescriptor::MyProc());
    int  ioProcNum(ParallelDescriptor::IOProcessorNumber());

    Box bControl(IntVect(0,0,0), IntVect(63,63,63));
    Box bBcast;

    Array<int> aSB;
    if(bIOP) {
      aSB = BoxLib::SerializeBox(bControl);
    }
    
    BoxLib::BroadcastArray(aSB, myProc, ioProcNum, ParallelDescriptor::CommunicatorAll());

    ParallelDescriptor::Barrier();
    BoxLib::USleep(myProc/10.0);
    for(int i(0); i < aSB.size(); ++i) {
      cout << myProc << "::  aSB[" << i << "] = " << aSB[i] << endl;
    }

    Box uSB = BoxLib::UnSerializeBox(aSB);
    ParallelDescriptor::Barrier();
    BoxLib::USleep(myProc/10.0);
    cout << myProc << "::  uSB = " << uSB << endl;

    ParallelDescriptor::Barrier();
    BoxArray baControl(bControl);
    baControl.maxSize(32);
    BoxArray baBcast;
    if(bIOP) {
      baBcast = baControl;
      cout << "baBcast = " << baBcast << endl;
    }

    Array<int> aBASerial;
    if(bIOP) {
      aBASerial = BoxLib::SerializeBoxArray(baControl);
      for(int i(0); i < aSB.size(); ++i) {
        cout << myProc << "::  aBASerial[" << i << "] = " << aBASerial[i] << endl;
      }
    }
    ParallelDescriptor::Barrier();

    BoxLib::BroadcastBoxArray(baBcast, myProc, ioProcNum,
                              ParallelDescriptor::CommunicatorAll());

    ParallelDescriptor::Barrier();
    BoxLib::USleep(myProc/10.0);
    cout << myProc << "::  baBcast = " << baBcast << endl;
    if(baBcast != baControl) {
      cout << myProc << "::  **** Error:  bad BoxArrayBroadcast:  baBcast baControl = "
           << baBcast << "  " << baControl << endl;
    }

/*
    if( ! fba.ok()) {
        BoxLib::Error("BoxArray is not OK");
    }

    if( ! fba.isDisjoint()) {
      BoxLib::Error("BoxArray is not disjoint");
    }

    fba.maxSize(32);

    if(bIOP) {
      std::cout << "number of grids in fba: " << fba.size() << '\n';
    }
*/



    BoxLib::Finalize();

    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
