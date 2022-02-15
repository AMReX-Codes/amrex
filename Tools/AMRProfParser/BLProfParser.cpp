// ----------------------------------------------------------------------
//  BLProfParser.cpp
// ----------------------------------------------------------------------
#include <iostream>

using std::cout;
using std::endl;

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>

using namespace amrex;

extern void PrintProfParserBatchUsage(std::ostream &os);
extern bool ProfParserBatchFunctions(int argc, char *argv[], bool runDefault,
                                     bool &bParserProf);


// ----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv, false);

  BL_PROFILE_VAR("main()", pmain);

  bool bIOP(ParallelDescriptor::IOProcessor());

  if(argc < 2) {
    if(bIOP) {
      cout << '\n';
      cout << "Usage:  " << argv[0] << "  [options]  profddirname\n";
      PrintProfParserBatchUsage(cout);
      cout << endl;
    }
    BL_PROFILE_VAR_STOP(pmain);
    amrex::BLProfiler::SetNoOutput();
    amrex::Finalize();
    return(-1);
  }


  bool runDefault(true), bParserProf(false);

  ProfParserBatchFunctions(argc, argv, runDefault, bParserProf);


  BL_PROFILE_VAR_STOP(pmain);

  if( ! bParserProf) {
    amrex::BLProfiler::SetNoOutput();
  }
  if(bIOP) {
    cout << endl;
  }

  amrex::Finalize();

  return(0);
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
