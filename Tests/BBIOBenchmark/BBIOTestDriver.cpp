// -------------------------------------------------------------
// BBIOTestDriver.cpp
// -------------------------------------------------------------
#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
using std::ios;

#include <unistd.h>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>

using std::cout;
using std::cerr;
using std::endl;

using namespace amrex;

void TestWriteNFiles(int nfiles, int nMB, bool raninit, bool mb2);
void TestReadNFiles(int nfiles, int nMB, bool raninit, bool mb2);
void SetDirName(const std::string &dirname);


// -------------------------------------------------------------
static void PrintUsage(const char *progName) {
  if(ParallelDescriptor::IOProcessor()) {
    cout << '\n';
    cout << "Usage:" << '\n';
    cout << progName << '\n';
    cout << "   [nfiles = nfiles]" << '\n';
    cout << "   [nMB = nMB]" << '\n';
    cout << "   [nsleep = nsleep]" << '\n';
    cout << "   [ntimes = ntimes]" << '\n';
    cout << "   [raninit = tf]" << '\n';
    cout << "   [mb2    = tf]" << '\n';
    cout << "   [dirName = dirname]" << '\n';
    cout << '\n';
    cout << "Running with default values." << '\n';
    cout << '\n';
  }
}


// -------------------------------------------------------------
int main(int argc, char *argv[]) {

  amrex::Initialize(argc,argv);
  VisMF::Initialize();

  if(argc == 1) {
    PrintUsage(argv[0]);
  }

  ParmParse pp;

  int myproc(ParallelDescriptor::MyProc());
  int nprocs(ParallelDescriptor::NProcs());
  int nsleep(0), nfiles(std::min(nprocs, 128));  // limit default to max of 128
  int nMB(100), ntimes(1);
  bool raninit(false), mb2(false);
  std::string dirName(".");

  pp.query("nfiles", nfiles);
  nfiles = std::max(1, std::min(nfiles, nprocs));

  pp.query("nMB", nMB);
  nMB = std::max(1, nMB);

  pp.query("ntimes", ntimes);
  ntimes = std::max(1, ntimes);

  pp.query("raninit", raninit);
  pp.query("mb2", mb2);

  pp.query("dirName", dirName);
  std::cout << "dirName = " << dirName << std::endl;
  SetDirName(dirName);

  if(ParallelDescriptor::IOProcessor()) {
    cout << endl;
    cout << "**************************************************" << endl;
    cout << "nprocs = " << nprocs << endl;
    cout << "nfiles = " << nfiles << endl;
    cout << "nMB    = " << nMB    << endl;
    cout << "ntimes = " << ntimes << endl;
    cout << "raninit = " << raninit << endl;
    cout << "mb2 = " << mb2 << endl;
  }

  pp.query("nsleep", nsleep);
  if(nsleep > 0) {  // test the timer
    double timerTimeStart = ParallelDescriptor::second();
    sleep(nsleep);  // for attaching a debugger or testing the timer
    double timerTime = ParallelDescriptor::second() - timerTimeStart;
    cout << "  ----- " << myproc << " :  " << "Sleep time = "
         << timerTime << "  (should be " << nsleep << " seconds)" << endl;
  }

  ParallelDescriptor::Barrier();

  for(int itimes(0); itimes < ntimes; ++itimes) {
    if(ParallelDescriptor::IOProcessor()) {
      cout << endl << "--------------------------------------------------" << endl;
      cout << "Testing NFiles Write" << endl;
    }

    TestWriteNFiles(nfiles, nMB, raninit, mb2);

    if(ParallelDescriptor::IOProcessor()) {
      cout << "==================================================" << endl;
      cout << endl;
    }
  }

  for(int itimes(0); itimes < ntimes; ++itimes) {
    if(ParallelDescriptor::IOProcessor()) {
      cout << endl << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      //cout << "Testing NFiles Read" << endl;
      cout << "Testing NFiles Write again" << endl;
    }

    TestWriteNFiles(nfiles, nMB, raninit, mb2);

    if(ParallelDescriptor::IOProcessor()) {
      cout << "##################################################" << endl;
      cout << endl;
    }
  }


  amrex::Finalize();
  return 0;
}
// -------------------------------------------------------------
// -------------------------------------------------------------
