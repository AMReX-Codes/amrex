// -------------------------------------------------------------
// IOTestDriver.cpp
// -------------------------------------------------------------
#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
using std::ios;

#include <unistd.h>

#include <ParallelDescriptor.H>
#include <Utility.H>
#include <ParmParse.H>
#include <MultiFab.H>
#include <VisMF.H>

using std::cout;
using std::cerr;
using std::endl;


void TestWriteNFiles(int nfiles, int maxgrid, int ncomps, int nboxes,
                     bool raninit, bool mb2);
void TestReadMF();


// -------------------------------------------------------------
static void PrintUsage(const char *progName) {
  if(ParallelDescriptor::IOProcessor()) {
    cout << '\n';
    cout << "Usage:" << '\n';
    cout << progName << '\n';
    cout << "   [nfiles = nfiles]" << '\n';
    cout << "   [maxgrid = maxgrid]" << '\n';
    cout << "   [ncomps = ncomps]" << '\n';
    cout << "   [nboxes = nboxes]" << '\n';
    cout << "   [nsleep = nsleep]" << '\n';
    cout << "   [ntimes = ntimes]" << '\n';
    cout << "   [raninit = tf]" << '\n';
    cout << "   [mb2    = tf]" << '\n';
    cout << '\n';
    cout << "Running with default values." << '\n';
    cout << '\n';
  }
}


// -------------------------------------------------------------
int main(int argc, char *argv[]) {

  BoxLib::Initialize(argc,argv);
  VisMF::Initialize();

  if(argc == 1) {
    PrintUsage(argv[0]);
  }

  ParmParse pp;

  int myproc(ParallelDescriptor::MyProc());
  int nprocs(ParallelDescriptor::NProcs());
  int nsleep(0), nfiles(std::min(nprocs, 128));  // limit default to max of 128
  int maxgrid(32), ncomps(4), nboxes(nprocs), ntimes(1);
  bool raninit(false), mb2(false);

  pp.query("nfiles", nfiles);
  nfiles = std::max(1, std::min(nfiles, nprocs));

  pp.query("maxgrid", maxgrid);
  maxgrid = std::max(4, std::min(maxgrid, 256));

  pp.query("ncomps", ncomps);
  ncomps = std::max(1, std::min(ncomps, 256));

  pp.query("nboxes", nboxes);
  nboxes = std::max(1, nboxes);

  pp.query("ntimes", ntimes);
  ntimes = std::max(1, ntimes);

  pp.query("raninit", raninit);
  pp.query("mb2", mb2);

  if(ParallelDescriptor::IOProcessor()) {
    cout << endl;
    cout << "**************************************************" << endl;
    cout << "nprocs = " << nprocs << endl;
    cout << "nfiles = " << nfiles << endl;
    cout << "maxgrid = " << maxgrid << endl;
    cout << "ncomps = " << ncomps << endl;
    cout << "nboxes = " << nboxes << endl;
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

    TestWriteNFiles(nfiles, maxgrid, ncomps, nboxes, raninit, mb2);

    if(ParallelDescriptor::IOProcessor()) {
      cout << "==================================================" << endl;
      cout << endl;
    }
  }

  for(int itimes(0); itimes < ntimes; ++itimes) {
    if(ParallelDescriptor::IOProcessor()) {
      cout << endl << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      cout << "Testing MF Read" << endl;
    }

    TestReadMF();

    if(ParallelDescriptor::IOProcessor()) {
      cout << "##################################################" << endl;
      cout << endl;
    }
  }


  BoxLib::Finalize();
  return 0;
}
// -------------------------------------------------------------
// -------------------------------------------------------------
