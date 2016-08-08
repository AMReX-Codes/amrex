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
#include <FabConv.H>

using std::cout;
using std::cerr;
using std::endl;


void DirectoryTests();
void FileTests();
void TestWriteNFiles(int nfiles, int maxgrid, int ncomps, int nboxes,
                     bool raninit, bool mb2);
void TestWriteNFilesRawNative(int nfiles, int maxgrid, int ncomps,
                              int nboxes, bool raninit, bool mb2,
			      VisMF::Header::Version writeMinMax,
			      bool groupsets, bool setbuf);
void TestReadMF(const std::string &mfName);
void NFileTests(int nOutFiles, const std::string &filePrefix);


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
    cout << "   [rbuffsize = rbs]" << '\n';
    cout << "   [wbuffsize = wbs]" << '\n';
    cout << "   [writeminmax = wmm]" << '\n';
    cout << "   [writefaminmax = wfamm]" << '\n';
    cout << "   [groupsets = groupsets]" << '\n';
    cout << "   [setbuf = setbuf]" << '\n';
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
  int rbs(8192), wbs(8192);
  bool raninit(false), mb2(false);
  bool writeminmax(false), writefaminmax(false);
  bool groupsets(false), setbuf(true);
  VisMF::Header::Version hVersion;

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

  pp.query("writeminmax", writeminmax);
  pp.query("writefaminmax", writefaminmax);
  BL_ASSERT( ! (writeminmax && writefaminmax));
  if(writeminmax) {
    hVersion = VisMF::Header::RawNativeMinMax_v1;
  } else if(writefaminmax) {
    hVersion = VisMF::Header::RawNativeFAMinMax_v1;
  } else {
    hVersion = VisMF::Header::RawNative_v1;
  }

  pp.query("groupsets", groupsets);
  pp.query("setbuf", setbuf);

  pp.query("rbuffsize", rbs);
  pp.query("wbuffsize", wbs);
  RealDescriptor::SetReadBufferSize(rbs);
  RealDescriptor::SetWriteBufferSize(wbs);

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
    cout << "rbuffsize = " << rbs << endl;
    cout << "wbuffsize = " << wbs << endl;
    cout << "writeminmax = " << writeminmax << endl;
    cout << "writefaminmax = " << writefaminmax << endl;
    cout << "groupsets = " << groupsets << endl;
    cout << "setbuf = " << setbuf << endl;
    cout << endl;
    cout << "sizeof(int) = " << sizeof(int) << endl;
    cout << "sizeof(size_t) = " << sizeof(size_t) << endl;
    cout << "sizeof(long) = " << sizeof(long) << endl;
    cout << "sizeof(long long) = " << sizeof(long long) << endl;
    cout << "sizeof(std::streampos) = " << sizeof(std::streampos) << endl;
    cout << "sizeof(std::streamoff) = " << sizeof(std::streamoff) << endl;
    cout << "sizeof(std::streamsize) = " << sizeof(std::streamsize) << endl;
    cout << endl;
    cout << "std::numeric_limits<int>::min()  = " << std::numeric_limits<int>::min() << endl;
    cout << "std::numeric_limits<int>::max()  = " << std::numeric_limits<int>::max() << endl;
    cout << "std::numeric_limits<Real>::min() = " << std::numeric_limits<Real>::min() << endl;
    cout << "std::numeric_limits<Real>::max() = " << std::numeric_limits<Real>::max() << endl;
    cout << endl;
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

/*
  for(int itimes(0); itimes < ntimes; ++itimes) {
    if(ParallelDescriptor::IOProcessor()) {
      cout << endl << "--------------------------------------------------" << endl;
      cout << "Testing NFile Operations" << endl;
    }

    int nOutFiles(5);
    std::string filePrefix("FiveFiles");
    NFileTests(nOutFiles, filePrefix);

    if(ParallelDescriptor::IOProcessor()) {
      cout << "==================================================" << endl;
      cout << endl;
    }
  }
*/


/*
  for(int itimes(0); itimes < ntimes; ++itimes) {
    if(ParallelDescriptor::IOProcessor()) {
      cout << endl << "--------------------------------------------------" << endl;
      cout << "Testing File Operations" << endl;
    }

    FileTests();

    if(ParallelDescriptor::IOProcessor()) {
      cout << "==================================================" << endl;
      cout << endl;
    }
  }
*/


/*
  for(int itimes(0); itimes < ntimes; ++itimes) {
    if(ParallelDescriptor::IOProcessor()) {
      cout << endl << "--------------------------------------------------" << endl;
      cout << "Testing Directory Operations" << endl;
    }

    DirectoryTests();

    if(ParallelDescriptor::IOProcessor()) {
      cout << "==================================================" << endl;
      cout << endl;
    }
  }
*/


/*
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
*/


/*
  for(int itimes(0); itimes < ntimes; ++itimes) {
    if(ParallelDescriptor::IOProcessor()) {
      cout << endl << "--------------------------------------------------" << endl;
      cout << "Testing NFiles Raw Native Write" << endl;
    }

    TestWriteNFilesRawNative(nfiles, maxgrid, ncomps, nboxes, raninit, mb2,
                             hVersion, groupsets, setbuf);

    if(ParallelDescriptor::IOProcessor()) {
      cout << "==================================================" << endl;
      cout << endl;
    }
  }
*/


  for(int itimes(0); itimes < ntimes; ++itimes) {
    if(ParallelDescriptor::IOProcessor()) {
      cout << endl << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      cout << "Testing MF Read" << endl;
    }

    TestReadMF("TestMF");
    TestReadMF("TestMFRawNative");
    TestReadMF("TestMFRawNativeMinMax");
    TestReadMF("TestMFRawNativeFAMinMax");

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
