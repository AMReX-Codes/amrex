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
                     bool raninit, bool mb2,
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
    cout << "   [nfiles            = nfiles   ]" << '\n';
    cout << "   [maxgrid           = maxgrid  ]" << '\n';
    cout << "   [ncomps            = ncomps   ]" << '\n';
    cout << "   [nboxes            = nboxes   ]" << '\n';
    cout << "   [nsleep            = nsleep   ]" << '\n';
    cout << "   [ntimes            = ntimes   ]" << '\n';
    cout << "   [raninit           = tf       ]" << '\n';
    cout << "   [mb2               = tf       ]" << '\n';
    cout << "   [rbuffsize         = rbs      ]" << '\n';
    cout << "   [wbuffsize         = wbs      ]" << '\n';
    cout << "   [groupsets         = groupsets]" << '\n';
    cout << "   [setbuf            = setbuf   ]" << '\n';
    cout << "   [nfileitertest     = tf       ]" << '\n';
    cout << "   [filetests         = tf       ]" << '\n';
    cout << "   [dirtests          = tf       ]" << '\n';
    cout << "   [testwritenfiles   = versions ]" << '\n';
    cout << "   [testreadmf        = tf       ]" << '\n';
    cout << "   [readFANames       = fanames  ]" << '\n';
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
  bool groupSets(false), setBuf(true);
  bool nfileitertest(false), filetests(false), dirtests(false);
  bool testreadmf(false);
  Array<int> testWriteNFilesVersions;
  Array<std::string> readFANames;


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

  int nWNFTests(pp.countval("testwritenfiles"));
  if(nWNFTests > 0) {
    pp.getarr("testwritenfiles", testWriteNFilesVersions, 0, nWNFTests);
  }

  pp.query("groupsets", groupSets);
  pp.query("setbuf", setBuf);

  pp.query("rbuffsize", rbs);
  pp.query("wbuffsize", wbs);
  RealDescriptor::SetReadBufferSize(rbs);
  RealDescriptor::SetWriteBufferSize(wbs);

  pp.query("nfileitertest", nfileitertest);
  pp.query("filetests", filetests);
  pp.query("dirtests", dirtests);
  pp.query("testreadmf", testreadmf);
  int nNames(pp.countval("readfanames"));
  if(nNames > 0) {
    pp.getarr("readfanames", readFANames, 0, nNames);
  }


  if(ParallelDescriptor::IOProcessor()) {
    cout << '\n';
    cout << "**************************************************" << '\n';
    cout << "nprocs            = " << nprocs << '\n';
    cout << "nfiles            = " << nfiles << '\n';
    cout << "maxgrid           = " << maxgrid << '\n';
    cout << "ncomps            = " << ncomps << '\n';
    cout << "nboxes            = " << nboxes << '\n';
    cout << "ntimes            = " << ntimes << '\n';
    cout << "raninit           = " << raninit << '\n';
    cout << "mb2               = " << mb2 << '\n';
    cout << "rbuffsize         = " << rbs << '\n';
    cout << "wbuffsize         = " << wbs << '\n';
    cout << "groupsets         = " << groupSets << '\n';
    cout << "setbuf            = " << setBuf << '\n';
    cout << "nfileitertest     = " << nfileitertest << '\n';
    cout << "filetests         = " << filetests << '\n';
    cout << "dirtests          = " << dirtests << '\n';
    cout << "testreadmf        = " << testreadmf << '\n';
    for(int i(0); i < testWriteNFilesVersions.size(); ++i) {
      cout << "testWriteNFilesVersions[" << i << "]    = " << testWriteNFilesVersions[i] << '\n';
    }
    for(int i(0); i < readFANames.size(); ++i) {
      cout << "readFANames[" << i << "]    = " << readFANames[i] << '\n';
    }

    cout << '\n';
    cout << "sizeof(int) = " << sizeof(int) << '\n';
    cout << "sizeof(size_t) = " << sizeof(size_t) << '\n';
    cout << "sizeof(long) = " << sizeof(long) << '\n';
    cout << "sizeof(long long) = " << sizeof(long long) << '\n';
    cout << "sizeof(std::streampos) = " << sizeof(std::streampos) << '\n';
    cout << "sizeof(std::streamoff) = " << sizeof(std::streamoff) << '\n';
    cout << "sizeof(std::streamsize) = " << sizeof(std::streamsize) << '\n';
    cout << '\n';
    cout << "std::numeric_limits<int>::min()  = " << std::numeric_limits<int>::min() << '\n';
    cout << "std::numeric_limits<int>::max()  = " << std::numeric_limits<int>::max() << '\n';
    cout << "std::numeric_limits<Real>::min() = " << std::numeric_limits<Real>::min() << '\n';
    cout << "std::numeric_limits<Real>::max() = " << std::numeric_limits<Real>::max() << '\n';
    cout << "***************************************************" << '\n';
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

  if(nfileitertest) {
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
  }



  if(filetests) {
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
  }




  if(dirtests) {
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
  }



  for(int v(0); v < testWriteNFilesVersions.size(); ++v) {
    if(ParallelDescriptor::IOProcessor()) {
      cout << "testWriteNFilesVersions[" << v << "] = " << testWriteNFilesVersions[v] << std::endl;
    }
    VisMF::Header::Version hVersion;
    switch(testWriteNFilesVersions[v]) {
      case 1:
        hVersion = VisMF::Header::Version_v1;
      break;
      case 2:
        hVersion = VisMF::Header::NoFabHeader_v1;
      break;
      case 3:
        hVersion = VisMF::Header::NoFabHeaderMinMax_v1;
      break;
      case 4:
        hVersion = VisMF::Header::NoFabHeaderFAMinMax_v1;
      break;
      default:
        BoxLib::Abort("**** Error:  bad hVersion.");
      }

    for(int itimes(0); itimes < ntimes; ++itimes) {
      if(ParallelDescriptor::IOProcessor()) {
        cout << endl << "--------------------------------------------------" << endl;
        cout << "Testing NFiles Write:  version = " << hVersion << endl;
      }

      TestWriteNFiles(nfiles, maxgrid, ncomps, nboxes, raninit, mb2,
                      hVersion, groupSets, setBuf);

      if(ParallelDescriptor::IOProcessor()) {
        cout << "==================================================" << endl;
        cout << endl;
      }
    }
  }



  if(testreadmf) {
    for(int itimes(0); itimes < ntimes; ++itimes) {
      if(ParallelDescriptor::IOProcessor()) {
        cout << endl << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Testing MF Read" << endl;
      }

      for(int i(0); i < readFANames.size(); ++i) {
        TestReadMF(readFANames[i]);
      }

      if(ParallelDescriptor::IOProcessor()) {
        cout << "##################################################" << endl;
        cout << endl;
      }
    }
  }



  BoxLib::Finalize();
  return 0;
}
// -------------------------------------------------------------
// -------------------------------------------------------------
