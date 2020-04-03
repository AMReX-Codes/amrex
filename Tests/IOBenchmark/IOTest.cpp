// -------------------------------------------------------------
// IOTest.cpp
// -------------------------------------------------------------
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_NFiles.H>

#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cerrno>
#include <deque>

#include <unistd.h>
#include <string.h>
#include <sys/stat.h>

using std::cout;
using std::endl;
using std::ends;
using std::ofstream;
using std::streamoff;

using namespace amrex;

const int XDIR(0);
const int YDIR(1);
const int ZDIR(2);
Real bytesPerMB(1.0e+06);
const bool verboseDir(true);


// -------------------------------------------------------------
void DirectoryTests() {
    int ndirs(256), nlevels(4);

    if(ParallelDescriptor::IOProcessor()) { 
      errno = 0;
      mkdir("testdir", 0755);
      std::cout << "_here 0:  errno = " << strerror(errno) << std::endl;
      errno = 0;
      rmdir("testdir");
      std::cout << "_here 1:  errno = " << strerror(errno) << std::endl;
      errno = 0;
      mkdir("testnest/n0/n1", 0755);
      std::cout << "_here 2:  errno = " << strerror(errno) << std::endl;
      errno = 0;
    }

    BL_PROFILE_VAR("mkdirs", mkdirs);
    for(int i(0); i < ndirs; ++i) {
      std::stringstream dirname;
      dirname << "dir" << i;
      if(ParallelDescriptor::IOProcessor()) {
        if( ! amrex::UtilCreateDirectory(dirname.str(), 0755, verboseDir)) {
          amrex::CreateDirectoryFailed(dirname.str());
        }
        for(int level(0); level < nlevels; ++level) {
          std::stringstream dirname;
          dirname << "dir" << i << "/Level_" << level;
          if( ! amrex::UtilCreateDirectory(dirname.str(), 0755, verboseDir)) {
            amrex::CreateDirectoryFailed(dirname.str());
          }
        }
      }
    }
    ParallelDescriptor::Barrier("waitfordir");
    BL_PROFILE_VAR_STOP(mkdirs);

    BL_PROFILE_VAR("renamedirs", renamedirs);
    for(int i(0); i < ndirs; ++i) {
      if(ParallelDescriptor::IOProcessor()) {
        std::stringstream dirname;
        dirname << "dir" << i;
        std::string newdirname;
        newdirname = dirname.str() + ".old";
        std::rename(dirname.str().c_str(), newdirname.c_str());
      }
    }
    ParallelDescriptor::Barrier("renamedirs");
    BL_PROFILE_VAR_STOP(renamedirs);
}


// -------------------------------------------------------------
void NFileTests(int nOutFiles, const std::string &filePrefix) {
  int myProc(ParallelDescriptor::MyProc());
  Vector<int> data(32);

  for(int i(0); i < data.size(); ++i) {
    data[i] = (100 * myProc) + i;
  }

  bool groupSets(false), setBuf(true);
  for(NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf); nfi.ReadyToWrite(); ++nfi) {
    nfi.Stream().write((const char *) data.dataPtr(), data.size() * sizeof(int));
  }
}



// -------------------------------------------------------------
void FileTests() {
  Vector<int> myInts(4096 * 4096);
  for(int i(0); i < myInts.size(); ++i) {
    myInts[i] = i;
  }

  std::fstream myFile;

  BL_PROFILE_VAR("makeafile", makeafile);
  myFile.open("myFile", std::ios::out|std::ios::trunc|std::ios::binary);
  myFile.write((const char *) myInts.dataPtr(), myInts.size() * sizeof(int));
  myFile.close();
  BL_PROFILE_VAR_STOP(makeafile);

  BL_PROFILE_VAR_NS("seektests", seektests);
  myFile.open("myFile", std::ios::in|std::ios::binary);
  myFile.seekg(0, std::ios::end);
  myFile.seekg(0, std::ios::beg);
  for(int i(0); i < myInts.size()/10; ++i) {
    BL_PROFILE_VAR_START(seektests);
    myFile.seekg(1, std::ios::cur);
    BL_PROFILE_VAR_STOP(seektests);
  }
  myFile.close();


  std::string dirname("/home/vince/Development/BoxLib/Tests/IOBenchmark/a/b/c/d");
  if(ParallelDescriptor::IOProcessor()) {
    if( ! amrex::UtilCreateDirectory(dirname, 0755, verboseDir)) {
      amrex::CreateDirectoryFailed(dirname);
    }
  }
  std::string rdirname("relative/e/f/g");
  if(ParallelDescriptor::IOProcessor()) {
    if( ! amrex::UtilCreateDirectory(rdirname, 0755, verboseDir)) {
      amrex::CreateDirectoryFailed(rdirname);
    }
  }
  std::string nsdirname("noslash");
  if(ParallelDescriptor::IOProcessor()) {
    if( ! amrex::UtilCreateDirectory(nsdirname, 0755, verboseDir)) {
      amrex::CreateDirectoryFailed(nsdirname);
    }
  }

}


// -------------------------------------------------------------
BoxArray MakeBoxArray(int maxgrid,  int nboxes) {
#if (BL_SPACEDIM == 2)
  IntVect ivlo(0, 0);
  IntVect ivhi(maxgrid - 1, maxgrid - 1);
#else
  IntVect ivlo(0, 0, 0);
  IntVect ivhi(maxgrid - 1, maxgrid - 1, maxgrid - 1);
#endif
  int iSide(pow(static_cast<Real>(nboxes), 1.0/3.0));
  Box tempBox(ivlo, ivhi);
  BoxArray bArray(nboxes);
  int ix(0), iy(0), iz(0);
  for(int ibox(0); ibox < nboxes; ++ibox) {
    Box sBox(tempBox);
    sBox.shift(XDIR, ix * maxgrid);
    sBox.shift(YDIR, iy * maxgrid);
#if (BL_SPACEDIM == 3)
    sBox.shift(ZDIR, iz * maxgrid);
#endif
    bArray.set(ibox, sBox);
    ++ix;
    if(ix > iSide) {
      ix = 0;
      ++iy;
    }
    if(iy > iSide) {
      iy = 0;
      ++iz;
    }
  }
  return bArray;
}


// -------------------------------------------------------------
void TestWriteNFiles(int nfiles, int maxgrid, int ncomps, int nboxes,
                     bool raninit, bool mb2,
		     VisMF::Header::Version whichVersion,
		     bool groupSets, bool setBuf,
		     bool useDSS, int nMultiFabs,
		     bool checkmf, const std::string &dirName)
{
  VisMF::SetNOutFiles(nfiles);
  VisMF::SetGroupSets(groupSets);
  VisMF::SetSetBuf(setBuf);
  VisMF::SetUseDynamicSetSelection(useDSS);
  if(mb2) {
    bytesPerMB = pow(2.0, 20);
  }

  bool useDir( ! dirName.empty());
  Vector<std::string> pathNames(nMultiFabs);
  if(useDir) {
    // ---- make the directory and nMultiFabs Level_n directories
    for(int nmf(0); nmf < nMultiFabs; ++nmf) {
      std::stringstream path;
      path << dirName << "/Level_" << nmf << "/";
      pathNames[nmf] = path.str();
    }
    if(ParallelDescriptor::IOProcessor()) {
      const bool verboseDir(false);
      if( ! amrex::UtilCreateDirectory(dirName, 0755, verboseDir)) {
        amrex::CreateDirectoryFailed(dirName);
      }
      for(int nmf(0); nmf < nMultiFabs; ++nmf) {
        if( ! amrex::UtilCreateDirectory(pathNames[nmf], 0755, verboseDir)) {
          amrex::CreateDirectoryFailed(pathNames[nmf]);
        }
      }
    }
    ParallelDescriptor::Barrier("waitfordirName");
  }

  BoxArray bArray(MakeBoxArray(maxgrid, nboxes));
  if(ParallelDescriptor::IOProcessor()) {
    cout << "  Timings for writing to " << nfiles << " files with version:  "
         << whichVersion << endl;
  }

  std::string mfName;
  switch(whichVersion) {
    case VisMF::Header::Version_v1:
      mfName = "TestMF";
    break;
    case VisMF::Header::NoFabHeader_v1:
      mfName = "TestMFNoFabHeader";
    break;
    case VisMF::Header::NoFabHeaderMinMax_v1:
      mfName = "TestMFNoFabHeaderMinMax";
    break;
    case VisMF::Header::NoFabHeaderFAMinMax_v1:
      mfName = "TestMFNoFabHeaderFAMinMax";
    break;
    default:
      amrex::Abort("**** Error in TestWriteNFiles:  bad version.");
  }

  // ---- make the MultiFabs
  Vector<std::string> mfNames(nMultiFabs);
  Vector<MultiFab *> multifabs(nMultiFabs);
  DistributionMapping dmap{bArray};
  for(int nmf(0); nmf < nMultiFabs; ++nmf) {
    std::stringstream suffix;
    suffix << "_" << nmf;
    if(useDir) {
      mfNames[nmf] = pathNames[nmf] + mfName + suffix.str();
    } else {
      mfNames[nmf] = mfName + suffix.str();
    }
    VisMF::RemoveFiles(mfNames[nmf], false);  // ---- not verbose

    multifabs[nmf] = new MultiFab(bArray, dmap, ncomps, 0);

    for(MFIter mfiset(*(multifabs[nmf])); mfiset.isValid(); ++mfiset) {
      for(int invar(0); invar < ncomps; ++invar) {
        if(raninit) {
          Real *dp = (*multifabs[nmf])[mfiset].dataPtr(invar);
	  for(int i(0); i < (*multifabs[nmf])[mfiset].box().numPts(); ++i) {
	    dp[i] = amrex::Random() + (1.0 + static_cast<Real> (invar));
	  }
        } else {
          (*multifabs[nmf])[mfiset].setVal<RunOn::Host>((100.0 * mfiset.index()) + invar +
	                                (static_cast<Real> (nmf) / 100.0), invar);
        }
      }
    }
  }


  long totalBytesWritten(0);


  VisMF::Header::Version currentVersion(VisMF::GetHeaderVersion());
  VisMF::SetHeaderVersion(whichVersion);

  ParallelDescriptor::Barrier("TestWriteNFiles:BeforeWrite");
  double wallTimeStart(ParallelDescriptor::second());

  for(int nmf(0); nmf < nMultiFabs; ++nmf) {
    totalBytesWritten += VisMF::Write(*multifabs[nmf], mfNames[nmf]);
  }
  double wallTime(ParallelDescriptor::second() - wallTimeStart);

  ParallelDescriptor::Barrier("TestWriteNFiles:AfterWrite");

  double wallTimeMax(wallTime);
  double wallTimeMin(wallTime);

  ParallelDescriptor::ReduceLongSum(totalBytesWritten, ParallelDescriptor::IOProcessorNumber());
  ParallelDescriptor::ReduceRealMin(wallTimeMin, ParallelDescriptor::IOProcessorNumber());
  ParallelDescriptor::ReduceRealMax(wallTimeMax, ParallelDescriptor::IOProcessorNumber());
  Real megabytes((static_cast<Real> (totalBytesWritten)) / bytesPerMB);

  if(ParallelDescriptor::IOProcessor()) {
    cout << std::setprecision(5);
    cout << "------------------------------------------" << endl;
    cout << "  Total megabytes       = " << megabytes << endl;
    cout << "  Write:  Megabytes/sec = " << megabytes/wallTimeMax << endl;
    cout << "  Wall clock time       = " << wallTimeMax << " s." << endl;
    cout << "  Min wall clock time   = " << wallTimeMin << " s." << endl;
    cout << "  Max wall clock time   = " << wallTimeMax << " s." << endl;
    cout << "------------------------------------------" << endl;
  }

  for(int nmf(0); nmf < nMultiFabs; ++nmf) {
    delete multifabs[nmf];
  }


  if(checkmf) {
    ParallelDescriptor::Barrier("TestWriteNFiles:checkmf");
    wallTime = ParallelDescriptor::second();

    bool isOk(true);
    for(int nmf(0); nmf < nMultiFabs; ++nmf) {
      isOk &= VisMF::Check(mfNames[nmf]);
    }
    wallTimeMax = ParallelDescriptor::second() - wallTime;
    ParallelDescriptor::ReduceRealMax(wallTimeMax, ParallelDescriptor::IOProcessorNumber());
    if(ParallelDescriptor::IOProcessor()) {
      cout << std::setprecision(5);
      cout << "------------------------------------------" << endl;
      cout << "VisMF::Check():  time = " << wallTimeMax << " s." << endl;
      if(isOk) {
        cout << "VisMF::Check():  multifab is ok." << endl;
      } else {
        cout << "**** Error:  VisMF::Check():  multifab is not ok." << endl;
      }
      cout << "------------------------------------------" << endl;
    }
  }

  VisMF::SetHeaderVersion(currentVersion);  // ---- set back to previous version
}


// -------------------------------------------------------------
void TestReadMF(const std::string &mfName, bool useSyncReads,
                int nMultiFabs, const std::string &dirName)
{
  bool useDir( ! dirName.empty());
  Vector<std::string> pathNames(nMultiFabs);
  Vector<std::string> mfNames(nMultiFabs);
  Vector<MultiFab *> multifabs(nMultiFabs);

  for(int nmf(0); nmf < nMultiFabs; ++nmf) {
    std::stringstream suffix;
    suffix << "_" << nmf;
    if(useDir) {
      std::stringstream path;
      path << dirName << "/Level_" << nmf << "/";
      pathNames[nmf] = path.str();
      mfNames[nmf] = pathNames[nmf] + mfName + suffix.str();
    } else {
      mfNames[nmf] = mfName + suffix.str();
    }
    multifabs[nmf] = new MultiFab;
  }

  VisMF::SetUseSynchronousReads(useSyncReads);
  VisMF::CloseAllStreams(); 

  ParallelDescriptor::Barrier("TestReadMF:BeforeRead");
  double wallTimeStart(ParallelDescriptor::second());

  Vector<Vector<char> > faHeaders(nMultiFabs);
  for(int nmf(0); nmf < nMultiFabs; ++nmf) {
    std::string faHName(mfNames[nmf] + "_H");
    bool bExitOnError(false);  // ---- dont exit if this file does not exist
    ParallelDescriptor::ReadAndBcastFile(faHName, faHeaders[nmf], bExitOnError);
  }
  VisMF::Read(*multifabs[0], mfNames[0], faHeaders[0].dataPtr(), 0); 
  const BoxArray& ba = multifabs[0]->boxArray();
  const DistributionMapping& dm = multifabs[0]->DistributionMap();
  const int ncomps = multifabs[0]->nComp();
  const int ng = multifabs[0]->nGrow();
  for(int nmf(1); nmf < nMultiFabs; ++nmf) {
      multifabs[nmf]->define(ba,dm,ncomps,ng);
      VisMF::Read(*multifabs[nmf], mfNames[nmf], faHeaders[nmf].dataPtr(), nmf); 
  }

  double wallTime(ParallelDescriptor::second() - wallTimeStart);

  ParallelDescriptor::Barrier("TestReadMF:AfterRead");

  for(int nmf(0); nmf < nMultiFabs; ++nmf) {
    for(int i(0); i < multifabs[nmf]->nComp(); ++i) {
      Real mfMin = multifabs[nmf]->min(i);
      Real mfMax = multifabs[nmf]->max(i);
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "MMMMMMMM:  i mfMin mfMax = " << i << "  " << mfMin << "  " << mfMax << std::endl;
      }
    }
  }

  double wallTimeMax(wallTime);
  double wallTimeMin(wallTime);

  ParallelDescriptor::ReduceRealMin(wallTimeMin);
  ParallelDescriptor::ReduceRealMax(wallTimeMax);

  long totalNBytes(0);

  for(int nmf(0); nmf < nMultiFabs; ++nmf) {
    for(MFIter mfi(*multifabs[nmf]); mfi.isValid(); ++mfi) {
      totalNBytes += (*multifabs[nmf])[mfi].nBytes();
    }
    delete multifabs[nmf];
  }
  ParallelDescriptor::ReduceLongSum(totalNBytes);

  Real megabytes((static_cast<Real> (totalNBytes)) / bytesPerMB);

  if(ParallelDescriptor::IOProcessor()) {
    cout << std::setprecision(5);
    cout << "------------------------------------------" << endl;
    cout << "  Total megabytes = " << megabytes << endl;
    cout << "  Read:  Megabytes/sec   = " << megabytes/wallTimeMax << endl;
    cout << "  Wall clock time = " << wallTimeMax << endl;
    cout << "  Min wall clock time = " << wallTimeMin << endl;
    cout << "  Max wall clock time = " << wallTimeMax << endl;
    cout << "------------------------------------------" << endl;
  }
}



// -------------------------------------------------------------
void DSSNFileTests(int noutfiles, const std::string &filePrefixIn,
                   bool useIter)
{
#ifdef BL_USE_MPI
  bool groupSets(false), setBuf(true);
  std::string filePrefix(filePrefixIn);

  if(useIter) {
    int myProc(ParallelDescriptor::MyProc());
    Vector<int> data(10240);

    for(int i(0); i < data.size(); ++i) {
      data[i] = (100 * myProc) + i;
    }

    NFilesIter nfi(noutfiles, filePrefix, groupSets, setBuf);
    nfi.SetDynamic(-1);
    for( ; nfi.ReadyToWrite(); ++nfi) {
      nfi.Stream().write((const char *) data.dataPtr(), data.size() * sizeof(int));
    }
  }

  filePrefix += "_Check";

  int myProc(ParallelDescriptor::MyProc());
  int nProcs    = ParallelDescriptor::NProcs();
  int nOutFiles = NFilesIter::ActualNFiles(noutfiles);
  int mySetPosition = NFilesIter::WhichSetPosition(myProc, nProcs, nOutFiles, groupSets);
  Vector<int> data(10240);
  int deciderProc(nProcs - 1), coordinatorProc(-1);
  int deciderTag(ParallelDescriptor::SeqNum());
  int coordinatorTag(ParallelDescriptor::SeqNum());
  int doneTag(ParallelDescriptor::SeqNum());
  int writeTag(ParallelDescriptor::SeqNum());
  bool finishedWriting(false);
  ParallelDescriptor::Message rmess;
  int remainingWriters(nProcs);

  for(int i(0); i < data.size(); ++i) {
    data[i] = (100 * myProc) + i;
  }

  NFilesIter::CheckNFiles(nProcs, nOutFiles, false);

  int nSetZeros(0), nonZeroDeciderProc(-1);
  for(int i(0); i < nProcs; ++i) {
    // ---- count zero set positions  and find an alternate decider
    if(NFilesIter::WhichSetPosition(i, nProcs, nOutFiles, groupSets) == 0) {
      ++nSetZeros;
    } else {
      nonZeroDeciderProc = i;  // ---- this will end up with the last value
    }
  }

  if(NFilesIter::WhichSetPosition(deciderProc, nProcs, nOutFiles, groupSets) == 0) {
    deciderProc = nonZeroDeciderProc;
  }


    if(mySetPosition == 0) {    // ---- write data
      int fileNumber(NFilesIter::FileNumber(nOutFiles, myProc, groupSets));
      std::ofstream csFile;
      std::string FullName(amrex::Concatenate(filePrefix, fileNumber, 5));
      csFile.open(FullName.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
      if( ! csFile.good()) { amrex::FileOpenFailed(FullName); }
      // ----------------------------- write to file here
      csFile.write((const char *) data.dataPtr(), data.size() * sizeof(int));
      // ----------------------------- end write to file here
      csFile.flush();
      csFile.close();
      finishedWriting = true;

      // ---- tell the decider we are done
      ParallelDescriptor::Send(&myProc, 1, deciderProc, deciderTag);

      // ---- wait to find out who will coordinate
      ParallelDescriptor::Recv(&coordinatorProc, 1, deciderProc, coordinatorTag);

      if(myProc == coordinatorProc) {
	Vector<std::deque<int> > procsToWrite(nOutFiles);    // ---- [fileNumber](procsToWriteToFileNumber)
	// ---- populate with the static nfiles sets
	for(int i(0); i < nProcs; ++i) {
          int fileNumber(NFilesIter::FileNumber(nOutFiles, i, groupSets));
          int procSet(NFilesIter::WhichSetPosition(i, nProcs, nOutFiles, groupSets));
	  if(procSet == 0) {    // ---- set 0 procs have already written their data
	    --remainingWriters;
	  }
	  if(procSet != 0) {
	    procsToWrite[fileNumber].push_back(i);
	  }
	}

        // ---- signal each remaining processor when to write and to which file
	std::set<int> availableFileNumbers;
	availableFileNumbers.insert(fileNumber);  // ---- the coordinators file number

	// ---- recv incoming available files
	while(remainingWriters > 0) {

	  int nextProcToWrite, nextFileNumberToWrite, nextFileNumberAvailable;
	  std::set<int>::iterator ait = availableFileNumbers.begin();
	  nextFileNumberToWrite = *ait;
	  availableFileNumbers.erase(nextFileNumberToWrite);

	  for(int nfn(0); nfn < procsToWrite.size(); ++nfn) {
	    // ---- start with the current next file number
	    // ---- get a proc from another file number if the queue is empty
	    int tempNFN((nextFileNumberToWrite + nfn) % procsToWrite.size());
	    if(procsToWrite[tempNFN].size() > 0) {
	      nextProcToWrite = procsToWrite[tempNFN].front();
	      procsToWrite[tempNFN].pop_front();
	      break;  // ---- found one
	    }
	  }

          ParallelDescriptor::Asend(&nextFileNumberToWrite, 1, nextProcToWrite, writeTag);

          ParallelDescriptor::Recv(&nextFileNumberAvailable, 1, MPI_ANY_SOURCE, doneTag);
	  availableFileNumbers.insert(nextFileNumberAvailable);
	  --remainingWriters;
	}

      } else {
        // ---- tell the coordinatorProc we are done writing
        ParallelDescriptor::Send(&fileNumber, 1, coordinatorProc, doneTag);
      }

    } else if(myProc == deciderProc) {  // ---- this proc decides who decides

      // ---- the first message received is the coordinator
      ParallelDescriptor::Recv(&coordinatorProc, 1, MPI_ANY_SOURCE, deciderTag);
      // ---- tell the coordinatorProc to start coordinating
      ParallelDescriptor::Asend(&coordinatorProc, 1, coordinatorProc, coordinatorTag);
      for(int i(0); i < nSetZeros - 1; ++i) {  // ---- tell the others who is coorinating
        int nonCoordinatorProc(-1);
        ParallelDescriptor::Recv(&nonCoordinatorProc, 1, MPI_ANY_SOURCE, deciderTag);
        ParallelDescriptor::Asend(&coordinatorProc, 1, nonCoordinatorProc, coordinatorTag);
      }
    }

    // ---- these are the rest of the procs who need to write
    if( ! finishedWriting) {  // ---- the deciderProc drops through to here
      int fileNumber;
      // ---- wait for signal to start writing
      rmess = ParallelDescriptor::Recv(&fileNumber, 1, MPI_ANY_SOURCE, writeTag);
      coordinatorProc = rmess.pid();
      std::string FullName(amrex::Concatenate(filePrefix, fileNumber, 5));

      std::ofstream csFile;
      csFile.open(FullName.c_str(), std::ios::out | std::ios::app | std::ios::binary);
      csFile.seekp(0, std::ios::end);   // set to eof
      if( ! csFile.good()) { amrex::FileOpenFailed(FullName); }
      // ----------------------------- write to file here
      csFile.write((const char *) data.dataPtr(), data.size() * sizeof(int));
      // ----------------------------- end write to file here
      csFile.flush();
      csFile.close();
      finishedWriting = true;

      // ---- signal we are finished
      ParallelDescriptor::Send(&fileNumber, 1, coordinatorProc, doneTag);
    }
ParallelDescriptor::Barrier();
#endif
}


// -------------------------------------------------------------
// -------------------------------------------------------------


