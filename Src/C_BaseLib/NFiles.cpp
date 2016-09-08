
#include <Utility.H>
#include <NFiles.H>


NFilesIter::NFilesIter(int noutfiles, const std::string &fileprefix,
                       bool groupsets, bool setBuf)
{
  isReading = false;
  groupSets = groupsets;
  myProc    = ParallelDescriptor::MyProc();
  nProcs    = ParallelDescriptor::NProcs();
  nOutFiles = ActualNFiles(noutfiles);
  nSets     = NSets(nProcs, nOutFiles);
  if(groupSets) {
    mySet     = myProc / nOutFiles;
  } else {
    mySet     = myProc % nSets;
  }
  int fileNumber(FileNumber(nOutFiles, myProc, groupSets));
  filePrefix    = fileprefix;
  fullFileName  = BoxLib::Concatenate(filePrefix, fileNumber, 5);

  finishedWriting = false;

  if(setBuf) {
    io_buffer.resize(VisMF::GetIOBufferSize());
    fileStream.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
  }

  bool checkNFiles(true);
  if(checkNFiles) {
    CheckNFiles(nProcs, nOutFiles, groupSets);
  }

}


NFilesIter::NFilesIter(const std::string &filename,
		       const Array<int> &readranks,
                       bool setBuf)
{
  isReading = true;
  myProc    = ParallelDescriptor::MyProc();
  nProcs    = ParallelDescriptor::NProcs();
  fullFileName  = filename;
  readRanks = readranks;
  myReadIndex = indexUndefined;
  for(int i(0); i < readRanks.size(); ++i) {
    if(myProc == readRanks[i]) {
      if(myReadIndex != indexUndefined) {
        BoxLib::Abort("**** Error in NFilesIter:  readRanks not unique.");
      }
      myReadIndex = i;
    }
  }

  if(myReadIndex == indexUndefined) {  // ---- nothing to read
    finishedReading = true;
    return;
  } else {
    finishedReading = false;
  }

  if(setBuf) {
    io_buffer.resize(VisMF::GetIOBufferSize());
    fileStream.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
  }

}


NFilesIter::~NFilesIter() {
}


bool NFilesIter::ReadyToWrite() {

  if(finishedWriting) {
    return false;
  }

  for(int iSet(0); iSet < nSets; ++iSet) {
    if(mySet == iSet) {
      if(iSet == 0) {   // ---- first set
	BoxLib::USleep(myProc);
	std::cout << myProc << ":: writing to:  " << fullFileName << std::endl;
        fileStream.open(fullFileName.c_str(),
                        std::ios::out | std::ios::trunc | std::ios::binary);
      } else {
	BoxLib::USleep(myProc);
	std::cout << myProc << ":: appending to:  " << fullFileName << std::endl;
        fileStream.open(fullFileName.c_str(),
                        std::ios::out | std::ios::app | std::ios::binary);
        fileStream.seekp(0, std::ios::end);   // ---- set to eof
      }
      if( ! fileStream.good()) {
        BoxLib::FileOpenFailed(fullFileName);
      }
      return true;
    }

    if(mySet == (iSet + 1)) {   // ---- next set waits
      int iBuff, waitForPID(-1), tag(-2);
      if(groupSets) {
        waitForPID = (myProc - nOutFiles);
        tag = (myProc % nOutFiles);
      } else {
        waitForPID = (myProc - 1);
        tag = (myProc / nSets);
      }
	BoxLib::USleep(myProc);
	std::cout << myProc << ":: waiting for:  " << waitForPID << std::endl;
      ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
    }
  }
  return false;
}


bool NFilesIter::ReadyToRead() {

  if(finishedReading) {
    return false;
  }

  if(myReadIndex != 0) {    // ---- wait for rank myReadIndex - 1
    int iBuff(-1), waitForPID(readRanks[myReadIndex - 1]);
    int tag(readRanks[0]);
    ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
  }

  fileStream.open(fullFileName.c_str(),
                  std::ios::in | std::ios::binary);
  if( ! fileStream.good()) {
    BoxLib::FileOpenFailed(fullFileName);
  }
  return true;
}


NFilesIter &NFilesIter::operator++() {
  if(isReading) {
    fileStream.close();

    if(myReadIndex < readRanks.size() - 1) {
      int iBuff(0), wakeUpPID(readRanks[myReadIndex + 1]);
      int tag(readRanks[0]);
      ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
    }
    finishedReading = true;

  } else {  // ---- writing

    fileStream.flush();
    fileStream.close();

    int iBuff(0), wakeUpPID(-1), tag(-2);
    if(groupSets) {
      wakeUpPID = (myProc + nOutFiles);
      tag = (myProc % nOutFiles);
    } else {
      wakeUpPID = (myProc + 1);
      tag = (myProc / nSets);
    }
    if(wakeUpPID < nProcs) {
      ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
    }
    finishedWriting = true;
  }

  return *this;
}


std::streampos NFilesIter::SeekPos() {
  return fileStream.tellp();
}


bool NFilesIter::CheckNFiles(int nProcs, int nOutFiles, bool groupSets)
{
  if(ParallelDescriptor::IOProcessor()) {
    std::set<int> fileNumbers;
    for(int i(0); i < nProcs; ++i) {
      fileNumbers.insert(FileNumber(nOutFiles, i, groupSets));
    }
    std::cout << "nOutFiles fileNumbers.size() = " << nOutFiles
              << "  " << fileNumbers.size() << std::endl;
    if(nOutFiles != fileNumbers.size()) {
      std::cout << "**** Different number of files." << std::endl;
      return false;
    }
  }
  return true;
}



