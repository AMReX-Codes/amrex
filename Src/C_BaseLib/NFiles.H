
#ifndef BL_NFILES_H
#define BL_NFILES_H

#include <string>
#include <fstream>

#include <ParallelDescriptor.H>
#include <VisMF.H>

//
// This class encapsulates writing to nfiles.
//
// an example:
// 
// for(NFilesIter nfi(nOutFiles, filePrefix); nfi.ReadyToWrite(); ++nfi) {
//   // ---- do your writing here
//   nfi.Stream().write((const char *) data.dataPtr(), nChars);
// }

//

class NFilesIter
{
  public:

    // ---- the nfiles will be named "filePrefix" + "00000"
    // ---- through "filePrefix" + "noutfiles-1"
    // ---- the range [1, nProcs] for noutfiles is enforced

    // ---- constructor for writing with static set selection
    NFilesIter(int noutfiles, const std::string &fileprefix,
               bool groupsets, bool setBuf);

    // ---- constructor for writing with dynamic set selection
    // ---- deciderproc defaults to nprocs - 1 if < 0
    void SetDynamic(int deciderproc = -1);

    // ---- constructor for reading
    NFilesIter(const std::string &fileName,
               const Array<int> &readRanks,
               bool setBuf = false);

    ~NFilesIter();

    bool ReadyToWrite();
    bool ReadyToRead();
    
    NFilesIter &operator++();

    std::fstream &Stream() { return fileStream; }

    // ---- get the current Stream()'s seek position
    std::streampos SeekPos();

    static int SetLength(int nProcs, int nOutFiles) {
      int anf(ActualNFiles(nOutFiles));
      return ((nProcs + (anf - 1)) / anf);
    }

    static int WhichSetPosition(int whichproc, int nprocs,
                                int noutfiles, bool groupsets)
    {
      int whichset;
      if(groupsets) {
        whichset = whichproc / noutfiles;
      } else {
        whichset = whichproc % SetLength(nprocs, noutfiles);
      }
      return whichset;
    }

    // ---- this returns the actual number of files used
    // ---- the range [1, nProcs] is enforced
    static int ActualNFiles(int nOutFiles)
    {
      return( std::max(1, std::min(ParallelDescriptor::NProcs(), nOutFiles)) );
    }

    // ---- this checks if nOutFiles equals the calculated number of files
    // ---- returns false if they do not match
    static bool CheckNFiles(int nProcs, int nOutFiles, bool groupSets);

    static int FileNumber(int nOutFiles, int whichProc, bool groupSets)
    {
      BL_ASSERT(whichProc >= 0 && whichProc < ParallelDescriptor::NProcs());
      int anf(ActualNFiles(nOutFiles));
      if(groupSets) {
        return(whichProc % anf);
      } else {
	int nProcs(ParallelDescriptor::NProcs());
        return(whichProc / SetLength(nProcs, anf));
      }
    }

    static std::string FileName(int nOutFiles,
                                const std::string &filePrefix,
				int whichProc, bool groupSets)
    {
      BL_ASSERT(whichProc >= 0 && whichProc < ParallelDescriptor::NProcs());
      return ( BoxLib::Concatenate(filePrefix,
                                   FileNumber(ActualNFiles(nOutFiles),
				              whichProc, groupSets), 5 ) );
    }

    static std::string FileName(int fileNumber,
                                const std::string &filePrefix)
    {
      return ( BoxLib::Concatenate(filePrefix, fileNumber, 5 ) );
    }

    // ---- this is the processor coorinating dynamic set selection
    int CoordinatorProc() const { return coordinatorProc; }

    // ---- these are the file numbers to which each rank wrote
    const Array<int> &FileNumbersWritten() const { return fileNumbersWritten; }


  private:

    int myProc;
    int nProcs;
    int nOutFiles;
    int nSets;
    int mySetPosition;
    std::string fullFileName;
    std::string filePrefix;
    std::fstream fileStream;
    VisMF::IO_Buffer io_buffer;
    bool finishedWriting;
    bool groupSets;
    bool isReading;
    bool finishedReading;
    Array<int> readRanks;
    Array<int> fileNumbersWritten;    // ---- [rank]
    int myReadIndex;
    bool useStaticSetSelection;
    int fileNumber, remainingWriters, deciderProc, coordinatorProc;
    int deciderTag, coordinatorTag, doneTag, writeTag;
    int nSetZeros, nonZeroDeciderProc;


    static const int indexUndefined = -1;

    NFilesIter();  // ---- disallow
};

#endif  /* BL_NFILES_H */
