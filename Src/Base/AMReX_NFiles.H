
#ifndef BL_NFILES_H
#define BL_NFILES_H

#include <string>
#include <fstream>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>

namespace amrex {


/**
* \brief This class encapsulates writing to nfiles.
*
* an example:
*
* for(NFilesIter nfi(nOutFiles, filePrefix); nfi.ReadyToWrite(); ++nfi) {
*   // do your writing here
*   nfi.Stream().write((const char *) data.dataPtr(), nChars);
* }
*/

class NFilesIter
{
  public:


    /**
    * \brief the nfiles will be named "filePrefix" + "00000"
    * through "filePrefix" + "noutfiles-1"
    * the range [1, nProcs] for noutfiles is enforced
    *
    * constructor for writing with static set selection
    *
    * \param noutfiles
    * \param &fileprefix
    * \param groupsets
    * \param setBuf
    */
    NFilesIter(int noutfiles, const std::string &fileprefix,
               bool groupsets, bool setBuf);


    /**
    * \brief call this to use dynamic set selection
    * deciderproc defaults to nprocs - 1 if \< 0
    *
    * \param deciderproc
    */
    void SetDynamic(int deciderproc = -1);
    bool GetDynamic() { return ( ! useStaticSetSelection); }


    /**
    * \brief call this to use a file per process for a sparse
    * set of writers.  ranksToWrite.size() will set noutfiles
    * rank ranksToWrite[n] will write to file_rank
    *
    * \param &ranksToWrite
    */
    void SetSparseFPP(const Vector<int> &ranksToWrite);
    bool GetSparseFPP() const { return useSparseFPP; }


    /**
    * \brief constructor for reading
    *
    * \param &fileName
    * \param &readRanks
    * \param setBuf
    */
    NFilesIter(const std::string &fileName,
               const Vector<int> &readRanks,
               bool setBuf = false);

    ~NFilesIter();


    /**
    * \brief if appendFirst is true, the first set for this
    * iterator will open the files in append mode
    *
    * \param appendFirst
    */
    bool ReadyToWrite(bool appendFirst = false);

    bool ReadyToRead();

    NFilesIter &operator++();

    std::fstream &Stream() { return fileStream; }


    /**
    * \brief get the current Stream()'s seek position
    */
    std::streampos SeekPos();

    static int LengthOfSet(int nProcs, int nOutFiles) {
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
        whichset = whichproc % LengthOfSet(nprocs, noutfiles);
      }
      return whichset;
    }


    /**
    * \brief this returns the actual number of files used
    * the range [1, nProcs] is enforced
    *
    * \param nOutFiles
    */
    static int ActualNFiles(int nOutFiles)
    {
      return( std::max(1, std::min(ParallelDescriptor::NProcs(), nOutFiles)) );
    }


    /**
    * \brief this checks if nOutFiles equals the calculated number of files
    * returns false if they do not match
    *
    * \param nProcs
    * \param nOutFiles
    * \param groupSets
    */
    static bool CheckNFiles(int nProcs, int nOutFiles, bool groupSets);

    static int FileNumber(int nOutFiles, int whichProc, bool groupSets)
    {
      BL_ASSERT(whichProc >= 0 && whichProc < ParallelDescriptor::NProcs());

      int anf(ActualNFiles(nOutFiles));
      if(groupSets) {
        return(whichProc % anf);
      } else {

	    int nProcs(ParallelDescriptor::NProcs());
        return(whichProc / LengthOfSet(nProcs, anf));
      }
    }

    const std::string &FileName() const { return fullFileName; }
    int FileNumber() const { return fileNumber; }

    static std::string FileName(int nOutFiles,
                                const std::string &filePrefix,
				int whichProc, bool groupSets)
    {
      BL_ASSERT(whichProc >= 0 && whichProc < ParallelDescriptor::NProcs());

      return ( amrex::Concatenate(filePrefix,
                                   FileNumber(ActualNFiles(nOutFiles),
				              whichProc, groupSets), minDigits ) );
    }

    static std::string FileName(int fileNumber,
                                const std::string &filePrefix)
    {
      return ( amrex::Concatenate(filePrefix, fileNumber, minDigits ) );
    }


    /**
    * \brief this is the processor coordinating dynamic set selection
    */
    int CoordinatorProc() const { return coordinatorProc; }


    /**
    * \brief these are the file numbers to which each rank wrote  [rank]
    * a rank only writes to one file
    */
    Vector<int> FileNumbersWritten();


    /**
    * \brief these are the order of ranks which wrote to each file
    * [filenumber][ranks in order they wrote to filenumber]
    */
    const Vector< Vector<int> > &FileNumbersWriteOrder() const { return fileNumbersWriteOrder; }

    void CleanUpMessages();

    static int  GetMinDigits()       { return minDigits; }

    static void SetMinDigits(int md) { minDigits = md;   }

  private:

    int myProc;
    int nProcs;
    int nOutFiles;
    int nSets;
    int mySetPosition;
    std::string fullFileName;
    std::string filePrefix;
    VisMF::IO_Buffer io_buffer;
    std::fstream fileStream;
    bool finishedWriting;
    bool groupSets;
    bool isReading;
    bool finishedReading;
    Vector<int> readRanks;
    Vector< Vector<int> > fileNumbersWriteOrder;        //!< [filenumber][ranks in order they wrote to filenumber]
    int myReadIndex;
    bool useStaticSetSelection;
    int fileNumber, remainingWriters, deciderProc, coordinatorProc;
    int deciderTag, coordinatorTag, doneTag, writeTag, stReadTag, stWriteTag;
    Vector<int> availableDeciders;
    Vector<int> setZeroProcs;
    bool useSparseFPP;
    Vector<int> sparseWritingRanks;
    int mySparseFileNumber;

    static int currentDeciderIndex;

    //! these were ignored by the decider procs and need to be cleaned up
    Vector<std::pair<int, int> > unreadMessages;    //!< [](tag, nmessages)

    static const int indexUndefined = -1;

    static int minDigits;        //!< for Concatenate

    NFilesIter();  //!< disallow
};

}

#endif  /* BL_NFILES_H */
