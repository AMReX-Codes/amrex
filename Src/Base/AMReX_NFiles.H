
#ifndef BL_NFILES_H
#define BL_NFILES_H
#include <AMReX_Config.H>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMFBuffer.H>

#include <fstream>
#include <string>

/**
 * \file AMReX_NFiles.H
 * \brief Helpers for writing/reading datasets across N output files.
 */

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
     * \brief Construct a writer iterator producing files `filePrefix00000`, …, `filePrefixNNNNN`.
     *
     * The requested file count is clamped to `[1, NProcs]`.
     *
     * \param noutfiles Desired number of output files.
     * \param fileprefix Prefix applied to every file name.
     * \param groupsets  If true, ranks are partitioned into contiguous groups; otherwise they interleave.
     * \param setBuf     If true, enable buffered I/O on the file stream.
     */
    NFilesIter(int noutfiles, std::string fileprefix,
               bool groupsets, bool setBuf);


    /**
     * \brief Enable dynamic set selection coordinated by \p deciderproc.
     *
     * The coordinator defaults to `NProcs-1` if \p deciderproc < 0.
     *
     * \param deciderproc Rank that decides which writers proceed.
     */
    void SetDynamic(int deciderproc = -1);
    [[nodiscard]] bool GetDynamic() const { return ( ! useStaticSetSelection); }


    /**
     * \brief Configure a sparse file-per-process mode.
     *
     * Each rank in \p ranksToWrite writes to its own file.
     *
     * \param ranksToWrite Ordered list of participating writer ranks.
     */
    void SetSparseFPP(const Vector<int> &ranksToWrite);
    [[nodiscard]] bool GetSparseFPP() const { return useSparseFPP; }


    /**
     * \brief Construct a reader for an existing multi-file dataset.
     *
     * \param fileName  Base filename/prefix.
     * \param readRanks Ranks that will read the data.
     * \param setBuf    If true, enable buffered I/O.
     */
    NFilesIter(std::string fileName,
               Vector<int> readRanks,
               bool setBuf = false);

    ~NFilesIter();

    NFilesIter (NFilesIter const&) = delete;
    NFilesIter (NFilesIter &&) = delete;
    NFilesIter& operator= (NFilesIter const&) = delete;
    NFilesIter& operator= (NFilesIter &&) = delete;

    /**
     * \brief Determine whether the current rank should write during this iteration.
     *
     * \param appendFirst Open the file in append mode for the first write set.
     * \return True if the caller may write; false otherwise.
     */
    bool ReadyToWrite(bool appendFirst = false);

    /**
     * \brief Determine whether the current rank should read during this iteration.
     * \return True if the caller may issue read operations; false if this rank is idle.
     */
    bool ReadyToRead();

    /**
     * \brief Advance to the next writer/reader set (pre-increment).
     * \return Reference to this iterator after advancement.
     */
    NFilesIter &operator++();

    //! Return the fstream associated with the current file set.
    std::fstream &Stream() { return fileStream; }


    /**
     * \brief Return the current seek position of Stream().
     *
     * \return Offset within the underlying file.
     */
    std::streampos SeekPos();

    /**
     * \brief Compute the maximum number of entries per set.
     *
     * \param nProcs    Total number of ranks.
     * \param nOutFiles Requested number of output files.
     * \return Size of each set such that every rank appears in exactly one set.
     */
    static int LengthOfSet(int nProcs, int nOutFiles) {
      int anf(ActualNFiles(nOutFiles));
      return ((nProcs + (anf - 1)) / anf);
    }

    /**
     * \brief Determine the set index occupied by a rank.
     *
     * \param whichproc Rank index to query.
     * \param nprocs    Total number of ranks.
     * \param noutfiles Requested number of output files.
     * \param groupsets If true, sets are contiguous; otherwise they interleave.
     * \return Zero-based set index holding \p whichproc.
     */
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
     * \brief Clamp the requested file count to `[1, NProcs]`.
     *
     * \param nOutFiles Requested number of output files.
     * \return Feasible number of files that will actually be used.
     */
    static int ActualNFiles(int nOutFiles)
    {
      return( std::max(1, std::min(ParallelDescriptor::NProcs(), nOutFiles)) );
    }


    /**
     * \brief Validate that \p nOutFiles matches the derived number of files for \p nProcs.
     *
     * \param nProcs    Number of MPI ranks participating.
     * \param nOutFiles Requested number of files.
     * \param groupSets Whether the writers are grouped contiguously.
     * \return True if the request is consistent, false otherwise.
     */
    static bool CheckNFiles(int nProcs, int nOutFiles, bool groupSets);

    /**
     * \brief Map a rank to its assigned file number.
     *
     * \param nOutFiles Requested number of files (will be clamped internally).
     * \param whichProc Rank id whose file number is queried.
     * \param groupSets Whether ranks are grouped contiguously when mapping to files.
     * \return Zero-based file number that \p whichProc uses.
     */
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

    [[nodiscard]] const std::string &FileName() const { return fullFileName; }
    [[nodiscard]] int FileNumber() const { return fileNumber; }

    /**
     * \brief Return the file name for the rank \p whichProc under the given policy.
     *
     * \param nOutFiles   Requested number of files.
     * \param filePrefix  File name prefix.
     * \param whichProc   Rank whose file name to build.
     * \param groupSets   Whether ranks are grouped contiguously.
     * \return Full file name including numeric suffix.
     */
    static std::string FileName(int nOutFiles,
                                const std::string &filePrefix,
                                int whichProc, bool groupSets)
    {
      BL_ASSERT(whichProc >= 0 && whichProc < ParallelDescriptor::NProcs());

      return ( amrex::Concatenate(filePrefix,
                                   FileNumber(ActualNFiles(nOutFiles),
                                              whichProc, groupSets), minDigits ) );
    }

    /**
     * \brief Return the file name corresponding to \p fileNumber.
     *
     * \param fileNumber Zero-based file index.
     * \param filePrefix File name prefix.
     * \return Full file name including numeric suffix.
     */
    static std::string FileName(int fileNumber,
                                const std::string &filePrefix)
    {
      return ( amrex::Concatenate(filePrefix, fileNumber, minDigits ) );
    }


    //! Return the rank coordinating dynamic set selection.
    [[nodiscard]] int CoordinatorProc() const { return coordinatorProc; }


    /**
     * \brief Return, for each rank, the file number it wrote to (at most one per rank).
     *
     * \return Vector indexed by rank id giving the file number used, or -1 if the rank never wrote.
     */
    Vector<int> FileNumbersWritten();


    /**
     * \brief For each file, return the ranks that wrote to it in chronological order.
     *
     * Indexed as `[file_number][writer_index]`.
     */
    [[nodiscard]] const Vector< Vector<int> > &FileNumbersWriteOrder() const { return fileNumbersWriteOrder; }

    //! Release any outstanding MPI control messages (dynamic mode only).
    void CleanUpMessages();

    //! Return the number of digits used when formatting file numbers.
    static int  GetMinDigits()       { return minDigits; }

    //! Override the number of digits used when formatting file numbers.
    static void SetMinDigits(int md) { minDigits = md;   }

  private:

    int myProc = -1;
    int nProcs = -1;
    int nOutFiles = -1;
    int nSets = -1;
    bool groupSets = false;
    int mySetPosition = -1;
    int fileNumber = -1;
    std::string filePrefix;
    std::string fullFileName;
    VisMFBuffer::IO_Buffer io_buffer;
    std::fstream fileStream;
    bool finishedWriting = false;
    bool isReading = false;
    bool finishedReading = false;
    Vector<int> readRanks;
    Vector< Vector<int> > fileNumbersWriteOrder;        //!< [filenumber][ranks in order they wrote to filenumber]
    int myReadIndex = -1;
    bool useStaticSetSelection = true;
    int remainingWriters = -1, deciderProc=-1, coordinatorProc = -1;
    int deciderTag = -1, coordinatorTag = -1, doneTag = -1, writeTag=-1;
    int stWriteTag = -1;
    int stReadTag = -1;
    Vector<int> availableDeciders;
    Vector<int> setZeroProcs;
    bool useSparseFPP = false;
    Vector<int> sparseWritingRanks;
    int mySparseFileNumber = -1;

    static int currentDeciderIndex;

    //! these were ignored by the decider procs and need to be cleaned up
    Vector<std::pair<int, int> > unreadMessages;    //!< [](tag, nmessages)

    static const int indexUndefined = -1;

    static AMREX_EXPORT int minDigits;        //!< for Concatenate
};

}

#endif  /* BL_NFILES_H */
