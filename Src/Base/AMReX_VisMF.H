
#ifndef AMREX_VISMF_H_
#define AMREX_VISMF_H_
#include <AMReX_Config.H>

#include <AMReX_AsyncOut.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_FabArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabConv.H>
#include <AMReX_NFiles.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMFBuffer.H>

#include <fstream>
#include <iostream>
#include <sstream>
#include <deque>
#include <map>
#include <numeric>
#include <string>
#include <type_traits>

/**
 * \file AMReX_VisMF.H
 * \brief Visualization MultiFab (VisMF) on-disk I/O management.
 */

namespace amrex {

class IArrayBox;

/**
* \ingroup amrex_io
* \brief File I/O for FabArray<FArrayBox>.
*  Wrapper class for reading/writing FabArray<FArrayBox> objects to disk in various "smart" ways.
*/
class VisMF
    : public VisMFBuffer
{
public:
    /**
    * \brief How we write out FabArray<FArrayBox>s.
    * These are deprecated, we always use NFiles.
    * For OneFilePerCPU, set NFiles to NProcs.
    */
    enum How { OneFilePerCPU, NFiles };
    /**
    * \brief Construct by reading in the on-disk VisMF of the specified name.
    * The FABs in the on-disk FabArray are read on demand unless
    * the entire FabArray is requested. The name here is the name of
    * the FabArray not the name of the on-disk files.
    */
    explicit VisMF (std::string fafab_name);
    ~VisMF () = default;
    VisMF (const VisMF&) = delete;
    VisMF (VisMF&&) = delete;
    VisMF& operator= (const VisMF&) = delete;
    VisMF& operator= (VisMF&&) = delete;
    //! A structure containing info regarding an on-disk FAB.
    struct FabOnDisk
    {
        //! The default constructor -- null out all fields.
        FabOnDisk () = default;
        //! Constructor that sets the two values.
        FabOnDisk (std::string name, Long offset);
        //! The two data values in a FabOnDisk structure.
        std::string m_name; //!< The name of file containing the FAB.
        Long m_head = 0;     //!< Offset to start of FAB in file.
    };
    //! An on-disk FabArray<FArrayBox> contains this info in a header file.
    struct Header
    {
        //! The versions of the FabArray<FArrayBox> Header code.
        enum Version {
            Undefined_v1           = 0,  //!< ---- undefined
            Version_v1             = 1,  //!< ---- auto converting version with headers
                                         //!< ---- for each fab in the data files and
                                         //!< ---- min and max values for each fab in the header
            NoFabHeader_v1         = 2,  //!< ---- no fab headers, no fab mins or maxes
            NoFabHeaderMinMax_v1   = 3,  //!< ---- no fab headers,
                                         //!< ---- min and max values for each fab in the header
            NoFabHeaderFAMinMax_v1 = 4   //!< ---- no fab headers, no fab mins or maxes,
                                         //!< ---- min and max values for each FabArray in the header
        };
        //! The default constructor.
        Header ();
        //! Construct from a FabArray<FArrayBox>.
        Header (const FabArray<FArrayBox>& mf, VisMF::How how, Version version = Version_v1,
                bool calcMinMax = true, MPI_Comm = ParallelDescriptor::Communicator());

        ~Header () = default;
        Header (Header&& rhs) noexcept = default;
        Header (Header const&) = delete;
        Header& operator= (Header const&) = delete;
        Header& operator= (Header &&) = delete;

        //! Calculate the min and max arrays
        void CalculateMinMax(const FabArray<FArrayBox>& mf,
                             int procToWrite = ParallelDescriptor::IOProcessorNumber(),
                             MPI_Comm = ParallelDescriptor::Communicator());
        //
        // The data.
        //
        int                  m_vers{VisMF::Header::Undefined_v1};  //!< The version of the Header.
        How                  m_how;   //!< How the MF was written to disk.
        int                  m_ncomp; //!< Number of components in MF.
        IntVect              m_ngrow; //!< The number of ghost cells in MF.
        BoxArray             m_ba;    //!< The BoxArray of the MF.
        Vector< FabOnDisk >   m_fod;   //!< FabOnDisk info for contained FABs.
        //
        // These are not defined if they are not in the header
        //
        Vector< Vector<Real> > m_min;   //!< The min()s of each component of FABs.  [findex][comp]
        Vector< Vector<Real> > m_max;   //!< The max()s of each component of FABs.  [findex][comp]
        Vector<Real>          m_famin; //!< The min()s of each component of the FabArray.  [comp]
        Vector<Real>          m_famax; //!< The max()s of each component of the FabArray.  [comp]
        RealDescriptor       m_writtenRD;
    };

    //! This structure is used to store the read order for each FabArray file
    struct FabReadLink
    {
        int rankToRead{-1};
        int faIndex{-1};
        Long fileOffset{-1};
        Box box;

        FabReadLink() = default;
        FabReadLink(int ranktoread, int faindex, Long fileoffset, const Box &b);
    };

    //! This structure is used to store file ifstreams that remain open
    struct PersistentIFStream
    {
        std::ifstream   *pstr{nullptr};
        std::streampos   currentPosition{0};
        bool             isOpen{false};
        VisMF::IO_Buffer ioBuffer;

        PersistentIFStream () = default;
        ~PersistentIFStream ();
        PersistentIFStream (PersistentIFStream const&) = delete;
        PersistentIFStream (PersistentIFStream &&) = delete;
        PersistentIFStream& operator= (PersistentIFStream const&) = delete;
        PersistentIFStream& operator= (PersistentIFStream &&) = delete;
    };

    /**
     * \brief Open (or reuse) the persistent stream associated with \p fileName.
     * \param fileName Fully qualified path to the FAB data file.
     * \return Pointer to an active ifstream positioned at the start of the file.
     */
    static std::ifstream *OpenStream(const std::string &fileName);
    /**
     * \brief Flush and optionally close the persistent stream for \p fileName.
     * \param fileName Fully qualified path to the FAB data file.
     * \param forceClose Close even if the stream is marked persistent.
     */
    static void CloseStream(const std::string &fileName, bool forceClose = false);
    //! Remove the persistent stream associated with \p fileName from the cache.
    static void DeleteStream(const std::string &fileName);
    //! Close every stream currently cached by VisMF.
    static void CloseAllStreams();
    //! Return true if the given header encodes data without per-FAB headers.
    static bool NoFabHeader(const VisMF::Header &hdr);

    //! The number of components in the on-disk FabArray<FArrayBox>.
    [[nodiscard]] int nComp () const;
    //! The grow factor of the on-disk FabArray<FArrayBox>.
    [[nodiscard]] int nGrow () const;
    //! Return the full ghost-cell vector (per dimension) of the on-disk data.
    [[nodiscard]] IntVect nGrowVect () const;
    //! # of FABs in the VisMF. Equal to # of Boxes in the BoxArray.
    [[nodiscard]] int size () const;
    //! The BoxArray of the on-disk FabArray<FArrayBox>.
    [[nodiscard]] const BoxArray& boxArray () const;
    //! The min of the FAB (in valid region) at specified index and component.
    [[nodiscard]] Real min (int fabIndex, int nComp) const;
    //! The min of the FabArray (in valid region) at specified component.
    [[nodiscard]] Real min (int nComp) const;
    //! The max of the FAB (in valid region) at specified index and component.
    [[nodiscard]] Real max (int fabIndex, int nComp) const;
    //! The max of the FabArray (in valid region) at specified component.
    [[nodiscard]] Real max (int nComp) const;

    /**
    * \brief The FAB at the specified index and component.
    *         Reads it from disk if necessary.
    *         This reads only the specified component.
    */
    [[nodiscard]] const FArrayBox& GetFab (int fabIndex, int compIndex) const;
    //! Delete()s the FAB at the specified index and component.
    void clear (int fabIndex, int compIndex);
    //! Delete()s the FAB at the specified index (all components).
    void clear (int fabIndex);
    //! Delete()s all the FABs.
    void clear ();
    /**
    * \brief Write a FabArray<FArrayBox> to disk in a "smart" way.
    * Returns the total number of bytes written on this processor.
    * If set_ghost is true, sets the ghost cells in the FabArray<FArrayBox> to
    * one-half the average of the min and max over the valid region
    * of each contained FAB.
    */
    /**
     * \brief Write \p mf to disk using the VisMF format.
     *
     * \param mf        FabArray to serialize.
     * \param name      Logical dataset name (determines directory layout).
     * \param how       File layout strategy (deprecated; use appropriate NFiles).
     * \param set_ghost Whether to fill ghost cells with the average of min/max before writing.
     * \return Total number of bytes written by this rank.
     */
    static Long Write (const FabArray<FArrayBox> &mf,
                       const std::string& name,
                       VisMF::How         how = NFiles,
                       bool               set_ghost = false);

    /**
     * \brief Launch an asynchronous write of \p mf (lvalue reference) on the AsyncOut queue.
     * \param mf                FabArray to serialize (must remain valid until the write completes).
     * \param mf_name           Logical dataset name.
     * \param valid_cells_only  If true, only valid regions are written.
     */
    static void AsyncWrite (const FabArray<FArrayBox>& mf, const std::string& mf_name,
                            bool valid_cells_only = false);
    /**
     * \brief Launch an asynchronous write of \p mf (rvalue reference) by moving ownership.
     * \param mf                FabArray to serialize (ownership transferred to AsyncOut).
     * \param mf_name           Logical dataset name.
     * \param valid_cells_only  If true, only valid regions are written.
     */
    static void AsyncWrite (FabArray<FArrayBox>&& mf, const std::string& mf_name,
                            bool valid_cells_only = false);

    /**
    * \brief Write only the header-file corresponding to FabArray<FArrayBox> to
    * disk without the corresponding FAB data. This writes BoxArray information
    * (which might still be needed by data post-processing tools such as yt)
    * when the FAB data is not needed. Returns the total number of bytes written
    * on this processor.
    */
    static Long WriteOnlyHeader (const FabArray<FArrayBox> & mf,
                                 const std::string         & mf_name,
                                 VisMF::How                  how = NFiles);
    //! this will remove nfiles associated with name and the header
    static void RemoveFiles(const std::string &name, bool verbose = false);

    /**
    * \brief Read a FabArray<FArrayBox> from disk written using
    * VisMF::Write().  If the FabArray<FArrayBox> fafab has been
    * fully defined, the BoxArray on the disk must match the BoxArray
    * in fafab.  If it is constructed with the default constructor,
    * the BoxArray on the disk will be used and a new
    * DistributionMapping will be made.  A pre-read FabArray header
    * can be passed in to avoid a read and broadcast.
    */
    /**
     * \brief Read a FabArray from disk into \p mf.
     *
     * \param mf               Destination FabArray (must be default-constructed or match the stored BoxArray).
     * \param name             Logical dataset name.
     * \param faHeader         Optional pre-read header buffer (pass nullptr to read from disk).
     * \param coordinatorProc  Rank coordinating the read (usually IOProcessor).
     * \param allow_empty_mf   Set to non-zero to permit empty multifabs.
     */
    static void Read (FabArray<FArrayBox> &mf,
                      const std::string &name,
                      const char *faHeader = nullptr,
                      int coordinatorProc = ParallelDescriptor::IOProcessorNumber(),
                      int allow_empty_mf = 0);

    //! Does FabArray exist?
    //! Return true if a VisMF dataset named \p name exists on disk.
    static bool Exist (const std::string &name);

    //! Read only the header of a FabArray, header will be resized here.
    //! Read only the header for dataset \p fafabName into \p header (resized as needed).
    static void ReadFAHeader (const std::string &fafabName,
                              Vector<char> &header);

    //! Check if the multifab is ok, false is returned if not ok
    //! Return true if all files belonging to \p name pass basic consistency checks.
    static bool Check (const std::string &name);
    //! Return the current file offset of the provided stream (in bytes).
    static Long FileOffset (std::ostream& os);
    //! Read the entire fab (all components).
    FArrayBox* readFAB (int idx, const std::string& mf_name);
    //! Read the specified fab component.
    FArrayBox* readFAB (int idx, int icomp);

    //! Return the number of output files currently used for VisMF writes.
    static int  GetNOutFiles ();
    //! Set the number of output files used for VisMF writes.
    static void SetNOutFiles (int noutfiles, MPI_Comm comm = ParallelDescriptor::Communicator());

    //! Return the number of concurrently open input streams.
    static int  GetMFFileInStreams ()  { return nMFFileInStreams; }
    //! Set the number of concurrently open input streams.
    static void SetMFFileInStreams (int nstreams, MPI_Comm comm = ParallelDescriptor::Communicator());

    //! Return the verbosity level used by VisMF.
    static int GetVerbose ()       { return verbose; }
    //! Set the verbosity level used by VisMF.
    static void SetVerbose (int v) { verbose = v; }

    //! Return the VisMF header version used for future writes.
    static VisMF::Header::Version GetHeaderVersion () { return currentVersion; }
    //! Override the VisMF header version used for future writes.
    static void SetHeaderVersion (VisMF::Header::Version version)
                                                   { currentVersion = version; }

    //! Return whether VisMF groups FABs into sets when writing.
    static bool GetGroupSets () { return groupSets; }
    //! Enable or disable grouping FABs into sets.
    static void SetGroupSets (bool groupsets) { groupSets = groupsets; }

    //! Return whether VisMF sets custom IO buffers.
    static bool GetSetBuf () { return setBuf; }
    //! Enable or disable custom IO buffers.
    static void SetSetBuf (bool setbuf) { setBuf = setbuf; }

    //! Return whether all ranks funnel reads through a single reader.
    static bool GetUseSingleRead () { return useSingleRead; }
    //! Enable or disable single-reader mode.
    static void SetUseSingleRead (bool usesingleread) { useSingleRead = usesingleread; }

    //! Return whether all ranks funnel writes through a single writer.
    static bool GetUseSingleWrite () { return useSingleWrite; }
    //! Enable or disable single-writer mode.
    static void SetUseSingleWrite (bool usesinglewrite) { useSingleWrite = usesinglewrite; }

    //! Return whether VisMF verifies stream positions when reading/writing.
    static bool GetCheckFilePositions () { return checkFilePositions; }
    //! Enable or disable stream-position checks.
    static void SetCheckFilePositions (bool cfp) { checkFilePositions = cfp; }

    //! Return whether persistent input streams are used.
    static bool GetUsePersistentIFStreams () { return usePersistentIFStreams; }
    //! Enable or disable persistent input streams.
    static void SetUsePersistentIFStreams (bool usepifs) { usePersistentIFStreams = usepifs; }

    //! Return whether synchronous reads are forced.
    static bool GetUseSynchronousReads () { return useSynchronousReads; }
    //! Enable or disable synchronous reads.
    static void SetUseSynchronousReads (bool usepsr) { useSynchronousReads = usepsr; }

    //! Return whether dynamic set selection is used for NFiles.
    static bool GetUseDynamicSetSelection () { return useDynamicSetSelection; }
    //! Enable or disable dynamic set selection.
    static void SetUseDynamicSetSelection (bool usedss) { useDynamicSetSelection = usedss; }

    //! Return whether VisMF skips flushing after writes.
    static bool GetNoFlushAfterWrite () { return noFlushAfterWrite; }
    //! Enable or disable flush suppression after writes.
    static void SetNoFlushAfterWrite (bool nfaw) { noFlushAfterWrite = nfaw; }

    //! Return whether VisMF enforces a barrier after finishing each level.
    static bool GetBarrierAfterLevel () { return barrierAfterLevel; }
    //! Enable or disable enforcing a barrier after each level is written.
    static void SetBarrierAfterLevel (bool bal) { barrierAfterLevel = bal; }

    //! Return the directory containing \p filename.
    static std::string DirName (const std::string& filename);
    //! Return the filename without directory components.
    static std::string BaseName (const std::string& filename);

    //! Initialize global VisMF state (streams, buffers, etc.).
    static void Initialize ();
    //! Finalize global VisMF state and release cached resources.
    static void Finalize ();

private:
    static FabOnDisk Write (const FArrayBox&   fab,
                            const std::string& filename,
                            std::ostream&      os,
                            Long&              bytes);

    static Long WriteHeaderDoit (const std::string &mf_name,
                                 VisMF::Header const &hdr);

    static Long WriteHeader (const std::string &mf_name,
                             VisMF::Header     &hdr,
                             int procToWrite = ParallelDescriptor::IOProcessorNumber(),
                             MPI_Comm comm = ParallelDescriptor::Communicator());

    //! fileNumbers must be passed in for dynamic set selection [proc]
    static void FindOffsets (const FabArray<FArrayBox> &mf,
                             const std::string &filePrefix,
                             VisMF::Header &hdr,
                             VisMF::Header::Version whichVersion,
                             NFilesIter &nfi,
                             MPI_Comm comm = ParallelDescriptor::Communicator());
    /**
    * \brief Make a new FAB from a fab in a FabArray<FArrayBox> on disk.
    * The returned *FAB will have either one component filled from
    * fafab[fabIndex][whichComp] or fafab[fabIndex].nComp() components.
    * whichComp == -1 means reads the whole FAB.
    * Otherwise read just that component.
    */
    static FArrayBox *readFAB (int                idx,
                               const std::string &mf_name,
                               const Header      &hdr,
                               int                whichComp = -1);
    //! Read the whole FAB into fafab[fabIndex]
    static void readFAB (FabArray<FArrayBox> &mf,
                         int                idx,
                         const std::string &mf_name,
                         const Header&      hdr);

    static void AsyncWriteDoit (const FabArray<FArrayBox>& mf, const std::string& mf_name,
                                bool is_rvalue, bool valid_cells_only);

    //! Name of the FabArray<FArrayBox>.
    std::string m_fafabname;
    //! The VisMF header as read from disk.
    Header m_hdr;
    //! We manage the FABs individually.
    mutable Vector< Vector<FArrayBox*> > m_pa;
    /**
    * \brief Persistent streams.  These open on demand and should
    * be closed when not needed with CloseAllStreams.
    * ~VisMF also closes them.  [filename, pifs]
    */
    static AMREX_EXPORT std::map<std::string, VisMF::PersistentIFStream> persistentIFStreams;
    //! The number of files to write for a FabArray<FArrayBox>.
    static AMREX_EXPORT int nOutFiles;
    static AMREX_EXPORT int nMFFileInStreams;

    static AMREX_EXPORT int verbose;
    static AMREX_EXPORT VisMF::Header::Version currentVersion;
    static AMREX_EXPORT bool groupSets;
    static AMREX_EXPORT bool setBuf;
    static AMREX_EXPORT bool useSingleRead;
    static AMREX_EXPORT bool useSingleWrite;
    static AMREX_EXPORT bool checkFilePositions;
    static AMREX_EXPORT bool usePersistentIFStreams;
    static AMREX_EXPORT bool useSynchronousReads;
    static AMREX_EXPORT bool useDynamicSetSelection;
    static AMREX_EXPORT bool allowSparseWrites;
    static AMREX_EXPORT bool noFlushAfterWrite;
    static AMREX_EXPORT bool barrierAfterLevel;
};

//! Write a FabOnDisk to an ostream in ASCII.
std::ostream& operator<< (std::ostream& os, const VisMF::FabOnDisk& fod);
//! Read a FabOnDisk from an istream.
std::istream& operator>> (std::istream& is, VisMF::FabOnDisk& fod);
//! Write an Vector<FabOnDisk> to an ostream in ASCII.
std::ostream& operator<< (std::ostream& os, const Vector<VisMF::FabOnDisk>& fa);
//! Read an Vector<FabOnDisk> from an istream.
std::istream& operator>> (std::istream& is, Vector<VisMF::FabOnDisk>& fa);
//! Write a VisMF::Header to an ostream in ASCII.
std::ostream& operator<< (std::ostream& os, const VisMF::Header& hd);
//! Read a VisMF::Header from an istream.
std::istream& operator>> (std::istream& is, VisMF::Header& hd);

/**
 * \brief Write `iMultiFab/FabArray<IArrayBox>`
 *
 *  This writes an `iMultiFab/FabArray<IArrayBox>` to files on disk, including
 *  a clear text file NAME_H and binary files NAME_D_00000 etc.
 *
 * \param fa is the iMultiFab to be written.
 * \param name is the base name for the files.
 */
template <typename FAB>
// This function does work for MultiFab, but let's disable it to avoid confusion.
std::enable_if_t<std::is_same_v<FAB,IArrayBox>>
Write (const FabArray<FAB>& fa, const std::string& name)
{
    BL_PROFILE("Write(FabArray)");
    AMREX_ASSERT(name.back() != '/');

    auto data_descriptor = FAB::getDataDescriptor();
    int data_bytes = data_descriptor->numBytes();

    bool useSparseFPP = false;
    const Vector<int> &pmap = fa.DistributionMap().ProcessorMap();
    std::set<int> procsWithData;
    Vector<int> procsWithDataVector;
    for(int i : pmap) {
        procsWithData.insert(i);
    }
    const int nOutFiles = VisMF::GetNOutFiles();
    if (static_cast<int>(procsWithData.size()) < nOutFiles) {
        useSparseFPP = true;
        for (auto i : procsWithData) {
            procsWithDataVector.push_back(i);
        }
    }

    std::string filePrefix = name + "_D_";

    NFilesIter nfi(nOutFiles, filePrefix, VisMF::GetGroupSets(), VisMF::GetSetBuf());

    if (useSparseFPP) {
        nfi.SetSparseFPP(procsWithDataVector);
    } else {
        nfi.SetDynamic();
    }

    const auto &fio = FAB::getFABio();

    for ( ; nfi.ReadyToWrite(); ++nfi) {
        for(MFIter mfi(fa); mfi.isValid(); ++mfi) {
            const FAB &fab = fa[mfi];
            {
                std::stringstream hss;
                fio.write_header(hss, fab, fab.nComp());
                auto hLength = static_cast<std::streamoff>(hss.tellp());
                auto tstr = hss.str();
                nfi.Stream().write(tstr.c_str(), hLength);
                nfi.Stream().flush();
            }
            auto const* fabdata = fab.dataPtr();
#ifdef AMREX_USE_GPU
            std::unique_ptr<FAB> hostfab;
            if (fab.arena()->isManaged() || fab.arena()->isDevice()) {
                hostfab = std::make_unique<FAB>(fab.box(), fab.nComp(), The_Pinned_Arena());
                Gpu::dtoh_memcpy_async(hostfab->dataPtr(), fab.dataPtr(),
                                       fab.size()*sizeof(typename FAB::value_type));
                Gpu::streamSynchronize();
                fabdata = hostfab->dataPtr();
            }
#endif
            Long writeDataItems = fab.box().numPts() * fa.nComp();
            Long writeDataSize = writeDataItems * data_bytes;
            nfi.Stream().write((char *) fabdata, writeDataSize);
            nfi.Stream().flush();
        }
    }

    int coordinatorProc = ParallelDescriptor::IOProcessorNumber();
    if (nfi.GetDynamic()) {
        coordinatorProc = nfi.CoordinatorProc();
    }

    if (ParallelDescriptor::MyProc() == coordinatorProc) {
        std::string header_file_name = name + "_H";
        VisMFBuffer::IO_Buffer io_buffer(VisMFBuffer::GetIOBufferSize());
        std::ofstream ofs;
        ofs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        ofs.open(header_file_name.c_str(), std::ios::out | std::ios::trunc);
        if (!ofs.good()) {
            amrex::FileOpenFailed(header_file_name);
        }

        ofs << "amrex::FabArray<" << FAB::getClassName() << "> v1.0\n";
        ofs << fa.nComp() << '\n';
        ofs << fa.nGrowVect() << '\n';
        fa.boxArray().writeOn(ofs);
        ofs << '\n';

        const DistributionMapping& dm = fa.DistributionMap();
        int nfabs = fa.boxArray().size();
        int nFiles = NFilesIter::ActualNFiles(nOutFiles);
        int nprocs = ParallelDescriptor::NProcs();

        Vector<Long> fabBytes(nfabs, 0);
        std::map<int, Vector<int> > rankBoxOrder;
        for (int i = 0; i < nfabs; ++i) {
            std::stringstream hss;
            FAB tmp(fa.fabbox(i), fa.nComp(), false);
            fio.write_header(hss, tmp, tmp.nComp());
            // Size includes header and data
            fabBytes[i] = static_cast<std::streamoff>(hss.tellp()) + tmp.size() * data_bytes;
            rankBoxOrder[dm[i]].push_back(i);
        }

        Vector<int> fileNumbers;
        if (nfi.GetDynamic()) {
            fileNumbers = nfi.FileNumbersWritten();
        } else if (nfi.GetSparseFPP()) {
            // if sparse, write to (file number = rank)
            fileNumbers.resize(nprocs);
            std::iota(fileNumbers.begin(), fileNumbers.end(), 0);
        } else {
            fileNumbers.resize(nprocs);
            for (int i = 0; i < nprocs; ++i) {
                fileNumbers[i] = NFilesIter::FileNumber(nFiles, i, VisMF::GetGroupSets());
            }
        }

        Vector<VisMF::FabOnDisk> fod(nfabs);

        const Vector< Vector<int> > &fileNumbersWriteOrder = nfi.FileNumbersWriteOrder();
        for (auto const& rv : fileNumbersWriteOrder) {
            Long currentOffset = 0;
            for (auto rank : rv) {
                auto rbo_it = rankBoxOrder.find(rank);
                if (rbo_it != rankBoxOrder.end()) {
                    Vector<int> const& index = rbo_it->second;
                    int whichFileNumber = fileNumbers[rank];
                    std::string const& whichFileName =
                        VisMF::BaseName(NFilesIter::FileName(whichFileNumber, filePrefix));
                    for (int i : index) {
                        fod[i].m_name = whichFileName;
                        fod[i].m_head = currentOffset;
                        currentOffset += fabBytes[i];
                    }
                }
            }
        }
        ofs << fod;
    }
}

/// \cond DOXYGEN_IGNORE
namespace detail {
template <typename FAB>
void read_fab (FAB& fab, VisMF::FabOnDisk const& fod, std::string const& name)
{
    std::string fullname = VisMF::DirName(name);
    fullname += fod.m_name;
    VisMFBuffer::IO_Buffer io_buffer(VisMFBuffer::GetIOBufferSize());
    std::ifstream ifs;
    ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    ifs.open(fullname.c_str(), std::ios::in | std::ios::binary);
    if (!ifs.good()) {
        amrex::FileOpenFailed(fullname);
    }
    ifs.seekg(fod.m_head, std::ios::beg);
    fab.readFrom(ifs);
}
}
/// \endcond

/**
 * \brief Read `iMultiFab/FabArray<IArrayBox>`
 *
 * This reads an `iMultiFab/FabArray<IArrayBox>` from disk.  If it has been
 * fully defined, the BoxArray on the disk must match the BoxArray in the
 * given `iMultiFab/FabArray<IArrayBox>` object.  If it is only constructed
 * with the default constructor, the BoxArray on the disk will be used and a
 * new DistributionMapping will be made.  When this function is used to
 * restart a calculation from checkpoint files, one should use a fully
 * defined `iMultiFab/FabArray<IArrayBox>` except for the first one in a
 * series of iMultiFab/MultiFab objects that share the same
 * BoxArray/DistributionMapping.  This will ensure that they share the same
 * BoxArray/DistributionMapping after restart.
 *
 * \param fa is the iMultiFab.
 * \param name is the base name for the files.
 */
template <typename FAB>
// This function does work for MultiFab, but let's disable it to avoid confusion.
std::enable_if_t<std::is_same_v<FAB,IArrayBox>>
Read (FabArray<FAB>& fa, const std::string& name)
{
    BL_PROFILE("Read(FabArray)");
    AMREX_ASSERT(name.back() != '/');

    BoxArray ba;
    int ncomp;
    IntVect ngrow;
    Vector<VisMF::FabOnDisk> fod;
    {
        std::string header_file_name = name + "_H";
        Vector<char> header_file_chars;
        ParallelDescriptor::ReadAndBcastFile(header_file_name, header_file_chars);
        std::string header_file_string(header_file_chars.data());
        std::stringstream ifs(header_file_string, std::istringstream::in);

        std::string type, version;
        ifs >> type >> version;
        AMREX_ASSERT(type == "amrex::FabArray<amrex::IArrayBox>" ||
                     type == "amrex::FabArray<amrex::FArrayBox>");
        ifs >> ncomp;
        ifs >> ngrow;
        ba.readFrom(ifs);
        ifs >> fod;
    }

    if (fa.empty()) {
        fa.define(ba, DistributionMapping{ba}, ncomp, ngrow);
    } else {
        AMREX_ASSERT(amrex::match(ba, fa.boxArray()));
    }

#ifdef AMREX_USE_MPI
    const int nopensperfile = VisMF::GetMFFileInStreams(); // # of concurrent readers per file
    const int myproc = ParallelDescriptor::MyProc();
    const int coordproc = ParallelDescriptor::IOProcessorNumber();

    int nreqs = 0;
    int allreadsindex = 0;
    std::map<std::string, int> filenames; // <filename, allreadsindex>

    const int nboxes = fa.size();
    const auto& dm = fa.DistributionMap();
    for (int i = 0; i < nboxes; ++i) {
        if (myproc == dm[i]) {
            ++nreqs;
        }
        if (myproc == coordproc) {
            std::string const& fname = fod[i].m_name;
            auto r =filenames.insert(std::make_pair(fname, allreadsindex));
            if (r.second) {
                ++allreadsindex;
            }
        }
    }

    const int readtag = ParallelDescriptor::SeqNum();
    const int donetag = ParallelDescriptor::SeqNum();

    if (myproc == coordproc) {
        std::multiset<int> availablefiles;  // [whichFile]  supports multiple reads/file
        Vector<std::map<int,std::map<Long,int> > > allreads; // [file]<proc,<seek,index>>

        const auto nfiles = static_cast<int>(filenames.size());
        for (int i = 0; i < nfiles; ++i) {
            for (int j = 0; j < nopensperfile; ++j) {
                availablefiles.insert(i);
            }
        }
        allreads.resize(nfiles);
        for (int i = 0; i < nboxes; ++i) {
            const auto whichproc = dm[i];
            const auto iseekpos = fod[i].m_head;
            std::string const& fname = fod[i].m_name;
            auto filenamesiter = filenames.find(fname);
            if (filenamesiter != filenames.end()) {
                const int fi = filenamesiter->second;
                allreads[fi][whichproc].insert(std::make_pair(iseekpos,i));
            } else {
                amrex::Error("Error in amrex::Read: filename not found "+fname);
            }
        }

        int totalioreqs = nboxes;
        int reqspending = 0;
        int iopfileindex;
        std::deque<int> iopreads;
        std::set<int> busyprocs;
        while (totalioreqs > 0) {
            auto afilesiter = availablefiles.begin();
            while (afilesiter != availablefiles.end()) {
                const int fi = *afilesiter;
                if (allreads[fi].empty()) {
                    availablefiles.erase(fi);
                    afilesiter = availablefiles.begin();
                    continue;
                }
                auto whichread = allreads[fi].begin();
                for ( ; whichread != allreads[fi].end(); ++whichread) {
                    const int tryproc = whichread->first;
                    if (busyprocs.find(tryproc) == busyprocs.end()) { // not busy
                        busyprocs.insert(tryproc);
                        Vector<int> vreads;
                        vreads.reserve(whichread->second.size());
                        for (auto const& kv : whichread->second) {
                            vreads.push_back(kv.second);
                        }
                        if (tryproc == coordproc) {
                            iopfileindex = fi;
                            for (auto x : vreads) {
                                iopreads.push_back(x);
                            }
                        } else {
                            ParallelDescriptor::Send(vreads, tryproc, readtag);
                            ++reqspending;
                        }
                        availablefiles.erase(afilesiter);
                        afilesiter = availablefiles.begin();
                        break;
                    }
                }
                if (whichread == allreads[fi].end()) {
                    ++afilesiter;
                } else {
                    allreads[fi].erase(whichread);
                }
            }

            while (!iopreads.empty()) {
                int i = iopreads.front();
                detail::read_fab(fa[i], fod[i], name);
                --totalioreqs;
                iopreads.pop_front();
                if (iopreads.empty()) {
                    availablefiles.insert(iopfileindex);
                    busyprocs.erase(coordproc);
                }
                int doneflag;
                MPI_Status status;
                ParallelDescriptor::IProbe(MPI_ANY_SOURCE, donetag, doneflag, status);
                if (doneflag) {
                    break;
                }
            }

            if (reqspending > 0) {
                Vector<int> idone(2);
                ParallelDescriptor::Message rmess = ParallelDescriptor::Recv(idone, MPI_ANY_SOURCE,
                                                                             donetag);
                const int i = idone[0];
                const int procdone = rmess.pid();
                totalioreqs -= idone[1];
                --reqspending;
                busyprocs.erase(procdone);
                std::string const& fname = fod[i].m_name;
                const int fi = filenames.find(fname)->second;
                availablefiles.insert(fi);
            }
        }
    } else {
        Vector<int> recreads(nreqs, -1);
        Vector<int> idone(2);
        while (nreqs > 0) {
            ParallelDescriptor::Message rmess = ParallelDescriptor::Recv(recreads, coordproc,
                                                                         readtag);
            const auto nrmess = static_cast<int>(rmess.count());
            for (int ir = 0; ir < nrmess; ++ir) {
                int i = recreads[ir];
                detail::read_fab(fa[i], fod[i], name);
            }
            nreqs -= nrmess;
            idone[0] = recreads[0];
            idone[1] = nrmess;
            ParallelDescriptor::Send(idone, coordproc, donetag);
        }
    }
#else
    for (MFIter mfi(fa); mfi.isValid(); ++mfi) {
        detail::read_fab(fa[mfi], fod[mfi.index()], name);
    }
#endif
}

}

#endif /*BL_VISMF_H*/
