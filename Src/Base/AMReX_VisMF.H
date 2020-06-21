
#ifndef BL_VISMF_H
#define BL_VISMF_H

#include <iosfwd>
#include <string>
#include <fstream>
#include <thread>
#include <utility>
#include <cstdint>
#include <queue>

#include <AMReX_REAL.H>
#include <AMReX_FabArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabConv.H>
#include <AMReX_AsyncOut.H>

namespace amrex {

class NFilesIter;
class MultiFab;

/**
* \brief File I/O for FabArray<FArrayBox>.
*  Wrapper class for reading/writing FabArray<FArrayBox> objects to disk in various "smart" ways.
*/

class VisMF
{
public:
    //! We try to do I/O with buffers of this size.
    enum { IO_Buffer_Size = 262144 * 8 };
    //! The type of a char buffer required by [p]setbuf().
#ifdef BL_SETBUF_SIGNED_CHAR
    typedef signed char Setbuf_Char_Type;
#else
    typedef char Setbuf_Char_Type;
#endif
    //! A simple character buffer for setbuf() usage.
    typedef Vector<Setbuf_Char_Type> IO_Buffer;
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
    explicit VisMF (const std::string& fafab_name);
    ~VisMF ();
    //! A structure containing info regarding an on-disk FAB.
    struct FabOnDisk
    {
        //! The default constructor -- null out all fields.
        FabOnDisk ();
        //! Constructor that sets the two values.
        FabOnDisk (const std::string& name, Long offset);
        //! The two data values in a FabOnDisk structure.
        std::string m_name; //!< The name of file containing the FAB.
        Long    m_head;     //!< Offset to start of FAB in file.
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
        Header (const FabArray<FArrayBox>& fafab, VisMF::How how, Version version = Version_v1,
		bool calcMinMax = true, MPI_Comm = ParallelDescriptor::Communicator());

        Header (Header&& rhs) noexcept = default;

        //! Calculate the min and max arrays
        void CalculateMinMax(const FabArray<FArrayBox>& fafab,
                             int procToWrite = ParallelDescriptor::IOProcessorNumber(),
                             MPI_Comm = ParallelDescriptor::Communicator());
        //
        // The data.
        //
        int                  m_vers;  //!< The version of the Header.
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
        int rankToRead;
        int faIndex;
        Long fileOffset;
        Box box;

        FabReadLink();
        FabReadLink(int ranktoread, int faindex, Long fileoffset, const Box &b);
    };

    //! This structure is used to store file ifstreams that remain open
    struct PersistentIFStream
    {
        std::ifstream   *pstr;
        std::streampos   currentPosition;
        bool             isOpen;
        VisMF::IO_Buffer ioBuffer;

        PersistentIFStream();
        ~PersistentIFStream();
    };

    /**
    * \brief Open the stream if it is not already open
    * Close the stream if not persistent or forced
    * Close all open streams
    */
    static std::ifstream *OpenStream(const std::string &fileName);
    static void CloseStream(const std::string &fileName, bool forceClose = false);
    static void DeleteStream(const std::string &fileName);
    static void CloseAllStreams();
    static bool NoFabHeader(const VisMF::Header &hdr);

    //! The number of components in the on-disk FabArray<FArrayBox>.
    int nComp () const;
    //! The grow factor of the on-disk FabArray<FArrayBox>.
    int nGrow () const;
    IntVect nGrowVect () const;
    //! # of FABs in the VisMF. Equal to # of Boxes in the BoxArray.
    int size () const;
    //! The BoxArray of the on-disk FabArray<FArrayBox>.
    const BoxArray& boxArray () const;
    //! The min of the FAB (in valid region) at specified index and component.
    Real min (int fabIndex, int nComp) const;
    //! The min of the FabArray (in valid region) at specified component.
    Real min (int nComp) const;
    //! The max of the FAB (in valid region) at specified index and component.
    Real max (int fabIndex, int nComp) const;
    //! The max of the FabArray (in valid region) at specified component.
    Real max (int nComp) const;

    /**
    * \brief The FAB at the specified index and component.
    *         Reads it from disk if necessary.
    *         This reads only the specified component.
    */
    const FArrayBox& GetFab (int fabIndex,
                             int compIndex) const;
    //! Delete()s the FAB at the specified index and component.
    void clear (int fabIndex,
                int compIndex);
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
    static Long Write (const FabArray<FArrayBox> &fafab,
                       const std::string& name,
                       VisMF::How         how = NFiles,
                       bool               set_ghost = false);

    static void AsyncWrite (const FabArray<FArrayBox>& mf, const std::string& mf_name,
                            bool valid_cells_only = false);
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
    static void Read (FabArray<FArrayBox> &fafab,
                      const std::string &name,
                      const char *faHeader = nullptr,
                      int coordinatorProc = ParallelDescriptor::IOProcessorNumber(),
                      int allow_empty_mf = 0);

    //! Does FabArray exist?
    static bool Exist (const std::string &name);

    //! Read only the header of a FabArray, header will be resized here.
    static void ReadFAHeader (const std::string &fafabName,
                              Vector<char> &header);

    //! Check if the multifab is ok, false is returned if not ok
    static bool Check (const std::string &name);
    //! The file offset of the passed ostream.
    static Long FileOffset (std::ostream& os);
    //! Read the entire fab (all components).
    FArrayBox* readFAB (int fabIndex, const std::string& fafabName);
    //! Read the specified fab component.
    FArrayBox* readFAB (int fabIndex, int icomp);

    static int  GetNOutFiles ();
    static void SetNOutFiles (int newoutfiles, MPI_Comm comm = ParallelDescriptor::Communicator());

    static int  GetMFFileInStreams ()  { return nMFFileInStreams; }
    static void SetMFFileInStreams (int nstreams, MPI_Comm comm = ParallelDescriptor::Communicator());

    static int GetVerbose ()       { return verbose; }
    static void SetVerbose (int v) { verbose = v; }

    static VisMF::Header::Version GetHeaderVersion () { return currentVersion; }
    static void SetHeaderVersion (VisMF::Header::Version version)
                                                   { currentVersion = version; }

    static bool GetGroupSets () { return groupSets; }
    static void SetGroupSets (bool groupsets) { groupSets = groupsets; }

    static bool GetSetBuf () { return setBuf; }
    static void SetSetBuf (bool setbuf) { setBuf = setbuf; }

    static bool GetUseSingleRead () { return useSingleRead; }
    static void SetUseSingleRead (bool usesingleread) { useSingleRead = usesingleread; }

    static bool GetUseSingleWrite () { return useSingleWrite; }
    static void SetUseSingleWrite (bool usesinglewrite) { useSingleWrite = usesinglewrite; }

    static bool GetCheckFilePositions () { return checkFilePositions; }
    static void SetCheckFilePositions (bool cfp) { checkFilePositions = cfp; }

    static bool GetUsePersistentIFStreams () { return usePersistentIFStreams; }
    static void SetUsePersistentIFStreams (bool usepifs) { usePersistentIFStreams = usepifs; }

    static bool GetUseSynchronousReads () { return useSynchronousReads; }
    static void SetUseSynchronousReads (bool usepsr) { useSynchronousReads = usepsr; }

    static bool GetUseDynamicSetSelection () { return useDynamicSetSelection; }
    static void SetUseDynamicSetSelection (bool usedss) { useDynamicSetSelection = usedss; }

    static Long GetIOBufferSize () { return ioBufferSize; }
    static void SetIOBufferSize (Long iobuffersize) {
      BL_ASSERT(iobuffersize > 0);
      ioBufferSize = iobuffersize;
    }

    static void Initialize ();
    static void Finalize ();

private:
    //
    // These are disallowed.
    //
    VisMF (const VisMF&);
    VisMF& operator= (const VisMF&);

    static FabOnDisk Write (const FArrayBox&   fab,
                            const std::string& filename,
                            std::ostream&      os,
                            Long&              bytes);

    static Long WriteHeaderDoit (const std::string &fafab_name,
                                 VisMF::Header const &hdr);

    static Long WriteHeader (const std::string &fafab_name,
                             VisMF::Header     &hdr,
                             int procToWrite = ParallelDescriptor::IOProcessorNumber(),
                             MPI_Comm comm = ParallelDescriptor::Communicator());

    //! fileNumbers must be passed in for dynamic set selection [proc]
    static void FindOffsets (const FabArray<FArrayBox> &fafab,
			     const std::string &fafab_name,
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
    static FArrayBox *readFAB (int                fabIndex,
                               const std::string &fafab_name,
                               const Header      &hdr,
                               int                whichComp = -1);
    //! Read the whole FAB into fafab[fabIndex]
    static void readFAB (FabArray<FArrayBox> &fafab,
                         int                fabIndex,
                         const std::string &fafab_name,
                         const Header&      hdr);

    static std::string DirName (const std::string& filename);

    static std::string BaseName (const std::string& filename);

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
    static std::map<std::string, VisMF::PersistentIFStream> persistentIFStreams;
    //! The number of files to write for a FabArray<FArrayBox>.
    static int nOutFiles;
    static int nMFFileInStreams;

    static int verbose;
    static VisMF::Header::Version currentVersion;
    static bool groupSets;
    static bool setBuf;
    static bool useSingleRead;
    static bool useSingleWrite;
    static bool checkFilePositions;
    static bool usePersistentIFStreams;
    static bool useSynchronousReads;
    static bool useDynamicSetSelection;
    static bool allowSparseWrites;

    static Long ioBufferSize;   //!< ---- the settable buffer size
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

}

#endif /*BL_VISMF_H*/
