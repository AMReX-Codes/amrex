
#ifndef BL_VISMF_H
#define BL_VISMF_H

#include <iosfwd>
#include <string>

#include <REAL.H>
#include <FabArray.H>
#include <FArrayBox.H>

//
// File I/O for FabArray<FArrayBox>
//
//  Wrapper class for reading/writing FabArray<FArrayBox>s to disk in various "smart" ways.
//

class VisMF
{
public:
    //
    // How we write out FabArray<FArrayBox>s.
    //
    enum How { OneFilePerCPU, NFiles };
    //
    // Construct by reading in the on-disk VisMF of the specified name.
    // The MF on-disk is read lazily. The name here is the name of
    // the MF not the name of the on-disk file.
    //
    explicit VisMF (const std::string& mf_name);
    //
    // A structure containing info regarding an on-disk FAB.
    //
    struct FabOnDisk
    {
        //
        // The default constructor -- null out all fields.
        //
        FabOnDisk ();
        //
        // Constructor that sets the two values.
        //
        FabOnDisk (const std::string& name, long offset);
        //
        // The two data values in a FabOnDisk structure.
        //
        std::string m_name; // The name of file containing the FAB.
        long    m_head;     // Offset to start of FAB in file.
    };
    //
    // An on-disk FabArray<FArrayBox> contains this info in a header file.
    //
    struct Header
    {
        //
        // The current version of the FabArray<FArrayBox> Header code.
        //
        enum { Version = 1 };
        //
        // The default constructor.
        //
        Header ();
        //
        // Construct from a FabArray<FArrayBox>.
        //
        Header (const FabArray<FArrayBox>& mf, How how);
        //
        // The data.
        //
        int                  m_vers;  // The version # of the Header.
        How                  m_how;   // How the MF was written to disk.
        int                  m_ncomp; // Number of components in MF.
        int                  m_ngrow; // The number of ghost cells in MF.
        BoxArray             m_ba;    // The BoxArray of the MF.
        Array< FabOnDisk >   m_fod;   // FabOnDisk info for contained FABs.
        Array< Array<Real> > m_min;   // The min()s of each component of FABs.
        Array< Array<Real> > m_max;   // The max()s of each component of FABs.
    };
    //
    // The number of components in the on-disk FabArray<FArrayBox>.
    //
    int nComp () const;
    //
    // The grow factor of the on-disk FabArray<FArrayBox>.
    //
    int nGrow () const;
    //
    // # of FABs in the VisMF. Equal to # of Boxes in the BoxArray.
    //
    int size () const;
    //
    // The BoxArray of the on-disk FabArray<FArrayBox>.
    //
    const BoxArray& boxArray () const;
    //
    // The min of the FAB (in valid region) at specified index and component.
    //
    Real min (int fabIndex,
              int nComp) const;
    //
    // The max of the FAB (in valid region) at specified index and component.
    //
    Real max (int fabIndex,
              int nComp) const;

    /* The FAB at the specified index and component.
               Reads it from disk if necessary.
               This reads only the specified component.
    */
    const FArrayBox& GetFab (int fabIndex,
                             int compIndex) const;
    //
    // Delete()s the FAB at the specified index and component.
    //
    void clear (int fabIndex,
                int compIndex);
    //
    // Delete()s the FAB at the specified index (all components).
    //
    void clear (int fabIndex);
    //
    // Delete()s all the FABs.
    //
    void clear ();
    //
    // Write a FabArray<FArrayBox> to disk in a "smart" way.
    // Returns the total number of bytes written on this processor.
    // If set_ghost is true, sets the ghost cells in the FabArray<FArrayBox> to
    // one-half the average of the min and max over the valid region
    // of each contained FAB.
    //
    static long Write (const FabArray<FArrayBox>&    mf,
                       const std::string& name,
                       VisMF::How         how = OneFilePerCPU,
                       bool               set_ghost = false);
    //
    // Read a FabArray<FArrayBox> from disk written using VisMF::Write().
    // The FabArray<FArrayBox> mf must have been defined using the default
    // constructor.
    //
    static void Read (FabArray<FArrayBox>&          mf,
                      const std::string& name);
    static void Check (const std::string& name);
    //
    // We try to do I/O with buffers of this size.
    //
    enum { IO_Buffer_Size = 40960 * 32 };
    //
    // The type of a char buffer required by [p]setbuf().
    //
#ifdef BL_SETBUF_SIGNED_CHAR
    typedef signed char Setbuf_Char_Type;
#else
    typedef char Setbuf_Char_Type;
#endif
    //
    // A simple character buffer for setbuf() usage.
    //
    typedef Array<Setbuf_Char_Type> IO_Buffer;
    //
    // The file offset of the passed ostream.
    //
    static long FileOffset (std::ostream& os);
    //
    // Read the entire fab (all components).
    //
    FArrayBox* readFAB (int                fabIndex,
                        const std::string& mfName);
    //
    // Read the specified fab component.
    //
    FArrayBox* readFAB (int fabIndex,
                        int ncomp);

    static void SetNOutFiles (int noutfiles);

    static void SetMFFileInStreams (int nstreams);

    static int GetNOutFiles ();
    static int GetVerbose ();
    static void SetVerbose (int verbose);

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
                            long&              bytes);

    static long WriteHeader (const std::string& mf_name,
                             VisMF::Header&     hdr);
    //
    // Read the fab.
    // If ncomp == -1 reads the whole FAB.
    // Otherwise read just that component.
    //
    static FArrayBox* readFAB (int                fabIndex,
                               const std::string& mf_name,
                               const Header&      hdr,
			       int                ncomp = -1);
    //
    // Read the whole FAB into mf[fabIndex]
    //
    static void readFAB (FabArray<FArrayBox>&          mf,
			 int                fabIndex,
			 const std::string& mf_name,
			 const Header&      hdr);

    static std::string DirName (const std::string& filename);

    static std::string BaseName (const std::string& filename);
    //
    // Name of the FabArray<FArrayBox>.
    //
    std::string m_mfname;
    //
    // The VisMF header as read from disk.
    //
    Header m_hdr;
    //
    // We manage the FABs individually.
    //
    mutable Array< Array<FArrayBox*> > m_pa;
    //
    // The number of files to write for a FabArray<FArrayBox>.
    //
    static int nOutFiles;
    static int nMFFileInStreams;

    static int verbose;
};
//
// Write a FabOnDisk to an ostream in ASCII.
//
std::ostream& operator<< (std::ostream& os, const VisMF::FabOnDisk& fod);
//
// Read a FabOnDisk from an istream.
//
std::istream& operator>> (std::istream& is, VisMF::FabOnDisk& fod);
//
// Write an Array<FabOnDisk> to an ostream in ASCII.
//
std::ostream& operator<< (std::ostream& os, const Array<VisMF::FabOnDisk>& fa);
//
// Read an Array<FabOnDisk> from an istream.
//
std::istream& operator>> (std::istream& is, Array<VisMF::FabOnDisk>& fa);
//
// Write a VisMF::Header to an ostream in ASCII.
//
std::ostream& operator<< (std::ostream& os, const VisMF::Header& hd);
//
// Read a VisMF::Header from an istream.
//
std::istream& operator>> (std::istream& is, VisMF::Header& hd);

#endif /*BL_VISMF_H*/
