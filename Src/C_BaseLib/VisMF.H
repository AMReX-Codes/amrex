#ifndef BL_VISMF_H
#define BL_VISMF_H

//
// $Id: VisMF.H,v 1.46 2001-07-19 16:57:37 lijewski Exp $
//

#include <MultiFab.H>

//
//@Man:
//@Memo: File I/O for MultiFabs.
/*@Doc:

  Wrapper class for reading/writing MultiFabs to disk in various "smart" ways.
*/

class VisMF
{
public:
    //
    //@ManDoc: How we write out MultiFabs.
    //
    enum How { OneFilePerCPU };

    /*@ManDoc: Construct by reading in the on-disk VisMF of the specified name.
               The MF on-disk is read lazily. The name here is the name of
               the MF not the name of the on-disk file.
    */
    explicit VisMF (const aString& mf_name);
    //
    //@ManDoc: A structure containing info regarding an on-disk FAB.
    //
    struct FabOnDisk
    {
        //
        //@ManDoc: The default constructor -- null out all fields.
        //
        FabOnDisk () {}
        //
        //@ManDoc: Constructor that sets the two values.
        //
        FabOnDisk (const aString& name, long offset)
            :
            m_name(name), m_head(offset) {}
        //
        // The two data values in a FabOnDisk structure.
        //
        aString m_name;  // The name of file containing the FAB.
        long    m_head;  // Offset to start of FAB in file.
        //
        // Prefixed to lines of FabOnDisk I/O.
        //
        static const aString Prefix;
    };
    //
    //@ManDoc: An on-disk MultiFab contains this info in a header file.
    //
    struct Header
    {
        //
        // The current version of the MultiFab Header code.
        //
        enum { Version = 1 };
        //
        //@ManDoc: The default constructor.
        //
        Header ();
        //
        //@ManDoc: Construct from a MultiFab.
        //
        Header (const MultiFab& mf, How how);
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
    //@ManDoc: The number of components in the on-disk MultiFab.
    //
    int nComp () const;
    //
    //@ManDoc: The grow factor of the on-disk MultiFab.
    //
    int nGrow () const;
    //
    //@ManDoc: # of FABs in the VisMF. Equal to # of Boxes in the BoxArray.
    //
    int length () const;
    //
    //@ManDoc: The BoxArray of the on-disk MultiFab.
    //
    const BoxArray& boxArray () const;
    //
    //@ManDoc: The min of the FAB (in valid region) at specified index and component.
    //
    Real min (int fabIndex,
              int nComp) const;
    //
    //@ManDoc: The max of the FAB (in valid region) at specified index and component.
    //
    Real max (int fabIndex,
              int nComp) const;

    /*@ManDoc: The FAB at the specified index and component.
               Reads it from disk if necessary.
               This reads only the specified component.
    */
    const FArrayBox& GetFab (int fabIndex,
                             int compIndex) const;
    //
    //@ManDoc: Delete()s the FAB at the specified index and component.
    //
    void clear (int fabIndex,
                int compIndex);
    //
    //@ManDoc: Delete()s the FAB at the specified index (all components).
    //
    void clear (int fabIndex);
    //
    //@ManDoc: Delete()s all the FABs.
    //
    void clear ();

    /*@ManDoc: Write a MultiFab to disk in a "smart" way.
               Returns the total number of bytes written on this processor.
               If set\_ghost is true, sets the ghost cells in the MultiFab to
               one-half the average of the min and max over the valid region
               of each contained FAB.
     */
    static long Write (const MultiFab& mf,
                       const aString&  name,
                       VisMF::How      how = OneFilePerCPU,
                       bool            set_ghost = false);

    /*@ManDoc: Read a MultiFab from disk written using VisMF::Write().
               The MultiFab `mf' must have been defined using the default
               constructor.
    */
    static void Read (MultiFab&      mf,
                      const aString& name);
    //
    //@ManDoc: We try to do I/O with buffers of this size.
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
    //@ManDoc: A simple character buffer for setbuf() usage.
    //
    typedef Array<Setbuf_Char_Type> IO_Buffer;
    //
    //@ManDoc: The file offset of the passed ostream.
    //
    static long FileOffset (std::ostream& os);
    //
    //@ManDoc: Read the entire fab (all components).
    //
    FArrayBox* readFAB (int            fabIndex,
                        const aString& mfName);
    //
    //@ManDoc: Read the specified fab component.
    //
    FArrayBox* readFAB (int fabIndex,
                        int ncomp);

private:
    //
    // These are disallowed.
    //
    VisMF (const VisMF&);
    VisMF& operator= (const VisMF&);

    static FabOnDisk Write (const FArrayBox& fab,
                            const aString&   filename,
                            std::ostream&    os,
                            long&            bytes);

    static long WriteHeader (const aString& mf_name,
                             VisMF::Header& hdr);

    //
    // Read the fab.
    // If ncomp == -1 reads the whole FAB.
    // Otherwise read just that component.
    //
    static FArrayBox* readFAB (int            fabIndex,
                               const aString& mf_name,
                               const Header&  hdr,
			       int            ncomp = -1);

    static aString DirName (const aString& filename);

    static aString BaseName (const aString& filename);
    //
    // The suffix appended to files containing FAB data.
    //
    static const aString FabFileSuffix;
    //
    // The suffix appended to files containing MultiFabs.
    //
    static const aString MultiFabHdrFileSuffix;
    //
    // Name of the MultiFab.
    //
    aString m_mfname;
    //
    // The VisMF header as read from disk.
    //
    Header m_hdr;
    //
    // We manage the FABs individually.
    //
    mutable Array< Array<FArrayBox*> > m_pa;
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
