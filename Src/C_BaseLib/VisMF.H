//BL_COPYRIGHT_NOTICE

#ifndef BL_VISMF_H
#define BL_VISMF_H

//
// $Id: VisMF.H,v 1.21 1997-11-14 00:59:55 lijewski Exp $
//

#include <MultiFab.H>

//
// A smart wrapper class for reading/writing MultiFabs to disk.
//

class VisMF
{
public:
    //
    // How we write out MultiFabs.
    //
    enum How { OneFilePerCPU, OneFilePerFab };
    //
    // Construct by reading in the on-disk VisMF of the specified name.
    //
    // The MF on-disk is read lazily.
    //
    // The name here is the name of the MF not the name of the on-disk file.
    //
    explicit VisMF (const aString& mf_name);
    //
    // A structure containing info regarding an on-disk FAB.
    //
    struct FabOnDisk
    {
        //
        // The default constructor -- null out all fields.
        //
        FabOnDisk () {}
        //
        // Constructor that sets the two values.
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
    // An on-disk MultiFab contains this info in a header file.
    //
    struct Header
    {
        //
        // The current version of the MultiFab Header code.
        //
        enum { Version = 1 };
        //
        // The default constructor -- does nothing.
        //
        Header ();
        //
        // Construct from a MultiFab -- fills in all but m_fod member.
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
    // The number of components in the on-disk MultiFab.
    //
    int nComp () const;
    //
    // The grow factor of the on-disk MultiFab.
    //
    int nGrow () const;
    //
    // The # of FABs in the VisMF.
    //
    // Equal to the # of Boxes in the BoxArray.
    //
    int length () const;
    //
    // The BoxArray of the on-disk MultiFab.
    //
    const BoxArray& boxArray () const;
    //
    // The min of the FAB (in valid region) at specified index and component.
    //
    Real min (int i, int nComp) const;
    //
    // The max of the FAB (in valid region) at specified index and component.
    //
    Real max (int i, int nComp) const;
    //
    // The FAB at the specified index.
    //
    // Reads it from disk if necessary.
    //
    const FArrayBox& operator[] (int i) const;
    //
    // Write a MultiFab to disk in a "smart" way.
    //
    static void Write (const MultiFab& mf,
                       const aString&  name,
                       VisMF::How      how);
    //
    // We try to do I/O with buffers of this size.
    //
    enum { IO_Buffer_Size = 40960 * 32 };
    //
    // The type of a char buffer required by setbuf().
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

private:
    //
    // These are disallowed.
    //
    VisMF (const VisMF&);
    VisMF& operator= (const VisMF&);
    //
    // Private helper functions.
    //
    static FabOnDisk Write (const FArrayBox& fab,
                            const aString&   filename,
                            ostream&         os);

    static void WriteHeader (const aString& mf_name,
                             VisMF::Header& hdr);

    static aString ToString (long l);

    void readFAB (int i) const;
    //
    // The suffix appended to files containing FAB data.
    //
    static const aString FabFileSuffix;
    //
    // The suffix appended to files containing MultiFabs.
    //
    static const aString MultiFabHdrFileSuffix;
    //
    // The VisMF header as read from disk.
    //
    Header m_hdr;
    //
    // We manage the FABs individually.
    //
    mutable PArray<FArrayBox> m_pa;
};

//
// Write a FabOnDisk to an ostream in ASCII.
//
ostream& operator<< (ostream& os, const VisMF::FabOnDisk& fod);

//
// Read a FabOnDisk from an istream.
//
istream& operator>> (istream& is, VisMF::FabOnDisk& fod);

//
// Write an Array<FabOnDisk> to an ostream in ASCII.
//
ostream& operator<< (ostream& os, const Array<VisMF::FabOnDisk>& fa);

//
// Read an Array<FabOnDisk> from an istream.
//
istream& operator>> (istream& is, Array<VisMF::FabOnDisk>& fa);

//
// Write a VisMF::Header to an ostream in ASCII.
//
ostream& operator<< (ostream& os, const VisMF::Header& hd);

//
// Read a VisMF::Header from an istream.
//
istream& operator>> (istream& is, VisMF::Header& hd);

//
// Inlines.
//

inline
int
VisMF::nComp () const
{
    return m_hdr.m_ncomp;
}

inline
int
VisMF::nGrow () const
{
    return m_hdr.m_ngrow;
}

inline
int
VisMF::length () const
{
    return m_hdr.m_ba.length();
}

inline
const BoxArray&
VisMF::boxArray () const
{
    return m_hdr.m_ba;
}

inline
Real
VisMF::min (int i, int nComp) const
{
    assert(0 <= i && i < m_hdr.m_ba.length());
    assert(0 <= nComp && nComp < m_hdr.m_ncomp);
    return m_hdr.m_min[i][nComp];
}

inline
Real
VisMF::max (int i, int nComp) const
{
    assert(0 <= i && i < m_hdr.m_ba.length());
    assert(0 <= nComp && nComp < m_hdr.m_ncomp);
    return m_hdr.m_max[i][nComp];
}

inline
const FArrayBox&
VisMF::operator[] (int i) const
{
    if (!m_pa.defined(i))
    {
        readFAB(i);
    }
    return m_pa[i];
}

#endif /*BL_VISMF_H*/
