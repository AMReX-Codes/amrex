//BL_COPYRIGHT_NOTICE

#ifndef BL_VISMF_H
#define BL_VISMF_H

//
// $Id: VisMF.H,v 1.7 1997-11-10 03:55:46 lijewski Exp $
//

#include <sys/types.h>
#include <unistd.h>

#ifdef BL_USE_NEW_HFILES
#include <iostream>
using std::ostream;
using std::istream;
#else
#include <iostream.h>
#endif

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
        // Construct from a MultiFab -- fills in all but m_fod member.
        //
        Header (const MultiFab& mf, VisMF::How how);
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
    // Write a MultiFab to disk in a "smart" way.
    //
    static void Write (const MultiFab& mf,
                       const aString&  name,
                       VisMF::How      how);

private:
    //
    // These are disallowed.
    //
    VisMF (const VisMF&);
    VisMF& operator= (const VisMF&);
    //
    // Writes a FAB to disk returning the FabOnDisk structure.
    //
    static FabOnDisk Write (const FArrayBox& fab,
                            const aString&   filename,
                            ostream&         os);

    static void WriteOneFilePerCPU (const MultiFab& mf,
                                    const aString&  mf_name);

    static void WriteOneFilePerFab (const MultiFab& mf,
                                    const aString&  mf_name);
    //
    // Returns the CPU number as a four character string.
    //
    static aString TheCpuNumber ();
    //
    // The suffix appended to files containing FAB data.
    //
    static const aString FabFileSuffix;
    //
    // The suffix appended to files containing MultiFabs.
    //
    static const aString MultiFabHdrFileSuffix;
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

#endif /*BL_VISMF_H*/
