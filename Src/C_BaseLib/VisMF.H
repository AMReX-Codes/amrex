//BL_COPYRIGHT_NOTICE

#ifndef BL_VISMF_H
#define BL_VISMF_H

//
// $Id: VisMF.H,v 1.2 1997-11-09 00:58:56 lijewski Exp $
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
        FabOnDisk () : m_head(0), m_data(0) {}
        //
        // Constructor that sets the file name field and  nulls out others.
        //
        FabOnDisk (const aString name) : m_name(name), m_head(0), m_data(0) {}

        aString m_name; // The name of file containing the FAB.
        off_t   m_head; // Offset of start of FAB header in file.
        off_t   m_data; // Offset of start of FAB data in file.
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
        // Construct from a MultiFab -- fills in all but m_fad member.
        //
        Header (const MultiFab& mf, VisMF::How how);
        //
        // Writes the Header in ASCII to the supplied ostream.
        //
        void writeOn (ostream& os) const;
        //
        // What's in the header.
        //
        int                  m_vers;  // The version # of the Header.
        How                  m_how;   // How the MF was written to disk.
        int                  m_ncomp; // Number of components in MF.
        BoxArray             m_ba;    // The BoxArray of the MF.
        Array< FabOnDisk >   m_fad;   // FabOnDisk info for contained FABs.
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
    // Writes a FAB to disk.
    //
    static FabOnDisk Write (const FArrayBox& fab,
                            const aString&   filename,
                            ostream&         os);

    static void WriteOneFilePerCPU (const MultiFab& mf,
                                    const aString&  mf_name);

    static void WriteOneFilePerFab (const MultiFab& mf,
                                    const aString&  mf_name);
};

//
// Write a FabOnDisk to an ostream in ASCII.
//

ostream& operator<< (ostream& os, const VisMF::FabOnDisk& fad);

//
// Write an Array<FabOnDisk> to an ostream in ASCII.
//

ostream& operator<< (ostream& os, const Array<VisMF::FabOnDisk>& fa);

#endif /*BL_VISMF_H*/
