
//
// $Id: VisMF.cpp,v 1.1 1997-11-08 23:24:17 lijewski Exp $
//

#include <VisMF.H>

VisMF::FabOnDisk
VisMF::Write (const FArrayBox& fab,
              const aString&   name,
              ostream&         os)
{
}

VisMF::Header::Header (const MultiFab& mf,
                       VisMF::How      how)
    :
    m_vers(VisMF::Header::Version),
    m_how(how),
    m_ncomp(mf.nComp()),
    m_ba(mf.boxArray()),
    m_fad(m_ba.length()),
    m_min(m_ba.length()),
    m_max(m_ba.length())
{
    //
    // Fill in m_fad ...
    //



    //
    // Fill in m_min and m_max ...
    //
}

void
VisMF::WriteOneFilePerCPU (const MultiFab& mf,
                           const aString&  name)
{
    VisMF::Header hdr(mf, VisMF::OneFilePerCPU);
}

void
VisMF::WriteOneFilePerFab (const MultiFab& mf,
                           const aString&  name)
{
    VisMF::Header hdr(mf, VisMF::OneFilePerFab);

    BoxLib::Error("Not implemented");
}

void
VisMF::Write (const MultiFab& mf,
              const aString&  name,
              VisMF::How      how)
{
    switch (how)
    {
    case OneFilePerCPU:
        WriteOneFilePerCPU(mf, name); break;
    case OneFilePerFab:
        WriteOneFilePerFab(mf, name); break;
    default:
        BoxLib::Error("Bad case in switch");
    }
}
