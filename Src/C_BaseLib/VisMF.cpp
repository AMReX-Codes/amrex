//BL_COPYRIGHT_NOTICE

//
// $Id: VisMF.cpp,v 1.3 1997-11-09 19:45:29 lijewski Exp $
//

#include <VisMF.H>

const aString VisMF::FabFileSuffix("_data_");
const aString VisMF::MultiFabHdrFileSuffix("_hdr_");
const aString VisMF::FabOnDisk::IOPrefix("FabOnDisk:");

ostream&
operator<< (ostream&                os,
            const VisMF::FabOnDisk& fad)
{
    const char SPC = ' ';

    os << VisMF::FabOnDisk::IOPrefix << SPC
       << fad.m_name                 << SPC
       << fad.m_head                 << SPC
       << fad.m_data;

    if (!os.good())
        BoxLib::Error("Write of VisMF::FabOnDisk failed");

    return os;
}

ostream&
operator<< (ostream&                       os,
            const Array<VisMF::FabOnDisk>& fa)
{
    for (int i = 0, N = fa.length(); i < N; i++)
    {
        os << fa[i] << '\n';
    }

    if (!os.good())
        BoxLib::Error("Write of Array<VisMF::FabOnDisk> failed");

    return os;
}

static
ostream&
operator<< (ostream&                    os,
            const Array< Array<Real> >& ar)
{
    for (int i = 0, N = ar.length(); i < N; i++)
    {
        for (int j = 0, M = ar[i].length(); j < M; j++)
        {
            os << ar[i][j] << ' ';
        }
        os << '\n';
    }

    if (!os.good())
        BoxLib::Error("Write of Array<Array<Real>> failed");

    return os;
}

void
VisMF::Header::writeOn (ostream& os) const
{
    os << m_vers     << '\n';
    os << int(m_how) << '\n';
    os << m_ncomp    << '\n';
    os << m_ba       << '\n';
    os << m_fad      << '\n';
    os << m_min      << '\n';
    os << m_max      << '\n';

    if (!os.good())
        BoxLib::Error("VisMF::Header::WriteOn() failed");
}

VisMF::FabOnDisk
VisMF::Write (const FArrayBox& fab,
              const aString&   filename,
              ostream&         os)
{
    VisMF::FabOnDisk fad(filename);

    //
    // TODO -- set the other two fields.
    //

    return fad;
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
    for (int i = 0, N = m_min.length(); i < N; i++)
    {
        m_min[i].resize(m_ncomp);
        m_max[i].resize(m_ncomp);
    }

    for (int i = 0, N = m_min.length(); i < N; i++)
    {
        const Box& subbox = m_ba[i];

        const FArrayBox& fab = mf[i];

        for (int j = 0; j < m_ncomp; j++)
        {
            //
            // TODO -- is using the subbox here correct?
            //
            m_min[i][j] = fab.min(subbox,j);
            m_max[i][j] = fab.max(subbox,j);
        }
    }
}

aString
VisMF::TheCpuNumber ()
{
    //
    // TODO -- finish this!!!
    //
    return aString("0000");
}

void
VisMF::WriteOneFilePerCPU (const MultiFab& mf,
                           const aString&  mf_name)
{
    //
    // TODO -- make this work in Parallel!!!
    //
    VisMF::Header hdr(mf, VisMF::OneFilePerCPU);

    aString fab_file_name = mf_name;

    fab_file_name += VisMF::FabFileSuffix;
    fab_file_name += VisMF::TheCpuNumber();

    {
        ofstream fab_file(fab_file_name.c_str());

        for (int i = 0, N = mf.length(); i < N; i++)
        {
            hdr.m_fad[i] = VisMF::Write(mf[i], fab_file_name, fab_file);
        }
    }
    //
    // The Header gets written last to its own file.
    //
    aString mf_hdr_file_name = mf_name;

    mf_hdr_file_name += VisMF::MultiFabHdrFileSuffix;

    ofstream mf_hdr_file(mf_hdr_file_name.c_str());

    hdr.writeOn(mf_hdr_file);
}

void
VisMF::WriteOneFilePerFab (const MultiFab& mf,
                           const aString&  mf_name)
{
    //
    // TODO -- make this work in Parallel!!!
    //
    VisMF::Header hdr(mf, VisMF::OneFilePerFab);

    BoxLib::Error("Not implemented");
}

void
VisMF::Write (const MultiFab& mf,
              const aString&  mf_name,
              VisMF::How      how)
{
    switch (how)
    {
    case OneFilePerCPU:
        WriteOneFilePerCPU(mf, mf_name); break;
    case OneFilePerFab:
        WriteOneFilePerFab(mf, mf_name); break;
    default:
        BoxLib::Error("Bad case in switch");
    }
}
