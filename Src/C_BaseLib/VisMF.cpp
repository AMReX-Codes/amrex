//BL_COPYRIGHT_NOTICE

//
// $Id: VisMF.cpp,v 1.14 1997-11-10 22:41:21 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <fstream>
using std::ios;
using std::ifstream;
using std::ofstream;
#else
#include <fstream.h>
#endif

#include <VisMF.H>

const aString VisMF::FabFileSuffix("_D_");

const aString VisMF::MultiFabHdrFileSuffix("_H");

const aString VisMF::FabOnDisk::Prefix("FOD:");

inline
void
GetTheChar (istream& is,
            char     ch)
{
    char c;
    is.get(c);
    assert(c == ch);
}

ostream&
operator<< (ostream&                os,
            const VisMF::FabOnDisk& fod)
{
    os << VisMF::FabOnDisk::Prefix << ' ' << fod.m_name << ' ' << fod.m_head;

    if (!os.good())
        BoxLib::Error("Write of VisMF::FabOnDisk failed");

    return os;
}

istream&
operator>> (istream&          is,
            VisMF::FabOnDisk& fod)
{
    aString str;
    is >> str;
    assert(str == VisMF::FabOnDisk::Prefix);

    GetTheChar(is, ' ');

    is >> fod.m_name;

    GetTheChar(is, ' ');

    is >> fod.m_head;

    if (!is.good())
        BoxLib::Error("Read of VisMF::FabOnDisk failed");

    return is;
}

ostream&
operator<< (ostream&                       os,
            const Array<VisMF::FabOnDisk>& fa)
{
    long i = 0, N = fa.length();

    os << N << '\n';

    for ( ; i < N; i++)
    {
        os << fa[i] << '\n';
    }

    if (!os.good())
        BoxLib::Error("Write of Array<VisMF::FabOnDisk> failed");

    return os;
}

istream&
operator>> (istream&                 is,
            Array<VisMF::FabOnDisk>& fa)
{
    long i = 0, N;

    is >> N;
    assert(N >= 0);

    fa.resize(N);

    for ( ; i < N; i++)
    {
        is >> fa[i];
        GetTheChar(is, '\n');
    }

    if (!is.good())
        BoxLib::Error("Read of Array<VisMF::FabOnDisk> failed");

    return is;
}

static
ostream&
operator<< (ostream&                    os,
            const Array< Array<Real> >& ar)
{
    long i = 0, N = ar.length(), M = (N == 0) ? 0 : ar[0].length();

    os << N << ',' << M << '\n';

    for ( ; i < N; i++)
    {
        assert(ar[i].length() == M);

        for (long j = 0; j < M; j++)
        {
            os << ar[i][j] << ',';
        }
        os << '\n';
    }

    if (!os.good())
        BoxLib::Error("Write of Array<Array<Real>> failed");

    return os;
}

static
istream&
operator>> (istream&              is,
            Array< Array<Real> >& ar)
{
    char ch;
    long i = 0, N, M;

    is >> N >> ch >> M;

    assert(N >= 0);
    assert(ch == ',');
    assert(M >= 0);

    GetTheChar(is, '\n');

    ar.resize(N);
    
    for ( ; i < N; i++)
    {
        ar[i].resize(M);

        for (long j = 0; j < M; j++)
        {
            is >> ar[i][j] >> ch;
            assert(ch == ',');
        }

        GetTheChar(is, '\n');
    }

    if (!is.good())
        BoxLib::Error("Read of Array<Array<Real>> failed");

    return is;
}

ostream&
operator<< (ostream&             os,
            const VisMF::Header& hd)
{
    os << hd.m_vers     << '\n';
    os << int(hd.m_how) << '\n';
    os << hd.m_ncomp    << '\n';
    os << hd.m_ngrow    << '\n';

    hd.m_ba.writeOn(os); os << '\n';

    os << hd.m_fod      << '\n';
    os << hd.m_min      << '\n';
    os << hd.m_max      << '\n';

    if (!os.good())
        BoxLib::Error("Write of VisMF::Header failed");

    return os;
}

istream&
operator>> (istream&       is,
            VisMF::Header& hd)
{
    is >> hd.m_vers;
    assert(hd.m_vers == VisMF::Header::Version);
    GetTheChar(is, '\n');

    int how;
    is >> how;
    switch (how)
    {
    case VisMF::OneFilePerCPU:
        hd.m_how = VisMF::OneFilePerCPU; break;
    case VisMF::OneFilePerFab:
        hd.m_how = VisMF::OneFilePerFab; break;
    default:
        BoxLib::Error("Bad case in switch");
    }
    GetTheChar(is, '\n');

    is >> hd.m_ncomp;
    assert(hd.m_ncomp >= 0);
    GetTheChar(is, '\n');

    is >> hd.m_ngrow;
    assert(hd.m_ngrow >= 0);
    GetTheChar(is, '\n');

    hd.m_ba = BoxArray(is);
    GetTheChar(is, '\n');

    is >> hd.m_fod;
    assert(hd.m_ba.length() == hd.m_fod.length());
    GetTheChar(is, '\n');

    is >> hd.m_min;
    GetTheChar(is, '\n');

    is >> hd.m_max;
    GetTheChar(is, '\n');

    assert(hd.m_ba.length() == hd.m_min.length());
    assert(hd.m_ba.length() == hd.m_max.length());

    if (!is.good())
        BoxLib::Error("Read of VisMF::Header failed");

    return is;
}

VisMF::FabOnDisk
VisMF::Write (const FArrayBox& fab,
              const aString&   filename,
              ostream&         os)
{
#ifdef __KCC
    VisMF::FabOnDisk fod(filename, os.tellp().offset());
#else
    VisMF::FabOnDisk fod(filename, os.tellp());
#endif

    fab.writeOn(os);

    return fod;
}

//
// This does not build a valid header.
//

VisMF::Header::Header ()
    :
    m_vers(0)
{}

VisMF::Header::Header (const MultiFab& mf,
                       VisMF::How      how)
    :
    m_vers(VisMF::Header::Version),
    m_how(how),
    m_ncomp(mf.nComp()),
    m_ngrow(mf.nGrow()),
    m_ba(mf.boxArray()),
    m_fod(m_ba.length()),
    m_min(m_ba.length()),
    m_max(m_ba.length())
{
    for (long i = 0, N = m_min.length(); i < N; i++)
    {
        m_min[i].resize(m_ncomp);
        m_max[i].resize(m_ncomp);

        const Box& valid_box = m_ba[i];

        for (long j = 0; j < m_ncomp; j++)
        {
            m_min[i][j] = mf[i].min(valid_box,j);
            m_max[i][j] = mf[i].max(valid_box,j);
        }
    }
}

void
VisMF::WriteHeader (const aString& mf_name,
                    VisMF::Header& hdr)
{
    //
    // When running in parallel only one processor should do this I/O.
    //
    if (ParallelDescriptor::IOProcessor())
    {
        //
        // TODO -- all headers must be passed to IOProcessor for reduction.
        //
        aString MFHdrFileName = mf_name;

        MFHdrFileName += VisMF::MultiFabHdrFileSuffix;

        ofstream MFHdrFile(MFHdrFileName.c_str());

        MFHdrFile << hdr;
    }
}

void
VisMF::WriteOneFilePerCPU (const MultiFab& mf,
                           const aString&  mf_name)
{
    VisMF::Header hdr(mf, VisMF::OneFilePerCPU);

    aString FabFileName = mf_name;

    FabFileName += VisMF::FabFileSuffix;
    FabFileName += VisMF::ToString(ParallelDescriptor::NProcs());

    {
        ofstream FabFile(FabFileName.c_str());

        if (!FabFile.good())
        {
            aString msg("Couldn't open file: ");
            msg += FabFileName;
            BoxLib::Error(msg.c_str());
        }

        for (ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
        {
            hdr.m_fod[mfi.index()] = VisMF::Write(mfi(), FabFileName, FabFile);
        }
    }

    VisMF::WriteHeader(mf_name, hdr);
}

//
// Returns the integer as an aString of at least four digits.
//

aString
VisMF::ToString (int i)
{
    char buf[32];

    sprintf(buf, "%04d", i);

    return aString(buf);
}

void
VisMF::WriteOneFilePerFab (const MultiFab& mf,
                           const aString&  mf_name)
{
    VisMF::Header hdr(mf, VisMF::OneFilePerFab);

    for (ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
    {
        aString FabFileName = mf_name;

        FabFileName += VisMF::FabFileSuffix;
        FabFileName += VisMF::ToString(mfi.index());

        ofstream FabFile(FabFileName.c_str());

        if (!FabFile.good())
        {
            aString msg("Couldn't open file: ");
            msg += FabFileName;
            BoxLib::Error(msg.c_str());
        }

        hdr.m_fod[mfi.index()] = VisMF::Write(mfi(), FabFileName, FabFile);
    }

    VisMF::WriteHeader(mf_name, hdr);
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

VisMF::VisMF (const aString& mf_name)
    :
    m_pa(PArrayManage)
{
    aString file = mf_name;

    file += VisMF::MultiFabHdrFileSuffix;

    ifstream ifs(file.c_str());

    if (!ifs.good())
    {
        aString msg("Couldn't open file: ");
        msg += file;
        BoxLib::Error(msg.c_str());
    }

    ifs >> m_hdr;

    m_pa.resize(m_hdr.m_ba.length());
}

void
VisMF::readFAB (int i) const
{
    Box fab_box = m_hdr.m_ba[i];

    if (m_hdr.m_ngrow)
    {
        fab_box.grow(m_hdr.m_ngrow);
    }

    FArrayBox* fab = new FArrayBox(fab_box, m_hdr.m_ncomp);
    if (fab == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);

    m_pa.set(i, fab);

    const aString& file = m_hdr.m_fod[i].m_name;

    ifstream ifs(file.c_str());

    if (!ifs.good())
    {
        aString msg("Couldn't open file: ");
        msg += file;
        BoxLib::Error(msg.c_str());
    }

    if (m_hdr.m_fod[i].m_head)
    {
        ifs.seekg(m_hdr.m_fod[i].m_head, ios::beg);
    }

    fab->readFrom(ifs);
}
