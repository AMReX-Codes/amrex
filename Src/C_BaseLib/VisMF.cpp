//BL_COPYRIGHT_NOTICE

//
// $Id: VisMF.cpp,v 1.22 1997-11-13 16:34:00 lijewski Exp $
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

typedef ParallelDescriptor PD;

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
        ParallelDescriptor::Abort("Write of VisMF::FabOnDisk failed");

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
        ParallelDescriptor::Abort("Read of VisMF::FabOnDisk failed");

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
        ParallelDescriptor::Abort("Write of Array<VisMF::FabOnDisk> failed");

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
        ParallelDescriptor::Abort("Read of Array<VisMF::FabOnDisk> failed");

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
        ParallelDescriptor::Abort("Write of Array<Array<Real>> failed");

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
        ParallelDescriptor::Abort("Read of Array<Array<Real>> failed");

    return is;
}

ostream&
operator<< (ostream&             os,
            const VisMF::Header& hd)
{
    //
    // Up the precision for the Reals in m_min and m_max.
    //
    int old_prec = os.precision(14);

    os << hd.m_vers     << '\n';
    os << int(hd.m_how) << '\n';
    os << hd.m_ncomp    << '\n';
    os << hd.m_ngrow    << '\n';

    hd.m_ba.writeOn(os); os << '\n';

    os << hd.m_fod      << '\n';
    os << hd.m_min      << '\n';
    os << hd.m_max      << '\n';

    os.precision(old_prec);

    if (!os.good())
        ParallelDescriptor::Abort("Write of VisMF::Header failed");

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
        ParallelDescriptor::Abort("Bad case in switch");
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
        ParallelDescriptor::Abort("Read of VisMF::Header failed");

    return is;
}

signed char*
VisMF::Large_IO_Buffer ()
{
    //
    // Used to implicitly manage allocation/deallocation of I/O buffer.
    //
    struct XXX
    {
        //
        // The only constructor
        //
        XXX (int size)
        {
            assert(size >= 0);
            if ((m_buffer = new signed char[size]) == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);
        }
        //
        // The destructor.
        //
        ~XXX () { delete [] m_buffer; m_buffer = 0; }
        //
        // The data we manage -- a buffer for I/O.
        //
        signed char* m_buffer;
    };

    static XXX IO_Buffer_Manager(VisMF::IO_Buffer_Size);

    assert(!(IO_Buffer_Manager.m_buffer == 0));

    return IO_Buffer_Manager.m_buffer;
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

//
// The more-or-less complete header only exists at PD::IOProcessor().
//

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
    //
    // Structure we use to pass min/max data to the IOProcessor.
    //
    // The variable length part of the message is 2*m_ncomp array of Reals.
    // The first part being the min values and the second the max values.
    //
    struct
    {
        int m_index; // Index of the FAB.
    }
    msg_hdr;

    int msg_hdr_size = sizeof(msg_hdr);
    PD::SetMessageHeaderSize(msg_hdr_size);
    //
    // Note that m_min and m_max are only calculated on CPU owning the fab.
    // We pass this data back to IOProcessor() so it sees the whole Header.
    //
    for (ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
    {
        int idx = mfi.index();

        m_min[idx].resize(m_ncomp);
        m_max[idx].resize(m_ncomp);

        const Box& valid_box = m_ba[idx];

        const FArrayBox& fab = mfi();

        assert(fab.box().contains(valid_box));

        for (long j = 0; j < m_ncomp; j++)
        {
            m_min[idx][j] = fab.min(valid_box,j);
            m_max[idx][j] = fab.max(valid_box,j);
        }

        if (!PD::IOProcessor())
        {
            msg_hdr.m_index = idx;

            Real* min_n_max = new Real[2 * m_ncomp];
            if (min_n_max == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);

            for (int i = 0; i < m_ncomp; i++)
            {
                min_n_max[i]           = m_min[idx][i];
                min_n_max[m_ncomp + i] = m_max[idx][i];
            }

            PD::SendData(PD::IOProcessor(),
                         &msg_hdr,
                         min_n_max,
                         2 * sizeof(Real) * m_ncomp);

            delete [] min_n_max;
        }
    }

    for (int len; PD::GetMessageHeader(len, &msg_hdr); )
    {
        assert(PD::IOProcessor());

        assert(len == 2 * sizeof(Real) * m_ncomp);

        m_min[msg_hdr.m_index].resize(m_ncomp);
        m_max[msg_hdr.m_index].resize(m_ncomp);

        Real* min_n_max = new Real[2 * m_ncomp];
        if (min_n_max == 0)
            BoxLib::OutOfMemory(__FILE__, __LINE__);

        PD::ReceiveData(min_n_max, len);

        for (int i = 0; i < m_ncomp; i++)
        {
            m_min[msg_hdr.m_index][i] = min_n_max[i];
            m_max[msg_hdr.m_index][i] = min_n_max[m_ncomp + i];

        }

        delete [] min_n_max;
    }
}

void
VisMF::WriteHeader (const aString& mf_name,
                    VisMF::Header& hdr)
{
    //
    // When running in parallel only one processor should do this I/O.
    //
    if (PD::IOProcessor())
    {
        aString MFHdrFileName = mf_name;

        MFHdrFileName += VisMF::MultiFabHdrFileSuffix;

        ofstream MFHdrFile;

#ifdef BL_T3E
        MFHdrFile.setbuf(VisMF::Large_IO_Buffer(), VisMF::IO_Buffer_Size);
        MFHdrFile.open(MFHdrFileName.c_str(), ios::out|ios::trunc);
#else
        MFHdrFile.open(MFHdrFileName.c_str(), ios::out|ios::trunc|ios::binary);
#endif
        MFHdrFile << hdr;
    }
}

//
// Returns the long as an aString of at least four digits.
//

aString
VisMF::ToString (long l)
{
    char buf[32];

    sprintf(buf, "%04ld", l);

    return aString(buf);
}

void
VisMF::Write (const MultiFab& mf,
              const aString&  mf_name,
              VisMF::How      how)
{
    //
    // Structure we use to pass FabOnDisk data to the IOProcessor.
    // The variable length part of the message is name of the on-disk FAB.
    //
    struct
    {
        int  m_index; // Index of the FAB.
        long m_head;  // Offset of the FAB in the file.
    }
    msg_hdr;

    VisMF::Header hdr(mf, VisMF::OneFilePerCPU);

    int msg_hdr_size = sizeof(msg_hdr);
    PD::SetMessageHeaderSize(msg_hdr_size);

    switch (how)
    {
    case OneFilePerCPU:
    {
        aString FabFileName = mf_name;

        FabFileName += VisMF::FabFileSuffix;
        FabFileName += VisMF::ToString(PD::MyProc());

        ofstream FabFile;

#ifdef BL_T3E
        FabFile.setbuf(VisMF::Large_IO_Buffer(), VisMF::IO_Buffer_Size);
        FabFile.open(FabFileName.c_str(), ios::out|ios::trunc);
#else
        FabFile.open(FabFileName.c_str(), ios::out|ios::trunc|ios::binary);
#endif
        if (!FabFile.good())
        {
            aString msg("Couldn't open file: ");
            msg += FabFileName;
            ParallelDescriptor::Abort(msg.c_str());
        }

        for (ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
        {
            hdr.m_fod[mfi.index()] = VisMF::Write(mfi(), FabFileName, FabFile);

            if (!PD::IOProcessor())
            {
                msg_hdr.m_index = mfi.index();
                msg_hdr.m_head  = hdr.m_fod[mfi.index()].m_head;

                PD::SendData(PD::IOProcessor(),
                             &msg_hdr,
                             FabFileName.c_str(),
                             FabFileName.length()+1); // Include NULL in MSG.
            }
        }
    }
    break;
    case OneFilePerFab:
    {
        for (ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
        {
            aString FabFileName = mf_name;

            FabFileName += VisMF::FabFileSuffix;
            FabFileName += VisMF::ToString(mfi.index());

            ofstream FabFile;

#ifdef BL_T3E
            FabFile.setbuf(VisMF::Large_IO_Buffer(), VisMF::IO_Buffer_Size);
            FabFile.open(FabFileName.c_str(), ios::out|ios::trunc);
#else
            FabFile.open(FabFileName.c_str(), ios::out|ios::trunc|ios::binary);
#endif
            if (!FabFile.good())
            {
                aString msg("Couldn't open file: ");
                msg += FabFileName;
                ParallelDescriptor::Abort(msg.c_str());
            }

            hdr.m_fod[mfi.index()] = VisMF::Write(mfi(), FabFileName, FabFile);

            assert(hdr.m_fod[mfi.index()].m_head == 0);

            if (!PD::IOProcessor())
            {
                msg_hdr.m_head  = 0;
                msg_hdr.m_index = mfi.index();

                PD::SendData(PD::IOProcessor(),
                             &msg_hdr,
                             FabFileName.c_str(),
                             FabFileName.length()+1); // Include NULL in MSG.
            }
        }
    }
    break;
    default:
        ParallelDescriptor::Abort("Bad case in switch");
    }

    for (int len; PD::GetMessageHeader(len, &msg_hdr); )
    {
        assert(PD::IOProcessor());

        char* fab_name = new char[len];
        if (fab_name == 0)
            BoxLib::OutOfMemory(__FILE__, __LINE__);

        PD::ReceiveData(fab_name, len);

        hdr.m_fod[msg_hdr.m_index].m_head = msg_hdr.m_head;
        hdr.m_fod[msg_hdr.m_index].m_name = fab_name;

        delete [] fab_name;
    }

    VisMF::WriteHeader(mf_name, hdr);
}

VisMF::VisMF (const aString& mf_name)
    :
    m_pa(PArrayManage)
{
    aString file = mf_name;

    file += VisMF::MultiFabHdrFileSuffix;

    ifstream ifs;

#ifdef BL_T3E
    ifs.setbuf(VisMF::Large_IO_Buffer(), VisMF::IO_Buffer_Size);
    ifs.open(file.c_str(), ios::in);
#else
    ifs.open(file.c_str(), ios::in|ios::binary);
#endif

    if (!ifs.good())
    {
        aString msg("Couldn't open file: ");
        msg += file;
        ParallelDescriptor::Abort(msg.c_str());
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

    ifstream ifs;

#ifdef BL_T3E
    ifs.setbuf(VisMF::Large_IO_Buffer(), VisMF::IO_Buffer_Size);
    ifs.open(file.c_str(), ios::in);
#else
    ifs.open(file.c_str(), ios::in|ios::binary);
#endif

    if (!ifs.good())
    {
        aString msg("Couldn't open file: ");
        msg += file;
        ParallelDescriptor::Abort(msg.c_str());
    }

    if (m_hdr.m_fod[i].m_head)
    {
        ifs.seekg(m_hdr.m_fod[i].m_head, ios::beg);
    }

    fab->readFrom(ifs);
}
