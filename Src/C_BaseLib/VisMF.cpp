//BL_COPYRIGHT_NOTICE

//
// $Id: VisMF.cpp,v 1.52 1998-07-07 19:43:29 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstdio>
#include <fstream>
using std::ios;
using std::ifstream;
using std::ofstream;
#else
#include <stdio.h>
#include <fstream.h>
#endif

//
// This MUST be defined if don't have pubsetbuf() in I/O Streams Library.
//
#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

#include <Utility.H>
#include <VisMF.H>
#include <Tracer.H>

#ifdef BL_USE_MPI
#include <mpi.h>
#endif

const aString VisMF::FabFileSuffix("_D_");

const aString VisMF::MultiFabHdrFileSuffix("_H");

const aString VisMF::FabOnDisk::Prefix("FabOnDisk:");

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

    is >> fod.m_name;
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

    ar.resize(N);
    
    for ( ; i < N; i++)
    {
        ar[i].resize(M);

        for (long j = 0; j < M; j++)
        {
            is >> ar[i][j] >> ch;
            assert(ch == ',');
        }
    }

    if (!is.good())
        BoxLib::Error("Read of Array<Array<Real>> failed");

    return is;
}

ostream&
operator<< (ostream&             os,
            const VisMF::Header& hd)
{
    //
    // Up the precision for the Reals in m_min and m_max.
    //
    int old_prec = os.precision(15);

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
        BoxLib::Error("Write of VisMF::Header failed");

    return os;
}

istream&
operator>> (istream&       is,
            VisMF::Header& hd)
{
    is >> hd.m_vers;
    assert(hd.m_vers == VisMF::Header::Version);

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

    is >> hd.m_ncomp;
    assert(hd.m_ncomp >= 0);

    is >> hd.m_ngrow;
    assert(hd.m_ngrow >= 0);

    hd.m_ba = BoxArray(is);

    is >> hd.m_fod;
    assert(hd.m_ba.length() == hd.m_fod.length());

    is >> hd.m_min;
    is >> hd.m_max;

    assert(hd.m_ba.length() == hd.m_min.length());
    assert(hd.m_ba.length() == hd.m_max.length());

    if (!is.good())
        BoxLib::Error("Read of VisMF::Header failed");

    return is;
}

aString
VisMF::BaseName (const aString& filename)
{
    assert(filename[filename.length() - 1] != '/');

    if (char* slash = strrchr(filename.c_str(), '/'))
    {
        //
        // Got at least one slash -- give'm the following tail.
        //
        return aString(slash + 1);
    }
    else
    {
        //
        // No leading directory portion to name.
        //
        return filename;
    }
}

aString
VisMF::DirName (const aString& filename)
{
    assert(filename[filename.length() - 1] != '/');

    static const aString TheNullString("");

    const char* str = filename.c_str();    

    if (char* slash = strrchr(str, '/'))
    {
        //
        // Got at least one slash -- give'm the dirname including last slash.
        //
        int len = (slash - str) + 1;

        char* buf = new char[len+1];

        strncpy(buf, str, len);

        buf[len] = 0; // Stringify

        aString dirname = buf;

        delete [] buf;

        return dirname;
    }
    else
    {
        //
        // No directory name here.
        //
        return TheNullString;
    }
}

VisMF::FabOnDisk
VisMF::Write (const FArrayBox& fab,
              const aString&   filename,
              ostream&         os,
              long&            bytes)
{
    VisMF::FabOnDisk fab_on_disk(filename, VisMF::FileOffset(os));

    fab.writeOn(os);
    //
    // Add in the number of bytes in the FAB including the FAB header.
    //
    bytes += (VisMF::FileOffset(os) - fab_on_disk.m_head);
    
    return fab_on_disk;
}

//
// This does not build a valid header.
//

VisMF::Header::Header ()
    :
    m_vers(0)
{}

//
// The more-or-less complete header only exists at IOProcessor().
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
    ParallelDescriptor::SetMessageHeaderSize(sizeof(int));
    //
    // Note that m_min and m_max are only calculated on CPU owning the fab.
    // We pass this data back to IOProcessor() so it sees the whole Header.
    //
    Array<Real> min_n_max(2 * m_ncomp);

    const int IoProc = ParallelDescriptor::IOProcessorNumber();

    for (ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
    {
        int idx = mfi.index();

        m_min[idx].resize(m_ncomp);
        m_max[idx].resize(m_ncomp);

        assert(mfi().box().contains(m_ba[idx]));

        for (long j = 0; j < m_ncomp; j++)
        {
            m_min[idx][j] = mfi().min(m_ba[idx],j);
            m_max[idx][j] = mfi().max(m_ba[idx],j);
        }

        if (!ParallelDescriptor::IOProcessor())
        {
            for (int i = 0; i < m_ncomp; i++)
            {
                min_n_max[i]         = m_min[idx][i];
                min_n_max[m_ncomp+i] = m_max[idx][i];
            }
#ifdef BL_USE_MPI
            int rc = MPI_Ssend(min_n_max.dataPtr(),
                               2 * m_ncomp,
                               mpi_data_type(min_n_max.dataPtr()),
                               IoProc,
                               //
                               // We use the index as the tag for uniqueness.
                               //
                               idx,
                               MPI_COMM_WORLD);
            if (!(rc == MPI_SUCCESS))
                ParallelDescriptor::Abort(rc);
#else
            ParallelDescriptor::SendData(IoProc,
                                         &idx,
                                         min_n_max.dataPtr(),
                                         2 * sizeof(Real) * m_ncomp);
#endif
        }
    }
#ifdef BL_USE_MPI
    if (ParallelDescriptor::IOProcessor())
    {
        MPI_Status        status;
        const Array<int>& procmap = mf.DistributionMap().ProcessorMap();

        for (int idx = 0, N = procmap.length(), rc = 0; idx < N; idx++)
        {
            if (!(procmap[idx] == IoProc))
            {
                if ((rc = MPI_Recv(min_n_max.dataPtr(),
                                   2 * m_ncomp,
                                   mpi_data_type(min_n_max.dataPtr()),
                                   //
                                   // We recieve one message for each
                                   // index that IoProc doesn't own.
                                   //
                                   procmap[idx],
                                   idx,
                                   MPI_COMM_WORLD,
                                   &status)) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);

                m_min[idx].resize(m_ncomp);
                m_max[idx].resize(m_ncomp);

                for (int i = 0; i < m_ncomp; i++)
                {
                    m_min[idx][i] = min_n_max[i];
                    m_max[idx][i] = min_n_max[m_ncomp+i];
                }
            }
        }
    }
#else
    for (int len, idx; ParallelDescriptor::GetMessageHeader(len, &idx); )
    {
        assert(len == 2 * sizeof(Real) * m_ncomp);

        m_min[idx].resize(m_ncomp);
        m_max[idx].resize(m_ncomp);

        ParallelDescriptor::ReceiveData(min_n_max.dataPtr(), len);

        for (int i = 0; i < m_ncomp; i++)
        {
            m_min[idx][i] = min_n_max[i];
            m_max[idx][i] = min_n_max[m_ncomp+i];
        }
    }
#endif
}

long
VisMF::WriteHeader (const aString& mf_name,
                    VisMF::Header& hdr)
{
    long bytes = 0;
    //
    // When running in parallel only one processor should do this I/O.
    //
    if (ParallelDescriptor::IOProcessor())
    {
        aString MFHdrFileName = mf_name;

        MFHdrFileName += VisMF::MultiFabHdrFileSuffix;

        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

        ofstream MFHdrFile;

        MFHdrFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.length());

        MFHdrFile.open(MFHdrFileName.c_str(), ios::out|ios::trunc);

        if (!MFHdrFile.good())
            Utility::FileOpenFailed(MFHdrFileName);

        MFHdrFile << hdr;
        //
        // Add in the number of bytes written out in the Header.
        //
        bytes += VisMF::FileOffset(MFHdrFile);
    }
    return bytes;
}

long
VisMF::Write (const MultiFab& mf,
              const aString&  mf_name,
              VisMF::How      how)
{
    assert(mf_name[mf_name.length() - 1] != '/');

    const int IoProc = ParallelDescriptor::IOProcessorNumber();

    long bytes = 0;
    //
    // Structure we use to pass FabOnDisk data to the IOProcessor.
    // The variable length part of the message is the name of the on-disk FAB.
    //
    struct
    {
        int  m_index; // Index of the FAB.
        long m_head;  // Offset of the FAB in the file.
    }
    msg_hdr;

    VisMF::Header hdr(mf, VisMF::OneFilePerCPU);

    ParallelDescriptor::SetMessageHeaderSize(sizeof(msg_hdr));

    char buf[sizeof(int) + 1];

#ifdef BL_USE_MPI
    const int NPBUF = 512;

    Array<char> packedbuffer(NPBUF);

    char* pbuf = packedbuffer.dataPtr();
#endif

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    switch (how)
    {
    case OneFilePerCPU:
    {
        aString FullFileName = mf_name;

        FullFileName += VisMF::FabFileSuffix;
        sprintf(buf, "%04ld", ParallelDescriptor::MyProc());
        FullFileName += buf;

        ofstream FabFile;

        FabFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.length());

        FabFile.open(FullFileName.c_str(), ios::out|ios::trunc);

        if (!FabFile.good())
            Utility::FileOpenFailed(FullFileName);

        aString TheBaseName = VisMF::BaseName(FullFileName);

        for (ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
        {
            hdr.m_fod[mfi.index()] = VisMF::Write(mfi(), TheBaseName, FabFile, bytes);

            if (!ParallelDescriptor::IOProcessor())
            {
                msg_hdr.m_index = mfi.index();
                msg_hdr.m_head  = hdr.m_fod[mfi.index()].m_head;
#ifdef BL_USE_MPI
                //
                // Pack FAB offset, length of FAB name, and name into buffer.
                //
                int rc = 0, pos = 0;

                if ((rc = MPI_Pack(&msg_hdr.m_head,
                                   1,
                                   MPI_LONG,
                                   pbuf,
                                   NPBUF,
                                   &pos,
                                   MPI_COMM_WORLD)) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);

                int len = TheBaseName.length()+1; // Includes the NULL.

                if ((rc = MPI_Pack(&len,
                                   1,
                                   MPI_INT,
                                   pbuf,
                                   NPBUF,
                                   &pos,
                                   MPI_COMM_WORLD)) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);

                if ((rc = MPI_Pack((void*)TheBaseName.c_str(),
                                   len,
                                   MPI_CHAR,
                                   pbuf,
                                   NPBUF,
                                   &pos,
                                   MPI_COMM_WORLD)) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);

                if ((rc = MPI_Ssend(pbuf,
                                    pos,
                                    MPI_PACKED,
                                    IoProc,
                                    //
                                    // We use index as the tag for uniqueness.
                                    //
                                    mfi.index(),
                                    MPI_COMM_WORLD)) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);
#else
                ParallelDescriptor::SendData(IoProc,
                                             &msg_hdr,
                                             TheBaseName.c_str(),
                                             //
                                             // Include NULL in MSG.
                                             //
                                             TheBaseName.length()+1);
#endif
            }
        }
        if (VisMF::FileOffset(FabFile) <= 0)
            Utility::UnlinkFile(FullFileName);
    }
    break;
    case OneFilePerFab:
    {
        for (ConstMultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
        {
            aString FullFileName = mf_name;

            FullFileName += VisMF::FabFileSuffix;
            sprintf(buf, "%04ld", mfi.index());
            FullFileName += buf;

            ofstream FabFile;

            FabFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.length());

            FabFile.open(FullFileName.c_str(), ios::out|ios::trunc);

            if (!FabFile.good())
                Utility::FileOpenFailed(FullFileName);

            aString TheBaseName = VisMF::BaseName(FullFileName);

            hdr.m_fod[mfi.index()] = VisMF::Write(mfi(), TheBaseName, FabFile, bytes);

            assert(hdr.m_fod[mfi.index()].m_head == 0);

            if (!ParallelDescriptor::IOProcessor())
            {
                msg_hdr.m_head  = 0;
                msg_hdr.m_index = mfi.index();
#ifdef BL_USE_MPI
                //
                // Pack FAB offset, length of name, and name into buffer.
                //
                int rc = 0, pos = 0;

                if ((rc = MPI_Pack(&msg_hdr.m_head,
                                   1,
                                   MPI_LONG,
                                   pbuf,
                                   NPBUF,
                                   &pos,
                                   MPI_COMM_WORLD)) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);

                int len = TheBaseName.length()+1; // Includes the NULL.

                if ((rc = MPI_Pack(&len,
                                   1,
                                   MPI_INT,
                                   pbuf,
                                   NPBUF,
                                   &pos,
                                   MPI_COMM_WORLD)) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);

                if ((rc = MPI_Pack((void*)TheBaseName.c_str(),
                                   len,
                                   MPI_CHAR,
                                   pbuf,
                                   NPBUF,
                                   &pos,
                                   MPI_COMM_WORLD)) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);

                if ((rc = MPI_Ssend(pbuf,
                                    pos,
                                    MPI_PACKED,
                                    IoProc,
                                    //
                                    // We use index as the tag for uniqueness.
                                    //
                                    mfi.index(),
                                    MPI_COMM_WORLD)) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);
#else
                ParallelDescriptor::SendData(IoProc,
                                             &msg_hdr,
                                             TheBaseName.c_str(),
                                             //
                                             // Include NULL in MSG.
                                             //
                                             TheBaseName.length()+1);
#endif
            }
        }
    }
    break;
    default:
        BoxLib::Error("Bad case in switch");
    }
#ifdef BL_USE_MPI
    if (ParallelDescriptor::IOProcessor())
    {
        MPI_Status        status;
        const Array<int>& procmap = mf.DistributionMap().ProcessorMap();

        for (int idx = 0, N = procmap.length(); idx < N; idx++)
        {
            if (!(procmap[idx] == IoProc))
            {
                int rc, len, pos = 0;

                if ((rc = MPI_Recv(pbuf,
                                   NPBUF,
                                   MPI_PACKED,
                                   //
                                   // We recieve one message for each
                                   // index that IoProc doesn't own.
                                   //
                                   procmap[idx],
                                   idx,
                                   MPI_COMM_WORLD,
                                   &status)) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);
                //
                // Unpack FAB offset, length of name, and name.
                //
                if ((rc = MPI_Unpack(pbuf,
                                     NPBUF,
                                     &pos,
                                     &msg_hdr.m_head,
                                     1,
                                     MPI_LONG,
                                     MPI_COMM_WORLD)) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);

                if ((rc = MPI_Unpack(pbuf,
                                     NPBUF,
                                     &pos,
                                     &len,
                                     1,
                                     MPI_INT,
                                     MPI_COMM_WORLD)) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);

                char* fab_name = new char[len];

                if ((rc = MPI_Unpack(pbuf,
                                     NPBUF,
                                     &pos,
                                     fab_name,
                                     len,
                                     MPI_CHAR,
                                     MPI_COMM_WORLD)) != MPI_SUCCESS)
                    ParallelDescriptor::Abort(rc);
                
                hdr.m_fod[idx].m_head = msg_hdr.m_head;
                hdr.m_fod[idx].m_name = fab_name;

                delete [] fab_name;
            }
        }
    }
#else
    for (int len; ParallelDescriptor::GetMessageHeader(len, &msg_hdr); )
    {
        assert(ParallelDescriptor::IOProcessor());

        char* fab_name = new char[len];

        ParallelDescriptor::ReceiveData(fab_name, len);

        hdr.m_fod[msg_hdr.m_index].m_head = msg_hdr.m_head;
        hdr.m_fod[msg_hdr.m_index].m_name = fab_name;

        delete [] fab_name;
    }
#endif

    bytes += VisMF::WriteHeader(mf_name, hdr);

    return bytes;
}

VisMF::VisMF (const aString& mf_name)
    :
    m_mfname(mf_name)
{
    aString FullHdrFileName = m_mfname;

    FullHdrFileName += VisMF::MultiFabHdrFileSuffix;

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    ifstream ifs;

    ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.length());

    ifs.open(FullHdrFileName.c_str(), ios::in);

    if (!ifs.good())
        Utility::FileOpenFailed(FullHdrFileName);

    ifs >> m_hdr;

    m_pa.resize(m_hdr.m_ncomp);

    for (int nComp = 0; nComp < m_pa.length(); ++nComp)
    {
        m_pa[nComp].resize(m_hdr.m_ba.length());

        for (int ii = 0; ii < m_pa[nComp].length(); ++ii)
        {
            m_pa[nComp][ii] = 0;
        }
    }
}

FArrayBox*
VisMF::readFAB (int                  idx,
                const aString&       mf_name,
                const VisMF::Header& hdr,
		int                  ncomp)
{
    Box fab_box = hdr.m_ba[idx];

    if (hdr.m_ngrow)
        fab_box.grow(hdr.m_ngrow);

    FArrayBox* fab = new FArrayBox(fab_box, ncomp == -1 ? hdr.m_ncomp : 1);

    aString FullFileName = VisMF::DirName(mf_name);

    FullFileName += hdr.m_fod[idx].m_name;
    
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    ifstream ifs;

    ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.length());

    ifs.open(FullFileName.c_str(), ios::in);

    if (!ifs.good())
        Utility::FileOpenFailed(FullFileName);

    if (hdr.m_fod[idx].m_head)
        ifs.seekg(hdr.m_fod[idx].m_head, ios::beg);

    if (ncomp == -1)
    {
        fab->readFrom(ifs);
    }
    else
    {
        fab->readFrom(ifs, ncomp);
    }

    return fab;
}

void
VisMF::Read (MultiFab&      mf,
             const aString& mf_name)
{
    VisMF::Header hdr;

    aString FullHdrFileName = mf_name;

    FullHdrFileName += VisMF::MultiFabHdrFileSuffix;
    {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

        ifstream ifs;

        ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.length());

        ifs.open(FullHdrFileName.c_str(), ios::in);

        if (!ifs.good())
            Utility::FileOpenFailed(FullHdrFileName);

        ifs >> hdr;
    }
    mf.define(hdr.m_ba, hdr.m_ncomp, hdr.m_ngrow, Fab_noallocate);

    for (MultiFabIterator mfi(mf); mfi.isValid(); ++mfi)
    {
        mf.setFab(mfi.index(), VisMF::readFAB(mfi.index(), mf_name, hdr));
    }

    assert(mf.ok());
}

void
VisMF::clear (int fabIndex)
{
    for (int ncomp = 0; ncomp < m_pa.length(); ++ncomp)
    {
        clear(ncomp, fabIndex);
    }
}

void
VisMF::clear ()
{
  for (int ncomp = 0; ncomp < m_pa.length(); ++ncomp)
  {
      for (int fabIndex = 0; fabIndex < m_pa[ncomp].length(); ++fabIndex)
      {
          clear(ncomp, fabIndex);
      }
  }
}
