//
// $Id: VisMF.cpp,v 1.92 2002-04-17 22:22:19 vince Exp $
//

#include <winstd.H>
#include <cstdio>
#include <fstream>

//
// This MUST be defined if don't have pubsetbuf() in I/O Streams Library.
//
#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

#include <ccse-mpi.H>
#include <Utility.H>
#include <VisMF.H>

const std::string VisMF::FabFileSuffix("_D_");
const std::string VisMF::MultiFabHdrFileSuffix("_H");
const std::string VisMF::FabOnDisk::Prefix("FabOnDisk:");

std::ostream&
operator<< (std::ostream&           os,
            const VisMF::FabOnDisk& fod)
{
    os << VisMF::FabOnDisk::Prefix << ' ' << fod.m_name << ' ' << fod.m_head;

    if (!os.good())
        BoxLib::Error("Write of VisMF::FabOnDisk failed");

    return os;
}

std::istream&
operator>> (std::istream&     is,
            VisMF::FabOnDisk& fod)
{
    std::string str;
    is >> str;

    BL_ASSERT(str == VisMF::FabOnDisk::Prefix);

    is >> fod.m_name;
    is >> fod.m_head;

    if (!is.good())
        BoxLib::Error("Read of VisMF::FabOnDisk failed");

    return is;
}

std::ostream&
operator<< (std::ostream&                  os,
            const Array<VisMF::FabOnDisk>& fa)
{
    long i = 0, N = fa.size();

    os << N << '\n';

    for ( ; i < N; i++)
    {
        os << fa[i] << '\n';
    }

    if (!os.good())
        BoxLib::Error("Write of Array<VisMF::FabOnDisk> failed");

    return os;
}

std::istream&
operator>> (std::istream&            is,
            Array<VisMF::FabOnDisk>& fa)
{
    long i = 0, N;

    is >> N;
    BL_ASSERT(N >= 0);

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
std::ostream&
operator<< (std::ostream&               os,
            const Array< Array<Real> >& ar)
{
    long i = 0, N = ar.size(), M = (N == 0) ? 0 : ar[0].size();

    os << N << ',' << M << '\n';

    for ( ; i < N; i++)
    {
        BL_ASSERT(ar[i].size() == M);

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
std::istream&
operator>> (std::istream&         is,
            Array< Array<Real> >& ar)
{
    char ch;
    long i = 0, N, M;

    is >> N >> ch >> M;

    BL_ASSERT(N >= 0);
    BL_ASSERT(ch == ',');
    BL_ASSERT(M >= 0);

    ar.resize(N);
    
    for ( ; i < N; i++)
    {
        ar[i].resize(M);

        for (long j = 0; j < M; j++)
        {
            is >> ar[i][j] >> ch;
            BL_ASSERT(ch == ',');
        }
    }

    if (!is.good())
        BoxLib::Error("Read of Array<Array<Real>> failed");

    return is;
}

std::ostream&
operator<< (std::ostream&        os,
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

std::istream&
operator>> (std::istream&  is,
            VisMF::Header& hd)
{
    is >> hd.m_vers;
    BL_ASSERT(hd.m_vers == VisMF::Header::Version);

    int how;
    is >> how;
    switch (how)
    {
    case VisMF::OneFilePerCPU:
        hd.m_how = VisMF::OneFilePerCPU; break;
    default:
        BoxLib::Error("Bad case in switch");
    }

    is >> hd.m_ncomp;
    BL_ASSERT(hd.m_ncomp >= 0);

    is >> hd.m_ngrow;
    BL_ASSERT(hd.m_ngrow >= 0);

    hd.m_ba.readFrom(is);

    is >> hd.m_fod;
    BL_ASSERT(hd.m_ba.size() == hd.m_fod.size());

    is >> hd.m_min;
    is >> hd.m_max;

    BL_ASSERT(hd.m_ba.size() == hd.m_min.size());
    BL_ASSERT(hd.m_ba.size() == hd.m_max.size());

    if (!is.good())
        BoxLib::Error("Read of VisMF::Header failed");

    return is;
}

VisMF::FabOnDisk::FabOnDisk () {}

VisMF::FabOnDisk::FabOnDisk (const std::string& name, long offset)
    :
    m_name(name),
    m_head(offset)
{}

int
VisMF::nComp () const
{
    return m_hdr.m_ncomp;
}

int
VisMF::nGrow () const
{
    return m_hdr.m_ngrow;
}

int
VisMF::size () const
{
    return m_hdr.m_ba.size();
}

const BoxArray&
VisMF::boxArray () const
{
    return m_hdr.m_ba;
}

Real
VisMF::min (int fabIndex,
            int nComp) const
{
    BL_ASSERT(0 <= fabIndex && fabIndex < m_hdr.m_ba.size());
    BL_ASSERT(0 <= nComp && nComp < m_hdr.m_ncomp);
    return m_hdr.m_min[fabIndex][nComp];
}

Real
VisMF::max (int fabIndex,
            int nComp) const
{
    BL_ASSERT(0 <= fabIndex && fabIndex < m_hdr.m_ba.size());
    BL_ASSERT(0 <= nComp && nComp < m_hdr.m_ncomp);
    return m_hdr.m_max[fabIndex][nComp];
}

const FArrayBox&
VisMF::GetFab (int fabIndex,
               int ncomp) const
{
    if (m_pa[ncomp][fabIndex] == 0)
    {
        m_pa[ncomp][fabIndex] = VisMF::readFAB(fabIndex,m_mfname,m_hdr,ncomp);
    }
    return *m_pa[ncomp][fabIndex];
}

void
VisMF::clear (int fabIndex,
              int compIndex)
{
    delete m_pa[compIndex][fabIndex];
}

long
VisMF::FileOffset (std::ostream& os)
{
    return
#if defined(__KCC) 
#if ((BL_KCC_MAJOR_VERSION >= 4) || (BL_KCC_MAJOR_VERSION == 3 && BL_KCC_MINOR_VERSION > 3) || (BL_KCC_MAJOR_VERSION == 0))
      os.tellp();
#else
    os.tellp().offset();
#endif
#else
    os.tellp();
#endif
}

FArrayBox*
VisMF::readFAB (int                idx,
                const std::string& mf_name)
{
    return VisMF::readFAB(idx,mf_name,m_hdr,-1);
}

FArrayBox*
VisMF::readFAB (int idx,
		int ncomp)
{
    return VisMF::readFAB(idx,m_mfname,m_hdr,ncomp);
}

std::string
VisMF::BaseName (const std::string& filename)
{
    BL_ASSERT(filename[filename.length() - 1] != '/');

    if (char* slash = strrchr(filename.c_str(), '/'))
    {
        //
        // Got at least one slash -- give'm the following tail.
        //
        return std::string(slash + 1);
    }
    else
    {
        //
        // No leading directory portion to name.
        //
        return filename;
    }
}

std::string
VisMF::DirName (const std::string& filename)
{
    BL_ASSERT(filename[filename.length() - 1] != '/');

    static const std::string TheNullString("");

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

        std::string dirname = buf;

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
VisMF::Write (const FArrayBox&   fab,
              const std::string& filename,
              std::ostream&      os,
              long&              bytes)
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
    m_fod(m_ba.size()),
    m_min(m_ba.size()),
    m_max(m_ba.size())
{
#ifdef BL_USE_MPI
    //
    // Note that m_min and m_max are only calculated on CPU owning the fab.
    // We pass this data back to IOProcessor() so it sees the whole Header.
    //
    const int SeqNo  = ParallelDescriptor::SeqNum();
    const int NProcs = ParallelDescriptor::NProcs();
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    int nFabs = 0;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        nFabs++;

        const int idx = mfi.index();

        m_min[idx].resize(m_ncomp);
        m_max[idx].resize(m_ncomp);

        BL_ASSERT(mf[mfi].box().contains(m_ba[idx]));

        for (long j = 0; j < m_ncomp; j++)
        {
            m_min[idx][j] = mf[mfi].min(m_ba[idx],j);
            m_max[idx][j] = mf[mfi].max(m_ba[idx],j);
        }
    }

    if (!ParallelDescriptor::IOProcessor())
    {
        if (nFabs)
        {
            Array<Real> senddata(2*m_ncomp*nFabs);

            int offset = 0;

            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                const int idx = mfi.index();

                for (int i = 0; i < m_ncomp; i++)
                {
                    senddata[offset+i]         = m_min[idx][i];
                    senddata[offset+m_ncomp+i] = m_max[idx][i];
                }

                offset += 2*m_ncomp;
            }

            BL_ASSERT(offset == 2*m_ncomp*nFabs);

            BL_MPI_REQUIRE( MPI_Send(senddata.dataPtr(),
                                     2*m_ncomp*nFabs,
                                     ParallelDescriptor::Mpi_typemap<Real>::type(),
                                     IOProc,
                                     SeqNo,
                                     ParallelDescriptor::Communicator()) );

            BL_ASSERT(offset == 2*m_ncomp*nFabs);
        }
    }
    else
    {
        const Array<int>& procmap = mf.DistributionMap().ProcessorMap();

        Array<int>           fabs(NProcs,0);
        Array<int>           indx(NProcs);
        Array<MPI_Request>   reqs(NProcs,MPI_REQUEST_NULL);
        Array<MPI_Status>    status(NProcs);
        Array< Array<Real> > data(NProcs);

        for (int i = 0, N = procmap.size(); i < N; i++)
            fabs[procmap[i]]++;

        fabs[IOProc] = 0;

        int NWaits = 0;

        for (int i = 0; i < NProcs; i++)
        {
            if (fabs[i])
            {
                NWaits++;

                data[i].resize(2*m_ncomp*fabs[i]);

                BL_MPI_REQUIRE( MPI_Irecv(data[i].dataPtr(),
                                          2*m_ncomp*fabs[i],
                                          ParallelDescriptor::Mpi_typemap<Real>::type(),
                                          i,
                                          SeqNo,
                                          ParallelDescriptor::Communicator(),
                                          &reqs[i]) );
            }
        }

        for (int completed; NWaits > 0; NWaits -= completed)
        {
            BL_MPI_REQUIRE( MPI_Waitsome(NProcs,
                                         reqs.dataPtr(),
                                         &completed,
                                         indx.dataPtr(),
                                         status.dataPtr()) );

            for (int k = 0; k < completed; k++)
            {
                int Ncpu = indx[k], offset = 0;

                for (int idx = 0, N = procmap.size(); idx < N; idx++)
                {
                    if (procmap[idx] == Ncpu)
                    {
                        m_min[idx].resize(m_ncomp);
                        m_max[idx].resize(m_ncomp);

                        for (int i = 0; i < m_ncomp; i++)
                        {
                            m_min[idx][i] = data[Ncpu][offset+i];
                            m_max[idx][i] = data[Ncpu][offset+m_ncomp+i];
                        }

                        offset += 2*m_ncomp;
                    }
                }

                BL_ASSERT(offset == 2*m_ncomp*fabs[Ncpu]);
            }
        }
    }
#else
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const int idx = mfi.index();

        m_min[idx].resize(m_ncomp);
        m_max[idx].resize(m_ncomp);

        BL_ASSERT(mf[mfi].box().contains(m_ba[idx]));

        for (long j = 0; j < m_ncomp; j++)
        {
            m_min[idx][j] = mf[mfi].min(m_ba[idx],j);
            m_max[idx][j] = mf[mfi].max(m_ba[idx],j);
        }
    }
#endif /*BL_USE_MPI*/
}

long
VisMF::WriteHeader (const std::string& mf_name,
                    VisMF::Header&     hdr)
{
    long bytes = 0;
    //
    // When running in parallel only one processor should do this I/O.
    //
    if (ParallelDescriptor::IOProcessor())
    {
        std::string MFHdrFileName = mf_name;

        MFHdrFileName += VisMF::MultiFabHdrFileSuffix;

        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

        std::ofstream MFHdrFile;

        MFHdrFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        MFHdrFile.open(MFHdrFileName.c_str(), std::ios::out|std::ios::trunc);

        if (!MFHdrFile.good())
            BoxLib::FileOpenFailed(MFHdrFileName);

        MFHdrFile << hdr;
        //
        // Add in the number of bytes written out in the Header.
        //
        bytes += VisMF::FileOffset(MFHdrFile);
    }
    return bytes;
}

long
VisMF::Write (const MultiFab&    mf,
              const std::string& mf_name,
              VisMF::How         how,
              bool               set_ghost)
{
    BL_ASSERT(mf_name[mf_name.length() - 1] != '/');

    const int MyProc = ParallelDescriptor::MyProc();

    long bytes = 0;

    VisMF::Header hdr(mf, how);

    if (set_ghost)
    {
        MultiFab* the_mf = const_cast<MultiFab*>(&mf);

        BL_ASSERT(!(the_mf == 0));
        BL_ASSERT(hdr.m_ba == mf.boxArray());
        BL_ASSERT(hdr.m_ncomp == mf.nComp());

        for (MFIter mfi(*the_mf); mfi.isValid(); ++mfi)
        {
            const int idx = mfi.index();

            for (int j = 0; j < hdr.m_ncomp; j++)
            {
                const Real val = (hdr.m_min[idx][j] + hdr.m_max[idx][j]) / 2;

                the_mf->get(mfi).setComplement(val, hdr.m_ba[idx], j, 1);
            }
        }
    }

    char buf[sizeof(int) + 1];

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::string FullFileName = mf_name;

    FullFileName += VisMF::FabFileSuffix;
    sprintf(buf, "%04d", MyProc);
    FullFileName += buf;

    std::ofstream FabFile;

    FabFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    FabFile.open(FullFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);

    if (!FabFile.good())
        BoxLib::FileOpenFailed(FullFileName);

    std::string basename = VisMF::BaseName(FullFileName);

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        hdr.m_fod[mfi.index()] = VisMF::Write(mf[mfi],basename,FabFile,bytes);

#ifdef BL_USE_MPI
    const int SeqNo  = ParallelDescriptor::SeqNum();
    const int NProcs = ParallelDescriptor::NProcs();
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    if (!ParallelDescriptor::IOProcessor())
    {
        int nFabs = 0, idx = 0;

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            nFabs++;

        if (nFabs)
        {
            Array<long> senddata(nFabs);

            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
                senddata[idx++] = hdr.m_fod[mfi.index()].m_head;

            BL_MPI_REQUIRE( MPI_Send(senddata.dataPtr(),
                                     nFabs,
                                     MPI_LONG,
                                     IOProc,
                                     SeqNo,
                                     ParallelDescriptor::Communicator()));
        }

        BL_ASSERT(idx == nFabs);
    }
    else
    {
        const Array<int>& procmap = mf.DistributionMap().ProcessorMap();

        Array<int>           fabs(NProcs,0);
        Array<int>           indx(NProcs);
        Array<MPI_Request>   reqs(NProcs,MPI_REQUEST_NULL);
        Array<MPI_Status>    status(NProcs);
        Array< Array<long> > data(NProcs);

        for (int i = 0, N = procmap.size(); i < N; i++)
            fabs[procmap[i]]++;

        fabs[IOProc] = 0;

        int NWaits = 0;

        for (int i = 0; i < NProcs; i++)
        {
            if (fabs[i])
            {
                NWaits++;

                data[i].resize(fabs[i]);

                BL_MPI_REQUIRE( MPI_Irecv(data[i].dataPtr(),
                                          fabs[i],
                                          MPI_LONG,
                                          i,
                                          SeqNo,
                                          ParallelDescriptor::Communicator(),
                                          &reqs[i]));
            }
        }

        for (int completed; NWaits > 0; NWaits -= completed)
        {
            BL_MPI_REQUIRE( MPI_Waitsome(NProcs,
                                         reqs.dataPtr(),
                                         &completed,
                                         indx.dataPtr(),
                                         status.dataPtr()));

            for (int k = 0; k < completed; k++)
            {
                int Ncpu = indx[k], offset = 0;

                for (int idx = 0, N = procmap.size(); idx < N; idx++)
                {
                    if (procmap[idx] == Ncpu)
                    {
                        hdr.m_fod[idx].m_head = data[Ncpu][offset++];

                        std::string name = mf_name;

                        name += VisMF::FabFileSuffix;
                        sprintf(buf, "%04d", Ncpu);
                        name += buf;

                        hdr.m_fod[idx].m_name = VisMF::BaseName(name);
                    }
                }

                BL_ASSERT(offset == fabs[Ncpu]);
            }
        }
    }
#endif /*BL_USE_MPI*/

    if (VisMF::FileOffset(FabFile) <= 0)
        BoxLib::UnlinkFile(FullFileName);

    bytes += VisMF::WriteHeader(mf_name, hdr);

    return bytes;
}

VisMF::VisMF (const std::string& mf_name)
    :
    m_mfname(mf_name)
{
    std::string FullHdrFileName = m_mfname;

    FullHdrFileName += VisMF::MultiFabHdrFileSuffix;

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ifstream ifs;

    ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    ifs.open(FullHdrFileName.c_str(), std::ios::in);

    if (!ifs.good())
        BoxLib::FileOpenFailed(FullHdrFileName);

    ifs >> m_hdr;

    m_pa.resize(m_hdr.m_ncomp);

    for (int nComp = 0; nComp < m_pa.size(); ++nComp)
    {
        m_pa[nComp].resize(m_hdr.m_ba.size());

        for (int ii = 0; ii < m_pa[nComp].size(); ++ii)
        {
            m_pa[nComp][ii] = 0;
        }
    }
}

FArrayBox*
VisMF::readFAB (int                  idx,
                const std::string&   mf_name,
                const VisMF::Header& hdr,
		int                  ncomp)
{
    Box fab_box = hdr.m_ba[idx];

    if (hdr.m_ngrow)
        fab_box.grow(hdr.m_ngrow);

    FArrayBox* fab = new FArrayBox(fab_box, ncomp == -1 ? hdr.m_ncomp : 1);

    std::string FullFileName = VisMF::DirName(mf_name);

    FullFileName += hdr.m_fod[idx].m_name;
    
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ifstream ifs;

    ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    ifs.open(FullFileName.c_str(), std::ios::in|std::ios::binary);

    if (!ifs.good())
        BoxLib::FileOpenFailed(FullFileName);

    if (hdr.m_fod[idx].m_head)
        ifs.seekg(hdr.m_fod[idx].m_head, std::ios::beg);

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
VisMF::Read (MultiFab&          mf,
             const std::string& mf_name)
{
    VisMF::Header hdr;

    std::string FullHdrFileName = mf_name;

    FullHdrFileName += VisMF::MultiFabHdrFileSuffix;
    {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

        std::ifstream ifs;

        ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        ifs.open(FullHdrFileName.c_str(), std::ios::in);

        if (!ifs.good())
            BoxLib::FileOpenFailed(FullHdrFileName);

        ifs >> hdr;
    }
    mf.define(hdr.m_ba, hdr.m_ncomp, hdr.m_ngrow, Fab_noallocate);

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        mf.setFab(mfi.index(), VisMF::readFAB(mfi.index(), mf_name, hdr));
    }

    BL_ASSERT(mf.ok());
}

void
VisMF::clear (int fabIndex)
{
    for (int ncomp = 0; ncomp < m_pa.size(); ++ncomp)
    {
        clear(ncomp, fabIndex);
    }
}

void
VisMF::clear ()
{
    for (int ncomp = 0; ncomp < m_pa.size(); ++ncomp)
    {
        for (int fabIndex = 0; fabIndex < m_pa[ncomp].size(); ++fabIndex)
        {
            clear(ncomp, fabIndex);
        }
    }
}

