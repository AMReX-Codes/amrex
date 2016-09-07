// TODO: need to work on read for upc++

#include <winstd.H>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <deque>
//
// This MUST be defined if don't have pubsetbuf() in I/O Streams Library.
//
#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

#include <ccse-mpi.H>
#include <Utility.H>
#include <VisMF.H>
#include <ParmParse.H>

static const char* TheMultiFabHdrFileSuffix = "_H";

static const char* TheFabOnDiskPrefix = "FabOnDisk:";

int VisMF::verbose = 1;

//
// Set these in Initialize().
//
int VisMF::nOutFiles(64);
int VisMF::nMFFileInStreams(1);

namespace
{
    bool initialized = false;
}

void
VisMF::Initialize ()
{
    if (initialized) return;
    //
    // Use the same defaults as in Amr.cpp.
    //
    VisMF::SetNOutFiles(nOutFiles);

    VisMF::SetMFFileInStreams(nMFFileInStreams);

    BoxLib::ExecOnFinalize(VisMF::Finalize);

    ParmParse pp("vismf");
    pp.query("v",verbose);

    initialized = true;
}

void
VisMF::Finalize ()
{
    initialized = false;
}

void
VisMF::SetNOutFiles (int noutfiles)
{
    nOutFiles = std::max(1, std::min(ParallelDescriptor::NProcs(), noutfiles));
}

void
VisMF::SetMFFileInStreams (int nstreams)
{
    nMFFileInStreams = std::max(1, std::min(ParallelDescriptor::NProcs(), nstreams));
}

int
VisMF::GetNOutFiles()
{
    return nOutFiles;
}

std::ostream&
operator<< (std::ostream&           os,
            const VisMF::FabOnDisk& fod)
{
    os << TheFabOnDiskPrefix << ' ' << fod.m_name << ' ' << fod.m_head;

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

    BL_ASSERT(str == TheFabOnDiskPrefix);

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
#ifdef BL_USE_FLOAT
    double dtemp;
#endif

    is >> N >> ch >> M;

    if ( N < 0 ) 
      BoxLib::Error("Expected a positive integer, N, got something else");
    if ( M < 0 ) 
      BoxLib::Error("Expected a positive integer, M, got something else");
    if ( ch != ',' ) 
      BoxLib::Error("Expected a ',' got something else");

    ar.resize(N);
    
    for ( ; i < N; i++)
    {
        ar[i].resize(M);

        for (long j = 0; j < M; j++)
        {
#ifdef BL_USE_FLOAT
            is >> dtemp >> ch;
            ar[i][j] = static_cast<Real>(dtemp);
#else
            is >> ar[i][j] >> ch;
#endif
	    if ( ch != ',' ) 
	      BoxLib::Error("Expected a ',' got something else");
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
    // Force it to be written in scientific notation to match fParallel code.
    //
    std::ios::fmtflags oflags = os.flags();
    os.setf(std::ios::floatfield, std::ios::scientific);
    int old_prec = os.precision(15);

    os << hd.m_vers     << '\n';
    os << int(hd.m_how) << '\n';
    os << hd.m_ncomp    << '\n';
    os << hd.m_ngrow    << '\n';

    hd.m_ba.writeOn(os); os << '\n';

    os << hd.m_fod      << '\n';
    os << hd.m_min      << '\n';
    os << hd.m_max      << '\n';

    os.flags(oflags);
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
    case VisMF::NFiles:
        hd.m_how = VisMF::NFiles; break;
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
    m_pa[compIndex][fabIndex] = 0;
}

long
VisMF::FileOffset (std::ostream& os)
{
    //
    // Set to the end of the file before doing the tellp().
    // This shouldn't be needed except for a bug we've found
    // on edison.  As long as it doesn't hurt anything we'll
    // go with it for now.  The reason it should be OK is that
    // all our open()s for writing in this file are in append
    // mode.  So tellp() should always be at the file end.
    //
    os.seekp(0, std::ios::end);

    return os.tellp();
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

    if (const char* slash = strrchr(filename.c_str(), '/'))
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

    if (const char* slash = strrchr(str, '/'))
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

VisMF::Header::Header (const FabArray<FArrayBox>& mf,
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
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    //
    // Calculate m_min and m_max on the CPU owning the fab.
    //
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

    Array<int> nmtags(ParallelDescriptor::NProcs(),0);
    Array<int> offset(ParallelDescriptor::NProcs(),0);

    const Array<int>& pmap = mf.DistributionMap().ProcessorMap();

    for (int i = 0, N = mf.size(); i < N; i++)
        nmtags[pmap[i]]++;

    for (int i = 0, N = nmtags.size(); i < N; i++)
        //
        // Each Fab corresponds to 2*m_ncomp Reals.
        //
        nmtags[i] *= 2*m_ncomp;

    for (int i = 1, N = offset.size(); i < N; i++)
        offset[i] = offset[i-1] + nmtags[i-1];

    Array<Real> senddata(nmtags[ParallelDescriptor::MyProc()]);

    if (senddata.empty())
        //
        // Can't let senddata be empty as senddata.dataPtr() will fail.
        //
        senddata.resize(1);

    int ioffset = 0;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const int idx = mfi.index();

        for (int i = 0; i < m_ncomp; i++)
        {
            senddata[ioffset+i]         = m_min[idx][i];
            senddata[ioffset+m_ncomp+i] = m_max[idx][i];
        }

        ioffset += 2*m_ncomp;
    }

    BL_ASSERT(ioffset == nmtags[ParallelDescriptor::MyProc()]);

    Array<Real> recvdata(mf.size()*2*m_ncomp);

    BL_COMM_PROFILE(BLProfiler::Gatherv, recvdata.size() * sizeof(Real),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());

    BL_MPI_REQUIRE( MPI_Gatherv(senddata.dataPtr(),
                                nmtags[ParallelDescriptor::MyProc()],
                                ParallelDescriptor::Mpi_typemap<Real>::type(),
                                recvdata.dataPtr(),
                                nmtags.dataPtr(),
                                offset.dataPtr(),
                                ParallelDescriptor::Mpi_typemap<Real>::type(),
                                IOProc,
                                ParallelDescriptor::Communicator()) );

    BL_COMM_PROFILE(BLProfiler::Gatherv, recvdata.size() * sizeof(Real),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    if (ParallelDescriptor::IOProcessor())
    {
        for (int i = 0, N = mf.size(); i < N; i++)
        {
            if (pmap[i] != IOProc)
            {
                m_min[i].resize(m_ncomp);
                m_max[i].resize(m_ncomp);
            }
        }

        for (int j = 0, N = mf.size(); j < N; j++)
        {
            if (pmap[j] != IOProc)
            {
                for (int k = 0; k < m_ncomp; k++)
                {
                    m_min[j][k] = recvdata[offset[pmap[j]]+k];
                    m_max[j][k] = recvdata[offset[pmap[j]]+k+m_ncomp];
                }

                offset[pmap[j]] += 2*m_ncomp;
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

#ifdef BL_FIXHEADERDENORMS
    if (ParallelDescriptor::IOProcessor())
    {
        for(int i = 0; i < m_min.size(); ++i)
        {
            for(int j = 0; j < m_min[i].size(); ++j)
            {
                if (std::abs(m_min[i][j]) < 1.0e-300)
                {
                    m_min[i][j] = 0.0;
                }
            }
        }
    }
#endif
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

        MFHdrFileName += TheMultiFabHdrFileSuffix;

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

        MFHdrFile.close();
    }
    return bytes;
}

long
VisMF::Write (const FabArray<FArrayBox>&    mf,
              const std::string& mf_name,
              VisMF::How         how,
              bool               set_ghost)
{
    BL_ASSERT(mf_name[mf_name.length() - 1] != '/');

    static const char* FabFileSuffix = "_D_";

    VisMF::Initialize();

    VisMF::Header hdr(mf, how);

    if (set_ghost)
    {
        FabArray<FArrayBox>* the_mf = const_cast<FabArray<FArrayBox>*>(&mf);

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

    long        bytes    = 0;
    const int   MyProc   = ParallelDescriptor::MyProc();
    const int   NProcs   = ParallelDescriptor::NProcs();
    const int   NSets    = (NProcs + (nOutFiles - 1)) / nOutFiles;
    const int   MySet    = MyProc/nOutFiles;
    std::string FullName = BoxLib::Concatenate(mf_name + FabFileSuffix, MyProc % nOutFiles, 4);

    const std::string BName = VisMF::BaseName(FullName);

    for (int iSet = 0; iSet < NSets; ++iSet)
    {
        if (MySet == iSet)
        {
            {
                VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    
                std::ofstream FabFile;

                FabFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

                if (iSet == 0)
                {
                    //
                    // First set.
                    //
                    FabFile.open(FullName.c_str(),
                                 std::ios::out|std::ios::trunc|std::ios::binary);
                }
                else
                {
                    FabFile.open(FullName.c_str(),
                                 std::ios::out|std::ios::app|std::ios::binary);
                    //
                    // Set to the end of the file.
                    //
                    FabFile.seekp(0, std::ios::end);
                }
                if ( ! FabFile.good()) {
                    BoxLib::FileOpenFailed(FullName);
		}

                for (MFIter mfi(mf); mfi.isValid(); ++mfi)
                {
                    hdr.m_fod[mfi.index()] = VisMF::Write(mf[mfi],BName,FabFile,bytes);
                }

                FabFile.flush();

                FabFile.close();
            }

            int iBuff     = 0;
            int wakeUpPID = (MyProc + nOutFiles);
            int tag       = (MyProc % nOutFiles);
            if (wakeUpPID < NProcs) {
                ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
	    }
        }
        if (MySet == (iSet + 1))
        {
            //
            // Next set waits.
            //
            int iBuff;
            int waitForPID = (MyProc - nOutFiles);
            int tag        = (MyProc % nOutFiles);
            ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
        }
    }

#ifdef BL_USE_MPI
    ParallelDescriptor::Barrier("VisMF::Write");

    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Array<int> nmtags(NProcs,0);
    Array<int> offset(NProcs,0);

    const Array<int>& pmap = mf.DistributionMap().ProcessorMap();

    for (int i = 0, N = mf.size(); i < N; i++)
        nmtags[pmap[i]]++;

    for (int i = 1, N = offset.size(); i < N; i++)
        offset[i] = offset[i-1] + nmtags[i-1];

    Array<long> senddata(nmtags[ParallelDescriptor::MyProc()]);

    if (senddata.empty())
        //
        // Can't let senddata be empty as senddata.dataPtr() will fail.
        //
        senddata.resize(1);

    int ioffset = 0;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        senddata[ioffset++] = hdr.m_fod[mfi.index()].m_head;

    BL_ASSERT(ioffset == nmtags[ParallelDescriptor::MyProc()]);

    Array<long> recvdata(mf.size());

    BL_COMM_PROFILE(BLProfiler::Gatherv, recvdata.size() * sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());

    BL_MPI_REQUIRE( MPI_Gatherv(senddata.dataPtr(),
                                nmtags[ParallelDescriptor::MyProc()],
                                ParallelDescriptor::Mpi_typemap<long>::type(),
                                recvdata.dataPtr(),
                                nmtags.dataPtr(),
                                offset.dataPtr(),
                                ParallelDescriptor::Mpi_typemap<long>::type(),
                                IOProc,
                                ParallelDescriptor::Communicator()) );

    BL_COMM_PROFILE(BLProfiler::Gatherv, recvdata.size() * sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    if (ParallelDescriptor::IOProcessor())
    {
        Array<int> cnt(NProcs,0);

        for (int j = 0, N = mf.size(); j < N; ++j)
        {
            const int i = pmap[j];

            hdr.m_fod[j].m_head = recvdata[offset[i]+cnt[i]];

            std::string name = BoxLib::Concatenate(mf_name + FabFileSuffix, i % nOutFiles, 4);

            hdr.m_fod[j].m_name = VisMF::BaseName(name);

            cnt[i]++;
        }
    }
#endif /*BL_USE_MPI*/

    bytes += VisMF::WriteHeader(mf_name, hdr);

    return bytes;
}

VisMF::VisMF (const std::string& mf_name)
    :
    m_mfname(mf_name)
{
    std::string FullHdrFileName = m_mfname;

    FullHdrFileName += TheMultiFabHdrFileSuffix;

    Array<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(FullHdrFileName, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream ifs(fileCharPtrString, std::istringstream::in);

    ifs >> m_hdr;

    m_pa.resize(m_hdr.m_ncomp);

    for (int nComp = 0; nComp < m_pa.size(); ++nComp)
    {
        m_pa[nComp].resize(m_hdr.m_ba.size());

        for (int ii = 0, N = m_pa[nComp].size(); ii < N; ++ii)
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

    std::string FullName = VisMF::DirName(mf_name);

    FullName += hdr.m_fod[idx].m_name;
    
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ifstream ifs;

    ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    ifs.open(FullName.c_str(), std::ios::in|std::ios::binary);

    if (!ifs.good())
        BoxLib::FileOpenFailed(FullName);

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

    ifs.close();

    return fab;
}

void
VisMF::readFAB (FabArray<FArrayBox>&            mf,
		int                  idx,
                const std::string&   mf_name,
                const VisMF::Header& hdr)
{
    FArrayBox& fab = mf[idx];

    std::string FullName = VisMF::DirName(mf_name);

    FullName += hdr.m_fod[idx].m_name;
    
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ifstream ifs;

    ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    ifs.open(FullName.c_str(), std::ios::in|std::ios::binary);

    if (!ifs.good())
        BoxLib::FileOpenFailed(FullName);

    if (hdr.m_fod[idx].m_head)
        ifs.seekg(hdr.m_fod[idx].m_head, std::ios::beg);

    fab.readFrom(ifs);

    ifs.close();
}

void
VisMF::Read (FabArray<FArrayBox>&          mf,
             const std::string& mf_name)
{
  BL_PROFILE("VisMF::Read()");

  VisMF::Initialize();

  if (verbose && ParallelDescriptor::IOProcessor())
  {
      std::cout << "VisMF::Read:  about to read:  " << mf_name << std::endl;
  }

#ifdef BL_VISMF_MSGCHECK
  ParallelDescriptor::Barrier("VisMF::Read::MSGCHECK_0");
  {
      MPI_Status mwstatus;
      int mwflag(0);
      bool bextra(false);
      ParallelDescriptor::IProbe(MPI_ANY_SOURCE, MPI_ANY_TAG, mwflag, mwstatus);
      while(mwflag) {
        bextra = true;
        std::cout << "### ### ### :  EXTRA MESSAGE BEFORE:  myproc = "
	          << ParallelDescriptor::MyProc() << '\n' << std::flush;
        ParallelDescriptor::Message rmess;
        std::vector<int> cread(mwstatus.count);
        rmess = ParallelDescriptor::Recv(cread, mwstatus.MPI_SOURCE, mwstatus.MPI_TAG);
        ParallelDescriptor::IProbe(MPI_ANY_SOURCE, MPI_ANY_TAG, mwflag, mwstatus);
      }
      if (bextra) {
        //BoxLib::Abort("EXTRA MESSAGES BEFORE");
      }
  }
  ParallelDescriptor::Barrier("VisMF::Read::MSGCHECK_1");
#endif

    VisMF::Header hdr;

    std::string FullHdrFileName = mf_name;

    FullHdrFileName += TheMultiFabHdrFileSuffix;

    Real hEndTime, hStartTime;

    {
        hStartTime = ParallelDescriptor::second();
        Array<char> fileCharPtr;
        ParallelDescriptor::ReadAndBcastFile(FullHdrFileName, fileCharPtr);
        std::string fileCharPtrString(fileCharPtr.dataPtr());
        std::istringstream ifs(fileCharPtrString, std::istringstream::in);
        hEndTime = ParallelDescriptor::second();

        ifs >> hdr;
    }
    mf.define(hdr.m_ba, hdr.m_ncomp, hdr.m_ngrow, Fab_allocate);

#ifdef BL_USE_MPI
    //
    // Here we limit the number of open files when reading a multifab.
    //
    Real startTime(ParallelDescriptor::second());
    static Real totalTime(0.0);
    int nReqs(0), ioProcNum(ParallelDescriptor::IOProcessorNumber());
    int myProc(ParallelDescriptor::MyProc());
    int nBoxes(hdr.m_ba.size());
    int totalIOReqs(nBoxes), nFiles(-1);
    std::vector<int> iDone(2);
    const int iDoneIndex(0), iDoneCount(1);
    std::set<int> busyProcs;  // [whichProc]
    std::map<std::string, int> fileNames;  // <filename, allreadsindex>
    std::multiset<int> availableFiles;  // [whichFile]  supports multiple reads/file
    int nOpensPerFile(nMFFileInStreams), allReadsIndex(0), messTotal(0);
    ParallelDescriptor::Message rmess;
    Array<std::map<int,std::map<int,int> > > allReads; // [file]<proc,<seek,index>>

    for(int i(0); i < nBoxes; ++i) {   // count the files
      int whichProc(mf.DistributionMap()[i]);
      if(whichProc == myProc) {
        ++nReqs;
      }
      if(ParallelDescriptor::IOProcessor()) {
        std::string fname(hdr.m_fod[i].m_name);
	if(fileNames.insert(std::pair<std::string,int>(fname,allReadsIndex)).second)
	{
	  ++allReadsIndex;
	}
      }
    }

    if(ParallelDescriptor::IOProcessor()) {    // fill availableFiles
      nFiles = fileNames.size();
      for(int i(0); i < nFiles; ++i) {
        for(int nOpens(0); nOpens < nOpensPerFile; ++nOpens) {
          availableFiles.insert(i);
        }
      }
      allReads.resize(nFiles);
      int whichProc, iSeekPos;
      std::map<std::string, int>::iterator fileNamesIter;
      for(int i(0); i < nBoxes; ++i) {   // fill allReads maps
        whichProc = mf.DistributionMap()[i];
        iSeekPos = hdr.m_fod[i].m_head;
        std::string fname(hdr.m_fod[i].m_name);
	fileNamesIter = fileNames.find(fname);
	if(fileNamesIter != fileNames.end()) {
	  int findex(fileNames.find(fname)->second);
	  allReads[findex][whichProc].insert(std::pair<int, int>(iSeekPos, i));
	} else {
	  std::cout << "**** Error:  filename not found = " << fname << std::endl;
	  BoxLib::Abort();
	}
      }
    }

    if(ParallelDescriptor::IOProcessor()) {  // manage the file locks
      int reqsPending(0), iopFileIndex;
      std::deque<int> iopReads;
      MPI_Status status;
      int doneFlag;
      while(totalIOReqs > 0) {
	std::vector<int> vReads;
        std::multiset<int>::iterator aFilesIter;
        aFilesIter = availableFiles.begin();
        while(aFilesIter != availableFiles.end()) {  // handle available files
	  int arIndex(*aFilesIter);
	  if(allReads[arIndex].empty()) {
            availableFiles.erase(arIndex);
            aFilesIter = availableFiles.begin();
	    continue;
	  }
          std::map<int,std::map<int,int> >::iterator whichRead;
	  for(whichRead = allReads[arIndex].begin();
	      whichRead != allReads[arIndex].end(); ++whichRead)
	  {
	    int tryProc(whichRead->first);
	    if(busyProcs.find(tryProc) == busyProcs.end()) {  // tryProc not busy
	      busyProcs.insert(tryProc);
	      int nReads(whichRead->second.size());
	      int ir(0);
	      vReads.resize(nReads);
              std::map<int,int>::iterator imiter;
	      for(imiter = whichRead->second.begin();
	          imiter != whichRead->second.end(); ++imiter)
	      {
	        vReads[ir] = imiter->second;  // the mfindex
		++ir;
	      }
	      if(tryProc == ioProcNum) {
		iopFileIndex = arIndex;
		for(int irr(0); irr < nReads; ++irr) {
	          iopReads.push_back(vReads[irr]);
		}
	      } else {
	        ParallelDescriptor::Send(vReads, tryProc, vReads[0]);
		++messTotal;
		++reqsPending;
	      }
              availableFiles.erase(aFilesIter);
              aFilesIter = availableFiles.begin();
	      break;
	    }
	  }
	  if(whichRead == allReads[arIndex].end()) {
	    ++aFilesIter;
	  } else {
	    allReads[arIndex].erase(whichRead);
	  }
        }  // end while(aFilesIter...)

	while( ! iopReads.empty()) {
	  int index = iopReads.front();
	  VisMF::readFAB(mf,index, mf_name, hdr);
	  --totalIOReqs;
	  iopReads.pop_front();
	  if(iopReads.empty()) {
	    availableFiles.insert(iopFileIndex);
	    busyProcs.erase(ioProcNum);
	  }
	  ParallelDescriptor::IProbe(MPI_ANY_SOURCE, MPI_ANY_TAG,
	                             doneFlag, status);
	  if(doneFlag) {
	    break;
	  }
	}

	if(reqsPending > 0) {
          rmess = ParallelDescriptor::Recv(iDone, MPI_ANY_SOURCE, MPI_ANY_TAG);

	  int index(iDone[iDoneIndex]);
	  int procDone(rmess.pid());
	  totalIOReqs -= iDone[iDoneCount];
	  --reqsPending;
	  busyProcs.erase(procDone);
          std::string fname(hdr.m_fod[index].m_name);
	  int fileIndex(fileNames.find(fname)->second);
	  availableFiles.insert(fileIndex);
	}

      }  // end while(totalIOReqs...)

    } else {  /// all other procs
      std::vector<int> recReads(nReqs);
      while(nReqs > 0) {
        rmess = ParallelDescriptor::Recv(recReads, ioProcNum, MPI_ANY_TAG);
        for(int ir(0); ir < rmess.count(); ++ir) {
	  int mfIndex(recReads[ir]);
	  VisMF::readFAB(mf,mfIndex, mf_name, hdr);
	}
        nReqs -= rmess.count();
	iDone[iDoneIndex] = recReads[0];
	iDone[iDoneCount] = rmess.count();
        ParallelDescriptor::Send(iDone, ioProcNum, recReads[0]);
      }
    }


#ifdef BL_VISMF_MSGCHECK
  ParallelDescriptor::Barrier("VisMF::Read::MSGCHECK_2");
  {
      MPI_Status mwstatus;
      int mwflag(0);
      bool bextra(false);
      ParallelDescriptor::IProbe(MPI_ANY_SOURCE, MPI_ANY_TAG, mwflag, mwstatus);
      while(mwflag) {
        bextra = true;
        std::cout << "### ### ### :  EXTRA MESSAGE AFTER:  myproc = "
	          << ParallelDescriptor::MyProc() << '\n' << std::flush;
        ParallelDescriptor::Message rmess;
        std::vector<int> cread(mwstatus.count);
        rmess = ParallelDescriptor::Recv(cread, mwstatus.MPI_SOURCE, mwstatus.MPI_TAG);
        ParallelDescriptor::IProbe(MPI_ANY_SOURCE, MPI_ANY_TAG, mwflag, mwstatus);
      }
      if(bextra) {
        //BoxLib::Abort("EXTRA MESSAGES AFTER");
      }
  }
#endif
    ParallelDescriptor::Barrier("VisMF::Read");

    if (ParallelDescriptor::IOProcessor() && false)
    {
      Real mfReadTime = ParallelDescriptor::second() - startTime;
      totalTime += mfReadTime;
      std::cout << "MFRead:::  nBoxes = "
                << nBoxes
                << "  nMessages = "
                << messTotal << '\n';
      std::cout << "MFRead:::  hTime = " << (hEndTime - hStartTime) << '\n';
      std::cout << "MFRead:::  mfReadTime = "
                << mfReadTime
                << "  totalTime = "
                << totalTime << std::endl;
    }
#else
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
	VisMF::readFAB(mf,mfi.index(), mf_name, hdr);
    }
#endif

    BL_ASSERT(mf.ok());
}



void
VisMF::Check (const std::string& mf_name)
{
  BL_PROFILE("VisMF::Check()");

  VisMF::Initialize();

  if(ParallelDescriptor::IOProcessor()) {
      std::cout << "VisMF::Check:  about to check:  " << mf_name << std::endl;

    char c;
    int nBadFabs(0);
    VisMF::Header hdr;
    std::string FullHdrFileName(mf_name);
    FullHdrFileName += TheMultiFabHdrFileSuffix;

    std::ifstream ifs(FullHdrFileName.c_str());

    ifs >> hdr;

    ifs.close();

    std::cout << "hdr.boxarray size =  " << hdr.m_ba.size() << std::endl;
    std::cout << "mf.ncomp =  " << hdr.m_ncomp << std::endl;
    std::cout << "number of fabs on disk =  " << hdr.m_fod.size() << std::endl;
    std::cout << "DirName = " << DirName(mf_name) << std::endl;
    std::cout << "mf_name = " << mf_name << std::endl;
    std::cout << "FullHdrFileName = " << FullHdrFileName << std::endl;

    // check that the string FAB is where it should be
    for(int i(0); i < hdr.m_fod.size(); ++i) {
      bool badFab(false);
      FabOnDisk &fod = hdr.m_fod[i];
      std::string FullName(VisMF::DirName(mf_name));
      FullName += fod.m_name;
      std::ifstream ifs;
      ifs.open(FullName.c_str(), std::ios::in|std::ios::binary);

      if( ! ifs.good()) {
        std::cout << "**** Error:  could not open file:  " << FullName << std::endl;
	continue;
      }

      ifs.seekg(fod.m_head, std::ios::beg);

      ifs >> c;
      if(c != 'F') {
        badFab = true;
      }
      ifs >> c;
      if(c != 'A') {
        badFab = true;
      }
      ifs >> c;
      if(c != 'B') {
        badFab = true;
      }
      if(badFab) {
	++nBadFabs;
        std::cout << "**** Error in file:  " << FullName << "  Bad Fab at index = "
	          << i << "   at seekpos = " << fod.m_head << std::endl;
      }
      ifs.close();

    }
    if(nBadFabs) {
      std::cout << "Total Bad Fabs = " << nBadFabs << std::endl;
    } else {
      std::cout << "No Bad Fabs." << std::endl;
    }

  }

}




void
VisMF::clear (int fabIndex)
{
    for (int ncomp = 0, N = m_pa.size(); ncomp < N; ++ncomp)
    {
        clear(ncomp, fabIndex);
    }
}

void
VisMF::clear ()
{
    for (int ncomp = 0, N = m_pa.size(); ncomp < N; ++ncomp)
    {
        for (int fabIndex = 0, M = m_pa[ncomp].size(); fabIndex < M; ++fabIndex)
        {
            clear(ncomp, fabIndex);
        }
    }
}
