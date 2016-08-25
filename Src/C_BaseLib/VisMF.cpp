// TODO: need to work on read for upc++

#include <winstd.H>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <deque>
#include <cerrno>
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
#include <NFiles.H>
#include <FPC.H>
#include <FabConv.H>

static const char *TheMultiFabHdrFileSuffix = "_H";
static const char *FabFileSuffix = "_D_";
static const char *TheFabOnDiskPrefix = "FabOnDisk:";

int VisMF::verbose(1);
VisMF::Header::Version VisMF::currentVersion(VisMF::Header::Version_v1);
bool VisMF::groupSets(false);
bool VisMF::setBuf(true);
bool VisMF::useSingleRead(false);
bool VisMF::useSingleWrite(false);

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
    BL_PROFILE("VisMF::Initialize");

    if(initialized) {
      return;
    }
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

    if( ! os.good()) {
        BoxLib::Error("Write of VisMF::FabOnDisk failed");
    }

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

    if( ! is.good()) {
        BoxLib::Error("Read of VisMF::FabOnDisk failed");
    }

    return is;
}

std::ostream&
operator<< (std::ostream&                  os,
            const Array<VisMF::FabOnDisk>& fa)
{
    long i(0), N(fa.size());

    os << N << '\n';

    for( ; i < N; ++i) {
        os << fa[i] << '\n';
    }

    if( ! os.good()) {
        BoxLib::Error("Write of Array<VisMF::FabOnDisk> failed");
    }

    return os;
}

std::istream&
operator>> (std::istream&            is,
            Array<VisMF::FabOnDisk>& fa)
{
    long i(0), N;

    is >> N;
    BL_ASSERT(N >= 0);

    fa.resize(N);

    for ( ; i < N; ++i) {
        is >> fa[i];
    }

    if( ! is.good()) {
        BoxLib::Error("Read of Array<VisMF::FabOnDisk> failed");
    }

    return is;
}

static
std::ostream&
operator<< (std::ostream&               os,
            const Array< Array<Real> >& ar)
{
    long i(0), N(ar.size()), M = (N == 0) ? 0 : ar[0].size();

    os << N << ',' << M << '\n';

    for( ; i < N; ++i) {
        BL_ASSERT(ar[i].size() == M);

        for(long j(0); j < M; ++j) {
            os << ar[i][j] << ',';
        }
        os << '\n';
    }

    if( ! os.good()) {
        BoxLib::Error("Write of Array<Array<Real>> failed");
    }

    return os;
}

static
std::istream&
operator>> (std::istream&         is,
            Array< Array<Real> >& ar)
{
    char ch;
    long i(0), N, M;
#ifdef BL_USE_FLOAT
    double dtemp;
#endif

    is >> N >> ch >> M;

    if( N < 0 ) {
      BoxLib::Error("Expected a positive integer, N, got something else");
    }
    if( M < 0 ) {
      BoxLib::Error("Expected a positive integer, M, got something else");
    }
    if( ch != ',' ) {
      BoxLib::Error("Expected a ',' got something else");
    }

    ar.resize(N);
    
    for( ; i < N; ++i) {
        ar[i].resize(M);

        for(long j = 0; j < M; ++j) {
#ifdef BL_USE_FLOAT
            is >> dtemp >> ch;
            ar[i][j] = static_cast<Real>(dtemp);
#else
            is >> ar[i][j] >> ch;
#endif
	    if( ch != ',' ) {
	      BoxLib::Error("Expected a ',' got something else");
	    }
        }
    }

    if( ! is.good()) {
        BoxLib::Error("Read of Array<Array<Real>> failed");
    }

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
    int oldPrec(os.precision(16));

    os << hd.m_vers     << '\n';
    os << int(hd.m_how) << '\n';
    os << hd.m_ncomp    << '\n';
    os << hd.m_ngrow    << '\n';

    hd.m_ba.writeOn(os); os << '\n';

    os << hd.m_fod      << '\n';

    if(hd.m_vers == VisMF::Header::Version_v1 ||
       hd.m_vers == VisMF::Header::NoFabHeaderMinMax_v1)
    {
      os << hd.m_min      << '\n';
      os << hd.m_max      << '\n';
    }

    if(hd.m_vers == VisMF::Header::NoFabHeaderFAMinMax_v1) {
      BL_ASSERT(hd.m_famin.size() == hd.m_ncomp);
      BL_ASSERT(hd.m_famin.size() == hd.m_famax.size());
      for(int i(0); i < hd.m_famin.size(); ++i) {
        os << hd.m_famin[i] << ',';
      }
      os << '\n';
      for(int i(0); i < hd.m_famax.size(); ++i) {
        os << hd.m_famax[i] << ',';
      }
      os << '\n';
    }

    if(hd.m_vers == VisMF::Header::NoFabHeader_v1       ||
       hd.m_vers == VisMF::Header::NoFabHeaderMinMax_v1 ||
       hd.m_vers == VisMF::Header::NoFabHeaderFAMinMax_v1)
    {
      os << FPC::NativeRealDescriptor() << '\n';
    }

    os.flags(oflags);
    os.precision(oldPrec);

    if( ! os.good()) {
        BoxLib::Error("Write of VisMF::Header failed");
    }

    return os;
}

std::istream&
operator>> (std::istream&  is,
            VisMF::Header& hd)
{
    is >> hd.m_vers;
    BL_ASSERT(hd.m_vers != VisMF::Header::Undefined_v1);

    int how;
    is >> how;
    switch(how) {
      case VisMF::OneFilePerCPU:
        hd.m_how = VisMF::OneFilePerCPU;
      break;
      case VisMF::NFiles:
        hd.m_how = VisMF::NFiles;
      break;
      default:
        BoxLib::Error("Bad case in VisMF::Header.m_how switch");
    }

    is >> hd.m_ncomp;
    BL_ASSERT(hd.m_ncomp >= 0);

    is >> hd.m_ngrow;
    BL_ASSERT(hd.m_ngrow >= 0);

    hd.m_ba.readFrom(is);

    is >> hd.m_fod;
    BL_ASSERT(hd.m_ba.size() == hd.m_fod.size());

    if(hd.m_vers == VisMF::Header::Version_v1 ||
       hd.m_vers == VisMF::Header::NoFabHeaderMinMax_v1)
    {
      is >> hd.m_min;
      is >> hd.m_max;
      BL_ASSERT(hd.m_ba.size() == hd.m_min.size());
      BL_ASSERT(hd.m_ba.size() == hd.m_max.size());
    }

    if(hd.m_vers == VisMF::Header::NoFabHeaderFAMinMax_v1) {
      char ch;
      hd.m_famin.resize(hd.m_ncomp);
      hd.m_famax.resize(hd.m_ncomp);
      for(int i(0); i < hd.m_famin.size(); ++i) {
        is >> hd.m_famin[i] >> ch;
	if( ch != ',' ) {
	  BoxLib::Error("Expected a ',' when reading hd.m_famin");
	}
      }
      for(int i(0); i < hd.m_famax.size(); ++i) {
        is >> hd.m_famax[i] >> ch;
	if( ch != ',' ) {
	  BoxLib::Error("Expected a ',' when reading hd.m_famax");
	}
      }
    }
    if(hd.m_vers == VisMF::Header::NoFabHeader_v1       ||
       hd.m_vers == VisMF::Header::NoFabHeaderMinMax_v1 ||
       hd.m_vers == VisMF::Header::NoFabHeaderFAMinMax_v1)
    {
      RealDescriptor writtenRD;
      is >> writtenRD;
      if(writtenRD != FPC::NativeRealDescriptor()) {
        BoxLib::Abort("**** Error:  NoFabHeaderFormat read different RealDescriptor.");
      }
    }


    if( ! is.good()) {
        BoxLib::Error("Read of VisMF::Header failed");
    }

    return is;
}

VisMF::FabOnDisk::FabOnDisk () {}

VisMF::FabOnDisk::FabOnDisk (const std::string& name, long offset)
    :
    m_name(name),
    m_head(offset)
{}


VisMF::FabReadLink::FabReadLink()
    :
    rankToRead(-1),
    faIndex(-1),
    fileOffset(-1)
{ }

VisMF::FabReadLink::FabReadLink(int ranktoread, int faindex, long fileoffset)
    :
    rankToRead(ranktoread),
    faIndex(faindex),
    fileOffset(fileoffset)
{ }

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

    if(m_hdr.m_min.size() == 0) {  // ---- these were not in the header
      return std::numeric_limits<int>::max();
    }

    return m_hdr.m_min[fabIndex][nComp];
}

Real
VisMF::min (int nComp) const
{
    BL_ASSERT(0 <= nComp && nComp < m_hdr.m_ncomp);

    if(m_hdr.m_famin.size() == 0) {  // ---- these were not in the header
      return std::numeric_limits<int>::max();
    }

    return m_hdr.m_famin[nComp];
}

Real
VisMF::max (int fabIndex,
            int nComp) const
{
    BL_ASSERT(0 <= fabIndex && fabIndex < m_hdr.m_ba.size());
    BL_ASSERT(0 <= nComp && nComp < m_hdr.m_ncomp);

    if(m_hdr.m_max.size() == 0) {  // ---- these were not in the header
      return -std::numeric_limits<int>::max();
    }

    return m_hdr.m_max[fabIndex][nComp];
}

Real
VisMF::max (int nComp) const
{
    BL_ASSERT(0 <= nComp && nComp < m_hdr.m_ncomp);

    if(m_hdr.m_famax.size() == 0) {  // ---- these were not in the header
      return -std::numeric_limits<int>::max();
    }

    return m_hdr.m_famax[nComp];
}

const FArrayBox&
VisMF::GetFab (int fabIndex,
               int ncomp) const
{
    if(m_pa[ncomp][fabIndex] == 0) {
        m_pa[ncomp][fabIndex] = VisMF::readFAB(fabIndex,m_fafabname,m_hdr,ncomp);
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
    return VisMF::readFAB(idx,m_fafabname,m_hdr,ncomp);
}

std::string
VisMF::BaseName (const std::string& filename)
{
    BL_ASSERT(filename[filename.length() - 1] != '/');

    if(const char *slash = strrchr(filename.c_str(), '/')) {
        //
        // Got at least one slash -- return the following tail.
        //
        return std::string(slash + 1);
    } else {
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

    const char *str = filename.c_str();    

    if(const char *slash = strrchr(str, '/')) {
        //
        // Got at least one slash -- return the dirname including last slash.
        //
        int len((slash - str) + 1);

        char *buf = new char[len+1];

        strncpy(buf, str, len);

        buf[len] = 0;   // Stringify

        std::string dirname = buf;

        delete [] buf;

        return dirname;
    } else {
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
    BL_PROFILE("VisMF::Write_fab");
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
    m_vers(VisMF::Header::Undefined_v1)
{}

//
// The more-or-less complete header only exists at IOProcessor().
//

VisMF::Header::Header (const FabArray<FArrayBox>& mf,
                       VisMF::How      how,
		       Version version)
    :
    m_vers(version),
    m_how(how),
    m_ncomp(mf.nComp()),
    m_ngrow(mf.nGrow()),
    m_ba(mf.boxArray()),
    m_fod(m_ba.size())
{
    BL_PROFILE("VisMF::Header");

    if(version == NoFabHeader_v1) {
      m_min.clear();
      m_max.clear();
      m_famin.clear();
      m_famax.clear();
      return;
    }

    if(version == NoFabHeaderFAMinMax_v1) {
      // ---- calculate FabArray min max values only
      m_min.clear();
      m_max.clear();
      m_famin.resize(m_ncomp,  std::numeric_limits<Real>::max());
      m_famax.resize(m_ncomp, -std::numeric_limits<Real>::max());

      for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const int idx = mfi.index();
        for(int i(0); i < m_ncomp; ++i) {
          m_famin[i] = std::min(m_famin[i], mf[mfi].min(m_ba[idx],i));
          m_famax[i] = std::max(m_famax[i], mf[mfi].max(m_ba[idx],i));
        }
      }
      ParallelDescriptor::ReduceRealMin(m_famin.dataPtr(), m_famin.size());
      ParallelDescriptor::ReduceRealMax(m_famax.dataPtr(), m_famax.size());

      return;
    }

    CalculateMinMax(mf);
}


void
VisMF::Header::CalculateMinMax (const FabArray<FArrayBox>& mf)
{
    m_min.resize(m_ba.size());
    m_max.resize(m_ba.size());

#ifdef BL_USE_MPI
    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();
    //
    // Calculate m_min and m_max on the CPU owning the fab.
    //
    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const int idx = mfi.index();

        m_min[idx].resize(m_ncomp);
        m_max[idx].resize(m_ncomp);

        BL_ASSERT(mf[mfi].box().contains(m_ba[idx]));

        for(long j(0); j < m_ncomp; ++j) {
            m_min[idx][j] = mf[mfi].min(m_ba[idx],j);
            m_max[idx][j] = mf[mfi].max(m_ba[idx],j);
        }
    }

    Array<int> nmtags(ParallelDescriptor::NProcs(),0);
    Array<int> offset(ParallelDescriptor::NProcs(),0);

    const Array<int>& pmap = mf.DistributionMap().ProcessorMap();

    for(int i(0), N = mf.size(); i < N; ++i) {
        ++nmtags[pmap[i]];
    }

    for(int i(0), N(nmtags.size()); i < N; ++i) {
        //
        // Each Fab corresponds to 2*m_ncomp Reals.
        //
        nmtags[i] *= 2*m_ncomp;
    }

    for(int i(1), N(offset.size()); i < N; ++i) {
        offset[i] = offset[i-1] + nmtags[i-1];
    }

    Array<Real> senddata(nmtags[ParallelDescriptor::MyProc()]);

    if(senddata.empty()) {
        //
        // Can't let senddata be empty as senddata.dataPtr() will fail.
        //
        senddata.resize(1);
    }

    int ioffset = 0;

    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const int idx = mfi.index();
        for(int i(0); i < m_ncomp; ++i) {
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
                                IOProcNumber,
                                ParallelDescriptor::Communicator()) );

    BL_COMM_PROFILE(BLProfiler::Gatherv, recvdata.size() * sizeof(Real),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    if(ParallelDescriptor::IOProcessor()) {
        for(int i(0), N(mf.size()); i < N; ++i) {
            if(pmap[i] != IOProcNumber) {
                m_min[i].resize(m_ncomp);
                m_max[i].resize(m_ncomp);
            }
        }

        for(int j(0), N(mf.size()); j < N; ++j) {
            if(pmap[j] != IOProcNumber) {
                for(int k(0); k < m_ncomp; ++k) {
                    m_min[j][k] = recvdata[offset[pmap[j]]+k];
                    m_max[j][k] = recvdata[offset[pmap[j]]+k+m_ncomp];
                }

                offset[pmap[j]] += 2*m_ncomp;
            }
        }
    }
#else
    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const int idx = mfi.index();

        m_min[idx].resize(m_ncomp);
        m_max[idx].resize(m_ncomp);

        BL_ASSERT(mf[mfi].box().contains(m_ba[idx]));

        for(long j(0); j < m_ncomp; ++j) {
            m_min[idx][j] = mf[mfi].min(m_ba[idx],j);
            m_max[idx][j] = mf[mfi].max(m_ba[idx],j);
        }
    }
#endif /*BL_USE_MPI*/

#ifdef BL_FIXHEADERDENORMS
    if(ParallelDescriptor::IOProcessor()) {
        for(int i(0); i < m_min.size(); ++i) {
            for(int j(0); j < m_min[i].size(); ++j) {
                if(std::abs(m_min[i][j]) < 1.0e-300) {
                    m_min[i][j] = 0.0;
                }
            }
        }
        for(int i(0); i < m_max.size(); ++i) {
            for(int j(0); j < m_max[i].size(); ++j) {
                if(std::abs(m_max[i][j]) < 1.0e-300) {
                    m_max[i][j] = 0.0;
                }
            }
        }
    }
#endif

    // ---- calculate fabarray min max values
    m_famin.resize(m_ncomp);
    m_famax.resize(m_ncomp);
    for(int i(0); i < m_ncomp; ++i) {
      m_famin[i] =  std::numeric_limits<Real>::max();
      m_famax[i] = -std::numeric_limits<Real>::max();
      for(int j(0); j < m_min[i].size(); ++j) {
        m_famin[i] = std::min(m_famin[i], m_min[i][j]);
      }
      for(int j(0); j < m_max[i].size(); ++j) {
        m_famax[i] = std::max(m_famax[i], m_max[i][j]);
      }
    }
}

long
VisMF::WriteHeader (const std::string& mf_name,
                    VisMF::Header&     hdr)
{
    BL_PROFILE("VisMF::WriteHeader");
    long bytesWritten(0);

    if(ParallelDescriptor::IOProcessor()) {
        std::string MFHdrFileName(mf_name);

        MFHdrFileName += TheMultiFabHdrFileSuffix;

        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

        std::ofstream MFHdrFile;

        MFHdrFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        MFHdrFile.open(MFHdrFileName.c_str(), std::ios::out|std::ios::trunc);

        if( ! MFHdrFile.good()) {
            BoxLib::FileOpenFailed(MFHdrFileName);
	}

        MFHdrFile << hdr;
        //
        // Add in the number of bytes written out in the Header.
        //
        bytesWritten += VisMF::FileOffset(MFHdrFile);

        MFHdrFile.close();
    }
    return bytesWritten;
}

long
VisMF::Write (const FabArray<FArrayBox>&    mf,
              const std::string& mf_name,
              VisMF::How         how,
              bool               set_ghost)
{
    BL_PROFILE("VisMF::Write_FabArray");
    BL_ASSERT(mf_name[mf_name.length() - 1] != '/');

    // ---- add stream retry
    // ---- add stream buffer (to nfiles)

    VisMF::Initialize();

    if(set_ghost) {
        FabArray<FArrayBox>* the_mf = const_cast<FabArray<FArrayBox>*>(&mf);

        for(MFIter mfi(*the_mf); mfi.isValid(); ++mfi) {
            const int idx(mfi.index());

            for(int j(0); j < mf.nComp(); ++j) {
                const Real valMin(mf[mfi].min(mf.box(idx), j));
                const Real valMax(mf[mfi].max(mf.box(idx), j));
                const Real val((valMin + valMax) / 2.0);

                the_mf->get(mfi).setComplement(val, mf.box(idx), j, 1);
            }
        }
    }

    VisMF::Header hdr(mf, how, currentVersion);

    long bytesWritten(0);
    const int myProc(ParallelDescriptor::MyProc());

    std::string filePrefix(mf_name + FabFileSuffix);

    if(currentVersion == VisMF::Header::Version_v1) {  // ---- write the old way

      for(NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf); nfi.ReadyToWrite(); ++nfi) {
        const std::string fName(NFilesIter::FileName(nOutFiles, filePrefix, myProc, groupSets));
        for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
          hdr.m_fod[mfi.index()] = VisMF::Write(mf[mfi], fName,nfi.Stream(), bytesWritten);
        }
      }

    } else {     // ---- write the fab data directly

      for(NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf); nfi.ReadyToWrite(); ++nfi) {
	if(VisMF::useSingleWrite) {
          for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
            bytesWritten += mf[mfi].nBytes();
          }

	  char *allFabData;
	  bool goodAlloc(true);
	  allFabData = new(std::nothrow) char[bytesWritten];
	  if(allFabData == nullptr) {
	    goodAlloc = false;
	  }

	  if(goodAlloc) {
	    long writePosition(0);
            for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
              const FArrayBox &fab = mf[mfi];
	      char *afPtr = allFabData + writePosition;
	      memcpy(afPtr, fab.dataPtr(), fab.nBytes());
              writePosition += fab.nBytes();
            }
            nfi.Stream().write(allFabData, bytesWritten);
	    delete [] allFabData;

	  } else {    // ---- cannot use one write

            for(MFIter mfi(mf); mfi.isValid(); ++mfi) {    // ---- write the fab data directly
              const FArrayBox &fab = mf[mfi];
              nfi.Stream().write((char *) fab.dataPtr(), fab.nBytes());
              bytesWritten += fab.nBytes();
            }
	  }

	} else {

          for(MFIter mfi(mf); mfi.isValid(); ++mfi) {    // ---- write the fab data directly
            const FArrayBox &fab = mf[mfi];
            nfi.Stream().write((char *) fab.dataPtr(), fab.nBytes());
            bytesWritten += fab.nBytes();
          }
	}
      }

    }

    VisMF::FindOffsets(mf, filePrefix, hdr, groupSets, currentVersion);

    bytesWritten += VisMF::WriteHeader(mf_name, hdr);

    return bytesWritten;
}


void
VisMF::FindOffsets (const FabArray<FArrayBox> &mf,
		    const std::string &filePrefix,
                    VisMF::Header &hdr,
		    bool groupSets,
		    VisMF::Header::Version whichVersion)
{
    BL_PROFILE("VisMF::FindOffsets");

    if(FArrayBox::getFormat() == FABio::FAB_ASCII ||
       FArrayBox::getFormat() == FABio::FAB_8BIT)
    {
#ifdef BL_USE_MPI
    const int IOProcNumber(ParallelDescriptor::IOProcessorNumber());
    const int nProcs(ParallelDescriptor::NProcs());

    Array<int> nmtags(nProcs,0);
    Array<int> offset(nProcs,0);

    const Array<int> &pmap = mf.DistributionMap().ProcessorMap();

    for(int i(0), N(mf.size()); i < N; ++i) {
        ++nmtags[pmap[i]];
    }

    for(int i(1), N(offset.size()); i < N; ++i) {
        offset[i] = offset[i-1] + nmtags[i-1];
    }

    Array<long> senddata(nmtags[ParallelDescriptor::MyProc()]);

    if(senddata.empty()) {
      // Can't let senddata be empty as senddata.dataPtr() will fail.
      senddata.resize(1);
    }

    int ioffset(0);

    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
      senddata[ioffset++] = hdr.m_fod[mfi.index()].m_head;
    }

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
                                IOProcNumber,
                                ParallelDescriptor::Communicator()) );

    BL_COMM_PROFILE(BLProfiler::Gatherv, recvdata.size() * sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    if(ParallelDescriptor::IOProcessor()) {
        Array<int> cnt(nProcs,0);

        for(int j(0), N(mf.size()); j < N; ++j) {
            const int i(pmap[j]);
            hdr.m_fod[j].m_head = recvdata[offset[i]+cnt[i]];

            const std::string name(NFilesIter::FileName(nOutFiles, filePrefix, i, groupSets));

            hdr.m_fod[j].m_name = VisMF::BaseName(name);

            ++cnt[i];
        }
    }
#endif /*BL_USE_MPI*/


    } else {  // ---- calculate offsets
      RealDescriptor *whichRD;
      if(FArrayBox::getFormat() == FABio::FAB_NATIVE) {
        whichRD = FPC::NativeRealDescriptor().clone();
      } else if(FArrayBox::getFormat() == FABio::FAB_NATIVE_32) {
        whichRD = FPC::Native32RealDescriptor().clone();
      } else if(FArrayBox::getFormat() == FABio::FAB_IEEE_32) {
        whichRD = FPC::NativeRealDescriptor().clone();
      }
      const  FABio &fio = FArrayBox::getFABio();
      int whichRDBytes(whichRD->numBytes());
      int nComps(mf.nComp());

      if(ParallelDescriptor::IOProcessor()) {   // ---- calculate offsets
	const BoxArray &mfBA = mf.boxArray();
	const DistributionMapping &mfDM = mf.DistributionMap();
	Array<long> fabHeaderBytes(mfBA.size(), 0);
	int nFiles(NFilesIter::ActualNFiles(nOutFiles));
	int whichFileNumber(-1), whichProc(-1);
	std::string whichFileName;
	Array<std::streampos> currentOffset(nFiles, 0);

        if(hdr.m_vers == VisMF::Header::Version_v1) {
	  for(int i(0); i < mfBA.size(); ++i) {
            std::stringstream hss;
	    FArrayBox tempFab(mf.fabbox(i), nComps, false);  // ---- no alloc
            fio.write_header(hss, tempFab, tempFab.nComp());
	    fabHeaderBytes[i] = hss.tellp();
	  }
	}

	std::map<int, Array<int> > rankBoxOrder;  // ---- [rank, index]

	for(int i(0); i < mfBA.size(); ++i) {
	  rankBoxOrder[mfDM[i]].push_back(i);
	}

	std::map<int, Array<int> >::iterator rboIter;
	for(rboIter = rankBoxOrder.begin(); rboIter != rankBoxOrder.end(); ++rboIter) {
	  Array<int> &index = rboIter->second;
	  whichProc = rboIter->first;
	  whichFileNumber = NFilesIter::FileNumber(nFiles, whichProc, groupSets);
	  whichFileName   = NFilesIter::FileName(nFiles, filePrefix, whichProc, groupSets);
	  for(int i(0); i < index.size(); ++i) {
	    hdr.m_fod[index[i]].m_name = VisMF::BaseName(whichFileName);
	    hdr.m_fod[index[i]].m_head = currentOffset[whichFileNumber];
	    //currentOffset[whichFileNumber] += mfBA[index[i]].numPts() * nComps * whichRDBytes
	    currentOffset[whichFileNumber] += mf.fabbox(index[i]).numPts() * nComps * whichRDBytes
	                                      + fabHeaderBytes[index[i]];
	  }
	}
      }
    
      delete whichRD;
    }
}


void
VisMF::RemoveFiles(const std::string &mf_name, bool verbose)
{
    if(ParallelDescriptor::IOProcessor()) {
      std::string MFHdrFileName(mf_name + TheMultiFabHdrFileSuffix);
      if(verbose) {
        std::cout << "---- removing:  " << MFHdrFileName << std::endl;
      }
      int retVal(std::remove(MFHdrFileName.c_str()));
      if(verbose) {
        if(retVal != 0) {
          std::cout << "---- error removing:  " << MFHdrFileName << "  errno = "
	            << strerror(errno) << std::endl;
        }
      }
      for(int ip(0); ip < nOutFiles; ++ip) {
        std::string fileName(NFilesIter::FileName(nOutFiles, mf_name + FabFileSuffix, ip, true));
        if(verbose) {
          std::cout << "---- removing:  " << fileName << std::endl;
	}
        int retVal(std::remove(fileName.c_str()));
        if(verbose) {
          if(retVal != 0) {
            std::cout << "---- error removing:  " << fileName << "  errno = "
	              << strerror(errno) << std::endl;
          }
	}
      }
    }
}


VisMF::VisMF (const std::string &fafab_name)
    :
    m_fafabname(fafab_name)
{
    std::string FullHdrFileName(m_fafabname);

    FullHdrFileName += TheMultiFabHdrFileSuffix;

    Array<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(FullHdrFileName, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream ifs(fileCharPtrString, std::istringstream::in);

    ifs >> m_hdr;

    m_pa.resize(m_hdr.m_ncomp);

    for(int nComp(0); nComp < m_pa.size(); ++nComp) {
        m_pa[nComp].resize(m_hdr.m_ba.size());

        for(int ii(0), N(m_pa[nComp].size()); ii < N; ++ii) {
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
    BL_PROFILE("VisMF::readFAB_idx");
    Box fab_box = hdr.m_ba[idx];

    if(hdr.m_ngrow) {
        fab_box.grow(hdr.m_ngrow);
    }

    FArrayBox *fab = new FArrayBox(fab_box, ncomp == -1 ? hdr.m_ncomp : 1);

    std::string FullName(VisMF::DirName(mf_name));

    FullName += hdr.m_fod[idx].m_name;
    
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ifstream ifs;

    ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    ifs.open(FullName.c_str(), std::ios::in|std::ios::binary);

    if( ! ifs.good()) {
        BoxLib::FileOpenFailed(FullName);
    }

    if(hdr.m_fod[idx].m_head) {
        ifs.seekg(hdr.m_fod[idx].m_head, std::ios::beg);
    }

    if(hdr.m_vers == Header::Version_v1) {
      if(ncomp == -1) {
        fab->readFrom(ifs);
      } else {
        fab->readFrom(ifs, ncomp);
      }
    } else {
      if(ncomp == -1) {
        ifs.read((char *) fab->dataPtr(), fab->nBytes());
      } else {
	const RealDescriptor &nrd = FPC::NativeRealDescriptor();
        long bytesPerComp(fab->box().numPts() * nrd.numBytes());
        ifs.seekg(bytesPerComp * ncomp, std::ios::cur);
        ifs.read((char *) fab->dataPtr(), bytesPerComp);
      }
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
    BL_PROFILE("VisMF::readFAB_mf");
    FArrayBox& fab = mf[idx];

    std::string FullName(VisMF::DirName(mf_name));

    FullName += hdr.m_fod[idx].m_name;
    
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ifstream ifs;

    ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    ifs.open(FullName.c_str(), std::ios::in|std::ios::binary);

    if( ! ifs.good()) {
        BoxLib::FileOpenFailed(FullName);
    }

    if(hdr.m_fod[idx].m_head) {
        ifs.seekg(hdr.m_fod[idx].m_head, std::ios::beg);
    }
    fab.readFrom(ifs);

    ifs.close();
}


void
VisMF::readFABNoFabHeader (FabArray<FArrayBox>&            mf,
		           int                  idx,
                           const std::string&   mf_name,
                           const VisMF::Header& hdr)
{
    BL_PROFILE("VisMF::readFABNoFabHeader_mf");
    FArrayBox &fab = mf[idx];

    std::string FullName(VisMF::DirName(mf_name));

    FullName += hdr.m_fod[idx].m_name;
    
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ifstream ifs;

    ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    ifs.open(FullName.c_str(), std::ios::in|std::ios::binary);

    if( ! ifs.good()) {
        BoxLib::FileOpenFailed(FullName);
    }

    if(hdr.m_fod[idx].m_head) {
        ifs.seekg(hdr.m_fod[idx].m_head, std::ios::beg);
    }

    ifs.read((char *) fab.dataPtr(), fab.nBytes());

    //// test ifs.gcount and tellg vs count

    ifs.close();
}


void
VisMF::Read (FabArray<FArrayBox> &mf,
             const std::string   &mf_name,
	     const char *faHeader)
{
    BL_PROFILE("VisMF::Read()");

    VisMF::Header hdr;
    Real hEndTime, hStartTime, faCopyTime(0.0);
    Real startTime(ParallelDescriptor::second());
    int nOpensPerFile(nMFFileInStreams), messTotal(0);
    static Real totalTime(0.0);

    if(verbose && ParallelDescriptor::IOProcessor()) {
      std::cout << "VisMF::Read:  about to read:  " << mf_name << std::endl;
    }

    VisMF::Initialize();

    std::string FullHdrFileName(mf_name + TheMultiFabHdrFileSuffix);
    {
        hStartTime = ParallelDescriptor::second();
        std::string fileCharPtrString;
	if(faHeader == 0) {
          Array<char> fileCharPtr;
          ParallelDescriptor::ReadAndBcastFile(FullHdrFileName, fileCharPtr);
          fileCharPtrString = fileCharPtr.dataPtr();
	} else {
          fileCharPtrString = faHeader;
	}
        std::istringstream ifs(fileCharPtrString, std::istringstream::in);

        ifs >> hdr;

        hEndTime = ParallelDescriptor::second();
    }

    bool noFabHeader(false);
    if(hdr.m_vers == VisMF::Header::NoFabHeader_v1       ||
       hdr.m_vers == VisMF::Header::NoFabHeaderMinMax_v1 ||
       hdr.m_vers == VisMF::Header::NoFabHeaderFAMinMax_v1)
    {
      noFabHeader = true;
    }


    mf.define(hdr.m_ba, hdr.m_ncomp, hdr.m_ngrow, Fab_allocate);

#ifdef BL_USE_MPI

if(noFabHeader) {

    // ---- Create an ordered map of which processors read which
    // ---- Fabs in each file
    
    std::map<std::string, Array<FabReadLink> > FileReadChains;        // ---- [filename, chain]
    std::map<std::string, Array<FabReadLink> > FileReadChainsSorted;  // ---- [filename, chain]

    int nBoxes(hdr.m_ba.size());
    for(int i(0); i < nBoxes; ++i) {   // ---- create the map
      int whichProc(mf.DistributionMap()[i]);
      std::string fname(hdr.m_fod[i].m_name);
      FileReadChains[fname].push_back(FabReadLink(whichProc, i, hdr.m_fod[i].m_head));
    }

    // ---- This code is only for reading in file order
    bool inFileOrder(true);
    std::map<std::string, Array<FabReadLink> >::iterator frcIter;
    std::string readFileName;
    std::map<std::string, std::set<int> > readFileRanks;  // ---- [filename, ranks]

    int indexFileOrder(0);
    FabArray<FArrayBox> fafabFileOrder;
    BoxArray baFileOrder(hdr.m_ba.size());
    Array<int> ranksFileOrder(mf.DistributionMap().size());
    ranksFileOrder[ranksFileOrder.size() - 1] = ParallelDescriptor::MyProc();

    for(frcIter = FileReadChains.begin(); frcIter != FileReadChains.end(); ++frcIter) {
      // ---- sort to rank ordering and compare with original
      const std::string &fileName = frcIter->first;
      Array<FabReadLink> &frc = frcIter->second;
      Array<FabReadLink> frcSorted = frcIter->second;  // ---- make a copy
      std::stable_sort(frcSorted.begin(), frcSorted.end(),
	               [] (const FabReadLink &a, const FabReadLink &b)
	                    { return a.rankToRead < b.rankToRead; } );

      for(int i(0); i < frc.size(); ++i) {
        if(frc[i].rankToRead != frcSorted[i].rankToRead ||
           frc[i].fileOffset != frcSorted[i].fileOffset)
	{
	  inFileOrder = false;
	}
        readFileRanks[fileName].insert(frcSorted[i].rankToRead);

	ranksFileOrder[indexFileOrder] = frcSorted[i].rankToRead;
	baFileOrder.set(indexFileOrder, hdr.m_ba[frc[i].faIndex]);

        FileReadChainsSorted[fileName].push_back(FabReadLink(frcSorted[i].rankToRead,
	                                                     indexFileOrder,
							     frcSorted[i].fileOffset));
	++indexFileOrder;
      }
    }

    DistributionMapping dmFileOrder(ranksFileOrder);

    if(inFileOrder) {
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "OOOOOOOO:  inFileOrder" << std::endl;
      }
    } else {
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "OOOOOOOO:  not inFileOrder" << std::endl;
      }
      // ---- make a temporary fabarray in file order
      fafabFileOrder.define(baFileOrder, hdr.m_ncomp, hdr.m_ngrow, dmFileOrder, Fab_allocate);
    }

    FabArray<FArrayBox> &whichFA = inFileOrder ? mf : fafabFileOrder;

    // ---- check that a rank only needs to read one file
    std::map<std::string, std::set<int> >::iterator rfrIter;
    std::set<int>::iterator setIter;

    for(rfrIter = readFileRanks.begin(); rfrIter != readFileRanks.end(); ++rfrIter) {
      std::set<int> &rfrSplitSet = rfrIter->second;
      if(rfrSplitSet.size() == 0) {
        continue;
      }
      // ---- split the set into nstreams sets
      int ssSize(rfrSplitSet.size());
      int nStreams(std::min(ssSize, nOpensPerFile));
      int ranksPerStream(ssSize / nStreams); // ---- plus some remainder... 
      Array<std::set<int> > streamSets(nStreams);
      int sIndex(0), sCount(0);
      for(setIter = rfrSplitSet.begin(); setIter != rfrSplitSet.end(); ++setIter) {
        streamSets[sIndex].insert(*setIter);
	if(++sCount >= ranksPerStream) {
	  sCount = 0;
	  ++sIndex;
	  sIndex = std::min(sIndex, streamSets.size() - 1);
	}
      }

      for(int iSet(0); iSet < streamSets.size(); ++iSet) {
        Array<int> readRanks;
        std::set<int> &rfrSet = streamSets[iSet];
        for(setIter = rfrSet.begin(); setIter != rfrSet.end(); ++setIter) {
          readRanks.push_back(*setIter);
        }

        const int myProc(ParallelDescriptor::MyProc());
        if(rfrSet.find(myProc) != rfrSet.end()) {  // ---- myProc needs to read this file
          const std::string &fileName = rfrIter->first;
	  std::string fullFileName(VisMF::DirName(mf_name) + fileName);
	  if(inFileOrder) {
	    frcIter = FileReadChains.find(fileName);
	  } else {
	    frcIter = FileReadChainsSorted.find(fileName);
	  }
          Array<FabReadLink> &frc = frcIter->second;
          for(NFilesIter nfi(fullFileName, readRanks); nfi.ReadyToRead(); ++nfi) {
	    if(VisMF::useSingleRead) {
	      // ---- confirm the data is contiguous in the stream
	      long firstOffset(-1);
	      for(int i(0); i < frc.size(); ++i) {
	        if(myProc == frc[i].rankToRead) {
		  firstOffset = frc[i].fileOffset;
		  break;
		}
	      }

	      bool dataIsContiguous(true);
	      long currentOffset(firstOffset);
	      long bytesToRead(0);

	      for(int i(0); i < frc.size(); ++i) {
	        if(myProc == frc[i].rankToRead) {
	          if(currentOffset != frc[i].fileOffset) {
                    dataIsContiguous = false;
	          } else {
	            FArrayBox &fab = whichFA[frc[i].faIndex];
                    currentOffset += fab.nBytes();
                    bytesToRead   += fab.nBytes();
		  }
	        }
	      }
	      char *allFabData;
	      bool goodAlloc(true);
	      if(dataIsContiguous) {
	        allFabData = new(std::nothrow) char[bytesToRead];
		if(allFabData == nullptr) {
		  goodAlloc = false;
		}
	      }
	      if(goodAlloc) {
                nfi.Stream().seekp(firstOffset, std::ios::beg);
                nfi.Stream().read(allFabData, bytesToRead);

		currentOffset = 0;  // ---- this is now relative to allFabData

	        for(int i(0); i < frc.size(); ++i) {
	          if(myProc == frc[i].rankToRead) {
		    char *afPtr = allFabData + currentOffset;
	            FArrayBox &fab = whichFA[frc[i].faIndex];
                    memcpy(fab.dataPtr(), afPtr, fab.nBytes());
                    currentOffset += fab.nBytes();
	          }
	        }
		delete [] allFabData;

	      } else {          // ---- cannot use one read
	        for(int i(0); i < frc.size(); ++i) {
	          if(myProc == frc[i].rankToRead) {
	            if(nfi.SeekPos() != frc[i].fileOffset) {
                      nfi.Stream().seekp(frc[i].fileOffset, std::ios::beg);
	            }
	            FArrayBox &fab = whichFA[frc[i].faIndex];
                    nfi.Stream().read((char *) fab.dataPtr(), fab.nBytes());
	          }
	        }
	      }

	    } else {
	      for(int i(0); i < frc.size(); ++i) {
	        if(myProc == frc[i].rankToRead) {
	          if(nfi.SeekPos() != frc[i].fileOffset) {
                    nfi.Stream().seekp(frc[i].fileOffset, std::ios::beg);
	          }
	          FArrayBox &fab = whichFA[frc[i].faIndex];
                  nfi.Stream().read((char *) fab.dataPtr(), fab.nBytes());
	        }
	      }
	    }
          }
        }

      }
    }

    if( ! inFileOrder) {
      faCopyTime = ParallelDescriptor::second();
      mf.copy(fafabFileOrder);
      faCopyTime = ParallelDescriptor::second() - faCopyTime;
    }

} else {

    //
    // Here we limit the number of open files when reading a multifab.
    //
    int nReqs(0), ioProcNum(ParallelDescriptor::IOProcessorNumber());
    int myProc(ParallelDescriptor::MyProc());
    int nBoxes(hdr.m_ba.size());
    int totalIOReqs(nBoxes), nFiles(-1);
    std::vector<int> iDone(2);
    const int iDoneIndex(0), iDoneCount(1);
    std::set<int> busyProcs;  // [whichProc]
    std::map<std::string, int> fileNames;  // <filename, allreadsindex>
    std::multiset<int> availableFiles;  // [whichFile]  supports multiple reads/file
    int allReadsIndex(0);
    ParallelDescriptor::Message rmess;
    Array<std::map<int,std::map<int,int> > > allReads; // [file]<proc,<seek,index>>
    Array<std::ifstream *> dataStreams;        // ---- persistent streams
    std::map<std::string, int> dataFileNames;  // ---- [filename, stream index]


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
	  int index(iopReads.front());
	  if(noFabHeader) {
	    VisMF::readFABNoFabHeader(mf,index, mf_name, hdr);
	  } else {
	    VisMF::readFAB(mf,index, mf_name, hdr);
	  }
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
	  if(noFabHeader) {
	    VisMF::readFABNoFabHeader(mf,mfIndex, mf_name, hdr);
	  } else {
	    VisMF::readFAB(mf,mfIndex, mf_name, hdr);
	  }
	}
        nReqs -= rmess.count();
	iDone[iDoneIndex] = recReads[0];
	iDone[iDoneCount] = rmess.count();
        ParallelDescriptor::Send(iDone, ioProcNum, recReads[0]);
      }
    }


    //ParallelDescriptor::Barrier("VisMF::Read");
}

    bool reportMFReadStats(true);
    if(ParallelDescriptor::IOProcessor() && reportMFReadStats) {
      Real mfReadTime = ParallelDescriptor::second() - startTime;
      totalTime += mfReadTime;
      std::cout << "FARead ::  nBoxes = " << hdr.m_ba.size()
                << "  nMessages = " << messTotal << '\n';
      std::cout << "FARead ::  hTime = " << (hEndTime - hStartTime) << '\n';
      std::cout << "FARead ::  faCopyTime = " << faCopyTime << '\n';
      std::cout << "FARead ::  mfReadTime = " << mfReadTime
                << "  totalTime = " << totalTime << std::endl;
    }

#else
    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
	if(noFabHeader) {
	  VisMF::readFABNoFabHeader(mf,mfi.index(), mf_name, hdr);
	} else {
	  VisMF::readFAB(mf,mfi.index(), mf_name, hdr);
	}
    }

#endif

    BL_ASSERT(mf.ok());
}


void
VisMF::ReadFAHeader (const std::string &fafabName,
	             Array<char> &faHeader)
{
    BL_PROFILE("VisMF::ReadFAHeader()");

    VisMF::Initialize();

    std::string FullHdrFileName(fafabName + TheMultiFabHdrFileSuffix);
    ParallelDescriptor::ReadAndBcastFile(FullHdrFileName, faHeader);
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
    for(int ncomp(0), N(m_pa.size()); ncomp < N; ++ncomp) {
        clear(ncomp, fabIndex);
    }
}

void
VisMF::clear ()
{
    for(int ncomp(0), N(m_pa.size()); ncomp < N; ++ncomp) {
        for(int fabIndex(0), M(m_pa[ncomp].size()); fabIndex < M; ++fabIndex) {
            clear(ncomp, fabIndex);
        }
    }
}
