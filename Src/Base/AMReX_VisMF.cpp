// TODO: need to work on read for upc++

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <deque>
#include <cerrno>

#include <AMReX_ccse-mpi.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_NFiles.H>
#include <AMReX_FPC.H>

namespace amrex {

static const char *TheMultiFabHdrFileSuffix = "_H";
static const char *FabFileSuffix = "_D_";
static const char *TheFabOnDiskPrefix = "FabOnDisk:";

std::map<std::string, VisMF::PersistentIFStream> VisMF::persistentIFStreams;

int VisMF::verbose(0);
VisMF::Header::Version VisMF::currentVersion(VisMF::Header::Version_v1);
bool VisMF::groupSets(false);
bool VisMF::setBuf(true);
bool VisMF::useSingleRead(false);
bool VisMF::useSingleWrite(false);
bool VisMF::checkFilePositions(false);
bool VisMF::usePersistentIFStreams(false);
bool VisMF::useSynchronousReads(false);
bool VisMF::useDynamicSetSelection(true);
bool VisMF::allowSparseWrites(true);

long VisMF::ioBufferSize(VisMF::IO_Buffer_Size);


//
// Set these in Initialize().
//
int VisMF::nOutFiles(256);
int VisMF::nMFFileInStreams(4);

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

    amrex::ExecOnFinalize(VisMF::Finalize);

    ParmParse pp("vismf");
    pp.query("v",verbose);

    int headerVersion(currentVersion);
    pp.query("headerversion", headerVersion);
    if(headerVersion != currentVersion) {
      currentVersion = static_cast<VisMF::Header::Version> (headerVersion);
    }

    pp.query("groupsets", groupSets);
    pp.query("setbuf", setBuf);
    pp.query("usesingleread", useSingleRead);
    pp.query("usesinglewrite", useSingleWrite);
    pp.query("checkfilepositions", checkFilePositions);
    pp.query("usepersistentifstreams", usePersistentIFStreams);
    pp.query("usesynchronousreads", useSynchronousReads);
    pp.query("usedynamicsetselection", useDynamicSetSelection);
    pp.query("iobuffersize", ioBufferSize);
    pp.query("allowsparsewrites", allowSparseWrites);

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
        amrex::Error("Write of VisMF::FabOnDisk failed");
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
        amrex::Error("Read of VisMF::FabOnDisk failed");
    }

    return is;
}

std::ostream&
operator<< (std::ostream&                  os,
            const Vector<VisMF::FabOnDisk>& fa)
{
    long i(0), N(fa.size());

    os << N << '\n';

    for( ; i < N; ++i) {
        os << fa[i] << '\n';
    }

    if( ! os.good()) {
        amrex::Error("Write of Vector<VisMF::FabOnDisk> failed");
    }

    return os;
}

std::istream&
operator>> (std::istream&            is,
            Vector<VisMF::FabOnDisk>& fa)
{
    long i(0), N;

    is >> N;
    BL_ASSERT(N >= 0);

    fa.resize(N);

    for ( ; i < N; ++i) {
        is >> fa[i];
    }

    if( ! is.good()) {
        amrex::Error("Read of Vector<VisMF::FabOnDisk> failed");
    }

    return is;
}

static
std::ostream&
operator<< (std::ostream&               os,
            const Vector< Vector<Real> >& ar)
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
        amrex::Error("Write of Vector<Vector<Real>> failed");
    }

    return os;
}

static
std::istream&
operator>> (std::istream&         is,
            Vector< Vector<Real> >& ar)
{
    char ch;
    long i(0), N, M;
#ifdef BL_USE_FLOAT
    double dtemp;
#endif

    is >> N >> ch >> M;

    if( N < 0 ) {
      amrex::Error("Expected a positive integer, N, got something else");
    }
    if( M < 0 ) {
      amrex::Error("Expected a positive integer, M, got something else");
    }
    if( ch != ',' ) {
      amrex::Error("Expected a ',' got something else");
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
	      amrex::Error("Expected a ',' got something else");
	    }
        }
    }

    if( ! is.good()) {
        amrex::Error("Read of Vector<Vector<Real>> failed");
    }

    return is;
}

std::ostream&
operator<< (std::ostream        &os,
            const VisMF::Header &hd)
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
      if(FArrayBox::getFormat() == FABio::FAB_NATIVE) {
        os << FPC::NativeRealDescriptor() << '\n';
      } else if(FArrayBox::getFormat() == FABio::FAB_NATIVE_32) {
        os << FPC::Native32RealDescriptor() << '\n';
      } else if(FArrayBox::getFormat() == FABio::FAB_IEEE_32) {
        os << FPC::Ieee32NormalRealDescriptor() << '\n';
      }
    }

    os.flags(oflags);
    os.precision(oldPrec);

    if( ! os.good()) {
        amrex::Error("Write of VisMF::Header failed");
    }

    return os;
}

std::istream&
operator>> (std::istream  &is,
            VisMF::Header &hd)
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
        amrex::Error("Bad case in VisMF::Header.m_how switch");
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
	  amrex::Error("Expected a ',' when reading hd.m_famin");
	}
      }
      for(int i(0); i < hd.m_famax.size(); ++i) {
        is >> hd.m_famax[i] >> ch;
	if( ch != ',' ) {
	  amrex::Error("Expected a ',' when reading hd.m_famax");
	}
      }
    }
    if(hd.m_vers == VisMF::Header::NoFabHeader_v1       ||
       hd.m_vers == VisMF::Header::NoFabHeaderMinMax_v1 ||
       hd.m_vers == VisMF::Header::NoFabHeaderFAMinMax_v1)
    {
      is >> hd.m_writtenRD;
    }


    if( ! is.good()) {
        amrex::Error("Read of VisMF::Header failed");
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

VisMF::FabReadLink::FabReadLink(int ranktoread, int faindex, long fileoffset,
                                const Box &b)
    :
    rankToRead(ranktoread),
    faIndex(faindex),
    fileOffset(fileoffset),
    box(b)
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
            int nc) const
{
    BL_ASSERT(0 <= fabIndex && fabIndex < m_hdr.m_ba.size());
    BL_ASSERT(0 <= nc && nc < m_hdr.m_ncomp);

    if(m_hdr.m_min.size() == 0) {  // ---- these were not in the header
      return std::numeric_limits<int>::max();
    }

    return m_hdr.m_min[fabIndex][nc];
}

Real
VisMF::min (int nc) const
{
    BL_ASSERT(0 <= nc && nc < m_hdr.m_ncomp);

    if(m_hdr.m_famin.size() == 0) {  // ---- these were not in the header
      return std::numeric_limits<int>::max();
    }

    return m_hdr.m_famin[nc];
}

Real
VisMF::max (int fabIndex,
            int nc) const
{
    BL_ASSERT(0 <= fabIndex && fabIndex < m_hdr.m_ba.size());
    BL_ASSERT(0 <= nc && nc < m_hdr.m_ncomp);

    if(m_hdr.m_max.size() == 0) {  // ---- these were not in the header
      return -std::numeric_limits<int>::max();
    }

    return m_hdr.m_max[fabIndex][nc];
}

Real
VisMF::max (int nc) const
{
    BL_ASSERT(0 <= nc && nc < m_hdr.m_ncomp);

    if(m_hdr.m_famax.size() == 0) {  // ---- these were not in the header
      return -std::numeric_limits<int>::max();
    }

    return m_hdr.m_famax[nc];
}

const FArrayBox&
VisMF::GetFab (int fabIndex,
               int ncomp) const
{
    if(m_pa[ncomp][fabIndex] == 0) {
        m_pa[ncomp][fabIndex] = VisMF::readFAB(fabIndex, m_fafabname, m_hdr, ncomp);
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
    return VisMF::readFAB(idx, mf_name, m_hdr, -1);
}

FArrayBox*
VisMF::readFAB (int idx,
		int ncomp)
{
    return VisMF::readFAB(idx, m_fafabname, m_hdr, ncomp);
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
		       Version version,
		       bool calcMinMax)
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

    if(calcMinMax) {
      CalculateMinMax(mf);
    }
}


void
VisMF::Header::CalculateMinMax (const FabArray<FArrayBox>& mf,
                                int procToWrite)
{
    BL_PROFILE("VisMF::CalculateMinMax");

    m_min.resize(m_ba.size());
    m_max.resize(m_ba.size());

#ifdef BL_USE_MPI
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

    Vector<int> nmtags(ParallelDescriptor::NProcs(), 0);
    Vector<int> offset(ParallelDescriptor::NProcs(), 0);

    const Vector<int> &pmap = mf.DistributionMap().ProcessorMap();

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

    Vector<Real> senddata(nmtags[ParallelDescriptor::MyProc()]);

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

    Vector<Real> recvdata(mf.size()*2*m_ncomp);

    BL_COMM_PROFILE(BLProfiler::Gatherv, recvdata.size() * sizeof(Real),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());

    BL_MPI_REQUIRE( MPI_Gatherv(senddata.dataPtr(),
                                nmtags[ParallelDescriptor::MyProc()],
                                ParallelDescriptor::Mpi_typemap<Real>::type(),
                                recvdata.dataPtr(),
                                nmtags.dataPtr(),
                                offset.dataPtr(),
                                ParallelDescriptor::Mpi_typemap<Real>::type(),
                                procToWrite,
                                ParallelDescriptor::Communicator()) );

    BL_COMM_PROFILE(BLProfiler::Gatherv, recvdata.size() * sizeof(Real),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    if(ParallelDescriptor::MyProc() == procToWrite) {
        for(int i(0), N(mf.size()); i < N; ++i) {
            if(pmap[i] != procToWrite) {
                m_min[i].resize(m_ncomp);
                m_max[i].resize(m_ncomp);
            }
        }

        for(int j(0), N(mf.size()); j < N; ++j) {
            if(pmap[j] != procToWrite) {
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
    if(ParallelDescriptor::MyProc() == procToWrite) {
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
    for(int comp(0); comp < m_ncomp; ++comp) {
      m_famin[comp] =  std::numeric_limits<Real>::max();
      m_famax[comp] = -std::numeric_limits<Real>::max();
    }

    for(int ibox(0); ibox < m_min.size(); ++ibox) {
      for(int comp(0); comp < m_min[ibox].size(); ++comp) {
        m_famin[comp] = std::min(m_famin[comp], m_min[ibox][comp]);
      }
    }
    for(int ibox(0); ibox < m_max.size(); ++ibox) {
      for(int comp(0); comp < m_max[ibox].size(); ++comp) {
        m_famax[comp] = std::max(m_famax[comp], m_max[ibox][comp]);
      }
    }
}


long
VisMF::WriteHeader (const std::string &mf_name,
                    VisMF::Header     &hdr,
		    int                procToWrite)
{
    BL_PROFILE("VisMF::WriteHeader");
    long bytesWritten(0);

    if(ParallelDescriptor::MyProc() == procToWrite) {
        std::string MFHdrFileName(mf_name);

        MFHdrFileName += TheMultiFabHdrFileSuffix;

        VisMF::IO_Buffer io_buffer(ioBufferSize);

        std::ofstream MFHdrFile;

        MFHdrFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        MFHdrFile.open(MFHdrFileName.c_str(), std::ios::out | std::ios::trunc);

        if( ! MFHdrFile.good()) {
            amrex::FileOpenFailed(MFHdrFileName);
	}

        MFHdrFile << hdr;

        //
        // Add in the number of bytes written out in the Header.
        //
        bytesWritten += VisMF::FileOffset(MFHdrFile);

        MFHdrFile.flush();
        MFHdrFile.close();

	if(checkFilePositions) {
          std::stringstream hss;
	  hss << hdr;
	  if(hss.tellp() != bytesWritten) {
	    std::cerr << "**** tellp error: hss.tellp() != bytesWritten :  "
	              << hss.tellp() << "  " << bytesWritten << std::endl;
	  }
	}
	
    }
    return bytesWritten;
}


long
VisMF::Write (const FabArray<FArrayBox>&    mf,
              const std::string& mf_name,
              VisMF::How         how,
              bool               set_ghost)
{
    BL_PROFILE("VisMF::Write(FabArray)");
    BL_ASSERT(mf_name[mf_name.length() - 1] != '/');
    BL_ASSERT(currentVersion != VisMF::Header::Undefined_v1);

    // ---- add stream retry
    // ---- add stream buffer (to nfiles)
    RealDescriptor *whichRD;
    if(FArrayBox::getFormat() == FABio::FAB_NATIVE) {
      whichRD = FPC::NativeRealDescriptor().clone();
    } else if(FArrayBox::getFormat() == FABio::FAB_NATIVE_32) {
      whichRD = FPC::Native32RealDescriptor().clone();
    } else if(FArrayBox::getFormat() == FABio::FAB_IEEE_32) {
      whichRD = FPC::Ieee32NormalRealDescriptor().clone();
    }
    bool doConvert(*whichRD != FPC::NativeRealDescriptor());

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

    // ---- check if mf has sparse data
    bool useSparseFPP(false);
    const Vector<int> &pmap = mf.DistributionMap().ProcessorMap();
    std::set<int> procsWithData;
    Vector<int> procsWithDataVector;
    for(int i(0); i < pmap.size(); ++i) {
      procsWithData.insert(pmap[i]);
    }
    if(allowSparseWrites && (procsWithData.size() < nOutFiles)) {
      useSparseFPP = true;
//      amrex::Print() << "SSSSSSSS:  in VisMF::Write:  useSparseFPP for:  " << mf_name << '\n';
      for(std::set<int>::iterator it = procsWithData.begin(); it != procsWithData.end(); ++it) {
        procsWithDataVector.push_back(*it);
      }
    }

    int coordinatorProc(ParallelDescriptor::IOProcessorNumber());
    long bytesWritten(0);
    bool calcMinMax(false);
    VisMF::Header hdr(mf, how, currentVersion, calcMinMax);

    std::string filePrefix(mf_name + FabFileSuffix);

    NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf);

    bool oldHeader(currentVersion == VisMF::Header::Version_v1);

      if(useSparseFPP) {
        nfi.SetSparseFPP(procsWithDataVector);
      } else if(useDynamicSetSelection) {
        nfi.SetDynamic();
      }
      for( ; nfi.ReadyToWrite(); ++nfi) {
	  // ---- find the total number of bytes including fab headers if needed
          const FABio &fio = FArrayBox::getFABio();
          int whichRDBytes(whichRD->numBytes()), nFABs(0);
          long writeDataItems(0), writeDataSize(0);
          for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
	    const FArrayBox &fab = mf[mfi];
	    if(oldHeader) {
	      std::stringstream hss;
	      fio.write_header(hss, fab, fab.nComp());
	      bytesWritten += hss.tellp();
	    }
	    bytesWritten += fab.box().numPts() * mf.nComp() * whichRDBytes;
	    ++nFABs;
	  }
	  char *allFabData(nullptr);
	  bool canCombineFABs(false);
	  if((nFABs > 1 || doConvert) && VisMF::useSingleWrite) {
	    allFabData = new(std::nothrow) char[bytesWritten];
	  }    // ---- else { no need to make a copy for one fab }
	  if(allFabData == nullptr) {
	    canCombineFABs = false;
	  } else {
	    canCombineFABs = true;
	  }

	  if(canCombineFABs) {
            long writePosition(0);
            for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
              int hLength(0);
              const FArrayBox &fab = mf[mfi];
	      writeDataItems = fab.box().numPts() * mf.nComp();
	      writeDataSize = writeDataItems * whichRDBytes;
	      char *afPtr = allFabData + writePosition;
	      if(oldHeader) {
	        std::stringstream hss;
	        fio.write_header(hss, fab, fab.nComp());
	        hLength = hss.tellp();
	        memcpy(afPtr, hss.str().c_str(), hLength);  // ---- the fab header
	      }
	      if(doConvert) {
	        RealDescriptor::convertFromNativeFormat(static_cast<void *> (afPtr + hLength),
		                                        writeDataItems,
		                                        fab.dataPtr(), *whichRD);
	      } else {    // ---- copy from the fab
	        memcpy(afPtr + hLength, fab.dataPtr(), writeDataSize);
	      }
              writePosition += hLength + writeDataSize;
            }
            nfi.Stream().write(allFabData, bytesWritten);
            nfi.Stream().flush();
	    delete [] allFabData;

	  } else {    // ---- write fabs individually
            for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
              int hLength(0);
              const FArrayBox &fab = mf[mfi];
	      writeDataItems = fab.box().numPts() * mf.nComp();
	      writeDataSize = writeDataItems * whichRDBytes;
	      if(oldHeader) {
	        std::stringstream hss;
	        fio.write_header(hss, fab, fab.nComp());
	        hLength = hss.tellp();
                nfi.Stream().write(hss.str().c_str(), hLength);    // ---- the fab header
                nfi.Stream().flush();
	      }
	      if(doConvert) {
	        char *cDataPtr = new char[writeDataSize];
	        RealDescriptor::convertFromNativeFormat(static_cast<void *> (cDataPtr),
		                                        writeDataItems,
		                                        fab.dataPtr(), *whichRD);
                nfi.Stream().write(cDataPtr, writeDataSize);
                nfi.Stream().flush();
	        delete [] cDataPtr;
	      } else {    // ---- copy from the fab
                nfi.Stream().write((char *) fab.dataPtr(), writeDataSize);
                nfi.Stream().flush();
	      }
            }
	  }
      }


    if(nfi.GetDynamic()) {
      coordinatorProc = nfi.CoordinatorProc();
    }

    if(currentVersion == VisMF::Header::Version_v1 ||
       currentVersion == VisMF::Header::NoFabHeaderMinMax_v1)
    {
      hdr.CalculateMinMax(mf, coordinatorProc);
    }

    VisMF::FindOffsets(mf, filePrefix, hdr, groupSets, currentVersion, nfi);

    bytesWritten += VisMF::WriteHeader(mf_name, hdr, coordinatorProc);

    delete whichRD;

    return bytesWritten;
}


void
VisMF::FindOffsets (const FabArray<FArrayBox> &mf,
		    const std::string &filePrefix,
                    VisMF::Header &hdr,
		    bool groupSets,
		    VisMF::Header::Version whichVersion,
		    NFilesIter &nfi)
{
    BL_PROFILE("VisMF::FindOffsets");

    const int myProc(ParallelDescriptor::MyProc());
    const int nProcs(ParallelDescriptor::NProcs());
    int coordinatorProc(ParallelDescriptor::IOProcessorNumber());
    if(nfi.GetDynamic()) {
      coordinatorProc = nfi.CoordinatorProc();
    }

    if(FArrayBox::getFormat() == FABio::FAB_ASCII ||
       FArrayBox::getFormat() == FABio::FAB_8BIT)
    {

#ifdef BL_USE_MPI
    Vector<int> nmtags(nProcs,0);
    Vector<int> offset(nProcs,0);

    const Vector<int> &pmap = mf.DistributionMap().ProcessorMap();

    for(int i(0), N(mf.size()); i < N; ++i) {
        ++nmtags[pmap[i]];
    }

    for(int i(1), N(offset.size()); i < N; ++i) {
        offset[i] = offset[i-1] + nmtags[i-1];
    }

    Vector<long> senddata(nmtags[myProc]);

    if(senddata.empty()) {
      // Can't let senddata be empty as senddata.dataPtr() will fail.
      senddata.resize(1);
    }

    int ioffset(0);

    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
      senddata[ioffset++] = hdr.m_fod[mfi.index()].m_head;
    }

    BL_ASSERT(ioffset == nmtags[myProc]);

    Vector<long> recvdata(mf.size());

    BL_COMM_PROFILE(BLProfiler::Gatherv, recvdata.size() * sizeof(long),
                    myProc, BLProfiler::BeforeCall());

    BL_MPI_REQUIRE( MPI_Gatherv(senddata.dataPtr(),
                                nmtags[myProc],
                                ParallelDescriptor::Mpi_typemap<long>::type(),
                                recvdata.dataPtr(),
                                nmtags.dataPtr(),
                                offset.dataPtr(),
                                ParallelDescriptor::Mpi_typemap<long>::type(),
                                coordinatorProc,
                                ParallelDescriptor::Communicator()) );

    BL_COMM_PROFILE(BLProfiler::Gatherv, recvdata.size() * sizeof(long),
                    myProc, BLProfiler::AfterCall());

    if(myProc == coordinatorProc) {
        Vector<int> cnt(nProcs,0);

        for(int j(0), N(mf.size()); j < N; ++j) {
            const int i(pmap[j]);
            hdr.m_fod[j].m_head = recvdata[offset[i]+cnt[i]];

            const std::string name(NFilesIter::FileName(nOutFiles, filePrefix, i, groupSets));

            hdr.m_fod[j].m_name = VisMF::BaseName(name);

            ++cnt[i];
        }
    }
#endif /*BL_USE_MPI*/


    } else {    // ---- calculate offsets

      RealDescriptor *whichRD;
      if(FArrayBox::getFormat() == FABio::FAB_NATIVE) {
        whichRD = FPC::NativeRealDescriptor().clone();
      } else if(FArrayBox::getFormat() == FABio::FAB_NATIVE_32) {
        whichRD = FPC::Native32RealDescriptor().clone();
      } else if(FArrayBox::getFormat() == FABio::FAB_IEEE_32) {
        whichRD = FPC::Ieee32NormalRealDescriptor().clone();
      }
      const FABio &fio = FArrayBox::getFABio();
      int whichRDBytes(whichRD->numBytes());
      int nComps(mf.nComp());

      if(myProc == coordinatorProc) {   // ---- calculate offsets
	const BoxArray &mfBA = mf.boxArray();
	const DistributionMapping &mfDM = mf.DistributionMap();
	Vector<long> fabHeaderBytes(mfBA.size(), 0);
	int nFiles(NFilesIter::ActualNFiles(nOutFiles));
	int whichFileNumber(-1);
	std::string whichFileName;
	Vector<long> currentOffset(nProcs, 0L);

        if(hdr.m_vers == VisMF::Header::Version_v1) {
	  // ---- find the length of the fab header instead of asking the file system
	  for(int i(0); i < mfBA.size(); ++i) {
            std::stringstream hss;
	    FArrayBox tempFab(mf.fabbox(i), nComps, false);  // ---- no alloc
            fio.write_header(hss, tempFab, tempFab.nComp());
	    fabHeaderBytes[i] = hss.tellp();
	  }
	}

	std::map<int, Vector<int> > rankBoxOrder;  // ---- [rank, boxarray index array]
	for(int i(0); i < mfBA.size(); ++i) {
	  rankBoxOrder[mfDM[i]].push_back(i);
	}

	Vector<int> fileNumbers;
        if(nfi.GetDynamic()) {
	  fileNumbers = nfi.FileNumbersWritten();
        }
         else if(nfi.GetSparseFPP()) {        // if sparse, write to (file number = rank)
 	  fileNumbers.resize(nProcs);
	  for(int i(0); i < nProcs; ++i) {
	    fileNumbers[i] = i;
          }
        }
        else {
	  fileNumbers.resize(nProcs);
	  for(int i(0); i < nProcs; ++i) {
	    fileNumbers[i] = NFilesIter::FileNumber(nFiles, i, groupSets);
	  }
	}

	const Vector< Vector<int> > &fileNumbersWriteOrder = nfi.FileNumbersWriteOrder();

	for(int fn(0); fn < fileNumbersWriteOrder.size(); ++fn) {
	  for(int ri(0); ri < fileNumbersWriteOrder[fn].size(); ++ri) {
	    int rank(fileNumbersWriteOrder[fn][ri]);
	    std::map<int, Vector<int> >::iterator rboIter = rankBoxOrder.find(rank);

            if(rboIter != rankBoxOrder.end()) {
              Vector<int> &index = rboIter->second;
	      whichFileNumber = fileNumbers[rank];
	      whichFileName   = VisMF::BaseName(NFilesIter::FileName(whichFileNumber, filePrefix));

	      for(int i(0); i < index.size(); ++i) {
                 hdr.m_fod[index[i]].m_name = whichFileName;
                 hdr.m_fod[index[i]].m_head = currentOffset[whichFileNumber];
                 currentOffset[whichFileNumber] += mf.fabbox(index[i]).numPts() * nComps * whichRDBytes
	                                           + fabHeaderBytes[index[i]];
              }
            } 
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
        int rv(std::remove(fileName.c_str()));
        if(verbose) {
          if(rv != 0) {
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

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(FullHdrFileName, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream infs(fileCharPtrString, std::istringstream::in);

    infs >> m_hdr;

    m_pa.resize(m_hdr.m_ncomp);

    for(int n(0); n < m_pa.size(); ++n) {
        m_pa[n].resize(m_hdr.m_ba.size());

        for(int ii(0), N(m_pa[n].size()); ii < N; ++ii) {
            m_pa[n][ii] = 0;
        }
    }
}


VisMF::~VisMF ()
{
}


FArrayBox*
VisMF::readFAB (int                  idx,
                const std::string   &mf_name,
                const VisMF::Header &hdr,
		int                  whichComp)
{
    BL_PROFILE("VisMF::readFAB_idx");
    Box fab_box(hdr.m_ba[idx]);
    if(hdr.m_ngrow) {
        fab_box.grow(hdr.m_ngrow);
    }

    FArrayBox *fab = new FArrayBox(fab_box, whichComp == -1 ? hdr.m_ncomp : 1);

    std::string FullName(VisMF::DirName(mf_name));
    FullName += hdr.m_fod[idx].m_name;

    std::ifstream *infs = VisMF::OpenStream(FullName);
    infs->seekg(hdr.m_fod[idx].m_head, std::ios::beg);

    if(hdr.m_vers == Header::Version_v1) {
      if(whichComp == -1) {    // ---- read all components
        fab->readFrom(*infs);
      } else {
        fab->readFrom(*infs, whichComp);
      }
    } else {
      if(whichComp == -1) {    // ---- read all components
	if(hdr.m_writtenRD == FPC::NativeRealDescriptor()) {
          infs->read((char *) fab->dataPtr(), fab->nBytes());
	} else {
          long readDataItems(fab->box().numPts() * fab->nComp());
          RealDescriptor::convertToNativeFormat(fab->dataPtr(), readDataItems,
	                                        *infs, hdr.m_writtenRD);
	}

      } else {
        long bytesPerComp(fab->box().numPts() * hdr.m_writtenRD.numBytes());
        infs->seekg(bytesPerComp * whichComp, std::ios::cur);
	if(hdr.m_writtenRD == FPC::NativeRealDescriptor()) {
          infs->read((char *) fab->dataPtr(), bytesPerComp);
	} else {
          long readDataItems(fab->box().numPts());  // ---- one component only
          RealDescriptor::convertToNativeFormat(fab->dataPtr(), readDataItems,
	                                        *infs, hdr.m_writtenRD);
	}
      }
    }

    VisMF::CloseStream(FullName);

    return fab;
}


void
VisMF::readFAB (FabArray<FArrayBox> &mf,
		int                  idx,
                const std::string&   mf_name,
                const VisMF::Header& hdr)
{
    BL_PROFILE("VisMF::readFAB_mf");
    FArrayBox &fab = mf[idx];

    std::string FullName(VisMF::DirName(mf_name));
    FullName += hdr.m_fod[idx].m_name;

    std::ifstream *infs = VisMF::OpenStream(FullName);
    infs->seekg(hdr.m_fod[idx].m_head, std::ios::beg);

    if(NoFabHeader(hdr)) {
      if(hdr.m_writtenRD == FPC::NativeRealDescriptor()) {
        infs->read((char *) fab.dataPtr(), fab.nBytes());
      } else {
        long readDataItems(fab.box().numPts() * fab.nComp());
        RealDescriptor::convertToNativeFormat(fab.dataPtr(), readDataItems,
	                                      *infs, hdr.m_writtenRD);
      }
    } else {
      fab.readFrom(*infs);
    }

    VisMF::CloseStream(FullName);
}


void
VisMF::Read (FabArray<FArrayBox> &mf,
             const std::string   &mf_name,
	     const char *faHeader,
	     int coordinatorProc)
{
    BL_PROFILE("VisMF::Read()");

    VisMF::Header hdr;
    Real hEndTime, hStartTime, faCopyTime(0.0);
    Real startTime(ParallelDescriptor::second());
    static Real totalTime(0.0);
    int myProc(ParallelDescriptor::MyProc());
    int messTotal(0);

    if(verbose && myProc == coordinatorProc) {
      std::cout << myProc << "::VisMF::Read:  about to read:  " << mf_name << std::endl;
    }

    std::string FullHdrFileName(mf_name + TheMultiFabHdrFileSuffix);

    {
        hStartTime = ParallelDescriptor::second();
        std::string fileCharPtrString;
	if(faHeader == nullptr) {
          Vector<char> fileCharPtr;
          ParallelDescriptor::ReadAndBcastFile(FullHdrFileName, fileCharPtr);
          fileCharPtrString = fileCharPtr.dataPtr();
	} else {
          fileCharPtrString = faHeader;
	}
        std::istringstream infs(fileCharPtrString, std::istringstream::in);

        infs >> hdr;

        hEndTime = ParallelDescriptor::second();
    }

    if (mf.empty()) {
	DistributionMapping dm(hdr.m_ba);
	mf.define(hdr.m_ba, dm, hdr.m_ncomp, hdr.m_ngrow, MFInfo(), FArrayBoxFactory());
    } else {
	BL_ASSERT(amrex::match(hdr.m_ba,mf.boxArray()));
    }

#ifdef BL_USE_MPI

  // ---- This limits the number of concurrent readers per file.
  int nOpensPerFile(nMFFileInStreams);
  int nProcs(ParallelDescriptor::NProcs());
  bool noFabHeader(NoFabHeader(hdr));

  if(noFabHeader && useSynchronousReads) {

    // ---- This code is only for reading in file order
    bool doConvert(hdr.m_writtenRD != FPC::NativeRealDescriptor());

    // ---- Create an ordered map of which processors read which
    // ---- Fabs in each file
    
    std::map<std::string, Vector<FabReadLink> > FileReadChains;        // ---- [filename, chain]
    std::map<std::string, std::set<int> > readFileRanks;              // ---- [filename, ranks]

    int nBoxes(hdr.m_ba.size());
    for(int i(0); i < nBoxes; ++i) {   // ---- create the map
      int undefined(-1);
      std::string fname(hdr.m_fod[i].m_name);
      FileReadChains[fname].push_back(FabReadLink(undefined, undefined, hdr.m_fod[i].m_head, hdr.m_ba[i]));
    }

    std::map<std::string, Vector<FabReadLink> >::iterator frcIter;

    int indexFileOrder(0);
    int currentRank(0);
    FabArray<FArrayBox> fafabFileOrder;
    BoxArray baFileOrder(hdr.m_ba.size());

    Vector<int> ranksFileOrder(mf.DistributionMap().size(), -1);

    Vector<int> nRanksPerFile(FileReadChains.size());
    amrex::NItemsPerBin(nProcs, nRanksPerFile);
    int currentFileIndex(0);

    for(frcIter = FileReadChains.begin(); frcIter != FileReadChains.end(); ++frcIter) {
      const std::string &fileName = frcIter->first;
      Vector<FabReadLink> &frc = frcIter->second;
      // ---- sort by offset
      std::sort(frc.begin(), frc.end(), [] (const FabReadLink &a, const FabReadLink &b)
	                                      { return a.fileOffset < b.fileOffset; } );

      Vector<int> nBoxesPerRank(nRanksPerFile[currentFileIndex]);
      amrex::NItemsPerBin(frc.size(), nBoxesPerRank);
      int frcIndex(0);

      for(int nbpr(0); nbpr < nBoxesPerRank.size(); ++nbpr) {
        for(int nb(0); nb < nBoxesPerRank[nbpr]; ++nb) {

	  baFileOrder.set(indexFileOrder, frc[frcIndex].box);
	  ranksFileOrder[indexFileOrder] = currentRank;
	  frc[frcIndex].rankToRead = currentRank;
	  frc[frcIndex].faIndex    = indexFileOrder;
	  readFileRanks[fileName].insert(currentRank);

	  ++frcIndex;
	  ++indexFileOrder;
        }
        ++currentRank;
        currentRank = std::min(currentRank, nProcs - 1);
      }
      ++currentFileIndex;
    }

    DistributionMapping dmFileOrder(std::move(ranksFileOrder));

    bool inFileOrder(mf.DistributionMap() == dmFileOrder && mf.boxArray() == baFileOrder);
    if(inFileOrder) {
      if(myProc == coordinatorProc && verbose) {
        std::cout << "VisMF::Read:  inFileOrder" << std::endl;
      }
    } else {
      if(myProc == coordinatorProc && verbose) {
        std::cout << "VisMF::Read:  not inFileOrder" << std::endl;
      }
      // ---- make a temporary fabarray in file order
      fafabFileOrder.define(baFileOrder, dmFileOrder, hdr.m_ncomp, hdr.m_ngrow, MFInfo(), mf.Factory());
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
      Vector<std::set<int> > streamSets(nStreams);
      int sIndex(0), sCount(0);
      for(setIter = rfrSplitSet.begin(); setIter != rfrSplitSet.end(); ++setIter) {
        streamSets[sIndex].insert(*setIter);
	if(++sCount >= ranksPerStream) {
	  sCount = 0;
	  ++sIndex;
	  sIndex = std::min<int>(sIndex, streamSets.size() - 1);
	}
      }

      for(int iSet(0); iSet < streamSets.size(); ++iSet) {
        Vector<int> readRanks;
        std::set<int> &rfrSet = streamSets[iSet];
        for(setIter = rfrSet.begin(); setIter != rfrSet.end(); ++setIter) {
          readRanks.push_back(*setIter);
        }

        if(rfrSet.find(myProc) != rfrSet.end()) {  // ---- myProc needs to read this file
          const std::string &fileName = rfrIter->first;
	  std::string fullFileName(VisMF::DirName(mf_name) + fileName);
	  frcIter = FileReadChains.find(fileName);
	  BL_ASSERT(frcIter != FileReadChains.end());
          Vector<FabReadLink> &frc = frcIter->second;
          for(NFilesIter nfi(fullFileName, readRanks); nfi.ReadyToRead(); ++nfi) {

	      // ---- confirm the data is contiguous in the stream
	      long firstOffset(-1);
	      for(int i(0); i < frc.size(); ++i) {
	        if(myProc == frc[i].rankToRead) {
		  firstOffset = frc[i].fileOffset;
		  break;
		}
	      }

	      bool dataIsContiguous(true);
	      long currentOffset(firstOffset), bytesToRead(0);
	      int nFABs(0);

	      for(int i(0); i < frc.size(); ++i) {
	        if(myProc == frc[i].rankToRead) {
	          if(currentOffset != frc[i].fileOffset) {
                    dataIsContiguous = false;
	          } else {
	            FArrayBox &fab = whichFA[frc[i].faIndex];
		    long fabBytesToRead(fab.box().numPts() * fab.nComp() * hdr.m_writtenRD.numBytes());
                    currentOffset += fabBytesToRead;
                    bytesToRead   += fabBytesToRead;
		    ++nFABs;
		  }
	        }
	      }
	      char *allFabData;
	      bool canCombineFABs(false);
	      if(nFABs > 1 && dataIsContiguous && VisMF::useSingleRead) {
	        allFabData = new(std::nothrow) char[bytesToRead];
		if(allFabData == nullptr) {
		  canCombineFABs = false;
		} else {
		  canCombineFABs = true;
		}
	      }
	      if(canCombineFABs) {
                nfi.Stream().seekp(firstOffset, std::ios::beg);
                nfi.Stream().read(allFabData, bytesToRead);

		currentOffset = 0;  // ---- this is now relative to allFabData

	        for(int i(0); i < frc.size(); ++i) {
	          if(myProc == frc[i].rankToRead) {
		    char *afPtr = allFabData + currentOffset;
	            FArrayBox &fab = whichFA[frc[i].faIndex];
		    long readDataItems(fab.box().numPts() * fab.nComp());
		    if(doConvert) {
		      RealDescriptor::convertToNativeFormat(fab.dataPtr(), readDataItems,
		                                            afPtr, hdr.m_writtenRD);
		    } else {
                      memcpy(fab.dataPtr(), afPtr, fab.nBytes());
		    }
                    currentOffset += readDataItems * hdr.m_writtenRD.numBytes();
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
		    long readDataItems(fab.box().numPts() * fab.nComp());
		    if(doConvert) {
		      RealDescriptor::convertToNativeFormat(fab.dataPtr(), readDataItems,
		                                            nfi.Stream(), hdr.m_writtenRD);
		    } else {
                      nfi.Stream().read((char *) fab.dataPtr(), fab.nBytes());
		    }
	          }
	        }
	      }

          }    // ---- end NFilesIter
        }

      }
    }

    if( ! inFileOrder) {
      faCopyTime = ParallelDescriptor::second();
      mf.copy(fafabFileOrder);
      faCopyTime = ParallelDescriptor::second() - faCopyTime;
    }

  } else {    // ---- (noFabHeader && useSynchronousReads) == false

    int nReqs(0), ioProcNum(coordinatorProc);
    int nBoxes(hdr.m_ba.size());
    int totalIOReqs(nBoxes), nFiles(-1);
    std::vector<int> iDone(2);
    const int iDoneIndex(0), iDoneCount(1);
    std::set<int> busyProcs;  // [whichProc]
    std::map<std::string, int> fileNames;  // <filename, allreadsindex>
    std::multiset<int> availableFiles;  // [whichFile]  supports multiple reads/file
    int allReadsIndex(0);
    ParallelDescriptor::Message rmess;
    Vector<std::map<int,std::map<long,int> > > allReads; // [file]<proc,<seek,index>>


    for(int i(0); i < nBoxes; ++i) {   // count the files
      int whichProc(mf.DistributionMap()[i]);
      if(whichProc == myProc) {
        ++nReqs;
      }
      if(myProc == coordinatorProc) {
        std::string fname(hdr.m_fod[i].m_name);
	if(fileNames.insert(std::pair<std::string,int>(fname,allReadsIndex)).second)
	{
	  ++allReadsIndex;
	}
      }
    }

    if(myProc == coordinatorProc) {    // fill availableFiles
      nFiles = fileNames.size();
      for(int i(0); i < nFiles; ++i) {
        for(int nOpens(0); nOpens < nOpensPerFile; ++nOpens) {
          availableFiles.insert(i);
        }
      }
      allReads.resize(nFiles);
      int whichProc;
      long iSeekPos;
      std::map<std::string, int>::iterator fileNamesIter;
      for(int i(0); i < nBoxes; ++i) {   // fill allReads maps
        whichProc = mf.DistributionMap()[i];
        iSeekPos = hdr.m_fod[i].m_head;
        std::string fname(hdr.m_fod[i].m_name);
	fileNamesIter = fileNames.find(fname);
	if(fileNamesIter != fileNames.end()) {
	  int findex(fileNames.find(fname)->second);
	  allReads[findex][whichProc].insert(std::pair<long, int>(iSeekPos, i));
	} else {
	  std::cout << "**** Error:  filename not found = " << fname << std::endl;
	  amrex::Abort("**** Error in VisMF::Read");
	}
      }
    }

    int readTag(ParallelDescriptor::SeqNum());
    int doneTag(ParallelDescriptor::SeqNum());

    if(myProc == coordinatorProc) {  // manage the file locks
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
          std::map<int,std::map<long,int> >::iterator whichRead;
	  for(whichRead = allReads[arIndex].begin();
	      whichRead != allReads[arIndex].end(); ++whichRead)
	  {
	    int tryProc(whichRead->first);
	    if(busyProcs.find(tryProc) == busyProcs.end()) {  // tryProc not busy
	      busyProcs.insert(tryProc);
	      int nReads(whichRead->second.size());
	      int ir(0);
	      vReads.resize(nReads);
              std::map<long,int>::iterator imiter;
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
	        ParallelDescriptor::Send(vReads, tryProc, readTag);
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
	  VisMF::readFAB(mf,index, mf_name, hdr);
	  --totalIOReqs;
	  iopReads.pop_front();
	  if(iopReads.empty()) {
	    availableFiles.insert(iopFileIndex);
	    busyProcs.erase(ioProcNum);
	  }
	  ParallelDescriptor::IProbe(MPI_ANY_SOURCE, doneTag, doneFlag, status);
	  if(doneFlag) {
	    break;
	  }
	}

	if(reqsPending > 0) {
          rmess = ParallelDescriptor::Recv(iDone, MPI_ANY_SOURCE, doneTag);

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
      std::vector<int> recReads(nReqs, -1);
      while(nReqs > 0) {
        rmess = ParallelDescriptor::Recv(recReads, ioProcNum, readTag);
        for(int ir(0); ir < static_cast<int>(rmess.count()); ++ir) {
	  int mfIndex(recReads[ir]);
	  VisMF::readFAB(mf,mfIndex, mf_name, hdr);
	}
        nReqs -= rmess.count();
	iDone[iDoneIndex] = recReads[0];
	iDone[iDoneCount] = rmess.count();
        ParallelDescriptor::Send(iDone, ioProcNum, doneTag);
      }
    }

  }

#else
    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
      VisMF::readFAB(mf,mfi.index(), mf_name, hdr);
    }
#endif

    if(VisMF::GetUsePersistentIFStreams()) {
      for(int idx(0); idx < hdr.m_fod.size(); ++idx) {
        std::string FullName(VisMF::DirName(mf_name));
        FullName += hdr.m_fod[idx].m_name;
        VisMF::DeleteStream(FullName);
      }
    }

    if(myProc == coordinatorProc && verbose) {
      Real mfReadTime = ParallelDescriptor::second() - startTime;
      totalTime += mfReadTime;
      std::cout << "FARead ::  nBoxes = " << hdr.m_ba.size()
                << "  nMessages = " << messTotal << '\n';
      std::cout << "FARead ::  hTime = " << (hEndTime - hStartTime) << '\n';
      std::cout << "FARead ::  faCopyTime = " << faCopyTime << '\n';
      std::cout << "FARead ::  mfReadTime = " << mfReadTime
                << "  totalTime = " << totalTime << std::endl;
    }

    BL_ASSERT(mf.ok());
}


bool
VisMF::Exist (const std::string& mf_name)
{
    std::string FullHdrFileName(mf_name + TheMultiFabHdrFileSuffix);
    int exist;
    if (ParallelDescriptor::IOProcessor()) {
        std::ifstream iss;
        iss.open(FullHdrFileName.c_str(), std::ios::in);
        exist = iss.good();
    }
    ParallelDescriptor::Bcast(&exist, 1, ParallelDescriptor::IOProcessorNumber());
    return exist;
}

void
VisMF::ReadFAHeader (const std::string &fafabName,
	             Vector<char> &faHeader)
{
    BL_PROFILE("VisMF::ReadFAHeader()");

    std::string FullHdrFileName(fafabName + TheMultiFabHdrFileSuffix);
    ParallelDescriptor::ReadAndBcastFile(FullHdrFileName, faHeader);
}


bool
VisMF::Check (const std::string& mf_name)
{
  BL_PROFILE("VisMF::Check()");

  int isOk(true);  // ---- int to broadcast
  int v1(true);

  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "---------------- VisMF::Check:  about to check:  " << mf_name << std::endl;

    char c;
    int nBadFabs(0);
    VisMF::Header hdr;
    std::string FullHdrFileName(mf_name);
    FullHdrFileName += TheMultiFabHdrFileSuffix;

    {
        std::ifstream ifs(FullHdrFileName.c_str());
        ifs >> hdr;
        ifs.close();
    }

    std::cout << "hdr.version =  " << hdr.m_vers << std::endl;
    std::cout << "hdr.boxarray size =  " << hdr.m_ba.size() << std::endl;
    std::cout << "mf.ncomp =  " << hdr.m_ncomp << std::endl;
    std::cout << "number of fabs on disk =  " << hdr.m_fod.size() << std::endl;
    std::cout << "DirName = " << DirName(mf_name) << std::endl;
    std::cout << "mf_name = " << mf_name << std::endl;
    std::cout << "FullHdrFileName = " << FullHdrFileName << std::endl;

    if(hdr.m_vers != VisMF::Header::Version_v1) {
     v1 = false;
     std::cout << "**** VisMF::Check currently only supports Version_v1." << std::endl;
    } else {

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
	          << i << "  seekpos = " << fod.m_head << "  box = " << hdr.m_ba[i]
	          << std::endl;
      }
      ifs.close();

    }
    if(nBadFabs) {
      std::cout << "Total Bad Fabs = " << nBadFabs << std::endl;
      isOk = false;
    } else {
      std::cout << "No Bad Fabs." << std::endl;
      isOk = true;
    }
    }

  }
  ParallelDescriptor::Bcast(&isOk, 1, ParallelDescriptor::IOProcessorNumber());
  ParallelDescriptor::Bcast(&v1,  1, ParallelDescriptor::IOProcessorNumber());

  return isOk;

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


bool VisMF::NoFabHeader(const VisMF::Header &hdr) {
  if(hdr.m_vers == VisMF::Header::NoFabHeader_v1       ||
    hdr.m_vers == VisMF::Header::NoFabHeaderMinMax_v1 ||
    hdr.m_vers == VisMF::Header::NoFabHeaderFAMinMax_v1)
  {
    return true;
  }
  return false;
}


VisMF::PersistentIFStream::PersistentIFStream()
    :
    pstr(0),
    currentPosition(0),
    isOpen(false)
{ }


VisMF::PersistentIFStream::~PersistentIFStream()
{
  if(isOpen) {
    pstr->close();
    delete pstr;
    pstr = 0;
    isOpen = false;
  }
}


std::ifstream *VisMF::OpenStream(const std::string &fileName) {
  VisMF::PersistentIFStream &pifs = VisMF::persistentIFStreams[fileName];
  if( ! pifs.isOpen) {
    pifs.pstr = new std::ifstream;
    pifs.pstr->open(fileName.c_str(), std::ios::in | std::ios::binary);
    if( ! pifs.pstr->good()) {
      delete pifs.pstr;
      amrex::FileOpenFailed(fileName);
    }
    pifs.isOpen = true;
    pifs.currentPosition = 0;
    if(setBuf) {
      pifs.ioBuffer.resize(ioBufferSize);
      pifs.pstr->rdbuf()->pubsetbuf(pifs.ioBuffer.dataPtr(), pifs.ioBuffer.size());
    }
  }

  return pifs.pstr;
}


void VisMF::CloseStream(const std::string &fileName, bool forceClose)
{
  if(usePersistentIFStreams && ! forceClose) {
    return;
  }

  VisMF::PersistentIFStream &pifs = VisMF::persistentIFStreams[fileName];
  if(pifs.isOpen) {
    pifs.pstr->close();
    delete pifs.pstr;
    pifs.pstr = 0;
    pifs.isOpen = false;
  }
  pifs.ioBuffer.clear();
}


void VisMF::DeleteStream(const std::string &fileName)
{
  if(usePersistentIFStreams) {
    auto psIter = VisMF::persistentIFStreams.find(fileName);
    if(psIter != VisMF::persistentIFStreams.end()) {
      VisMF::persistentIFStreams.erase(psIter);
    }
  }
}


void VisMF::CloseAllStreams() {
  VisMF::persistentIFStreams.clear();
}

}
