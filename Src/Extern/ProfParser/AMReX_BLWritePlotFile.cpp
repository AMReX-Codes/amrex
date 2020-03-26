// --------------------------------------------------------------------
// AMReX_BLWritePlotFile.cpp
// --------------------------------------------------------------------
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <unistd.h>

#include <AMReX_BLWritePlotFile.H>
#include <AMReX_FPC.H>

namespace amrex{

// --------------------------------------------------------------------
void WritePlotfile(const std::string         &pfversion,
                   const Vector<MultiFab>     &data,
                   const Real                 time,
                   const Vector<Real>         &probLo,
                   const Vector<Real>         &probHi,
                   const Vector<int>          &refRatio,
                   const Vector<Box>          &probDomain,
                   const Vector<Vector<Real> > &dxLevel,
                   const int                  coordSys,
                   const std::string         &oFile,
                   const Vector<std::string>  &names,
                   const bool                 verbose,
                   const bool                 isCartGrid,
                   const Real                *vfeps,
		   const int                 *levelSteps)
{
    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(oFile,0755)) {
         amrex::CreateDirectoryFailed(oFile);
      }
    }
    // Force other processors to wait untill directory is built.
    ParallelDescriptor::Barrier();

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream os;
    const int finestLevel(data.size() - 1);

    if(ParallelDescriptor::IOProcessor()) {
        
    std::string oFileHeader(oFile);
    oFileHeader += "/Header";
    
    //os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    
    if(verbose && ParallelDescriptor::IOProcessor()) {
      std::cout << "Opening file = " << oFileHeader << '\n';
    }
    
    os.open(oFileHeader.c_str(), std::ios::out|std::ios::binary);
    
    if(os.fail()) {
      amrex::FileOpenFailed(oFileHeader);
    }
    //
    // Start writing plotfile.
    //
    os << pfversion << '\n';
    int n_var = data[0].nComp();
    os << n_var << '\n';
    for (int n = 0; n < n_var; n++) os << names[n] << '\n';
    os << BL_SPACEDIM << '\n';
    os << std::setprecision(30) << time << '\n';
    os << finestLevel << '\n';
    for (int i = 0; i < BL_SPACEDIM; i++) os << probLo[i] << ' ';
    os << '\n';
    for (int i = 0; i < BL_SPACEDIM; i++) os << probHi[i] << ' ';
    os << '\n';
    for (int i = 0; i < finestLevel; i++) os << refRatio[i] << ' ';
    os << '\n';
    for (int i = 0; i <= finestLevel; i++) os << probDomain[i] << ' ';
    os << '\n';
    if(levelSteps != 0) {
      for (int i = 0; i <= finestLevel; i++) os << levelSteps[i] << ' ';
    } else {
      for (int i = 0; i <= finestLevel; i++) os << 0 << ' ';
    }
    os << '\n';
    for(int i = 0; i <= finestLevel; i++) {
      for(int k = 0; k < BL_SPACEDIM; k++) {
            os << dxLevel[i][k] << ' ';
      }
      os << '\n';
    }
    if(isCartGrid) {
      for(int i(0); i <= finestLevel; i++) {
        os << vfeps[i] << ' ';
      }
      os << '\n';
    }
    os << coordSys << '\n';
    os << 0 << '\n';        // --------------- The bndry data width.

    }

    //
    // Write out level by level.
    //
    for(int iLevel(0); iLevel <= finestLevel; ++iLevel) {
        // Write state data.
        const BoxArray &ba = data[iLevel].boxArray();
        int nGrids = ba.size();
        char buf[64];
        sprintf(buf, "Level_%d", iLevel);
        
        if(ParallelDescriptor::IOProcessor()) {
            os << iLevel << ' ' << nGrids << ' ' << time << '\n';
            if(levelSteps != 0) {
              os << levelSteps[iLevel] << '\n';
	    } else {
              os << 0 << '\n';
	    }
            
            for(int i(0); i < nGrids; ++i) {
              const Box &b = ba[i];
              for(int n(0); n < BL_SPACEDIM; ++n) {
                Real glo = b.smallEnd()[n] * dxLevel[iLevel][n];
                Real ghi = (b.bigEnd()[n]+1) * dxLevel[iLevel][n];
                os << glo << ' ' << ghi << '\n';
              }
            }
            // Build the directory to hold the MultiFabs at this level.
            std::string Level(oFile);
            Level += '/';
            Level += buf;
            
            if( ! amrex::UtilCreateDirectory(Level, 0755)) {
              amrex::CreateDirectoryFailed(Level);
	    }
        }
        // Force other processors to wait till directory is built.
        ParallelDescriptor::Barrier();
        // Now build the full relative pathname of the MultiFab.
        static const std::string MultiFabBaseName("MultiFab");
        
        std::string PathName(oFile);
        PathName += '/';
        PathName += buf;
        PathName += '/';
        PathName += MultiFabBaseName;
        
        if(ParallelDescriptor::IOProcessor()) {
            //
            // The full name relative to the Header file.
            //
            std::string RelativePathName(buf);
            RelativePathName += '/';
            RelativePathName += MultiFabBaseName;
            os << RelativePathName << '\n';
        }
        VisMF::Write(data[iLevel], PathName);
    }
    
    if(ParallelDescriptor::IOProcessor()) {
      os.close();
    }
}



// --------------------------------------------------------------------
std::string VisMFBaseName(const std::string& filename) {
    BL_ASSERT(filename[filename.length() - 1] != '/');
    if(const char *slash = strrchr(filename.c_str(), '/')) {
        return std::string(slash + 1);
    } else {
        return filename;
    }
}


// --------------------------------------------------------------------
void Write2DBoxFrom3D(const Box &box, std::ostream &os, int whichPlane) {
  os << '(';
  switch(whichPlane) {
    case 0:
      os << '(' << box.smallEnd(1) << ',' << box.smallEnd(2) << ')' << ' '
         << '(' << box.bigEnd(1)   << ',' << box.bigEnd(2)   << ')' << ' '
         << '(' << box.type(1)     << ',' << box.type(2)     << ')';
    break;

    case 1:
      os << '(' << box.smallEnd(0) << ',' << box.smallEnd(2) << ')' << ' '
         << '(' << box.bigEnd(0)   << ',' << box.bigEnd(2)   << ')' << ' '
         << '(' << box.type(0)     << ',' << box.type(2)     << ')';
    break;

    case 2:
      os << '(' << box.smallEnd(0) << ',' << box.smallEnd(1) << ')' << ' '
         << '(' << box.bigEnd(0)   << ',' << box.bigEnd(1)   << ')' << ' '
         << '(' << box.type(0)     << ',' << box.type(1)     << ')';
    break;

    case 3:  // for testing
      os << '(' << box.smallEnd(0) << ',' << box.smallEnd(1) << ',' << box.smallEnd(2) << ')' << ' '
         << '(' << box.bigEnd(0)   << ',' << box.bigEnd(1)   << ',' << box.bigEnd(2)   << ')' << ' '
         << '(' << box.type(0)     << ',' << box.type(1)     << ',' << box.type(2)     << ')';
    break;
  }
  os << ')';

}


#include <AMReX_FabConv.H>
#include <AMReX_FPC.H>
// --------------------------------------------------------------------
VisMF::FabOnDisk VisMFWrite(const FArrayBox &fabIn, const std::string &filename,
                            std::ostream &os, long &bytes, int whichPlane)
{
    FArrayBox fab;
    if(fabIn.box().length(whichPlane) > 1) {
      Box b(fabIn.box());
      b.setBig(whichPlane, b.smallEnd(whichPlane));
      fab.resize(b, fabIn.nComp());
    } else {
      fab.resize(fabIn.box(), fabIn.nComp());
    }
    fab.copy<RunOn::Host>(fabIn);


    VisMF::FabOnDisk fab_on_disk(filename, VisMF::FileOffset(os));
// ====================
    //fab.writeOn(os);
    // ====================
    //fabio->write_header(os, *this, num_comp);
    //os << "FAB " << *rd;
    //os << fab.box() << ' ' << fab.nComp() << '\n';
    // ====================

    // ====================
    //fabio->write(os, *this, 0, fab.nComp());
    const long base_siz  = fab.box().numPts();
    const Real *comp_ptr = fab.dataPtr(0);
    const long siz       = base_siz*fab.nComp();
    const int *ord = FPC::reverse_double_order;
    RealDescriptor *rd = new RealDescriptor(FPC::ieee_double, ord, 8);
      os << "FAB " << *rd;
      //os << fab.box();
      Write2DBoxFrom3D(fab.box(), os, whichPlane);
      os << ' ' << fab.nComp() << '\n';
    RealDescriptor::convertFromNativeFormat(os, siz, comp_ptr, *rd);
    //os.write((char *) comp_ptr, siz*sizeof(double));
    // ====================
// ====================
    bytes += (VisMF::FileOffset(os) - fab_on_disk.m_head);
    return fab_on_disk;
}


// --------------------------------------------------------------------
static std::ostream &operator<<(std::ostream &os,
                                const Vector< Vector<Real> > &ar)
{
    long i = 0, N = ar.size(), M = (N == 0) ? 0 : ar[0].size();
    os << N << ',' << M << '\n';
    for( ; i < N; i++) {
      BL_ASSERT(ar[i].size() == M);
      for(long j = 0; j < M; j++) {
        os << ar[i][j] << ',';
      }
      os << '\n';
    }
    if( ! os.good()) {
      Error("Write of Vector<Vector<Real>> failed");
    }
    return os;
}


// --------------------------------------------------------------------
long VisMFWriteHeader(const std::string &mf_name, VisMF::Header &hdr,
                      int whichPlane)
{
  long bytes(0);
  std::string MFHdrFileName(mf_name);
  MFHdrFileName += "_H";
  VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
  std::ofstream MFHdrFile;
  MFHdrFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
  MFHdrFile.open(MFHdrFileName.c_str(), std::ios::out|std::ios::trunc);
  if( ! MFHdrFile.good()) {
    FileOpenFailed(MFHdrFileName);
  }
  // ===============================
  //MFHdrFile << hdr;
    //std::ios::fmtflags oflags = MFHdrFile.flags();
    MFHdrFile.setf(std::ios::floatfield, std::ios::scientific);
    //int old_prec = MFHdrFile.precision(15);

  MFHdrFile << hdr.m_vers << '\n';
  MFHdrFile << int(hdr.m_how) << '\n';
  MFHdrFile << hdr.m_ncomp    << '\n';
  MFHdrFile << hdr.m_ngrow    << '\n'; ;

  MFHdrFile << '(' << hdr.m_ba.size() << ' ' << 0 << '\n';

  for(int i(0); i < hdr.m_ba.size(); ++i) {
    Write2DBoxFrom3D(hdr.m_ba[i], MFHdrFile, whichPlane);
    MFHdrFile << '\n';
  }

  MFHdrFile << ')' << '\n';


  MFHdrFile << hdr.m_fod << '\n';
  MFHdrFile << hdr.m_min << '\n';
  MFHdrFile << hdr.m_max << '\n';

  // ===============================
  bytes += VisMF::FileOffset(MFHdrFile);
  MFHdrFile.close();
  return bytes;
}


// --------------------------------------------------------------------
void WritePlotfile2DFrom3D(const std::string &pfversion,
                           const Vector<MultiFab>     &data,
                           const Real                 time,
                           const Vector<Real>         &probLo,
                           const Vector<Real>         &probHi,
                           const Vector<int>          &refRatio,
                           const Vector<Box>          &probDomain,
                           const Vector<Vector<Real> > &dxLevel,
                           const int                  coordSys,
                           const std::string         &oFile,
                           const Vector<std::string>  &names,
                           const bool                 verbose,
		           const bool                 isCartGrid,
		           const Real                *vfeps,
		           const int                 *levelSteps)
{

    int whichPlane;
    if(probDomain[0].length(0) == 1) {
      whichPlane = 0;
    } else if(probDomain[0].length(1) == 1) {
      whichPlane = 1;
    } else if(probDomain[0].length(2) == 1) {
      whichPlane = 2;
    } else {
      Abort("Error:  no length = 1 plane.");
    }

    if(ParallelDescriptor::IOProcessor()) {
      if( ! UtilCreateDirectory(oFile,0755)) {
         CreateDirectoryFailed(oFile);
      }
    }
    // Force other processors to wait till directory is built.
    ParallelDescriptor::Barrier();
    
    std::string oFileHeader(oFile);
    oFileHeader += "/Header";
    
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    
    std::ofstream os;
    
    //os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    
    if(verbose && ParallelDescriptor::IOProcessor()) {
      std::cout << "Opening file = " << oFileHeader << '\n';
    }
    
    os.open(oFileHeader.c_str(), std::ios::out|std::ios::binary);
    
    if(os.fail()) {
      FileOpenFailed(oFileHeader);
    }
    //
    // Start writing plotfile.
    //
    int spacedim(2);
    int bl_spacedim(BL_SPACEDIM);
    os << pfversion << '\n';
    int n_var = data[0].nComp();
    os << n_var << '\n';
    for (int n = 0; n < n_var; n++) os << names[n] << '\n';
    os << spacedim << '\n';
    os << std::setprecision(30) << time << '\n';
    const int finestLevel = data.size() - 1;
    os << finestLevel << '\n';
    for (int i = 0; i < bl_spacedim; i++) {
      if(i != whichPlane) {
        os << probLo[i] << ' ';
      }
    }
    os << '\n';
    for (int i = 0; i < bl_spacedim; i++) {
      if(i != whichPlane) {
        os << probHi[i] << ' ';
      }
    }
    os << '\n';
    for (int i = 0; i < finestLevel; i++) os << refRatio[i] << ' ';
    os << '\n';
    for (int i = 0; i <= finestLevel; i++) {
      //os << probDomain[i] << ' ';
      Write2DBoxFrom3D(probDomain[i], os, whichPlane);
      os << ' ';
    }
    os << '\n';
    if(levelSteps != 0) {
      for (int i = 0; i <= finestLevel; i++) os << levelSteps[i] << ' ';
    } else {
      for (int i = 0; i <= finestLevel; i++) os << 0 << ' ';
    }
    os << '\n';
    for(int i = 0; i <= finestLevel; i++) {
      for(int k = 0; k < bl_spacedim; k++) {
        if(k != whichPlane) {
          os << dxLevel[i][k] << ' ';
	}
      }
      os << '\n';
    }
    if(isCartGrid) {
      for(int i(0); i <= finestLevel; i++) {
        os << vfeps[i] << ' ';
      }
      os << '\n';
    }
    os << coordSys << '\n';
    os << 0 << '\n';                  // --------------- The bndry data width.
    //
    // Write out level by level.
    //
    for(int iLevel(0); iLevel <= finestLevel; ++iLevel) {
        //
        // Write state data.
        //
        const BoxArray &ba = data[iLevel].boxArray();
        int nGrids = ba.size();
        char buf[64];
        sprintf(buf, "Level_%d", iLevel);
        
        if(ParallelDescriptor::IOProcessor()) {
            os << iLevel << ' ' << nGrids << ' ' << time << '\n';
            if(levelSteps != 0) {
              os << levelSteps[iLevel] << '\n';
	    } else {
              os << 0 << '\n';
	    }
            
            for(int i(0); i < nGrids; ++i) {
              const Box &b = ba[i];
              for(int n(0); n < bl_spacedim; ++n) {
                if(n != whichPlane) {
                  Real glo = b.smallEnd()[n] * dxLevel[iLevel][n];
                  Real ghi = (b.bigEnd()[n]+1) * dxLevel[iLevel][n];
                  os << glo << ' ' << ghi << '\n';
		}
              }
            }
            //
            // Build the directory to hold the MultiFabs at this level.
            //
            std::string Level(oFile);
            Level += '/';
            Level += buf;
            
            if( ! UtilCreateDirectory(Level, 0755)) {
              CreateDirectoryFailed(Level);
	    }
        }
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        //
        // Now build the full relative pathname of the MultiFab.
        //
        static const std::string MultiFabBaseName("MultiFab");
        
        std::string PathName(oFile);
        PathName += '/';
        PathName += buf;
        PathName += '/';
        PathName += MultiFabBaseName;
        
        if(ParallelDescriptor::IOProcessor()) {
            // The full name relative to the Header file.
            std::string RelativePathName(buf);
            RelativePathName += '/';
            RelativePathName += MultiFabBaseName;
            os << RelativePathName << '\n';
        }
// ======================================
        //VisMF::Write(data[iLevel], PathName);

        string mf_name(PathName);
        VisMF::How how(VisMF::OneFilePerCPU);
        const MultiFab &mf = data[iLevel];

    BL_ASSERT(mf_name[mf_name.length() - 1] != '/');

    static const char* FabFileSuffix = "_D_";

    VisMF::Initialize();
    VisMF::Header hdr(mf, how);

    long        bytes    = 0;
    std::string FullName = Concatenate(mf_name + FabFileSuffix, 0, 4);

    const std::string BName = VisMFBaseName(FullName);

    std::ofstream FabFile;
    FabFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    FabFile.open(FullName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);
    if (!FabFile.good())
      FileOpenFailed(FullName);

    for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
      hdr.m_fod[mfi.index()] = VisMFWrite(mf[mfi],BName,FabFile,bytes, whichPlane);
    }
    FabFile.flush();
    FabFile.close();

    bytes += VisMFWriteHeader(mf_name, hdr, whichPlane);

// ======================================

    }
    
    os.close();
}
// --------------------------------------------------------------------
// --------------------------------------------------------------------
}
