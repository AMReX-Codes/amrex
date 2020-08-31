// ----------------------------------------------------------------------
//  AMReX_BLProfUtilities.cpp
// ----------------------------------------------------------------------
#ifndef BL_BLPROFUTILITIES_CPP
#define BL_BLPROFUTILITIES_CPP

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <cstring>
#include <list>
#include <vector>
#include <typeinfo>
#include <limits>
#include <algorithm>
#include <sys/stat.h>

using std::cout;
using std::endl;
using std::flush;
using std::string;

#include <AMReX_SPACE.H>
#include <AMReX.H>
#include <AMReX_BLProfStats.H>
#include <AMReX_CommProfStats.H>
#include <AMReX_RegionsProfStats.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_AVGDOWN_F.H>
#include <AMReX_BLProfUtilities.H>

const int maxFabXDim(1024);

#define SHOWVAL(val) { cout << #val << " = " << val << endl; }

// --------------------------------------------------------------------
std::string amrex::SanitizeName(const std::string &sname) {
  std::string s(sname);
  std::size_t found;
  Vector<char> replaceChars;
  replaceChars.push_back(' ');
  replaceChars.push_back(':');
  replaceChars.push_back('(');
  replaceChars.push_back(')');
  replaceChars.push_back('[');
  replaceChars.push_back(']');
  replaceChars.push_back('}');
  replaceChars.push_back('{');
  replaceChars.push_back('<');
  replaceChars.push_back('>');
  replaceChars.push_back(',');
  for(int c(0); c < replaceChars.size(); ++c) {
    while((found = s.find(replaceChars[c])) != std::string::npos) {
      s.replace(found, 1, "_");
    }
  }
  return s;
}

// --------------------------------------------------------------------
// ---- this resolves only a simple case
void amrex::SimpleRemoveOverlap(BoxArray &ba) {
  for(int dim(0); dim < BL_SPACEDIM; ++dim) {
    for(int i(0); i < ba.size(); ++i) {
      Box b(ba[i]);
      for(int n(0); n < ba.size(); ++n) {
        if(i != n) {
          const Box &bn = ba[n];
	  if(b.bigEnd(dim) == bn.smallEnd(dim)) {
	    b.growHi(dim, - 1);
	    ba.set(i, b);
	  }
        }
      }
    }
  }
}


// --------------------------------------------------------------------
void amrex::avgDown_doit(const FArrayBox &fine_fab, FArrayBox &crse_fab,
                         const Box &ovlp, int scomp, int dcomp, int ncomp,
		         Vector<int> &ratio)
{
    const int  *ovlo   = ovlp.loVect();
    const int  *ovhi   = ovlp.hiVect();
    const int  *flo    = fine_fab.loVect();
    const int  *fhi    = fine_fab.hiVect();
    const Real *f_dat  = fine_fab.dataPtr(scomp);
    const int  *clo    = crse_fab.loVect();
    const int  *chi    = crse_fab.hiVect();
    Real       *c_dat  = crse_fab.dataPtr(dcomp);
    const int  *rr     = ratio.dataPtr();

#if (BL_SPACEDIM == 2)
    FORT_MAXVAL_AVGDOWN(c_dat, AMREX_ARLIM(clo), AMREX_ARLIM(chi), &ncomp,
                    f_dat, AMREX_ARLIM(flo), AMREX_ARLIM(fhi),
                    ovlo, ovhi, rr);
#else
    FORT_CV_AVGDOWN(c_dat, AMREX_ARLIM(clo), AMREX_ARLIM(chi), &ncomp,
                    f_dat, AMREX_ARLIM(flo), AMREX_ARLIM(fhi),
                    ovlo, ovhi, rr);
#endif
}


// --------------------------------------------------------------------
amrex::Box amrex::FixCoarseBoxSize(const Box &fineBox, int rr) {
  // -- we need this because these datasets are not always properly blocked
  Box coarseBox(amrex::coarsen(fineBox, rr));
  for(int dim(0); dim < BL_SPACEDIM; ++dim) {
    while((coarseBox.smallEnd(dim) * rr) < fineBox.smallEnd(dim)) {
      coarseBox.growLo(dim, -1);
    }
    while((coarseBox.bigEnd(dim) * rr) > (fineBox.bigEnd(dim) - rr)) {
      coarseBox.growHi(dim, -1);
    }
  }
  return coarseBox;
}


// --------------------------------------------------------------------
void amrex::avgDown(MultiFab &S_crse, MultiFab &S_fine, int scomp, int dcomp,
                    int ncomp, Vector<int> &ratio)
{
  BL_PROFILE("avgDown()");
    const BoxArray &fgrids = S_fine.boxArray();
    // ---- coarsen on processors owning the fine data.
    BoxArray crse_S_fine_BA(fgrids.size());

    for(int i = 0; i < fgrids.size(); ++i) {
      Box fixedCoarseBox(FixCoarseBoxSize(fgrids[i], ratio[0]));
      crse_S_fine_BA.set(i, fixedCoarseBox);
    }

    const DistributionMapping &S_fine_DM = S_fine.DistributionMap();
    MultiFab crse_S_fine(crse_S_fine_BA, S_fine_DM, ncomp, 0);

    for(MFIter mfi(S_fine); mfi.isValid(); ++mfi) {
      const int i(mfi.index());
      avgDown_doit(S_fine[i], crse_S_fine[i], crse_S_fine_BA[i],
	           scomp, 0, ncomp, ratio);
    }
    S_crse.copy(crse_S_fine, 0, dcomp, ncomp);
}


// ----------------------------------------------------------------------
void amrex::PrintTimeRangeList(const std::list<RegionsProfStats::TimeRange> &trList)
{
  std::list<RegionsProfStats::TimeRange>::const_iterator it;
  for(it = trList.begin(); it != trList.end(); ++it) {
    if(it->startTime > it->stopTime) {
      cout << "  [Time Range Uninitialized] " << endl;
    } else {
      cout << "  " << *it;
    }
  }
  cout << endl;
}


// ----------------------------------------------------------------------
void amrex::RedistFiles() {
    amrex::Abort("**** FixMe:  RedistFiles()");
    bool bIOP(ParallelDescriptor::IOProcessor());
    int myProc(ParallelDescriptor::MyProc());
    int nProcs(ParallelDescriptor::NProcs());
    const Vector<string> &headerFileNames = CommProfStats::GetHeaderFileNames();
    if(nProcs != headerFileNames.size()) {
      if(bIOP) {
        cout << "**** Error:  run with nprocs = " << headerFileNames.size() << endl;
      }
      amrex::Abort();
    }
    std::string cpdir("bl_comm_prof");
    std::string cpredistdir("bl_comm_prof_redist");
    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(cpredistdir, 0755)) {
        amrex::CreateDirectoryFailed(cpredistdir);
      }
    }
    ParallelDescriptor::Barrier("BLProfParser::runRedist::waitfordir");

    std::string hFile(headerFileNames[myProc]);
    std::string dFile(hFile);
    dFile.replace(hFile.find("_H_"), 3, "_D_");

    std::string readHFile(cpdir + '/' + hFile);
    std::string readDFile(cpdir + '/' + dFile);
    std::string writeHFile(cpredistdir + '/' + hFile);
    std::string writeDFile(cpredistdir + '/' + dFile);

    std::ifstream iss;
    std::ofstream oss;
    if(bIOP) {  // main header file
      std::string readCPHFile(cpdir + '/' + "bl_comm_prof_H");
      std::string writeCPHFile(cpredistdir + '/' + "bl_comm_prof_H");
      iss.open(readCPHFile.c_str(), std::ios::in);
      iss.seekg(0, std::ios::end);
      long fileLength(iss.tellg());
      char* charBuf = (char*) malloc(fileLength);
      iss.seekg(0, std::ios::beg);
      iss.read(charBuf, fileLength);
      iss.close();

      oss.open(writeCPHFile.c_str(), std::ios::out);
      oss.write(charBuf, fileLength);
      oss.close();
      free(charBuf);
    }

    {  // header files
      iss.open(readHFile.c_str(), std::ios::in);
      iss.seekg(0, std::ios::end);
      long fileLength(iss.tellg());
      char* charBuf = (char*) malloc(fileLength);
      iss.seekg(0, std::ios::beg);
      iss.read(charBuf, fileLength);
      iss.close();

      oss.open(writeHFile.c_str(), std::ios::out);
      oss.write(charBuf, fileLength);
      oss.close();
      free(charBuf);
    }

    {  // data files
      iss.open(readDFile.c_str(), std::ios::in);
      iss.seekg(0, std::ios::end);
      long fileLength(iss.tellg());
      char* charBuf = (char*) malloc(fileLength);
      iss.seekg(0, std::ios::beg);
      iss.read(charBuf, fileLength);
      iss.close();

      oss.open(writeDFile.c_str(), std::ios::out);
      oss.write(charBuf, fileLength);
      oss.close();
      free(charBuf);
    }
}

// ----------------------------------------------------------------------
int amrex::NHops(const Box &tbox, const IntVect &ivfrom, const IntVect &ivto) {
  int nhops(0);
  for(int d(0); d < BL_SPACEDIM; ++d) {
    int bl(tbox.length(d));
    int ivl(std::min(ivfrom[d], ivto[d]));
    int ivh(std::max(ivfrom[d], ivto[d]));
    int dist(std::min(ivh - ivl, ivl + bl - ivh));
    nhops += dist;
  }
  return nhops;
}


// ----------------------------------------------------------------------
void amrex::Write2DFab(const string &filenameprefix, const int xdim, const int ydim,
                       const double *data)
{
  string filename(filenameprefix + ".2d.fab");
  std::ofstream bfabout(filename.c_str());
  bfabout << "FAB ((8, (64 11 52 0 1 12 0 1023)),(8, (8 7 6 5 4 3 2 1)))((0,0) ("
          << xdim - 1 << "," << ydim - 1 << ") (0,0)) 1" << '\n';
  bfabout.write((char *) data, xdim * ydim * sizeof(double));
  bfabout.close();
}


// ----------------------------------------------------------------------
void amrex::Write2DText(const string &filenameprefix, const int xdim, const int ydim,
                        const double *data)
{
  string filename(filenameprefix + ".2d.txt");
  std::ofstream bfabout(filename.c_str());
  for(int j(0); j < ydim; ++j) {
    for(int i(0); i < xdim; ++i) {
      int index(((ydim - 1 - j) * xdim) + i);
      if(index < 0 || index >= xdim * ydim) {
        cout << "********* bad index = " << index << endl;
        exit(-4);
      }
      bfabout << data[index] << '\t';
    }
    bfabout << endl;
  }
  bfabout.close();
}


// ----------------------------------------------------------------------
void amrex::Write3DFab(const string &filenameprefix, const int xdim, const int ydim,
                const int zdim, const double *data)
{
  string filename(filenameprefix + ".3d.fab");
  std::ofstream bfabout(filename.c_str());
  bfabout << "FAB ((8, (64 11 52 0 1 12 0 1023)),(8, (8 7 6 5 4 3 2 1)))((0,0,0) ("
          << xdim - 1 << "," << ydim - 1 << "," << zdim - 1 << ") (0,0,0)) 1" << '\n';
  bfabout.write((char *) data, xdim * ydim * zdim * sizeof(double));
  bfabout.close();
}


// ----------------------------------------------------------------------
void amrex::WriteFab(const string &filenameprefix, const int xdim, const int ydim,
                     const double *data)
{
  if(xdim <= maxFabXDim) {

    Write2DFab(filenameprefix, xdim, ydim, data);

  } else {     // wrap the 2d fab into a 3d one

    int index(-1);
    int xdim3d(maxFabXDim), ydim3d(ydim), zdim3d(xdim / maxFabXDim);
    if(xdim % maxFabXDim) {
      ++zdim3d;
    }
    Vector<double> data3d(xdim3d * ydim3d * zdim3d, (Real) BLProfiler::InvalidCFT);

    for(int iy(0); iy < ydim; ++iy) {
      for(int ix(0); ix < xdim; ++ix) {
        index = (ix%xdim3d) + (iy*xdim3d) + ((ix/xdim3d)*xdim3d*ydim3d);
        data3d[index] = data[ix + iy*xdim];
      }
    }
    Write3DFab(filenameprefix, xdim3d, ydim3d, zdim3d, &data3d[0]);
  }
}
// ----------------------------------------------------------------------
long amrex::FileSize(const std::string &filename) {
  struct stat statBuff;
  int rCode(::stat(filename.c_str(), &statBuff));
  return(rCode == 0 ? statBuff.st_size : -1);
}


// ----------------------------------------------------------------------
void amrex::MakeFuncPctTimesMF(const Vector<Vector<BLProfStats::FuncStat> > &funcStats,
                               const Vector<std::string> &blpFNames,
		               const std::map<std::string, BLProfiler::ProfStats> &mProfStats,
			       Real runTime, int dataNProcs)
{
#ifdef BL_TRACE_PROFILING
#if (BL_SPACEDIM == 2)
  Vector<std::pair<Real, int> > funcPctTimes(funcStats.size());
  for(int fnum(0); fnum < funcStats.size(); ++fnum) {
    const std::string &fName(blpFNames[fnum]);
    std::map<std::string, BLProfiler::ProfStats>::const_iterator mpsit = mProfStats.find(fName);
    if(mpsit != mProfStats.end()) {
      const std::string &fNamePS(mpsit->first);
      const BLProfiler::ProfStats &ps = mpsit->second;
      if(fNamePS == fName && runTime > 0.0) {
        Real percent((ps.totalTime / dataNProcs) / runTime);
        funcPctTimes[fnum].first  = percent;
        funcPctTimes[fnum].second = fnum;
      } else {
        cout << "**** Error MakeFuncPctTimesMF:  runTime <= 0.0 || fNamePS != fname:  " << runTime << "  "
             << fNamePS << "  " << fName << endl;
      }
    }
  }

  // make a fab with function time data per process
  // needs to be parallelized and plotfiled
  if(funcPctTimes.size() > 0) {
    std::sort(funcPctTimes.begin(), funcPctTimes.end(), BLProfiler::fTTComp());
    Box fptBox(IntVect(0,0), IntVect(dataNProcs - 1, funcPctTimes.size() - 1));
    FArrayBox fptFab(fptBox, 1);
    fptFab.setVal(0.0);
    Real *dptr = fptFab.dataPtr(0);
    for(int fnum(0); fnum < funcStats.size(); ++fnum) {
      int iy(funcPctTimes[fnum].second);
      int xsize(funcStats[iy].size());
      for(int p(0); p < xsize; ++p) {
        int index(((funcStats.size() - 1 - fnum) * xsize) + p);
        const BLProfStats::FuncStat &fs = funcStats[iy][p];
        dptr[index] = fs.totalTime;
      }
    }
    std::ofstream tfout("funcPctTimes.fab");
    fptFab.writeOn(tfout);
    tfout.close();
  }
#endif
#endif
}


// ----------------------------------------------------------------------
void amrex::CollectMProfStats(std::map<std::string, BLProfiler::ProfStats> &mProfStats,
                              const Vector<Vector<BLProfStats::FuncStat> > &funcStats,
                              const Vector<std::string> &fNames,
		              Real runTime, int whichProc)
{
  for(int fnum(0); fnum < funcStats.size(); ++fnum) {
    BLProfiler::ProfStats ps;
    ps.minTime    = runTime;
    ps.maxTime    = 0.0;
    for(int p(0); p < funcStats[fnum].size(); ++p) {
      const BLProfStats::FuncStat &fs = funcStats[fnum][p];
      if(p == whichProc) {  // ---- store ncalls for whichProc only
        ps.nCalls = fs.nCalls;
      }
      ps.totalTime += fs.totalTime;
      ps.minTime    = std::min(ps.minTime, fs.totalTime);
      ps.maxTime    = std::max(ps.maxTime, fs.totalTime);
    }
    Real numProcs(funcStats[fnum].size());
    if(numProcs > 0.0) {
      ps.avgTime = ps.totalTime / numProcs;
    } else {
      ps.avgTime = 0.0;
    }
    Real variance(0.0);
    for(int p(0); p < funcStats[fnum].size(); ++p) {
      const BLProfStats::FuncStat &fs = funcStats[fnum][p];
      variance += (fs.totalTime - ps.avgTime) * (fs.totalTime - ps.avgTime);
    }
    if(funcStats[fnum].size() > 0) {
      ps.variance = variance / static_cast<Real> (funcStats[fnum].size());
    } else {
      ps.variance = 0.0;
    }

    const std::string &fName(fNames[fnum]);
    mProfStats.insert(std::pair<std::string, BLProfiler::ProfStats>(fName, ps));
  }
}


// ----------------------------------------------------------------------
void amrex::GraphTopPct(const std::map<std::string, BLProfiler::ProfStats> &mProfStats,
                        const Vector<Vector<BLProfStats::FuncStat> > &funcStats,
                        const Vector<std::string> &fNames,
		        Real runTime, int dataNProcs, Real gPercent)
{
  for(int fnum(0); fnum < funcStats.size(); ++fnum) {
    const std::string &fName(fNames[fnum]);
    auto mpsit = mProfStats.find(fName);
    if(mpsit != mProfStats.end()) {
      const std::string &fNamePS(mpsit->first);
      const BLProfiler::ProfStats &ps = mpsit->second;
      if(fNamePS == fName && runTime > 0.0) {
        Real percent((ps.totalTime / dataNProcs) / runTime);
        if(percent >= gPercent) {
          const std::string &xgrName(SanitizeName(fName) + ".xgr");
          cout << "xxxxxxxxxxxxxxxxx graphing:  " << xgrName << "  " << percent << endl;
          std::ofstream xgout(xgrName.c_str());
          for(int p(0); p < funcStats[fnum].size(); ++p) {
            const BLProfStats::FuncStat &fs = funcStats[fnum][p];
            xgout << p << " " << fs.totalTime << '\n';
          }
          xgout.close();
        }
      } else {
        cout << "**** Error GraphTopPct:  runTime <= 0.0 || fNamePS != fname:  " << runTime << "  "
             << fNamePS << "  " << fName << endl;
      }
    }
  }
}
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
#endif
