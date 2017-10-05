// ----------------------------------------------------------------------
//  AMReX_BLProfUtilities.H
// ----------------------------------------------------------------------
#ifndef BL_BLPROFUTILITIES_H
#define BL_BLPROFUTILITIES_H

//#include <cstdlib>
//#include <cstdio>
//#include <iostream>
#include <string>
#include <cstring>
//#include <list>
//#include <vector>
#include <typeinfo>
//#include <limits>
//#include <algorithm>
//#include <sys/stat.h>

//using std::cout;
//using std::endl;
//using std::flush;
//using std::string;

//#include <AMReX_SPACE.H>
//#include <AMReX.H>
#include <AMReX_BLProfStats.H>
#include <AMReX_CommProfStats.H>
#include <AMReX_RegionsProfStats.H>
//#include <AMReX_Utility.H>
//#include <AMReX_ParallelDescriptor.H>
//#include <AMReX_Array.H>
#include <AMReX_Vector.H>
//#include <AMReX_IntVect.H>
//#include <AMReX_Box.H>
//#include <AMReX_FArrayBox.H>
//#include <AMReX_MultiFab.H>
//#include <AMReX_AVGDOWN_F.H>

namespace amrex{

// --------------------------------------------------------------------
std::string SanitizeName(const std::string &sname); 
void SimpleRemoveOverlap(BoxArray &ba);
void avgDown_doit(const FArrayBox &fine_fab, FArrayBox &crse_fab,
                  const Box &ovlp, int scomp, int dcomp, int ncomp,
		  Vector<int> &ratio);
Box FixCoarseBoxSize(const Box &fineBox, int rr);
void avgDown(MultiFab &S_crse, MultiFab &S_fine, int scomp, int dcomp,
             int ncomp, Vector<int> &ratio);
void PrintTimeRangeList(const std::list<RegionsProfStats::TimeRange> &trList);
void RedistFiles();
int NHops(const Box &tbox, const IntVect &ivfrom, const IntVect &ivto); 
void Write2DFab(const string &filenameprefix, const int xdim, const int ydim,
                const double *data);
void Write2DText(const string &filenameprefix, const int xdim, const int ydim,
                 const double *data);
void Write3DFab(const string &filenameprefix, const int xdim, const int ydim,
                const int zdim, const double *data);
void WriteFab(const string &filenameprefix, const int xdim, const int ydim,
              const double *data);
long FileSize(const std::string &filename);
void MakeFuncPctTimesMF(const Vector<Vector<BLProfStats::FuncStat> > &funcStats,
                        const Vector<std::string> &blpFNames,
			const std::map<std::string, BLProfiler::ProfStats> &mProfStats,
			Real runTime, int dataNProcs);
void CollectMProfStats(std::map<std::string, BLProfiler::ProfStats> &mProfStats,
                       const Vector<Vector<BLProfStats::FuncStat> > &funcStats,
                       const Vector<std::string> &fNames,
		       Real runTime, int whichProc);
void GraphTopPct(const std::map<std::string, BLProfiler::ProfStats> &mProfStats,
                 const Vector<Vector<BLProfStats::FuncStat> > &funcStats,
                 const Vector<std::string> &fNames,
		 Real runTime, int dataNProcs, Real gPercent);
}
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
#endif
