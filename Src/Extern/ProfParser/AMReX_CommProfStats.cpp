// ----------------------------------------------------------------------
//  CommProfStats.cpp
// ----------------------------------------------------------------------
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <limits>
#include <iomanip>
#include <set>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::flush;
using std::string;
using std::ifstream;
using std::map;
using std::unordered_map;
using std::unordered_multimap;
using std::vector;
using std::pair;

#include <AMReX_CommProfStats.H>
#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfUtilities.H>

using namespace amrex;

#ifdef _OPENMP
#include <omp.h>
#endif

#define SHOWVAL(val) { cout << #val << " = " << val << endl; }

const int XDIR(0);
const int YDIR(1);
const int ZDIR(2);

//Vector<Vector<Real> > barrierExitTimes;  // [proc, bnum]
//Vector<Vector<Real> > barrierSkewTimes;
bool bExitTimesDone(false);
bool bSkewTimesDone(false);

int nTopPts(0);
int maxTopNodeNum(0);

Box topoBox;
std::map<int, IntVect> pNumTopIVMap;  // [procNumber, topological iv position]
bool bTopoFabInited(false);
Vector<int> procNodeNumber;

bool CommProfStats::bInitDataBlocks(true);

Vector<int> CommProfStats::rankFromProx;  // [prox]
Vector<int> CommProfStats::proxFromRank;  // [rank]
bool CommProfStats::bProxMapOK(false);
bool CommProfStats::persistentStreams(true);

int CommProfStats::cpVersion(-1);
int CommProfStats::csSize(-1);
int CommProfStats::finestLevel(-1);
int CommProfStats::maxLevel(-1);
Vector<IntVect> CommProfStats::calcRefRatios;
Vector<Box> CommProfStats::probDomain;
Vector<std::string> CommProfStats::commHeaderFileNames;
std::map<std::string, int> CommProfStats::commDataFileNames;
Vector<std::ifstream *> CommProfStats::commDataStreams;


extern int amrex::NHops(const Box &tbox, const IntVect &ivfrom, const IntVect &ivto);
extern void amrex::Write2DFab(const string &filenameprefix, const int xdim, const int ydim,
                              const double *data);
extern void amrex::Write2DText(const string &filenameprefix, const int xdim, const int ydim,
                               const double *data);
extern void amrex::Write3DFab(const string &filenameprefix, const int xdim, const int ydim,
                              const int zdim, const double *data);
extern void amrex::WriteFab(const string &filenameprefix, const int xdim, const int ydim,
                            const double *data);


// ----------------------------------------------------------------------
IntVect TopIVFromProcNumber(const int procnum) {
  IntVect iv;
  std::map<int, IntVect>::iterator it = pNumTopIVMap.find(procnum);
  if(it == pNumTopIVMap.end()) {
      std::cerr << "**** Error in TopIVFromProcNumber:  procnum not found:  "
                << procnum << std::endl;
  } else {
    iv = it->second;
    if(iv != pNumTopIVMap[procnum]) {
      std::cerr << "**** Error in TopIVFromProcNumber:  procnum not matched:  "
                << procnum << std::endl;
    }
  }
  return iv;
}


// ----------------------------------------------------------------------
CommProfStats::CommProfStats()
  :  currentDataBlock(0)
{
  TopoMap.resize(1);

  dataBlocks.reserve(2048);
}


// ----------------------------------------------------------------------
CommProfStats::~CommProfStats() {
}


// ----------------------------------------------------------------------
void CommProfStats::InitDataFileNames(const Vector<std::string> &hfn) {
  for(int i(0); i < hfn.size(); ++i) {
    std::string dFileName(hfn[i]);
    dFileName.replace(dFileName.find("_H_"), 3, "_D_");
    commDataFileNames.insert(std::pair<std::string, int>(dFileName, i));
  }
}


// ----------------------------------------------------------------------
void CommProfStats::InitCommDataBlock(const int proc, const long ncommstats,
                                      const std::string &filename,
				      const long seekpos,
				      const std::string &procname,
				      const int nodenumber)
{
  int streamindex;
  std::map<std::string, int>::iterator it =  commDataFileNames.find(filename);
  if(it == commDataFileNames.end()) {
    streamindex = commDataFileNames.size();
    commDataFileNames.insert(std::pair<std::string, int>(filename, streamindex));
  } else {
    streamindex = it->second;
  }
  if(bInitDataBlocks) {
    currentProc = proc;
    dataBlocks.push_back(DataBlock(proc, ncommstats, filename, seekpos,
                                   procname, nodenumber, streamindex));
    currentDataBlock = dataBlocks.size() - 1;
  }
}


// ----------------------------------------------------------------------
int CommProfStats::AfterBarrier(const int /*proc*/, const double /*t*/) {
/*
  if(bExitTimesDone) {
    if(barrierExitTimes[proc][0] > t) {
      return 0;
    }
    for(int i(0); i < barrierExitTimes[proc].size() - 1; ++i) {
      if(barrierExitTimes[proc][i + 1] > t) {
        return i;
      }
    }
    return (barrierExitTimes[proc].size() - 1);
  } else {
    cerr << "**** Error in AfterBarrier:  bExitTimesDone == false" << endl;
    return -1;
  }
*/
    return -1;
}


// ----------------------------------------------------------------------
void CommProfStats::AddBarrier(long bnum, const std::string &bname, long index)
{
  int bNameNumber(-1);
  std::map<int, std::string>::iterator it = barrierNumbersToNames.find(bnum);
  if(it == barrierNumbersToNames.end()) {
    bNameNumber = barrierNames.size();
    barrierNames.push_back(bname);
    barrierNumbersToNames.insert(std::pair<int, std::string>(bNameNumber, bname));
  } else {
    bNameNumber = bnum;
  }
  if(dataBlocks[currentDataBlock].barriers.size() <= bnum) {
    dataBlocks[currentDataBlock].barriers.resize(bnum + 1);
    dataBlocks[currentDataBlock].barriers[bnum] = BarrierEntry(bnum, bname,index, bNameNumber);
  }
}


// ----------------------------------------------------------------------
void CommProfStats::AddReduction(const long rnum, const long index) {
  dataBlocks[currentDataBlock].reductions.push_back(ReductionEntry(rnum, index));
}


// ----------------------------------------------------------------------
void CommProfStats::AddNameTagName(const string &name) {
  if(std::find(nameTagNames.begin(), nameTagNames.end(), name) == nameTagNames.end()) {
    nameTagNames.push_back(name);
  }
}


// ----------------------------------------------------------------------
void CommProfStats::AddNameTag(const long ntindex, const long seekindex) {
  dataBlocks[currentDataBlock].nameTags.push_back(NameTagEntry(ntindex, seekindex));
}


// ----------------------------------------------------------------------
void CommProfStats::AddTagRange(const long tmin, const long tmax) {
  tagMin = tmin;
  tagMax = tmax;
}


// ----------------------------------------------------------------------
void CommProfStats::AddTimeMinMax(const double tmin, const double tmax) {
  dataBlocks[currentDataBlock].timeMin = tmin;
  dataBlocks[currentDataBlock].timeMax = tmax;
}


// ----------------------------------------------------------------------
void CommProfStats::AddTimerTime(const double tt) {
  dataBlocks[currentDataBlock].timerTime = tt;
  //++nTimerTimes;
}


// ----------------------------------------------------------------------
void CommProfStats::AddGridLevel(const int /*level*/, const int /*ngrids*/) {
}


// ----------------------------------------------------------------------
void CommProfStats::AddGrid3D(int /*level*/, int /*xlo*/, int /*ylo*/, int /*zlo*/,
			      int /*xhi*/, int /*yhi*/, int /*zhi*/,
			      int /*xc*/,  int /*yc*/,  int /*zc*/,
			      int xn,  int yn,  int zn, int proc)
{
  long nPoints(xn * yn * zn);
  std::map<int, long>::iterator it = glMap.find(proc);
  if(it == glMap.end()) {
    glMap.insert(std::make_pair(proc, nPoints));
  } else {
    glMap[proc] += nPoints;
  }
  int nP(nPoints);
  std::map<int, int>::iterator its = glSizeMap.find(nP);
  if(its == glSizeMap.end()) {
    glSizeMap.insert(std::make_pair(nP, 1));
  } else {
    glSizeMap[nP] += 1;
  }
}


// ----------------------------------------------------------------------
void CommProfStats::AddRefRatio(const int level, const IntVect &rr) {
  if((level + 1) > calcRefRatios.size()) {
    calcRefRatios.resize(level + 1);
  }
  calcRefRatios[level] = rr;
}


// ----------------------------------------------------------------------
void CommProfStats::AddProbDomain(const int level, const Box &pd) {
  if((level + 1) > probDomain.size()) {
    probDomain.resize(level + 1);
  }
  probDomain[level] = pd;
}


// ----------------------------------------------------------------------
void CommProfStats::AddTopoCoord(const int nid, const int node,
                                 const int tx, const int ty, const int tz,
                                 const bool /*servicenode*/)
{
#if (BL_SPACEDIM == 2)
  amrex::ignore_unused(nid, node, tx, ty, tz);
  cout << "**** Error:  CommProfStats::AddTopoCoord not supported for 2D" << endl;
#else
  cout << "TopoMap.size() = " << TopoMap.size() << " nid node = " << nid << "  " << node << endl;
  TopoMap[node].insert(std::pair<int, IntVect>(nid, IntVect(AMREX_D_DECL(tx, ty, tz))));
  ++nTopPts;
  maxTopNodeNum = std::max(maxTopNodeNum, node);
#endif
}


// ----------------------------------------------------------------------
void CommProfStats::AddCommHeaderFileName(const string &hfn) {
  if(std::find(commHeaderFileNames.begin(),
     commHeaderFileNames.end(), hfn) == commHeaderFileNames.end())
  {
    commHeaderFileNames.push_back(hfn);
  }
}


// ----------------------------------------------------------------------
void CommProfStats::InitProxMap() {
  string filename("RankProxOrder.txt");
  ifstream rpo(filename.c_str());
  if( ! rpo.good()) {
    bProxMapOK = false;
    if(ParallelDescriptor::IOProcessor()) {
      cout << "**** Error in CommProfStats::InitProxMap:  cannot open file." << endl;
    }
  } else {
    int nprocs, r, p;
    rpo >> nprocs;
    cout << "CommProfStats::InitProxMap:  nprocs = " << nprocs << endl;
    rankFromProx.resize(nprocs);
    proxFromRank.resize(nprocs);
    for(int i(0); i < nprocs; ++i) {
      rpo >> r >> p;
      rankFromProx[p] = r;
      proxFromRank[r] = p;
    }
    rpo.close();
    bProxMapOK = true;
  }
}


// ----------------------------------------------------------------------
void CommProfStats::WriteTopoFab() {
#if (BL_SPACEDIM == 2)
  cout << "**** Error:  CommProfStats::WriteTopoFab not supported for 2D" << endl;
#else
  IntVect ivmin(AMREX_D_DECL(100000, 100000, 100000));
  IntVect ivmax(AMREX_D_DECL(-100000, -100000, -100000));
  std::map<int, IntVect>::iterator it;
  for(int i(0); i < TopoMap.size(); ++i) {
    for(it = TopoMap[i].begin(); it != TopoMap[i].end(); ++it) {
      ivmin.min(it->second);
      ivmax.max(it->second);
    }
  }
  Box tDomain(ivmin, ivmax);
  cout << "tDomain = " << tDomain << "  npts = " << tDomain.numPts()
       << "  nTopPts = " << nTopPts << "  maxTopNodeNum = " << maxTopNodeNum << endl;
  FArrayBox tFab(tDomain, maxTopNodeNum + 1);
  tFab.setVal<RunOn::Host>(-1);
  for(int i(0); i < TopoMap.size(); ++i) {
    for(it = TopoMap[i].begin(); it != TopoMap[i].end(); ++it) {
      tFab(it->second, i) = it->first;
    }
  }
  std::ofstream tfout("topolcoords.3d.fab");
  tFab.writeOn(tfout);
  tfout.close();

  Box tBox(tFab.box());
  cout << "tBox = " << tBox << "  ncomp = " << tFab.nComp() << endl;
  topoBox = tBox;

  for(int nc(0); nc < tFab.nComp(); ++nc) {
    for(IntVect iv(tBox.smallEnd()); iv <= tBox.bigEnd(); tBox.next(iv)) {
      int pnum(tFab(iv, nc));
      if(pnum >= 0) {
        //std::cout << ">>>> iv pnum = " << iv << "  " << pnum << std::endl;
        pNumTopIVMap.insert(std::pair<int, IntVect>(pnum, iv));
        //topIVpNumMM.insert(std::pair<IntVect, int>(iv, pnum));
      }
    }
  }
  bTopoFabInited = true;
#endif
}


// ----------------------------------------------------------------------
void CommProfStats::OpenAllStreams(const std::string &dirname) {
  BL_PROFILE_VAR("CommProfStats::OpenAllStreams", cpsopenallstreams);
  commDataStreams.resize(commDataFileNames.size());
  int dsIndex(0);
  for(std::map<std::string, int>::iterator it = commDataFileNames.begin();
      it != commDataFileNames.end(); ++it)
  {
    std::string fullFileName(dirname + '/' + it->first);
    commDataStreams[dsIndex] = new std::ifstream(fullFileName.c_str());

    if (commDataStreams[dsIndex]->fail())
    {
      cout << "****commDataStreams failed. Continuing without persistent streams." << std::endl;
      persistentStreams = false;
      CloseAllStreams();
      break;
    }

    ++dsIndex;
  }
  BL_PROFILE_VAR_STOP(cpsopenallstreams);
}


// ----------------------------------------------------------------------
void CommProfStats::CloseAllStreams() {
  BL_PROFILE_VAR("CommProfStats::CloseAllStreams", cpsclosellstreams);
  for(int i(0); i < commDataStreams.size(); ++i) {
    if (commDataStreams[i] != nullptr)
    {
      if (commDataStreams[i]->is_open())
      {
        commDataStreams[i]->close();
      }
      delete commDataStreams[i];
      commDataStreams[i] = nullptr;
    }
  }
  BL_PROFILE_VAR_STOP(cpsclosellstreams);
}


// ----------------------------------------------------------------------
void CommProfStats::ReadCommStats(DataBlock &dBlock) {
  if(dBlock.vCommStats.size() != dBlock.size) {
    dBlock.vCommStats.resize(dBlock.size);
  }
  std::string fullFileName(dirName + '/' + dBlock.fileName);
  BL_PROFILE_VAR("OpenStream", openstream);
  std::ifstream instr(fullFileName.c_str());
  BL_PROFILE_VAR_STOP(openstream);
  long dataSize(dBlock.size * csSize);
  instr.seekg(dBlock.seekpos);
  instr.read(reinterpret_cast<char *>(dBlock.vCommStats.dataPtr()), dataSize);
  instr.close();
}


// ----------------------------------------------------------------------
void CommProfStats::ReadCommStatsNoOpen(DataBlock &dBlock) {
  if(dBlock.vCommStats.size() != dBlock.size) {
    dBlock.vCommStats.resize(dBlock.size);
  }
  std::ifstream *instr = commDataStreams[dBlock.streamIndex];
  long dataSize(dBlock.size * csSize);
  instr->seekg(dBlock.seekpos);
  instr->read(reinterpret_cast<char *>(dBlock.vCommStats.dataPtr()), dataSize);
}


// ----------------------------------------------------------------------
bool CommProfStats::ReadCommStats(DataBlock &dBlock, const int nmessages) {
  int leftToRead(dBlock.size - dBlock.readoffset);
  int readSize(std::min(leftToRead, nmessages));
  int readPos(dBlock.seekpos + dBlock.readoffset * csSize);
  if(dBlock.vCommStats.size() != readSize) {
    dBlock.vCommStats.resize(readSize);
  }
  std::string fullFileName(dirName + '/' + dBlock.fileName);

  std::ifstream instr(fullFileName.c_str());
  int dataSize(readSize * csSize);
  instr.seekg(readPos);
  instr.read((char *) dBlock.vCommStats.dataPtr(), dataSize);
  instr.close();

  dBlock.readoffset += readSize;
  return(dBlock.readoffset < dBlock.size);
}


// ----------------------------------------------------------------------
void CommProfStats::ClearCommStats(DataBlock &dBlock) {
  dBlock.barriers.clear();
  dBlock.reductions.clear();
  dBlock.nameTags.clear();
  dBlock.vCommStats.clear();

  Vector<BarrierEntry>().swap(dBlock.barriers);           // delete memory
  Vector<ReductionEntry>().swap(dBlock.reductions);
  Vector<NameTagEntry>().swap(dBlock.nameTags);
  Vector<BLProfiler::CommStats>().swap(dBlock.vCommStats);
}


// ----------------------------------------------------------------------
void CommProfStats::CheckCommData(Vector<long> &nBMin, Vector<long> &nBMax,
                                  Vector<long> &nRMin, Vector<long> &nRMax)
{
  int myProc(ParallelDescriptor::MyProc());
  cout << myProc << ":  " << "------------------------------------ checking comm data." << endl;
  SHOWVAL(dataNProcs);
  SHOWVAL(csSize);
  SHOWVAL(dataBlocks.size());
  cout << myProc << ":  " << "----" << endl;

  Vector<long> nBarriers(dataNProcs, 0);
  Vector<long> nReductions(dataNProcs, 0);
  if(nBMin.size() != dataNProcs) {
    nBMin.resize(dataNProcs, std::numeric_limits<long>::max());
  }
  if(nBMax.size() != dataNProcs) {
    nBMax.resize(dataNProcs, std::numeric_limits<long>::min());
  }
  if(nRMin.size() != dataNProcs) {
    nRMin.resize(dataNProcs, std::numeric_limits<long>::max());
  }
  if(nRMax.size() != dataNProcs) {
    nRMax.resize(dataNProcs, std::numeric_limits<long>::min());
  }

  int idb;

  for(idb = 0; idb < dataBlocks.size(); ++idb) {  // ---- go through dataBlocks
    DataBlock &dBlock = dataBlocks[idb];
    if(verbose) {
      cout << myProc << ":  " << "CommProfProc  " << dBlock.proc << "  nCommStats  "
           << dBlock.size << "  datafile  " << dBlock.fileName << "  seekpos  "
	   << dBlock.seekpos << endl;
      cout << myProc << ":  " << "barriers.size() =  " << dBlock.barriers.size()
           << "  " << dBlock.proc << endl;
    }

    ReadCommStats(dBlock);

    // --------------------- check barrier integrity
    for(int idbb(0); idbb < dBlock.barriers.size(); ++idbb) {
      BarrierEntry &be = dBlock.barriers[idbb];
      int bNumber(be.number);
      std::string bName(be.name);
      int bIndex(be.seekIndex);
      int nCS(dBlock.vCommStats.size());

      if(bIndex > nCS-1 || bIndex+1 > nCS-1) {
        cerr << "**** Error:  bad bIndex:  " << bIndex << "  " << nCS << endl;
	continue;
      }
      BLProfiler::CommStats &cs = dBlock.vCommStats[bIndex];
      BLProfiler::CommStats &csNext = dBlock.vCommStats[bIndex + 1];

      if(cs.cfType == BLProfiler::Barrier && csNext.cfType == BLProfiler::Barrier &&
         cs.tag == csNext.tag  &&     // these are the barrier numbers
	 cs.tag == bNumber )
      {
        if(cs.commpid != BLProfiler::BeforeCall() || csNext.commpid != BLProfiler::AfterCall()) {
          cerr << "**** Error:  bad Barrier before, after." << endl;
          cerr << BLProfiler::CommStats::CFTToString(cs.cfType) << "   "
               << BLProfiler::CommStats::CFTToString(csNext.cfType) << "   "
               << (cs.commpid) << "   " << (csNext.commpid) << "   "
               << (cs.tag) << "   " << (csNext.tag) << "   "
               << (cs.timeStamp) << "   " << (csNext.timeStamp) << "   "
	       << endl;
	}
      } else {
        cerr << "**** Error:  bad Barriers." << endl;
          cerr << BLProfiler::CommStats::CFTToString(cs.cfType) << "   "
               << BLProfiler::CommStats::CFTToString(csNext.cfType) << "   "
               << (cs.commpid) << "   " << (csNext.commpid) << "   "
               << (cs.tag) << "   " << (csNext.tag) << "   "
               << (cs.timeStamp) << "   " << (csNext.timeStamp) << "   "
	       << bName << "  " << bNumber
	       << endl;
      }
    }

    // --------------------- check reduction integrity
    for(int idbb(0); idbb < dBlock.reductions.size(); ++idbb) {
      ReductionEntry &re = dBlock.reductions[idbb];
      //int rNumber(re.number);
      int bIndex(re.seekIndex);
      int nCS(dBlock.vCommStats.size());

      if(bIndex > nCS-1 || bIndex+1 > nCS-1) {
        cerr << "**** Error:  bad reduction bIndex:  " << bIndex << "  " << nCS << endl;
	continue;
      }

      BLProfiler::CommStats &cs = dBlock.vCommStats[bIndex];
      BLProfiler::CommStats &csNext = dBlock.vCommStats[bIndex + 1];

      if(cs.cfType == csNext.cfType &&  // need to check each type
         cs.tag == csNext.tag)  // these are the reductions numbers
      {
        if(cs.commpid != BLProfiler::BeforeCall() || csNext.commpid != BLProfiler::AfterCall()) {
          cerr << "**** Error:  bad Reduction before, after." << endl;
          cerr << BLProfiler::CommStats::CFTToString(cs.cfType) << "   "
               << BLProfiler::CommStats::CFTToString(csNext.cfType) << "   "
               << (cs.size) << "   " << (csNext.size) << "   "
               << (cs.commpid) << "   " << (csNext.commpid) << "   "
               << (cs.tag) << "   " << (csNext.tag) << "   "
               << (cs.timeStamp) << "   " << (csNext.timeStamp) << "   "
	       << endl;
	}
      } else {
        cerr << "**** Error:  bad Reductions." << endl;
          cerr << BLProfiler::CommStats::CFTToString(cs.cfType) << "   "
               << BLProfiler::CommStats::CFTToString(csNext.cfType) << "   "
               << (cs.size) << "   " << (csNext.size) << "   "
               << (cs.commpid) << "   " << (csNext.commpid) << "   "
               << (cs.tag) << "   " << (csNext.tag) << "   "
               << (cs.timeStamp) << "   " << (csNext.timeStamp) << "   "
	       << endl;
      }
    }

    for(int idbb(0); idbb < dBlock.barriers.size(); ++idbb) {
      BarrierEntry &be = dBlock.barriers[idbb];
      nBarriers[dBlock.proc] = std::max(nBarriers[dBlock.proc], be.number);
    }
    for(int idbb(0); idbb < dBlock.reductions.size(); ++idbb) {
      ReductionEntry &re = dBlock.reductions[idbb];
      nReductions[dBlock.proc] = std::max(nReductions[dBlock.proc], re.number);
    }

    ClearCommStats(dBlock);
  }

  //nBMin = std::min(nBMin, (*std::min_element(nBarriers.begin(), nBarriers.end())));
  //nBMax = std::max(nBMax, (*std::max_element(nBarriers.begin(), nBarriers.end())));
  //cout << "nBarriers.minmax = " << nBMin << "  " << nBMax << endl;
  //if(nBMin != nBMax) {
    //cerr << "**** Error:  different number of barriers per processor." << endl;
  //}

  //nRMin = std::min(nRMin, (*std::min_element(nReductions.begin(), nReductions.end())));
  //nRMax = std::max(nRMax, (*std::max_element(nReductions.begin(), nReductions.end())));
  //cout << "nReductions.minmax = " << nRMin << "  " << nRMax << endl;
  //if(nRMin != nRMax) {
    //cerr << "**** Error:  different number of reductions per processor." << endl;
  //}

}


// ----------------------------------------------------------------------
void CommProfStats::FillSendFAB(long &totalSends, long &totalSentData,
				Vector<long> &totalSendsPerProc,
				Vector<long> &totalSentDataPerProc,
				FArrayBox &sendFAB, bool proxmap)
{
  BL_PROFILE("CommProfStats::FillSendFAB");

#if (BL_SPACEDIM == 2)

  if(proxmap) {

  for(int idb(0); idb < dataBlocks.size(); ++idb) {    // ---- go through dataBlocks
    Real *tsp  = sendFAB.dataPtr(0);
    Real *tsdp = sendFAB.dataPtr(1);
    Box dataBox(sendFAB.box());
    int smallX(dataBox.smallEnd(XDIR)), bigX(dataBox.bigEnd(XDIR));
    int smallY(dataBox.smallEnd(YDIR)), bigY(dataBox.bigEnd(YDIR));
    int xlen(dataBox.length(XDIR)), rankFrom, rankTo, proxFrom, proxTo, proc;

    DataBlock &dBlock = dataBlocks[idb];
    rankFrom = dBlock.proc;
    proxFrom = proxFromRank[rankFrom];
    proc = dBlock.proc;

    if(proxFrom >= smallX && proxFrom <= bigX) {    // ---- within from proc range
      int index, offsetX(proxFrom - smallX);

      BL_PROFILE_VAR("FillSendFABIO", fillsendfabio);

      if (persistentStreams){
        ReadCommStatsNoOpen(dBlock);
      } else {
        ReadCommStats(dBlock);
      }

      BL_PROFILE_VAR_STOP(fillsendfabio);

      for(int i(0); i < dBlock.vCommStats.size(); ++i) {    // ---- find sends and sum
        BLProfiler::CommStats &cs = dBlock.vCommStats[i];
        if(IsSend(cs.cfType)) {
	  if(cs.size != BLProfiler::AfterCall()) {
	    rankTo = cs.commpid;
	    proxTo = proxFromRank[rankTo];
	    if(proxTo >= smallY && proxTo <= bigY) {    // ---- within to proc range
	      if(InTimeRange(proc, cs.timeStamp)) {
	        ++totalSends;
                totalSentData += cs.size;
	        ++totalSendsPerProc[proxFrom];
	        totalSentDataPerProc[proxFrom] += cs.size;
	        index = (offsetX) + (xlen * (proxTo - smallY));
	        tsp[index]  += 1.0;
	        tsdp[index] += cs.size;
	      }
	    }
	  }
        }
      }
      ClearCommStats(dBlock);
    }
  }


  } else {


  for(int idb(0); idb < dataBlocks.size(); ++idb) {    // ---- go through dataBlocks
    //cout << idb + 1 << "  " << flush;
    //if((idb + 1) % 20 == 0) {
      //cout << endl;
    //}
    Real *tsp  = sendFAB.dataPtr(0);
    Real *tsdp = sendFAB.dataPtr(1);
    Box dataBox(sendFAB.box());
    int smallX(dataBox.smallEnd(XDIR)), bigX(dataBox.bigEnd(XDIR));
    int smallY(dataBox.smallEnd(YDIR)), bigY(dataBox.bigEnd(YDIR));
    int xlen(dataBox.length(XDIR)), proc;

    DataBlock &dBlock = dataBlocks[idb];
    proc = dBlock.proc;

    if(proc >= smallX && proc <= bigX) {    // ---- within from proc range
      int index, offsetX(proc - smallX);

      BL_PROFILE_VAR("FillSendFABIO", fillsendfabio);

      if (persistentStreams){
        ReadCommStatsNoOpen(dBlock);
      } else {
        ReadCommStats(dBlock);
      }

      BL_PROFILE_VAR_STOP(fillsendfabio);

      for(int i(0); i < dBlock.vCommStats.size(); ++i) {    // ---- find sends and sum
        BLProfiler::CommStats &cs = dBlock.vCommStats[i];
        if(IsSend(cs.cfType)) {
	  if(cs.size != BLProfiler::AfterCall()) {
	    if(cs.commpid >= smallY && cs.commpid <= bigY) {    // ---- within to proc range
	      if(InTimeRange(proc, cs.timeStamp)) {
	        ++totalSends;
                totalSentData += cs.size;
	        ++totalSendsPerProc[proc];
	        totalSentDataPerProc[proc] += cs.size;
	        index = (offsetX) + (xlen * (cs.commpid - smallY));
	        tsp[index]  += 1.0;
	        tsdp[index] += cs.size;
	      }
	    }
	  }
        }
      }
      ClearCommStats(dBlock);
    }
  }

  }

#else
  cerr << "**** Error in ReportSendsFABS:  must compile with DIM = 2." << endl;
#endif
}


// ----------------------------------------------------------------------
void CommProfStats::ReportSyncPointDataSetup(long &nBMax, long &nRMax)
{
  Vector<long> nBarriers(dataNProcs, 0L);
  Vector<long> nReductions(dataNProcs, 0L);
  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
    int proc(dBlock.proc);

    for(int idbb(0); idbb < dBlock.barriers.size(); ++idbb) {
      BarrierEntry &be = dBlock.barriers[idbb];
      nBarriers[proc] = std::max(nBarriers[proc], be.number);
    }
    for(int idbb(0); idbb < dBlock.reductions.size(); ++idbb) {
      ReductionEntry &re = dBlock.reductions[idbb];
      nReductions[proc] = std::max(nReductions[proc], re.number);
    }
  }
  //nBMin = *std::min_element(nBarriers.begin(), nBarriers.end());
  nBMax = *std::max_element(nBarriers.begin(), nBarriers.end());
  //nRMin = *std::min_element(nReductions.begin(), nReductions.end());
  nRMax = *std::max_element(nReductions.begin(), nReductions.end());
}


// ----------------------------------------------------------------------
void CommProfStats::ReportSyncPointData(Vector<Vector<Real> > &barrierExitTimes,
                                        Vector<Vector<Real> > &barrierWaitTimes,
                                        Vector<Vector<Real> > &reductionWaitTimes,
					bool bDoReductions)
{
  procNodeNumber.resize(dataNProcs);

  for(int idb(0); idb < dataBlocks.size(); ++idb) {  // ---- go through dataBlocks
    //cout << idb + 1 << "  " << flush;
    //if((idb + 1) % 20 == 0) {
      //cout << endl;
    //}
    DataBlock &dBlock = dataBlocks[idb];

    if (persistentStreams){
      ReadCommStatsNoOpen(dBlock);
    } else {
      ReadCommStats(dBlock);
    }

    // ------------------------------------------------ collect barrier timings
    for(int i(0); i < dBlock.barriers.size(); ++i) {
      BarrierEntry &be = dBlock.barriers[i];
      if(be.seekIndex < 0) {  // ---- skip these unused entries, results from flushing
        //cout << "***************** be.seekIndex < 0 :: = " << be.seekIndex
             //<< "  i = " << i << "  dBlock.barriers.size() = " << dBlock.barriers.size()
             //<< "  be.name = " << be.name << endl;
        continue;
      }
      BLProfiler::CommStats &cs = dBlock.vCommStats[be.seekIndex];
      BLProfiler::CommStats &csNext = dBlock.vCommStats[be.seekIndex + 1];

      if(InTimeRange(dBlock.proc, cs.timeStamp)) {
        barrierWaitTimes[dBlock.proc][be.number] = csNext.timeStamp - cs.timeStamp;
        barrierExitTimes[dBlock.proc][be.number] = csNext.timeStamp;
      }

      //double zeroBarrierTime = barrierExitTimes[0][be.number];
      //for(int isk(0); isk < barrierSkewTimes.size(); ++isk) {
        //barrierSkewTimes[isk][be.number] = barrierExitTimes[isk][be.number] - zeroBarrierTime;
      //}
    }
    bExitTimesDone = true;
    bSkewTimesDone = true;

    // ------------------------------------------------ collect reduction wait times
    for(int i(0); i < dBlock.reductions.size(); ++i) {
      ReductionEntry &re = dBlock.reductions[i];
      BLProfiler::CommStats &cs = dBlock.vCommStats[re.seekIndex];
      BLProfiler::CommStats &csNext = dBlock.vCommStats[re.seekIndex + 1];

      if(bDoReductions) {
        if(InTimeRange(dBlock.proc, cs.timeStamp)) {
          reductionWaitTimes[dBlock.proc][re.number] = csNext.timeStamp - cs.timeStamp;
        }
      }
    }

    procNodeNumber[dBlock.proc] = dBlock.nodeNumber;

    ClearCommStats(dBlock);
  }

  // ------------------------------------------------ proc node numbers
  std::ofstream pnnout("procNodeNumber.xgr");
  for(int ip(0); ip < dataNProcs; ++ip) {
    pnnout << ip << " " << procNodeNumber[ip] << '\n';
  }
  pnnout.close();
}


// ----------------------------------------------------------------------
void CommProfStats::ReportStats(long &totalSentData, long &totalNCommStats,
                                Vector<long> &totalFuncCalls,
				int bytesPerSlot, Vector<long> &msgSizes,
				int &minMsgSize, int &maxMsgSize,
                                Real &timeMin, Real &timeMax, Real &timerTime,
				Vector<int> &rankNodeNumbers)
{
  amrex::ignore_unused(timerTime);

  amrex::Print(Print::AllProcs) << ParallelDescriptor::MyProc() << "::Processing "
                                << dataBlocks.size() << " data blocks:" << endl;

  int nMsgSizes(msgSizes.size());
  int slot(0), highSlot(nMsgSizes - 1);

  for(int idb(0); idb < dataBlocks.size(); ++idb) {  // ---- go through dataBlocks
    if(verbose) {
      //cout << idb + 1 << "  " << flush;
      //if((idb + 1) % 20 == 0) {
        //cout << endl;
      //}
    }
    DataBlock &dBlock = dataBlocks[idb];
    if (persistentStreams){
      ReadCommStatsNoOpen(dBlock);
    } else {
      ReadCommStats(dBlock);
    }

    rankNodeNumbers[dBlock.proc] = dBlock.nodeNumber;

    totalNCommStats += dBlock.size;
    for(int i(0); i < dBlock.vCommStats.size(); ++i) {  // ------- sum sent data
      BLProfiler::CommStats &cs = dBlock.vCommStats[i];
      if(IsSend(cs.cfType)) {
	if(cs.size != BLProfiler::AfterCall()) {
	  if(InTimeRange(dBlock.proc, cs.timeStamp)) {
            totalSentData += cs.size;
	    slot = std::min(cs.size/bytesPerSlot, highSlot);
	    ++msgSizes[slot];
	    minMsgSize = std::min(cs.size, minMsgSize);
	    maxMsgSize = std::max(cs.size, maxMsgSize);
	  }
	}
      }
    }

    for(int i(0); i < dBlock.vCommStats.size(); ++i) {  // ----- sum function calls
      BLProfiler::CommStats &cs = dBlock.vCommStats[i];
      if((cs.size > -1 && cs.cfType != BLProfiler::Waitsome) ||
         (cs.size == BLProfiler::BeforeCall() && cs.cfType == BLProfiler::Waitsome))
      {
	if(InTimeRange(dBlock.proc, cs.timeStamp)) {
	  if(cs.cfType >= 0 && cs.cfType < totalFuncCalls.size()) {
            ++totalFuncCalls[cs.cfType];
	  } else {
	    std::cout << "--------:: totalFuncCalls.size() cs.cfType = " << totalFuncCalls.size()
                      << "  " << cs.cfType << std::endl;
	  }
	}
      }
    }

    // ------------------------------------------------ find minmax times
    timeMin = std::min(timeMin, dBlock.timeMin);
    timeMax = std::max(timeMax, dBlock.timeMax);

    // this sums timerTimes for all datablocks for this mpi process (not dataNProcs)
    //timerTime += dBlock.timerTime;

    ClearCommStats(dBlock);
  }
  //if(nTimerTimes > 0) {
    //timerTime /= nTimerTimes;
  //}
  amrex::Print(Print::AllProcs) << ParallelDescriptor::MyProc() << "::done." << '\n';
}


// ----------------------------------------------------------------------
void CommProfStats::FindTimeRange(BLProfStats::TimeRange& tr) {
  for(int idb(0); idb < dataBlocks.size(); ++idb) {  // ---- go through dataBlocks
    DataBlock &dBlock = dataBlocks[idb];
    tr.startTime = std::min(tr.startTime, dBlock.timeMin);
    tr.stopTime  = std::max(tr.stopTime,  dBlock.timeMax);
  }
}

// ----------------------------------------------------------------------
void CommProfStats::TimelineFAB(FArrayBox &timelineFAB, const Box &probDomain,
                                const BLProfStats::TimeRange tr,
                                const int rankMin, const int rankMax,
			        const int rankStride,
				const Real ntnMultiplier, const Vector<Real> &ntnNumbers,
				const Real bnMultiplier, const Vector<Real> &bnNumbers)
{
  BL_PROFILE("CommProfStats::TimelineFAB()");

  amrex::ignore_unused(rankMin, rankMax, ntnMultiplier, bnMultiplier);

  Real tlo = tr.startTime;
  Real thi = tr.stopTime;
  Real timeRangeAll(thi - tlo);
  //Real ooTimeRangeAll(1.0 / timeRangeAll);
  Real dt(timeRangeAll / probDomain.length(XDIR));
  Real fabTimeLo(tlo + ((timelineFAB.box().smallEnd(XDIR) - probDomain.smallEnd(XDIR)) * dt));
  Real fabTimeHi(fabTimeLo + (timelineFAB.box().length(XDIR) * dt));
  Real timeRangeFab(fabTimeHi - fabTimeLo);
  Real ooTimeRangeFab(1.0 / timeRangeFab);

  long index(-1), xi(-1);
  //int nTimeSlotsTotal(probDomain.length(XDIR));
  int nTimeSlotsFab(timelineFAB.box().length(XDIR));
  //int nRanksTotal(probDomain.length(YDIR));
  int fabRankLo((timelineFAB.box().smallEnd(YDIR) - probDomain.smallEnd(YDIR)) * rankStride);
  int fabRankHi(fabRankLo + (timelineFAB.box().length(YDIR) - 1) * rankStride);

  for(int idb(0); idb < dataBlocks.size(); ++idb) {  // ---- go through dataBlocks
    DataBlock &dBlock = dataBlocks[idb];

    if(dBlock.timeMin > fabTimeHi || dBlock.timeMax < fabTimeLo) {  // entire block outside time range
      continue;
    }
    int proc(dBlock.proc);
    if(proc > fabRankHi || proc < fabRankLo || proc % rankStride != 0) {  // block outside rank range
      continue;
    }

    if (persistentStreams){
      ReadCommStatsNoOpen(dBlock);
    } else {
      ReadCommStats(dBlock);
    }
    Real *timeline(timelineFAB.dataPtr(0));
    Real *mpiCount(timelineFAB.dataPtr(1));

    long prevIndex(0);
    for(long i(0); i < dBlock.vCommStats.size(); ++i) {
      BLProfiler::CommStats &cs = dBlock.vCommStats[i];
      Real ts(cs.timeStamp);
      if((ts <= fabTimeHi && ts >= fabTimeLo) && InTimeRange(dBlock.proc, ts)) {  // within time range
        xi = long( nTimeSlotsFab * ((ts - fabTimeLo) * ooTimeRangeFab) );
	if(xi == nTimeSlotsFab) {
	  --xi;
	}
	if(xi < 0 || xi >= nTimeSlotsFab) {
	  SHOWVAL(xi)
	  SHOWVAL(nTimeSlotsFab)
	  SHOWVAL(timeRangeFab)
	  SHOWVAL(ts)
	  amrex::Abort("xi out of range.");
	}
        index = (((proc - fabRankLo) / rankStride) * nTimeSlotsFab) + xi;
	if(index < 0 || index >= timelineFAB.box().numPts()) {
	  SHOWVAL(index)
	  SHOWVAL(timelineFAB.box().size())
	  amrex::Abort("index out of range.");
	}

        timeline[index] = cs.cfType;
        mpiCount[index] += 1.0;

	Real ntnMult(0.0), bnMult(0.0);
	if(cs.cfType == BLProfiler::NameTag) {  // ---- add encoded value for the name tag name
	  ntnMult = ntnNumbers[cs.tag];
          timeline[index] += ntnMult;
	}
	if(cs.cfType == BLProfiler::Barrier) {  // ---- add encoded value for the barrier name
	  if(cs.tag >=  dBlock.barriers.size()) {
	    cout << "******** TimelineFAB::0" << endl;
	    SHOWVAL(i);
	    SHOWVAL(cs.tag);
	    SHOWVAL(dBlock.barriers.size());
	    for(int ib(0); ib < dBlock.barriers.size(); ++ib) {
	      cout << "be[" << ib << "] = " << dBlock.barriers[ib] << endl;
	    }
	    amrex::Abort("--------- bad barrier.");
	  }
	  BarrierEntry &be = dBlock.barriers[cs.tag];
	  bnMult = bnNumbers[be.bNameNumber];
          timeline[index] += bnMult;
	}
        // now fill in gaps
        if(i > 0 && xi < nTimeSlotsFab) {
          BLProfiler::CommStats &csPrev = dBlock.vCommStats[i-1];
	  if(cs.cfType == csPrev.cfType
	     &&
	     (
	       (
	         csPrev.size == BLProfiler::BeforeCall()
		 ||
	         cs.size     == BLProfiler::AfterCall()
	       )
	       ||
	       (
	         csPrev.tag == BLProfiler::BeforeCall()
		 &&
	         cs.tag     == BLProfiler::AfterCall()
	       )

	     )
	    )
	  {
            Real prevTs(csPrev.timeStamp);
	    if(prevTs < fabTimeLo) {
	      cout << "::::  prevTs fabTimeLo = " << prevTs << "  " << fabTimeLo << endl;
	      prevTs = fabTimeLo;
	    }
            long prevXi = long(nTimeSlotsFab * ((prevTs - fabTimeLo) * ooTimeRangeFab));
            prevIndex = (((proc - fabRankLo) / rankStride) * nTimeSlotsFab) + prevXi;
	    for(long idx(prevIndex); idx < index; ++idx) {
	      if(idx < 0 || idx >= timelineFAB.box().numPts() || idx > index) {
	        SHOWVAL(proc)
	        SHOWVAL(ts)
	        SHOWVAL(prevTs)
	        SHOWVAL(fabTimeLo)
	        SHOWVAL(prevXi)
	        SHOWVAL(idx)
	        SHOWVAL(index)
	        SHOWVAL(prevIndex)
	        SHOWVAL(timelineFAB.box().size())
	        //amrex::Abort("idx out of range.");
	        amrex::Print() << "CommProfStats::TimelineFAB::idx out of range." << std::endl;
		continue;
	      }
              timeline[idx] = cs.cfType;
              mpiCount[idx] += 1.0;
	      if(cs.cfType == BLProfiler::NameTag) {
                timeline[idx] += ntnMult;
	      }
	      if(cs.cfType == BLProfiler::Barrier) {
                timeline[idx] += bnMult;
	      }
	    }
	  }
        }
        prevIndex = index;
      }
    }

    ClearCommStats(dBlock);
  }
}


// ----------------------------------------------------------------------
bool CommProfStats::IsSend(const BLProfiler::CommFuncType &cft) {
  return(cft == BLProfiler::AsendTsii  ||
         cft == BLProfiler::AsendTsiiM ||
         cft == BLProfiler::AsendvTii  ||
	 cft == BLProfiler::SendTsii   ||
	 cft == BLProfiler::SendvTii);
}


// ----------------------------------------------------------------------
bool CommProfStats::IsRecv(const BLProfiler::CommFuncType &cft) {
  return(cft == BLProfiler::ArecvTsii  ||
         cft == BLProfiler::ArecvTsiiM ||
         cft == BLProfiler::ArecvTii   ||
         cft == BLProfiler::ArecvvTii  ||
         cft == BLProfiler::RecvTsii   ||
         cft == BLProfiler::RecvvTii   ||
         cft == BLProfiler::Waitsome);
}


// ----------------------------------------------------------------------
bool CommProfStats::IsBlockingRecv(const BLProfiler::CommFuncType &cft) {
  return(cft == BLProfiler::RecvTsii   ||
         cft == BLProfiler::RecvvTii   ||
         cft == BLProfiler::Waitsome);
}


// ----------------------------------------------------------------------
void CommProfStats::SendRecvList(std::multimap<Real, SendRecvPairUnpaired> &srMMap)
{
  BL_PROFILE("SendRecvList");

  for(int idb(0); idb < dataBlocks.size(); ++idb) {    // ---- go through dataBlocks
    DataBlock &dBlock = dataBlocks[idb];
    int proc(dBlock.proc);

    ReadCommStats(dBlock);

    for(int i(0); i < dBlock.vCommStats.size(); ++i) {    // ---- find sends and recvs
      BLProfiler::CommStats &cs = dBlock.vCommStats[i];
      if(IsSend(cs.cfType)) {
        if(cs.size != BLProfiler::AfterCall()) {
          if(InTimeRange(proc, cs.timeStamp)) {
            srMMap.insert(std::pair<Real, SendRecvPairUnpaired>(cs.timeStamp,
	                         SendRecvPairUnpaired(cs.cfType, proc, cs.commpid,
                                                      cs.size, cs.tag, cs.timeStamp,
					              -1000.0, 0.0)));
          }
        }
      }
      if(IsRecv(cs.cfType)) {
        if(cs.size != BLProfiler::AfterCall()) {
          if(InTimeRange(proc, cs.timeStamp)) {
            srMMap.insert(std::pair<Real, SendRecvPairUnpaired>(cs.timeStamp,
	                          SendRecvPairUnpaired(cs.cfType, cs.commpid, proc,
                                                       cs.size, cs.tag, -2000.0,
						       cs.timeStamp, 0.0)));
          }
        }
      }
    }
    ClearCommStats(dBlock);
  }
}


// ----------------------------------------------------------------------
void CommProfStats::SendRecvData(const std::string &filenameprefix,
                                  const double tlo, const double thi)
{
  amrex::ignore_unused(filenameprefix, tlo, thi);

  double dstart(amrex::ParallelDescriptor::second());
  Real timeMin(std::numeric_limits<Real>::max());
  Real timeMax(-std::numeric_limits<Real>::max());

  Vector<Vector<Real> > sendCallTimes(dataNProcs);  // [proc, bnum]
  Vector<Vector<Real> > recvCallTimes(dataNProcs);  // [proc, bnum]
  Vector<Real> dataSentPerProc(dataNProcs * dataNProcs, 0.0);  // [fromproc, toproc]
  Vector<unordered_map<int, unordered_multimap<long, SendRecvPairUnpaired> > > unpairedMessages;
                                                                       // [tag,[hash,msges]]
  Vector<Vector<SendRecvPair> > pairedMessages;
  //unpairedMessages.rehash(minComps);

  //double st(0.0);
  double s0(0.0), s1(0.0), s2(0.0), s3(0.0), rtime(0.0), ptime(0.0);
  long maxupmsize(0), upmsize(0), maxupmmapsize(0);
  Vector<long> upmsizeV;
  int maxTag(-1), maxNMatches(-1);
  //bool matches(false);
  long mTotal(0), mNeg(0);
  //float maxLF0(0.0), maxLF1(0.0);
  Real filterTLo(00.0), filterTHi(30000.0);

#ifdef _OPENMP
  int nReads(4);
  Vector<omp_lock_t> locks(nReads);
  for(int i(0); i < locks.size(); ++i) {
    omp_init_lock(&(locks[i]));
  }
  int nThreads(omp_get_max_threads());
  int myThread(omp_get_thread_num());
#else
  int nThreads(1);
  int myThread(0);
#endif

  procNodeNumber.resize(dataNProcs);

  Vector<DataBlock> dBlockV(nThreads);

  unpairedMessages.resize(nThreads);
  pairedMessages.resize(nThreads);
  for(int i(0); i < pairedMessages.size(); ++i) {
    pairedMessages[i].reserve(2000000);
  }

  int anyDataLeft(dataBlocks.size());
  vector<bool> dataLeft(dataBlocks.size(), true);

  // ---- this part orders the indicies by time, then processor number
  Vector<int> idbIndex;
  bool resort(true);
  if(resort) {
    map<int, multimap<int, int> > idbMM;    // [time, [proc, index]]
    map<int, multimap<int, int> >::iterator idb_m_MMiter;
    multimap<int, int>::iterator idbMMiter;
    //cout << "==============" << endl;
    int countI(0);
    for(int idb(0); idb < dataBlocks.size(); ++idb) {
      //cout << "idb dB.proc = " << idb << "  " << dataBlocks[idb].proc << endl;
      //cout << idb << "  " << dataBlocks[idb].proc << endl;
      int t0 = int(dataBlocks[idb].timeMin * 1.0);
      idbMM[t0].insert(std::pair<int, int>(dataBlocks[idb].proc, idb));
      ++countI;
    }
    //cout << "%%%%%%%%%%%%%%" << endl;
    int count(0);
    for(idb_m_MMiter = idbMM.begin(); idb_m_MMiter != idbMM.end(); ++idb_m_MMiter) {
      multimap<int, int> &mm = idb_m_MMiter->second;
      //cout << "---- " << idb_m_MMiter->first << endl;
      for(idbMMiter = mm.begin(); idbMMiter != mm.end(); ++idbMMiter) {
        //cout << "idb idbMM.proc = " << idbMMiter->second << "  " << idbMMiter->first << endl;
        idbIndex.push_back(idbMMiter->second);
        ++count;
      }
    }
    //cout << "==============  count countI = " << count << "  " << countI << endl;
    //cout << "**************  idbIndex.size() = " << idbIndex.size() << endl;
    //for(idb = 0; idb < idbIndex.size(); ++idb) {
      //cout << idb << "  " << idbIndex[idb] << "  " << dataBlocks[idbIndex[idb]].timeMin
           //<< "  " << dataBlocks[idbIndex[idb]].size << endl;
    //}
    //cout << "**************" << endl;
  } else {
    idbIndex.resize(dataBlocks.size());
    for(int idbII = 0; idbII < idbIndex.size(); ++idbII) {
      idbIndex[idbII] = idbII;
    }
  }


  while(anyDataLeft) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
  //for(idb = 0; idb < dataBlocks.size(); ++idb) {
  for(int idbII = 0; idbII < dataBlocks.size(); ++idbII) {
    int idb(idbIndex[idbII]);
    //int mySet(myThread % nReads);
    //int mySet(idb % nReads);

    DataBlock &dB = dataBlocks[idb];
    DataBlock &dBlock = dBlockV[myThread];
    dBlock.proc = dB.proc;
    dBlock.size = dB.size;
    dBlock.fileName = dB.fileName;
    dBlock.seekpos = dB.seekpos;
    dBlock.readoffset = dB.readoffset;
    //DataBlock dBlock(dB.proc, dB.size, dB.fileName, dB.seekpos);

    procNodeNumber[dBlock.proc] = dB.nodeNumber;

    timeMin = std::min(timeMin, dB.timeMin);
    timeMax = std::max(timeMax, dB.timeMax);


    Real pctDone(100.0 - 100.0 * ((Real)(dB.size - dB.readoffset) / (Real) dB.size));
    cout << "Processing data block " << idb + 1 << " : " << idbII + 1 << " / "
         << dataBlocks.size()
         << "   pctDone = " << pctDone << endl;

    if(dBlock.timeMin > filterTHi || dBlock.timeMax < filterTLo) {  // block outside range
      continue;
    }
    double rstart(amrex::ParallelDescriptor::second());
  //  omp_set_lock(&(locks[mySet]));    // sort of a semaphore
      //ReadCommStats(dBlock);

      //int readSize(64 + dBlock.size/2);
      //int readSize(dBlock.size);
      int readSize(1000000);

      if(dataLeft[idb]) {
        dataLeft[idb] = ReadCommStats(dBlock, readSize);
        if(dataLeft[idb] == false) {
	  --anyDataLeft;
	}
      } else {
        continue;
      }

      dB.readoffset = dBlock.readoffset;
  //  omp_unset_lock(&(locks[mySet]));
    rtime += amrex::ParallelDescriptor::second() - rstart;

    Real sttemp(0.0), rttemp(0.0), tttemp(0.0);
    int proc(dBlock.proc);
    //int prevIndex(0);

    for(int ics(0); ics < dBlock.vCommStats.size(); ++ics) {
      BLProfiler::CommStats &cs = dBlock.vCommStats[ics];
      Real ts(cs.timeStamp);

      if(ts > filterTHi || ts < filterTLo) {
        continue;
      }
      if(IsBlockingRecv(cs.cfType) && cs.size == BLProfiler::BeforeCall()) {
        continue;
      }
      if(IsSend(cs.cfType) && cs.size == BLProfiler::AfterCall()) {
        continue;
      }
      if(cs.tag == BLProfiler::NoTag()) {
        continue;
      }
      if(cs.tag < 1000) {
        continue;
      }

      if(IsSend(cs.cfType) || IsBlockingRecv(cs.cfType)) {
        double pstart(amrex::ParallelDescriptor::second());

double st0 = amrex::ParallelDescriptor::second();

        int fromProcTemp, toProcTemp;
	if(IsSend(cs.cfType)) {
	  fromProcTemp = proc;
	  toProcTemp   = cs.commpid;
	} else {
	  fromProcTemp = cs.commpid;
	  toProcTemp   = proc;
	}
	//long hlong(SendRecvPair::HashLong(fromProcTemp, toProcTemp, cs.size, cs.tag));
	long hlong = SendRecvPair::HashLong(fromProcTemp, toProcTemp, cs.tag);

	//bool foundMatch(false);
	int  nMatches(0);

	pair<unordered_multimap<long, SendRecvPairUnpaired>::iterator,
	     unordered_multimap<long, SendRecvPairUnpaired>::iterator> upmSRPERI;
	Vector<unordered_multimap<long, SendRecvPairUnpaired>::iterator> upmSRPMatchSave;

	maxTag = std::max(maxTag, cs.tag);

	unordered_multimap<long, SendRecvPairUnpaired> &upm = unpairedMessages[myThread][cs.tag];
	//maxupmmapsize = std::max(static_cast<long>(unpairedMessages[myThread].size()), maxupmmapsize);
	upmSRPERI = upm.equal_range(hlong);
	for(unordered_multimap<long, SendRecvPairUnpaired>::iterator upmsrit = upmSRPERI.first;
	    upmsrit != upmSRPERI.second; ++upmsrit)
	{
	  SendRecvPairUnpaired &srp = upmsrit->second;
	  if(srp.Matches(fromProcTemp, toProcTemp, cs.size, cs.tag)) {
	    if((IsSend(srp.unmatchedCFType) && IsBlockingRecv(cs.cfType)) || 
	       (IsSend(cs.cfType) && IsBlockingRecv(srp.unmatchedCFType)))
	    {
	      upmSRPMatchSave.push_back(upmsrit);
	      ++nMatches;
	      maxNMatches = std::max(maxNMatches, nMatches);
	    }
	  }
	}
s0 += amrex::ParallelDescriptor::second() - st0;

        if(ics+1 >= dBlock.vCommStats.size()) {
	  continue;
	}

	if(nMatches == 0) {      // ------------------------ unpaired
double st1 = amrex::ParallelDescriptor::second();
	  if(IsSend(cs.cfType)) {
	    sttemp = cs.timeStamp;
	    tttemp = sttemp;
	  } else {
	    if(ics+1 >= dBlock.vCommStats.size()) {
	      cout << "!!!!!!!!!!!!!!!!!!!!!!! cs on end (unpaired)." << endl;
	    }
            BLProfiler::CommStats &csNext = dBlock.vCommStats[ics+1];
	    rttemp = csNext.timeStamp;
	    tttemp = rttemp;
	  }
	  upm.insert(std::pair<long, SendRecvPairUnpaired>(hlong, 
	                                      SendRecvPairUnpaired(cs.cfType, fromProcTemp,
						toProcTemp, cs.size, cs.tag,
						sttemp, rttemp, tttemp)));
	  ++upmsize;
	  maxupmsize = std::max(upmsize, maxupmsize);
s1 += amrex::ParallelDescriptor::second() - st1;

	} else {                 // ------------------------ paired

double st2 = amrex::ParallelDescriptor::second();
	  unordered_multimap<long, SendRecvPairUnpaired>::iterator
	                                  earliestMatchSave = upmSRPMatchSave[0];
	  for(int i(1); i < upmSRPMatchSave.size(); ++i) {
	    if(upmSRPMatchSave[i]->second.totalTime < earliestMatchSave->second.totalTime) {
	      earliestMatchSave = upmSRPMatchSave[i];
	    }
	  }
	  SendRecvPairUnpaired &srpup = earliestMatchSave->second;
	  sttemp = srpup.sendTime;
	  rttemp = srpup.recvTime;
          if(IsSend(cs.cfType)) {
	    sttemp = cs.timeStamp;
	  } else {
	    if(ics+1 >= dBlock.vCommStats.size()) {
	      cout << "!!!!!!!!!!!!!!!!!!!!!!! cs on end (paired)." << endl;
	    }
            BLProfiler::CommStats &csNext = dBlock.vCommStats[ics + 1];
	    rttemp = csNext.timeStamp;
	  }
	  tttemp = rttemp - sttemp;
          if(tttemp < 0.0) {
            ++mNeg;
	  }
	  ++mTotal;
	  pairedMessages[myThread].push_back(SendRecvPair(srpup.fromProc, srpup.toProc,
	                                        srpup.dataSize, srpup.tag,
					        sttemp, rttemp, tttemp));
	  upm.erase(earliestMatchSave);
	  --upmsize;

s2 += amrex::ParallelDescriptor::second() - st2;
	}  // end else paired

    ptime += amrex::ParallelDescriptor::second() - pstart;

      }  // end if(IsSend...)

    }  // end for(ics...)
//st = amrex::ParallelDescriptor::second();
    //ClearCommStats(dBlock);
//s3 += amrex::ParallelDescriptor::second() - st;

    //SHOWVAL(unpairedMessages[myThread].bucket_count());
    //SHOWVAL(unpairedMessages[myThread].load_factor());
    SHOWVAL(upmsize);

    upmsizeV.push_back(upmsize);


  }  // end for(idb...)

  } // end while




#ifdef _OPENMP
  for(int i(0); i < locks.size(); ++i) {
    omp_destroy_lock(&(locks[i]));
  }
#endif

  cout << endl;
  SHOWVAL(maxTag);
  SHOWVAL(maxNMatches);
  SHOWVAL(mTotal);
  SHOWVAL(mNeg);
  double pNeg(static_cast<double>(mNeg) / static_cast<double>(mTotal));
  SHOWVAL(pNeg);


  SHOWVAL(upmsize);
  SHOWVAL(maxupmsize);
  SHOWVAL(maxupmmapsize);
  SHOWVAL(s0);
  SHOWVAL(s1);
  SHOWVAL(s2);
  SHOWVAL(s3);
  cout << "%%%%%%%%%%%%% SendRecvData pairing time = " << amrex::ParallelDescriptor::second() - dstart << endl;
  SHOWVAL(rtime);
  SHOWVAL(ptime);
  SHOWVAL(ptime-s0-s1-s2-s3);
  //SHOWVAL(maxLF0);
  //SHOWVAL(maxLF1);


//return;

  std::ofstream upmout("upmsize.xgr");
  for(int i(0); i < upmsizeV.size(); ++i) {
    upmout << i << " " << upmsizeV[i] << endl;
  }
  upmout.close();


// =================================================

  // ------------------------------------------------ print send times
  /*
  std::ofstream tsdout("sendtimes.xgr");
  for(int i(0); i < pairedMessages[0].size(); ++i) {
    SendRecvPair &srp = pairedMessages[i];
    tsdout << i << " " << srp.totalTime << endl;
  }
  tsdout.close();
  */

  {
    int index;
    Vector<double> sdataTT(dataNProcs * dataNProcs, 0.0);
    Vector<double> sdataDS(dataNProcs * dataNProcs, 0.0);
    Vector<double> sdataNC(dataNProcs * dataNProcs, 0.0);
    Vector<double> sdataNHops(dataNProcs * dataNProcs, 0.0);
  SHOWVAL(pairedMessages[0].size());
    cout << endl << "---------------------------------- send receive pairs" << endl;
    for(int i(0); i < pairedMessages[0].size(); ++i) {
      SendRecvPair &srp = pairedMessages[0][i];
      index = srp.fromProc + srp.toProc * dataNProcs;
      sdataTT[index] += srp.totalTime;
      sdataDS[index] += srp.dataSize;
      //cout << "srp:  i fromProc toProc dataSize = " << i << "  " <<  srp.fromProc
           //<< "  " << srp.toProc << "  " << srp.dataSize << "  " << srp.sendTime << endl;
      sdataNC[index] += 1.0;  // count

      if(bTopoFabInited) {
	IntVect ivfrom(TopIVFromProcNumber(procNodeNumber[srp.fromProc]));
	IntVect ivto(TopIVFromProcNumber(procNodeNumber[srp.toProc]));
	int nhops(NHops(topoBox, ivfrom, ivto));
        sdataNHops[index] += nhops;
      }

    }
    Write2DFab("srTotalTime", dataNProcs, dataNProcs, &sdataTT[0]);
    Write2DFab("srTotalData", dataNProcs, dataNProcs, &sdataDS[0]);
    Write2DText("srTotalData", dataNProcs, dataNProcs, &sdataDS[0]);
    Write2DFab("srTotalMsgs", dataNProcs, dataNProcs, &sdataNC[0]);
    Write2DText("srTotalMsgs", dataNProcs, dataNProcs, &sdataNC[0]);
    if(bTopoFabInited) {
      Write2DFab("srTotalNHops", dataNProcs, dataNProcs, &sdataNHops[0]);
    }

  cout << endl << "---------------------------------- total sent data per proc" << endl;
  /*
  for(int j(dataNProcs-1); j >= 0; --j) {
    for(int i(0); i < dataNProcs; ++i) {
      cout << sdataDS[i + j * dataNProcs] << "  ";
    }
    cout << '\n';
  }
  cout << endl;
  */
  }

  {
    cout << endl << "---------------------------------- total time time view" << endl;
    Real timeMinZ(0.0);
    Real timeMaxZ(timeMax + 0.001);
    Real timeMaxRange(timeMaxZ - timeMinZ);
    int ntimeslots = (int) timeMaxRange;
    ntimeslots = std::max(512, ntimeslots);
    //Real dt(timeMaxRange / (Real) ntimeslots);
    int index(-1), zi(-1);
    SHOWVAL(timeMinZ);
    SHOWVAL(timeMax);
    SHOWVAL(timeMaxZ);
    SHOWVAL(timeMaxRange);
    SHOWVAL(ntimeslots);
    SHOWVAL(dataNProcs);

    Vector<double> ttData(dataNProcs * dataNProcs * ntimeslots, 0.0);
    Vector<double> bwData(dataNProcs * dataNProcs * ntimeslots, 0.0);
    for(int i(0); i < pairedMessages[0].size(); ++i) {
      SendRecvPair &srp = pairedMessages[0][i];
      Real ts(srp.recvTime);
      zi = (ntimeslots - 1) - (int) ((ntimeslots - 1) * ((timeMaxRange - ts) / timeMaxRange));
      index = srp.fromProc + (dataNProcs * srp.toProc) + (dataNProcs * dataNProcs * zi);
      if(index >= ttData.size()) {
        SHOWVAL(ts);
        SHOWVAL(zi);
        SHOWVAL(index);
	abort();
      }
      ttData[index] += srp.totalTime;
      bwData[index] += srp.dataSize;
      //cout << "srp:  i fromProc toProc dataSize = " << i << "  " <<  srp.fromProc
           //<< "  " << srp.toProc << "  " << srp.dataSize << "  " << srp.sendTime << endl;
    }
    for(int i(0); i < bwData.size(); ++i) {
      if(ttData[i] > 0.0) {
        bwData[i] /= ttData[i];
      }
    }

    Write3DFab("SRP_TT", dataNProcs, dataNProcs, ntimeslots, &ttData[0]);
    Write3DFab("SRP_BW", dataNProcs, dataNProcs, ntimeslots, &bwData[0]);
  }


return;


#if 0

  {
    int index;
    Vector<double> sdataTT(dataNProcs * dataNProcs, 0.0);
    cout << endl << "---------------------------------- time data total calculated" << endl;
    for(int i(0); i < pairedMessages.size(); ++i) {
      SendRecvPair &srp = pairedMessages[i];
      int fp, tp;
      if(srp.fromProc > srp.toProc) {
        fp = srp.fromProc;
        tp = srp.toProc;
      } else {
        fp = srp.toProc;
        tp = srp.fromProc;
      }
      index = fp + tp * dataNProcs;
      sdataTT[index] += srp.totalTime;
    }
    Write2DFab("timeDataCalc", dataNProcs, dataNProcs, &sdataTT[0]);
  }

  // ------------------------------------------------ unskew
  if(bExitTimesDone) {
    for(int i(0); i < pairedMessages.size(); ++i) {
      SendRecvPair &srp = pairedMessages[i];
      //tsdout << i << " " << srp.totalTime << endl;
      //srp. sendTime recvTime fromProc toProc...
      int sB = AfterBarrier(srp.fromProc, srp.sendTime);
      int rB = AfterBarrier(srp.toProc,   srp.recvTime);
      if(sB != rB) {
        cout << "**** Error:  over barrier:  sB rB = " << sB << "  " << rB << endl;
      }
      if(i < 32) {
        cout << "unskew:  fP sT  tP rT  betS betR s-r = " << srp.fromProc << "  " << srp.sendTime
             << "  >  " << srp.toProc << "  " << srp.recvTime << "    "
	     << barrierExitTimes[srp.fromProc][sB]
	     << "  " << barrierExitTimes[srp.toProc][rB] << "    "
	     << barrierExitTimes[srp.fromProc][sB] - barrierExitTimes[srp.toProc][rB] << endl;
      }
      //srp.totalTime += barrierExitTimes[srp.fromProc][sB] - barrierExitTimes[srp.toProc][rB];
      srp.totalTime = (srp.recvTime - barrierExitTimes[srp.toProc][rB]) -
                      (srp.sendTime - barrierExitTimes[srp.fromProc][sB]);
    }
  }

  {
    int index;
    Vector<double> sdataTT(dataNProcs * dataNProcs, 0.0);
    cout << endl << "---------------------------------- send receive pairs" << endl;
    for(int i(0); i < pairedMessages.size(); ++i) {
      SendRecvPair &srp = pairedMessages[i];
      index = srp.fromProc + srp.toProc * dataNProcs;
      sdataTT[index] += srp.totalTime;
    }
    Write2DFab("srTotalTimeUnskewed", dataNProcs, dataNProcs, &sdataTT[0]);
  }

  {
    int index;
    Vector<double> sdataTT(dataNProcs * dataNProcs, 0.0);
    cout << endl << "---------------------------------- time data total calculated" << endl;
    for(int i(0); i < pairedMessages.size(); ++i) {
      SendRecvPair &srp = pairedMessages[i];
      int fp, tp;
      if(srp.fromProc > srp.toProc) {
        fp = srp.fromProc;
        tp = srp.toProc;
      } else {
        fp = srp.toProc;
        tp = srp.fromProc;
      }
      index = fp + tp * dataNProcs;
      sdataTT[index] += srp.totalTime;
    }
    Write2DFab("timeDataCalcUnskewed", dataNProcs, dataNProcs, &sdataTT[0]);
  }


//return;


  sendCallTimes.reserve(2000000);
  recvCallTimes.reserve(2000000);
  cout << "Processing " << dataBlocks.size() << " data blocks:" << endl;
  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    cout << idb + 1 << "  " << flush;
    if((idb + 1) % 20 == 0) {
      cout << endl;
    }
    DataBlock &dBlock = dataBlocks[idb];
    ReadCommStats(dBlock);

    int proc(dBlock.proc);
    int prevIndex(0);
    for(int i(0); i < dBlock.vCommStats.size(); ++i) {
      BLProfiler::CommStats &cs = dBlock.vCommStats[i];

      // time for send call to return
      if(IsSend(cs.cfType)) {
        index = proc + dataNProcs * cs.commpid;
        dataSentPerProc[index] += cs.size;
#ifdef DEBUG
        if(i+1 >=  dBlock.vCommStats.size()) {
          cerr << "**** Error in SendRecvData:  bad i." << endl;
          abort();
        }
#endif
        BLProfiler::CommStats &csNext = dBlock.vCommStats[i+1];
        sendCallTimes[proc].push_back(csNext.timeStamp - cs.timeStamp);
        ++i;
        continue;
      }

      // time for recv call to return
      if(IsRecv(cs.cfType)) {
#ifdef DEBUG
        if(i+1 >=  dBlock.vCommStats.size()) {
          cerr << "**** Error in SendRecvData:  bad i." << endl;
          abort();
        }
#endif
        BLProfiler::CommStats &csNext = dBlock.vCommStats[i+1];
        recvCallTimes[proc].push_back(csNext.timeStamp - cs.timeStamp);
        ++i;
        continue;
      }
    }

    ClearCommStats(dBlock);
  }
  cout << endl;

  // ------------------------------------------------ print send call times
  {
    int maxSends(0);
    for(int i(0); i < sendCallTimes.size(); ++i) {
      maxSends = std::max(maxSends, static_cast<int> (sendCallTimes[i].size()));
    }
    SHOWVAL(maxSends);
    Vector<double> sdata(maxSends * dataNProcs, 0.0);

    for(int ip(0); ip < dataNProcs; ++ip) {
      for(int is(0); is < sendCallTimes[ip].size(); ++is) {
        index = is + ip * maxSends;
        sdata[index] = sendCallTimes[ip][is];
      }
    }
    cout << endl << "writing send wait times." << endl;
    WriteFab("sendCallTimes", maxSends, dataNProcs, &sdata[0]);
  }


  // ------------------------------------------------ print recv call times
  {
    int maxRecvs(0);
    for(int i(0); i < recvCallTimes.size(); ++i) {
      maxRecvs = std::max(maxRecvs, static_cast<int> (recvCallTimes[i].size()));
    }
    SHOWVAL(maxRecvs);
    Vector<double> rdata(maxRecvs * dataNProcs, 0.0);

    for(int ip(0); ip < dataNProcs; ++ip) {
      for(int is(0); is < recvCallTimes[ip].size(); ++is) {
        index = is + ip * maxRecvs;
        rdata[index] = recvCallTimes[ip][is];
      }
    }
    cout << endl << "writing recv wait times." << endl;
    WriteFab("recvCallTimes", maxRecvs, dataNProcs, &rdata[0]);
  }

  Write2DFab("dataSentPerProc", dataNProcs, dataNProcs, &dataSentPerProc[0]);

  cout << "%%%%%%%%%%%%% SendRecvData time = "
       << amrex::ParallelDescriptor::second() - dstart << endl;


#endif

}


// ----------------------------------------------------------------------
void CommProfStats::InitEdisonTopoMF() {
  edisonNodeFab.resize(32 * 24 * 8, -1.0);
  edisonCPUFab.resize(64 * 24 * 96, -2.0);
  edisonNodeXYZ.resize(16);
  for(int ix(0); ix < edisonNodeXYZ.size(); ++ix) {
    edisonNodeXYZ[ix].resize(6);
    for(int iy(0); iy < edisonNodeXYZ[ix].size(); ++iy) {
      edisonNodeXYZ[ix][iy].resize(16);
    }
  }
}


// ----------------------------------------------------------------------
void CommProfStats::WriteEdisonTopoMF() {

  Write3DFab("edisonNode", 32, 24, 8, edisonNodeFab.dataPtr());

  int ip(0);
  std::ofstream pnnout("edisonNodeXYZ.xgr");
  for(int ix(0); ix < edisonNodeXYZ.size(); ++ix) {
    edisonNodeXYZ[ix].resize(6);
    for(int iy(0); iy < edisonNodeXYZ[ix].size(); ++iy) {
      for(int iz(0); iz < edisonNodeXYZ[ix][iy].size(); ++iz) {
        cout << "XYZ pid = " << ix << "  " << iy << "  " << iz << "  "
	     << edisonNodeXYZ[ix][iy][iz] << endl;
        pnnout << ip++ << " " << edisonNodeXYZ[ix][iy][iz] << '\n';
      }
    }
  }
  pnnout.close();
}


// ----------------------------------------------------------------------
void CommProfStats::AddEdisonPID(int X, int Y, int Z,
                                 int cab, int row, int cage, int slot,
                                 int cpu, int pid)
{
  amrex::ignore_unused(cage, slot);

  int ix, iy, iz, index;
  int ixGroup, ixCab, ixSlot, ixNode, iyCage, iySlot, izGroup, izNode;

  edisonNodeXYZ[X][Y][Z] = pid;

  ixGroup = X / 4;
  if(ixGroup != row) {
    cout << "ixGroup row = " << ixGroup << "  " << row << endl;
  }
  ixCab = cab % 2;
  ixSlot = Z / 8;
  ixNode = cpu % 2;
  ix = (8 * ixGroup) + (4 * ixCab) + (2 * ixSlot) + ixNode;

  iyCage = Y % 3;
  iySlot = Z % 8;
  iy = (8 * iyCage) + iySlot;

  izGroup = X % 4;
  izNode = cpu / 2;
  iz = (2 * izGroup) + izNode;

  index = (iz * 32 * 24) + (iy * 32) + ix;

  Real *ptr = edisonNodeFab.dataPtr();
  ptr[index] = pid;
  if(X == 11 && Y == 3) {
    cout << "AddEdisonPID::113:  pid index ix iy iz X Y Z = " << pid << "  " << index << "  "
         << ix << "  " << iy << "  " << iz << "  "
         << X << "  " << Y << "  " << Z << endl;
  }
  if(X == 3 && Y == 0) {
    cout << "AddEdisonPID::30:  pid index ix iy iz = " << pid << "  " << index << "  "
         << ix << "  " << iy << "  " << iz << endl;
  }
}




// ----------------------------------------------------------------------
std::ostream &operator<<(std::ostream &os, const CommProfStats::BarrierEntry &be)
{
  os << "BE:  num seek name bnn = " << be.number << ' ' << be.seekIndex
     << ' ' << be.name << ' ' << be.bNameNumber << endl;
  return os;
}


// ----------------------------------------------------------------------
CommProfStats::BarrierEntry &CommProfStats::BarrierEntry::operator=(const CommProfStats::BarrierEntry &be)
{
  number = be.number;
  name = be.name;
  seekIndex = be.seekIndex;
  bNameNumber = be.bNameNumber;
  return *this;
}
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
