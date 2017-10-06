// -------------------------------------------------------------
// Pieces.cpp
// -------------------------------------------------------------
#include <AMReX_Vector.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_BLProfiler.H>
#include <RegionsProfStats.H>
#include <CommProfStats.H>
#include <iostream>
#include <strstream>
#include <fstream>
#include <iomanip>

#include <unistd.h>

using std::cout;
using std::endl;
using std::ostrstream;
using std::ofstream;
using std::streamoff;

const int XDIR(0);
const int YDIR(1);
const int ZDIR(2);


extern void PrintTimeRangeList(const std::list<RegionsProfStats::TimeRange> &trList);


// ----------------------------------------------------------------------
void TestRemovePiece(std::list<RegionsProfStats::TimeRange> &removeFromHere,
                     const RegionsProfStats::TimeRange &pieceToRemove)
{
  RegionsProfStats::RemovePiece(removeFromHere, pieceToRemove);
  cout << "Removing " << pieceToRemove << endl;
  cout << "  removeFromHere =";
  PrintTimeRangeList(removeFromHere);
}


// ----------------------------------------------------------------------
void TestAddPiece(std::list<RegionsProfStats::TimeRange> &addToHere,
                     const RegionsProfStats::TimeRange &pieceToAdd)
{
  RegionsProfStats::AddPiece(addToHere, pieceToAdd);
  cout << "Adding " << pieceToAdd << endl;
  cout << "  addToHere =";
  PrintTimeRangeList(addToHere);
}


// ----------------------------------------------------------------------
void TestRangeIntersection(std::list<RegionsProfStats::TimeRange> &rangeList,
                           const RegionsProfStats::TimeRange &pieceToIntersect)
{
  std::list<RegionsProfStats::TimeRange> intersectList =
               RegionsProfStats::RangeIntersection(rangeList, pieceToIntersect);
  cout << "Intersecting " << pieceToIntersect << endl;
  cout << "  intersectList =";
  PrintTimeRangeList(intersectList);
}



// -------------------------------------------------------------
int main(int argc, char *argv[]) {

  amrex::Initialize(argc,argv);

  BLProfiler::SetBlProfDirName("blpp_prof");
  BL_PROFILE_VAR("main()", pmain);

  bool bIOP(ParallelDescriptor::IOProcessor());
  int myproc(ParallelDescriptor::MyProc());
  int nprocs(ParallelDescriptor::NProcs());

  if(bIOP) {
    cout << "**** sizeof(bool)        = " << sizeof(bool) << endl;
    cout << "**** sizeof(Real)        = " << sizeof(Real) << endl;
    cout << "**** sizeof(std::string) = " << sizeof(std::string) << endl;
    cout << "**** sizeof(BLProfiler)    = " << sizeof(BLProfiler) << endl;
    cout << "**** sizeof(BLProfiler::ProfStats)    = " << sizeof(BLProfiler::ProfStats) << endl;
    cout << "**** sizeof(BLProfiler::CallStats)    = " << sizeof(BLProfiler::CallStats) << endl;
    cout << "**** sizeof(BLProfiler::CommStats)    = " << sizeof(BLProfiler::CommStats) << endl;
    cout << "**** sizeof(BLProfStats)      = " << sizeof(BLProfStats) << endl;
    cout << "**** sizeof(CommProfStats)    = " << sizeof(CommProfStats) << endl;
    cout << "**** sizeof(RegionsProfStats) = " << sizeof(RegionsProfStats) << endl;
  }

  // test the RemovePiece function
  {
    cout << endl << "======================== testing RemovePiece" << endl;
    std::list<RegionsProfStats::TimeRange> removeFromHere;
    std::list<RegionsProfStats::TimeRange>::iterator it;
    removeFromHere.push_back(RegionsProfStats::TimeRange(0.0, 10.0));
    cout << "  removeFromHere =";
    PrintTimeRangeList(removeFromHere);

    TestRemovePiece(removeFromHere, RegionsProfStats::TimeRange(4.2, 6.4));
    TestRemovePiece(removeFromHere, RegionsProfStats::TimeRange(1.3, 2.7));
    TestRemovePiece(removeFromHere, RegionsProfStats::TimeRange(3.0, 3.9));
    TestRemovePiece(removeFromHere, RegionsProfStats::TimeRange(2.4, 3.2));
    TestRemovePiece(removeFromHere, RegionsProfStats::TimeRange(-10.5, -3.2));
    TestRemovePiece(removeFromHere, RegionsProfStats::TimeRange(-1.7, 0.9));
    TestRemovePiece(removeFromHere, RegionsProfStats::TimeRange(9.42, 10.2));

    cout << "============================================" << endl << endl;
  }

  // test the AddPiece function
  std::list<RegionsProfStats::TimeRange> addToHere;
  {
    cout << endl << "======================== testing AddPiece" << endl;
    std::list<RegionsProfStats::TimeRange>::iterator it;
    cout << "  addToHere =";
    PrintTimeRangeList(addToHere);

    TestAddPiece(addToHere, RegionsProfStats::TimeRange(4.2, 6.4));
    TestAddPiece(addToHere, RegionsProfStats::TimeRange(1.3, 2.7));
    TestAddPiece(addToHere, RegionsProfStats::TimeRange(3.0, 3.9));
    TestAddPiece(addToHere, RegionsProfStats::TimeRange(2.6, 2.8));

    cout << "============================================" << endl << endl;
    {
      std::list<RegionsProfStats::TimeRange> addToHere;
      RegionsProfStats::TimeRange trAll(0.0, 10.0);
      addToHere.push_back(trAll);
      cout << "  addToHere ="; PrintTimeRangeList(addToHere);
      TestAddPiece(addToHere, RegionsProfStats::TimeRange(4.2, 6.4));
      cout << "  addToHere ="; PrintTimeRangeList(addToHere);
    }
    cout << "============================================" << endl << endl;
  }

  // test the RangeIntersection function
  {
    cout << endl << "======================== testing RangeIntersection" << endl;
    std::list<RegionsProfStats::TimeRange> rangeList = addToHere;
    std::list<RegionsProfStats::TimeRange>::iterator it;
    cout << "  rangeList =";
    PrintTimeRangeList(rangeList);

    TestRangeIntersection(rangeList, RegionsProfStats::TimeRange(2.5, 3.1));

    cout << "============================================" << endl << endl;
  }

  // test the RangeIntersection function
  {
    cout << endl << "======================== testing sort TimeRange list" << endl;
    std::list<RegionsProfStats::TimeRange> trSorted;
    trSorted.push_back(RegionsProfStats::TimeRange(8.0, 10.0));
    trSorted.push_back(RegionsProfStats::TimeRange(2.0, 3.0));
    trSorted.push_back(RegionsProfStats::TimeRange(2.5, 5.6));
    trSorted.push_back(RegionsProfStats::TimeRange(4.0, 5.0));
    trSorted.push_back(RegionsProfStats::TimeRange(-1.0, 0.7));
    cout << "  trSorted (orig) ="; PrintTimeRangeList(trSorted);
    trSorted.sort(BLProfStats::TimeRangeCompare());
    cout << "  trSorted (sorted) ="; PrintTimeRangeList(trSorted);
    TestAddPiece(trSorted, RegionsProfStats::TimeRange(2.6, 2.8));
    cout << "  trSorted (sorted) ="; PrintTimeRangeList(trSorted);
    TestAddPiece(trSorted, RegionsProfStats::TimeRange(-2.6, 12.8));
    cout << "  trSorted (sorted) ="; PrintTimeRangeList(trSorted);
    cout << "============================================" << endl << endl;
  }

  double dz(0);
  cout << "dz = " << dz << endl;
  cout << std::setiosflags(std::ios::showpoint) << "dz = " << dz << endl;
  int di(1000000);
  cout << "di = " << di << endl;
  const std::locale oldLoc(cout.std::ios_base::getloc());
  cout.imbue(std::locale(""));
  cout << "di = " << di << endl;
  cout << std::setiosflags(std::ios::showpoint) << "dz = " << dz << endl;
  dz = di;
  cout << std::setiosflags(std::ios::showpoint) << "dz = " << dz << endl;
  cout.imbue(oldLoc);
  cout << "di = " << di << endl;



  BL_PROFILE_VAR_STOP(pmain);

  BLProfiler::SetNoOutput();

  amrex::Finalize();
  return 0;
}
// -------------------------------------------------------------
// -------------------------------------------------------------
