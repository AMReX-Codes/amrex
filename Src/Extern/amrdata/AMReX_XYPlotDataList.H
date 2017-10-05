// -------------------------------------------------------------------
// XYPlotDataList.H
// -------------------------------------------------------------------
#ifndef XYPLOTDATALIST_H
#define XYPLOTDATALIST_H

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_BLassert.H>
#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>

#include <list>
#include <string>
using std::string;
using std::list;

using namespace amrex;

namespace amrex {

class OrderedBoxes {
  public:
    OrderedBoxes(const int level, int idir, const amrex::Box &databox, const amrex::Box &finebox)
      : iLevel(level), iDir(idir), dataBox(databox), finestBox(finebox)
    { }

    int ILevel() const { return iLevel; }
    const amrex::Box &DataBox() const { return dataBox; }
    const amrex::Box &FinestBox() const { return finestBox; }
    bool operator<(const OrderedBoxes &ob) const {
      return(finestBox.smallEnd()[iDir] < ob.finestBox.smallEnd()[iDir]);
    }

  protected:
    int iLevel, iDir;
    amrex::Box dataBox, finestBox;

};


class XYPlotDataListLink {
  public:
    XYPlotDataListLink(const Real *d, int sX, int len)
      : xypdllLength(len), startXi(sX)
    { 
      BL_ASSERT(len != 0);
      endXi = sX + len;
      xypdllData = new Real[len];
      for(int i(0); i < len; ++i) {
        xypdllData[i] = d[i];
      }
    }

    ~XYPlotDataListLink() {
      delete [] xypdllData;
    }

    Real *XYPDLLData() const { return xypdllData; }
    int XYDPLLength() const { return xypdllLength; }
    int StartXi() const { return startXi; }
    int EndXi() const { return endXi; }
    bool operator<(const XYPlotDataListLink &xyp) const {
      return(startXi < xyp.startXi);
    }

  protected:
    Real *xypdllData;
    int xypdllLength, startXi, endXi;
};


class XYPlotDataList {
  public:
    XYPlotDataList(const string &derived, int minlevel, int maxlevel,
		   int gridlinein,
		   const amrex::Vector<int> &ratiolist,
		   const amrex::Vector<Real> &dx,
		   const amrex::Vector<char *> &intersectpoint,
		   Real startx = 0.0);
    XYPlotDataList(XYPlotDataList *src);
    ~XYPlotDataList();

    void AddFArrayBox(amrex::FArrayBox &fab, int direction, int level);

    void UpdateStats(void);
    void SetLevel(int new_level) {
      curLevel = (new_level < maxLevel) ? new_level : maxLevel;
    }
    const string &DerivedName() const { return xypdlDerived; }
    const XYPlotDataList *CopiedFrom() const { return copiedFrom; }
    void SetCopiedFrom(XYPlotDataList *xypdl) { copiedFrom = xypdl; }
    int MaxLevel()  const { return maxLevel; }
    int NumPoints() const { return numPoints[curLevel]; }
    int CurLevel()  const { return curLevel; }
    int Gridline()  const { return gridline; }
    double StartX() const { return startX; }
    double EndX()   const { return endX; }
    double XYPDLLoY(const int ilev) const { return xypdlLoY[ilev]; }
    double XYPDLHiY(const int ilev) const { return xypdlHiY[ilev]; }
    char *IntersectPoint(const int ilev) const { return intersectPoint[ilev]; }
    amrex::Vector<double> &XVal(int level) { return xypdlXVal[level]; }
    amrex::Vector<double> &YVal(int level) { return xypdlYVal[level]; }

  protected:
    amrex::Vector<list<XYPlotDataListLink *> > dataSets; 
    amrex::Vector<int> xypdlRatios;        // Used for insertion of data sets beginning
                                    // in the middle of a single cell [level]
    amrex::Vector<Real> dX;                // dX in the proper direction [level]
    amrex::Vector<char *> intersectPoint;  // Intersected value (in problem space)
                                    // for labelling plots [level]
    amrex::Vector<double> xypdlLoY, xypdlHiY;  // [level]
    amrex::Vector<amrex::Vector<double> > xypdlXVal, xypdlYVal;  // [level][data]

    amrex::Vector<int> numPoints;
    int  minLevel, maxLevel, curLevel, gridline, whichDir;
    bool updatedQ;
    double startX, endX, offsetX;
    string xypdlDerived;
    XYPlotDataList *copiedFrom;
    amrex::Vector<amrex::BoxList> fabBoxLists, fillBoxLists;
};
}
#endif
// -------------------------------------------------------------------
// -------------------------------------------------------------------
