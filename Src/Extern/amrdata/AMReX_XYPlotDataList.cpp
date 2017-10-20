// -------------------------------------------------------------------
// XYPlotDataList.cpp 
// -------------------------------------------------------------------
#include <AMReX_XYPlotDataList.H>
#include <AMReX_Utility.H>
#include <cfloat>
#include <limits>

using namespace amrex;


// -------------------------------------------------------------------
XYPlotDataList::XYPlotDataList(const string &derived, int minlevel,
                               int maxlevel, int gridlinein,
			       const Vector<int> &ratiolist,
			       const Vector<Real> &dx,
			       const Vector<char *> &intersectpoint,
			       Real offsetx)
  : dataSets(maxlevel + 1),
    xypdlRatios(ratiolist),
    dX(dx),
    intersectPoint(intersectpoint),
    xypdlLoY(maxlevel + 1),
    xypdlHiY(maxlevel + 1),
    xypdlXVal(maxlevel + 1),
    xypdlYVal(maxlevel + 1),
    numPoints(maxlevel + 1),
    minLevel(minlevel),
    maxLevel(maxlevel),
    gridline(gridlinein),
    offsetX(offsetx),
    xypdlDerived(derived),
    fabBoxLists(maxlevel + 1),
    fillBoxLists(maxlevel + 1)
{
  updatedQ = false;
  curLevel = 0;
  copiedFrom = NULL;
}


// -------------------------------------------------------------------
XYPlotDataList::XYPlotDataList(XYPlotDataList *src)
  : dataSets(src->dataSets),
    xypdlRatios(src->xypdlRatios),
    dX(src->dX),
    intersectPoint(src->intersectPoint),
    xypdlLoY(src->xypdlLoY),
    xypdlHiY(src->xypdlHiY),
    xypdlXVal(src->xypdlXVal),
    xypdlYVal(src->xypdlYVal),
    numPoints(src->numPoints),
    maxLevel(src->maxLevel),
    curLevel(src->curLevel),
    gridline(src->gridline),
    updatedQ(src->updatedQ),
    startX(src->startX),
    endX(src->endX),
    offsetX(src->offsetX),
    xypdlDerived(src->xypdlDerived)
{
  if(src->copiedFrom) {
    copiedFrom = src->copiedFrom;
  } else {
    copiedFrom = src;
  }
}


// -------------------------------------------------------------------
XYPlotDataList::~XYPlotDataList() {
  if(copiedFrom == NULL) {
    for(int ilev(0); ilev <= maxLevel; ++ilev) {
      delete intersectPoint[ilev];
      for(list<XYPlotDataListLink *>::iterator li = dataSets[ilev].begin();
          li != dataSets[ilev].end(); ++li)
      {
        delete (*li);
      }
    }
  }
}


// -------------------------------------------------------------------
void XYPlotDataList::AddFArrayBox(FArrayBox &fab, int whichdir, int level) {
  BL_ASSERT(level >= 0);
  BL_ASSERT(level <= maxLevel);

  whichDir = whichdir;
  if(fabBoxLists[level].isEmpty()) {  // set correct index type
    fabBoxLists[level].convert(fab.box().ixType());
  }
  fabBoxLists[level].push_back(fab.box());
  updatedQ = false;
  int istartx(fab.smallEnd()[whichdir]);
  XYPlotDataListLink *pdll = new XYPlotDataListLink(fab.dataPtr(), istartx,
                                                    fab.length()[whichdir]);
  list<XYPlotDataListLink *>::iterator li = dataSets[level].begin();
  if(li == dataSets[level].end()) {
    dataSets[level].push_back(pdll);
  } else {
    while((*li)->StartXi() < istartx && li != dataSets[level].end()) {
      ++li;
      if(li == dataSets[level].end()) {
        break;
      }
    }
    dataSets[level].insert(li, pdll);
  }
}


// -------------------------------------------------------------------
void XYPlotDataList::UpdateStats() {
  if(updatedQ) {
    return;
  }
  BL_ASSERT(dataSets[minLevel].empty() == false);

  int ilev, startxi, endxi;
  startxi = dataSets[minLevel].front()->StartXi();
  endxi   = dataSets[minLevel].back()->EndXi();
  startX  = offsetX + dX[minLevel] * (double) startxi;
  endX    = offsetX + dX[minLevel] * (double) endxi;

for(int iCurLevel(maxLevel); iCurLevel >= minLevel; --iCurLevel) {

  // back out the probDomain
  Box probDomain = fabBoxLists[minLevel].minimalBox();
  probDomain.refine(amrex::CRRBetweenLevels(minLevel, iCurLevel, xypdlRatios));
  for(int isd(0); isd < BL_SPACEDIM; ++isd) {
    if(isd != whichDir) {  // squish the pd
      if(fabBoxLists[iCurLevel].size() > 0) {
        probDomain.setSmall(isd, (*(fabBoxLists[iCurLevel].begin())).smallEnd(isd));
        probDomain.setBig(isd, (*(fabBoxLists[iCurLevel].begin())).bigEnd(isd));
      }
    }
  }

  Vector<BoxList> unfilledBoxLists(iCurLevel + 1);
  if(unfilledBoxLists[iCurLevel].isEmpty()) {  // convert to correct type
    unfilledBoxLists[iCurLevel].convert(probDomain.ixType());
  }
  unfilledBoxLists[iCurLevel].push_back(probDomain);

  for(ilev = iCurLevel; ilev >= minLevel; --ilev) {
    if(fabBoxLists[ilev].isEmpty()) {  // set correct index type
      fabBoxLists[ilev].convert(probDomain.ixType());
    }
    fillBoxLists[ilev].clear(); 
    fillBoxLists[ilev].convert(unfilledBoxLists[ilev].ixType());
    fillBoxLists[ilev].join(unfilledBoxLists[ilev]);
    fillBoxLists[ilev].intersect(fabBoxLists[ilev]);
    BoxList tempUnfilled(probDomain.ixType());
    for(BoxList::iterator bli = unfilledBoxLists[ilev].begin();
        bli != unfilledBoxLists[ilev].end(); ++bli)
    {
      tempUnfilled.join(amrex::complementIn(*bli, fabBoxLists[ilev]));
    }
    unfilledBoxLists[ilev].clear();
    unfilledBoxLists[ilev].join(tempUnfilled);

    if(ilev > minLevel) {
      unfilledBoxLists[ilev - 1].clear();
      unfilledBoxLists[ilev - 1].convert(probDomain.ixType());
      unfilledBoxLists[ilev - 1].join(unfilledBoxLists[ilev]);
      unfilledBoxLists[ilev - 1].coarsen(amrex::CRRBetweenLevels(ilev - 1,
                                         ilev, xypdlRatios));
      for(int isd(0); isd < BL_SPACEDIM; ++isd) {
        if(isd != whichDir) {  // squish the boxlist
          if(fabBoxLists[ilev - 1].size() > 0) {
            for(BoxList::iterator bli = unfilledBoxLists[ilev - 1].begin();
                bli != unfilledBoxLists[ilev - 1].end(); ++bli)
            {
             (*bli).setSmall(isd, (*(fabBoxLists[ilev - 1].begin())).smallEnd(isd));
             (*bli).setBig(isd, (*(fabBoxLists[ilev - 1].begin())).bigEnd(isd));
	    }
          }
        }
      }
    }
  }

  numPoints[iCurLevel] = 0;
  for(ilev = minLevel; ilev <= iCurLevel; ++ilev) {
    for(BoxList::iterator bli = fillBoxLists[ilev].begin();
        bli != fillBoxLists[ilev].end(); ++bli)
    {
      numPoints[iCurLevel] += (*bli).length(whichDir);
    }
  }
  xypdlXVal[iCurLevel].resize(numPoints[iCurLevel]);
  xypdlYVal[iCurLevel].resize(numPoints[iCurLevel]);

  list<OrderedBoxes> orderedBoxes;
  for(ilev = minLevel; ilev <= iCurLevel; ++ilev) {
    for(BoxList::iterator bli = fillBoxLists[ilev].begin();
        bli != fillBoxLists[ilev].end(); ++bli)
    {
      Box refinedBox(*bli);
      refinedBox.refine(amrex::CRRBetweenLevels(ilev, iCurLevel, xypdlRatios));
      orderedBoxes.push_back(OrderedBoxes(ilev, whichDir, *bli, refinedBox));
    }
  }
  orderedBoxes.sort();
  int xIndex(0);
  for(list<OrderedBoxes>::iterator obli = orderedBoxes.begin();
        obli != orderedBoxes.end(); ++obli)
  {
    for(int i((*obli).DataBox().smallEnd()[whichDir]);
        i <= (*obli).DataBox().bigEnd()[whichDir]; ++i)
    {
      int obLev((*obli).ILevel());
      Real xval((0.5 + i) * dX[obLev] + offsetX);
      xypdlXVal[iCurLevel][xIndex] = xval;
      for(list<XYPlotDataListLink *>::iterator li = (dataSets[obLev]).begin();
          li != (dataSets[obLev]).end(); ++li)
      {
        XYPlotDataListLink *xypd = *li;
	if(i >= xypd->StartXi() && i < xypd->EndXi()) {
	  Real yval(xypd->XYPDLLData()[i - xypd->StartXi()]);
          xypdlYVal[iCurLevel][xIndex] = yval;
	}
      }
      ++xIndex;
    }
  }

  xypdlLoY[iCurLevel] =  std::numeric_limits<Real>::max();
  xypdlHiY[iCurLevel] = -std::numeric_limits<Real>::max();
  for(ilev = minLevel; ilev <= iCurLevel; ++ilev) {
    for(int ii(0); ii < xypdlYVal[ilev].size(); ++ii) {
      xypdlLoY[iCurLevel] = std::min(xypdlLoY[iCurLevel], xypdlYVal[ilev][ii]);
      xypdlHiY[iCurLevel] = std::max(xypdlHiY[iCurLevel], xypdlYVal[ilev][ii]);
    }
  }

}  // end for(iCurLevel...)

  updatedQ = true;
}
// -------------------------------------------------------------------
// -------------------------------------------------------------------
