#include <iostream>
#include <iomanip>

#include "AMReX_NormalDerivative.H"
#include "AMReX_IFData.H"

//leave to default faulse.   Moving the coords works better
//but makes for weird convergence tests
bool LocalCoordMoveSwitch::s_turnOffMoveLocalCoords = true;

// empty constructor (dim == 1)
IFData<1>::IFData()
{
}

IFData<1>::IFData(const IFData<1>& a_IFData)
  :m_cornerSigns    (a_IFData.m_cornerSigns),
   m_intersection   (a_IFData.m_intersection),
   m_globalCoord    (a_IFData.m_globalCoord),
   m_cellCenterCoord(a_IFData.m_cellCenterCoord),
   m_parentCoord    (a_IFData.m_parentCoord),
   m_allVerticesIn  (a_IFData.m_allVerticesIn),
   m_allVerticesOut (a_IFData.m_allVerticesOut),
   m_allVerticesOn  (a_IFData.m_allVerticesOn),
   m_badNormal      (a_IFData.m_badNormal)
{
}

// Constructor from the implicit function
IFData<1>::IFData(const IFData<2> & a_2DIFData,
                  const int       & a_maxOrder,
                  const int       & a_idir,
                  const int       & a_hilo)
  :m_globalCoord(a_2DIFData.m_globalCoord,a_idir),
   m_cellCenterCoord(a_2DIFData.m_cellCenterCoord,a_idir),
   m_parentCoord(a_2DIFData.m_localCoord,a_idir)
{
  // we want the edge on the a_hilo side of the square with normal in the
  // a_idir direction
  IFData<2>localInfo = a_2DIFData;

  // This 2D edgeIndex locates the 1D edge in the edgeIntersection map
  IndexTM<int,2>twoDEdge;
  twoDEdge[0] = (a_idir + 1)%2;
  twoDEdge[1] = a_hilo;

  m_intersection = LARGEREALVAL;
  if (localInfo.m_intersections.find(twoDEdge) != localInfo.m_intersections.end())
  {
    m_intersection = localInfo.m_intersections[twoDEdge];
  }

  // This 2D vertex locates the hi and lo ends of the 1D segment in the
  // cornerSigns map
  IndexTM<int,2>loPt2D;
  loPt2D[(a_idir + 1)%2] = 0;
  loPt2D[a_idir] = a_hilo;

  IndexTM<int,2>hiPt2D;
  hiPt2D[(a_idir+ 1)%2] = 1;
  hiPt2D[a_idir] = a_hilo;

  if (localInfo.m_cornerSigns.find(loPt2D) != localInfo.m_cornerSigns.end())
  {
    m_cornerSigns[0] = localInfo.m_cornerSigns[loPt2D];
  }
  else
  {
    amrex::Abort("Lo endpoint not in Map");
  }

  if (localInfo.m_cornerSigns.find(hiPt2D) != localInfo.m_cornerSigns.end())
  {
    m_cornerSigns[1] = localInfo.m_cornerSigns[hiPt2D];
  }
  else
  {
    amrex::Abort("Hi endpoint not in Map");
  }

  // set bools
  m_allVerticesIn  = true;
  m_allVerticesOut = true;
  m_allVerticesOn  = true;

  if (m_cornerSigns[0] != ON || m_cornerSigns[1] != ON)
  {
    m_allVerticesOn  = false;
  }

  if (m_cornerSigns[0] == IN || m_cornerSigns[1] == IN)
  {
    m_allVerticesOut = false;
  }

  if (m_cornerSigns[0] == OUT || m_cornerSigns[1] == OUT)
  {
    m_allVerticesIn = false;
  }

  //there is no normal in one dimension. However, if m_badNormal = true at a lower dimension, then the higher dimension refines.
  m_badNormal = false;
}

// Destructor (dim == 1)
IFData<1>::~IFData()
{
}


// equals operator
void IFData<1>::operator=(const IFData & a_IFData)
{
  if (this != &a_IFData)
  {
    m_cornerSigns     = a_IFData.m_cornerSigns;
    m_intersection    = a_IFData.m_intersection;
    m_globalCoord     = a_IFData.m_globalCoord;
    m_cellCenterCoord = a_IFData.m_cellCenterCoord;
    m_parentCoord     = a_IFData.m_parentCoord;
    m_allVerticesIn   = a_IFData.m_allVerticesIn;
    m_allVerticesOut  = a_IFData.m_allVerticesOut;
    m_allVerticesOn   = a_IFData.m_allVerticesOn;
    m_badNormal   = a_IFData.m_badNormal;
  }
}

