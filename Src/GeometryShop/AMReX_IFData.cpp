#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include <iostream>
#include <iomanip>

#include "NormalDerivative.H"
#include "IFData.H"
//#include "CoordinateSystem.H"
#include "NamespaceHeader.H"

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
    MayDay::Abort("Lo endpoint not in Map");
  }

  if (localInfo.m_cornerSigns.find(hiPt2D) != localInfo.m_cornerSigns.end())
  {
    m_cornerSigns[1] = localInfo.m_cornerSigns[hiPt2D];
  }
  else
  {
    MayDay::Abort("Hi endpoint not in Map");
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

void IFData<1>::print(ostream& a_out) const
{
  string padding = "  ";
  for (int i = 0; i < GLOBALDIM - 1; i++)
  {
    padding += "  ";
  }

  typedef map<int,int> oneDCornerSigns;

  for (oneDCornerSigns::const_iterator it = m_cornerSigns.begin();
       it != m_cornerSigns.end(); ++it)
  {
    a_out << padding << "Vertex "
                     << "("
                     << it->first
                     << ") = "
                     << it->second
                     << "\n";
  }
  a_out << padding << "\n";

  a_out << padding << "m_allVerticesIn  = " << m_allVerticesIn  << "\n";
  a_out << padding << "m_allVerticesOut = " << m_allVerticesOut << "\n";
  a_out << padding << "m_allVerticesOn  = " << m_allVerticesOn  << "\n";
  a_out << padding << "m_badNormal      = " << m_badNormal      << "\n";
  a_out << padding << "\n";

  if (!m_allVerticesOut)
  {
    a_out << padding << "m_globalCoord     = " << m_globalCoord    ;
    a_out << padding << "m_cellCenterCoord = " << m_cellCenterCoord;
    a_out << padding << "m_parentCoord     = " << m_parentCoord    ;
  }
  else
  {
    a_out << padding << "All vertices out" << "\n";
  }
  a_out << padding << "\n";

  int lo = LARGEINTVAL;
  int hi = LARGEINTVAL;
  if (m_cornerSigns.find(0)!= m_cornerSigns.end())
  {
    lo = m_cornerSigns.find(0)->second;
  }
  else
  {
    MayDay::Abort("No lo in cornerSigns");
  }

  if (m_cornerSigns.find(1)!= m_cornerSigns.end())
  {
    hi = m_cornerSigns.find(1)->second;
  }
  else
  {
    MayDay::Abort("No hi in cornerSigns");
  }

  if (lo == OUT && hi == OUT)
  {
    a_out << padding << "Edge Covered" << "\n";
  }
  else if (lo==IN  && hi==IN)
  {
    a_out << padding << "Edge Uncovered" << "\n";
  }
  else if (lo==ON  || hi==ON)
  {
    a_out << padding << "Edge has vertex on the interface" << "\n";
  }
  else
  {
    if (m_intersection == LARGEREALVAL)
    {
      MayDay::Warning("--- No intersection pt");
    }
    else
    {
      std::ios::fmtflags origFlags = a_out.flags();
      int origWidth = a_out.width();
      int origPrecision = a_out.precision();

      a_out << padding << "Intersection pt is "
                       << setw(23)
                       << setprecision(16)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << m_intersection
                       << "\n";

      a_out.flags(origFlags);
      a_out.width(origWidth);
      a_out.precision(origPrecision);
    }
  }
  a_out << padding << "\n";
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

#include "NamespaceFooter.H"
