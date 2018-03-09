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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "MayDay.H"

#include "CutCellMoments.H"

#include "NamespaceHeader.H"

// Null constructor
CutCellMoments<1>::CutCellMoments()
{
}

// Copy constructor
CutCellMoments<1>::CutCellMoments(const CutCellMoments<1>& a_thisCutCellMoments)

  :m_moments(a_thisCutCellMoments.m_moments),
   m_IFData(a_thisCutCellMoments.m_IFData),
   m_numActiveBounds(0)
{
}

// This constructor is used in the recursion
CutCellMoments<1>::CutCellMoments(const IFData<1> & a_info)
  :m_IFData(a_info),
   m_numActiveBounds(0)
{
  m_badNormal = false;
}

// Destructor
CutCellMoments<1>::~CutCellMoments()
{
}

void CutCellMoments<1>::dump() const
{
  print(pout());
}

Real CutCellMoments<1>::changeMomentCoordinates(OneDMoments           & a_refinedMomentMap,
                                                const IndexTM<int,1>  & a_monomial,
                                                const IndexTM<Real,1> & a_refinedCenterDelta)
{
  Real moment = 0.0;
  int degree = 0;
  degree += a_monomial[0];

  // Loop over the possible values of r in the development of the Taylor series
  for (int r = 0; r <= degree; r++)
  {
    // Generate all the possible monomials of degree r and add up the moment
    IndexTM<int,1> derivative;
    derivative[0] = r;

    // Add up the relevant term of the refined moment map to the map
    Real coeff = 1.0;
    Real factorial = 1.0;

    for (int j = 0; j < derivative[0]; j++)
    {
      coeff *= a_monomial[0] - j;
      factorial *= j+1;
    }

    coeff *= pow(a_refinedCenterDelta[0],a_monomial[0]-derivative[0]);

    moment += coeff * a_refinedMomentMap[derivative] / factorial;
  }

  return moment;
}

void CutCellMoments<1>::changeMomentCoordinatesToCellCenter()
{
  // Move moments from parent coord to cell center coord
  IndexTM<Real,1> delta = m_IFData.m_cellCenterCoord.m_origin;
  delta                -= m_IFData.m_parentCoord    .m_origin;

  OneDMoments copyMoments = m_moments;
  for (OneDMoments::const_iterator it = copyMoments.begin();
       it != copyMoments.end(); ++it)
    {
      IndexTM<int,1> mono = it->first;
      m_moments[mono] = changeMomentCoordinates(copyMoments, mono, delta);
    }
}

void CutCellMoments<1>::changeMomentCoordinatesToParentCenter()
{
  // Move moments from cell center coord to parent coord
  IndexTM<Real,1> delta = m_IFData.m_parentCoord    .m_origin;
  delta                -= m_IFData.m_cellCenterCoord.m_origin;

  OneDMoments copyMoments = m_moments;
  for (OneDMoments::const_iterator it = copyMoments.begin();
       it != copyMoments.end(); ++it)
    {
      IndexTM<int,1> mono = it->first;
      m_moments[mono] = changeMomentCoordinates(copyMoments, mono, delta);
    }
}

void CutCellMoments<1>::initialize(CutCellMoments<1> & a_refinedCutCell)
{
  initializeMap(m_EBmoments,a_refinedCutCell.m_EBmoments);
  initializeMap(m_moments,a_refinedCutCell.m_moments);
}

void CutCellMoments<1>::initializeMap(OneDMoments & a_map1,
                                      OneDMoments & a_map2)
{
  for (OneDMoments::const_iterator it = a_map2.begin(); it != a_map2.end(); ++it)
  {
    a_map1[it->first] = 0.0;
  }
}

Real CutCellMoments<1>::getBdMoment(const IndexTM<int,1>  & a_mono,
                                    const IFData<2>       & a_IFData,
                                    const IndexTM<Real,1> & a_refinedCenterDelta,
                                    OneDMoments             a_fullCellMap)
{
  Real moment = LARGEREALVAL;

  if (a_IFData.m_allVerticesOut)
  {
    moment = 0.0;
  }
  else if (a_IFData.m_allVerticesIn)
  {
    Real loPt = -0.5*m_IFData.m_globalCoord.m_dx[0];
    Real hiPt =  0.5*m_IFData.m_globalCoord.m_dx[0];

    moment = pow(hiPt,a_mono[0] + 1) - pow(loPt,a_mono[0] +1);
    moment = changeMomentCoordinates(a_fullCellMap,a_mono,a_refinedCenterDelta);
  }
  else
  {
    moment = changeMomentCoordinates(m_moments,a_mono,a_refinedCenterDelta);
  }

  return moment;
}

Real CutCellMoments<1>::getBdEBMoment(const IndexTM<int,1>  & a_mono,
                                      const IFData<2>       & a_IFData,
                                      const IndexTM<Real,1> & a_refinedCenterDelta)
{
  Real EBmoment = 0.0;

  return EBmoment;
}

void CutCellMoments<1>::addBdMoments(CutCellMoments<1>     & a_coarseBdCutCell,
                                     const IFData<2>       & a_IFData,
                                     const int             & a_degreePmax,
                                     const bool            & a_useConstraints,
                                     const IndexTM<Real,1> & a_refinedCenterDelta,
                                     const IndexTM<int,1>  & a_localHilo)
{
  OneDMoments fullCellMap;

  if (a_IFData.m_allVerticesIn)
  {
    IndexTM<int,1> degree;

    for (int iDegree = 0; iDegree <= a_degreePmax; ++iDegree)
    {
      degree[0] = iDegree;
      Real loPt = -0.5*m_IFData.m_globalCoord.m_dx[0];
      Real hiPt =  0.5*m_IFData.m_globalCoord.m_dx[0];
      fullCellMap[degree] = pow(hiPt,degree[0] + 1) - pow(loPt,degree[0] +1);
      fullCellMap[degree] /= degree[0]+1;
    }
  }

  for (int iDegree = 0; iDegree <= a_degreePmax; ++iDegree)
  {
    IndexTM<int,1>degree;
    degree[0] = iDegree;
    Real bdMoment;
    if (a_IFData.m_allVerticesIn)
    {
      bdMoment = getBdMoment(degree,a_IFData,a_refinedCenterDelta,fullCellMap);
    }
    else
    {
      bdMoment = getBdMoment(degree,a_IFData,a_refinedCenterDelta);
    }

    a_coarseBdCutCell.m_moments[degree] += bdMoment;
  }
}

Real CutCellMoments<1>::getMoment(const IndexTM<int,1> & a_mono,
                                  const EBorVol        & a_EBorVOL) const
{
  return getMoment(a_mono);
}

// Method for reading moments
Real CutCellMoments<1>::getMoment(const IndexTM<int,1> & a_mono) const
{
  Real moment = LARGEREALVAL;

  // Find the vol in the map
  OneDMoments::const_iterator it = m_moments.find(a_mono);

  if (it != m_moments.end())
  {
    moment = it->second;
  }
  else
  {
    MayDay::Abort("No volume moment in m_moments");
  }

  return moment;
}

// Read geometric data and sanity check
Real CutCellMoments<1>::getVol(const EBorVol & a_EBorVol) const
{
  IndexTM<int,1> zero;
  zero[0] = 0;
  Real volume = getMoment(zero);

  return volume;
}

IndexTM<Real,1> CutCellMoments<1>::getCentroid(const EBorVol& a_EBorVol) const
{
  IndexTM<Real,1> centroid;
  IndexTM<int,1> first;
  first[0] = 1;

  Real volume = getVol(a_EBorVol);

  if (volume <= 0)
  {
    MayDay::Abort("Volume <= 0");
  }
  else
  {
    centroid[0] = getMoment(first)/volume;
  }

  return centroid;
}

bool CutCellMoments<1>::isCovered() const
{
  return m_IFData.m_allVerticesOut;
}

bool CutCellMoments<1>::isRegular() const
{
  return m_IFData.m_allVerticesIn;
}

void CutCellMoments<1>::print(ostream & a_out) const
{
  string padding = "  ";
  for (int i = 0; i < GLOBALDIM - 1; i++)
  {
    padding += "  ";
  }

  a_out << padding << "Dim = 1" << "\n";
  a_out << padding << "\n";

  for (OneDMoments::const_iterator it = m_moments.begin(); it != m_moments.end();++it)
  {
    std::ios::fmtflags origFlags = a_out.flags();
    int origWidth = a_out.width();
    int origPrecision = a_out.precision();

    a_out << padding << "Integral "
                     << it->first
                     << " = "
                     << setw(23)
                     << setprecision(16)
                     << setiosflags(ios::showpoint)
                     << setiosflags(ios::scientific)
                     << it->second
                     << "\n";

    a_out.flags(origFlags);
    a_out.width(origWidth);
    a_out.precision(origPrecision);
  }
  a_out << padding << "\n";

  a_out << padding << "IFData:" << "\n";
  a_out << m_IFData;
}

void CutCellMoments<1>::operator=(const CutCellMoments<1> & a_cutCellMoments)
{
  // Only copy if the two objects are distinct
  if (this != &a_cutCellMoments)
  {
    m_IFData  = a_cutCellMoments.m_IFData;
    m_moments = a_cutCellMoments.m_moments;
  }
}

#include "NamespaceFooter.H"
