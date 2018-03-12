


#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>


#include "AMReX_NoRefinement.H"
#include "AMReX_LSProblem.H"
#include "AMReX_MinimalCCCM.H"

namespace amrex
{

// Null constructor
  MinimalCCCM<1>::MinimalCCCM()
  {
  }

// Copy constructor
  MinimalCCCM<1>::MinimalCCCM(const MinimalCCCM<1>& a_thisMinimalCCCM)

    :m_cutCellMoments(a_thisMinimalCCCM.m_cutCellMoments)
  {
  }

// This constructor is used in the recursion
  MinimalCCCM<1>::MinimalCCCM(const IFData<1> & a_info)
    :m_cutCellMoments(a_info)
  {
  }

// Destructor
  MinimalCCCM<1>::~MinimalCCCM()
  {
  }

// Integrate along line segments aligned in a coordinate direction
  void MinimalCCCM<1>::computeMoments(const int              & a_orderPmax,
                                      const int              & a_degreePmax)
  {
    int lo = 0;
    int hi = 1;
    int loSign = m_cutCellMoments.m_IFData.m_cornerSigns[lo];
    int hiSign = m_cutCellMoments.m_IFData.m_cornerSigns[hi];

    // If entire edge out of the fluid, then moments = 0.0
    if (loSign <= ON && hiSign<= ON)
    {
      for (int iDegree = 0; iDegree <= a_degreePmax; ++iDegree)
      {
        // Definition of m_cutCellMoments.m_moments typedef requires that we define degree as
        // a oneTuple:
        IndexTM<int,1> degree;
        degree[0] = iDegree;
        m_cutCellMoments.m_moments[degree] = 0.0;
      }
    }
    else
    {
      // Assign loPt and hiPt in m_cutCellMoments.m_IFData.m_parentCoord system

      Real loPt = LARGEREALVAL;
      Real hiPt = LARGEREALVAL;

      // m_origin is an IndexTM<Real,1>, which implies we need to use [0] everywhere
      // m_intersection is undefined if hiSign >= ON  && loSign >= ON
      if (loSign >= ON)
      {
        loPt = m_cutCellMoments.m_IFData.m_parentCoord.convertDir(
          -0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0],
          m_cutCellMoments.m_IFData.m_cellCenterCoord,
          0);
      }
      else
      {
        BL_ASSERT(-0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0]
                  <= m_cutCellMoments.m_IFData.m_intersection
                  &&
                  m_cutCellMoments.m_IFData.m_intersection
                  <= 0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0]);

        loPt = m_cutCellMoments.m_IFData.m_parentCoord.convertDir(
          m_cutCellMoments.m_IFData.m_intersection,
          m_cutCellMoments.m_IFData.m_cellCenterCoord,
          0);
      }

      if (hiSign >= ON)
      {
        hiPt = m_cutCellMoments.m_IFData.m_parentCoord.convertDir(
          0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0],
          m_cutCellMoments.m_IFData.m_cellCenterCoord,
          0);
      }
      else
      {
        BL_ASSERT(-0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0]
                  <= m_cutCellMoments.m_IFData.m_intersection
                  &&
                  m_cutCellMoments.m_IFData.m_intersection
                  <= 0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[0]);

        hiPt = m_cutCellMoments.m_IFData.m_parentCoord.convertDir(
          m_cutCellMoments.m_IFData.m_intersection,
          m_cutCellMoments.m_IFData.m_cellCenterCoord,
          0);
      }

      // Integrate x^degree over the line segment[loPt,hiPt]
      computeMomentsUsingBinomial(loPt,hiPt,loSign,hiSign,a_degreePmax);
    }
  }

  void MinimalCCCM<1>::simpleComputeMoments(const Real & a_loPt,
                                            const Real & a_hiPt,
                                            const int  & a_degreePmax)
  {
    for (int iDegree = 0; iDegree <= a_degreePmax; ++iDegree)
    {
      //definition of m_cutCellMoments.m_moments typedef requires that we define degree thus:
      IndexTM<int,1>degree;
      degree[0] = iDegree;
      m_cutCellMoments.m_moments[degree] = POW(a_hiPt,iDegree + 1) - POW(a_loPt,iDegree +1);
      //    Real dxFactor = POW(m_cutCellMoments.m_IFData.m_globalCoord.m_dx[0],iDegree + 1);
      //m_cutCellMoments.m_moments[degree] *= dxFactor;
      m_cutCellMoments.m_moments[degree] /= (iDegree + 1);
    }
  }

  void MinimalCCCM<1>::computeMomentsUsingBinomial(const Real & a_loPt,
                                                   const Real & a_hiPt,
                                                   const int  & a_loSign,
                                                   const int  & a_hiSign,
                                                   const int  & a_degreePmax)
  {
    for (int iDegree = 0; iDegree <= a_degreePmax; ++iDegree)
    {
      Real epsilon = a_hiPt - a_loPt;

      // Definition of m_cutCellMoments.m_moments typedef requires that we define degree thus:
      IndexTM<int,1> degree;
      degree[0] = iDegree;
      m_cutCellMoments.m_moments[degree] = 0.0;

      // New method without substracting higher order terms
      for (int j = 1; j <= iDegree + 1; j++)
      {
        int bigger = j;
        int smaller = j;
        if (iDegree + 1 - j > j)
        {
          bigger = iDegree + 1 - j;
        }
        else
        {
          smaller = iDegree + 1 - j;
        }

        int numerator = 1;
        for (int i = bigger + 1; i <= iDegree + 1; ++i)
        {
          numerator *= i;
        }

        int denominator = 1;
        for (int i = 1; i <= smaller; ++i)
        {
          denominator *= i;
        }

        Real factor = numerator / denominator;
        if (a_loSign >= ON)
        {
          m_cutCellMoments.m_moments[degree] += factor * POW(a_loPt,iDegree + 1 - j) * POW(epsilon,j);
        }
        else if (a_hiSign >= ON)
        {
          m_cutCellMoments.m_moments[degree] -= factor * POW(a_hiPt,iDegree + 1 - j) * POW(epsilon,j) * POW(-1.0,j);
        }
      }

      //Real dxFactor = POW(m_cutCellMoments.m_IFData.m_globalCoord.m_dx[0],iDegree + 1);
      //m_cutCellMoments.m_moments[degree] *= dxFactor;
      m_cutCellMoments.m_moments[degree] /= iDegree + 1;
    }
  }

  void MinimalCCCM<1>::print(ostream & a_out) const
  {
    m_cutCellMoments.print(a_out);
  }

  void MinimalCCCM<1>::operator=(const MinimalCCCM<1> & a_MinimalCCCM)
  {
    // Only copy of the two objects are distinct
    if (this != &a_MinimalCCCM)
    {
      m_cutCellMoments = a_MinimalCCCM.m_cutCellMoments;
    }
  }
}

