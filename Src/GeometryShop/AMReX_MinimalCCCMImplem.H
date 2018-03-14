#ifndef _MINIMALCCCMIMPLEM_H_
#define _MINIMALCCCMIMPLEM_H_

#include "AMReX_NoRefinement.H"
#include "AMReX_LSProblem.H"
#include "AMReX_Factorial.H"


// Null constructor
template <int dim> MinimalCCCM<dim>::MinimalCCCM()
{
}

// Copy constructor
template <int dim> MinimalCCCM<dim>::MinimalCCCM(const MinimalCCCM<dim>& a_MinimalCCCM)
  :m_cutCellMoments(a_MinimalCCCM.m_cutCellMoments)
{
}

// Use this for initializing
template <int dim> MinimalCCCM<dim>::MinimalCCCM(const IFData<dim>& a_info)
  :m_cutCellMoments(a_info)
{
  m_cutCellMoments.m_bdCCOn = false;

  for (int hilo = 0; hilo < 2; ++hilo)
  {
    for (int idir = 0; idir < dim; ++idir)
    {
      // Identifier for which boundary cutCellMoment
      Iv2 bdId;
      bdId[BDID_DIR] = idir;
      bdId[BDID_HILO] = hilo;

      IFData<dim-1> reducedInfo(a_info,a_info.m_maxOrder+1,idir,hilo);

      CutCellMoments<dim-1>bdCutCellMoments(reducedInfo);

      m_cutCellMoments.m_bdCutCellMoments[bdId] = bdCutCellMoments;


      // Notice whether at least one lower dimensional cutCell is on the
      // interface
      if (reducedInfo.m_allVerticesOn)
      {
        m_cutCellMoments.m_bdCCOn = true;
      }
    }
  }
}

// Destructor
template <int dim> MinimalCCCM<dim>::~MinimalCCCM()
{
}


// Solve
template <int dim> void MinimalCCCM<dim>::computeMoments(const int                & a_orderPmax,
                                                         const int                & a_degreePmax)
{
  int zerothOrderOfAccuracy = 0;
  int highestDegree         = a_orderPmax + a_degreePmax;

  Vector<Real> RNorm(3);
  for (int i = 0; i < 3; i++)
  {
    RNorm[i] = LARGEREALVAL;
  }

  for (int i = 0; i <= highestDegree; i++)
  {
    m_cutCellMoments.m_residual.push_back(RNorm);
  }

  m_boundaryMomentsComputed = false;

  computeMomentsRecursively(zerothOrderOfAccuracy,
                            highestDegree);

  IndexTM<int,dim> refineInDir;

  // Move moments from local coord to parent coord
  IndexTM<Real,dim> delta = m_cutCellMoments.m_IFData.m_parentCoord.m_origin;
  delta                  -= m_cutCellMoments.m_IFData.m_localCoord .m_origin;

  PthMoment copyMoments = m_cutCellMoments.m_moments;
  for (typename PthMoment::const_iterator it = copyMoments.begin();
       it != copyMoments.end(); ++it)
  {
    IvDim mono = it->first;
    m_cutCellMoments.m_moments[mono] = m_cutCellMoments.changeMomentCoordinates(copyMoments, mono, delta);
  }

  PthMoment copyEBMoments = m_cutCellMoments.m_EBmoments;
  for (typename PthMoment::const_iterator it = copyEBMoments.begin(); it != copyEBMoments.end(); ++it)
  {
    IvDim mono = it->first;
    m_cutCellMoments.m_EBmoments[mono] = m_cutCellMoments.changeMomentCoordinates(copyEBMoments, mono, delta);
  }

  // Move bdCutCell moments from parent coord to cell center coord. From here on bdCutCell moments will
  // not be used in any least squares problem

  for (typename BdCutCellMoments::iterator it = m_cutCellMoments.m_bdCutCellMoments.begin(); it != m_cutCellMoments.m_bdCutCellMoments.end(); ++it)
  {
    it->second.changeMomentCoordinatesToCellCenter();
  }
}

template <int dim> void MinimalCCCM<dim>::computeMomentsRecursively(const int                & a_orderPmax,
                                                                    const int                & a_degreePmax)
{
  //pout() << "a_orderPmax  = " << a_orderPmax << endl;
  //pout() << "a_degreePmax = " << a_degreePmax << endl;
  //pout() << "m_cutCellMoments.m_IFData.m_maxOrder " << m_cutCellMoments.m_IFData.m_maxOrder  << endl;
  BL_ASSERT(m_cutCellMoments.m_IFData.m_maxOrder >= a_orderPmax);

  if (m_cutCellMoments.m_IFData.m_allVerticesOut)
  {
    for (int iOrder = 0; iOrder <= a_orderPmax; ++iOrder)
    {
      LSProblem<dim> lsp(iOrder + a_degreePmax, false);

      // Fill moments
      const PthMomentLoc monDegreeP = lsp.getMonomialLocMapDegreeP();
      for (typename PthMomentLoc::const_iterator it = monDegreeP.begin(); it != monDegreeP.end(); ++it)
      {
        m_cutCellMoments.m_EBmoments[it->first] = 0.0;
      }

      const PthMomentLoc monDegreePLess1 = lsp.getMonomialLocMapDegreePLess1();
      for (typename PthMomentLoc::const_iterator it = monDegreePLess1.begin(); it != monDegreePLess1.end(); ++it)
      {
        m_cutCellMoments.m_moments[it->first] = 0.0;
      }
    }
  }
  else if (m_cutCellMoments.m_IFData.m_allVerticesIn && !m_cutCellMoments.m_bdCCOn)
  {
    for (int iOrder = 0; iOrder <= a_orderPmax; ++iOrder)
    {
      LSProblem<dim> lsp(iOrder + a_degreePmax, false);

      // Fill moments of degree P and P-1
      const PthMomentLoc monDegreeP = lsp.getMonomialLocMapDegreeP();
      for (typename PthMomentLoc::const_iterator it = monDegreeP.begin(); it != monDegreeP.end(); ++it)
      {
        m_cutCellMoments.m_EBmoments[it->first] = 0.0;
      }

      const PthMomentLoc monDegreePLess1 = lsp.getMonomialLocMapDegreePLess1();
      for (typename PthMomentLoc::const_iterator it = monDegreePLess1.begin(); it != monDegreePLess1.end(); ++it)
      {
        m_cutCellMoments.m_moments[it->first] = m_cutCellMoments.fullCellQuadrature(it->first,m_cutCellMoments.m_IFData.m_parentCoord);
      }
    }
  }
  else
  {
    // Only compute the boundary moments if they haven't already been
    // computed (earlier in the recursion)
    if (!m_boundaryMomentsComputed)
    {
      for (typename BdCutCellMoments::iterator it = m_cutCellMoments.m_bdCutCellMoments.begin(); it != m_cutCellMoments.m_bdCutCellMoments.end(); ++it)
      {
        MinimalCCCM<dim-1> subProblem(it->second.m_IFData);
        subProblem.computeMoments(a_orderPmax,a_degreePmax+1);
        it->second = subProblem.m_cutCellMoments;
      }

      m_boundaryMomentsComputed = true;
    }

    // Make a LSProb
    IvDim zeroDerivative = IndexTM<int,dim>::Zero;
    LSProblem<dim> lsp(a_orderPmax,a_degreePmax, false,m_cutCellMoments.m_IFData.m_normalDerivatives[zeroDerivative]);

    Vector<Real> rhs = computeRhs(lsp,a_orderPmax);

    // Solve the problem and return residual
    int lsCode=lsp.invertNormalEq(rhs,m_cutCellMoments.m_residual[a_degreePmax]);
    if (lsCode != 0)
    {
//      amrex::pout() << "Geometry generation least squares problem failed with residual:"
//             << m_cutCellMoments.m_residual[a_degreePmax]<< endl;
//      lsp.print(pout());
//      amrex::pout () << "Problem occurred generating geometry for these cut cell moments: " << m_cutCellMoments << endl;
      amrex::Error("Geometry generation error.[2]");
    }


    // Fill moments
    const PthMomentLoc monDegreePLess1 = lsp.getMonomialLocMapDegreePLess1();
    for (typename PthMomentLoc::const_iterator it = monDegreePLess1.begin();
         it != monDegreePLess1.end(); ++it)
    {
      m_cutCellMoments.m_moments[it->first] = lsp.getUnknown(it->second + lsp.getNumberDegP());
    }

    const PthMomentLoc monDegreeP = lsp.getMonomialLocMapDegreeP();
    for (typename PthMomentLoc::const_iterator it = monDegreeP.begin(); it != monDegreeP.end(); ++it)
    {
      m_cutCellMoments.m_EBmoments[it->first] = lsp.getUnknown(it->second);
    }

  }

  if (a_degreePmax > 0)
  {
    computeMomentsRecursively(a_orderPmax + 1,
                              a_degreePmax - 1);
  }
}

template <int dim> Vector<Real> MinimalCCCM<dim>::computeRhs(LSProblem<dim> & a_lsp,
                                                             const int      & a_orderPmax)
{
  // Resize rhs
  int numEq = dim*a_lsp.getNumberDegP();
  Vector<Real> rhs(numEq);

  // For each moment iterate thru bd CutCellMoments incrementing same comp
  // of rhs
  const LocPthMoment& locMap = a_lsp.getLocMonomialMapDegreeP();
  for (typename LocPthMoment::const_iterator it = locMap.begin(); it != locMap.end(); ++it)
  {
    int jth = it->first;
    IvDim mono = it->second;

    int hiSide = 1;
    int loSide = 0;
    Iv2 bdId;

    for (int idir = 0; idir < dim; ++idir)
    {
      // Which lower dimensional monomial corresponds (mono,j)
      IndexTM<int,dim-1> mono1Less;
      for (int jdir = 0; jdir < dim; ++jdir)
      {
        if (jdir < idir)
        {
          mono1Less[jdir] = mono[jdir];
        }
        else if (jdir > idir)
        {
          mono1Less[jdir-1] = mono[jdir];
        }
      }

      bdId[0] = idir;
      bdId[1] = hiSide;

      Real hiMom = m_cutCellMoments.m_bdCutCellMoments[bdId].m_moments[mono1Less];

      bdId[1] = loSide;

      Real loMom = m_cutCellMoments.m_bdCutCellMoments[bdId].m_moments[mono1Less];
      int exponent = it->second[idir];

      Real loSideValue;
      Real hiSideValue;

      loSideValue = m_cutCellMoments.m_IFData.m_localCoord.convertDir(-0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[idir],
                                                                      m_cutCellMoments.m_IFData.m_cellCenterCoord,
                                                                      idir);

      hiSideValue = m_cutCellMoments.m_IFData.m_localCoord.convertDir( 0.5*m_cutCellMoments.m_IFData.m_cellCenterCoord.m_dx[idir],
                                                                       m_cutCellMoments.m_IFData.m_cellCenterCoord,
                                                                       idir);

      Real loFactor = std::pow(loSideValue,exponent);
      Real hiFactor = std::pow(hiSideValue,exponent);

      rhs[(dim*jth) + idir] = hiMom*hiFactor - loMom*loFactor;

      // Add the Taylor series terms
      for (int order = 1; order <= a_orderPmax; order++)
      {
        Vector<IvDim> taylorMonomials;

        generateMultiIndices(taylorMonomials,order);

        for (int i = 0; i < taylorMonomials.size(); i++)
        {
          const IvDim & taylorMonomial = taylorMonomials[i];

          IvDim totalMonomial = mono + taylorMonomial;

          if (m_cutCellMoments.m_EBmoments.find(totalMonomial) !=
              m_cutCellMoments.m_EBmoments.end())
          {
            Real normalDerivative = m_cutCellMoments.m_IFData.m_normalDerivatives[taylorMonomial][idir];
            Real fact = factorial(taylorMonomial);

            Real moment = m_cutCellMoments.m_EBmoments[totalMonomial];

            rhs[(dim*jth) + idir] += normalDerivative * moment / fact;
          }
          else
          {
            amrex::Error("Unable to find needed monomial for Taylor series");
          }
        }
      }
    }
  }

  return rhs;
}

// Operators
template <int dim> void MinimalCCCM<dim>::operator=(const MinimalCCCM<dim> & a_MinimalCCCM)
{
  // Only copy if the objects are distinct
  if (this != &a_MinimalCCCM)
  {
    m_cutCellMoments = a_MinimalCCCM.m_cutCellMoments;
  }
}


template <int dim> Real MinimalCCCM<dim>::factorial(const IvDim & a_multiIndex) const
{
  Real fact = 1.0;

  for (int i = 0; i < dim; i++)
  {
    for (int j = 2; j <= a_multiIndex[i]; j++)
    {
      fact *= j;
    }
  }

  return fact;
}


#endif
