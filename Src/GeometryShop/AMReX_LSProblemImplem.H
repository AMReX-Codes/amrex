#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LSPROBLEMIMPLEM_H_
#define _LSPROBLEMIMPLEM_H_

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
#include <string>
using std::endl;

#include "AMReX_ConstrainedLS.H"
#include "AMReX_LSquares.H"
#include "AMReX_MultiIndex.H"
#include "AMReX_IFData.H"
#include "AMReX_CutCellMoments.H"

namespace amrex
{

  template<int dim> int LSProblem<dim>::nChooseR(int a_n,
                                                 int a_r)
  {
    if (a_r == 0) return 1;
    int num = 1;
    int den = 1;
    for (int i = 1; i <= a_r; ++i)
    {
      num *= (a_n+1-i);
      den *= i;
    }
    return num/den;
  }

  template<int dim> void LSProblem<dim>::monoMaxMin(Real                   & a_maxVal,
                                                    Real                   & a_minVal,
                                                    const IndexTM<int,dim> & a_mono,
                                                    const IFData<dim>      & a_IFData)
  {
    a_maxVal = 0.0;
    a_minVal = 0.0;

    for (typename IFData<dim>::CornerSigns::const_iterator cornerIt = a_IFData.m_cornerSigns.begin();
         cornerIt != a_IFData.m_cornerSigns.end();
         ++cornerIt)
    {
      const typename IFData<dim>::Vertex & vertex = cornerIt->first;

      // represent the vertex as an RvDim in cell centered coordinates
      typename IFData<dim>::RvDim corner;
      for (int idir = 0; idir < dim; ++idir)
      {
        corner[idir] = vertex[idir] - 0.5;
        corner[idir] *= a_IFData.m_cellCenterCoord.m_dx[idir];
      }

      // compute coordinates of corner in local coordinates
      typename IFData<dim>::RvDim cornerCoord = a_IFData.m_localCoord.convert(corner,a_IFData.m_cellCenterCoord);

      // monomial value at the corner
      Real monoCorner = 1.0;
      for (int idir = 0; idir < dim ; ++idir)
      {
        monoCorner *= POW(cornerCoord[idir],a_mono[idir]);
      }

      if (monoCorner > a_maxVal)
      {
        a_maxVal = monoCorner;
      }

      if (monoCorner < a_minVal)
      {
        a_minVal = monoCorner;
      }
    }
  }

  template<int dim> void LSProblem<dim>::computeBounds(const IndexTM<Real,dim>   & a_dx,
                                                       const CutCellMoments<dim> & a_ccm)
  {
    m_lowerBound.assign(-HUGE, m_lowerBound.size());
    m_upperBound.assign( HUGE, m_upperBound.size());

    int nEBBounds = 0;
    if (m_degreeP == 0 && dim == 2) // || m_degreeP == 1)
    {
      nEBBounds = m_numP;       // A lower bound for each EB. Haven't figured out a useful u.b. yet
    }

    // int nVar=m_numP + m_numPLess1; // total variables in problem

    // This constraint is a lower bound on the zero monomial of the EB.
    // It is based on inverting integral(mono)<Max(mono)*integral(1)
    if (m_degreeP == 0 && dim == 2)
    {
      Real lobnd = 0.0; // minimum, but we can do better
      IvDim whichMono; // for tracking
      for (typename CutCellMoments<dim>::PthMoment::const_iterator it = a_ccm.m_EBmoments.begin();
           it != a_ccm.m_EBmoments.end();++it)
      {
        IvDim mono = it->first;
        Real maxMono = 0.0;
        Real minMono = 0.0;
        monoMaxMin(maxMono,minMono,mono,a_ccm.m_IFData);
        Real val = it->second;
        if (minMono < 0.0)
        {
          Real newLoBnd = val/minMono;
          if (newLoBnd > lobnd)
          {
            lobnd = newLoBnd;
            whichMono = mono;
          }
        }
        if (maxMono > 0.0)
        {
          Real newLoBnd = val/maxMono;
          if (newLoBnd > lobnd)
          {
            lobnd = newLoBnd;
            whichMono = mono;
          }
        }
      }

      // The constraints matrix fits CiT x + ci >=0
      // so the rows are variables and the columns are constraints (may be a FORTRAN holdover)
      IvDim mono = IvDim::Zero;

      int ndx = 0; // todo: guarantee there should be only one for (0,0)
      int variableNdx = ndx;
      m_lowerBound[variableNdx]=lobnd;
      m_upperBound[variableNdx]=HUGE;

    }

    // Now create constraints for the volume integrals.
    for (typename PthMomentLoc::const_iterator it = m_monoLocPLess1.begin();
         it != m_monoLocPLess1.end();++it)
    {
      IvDim mono = it->first;
      // The constraints matrix fits CiT x + ci >=0
      // so the rows are variables and the columns are constraints (may be a FORTRAN holdover)
      int ndx = it->second;
      int variableNdx = ndx + m_numP;
      Real lobnd, hibnd;

      // Compute the lower and upper bounds
      // based on quadrature of all the regions that contribute
      // negatively or positively.
      momentBounds(lobnd,hibnd,mono, a_ccm.m_IFData);



      bool isZeroDegree = (mono.sum() == 0);

      if (isZeroDegree )
      {
        Real zeroDegreeLoBnd = 0.0;

        // Refine the lower bound of the zero degree monomial
        // based by inverting the identity
        // integral(mono)/integral(zeromono) < Max(mono)
        // Note that we already have integral(mono) from previous steps
        IvDim whichMono;
        // Real valAtWhichMono=LARGEREALVAL;
        for (typename CutCellMoments<dim>::PthMoment::const_iterator it2 = a_ccm.m_moments.begin();
             it2 != a_ccm.m_moments.end();++it2)
        {
          IvDim mono2 = it2->first;
          Real val = it2->second;
          Real maxMono = 0.0;
          Real minMono = 0.0;
          monoMaxMin(maxMono,minMono,mono2,a_ccm.m_IFData);
          if (minMono < 0.0)
          {
            Real newLoBnd = val/minMono;
            if (newLoBnd > zeroDegreeLoBnd)
            {
              zeroDegreeLoBnd = newLoBnd;
              whichMono = mono2;
            }
          }
          if (maxMono > 0.0)
          {
            Real newLoBnd = val/maxMono;
            if (newLoBnd > zeroDegreeLoBnd)
            {
              zeroDegreeLoBnd = newLoBnd;
              whichMono = mono2;
            }
          }
        }
        lobnd = std::max(lobnd, zeroDegreeLoBnd);
      }


      m_lowerBound[variableNdx] = lobnd;
      m_upperBound[variableNdx] = hibnd;

    }


  }

// Calculate upper and lower bounds for a moment
  template<int dim> void LSProblem<dim>::momentBounds(Real              & a_lobnd,
                                                      Real              & a_hibnd,
                                                      const IvDim       & a_mono,
                                                      const IFData<dim> & a_IFData)
  {
    a_lobnd = 0.0;
    a_hibnd = 0.0;
    // Enumerate though every vertex, and integrate the contributions
    // from the (local) cell center to the vertex in local coordinates
    typename IFData<dim>::CornerSigns::const_iterator cornerIt;
    // pout() << a_mono << endl;
    for (cornerIt =  a_IFData.m_cornerSigns.begin();
         cornerIt != a_IFData.m_cornerSigns.end();
         ++cornerIt)
    {
      const typename IFData<dim>::Vertex & vertex = cornerIt->first;
      // represent the vertex as an RvDim in cell centered coordinates
      typename IFData<dim>::RvDim corner;
      for (int idir = 0; idir < dim; ++idir)
      {
        corner[idir] = vertex[idir] - 0.5;
        corner[idir] *= a_IFData.m_cellCenterCoord.m_dx[idir];
      }

      // compute coordinates of corner in local coordinates
      typename IFData<dim>::RvDim cornerCoord
        = a_IFData.m_localCoord.convert(corner,a_IFData.m_cellCenterCoord);

      Real partialCellMag = 1.0;
      for (int idir = 0; idir < dim ; ++idir)
      {
        Real pPlus1 =(Real)( a_mono[idir]+1);
        partialCellMag *= POW(cornerCoord[idir],pPlus1)/pPlus1;
        if (cornerCoord[idir] < 0.0) partialCellMag *=-1.0;
      }
      if (partialCellMag > 0.0) a_hibnd += partialCellMag;
      if (partialCellMag < 0.0) a_lobnd += partialCellMag;
    }
  }

  template<int dim> LSProblem<dim>::~LSProblem()
  {
    if (m_matrix != NULL)
    {
      // free m_matrix
      int numRows = dim*m_numP;
      int numCols = m_numP + m_numPLess1;
      freeArray(numRows,numCols,m_matrix);

    }
  }

// this constructor is used when just a list of monomials is wanted
  template<int dim> LSProblem<dim>::LSProblem(const int  & a_degreeP,
                                              const bool & a_useConstraints)
    :m_degreeP(a_degreeP),
     m_numActiveBounds(0),
     m_useConstraints(a_useConstraints)
  {
    fillMap(m_monoLocP,m_locMonoP,m_degreeP);
    fillMap(m_monoLocPLess1,m_locMonoPLess1,m_degreeP-1);
    m_numP      = numMonomials(m_degreeP);
    m_numPLess1 = numMonomials(m_degreeP-1);
    m_matrix = NULL;
  }

// constructor for solving a LS problem
  template<int dim> LSProblem<dim>::LSProblem(const int               & a_order,
                                              const int               & a_degreeP,
                                              const bool              & a_useConstraints,
                                              const IndexTM<Real,dim> & a_normal)
    :m_order(a_order),
     m_degreeP(a_degreeP),
     m_numActiveBounds(0),
     m_useConstraints(a_useConstraints),
     m_normal(a_normal)
  {
    fillMap(m_monoLocP,m_locMonoP,m_degreeP);
    fillMap(m_monoLocPLess1,m_locMonoPLess1,m_degreeP-1);
    m_numP      = numMonomials(m_degreeP);
    m_numPLess1 = numMonomials(m_degreeP-1);
    setMatrix();

    m_unknowns.resize(m_numP + m_numPLess1);
    m_rhs.resize(m_numP*dim);

    if (a_useConstraints)
    {
      m_lowerBound.resize(m_numP+m_numPLess1);
      m_upperBound.resize(m_numP+m_numPLess1);
    }
  }

//
  template<int dim> int LSProblem<dim>::invertNormalEq(const Vector<Real> & a_rhs,
                                                       Vector<Real>       & a_residual)
  {
    m_rhs = a_rhs;
    int retCode = -1;
    if (!m_useConstraints)
    {
      LSquares lsSolver;
      lsSolver.LeastSquares(m_matrix,m_unknowns,m_rhs);
      m_numActiveBounds = 0;
      retCode = 0;
    }
    else
    {
      ConstrainedLS cls;
      ConstrainedLS::LSResult result = cls.solveBoundConstrained(m_unknowns,
                                                                 m_matrix,
                                                                 m_rhs,
                                                                 m_lowerBound,
                                                                 m_upperBound);
      m_numActiveBounds = cls.numberActiveConstraints();
      retCode = (result == ConstrainedLS::SUCCESS) ? 0 : -1;
      switch(result)
      {
      case(ConstrainedLS::SUCCESS):
        break;
      case(ConstrainedLS::SINGULAR):
        pout() << "Singular LS problem. Linearly dependent columns (zero normal?)" << endl;
        break;
      case(ConstrainedLS::INCONSISTENT_BOUNDS):
        pout() << "Conflicting upper/lower bounds in LS problem. Zero or bad normal?" << endl;
        break;
      case(ConstrainedLS::UNDERDETERMINED):
        amrex::Error("Malformed LS problem. Underdetermined");
      case(ConstrainedLS::UNCONVERGED):
        pout() << "BVLS failed to converge properly" << endl;
        break;
      }

    }

    a_residual.resize(3);
    a_residual[0] = 0.0;
    a_residual[1] = 0.0;
    a_residual[2] = 0.0;
    Real maxRi = 0.0;
    for (int i = 0 ; i < m_numP*dim ; i++)
    {
      Real AtimeX = 0.0;
      for (int j = 0 ; j < m_numP + m_numPLess1 ; j++)
      {
        AtimeX += m_matrix[i][j] * m_unknowns[j];
      }
      Real ri = std::abs(AtimeX - m_rhs[i]);
      if (ri > maxRi)
      {
        a_residual[0] = ri;
        maxRi = ri;
      }
      a_residual[1] += ri;
      a_residual[2] += ri * ri;
    }
    a_residual[2] = sqrt(a_residual[2]);
    // pout() << "Residual : " << a_residual[1] << ":" << a_residual[2] << endl;
    return retCode;
  }

  template<int dim> int LSProblem<dim>::factorial(const int & a_n,
                                                  const int & a_m)
  {
    int retval = 1;
    if (a_n < 0 || a_m < 0)
    {
      amrex::Abort("Attempting n! for n < 0");
    }
    for (int i = a_m+1; i <= a_n; ++i)
    {
      retval *= i;
    }
    return retval;
  }

  template<int dim> int LSProblem<dim>::numMonomials(const int & a_monoDegree)
  {
    int retval = LARGEINTVAL;
    if (a_monoDegree == -1)
    {
      retval =  0;
    }
    else
    {
      int bigger;
      int smaller;
      if (dim- 1 > a_monoDegree)
      {
        bigger = dim-1;
        smaller = a_monoDegree;
      }
      else
      {
        smaller = dim-1;
        bigger = a_monoDegree;
      }
      int numerator = factorial(dim - 1 + a_monoDegree,bigger);
      int denominator = factorial(smaller);
      // pout() << "numerator = "   << numerator   << endl;
      // pout() << "denominator = " << denominator << endl;
      retval =  numerator/denominator;
    }
    return retval;
  }

// uses dimension & degree create matrix for overdetermined system
  template<int dim> void LSProblem<dim>::setMatrix()
  {
    // solving m_matrix[x] = b

    // initializes m_matrix to zeros

    int numRows = dim*m_numP;
    int numCols = m_numP + m_numPLess1;
    allocArray(numRows,numCols,m_matrix);

    // iterate through the list of mono of Degree P
    // pout() << " Number of mono of DegreeP = " << m_monoLocP.size() << endl;
    for (typename PthMomentLoc::const_iterator it = m_monoLocP.begin();
         it != m_monoLocP.end();++it)
    {
      for (int idir = 0; idir< dim; ++idir)
      {
        // this entry corresponds to the integral over the boundary
        int row = dim*(it->second) + idir;
        int pcol = it->second;
        m_matrix[row][pcol] = -m_normal[idir];

        // differentiate the mono and
        IvDim Dmono;
        IvDim mono = it->first;
        int coeff = LARGEINTVAL;
        // diff(it->first) = coeff*Dmono
        differentiate(coeff,Dmono,idir,mono);
        int pLess1Col= LARGEINTVAL;
        if (coeff == LARGEINTVAL)
        {
          amrex::Abort("problem wth differentiate");
        }
        if (coeff > 0)
        {
          // find which mono this is in the list
          if (m_monoLocPLess1.find(Dmono) != m_monoLocPLess1.end())
          {
            pLess1Col = m_monoLocPLess1[Dmono] + m_numP;
            m_matrix[row][pLess1Col] = coeff;
          }
          else
          {
            amrex::Abort("can't find derived mono");
          }
        }
      }
    }
  }

// differentiate a_mono w.r.t. x_idir. Answer = a_coeff*a_Dmono
  template<int dim> void LSProblem<dim>::differentiate(int         & a_coeff,
                                                       IvDim       & a_Dmono,
                                                       int         & a_idir,
                                                       const IvDim & a_mono)
  {
    if (a_mono[a_idir] > 0)
    {
      a_coeff          = a_mono[a_idir];
      a_Dmono          = a_mono;
      a_Dmono[a_idir] -= 1;
    }
    else if (a_mono[a_idir] == 0)
    {
      a_coeff = 0;
      for (int idir = 0; idir < dim; ++idir)
      {
        a_Dmono[idir] = LARGEINTVAL;
      }
    }
    else
    {
      amrex::Abort("Monomial has negative power");
    }
  }

  template<int dim> void LSProblem<dim>::setRhs(const Vector<Real> & a_rhs)
  {
  }

  template<int dim> void LSProblem<dim>::fillMap(PthMomentLoc & a_monoLoc,
                                                 LocPthMoment & a_locMono,
                                                 const int    & a_degree)
  {
    if (a_degree >= 0)
    {
      Vector<IvDim> monomials;

      generateMultiIndices(monomials,a_degree);

      for (int i = 0; i < monomials.size(); i++)
      {
        const IvDim & monomial = monomials[i];

        a_monoLoc[monomial] = i;
        a_locMono[i] = monomial;
      }
    }
  }

  template<int dim> void LSProblem<dim>::allocArray(const int & a_rows,
                                                    const int & a_cols,
                                                    Real**    & a_A)
  {
    a_A = new Real* [a_rows];

    for (int i = 0; i < a_rows; i++)
    {
      a_A[i] = new Real [a_cols];

      Real* scanA = a_A[i];
      for (int j = 0; j < a_cols; j++)
      {
        *(scanA++) = 0.0;
      }
    }
  }

  template<int dim> void LSProblem<dim>::freeArray(const int & a_rows,
                                                   const int & a_cols,
                                                   Real**    & a_A)
  {
    for (int i = 0; i < a_rows; i++)
    {
      delete[] a_A[i];
    }

    delete[] a_A;
  }

  template<int dim > void LSProblem<dim>::print(ostream & a_out) const
  {
    a_out << "Dim = " << dim << ", degree = " << m_degreeP << '\n';
    a_out << "m_monoLocP has " << m_monoLocP.size() << " elements" << '\n';

    for (typename PthMomentLoc::const_iterator it = m_monoLocP.begin();
         it != m_monoLocP.end();++it)
    {
      a_out << "Monomial =  " << it->first << ", Loc = " << it->second << '\n';
    }

    a_out << "Dim = " << dim << '\n';
    a_out << "m_locMonoP has " << m_locMonoP.size() << " elements" << '\n';

    for (typename LocPthMoment::const_iterator it = m_locMonoP.begin();
         it != m_locMonoP.end();++it)
    {
      a_out << "Loc =  " << it->first << ", Monomial = " << it->second << '\n';
    }
    // one degree lower
    a_out << "m_locMonoPLess1 has " << m_locMonoPLess1.size() << " elements" << '\n';
    for (typename PthMomentLoc::const_iterator it = m_monoLocPLess1.begin();
         it != m_monoLocPLess1.end();++it)
    {
      a_out << "Monomial =  " << it->first << ", Loc = " << it->second << '\n';
    }

    a_out << "m_locMonoPLess1 has " << m_locMonoPLess1.size() << " elements" << '\n';

    for (typename LocPthMoment::const_iterator it = m_locMonoPLess1.begin();
         it != m_locMonoPLess1.end();++it)
    {
      a_out << "Loc =  " << it->first << ", Monomial = " << it->second << '\n';
    }

    a_out << "Matrix and rhs for least squares problem of dim = " << dim << endl;
    outputMatrix();
    outputRhs();
    outputUnknowns();
    outputBounds();
  }

  template<int dim> ostream& operator<<(ostream        & a_out,
                                        LSProblem<dim> & a_lsProblem)
  {
    a_lsProblem.print(a_out);
    return a_out;
  }

  template<int dim> void LSProblem<dim>::outputMatrix() const
  {
    int rows = m_locMonoP.size()*dim;
    int cols = m_locMonoPLess1.size() + m_locMonoP.size();
    pout() << "numRows = " << rows << ", numCols = " << cols << endl;
    // pout() << "outputting " << name << endl;
    for (int i = 0; i < rows; i++)
    {
      for (int j = 0; j < cols; j++)
      {
        pout() << m_matrix[i][j] << " ";
      }
      pout() << endl;
    }
  }

  template<int dim> void LSProblem<dim>::outputRhs() const
  {
    pout() << "Outputting Rhs" << endl;
    for (int i = 0; i < m_rhs.size(); i++)
    {
      pout() << "m_rhs[" << i << "] = " << m_rhs[i] << endl;
    }
  }

  template<int dim> void LSProblem<dim>::outputUnknowns()const
  {
    pout() << "Outputting Unknowns" << endl;
    for (int i = 0; i < m_unknowns.size(); i++)
    {
      pout() << "m_unknowns[" << i << "] = " << m_unknowns[i] << endl;
    }
  }

  template<int dim> void LSProblem<dim>::outputBounds()const
  {
    if (m_useConstraints)
    {
      pout() << "Outputting Lower/Upper bounds" << endl;
      for (int i = 0; i < m_unknowns.size(); i++)
      {
        pout() << "m_lowerBound["  <<  i << "] = " << m_lowerBound[i]
               << " m_upperBound[" <<  i << "] = " << m_upperBound[i] << endl;
      }
    }
    else
    {
      pout() << "Problem is unconstrained" << endl;
    }
  }
}


#endif
