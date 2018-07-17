#include "AMReX_PolynomialIF.H"

namespace amrex
{


  ////
  PolynomialIF::
  PolynomialIF (const Vector<PolyTerm>& a_polynomial,
                const bool&             a_inside)
  {
    m_polynomial = a_polynomial;
    m_inside     = a_inside;
  }


  ////
  Real 
  PolynomialIF::
  value (const RealVect         & a_point,
         const Vector<PolyTerm> & a_polynomial) const
  {
    Real retval;

    int size = a_polynomial.size();

    // Evaluate the polynomial
    retval = 0.0;
    for (int iterm = 0; iterm < size; iterm++)
    {
      Real cur;

      cur = a_polynomial[iterm].coef;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
        cur *= pow(a_point[idir],a_polynomial[iterm].powers[idir]);
      }

      retval += cur;
    }

    // Change the sign to change inside to outside
    if (!m_inside)
    {
      retval = -retval;
    }

    return retval;
  }
  ///
  Real PolynomialIF::value(const RealVect& a_point) const
  {
    return value(a_point,m_polynomial);
  }
  ///
  BaseIF* PolynomialIF::newImplicitFunction() const
  {
    PolynomialIF* polynomialPtr = new PolynomialIF(m_polynomial,
                                                   m_inside);

    return static_cast<BaseIF*>(polynomialPtr);
  }



}
