#ifndef AMREX_EB2_IF_POLYNOMIAL_H_
#define AMREX_EB2_IF_POLYNOMIAL_H_

#include <AMReX_Array.H>
#include <AMReX_Vector.H>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {



/*******************************************************************************
 *                                                                             *
 * Represents one term in a general polynomial                                 *
 *                                                                             *
 *******************************************************************************/
class PolyTerm
{
public:

    // Coefficient of this polynomial term
    Real coef;

    // Powers of this polynomial term
    IntVect powers;
};



class PolynomialIF
{
public:

    // inside: is the fluid inside the ellipsoid?
    PolynomialIF (const Vector<PolyTerm> & a_polynomial, bool a_inside = true)
        : m_polynomial(a_polynomial), m_inside(a_inside),
          m_sign( a_inside ? 1.0 : -1.0 )
        {}

    ~PolynomialIF () {}

    PolynomialIF (const PolynomialIF& rhs) = default;
    PolynomialIF (PolynomialIF&& rhs) = default;
    PolynomialIF& operator= (const PolynomialIF& rhs) = delete;
    PolynomialIF& operator= (PolynomialIF&& rhs) = delete;

    Real operator() (const RealArray& p) const {
        Real retval = 0.0;

        // Evaluate the polynomial
        for (int iterm = 0; iterm < m_polynomial.size(); iterm++) {
            Real cur = m_polynomial[iterm].coef;

            for (int idir = 0; idir < SpaceDim; idir++)
                cur *= pow(p[idir], m_polynomial[iterm].powers[idir]);

            retval += cur;
        }

        // Change the sign to change inside to outside
        return m_sign*retval;
    }

protected:
    Vector<PolyTerm> m_polynomial;
    bool             m_inside;
    Real             m_sign;
};

}}

#endif
