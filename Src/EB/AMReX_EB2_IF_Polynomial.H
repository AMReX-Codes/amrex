#ifndef AMREX_EB2_IF_POLYNOMIAL_H_
#define AMREX_EB2_IF_POLYNOMIAL_H_

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {


/********************************************************************************
 *                                                                              *
 * Integer-based powers:                                                        *
 *     -- Do not use std::pow => we don't need non-integer powers               *
 *                                                                              *
 ********************************************************************************/

inline Real exp_int (Real base, int exponent)
{
    if (exponent < 0)   return 1./exp_int(base, -exponent);
    if (exponent == 0 ) return 1.;
    if (exponent == 1 ) return base;

    Real tmp = exp_int(base, exponent/2);

    if (exponent%2 == 0) return tmp*tmp;
    else return base*tmp*tmp;
}



/********************************************************************************
 *                                                                              *
 * Represents one term in a general polynomial                                  *
 *                                                                              *
 ********************************************************************************/
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
                cur *= exp_int(p[idir], m_polynomial[iterm].powers[idir]);

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
