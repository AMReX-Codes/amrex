#ifndef AMREX_EB2_IF_POLYNOMIAL_H_
#define AMREX_EB2_IF_POLYNOMIAL_H_
#include <AMReX_Config.H>

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <memory>
#include <cmath>

/**
 * \file AMReX_EB2_IF_Polynomial.H
 *
 * Support for implicit functions defined by multivariate polynomials.
 */

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex::EB2 {

/********************************************************************************
 *                                                                              *
 * Represents one term in a general polynomial                                  *
 *                                                                              *
 ********************************************************************************/
//! Represents one monomial in a polynomial implicit function.
class PolyTerm
{
public:

    //! Coefficient of this polynomial term
    Real coef;

    //! Powers of this polynomial term
    IntVect powers;
};

template <unsigned int N>
//! GPU-ready polynomial implicit function with a fixed number of terms.
class PolyIF
    : public GPUable
{
public:

    //! Store the polynomial terms and orientation (inside/outside).
    PolyIF (const GpuArray<PolyTerm,N> & a_polynomial, bool a_inside = true)
        : m_polynomial(a_polynomial),
          m_sign( a_inside ? 1.0_rt : -1.0_rt )
    {}

    //! Evaluate the polynomial at explicit coordinates.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        Real retval = 0.0_rt;
        for (auto const& term : m_polynomial) {
            retval += term.coef * AMREX_D_TERM(  std::pow(x, term.powers[0]),
                                               * std::pow(y, term.powers[1]),
                                               * std::pow(z, term.powers[2]));
        }
        return m_sign*retval;
    }

    //! Array-based helper overload.
    AMREX_FORCE_INLINE Real operator() (const RealArray& p) const noexcept {
        return this->operator()(AMREX_D_DECL(p[0],p[1],p[2]));
    }

protected:
    GpuArray<PolyTerm,N> m_polynomial;
    Real                 m_sign;
};


//! Host-only dynamically-sized polynomial implicit function.
class PolynomialIF  // No GPU support
{
public:

    //! Store the term vector and orientation (inside/outside).
    PolynomialIF (const Vector<PolyTerm> & a_polynomial, bool a_inside = true)
        : m_polynomial(a_polynomial),
          m_inside(a_inside),
          m_sign( a_inside ? 1.0_rt : -1.0_rt ),
          m_size(int(m_polynomial.size()))
    {}

    ~PolynomialIF () = default;
    PolynomialIF (const PolynomialIF& rhs) = default;
    PolynomialIF (PolynomialIF&& rhs) = default;
    PolynomialIF& operator= (const PolynomialIF& rhs) = delete;
    PolynomialIF& operator= (PolynomialIF&& rhs) = delete;

    [[nodiscard]] AMREX_FORCE_INLINE
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        Real retval = 0.0_rt;
        for (auto const& term : m_polynomial) {
            retval += term.coef * AMREX_D_TERM(  std::pow(x, term.powers[0]),
                                               * std::pow(y, term.powers[1]),
                                               * std::pow(z, term.powers[2]));
        }
        return m_sign*retval;
    }

    //! Array-based helper overload.
    [[nodiscard]] AMREX_FORCE_INLINE Real operator() (const RealArray& p) const noexcept {
        return this->operator()(AMREX_D_DECL(p[0],p[1],p[2]));
    }

protected:
    Vector<PolyTerm> m_polynomial;
    bool             m_inside;
    Real             m_sign;
    int              m_size;
};

}

#endif
