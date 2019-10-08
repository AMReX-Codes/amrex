#ifndef AMREX_EB2_IF_POLYNOMIAL_H_
#define AMREX_EB2_IF_POLYNOMIAL_H_

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <memory>
#include <cmath>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

/********************************************************************************
 *                                                                              *
 * Represents one term in a general polynomial                                  *
 *                                                                              *
 ********************************************************************************/
class PolyTerm
{
public:

    //! Coefficient of this polynomial term
    Real coef;

    //! Powers of this polynomial term
    IntVect powers;
};



class PolynomialIF
    : public GPUable
{
public:

    //! inside: is the fluid inside the ellipsoid?
    PolynomialIF (const Vector<PolyTerm> & a_polynomial, bool a_inside = true)
        : m_polynomial(a_polynomial),
          m_inside(a_inside),
          m_sign( a_inside ? 1.0 : -1.0 ),
          m_size(m_polynomial.size()),
          m_dp(static_cast<PolyTerm*>(The_Device_Arena()->alloc(sizeof(PolyTerm)*m_size))),
          m_sp(m_dp,Gpu::Deleter(The_Device_Arena()))
    {
        Gpu::htod_memcpy(m_dp, m_polynomial.data(), m_size*sizeof(PolyTerm));
    }

    ~PolynomialIF () {}

    PolynomialIF (const PolynomialIF& rhs) = default;
    PolynomialIF (PolynomialIF&& rhs) = default;
    PolynomialIF& operator= (const PolynomialIF& rhs) = delete;
    PolynomialIF& operator= (PolynomialIF&& rhs) = delete;

    AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
#ifdef AMREX_DEVICE_COMPILE
        PolyTerm const* AMREX_RESTRICT pt = m_dp;
#else
        PolyTerm const* AMREX_RESTRICT pt = m_polynomial.data();
#endif
        Real retval = 0.0;
        for (int iterm = 0; iterm < m_size; ++iterm) {
            retval += pt[iterm].coef * AMREX_D_TERM(  std::pow(x, pt[iterm].powers[0]),
                                                    * std::pow(y, pt[iterm].powers[1]),
                                                    * std::pow(z, pt[iterm].powers[2]));
        }
        return m_sign*retval;
    }

    inline Real operator() (const RealArray& p) const noexcept {
        return this->operator()(AMREX_D_DECL(p[0],p[1],p[2]));
    }

protected:
    Vector<PolyTerm> m_polynomial;
    bool             m_inside;
    Real             m_sign;
    int              m_size;
    PolyTerm*        m_dp;
    std::shared_ptr<PolyTerm> m_sp;
};

}}

#endif
