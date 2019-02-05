#ifndef AMREX_EB2_IF_BOX_H_
#define AMREX_EB2_IF_BOX_H_

#include <AMReX_Array.H>

#include <algorithm>
#include <limits>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

class BoxIF
{
public:

    //! inside: is the fluid inside the box?
    BoxIF (const RealArray& a_lo, const RealArray& a_hi, bool a_inside)
        : m_lo(a_lo), m_hi(a_hi), m_inside(a_inside),
          m_sign( a_inside ? 1.0 : -1.0 )
        {}

    ~BoxIF () {}

    BoxIF (const BoxIF& rhs) noexcept = default;
    BoxIF (BoxIF&& rhs) noexcept = default;
    BoxIF& operator= (const BoxIF& rhs) = delete;
    BoxIF& operator= (BoxIF&& rhs) = delete;

    Real operator() (const RealArray& p) const
    {
        Real r = std::numeric_limits<Real>::lowest();
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            r = std::max(r,   p[i] - m_hi[i]);
            r = std::max(r, -(p[i] - m_lo[i]));
        }
        return r*m_sign;
    }


protected:

    RealArray m_lo;
    RealArray m_hi;
    bool      m_inside;
    //
    Real      m_sign;
};

}}

#endif
