#ifndef AMREX_EB2_IF_CYLINDER_H_
#define AMREX_EB2_IF_CYLINDER_H_

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

#include <algorithm>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

class CylinderIF
    : GPUable
{
public:
    // inside: is the fluid inside the cylinder?

    //! Overloaded constructor for infinite cylinder IF: if the user does not
    //! provide a cylinder height, the IF corresponds to a cylinder that spans
    //! the whole domain.
    CylinderIF (Real a_radius, int a_direction,
                const RealArray& a_center, bool a_inside)
        : CylinderIF(a_radius, -1.0, a_direction, a_center, a_inside)
    { }

    //! Overloaded constructor for finite cylinder IF: if the user does specify a
    //! cylinder height (length), then the IF corresponds to a cylinder extending
    //! only to height/2 from it's center.
    CylinderIF (Real a_radius, Real a_height, int a_direction,
                const RealArray& a_center, bool a_inside)
        : m_radius(a_radius), m_height(a_height), m_direction(a_direction),
          m_center(makeXDim3(a_center)),
          m_sign(a_inside ? 1.0 : -1.0)
        {
            AMREX_ASSERT(m_direction < AMREX_SPACEDIM);
        }

    CylinderIF (const CylinderIF& rhs) noexcept = default;
    CylinderIF (CylinderIF&& rhs) noexcept = default;
    CylinderIF& operator= (const CylinderIF& rhs) = delete;
    CylinderIF& operator= (CylinderIF&& rhs) = delete;

    AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        XDim3 pos{AMREX_D_DECL(x-m_center.x,
                               y-m_center.y,
                               z-m_center.z)};
        Real d2 = 0.0;
        Real pdir;
        switch (m_direction) {
        case 0 :
        {
#if (AMREX_SPACEDIM == 3)
            d2 = pos.y*pos.y+pos.z*pos.z;
#elif (AMREX_SPACEDIM == 2)
            d2 = pos.y*pos.y;
#endif
            pdir = pos.x;
            break;
        }
        case 1:
        {
#if (AMREX_SPACEDIM == 3)
            d2 = pos.x*pos.x+pos.z*pos.z;
#elif (AMREX_SPACEDIM == 2)
            d2 = pos.x*pos.x;
#endif
            pdir = pos.y;
            break;
        }
        default:
        {
#if (AMREX_SPACEDIM == 3)
            d2 = pos.x*pos.x+pos.y*pos.y;
#endif
            pdir = pos.z;
            break;
        }
        }

        d2 -= m_radius*m_radius;

        if (m_height < 0.0) {
            return d2*m_sign;
        } else {
            Real rtop = ( pdir - 0.5*m_height);
            Real rbot = (-pdir - 0.5*m_height);
            Real r = amrex::max(d2,rtop,rbot);
            return r*m_sign;
        }
    }

    inline Real operator() (const RealArray& p) const noexcept
    {
        return this->operator() (AMREX_D_DECL(p[0], p[1], p[2]));
    }

protected:

    Real      m_radius;
    Real      m_height;
    int       m_direction;
    XDim3     m_center;
    //
    Real      m_sign;
};

}}

#endif
