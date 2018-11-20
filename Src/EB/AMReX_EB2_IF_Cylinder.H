#ifndef AMREX_EB2_IF_CYLINDER_H_
#define AMREX_EB2_IF_CYLINDER_H_

#include <AMReX_Array.H>

#include <algorithm>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

class CylinderIF
{
public:
    // inside: is the fluid inside the cylinder?

    // Overloaded constructor for infinite cylinder IF: if the user does not
    // provide a cylinder height, the IF corresponds to a cylinder that spans
    // the whole domain.
    CylinderIF (Real a_radius, int a_direction,
                const RealArray& a_center, bool a_inside)
        : CylinderIF(a_radius, -1.0, a_direction, a_center, a_inside)
    { }

    // Overloaded constructor for finite cylinder IF: if the user does specify a
    // cylinder height (length), then the IF corresponds to a cylinder extending
    // only to height/2 from it's center.
    CylinderIF (Real a_radius, Real a_height, int a_direction,
                const RealArray& a_center, bool a_inside)
        : m_radius(a_radius), m_height(a_height), m_direction(a_direction),
          m_center(a_center), m_inside(a_inside),
          m_radius2(a_radius*a_radius),
          m_halfheight(0.5*a_height),
          m_sign(a_inside ? 1.0 : -1.0)
        {
            AMREX_ASSERT(m_direction < AMREX_SPACEDIM);
        }

    ~CylinderIF () {}

    CylinderIF (const CylinderIF& rhs) noexcept = default;
    CylinderIF (CylinderIF&& rhs) noexcept = default;
    CylinderIF& operator= (const CylinderIF& rhs) = delete;
    CylinderIF& operator= (CylinderIF&& rhs) = delete;

    Real operator() (const RealArray& p) const
    {
        RealArray pos{AMREX_D_DECL(p[0]-m_center[0],
                                   p[1]-m_center[1],
                                   p[2]-m_center[2])};
        Real d2 = 0.0;
        switch (m_direction) {
        case 0 :
        {
#if (AMREX_SPACEDIM == 3)
            d2 = pos[1]*pos[1]+pos[2]*pos[2];
#elif (AMREX_SPACEDIM == 2)
            d2 = pos[1]*pos[1];
#endif
            break;
        }
        case 1:
        {
#if (AMREX_SPACEDIM == 3)
            d2 = pos[0]*pos[0]+pos[2]*pos[2];
#elif (AMREX_SPACEDIM == 2)
            d2 = pos[0]*pos[0];
#endif
            break;
        }
        default:
        {
#if (AMREX_SPACEDIM == 3)
            d2 = pos[0]*pos[0]+pos[1]*pos[1];
#endif
            break;
        }
        }

        d2 -= m_radius2;

        if (m_height < 0.0) {
            return d2*m_sign;
        } else {
            Real rtop = ( pos[m_direction] - m_halfheight);
            Real rbot = (-pos[m_direction] - m_halfheight);
            Real r = std::max(d2,std::max(rtop,rbot));
            return r*m_sign;
        }
    }


protected:

    Real      m_radius;
    Real      m_height;
    int       m_direction;
    RealArray m_center;
    bool      m_inside;
    //
    Real      m_radius2;
    Real      m_halfheight;
    Real      m_sign;
};

}}

#endif
