#ifndef AMREX_EB2_IF_CYLINDER_H_
#define AMREX_EB2_IF_CYLINDER_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

#include <algorithm>

/**
 * \file AMReX_EB2_IF_Cylinder.H
 *
 * Implicit function describing an axis-aligned solid or hollow cylinder.
 */

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex::EB2 {

//! Signed implicit function for cylinders along a coordinate axis.
class CylinderIF
    : GPUable
{
public:
    //! Build an infinite cylinder (extends across the domain in its axis direction).
    CylinderIF (Real a_radius, int a_direction,
                const RealArray& a_center, bool a_inside)
        : CylinderIF(a_radius, -1.0_rt, a_direction, a_center, a_inside)
    { }

    //! Build a finite cylinder with height `a_height`.
    CylinderIF (Real a_radius, Real a_height, int a_direction,
                const RealArray& a_center, bool a_inside)
        : m_radius(a_radius), m_height(a_height), m_direction(a_direction),
          m_center(makeXDim3(a_center)),
          m_sign(a_inside ? 1.0_rt : -1.0_rt)
        {
            AMREX_ASSERT(m_direction < AMREX_SPACEDIM);
        }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
#if (AMREX_SPACEDIM == 3)
        XDim3 pos{.x = x-m_center.x, .y = y-m_center.y, .z = z-m_center.z};
#else
        XDim3 pos{.x = x-m_center.x, .y = y-m_center.y, .z = 0.0_rt};
#endif
        Real d2 = 0.0_rt;
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
            d2 = pos.x*pos.x+pos.y*pos.y;
#if (AMREX_SPACEDIM == 3)
            pdir = pos.z;
#else
            pdir = 0.0_rt;
#endif
            break;
        }
        }

        d2 -= m_radius*m_radius;

        if (m_height < 0.0_rt) {
            return d2*m_sign;
        } else {
            Real rtop = ( pdir - 0.5_rt*m_height);
            Real rbot = (-pdir - 0.5_rt*m_height);
            Real r = amrex::max(d2,rtop,rbot);
            return r*m_sign;
        }
    }

    //! Evaluate the implicit function at a given coordinate.
    [[nodiscard]] AMREX_FORCE_INLINE Real operator() (const RealArray& p) const noexcept
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

}

#endif
