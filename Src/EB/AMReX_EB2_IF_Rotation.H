#ifndef AMREX_EB2_IF_ROTATION_H_
#define AMREX_EB2_IF_ROTATION_H_

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_Array.H>
#include <type_traits>
#include <cmath>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

template <class F>
class RotationIF
{
public:

    RotationIF (F&& a_f, Real angle, int dir)
        : m_f(std::move(a_f)),
          m_cos_angle(std::cos(angle)),
          m_sin_angle(std::sin(angle)),
          m_dir(dir)
        {}
    RotationIF (F const& a_f, Real angle, int dir)
        : m_f(a_f),
          m_cos_angle(std::cos(angle)),
          m_sin_angle(std::sin(angle)),
          m_dir(dir)
        {}

    RotationIF (const RotationIF& rhs) = default;
    RotationIF (RotationIF&& rhs) = default;
    RotationIF& operator= (const RotationIF& rhs) = delete;
    RotationIF& operator= (RotationIF&& rhs) = delete;

// Note that angle is measured in radians 
#if (AMREX_SPACEDIM==2)
    inline Real operator() (const RealArray& p) const noexcept
    {
        Real x =  p[0]*m_cos_angle + p[1]*m_sin_angle; 
        Real y = -p[0]*m_sin_angle + p[1]*m_cos_angle; 
        return m_f({x, y});
    }

    template <class U=F, class = typename std::enable_if<IsGPUable<U>::value>::type>
    AMREX_GPU_HOST_DEVICE inline
    Real operator() (Real x, Real y) const noexcept
    {
        return m_f( x*m_cos_angle + y*m_sin_angle,
                   -x*m_sin_angle + y*m_cos_angle);   
    }
#endif

#if (AMREX_SPACEDIM==3)
    inline Real operator() (const RealArray& p) const noexcept
    {
        switch(m_dir) {
        case(0):
   	{
            Real y =  p[1]*m_cos_angle + p[2]*m_sin_angle; 
            Real z = -p[1]*m_sin_angle + p[2]*m_cos_angle; 
            return m_f({p[0], y, z});
        }
        case(1):
        {
            Real x = p[0]*m_cos_angle - p[2]*m_sin_angle; 
            Real z = p[0]*m_sin_angle + p[2]*m_cos_angle; 
            return m_f({x, p[1], z});
	}
        default:
	{
            Real x =  p[0]*m_cos_angle + p[1]*m_sin_angle; 
            Real y = -p[0]*m_sin_angle + p[1]*m_cos_angle; 
            return m_f({x, y, p[2]});
	}
        }
    }

    template <class U=F, class = typename std::enable_if<IsGPUable<U>::value>::type>
    AMREX_GPU_HOST_DEVICE inline
    Real operator() (Real x, Real y, Real z) const noexcept
    {
        switch(m_dir) {
        case(0):
   	{
            return m_f(x,
                        y*m_cos_angle + z*m_sin_angle,
                       -y*m_sin_angle + z*m_cos_angle); 
        }
        case(1):
        {
            return m_f(x*m_cos_angle - z*m_sin_angle,
                       y,
                       x*m_sin_angle + z*m_cos_angle);
	}
        default:
	{
            return m_f( x*m_cos_angle + y*m_sin_angle,
                       -x*m_sin_angle + y*m_cos_angle,
                        z);
	}
        }
    }
#endif

protected:

    F m_f;
    Real m_cos_angle;
    Real m_sin_angle;
    int m_dir; 
};

template <class F>
struct IsGPUable<RotationIF<F>, typename std::enable_if<IsGPUable<F>::value>::type>
    : std::true_type {};

template <class F>
constexpr RotationIF<typename std::decay<F>::type>
rotate (F&&f, const Real angle, const int dir)
{
    return RotationIF<typename std::decay<F>::type>(std::forward<F>(f),angle, dir);
}

}}

#endif

