#ifndef AMREX_EB2_IF_ROTATION_H_
#define AMREX_EB2_IF_ROTATION_H_

#include <AMReX_Array.H>
#include <type_traits>
#include <cmath>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

template <class F>
class RotationIF
{
public:

    RotationIF (F&& a_f, const Real& angle, const int dir)
        : m_f(std::move(a_f)), m_angle(angle), m_dir(dir)
        {}
    RotationIF (F const& a_f, const Real& angle, const int dir)
        : m_f(a_f), m_angle(angle), m_dir(dir)
        {}

    ~RotationIF () {}

    RotationIF (const RotationIF& rhs) = default;
    RotationIF (RotationIF&& rhs) = default;
    RotationIF& operator= (const RotationIF& rhs) = delete;
    RotationIF& operator= (RotationIF&& rhs) = delete;

// Note that angle is measured in radians 
#if (AMREX_SPACEDIM==2)
    Real operator() (const RealArray& p) const
	    {
		Real x =  p[0]*std::cos(m_angle) + p[1]*std::sin(m_angle); 
		Real y = -p[0]*std::sin(m_angle) + p[1]*std::cos(m_angle); 
	        return m_f({x, y});
	    }

#endif

#if (AMREX_SPACEDIM==3)
    Real operator() (const RealArray& p) const
	{
   	if(m_dir==0){
		Real y =  p[1]*std::cos(m_angle) + p[2]*std::sin(m_angle); 
		Real z = -p[1]*std::sin(m_angle) + p[2]*std::cos(m_angle); 
	        return m_f({p[0], y, z});
	 }
	else if(m_dir==1){
		Real x = p[0]*std::cos(m_angle) - p[2]*std::sin(m_angle); 
		Real z = p[0]*std::sin(m_angle) + p[2]*std::cos(m_angle); 
	        return m_f({x, p[1], z});
	}
	else {  // if(m_dir==2){
		Real x =  p[0]*std::cos(m_angle) + p[1]*std::sin(m_angle); 
		Real y = -p[0]*std::sin(m_angle) + p[1]*std::cos(m_angle); 
	        return m_f({x, y, p[2]});
	}
	}
#endif

protected:

    F m_f;
    Real m_angle; 
    int m_dir; 
};

template <class F>
constexpr RotationIF<typename std::decay<F>::type>
rotate (F&&f, const Real& angle, const int dir)
{
    return RotationIF<typename std::decay<F>::type>(std::forward<F>(f),angle, dir);
}

}}

#endif

