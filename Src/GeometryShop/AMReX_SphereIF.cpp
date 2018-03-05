#include "AMReX_SphereIF.H"

namespace amrex
{

  SphereIF::
  SphereIF(const Real&     a_radius,
           const RealVect& a_center,
           const bool&     a_inside)
    
  {
    m_radius  = a_radius;
    m_radius2 = m_radius*m_radius;
    m_inside  = a_inside;
    m_center  = a_center;
  }

  Real
  SphereIF::
  value(const RealVect& a_point) const
  {
    RealVect dist = a_point - m_center;
    Real distance2 = dist.radSquared();
    Real retval = distance2 - m_radius2;
    // Change the sign to change inside to outside
    if (!m_inside)
      {
        retval = -retval;
      }

    return retval;
  }

  BaseIF* 
  SphereIF::
  newImplicitFunction() const
  {
    SphereIF* spherePtr = new SphereIF(m_radius,
                                       m_center,
                                       m_inside);

    return static_cast<BaseIF*>(spherePtr);
  }
  
}

