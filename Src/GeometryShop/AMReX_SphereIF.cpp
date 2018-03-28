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
  
    ///return the partial derivative at the point
  Real 
  SphereIF::
  derivative(const  IntVect& a_deriv,
             const RealVect& a_point) const
  {
    //negative derivs make no sense
    BL_ASSERT(a_deriv.min() >= 0);
    Real retval = 0;
    int maxDir = a_deriv.maxDir(false);
    int derivativeOrder = a_deriv.sum();

    if (derivativeOrder == 0)
    {
      retval = value(a_point);
    }
    else if (derivativeOrder == 1)
    {
      retval = 2.0*(a_point[maxDir] - m_center[maxDir]);
    }
    else if (derivativeOrder == 2)
    {
      if (a_deriv[maxDir] == 2)
      {
        // unmixed second partial = 2.0
        retval = 2.0;
      }
      else
      {
        // mixed partials = 0.0
        retval = 0.0;
      }
    }
    else
    {
      // higher partials = 0.0
      retval = 0.0;
    }

    // Change the sign to change inside to outside
    if (!m_inside && derivativeOrder > 0)
    {
      retval = -retval;
    }

    return retval;
  }
}
