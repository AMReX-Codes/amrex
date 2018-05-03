#include "AMReX_EllipsoidIF.H"

namespace amrex
{

  EllipsoidIF::EllipsoidIF(const RealVect& a_radii,
                           const RealVect& a_center,
                           const bool&     a_inside)
  {
    // Remember the parameters
    m_radii  = a_radii;
    m_center = a_center;
    m_inside = a_inside;

    // Precompute the radii squared
    m_radii2  = m_radii;
    m_radii2 *= m_radii;
  }

  EllipsoidIF::EllipsoidIF(const EllipsoidIF& a_inputIF)
  {
    // Remember the parameters
    m_radii  = a_inputIF.m_radii;
    m_center = a_inputIF.m_center;
    m_inside = a_inputIF.m_inside;

    // Precompute the radii squared
    m_radii2  = m_radii;
    m_radii2 *= m_radii;
  }

  EllipsoidIF::~EllipsoidIF()
  {
  }

  Real EllipsoidIF::value(const RealVect& a_point) const
  {
    Real retval;

    // Compute the equation of the ellipsoid
    Real sum;

    sum = 0.0;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real cur;
      cur = a_point[idir] - m_center[idir];

      sum += cur*cur / m_radii2[idir];
    }

    // The sum should be 1.0 on the surface of the ellipsoid
    retval = sum - 1.0;

    // Change the sign to change inside to outside
    if (!m_inside)
    {
      retval = -retval;
    }

    return retval;
  }

  Real 
  EllipsoidIF::
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
      retval = 2.0*(a_point[maxDir] - m_center[maxDir])/m_radii2[maxDir];
    }
    else if (derivativeOrder == 2)
    {
      if (a_deriv[maxDir] == 2)
      {
        // unmixed second partial = 2.0
        retval = 2.0/m_radii2[maxDir];
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

  BaseIF* EllipsoidIF::newImplicitFunction() const
  {
    EllipsoidIF* ellipsoidPtr = new EllipsoidIF(m_radii,
                                                m_center,
                                                m_inside);

    return static_cast<BaseIF*>(ellipsoidPtr);
  }
  
}

