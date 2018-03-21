#include "AMReX_PlaneIF.H"

namespace amrex
{
  PlaneIF::
  PlaneIF(const RealVect& a_normal,
          const RealVect& a_point,
          const bool&     a_inside)
      : m_normal(a_normal),
        m_point (a_point),
        m_inside(a_inside)
  {}

  Real 
  PlaneIF::
  value(const RealVect & a_point) const
  {
    Real retval = 0;
    for (int idir = 0 ; idir < SpaceDim; idir++)
      {
        retval += (a_point[idir] - m_point[idir]) * m_normal[idir];
      }
    if(m_inside)
      {
        retval = -retval;
      }
    return retval;
  }

  BaseIF* 
  PlaneIF::
  newImplicitFunction() const
  {
    PlaneIF* newPtr = new PlaneIF(m_normal,
                                     m_point,
                                     m_inside);


    return static_cast<BaseIF*>(newPtr);
  }

    ///return the partial derivative at the point
  Real 
  PlaneIF::
  derivative(const  IntVect& a_deriv,
             const RealVect& a_point) const
  {
    //negative derivs make no sense
    BL_ASSERT(a_deriv.min() >= 0);
    int order = a_deriv.sum();

    Real retval = 0;

    if (order == 0)
    {
      retval = value(a_point);
    }
    else if (order == 1)
    {
      bool doAbs = true;
      retval = m_normal[a_deriv.maxDir(doAbs)];

      if (m_inside)
      {
        retval = -retval;
      }
    }
    else
    {
      retval = 0.0;
    }
    if (m_inside && order > 0)
    {
      retval = -retval;
    }

    return retval;
  }
}
