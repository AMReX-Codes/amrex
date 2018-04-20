#include "AMReX_LatheIF.H"
#include "AMReX_UnionIF.H"

namespace amrex
{
  LatheIF::LatheIF(const BaseIF& a_impFunc1,
                   const bool&   a_inside)
  {
    // Make a copy of the implicit function and note only one was given
    m_impFunc1 = a_impFunc1.newImplicitFunction();

    // Save inside flag
    m_inside = a_inside;
  }


  LatheIF::~LatheIF()
  {
    delete m_impFunc1;
  }

  Real LatheIF::value(const RealVect& a_point) const
  {
    Real retval;
    Real x = a_point[0];
    Real y = a_point[1];

    // Get r value
    Real r , a;
    r= x*x;
    a= y*y;

    r = sqrt(r+a);

#if AMREX_SPACEDIM == 2
    RealVect coord(r,0.0);

    retval =  m_impFunc1->value(coord);
#elif AMREX_SPACEDIM == 3
    Real z = a_point[2];
    Real r1,z1;

    r1 = r;
    z1 = z;

    RealVect coord2(r1,z1,0.0);

    retval = m_impFunc1->value(coord2);
#endif

    // Change the sign to change inside to outside
    if (!m_inside)
    {
      retval = -retval;
    }

    return retval;
  }


  BaseIF* LatheIF::newImplicitFunction() const
  {
    LatheIF* lathePtr;

    lathePtr = new LatheIF(*m_impFunc1,m_inside);

    return static_cast<BaseIF*>(lathePtr);
  }

}
