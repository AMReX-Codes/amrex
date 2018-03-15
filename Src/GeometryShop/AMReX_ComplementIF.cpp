#include "AMReX_ComplementIF.H"

namespace amrex
{

  ComplementIF::ComplementIF(const BaseIF& a_impFunc)
  {
    m_impFunc = a_impFunc.newImplicitFunction();
  }


  ComplementIF::~ComplementIF()
  {
    delete m_impFunc;
  }


  Real ComplementIF::value(const RealVect& a_point) const
  {
    Real retval;

    // Implicit function value
    retval = m_impFunc->value(a_point);

    // Return the negative because  this is the complement 
    retval = -retval;

    return retval;
  }


  BaseIF* ComplementIF::newImplicitFunction() const
  {
    ComplementIF* complementPtr = new ComplementIF(*m_impFunc);

    return static_cast<BaseIF*>(complementPtr);
  }
}


