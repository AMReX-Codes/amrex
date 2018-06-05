#include "AMReX_IntersectionIF.H"

namespace amrex
{
  IntersectionIF::~IntersectionIF()
  {
    // Delete all the copies
    for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
    {
      if (m_impFuncs[ifunc] != NULL)
      {
        delete m_impFuncs[ifunc];
      }
    }
  }

  Real IntersectionIF::value(const RealVect& a_point) const
  {

    // Maximum of the implicit functions values
    Real retval;

    retval = -1.0;

    // Find the maximum value and return it
    if (m_numFuncs > 0)
    {
      retval = m_impFuncs[0]->value(a_point);

      for (int ifunc = 1; ifunc < m_numFuncs; ifunc++)
      {
        Real cur;

        cur = m_impFuncs[ifunc]->value(a_point);
        if (cur > retval)
        {
          retval = cur;
        }
      }
    }

    return retval;
  }


  Real 
  IntersectionIF::
  derivative(const  IntVect& a_deriv,
             const RealVect& a_point) const
  {
    // Maximum of the implicit functions values
    Real retval;

    retval = -1.0;

    // Find the maximum value and return it
    if (m_numFuncs > 0)
    {

      int whichfunc = 0;
      Real funcval = m_impFuncs[0]->value(a_point);

      for (int ifunc = 1; ifunc < m_numFuncs; ifunc++)
      {
        Real cur = m_impFuncs[ifunc]->value(a_point);
        if (cur > funcval)
        {
          funcval = cur;
          whichfunc = ifunc;
        }
      }
      retval = m_impFuncs[whichfunc]->derivative(a_deriv, a_point);
    }

    return retval;
  }

  BaseIF* IntersectionIF::newImplicitFunction() const
  {
    IntersectionIF* intersectionPtr = new IntersectionIF(m_impFuncs);

    return static_cast<BaseIF*>(intersectionPtr);
  }

  IntersectionIF::IntersectionIF(const Vector<BaseIF *>& a_impFuncs)
  {
    // Number of implicit function in intersection
    m_numFuncs = a_impFuncs.size();

    // Vector of implicit function pointers
    m_impFuncs.resize(m_numFuncs);

    // Make copies of the implicit functions

    for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
    {
      if (a_impFuncs[ifunc] == NULL)
      {
        m_impFuncs[ifunc] = NULL;
      }
      else
      {
        m_impFuncs[ifunc] = a_impFuncs[ifunc]->newImplicitFunction();
      }
    }
  }
}
