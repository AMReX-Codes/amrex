#include "AMReX_UnionIF.H"
#include <AMReX_Array.H>
#include <AMReX_Vector.H>

namespace amrex
{

  UnionIF::UnionIF(const Vector<BaseIF *>& a_impFuncs)
  {
    // Number of implicit function in union
    m_numFuncs = a_impFuncs.size();

    // vector of implicit function pointers
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


  UnionIF::~UnionIF()
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

  Real UnionIF::value(const RealVect& a_point) const
  {
    // Minimum of the implicit functions values
    Real retval;

    retval = 1.0;

    // Find the minimum value and return it
    if (m_numFuncs > 0)
    {
      retval = m_impFuncs[0]->value(a_point);

      for (int ifunc = 1; ifunc < m_numFuncs; ifunc++)
      {
        Real cur;

        cur = m_impFuncs[ifunc]->value(a_point);
        if (cur < retval)
        {
          retval = cur;
        }
      }
    }

    return retval;
  }

  BaseIF* UnionIF::newImplicitFunction() const
  {
    UnionIF* unionPtr = new UnionIF(m_impFuncs);

    return static_cast<BaseIF*>(unionPtr);
  }

}

