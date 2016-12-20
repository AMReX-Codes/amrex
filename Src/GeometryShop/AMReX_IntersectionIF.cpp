#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "IntersectionIF.H"

#include "NamespaceHeader.H"

IntersectionIF::IntersectionIF(const BaseIF& a_impFunc1,
                               const BaseIF& a_impFunc2)
{
  // Number of implicit function in intersection
  m_numFuncs = 2;

  // Vector of implicit function pointers
  m_impFuncs.resize(m_numFuncs);

  // Make copies of the implicit functions
  m_impFuncs[0] = a_impFunc1.newImplicitFunction();
  m_impFuncs[1] = a_impFunc2.newImplicitFunction();
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

IntersectionIF::IntersectionIF(const IntersectionIF& a_inputIF)
{
  // Number of implicit function in intersection
  m_numFuncs = a_inputIF.m_impFuncs.size();

  // Vector of implicit function pointers
  m_impFuncs.resize(m_numFuncs);



  // Make copies of the implicit functions
  for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
  {
    if (a_inputIF.m_impFuncs[ifunc] == NULL)
    {
      m_impFuncs[ifunc] = NULL;
    }
    else
    {
      m_impFuncs[ifunc] = a_inputIF.m_impFuncs[ifunc]->newImplicitFunction();
    }
  }
}

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

Real IntersectionIF::value(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  int closestIF = -1;
  findClosest(a_point,closestIF);

  if (closestIF == -1)
  {
    return -1.0;
  }
  else
  {
    return m_impFuncs[closestIF]->value(a_point);
  }
}

Real IntersectionIF::value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                           const IndexTM<Real,GLOBALDIM>& a_point) const
{
  int closestIF = -1;
  findClosest(a_point,closestIF);

  if (closestIF == -1)
  {
    if (a_partialDerivative.sum() == 0)
    {
      return -1.0;
    }
    else
    {
      return 0.0;
    }
  }
  else
  {
    return m_impFuncs[closestIF]->value(a_partialDerivative,a_point);
  }
}

BaseIF* IntersectionIF::newImplicitFunction() const
{
  IntersectionIF* intersectionPtr = new IntersectionIF(m_impFuncs);

  return static_cast<BaseIF*>(intersectionPtr);
}

void IntersectionIF::findClosest(const IndexTM<Real,GLOBALDIM> & a_point,
                                 int                           & a_closestIF) const
{
  Real retval = 0.0;

  if (m_numFuncs > 0)
    {
      retval = m_impFuncs[0]->value(a_point);
      a_closestIF = 0;

      for (int ifunc = 1; ifunc < m_numFuncs; ifunc++)
        {
          Real cur;
          cur = m_impFuncs[ifunc]->value(a_point);
          if (cur > retval)
            {
              retval = cur;
              a_closestIF = ifunc;
            }
        }
    }
}

bool IntersectionIF::fastIntersection(const RealVect& a_low,
                                      const RealVect& a_high) const
{
  for (int i=0; i<m_impFuncs.size(); i++)
    {
      if (!m_impFuncs[i]->fastIntersection(a_low, a_high)) return false;
    }
  return true;
}

GeometryService::InOut IntersectionIF::InsideOutside(const RealVect& a_low,
                                                     const RealVect& a_high) const
{
  CH_assert(fastIntersection(a_low, a_high));

  bool allRegular = true;
  for (int i = 0; i<m_numFuncs; ++i)
    {
      GeometryService::InOut r =  m_impFuncs[i]->InsideOutside(a_low, a_high);
      if (r == GeometryService::Covered) return r;
      if (r == GeometryService::Irregular) allRegular = false;
    }
  if (allRegular) return GeometryService::Regular;
  return GeometryService::Irregular;
}

#include "NamespaceFooter.H"
