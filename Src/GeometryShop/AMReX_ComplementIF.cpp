#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ComplementIF.H"

#include "NamespaceHeader.H"

ComplementIF::ComplementIF(const BaseIF& a_impFunc,
                           const bool&   a_complement)
{
  m_impFunc = a_impFunc.newImplicitFunction();
  m_complement = a_complement;
}

ComplementIF::ComplementIF(const ComplementIF& a_inputIF,
                           const bool&         a_complement)
{
  m_impFunc = a_inputIF.m_impFunc->newImplicitFunction();
  m_complement = a_complement;
}

ComplementIF::~ComplementIF()
{
  delete m_impFunc;
}

void ComplementIF::GetParams(bool& a_complement) const
{
  // Copy parameter information over
  a_complement = m_complement;
}

void ComplementIF::SetParams(const bool& a_complement)
{
  // Set parameter information
  m_complement = a_complement;
}

Real ComplementIF::value(const RealVect& a_point) const
{
  Real retval;

  // Implicit function value
  retval = m_impFunc->value(a_point);

  // Return the negative if complement is turned on (true)
  if (m_complement)
  {
    retval = -retval;
  }

  return retval;
}

Real ComplementIF::value(const IndexTM<Real,GLOBALDIM>& a_point) const
{

  Real retval;

  // Implicit function value
  retval = m_impFunc->value(a_point);

  // Return the negative if complement is turned on (true)
  if (m_complement)
  {
    retval = -retval;
  }

  return retval;

}

Real ComplementIF::value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                        const IndexTM<Real,GLOBALDIM>& a_point) const
{
  Real retval = m_impFunc->value(a_partialDerivative,a_point);

  // Return the negative if complement is turned on (true)
  if (m_complement)
  {
    retval = -retval;
  }

  return retval;
}

GeometryService::InOut ComplementIF::InsideOutside(const RealVect& a_low, const RealVect& a_high) const
{
  GeometryService::InOut r = m_impFunc->InsideOutside(a_low, a_high);
  if (!m_complement) return r;
  if (r == GeometryService::Regular) return GeometryService::Covered;
  if (r == GeometryService::Covered) return GeometryService::Regular;
  return GeometryService::Irregular;
}



BaseIF* ComplementIF::newImplicitFunction() const
{
  ComplementIF* complementPtr = new ComplementIF(*m_impFunc,m_complement);

  return static_cast<BaseIF*>(complementPtr);
}

#include "NamespaceFooter.H"

