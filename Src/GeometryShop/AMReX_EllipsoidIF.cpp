#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EllipsoidIF.H"

#include "NamespaceHeader.H"

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

void EllipsoidIF::GetParams(RealVect& a_radii,
                            RealVect& a_center,
                            bool&     a_inside) const
{
  // Copy parameter information over
  a_radii  = m_radii;
  a_center = m_center;
  a_inside = m_inside;
}

void EllipsoidIF::SetParams(const RealVect& a_radii,
                            const RealVect& a_center,
                            const bool&     a_inside)
{
  // Set parameter information
  m_radii  = a_radii;
  m_center = a_center;
  m_inside = a_inside;

  // Precompute the radii squared
  m_radii2  = m_radii;
  m_radii2 *= m_radii;
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

BaseIF* EllipsoidIF::newImplicitFunction() const
{
  EllipsoidIF* ellipsoidPtr = new EllipsoidIF(m_radii,
                                              m_center,
                                              m_inside);

  return static_cast<BaseIF*>(ellipsoidPtr);
}

#include "NamespaceFooter.H"
