#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include <cmath>

#include "BaseIF.H"
#include "HyperSphereIF.H"

#include "NamespaceHeader.H"

HyperSphereIF::HyperSphereIF(const Real                   & a_radius,
                             const IndexTM<Real,GLOBALDIM>& a_center,
                             const bool                   & a_inside)
{
  // Remember the parameters
  m_radius = a_radius;
  m_center = a_center;
  m_inside = a_inside;



  // Precompute the radius squared
  m_radius2 = m_radius * m_radius;
}

HyperSphereIF::HyperSphereIF(const HyperSphereIF& a_inputIF)
{
  // Remember the parameters
  m_radius = a_inputIF.m_radius;
  m_center = a_inputIF.m_center;
  m_inside = a_inputIF.m_inside;

  // Precompute the radius squared
  m_radius2 = m_radius * m_radius;
}

HyperSphereIF::~HyperSphereIF()
{
}

void HyperSphereIF::GetParams(Real                   & a_radius,
                              IndexTM<Real,GLOBALDIM>& a_center,
                              bool                   & a_inside) const
{
  // Copy parameter information over
  a_radius = m_radius;
  a_center = m_center;
  a_inside = m_inside;

}

void HyperSphereIF::SetParams(const Real                   & a_radius,
                              const IndexTM<Real,GLOBALDIM>& a_center,
                              const bool                   & a_inside)
{
  // Set parameter information
  m_radius = a_radius;
  m_center = a_center;
  m_inside = a_inside;

  // Precompute the radius squared
  m_radius2 = m_radius * m_radius;
}
Real HyperSphereIF::value(const RealVect & a_point) const
{
  IndexTM<Real,GLOBALDIM> pt;

  if (GLOBALDIM == 3 && SpaceDim == 2)
    {
      MayDay::Abort("HyperPlaneIF should be wrapped in ReferenceHeightIF when GLOBALDIM==3 and SpaceDim==2");
    }
  else
    {
      for (int idir = 0; idir < SpaceDim; ++idir)
        {
          pt[idir] = a_point[idir];
        }
    }
  return value(pt);
}

Real HyperSphereIF::value(const IndexTM<Real,GLOBALDIM> & a_point) const
{
  Real retval;

  // The distance squared for m_center to a_point
  Real distance2;

  // Compute the distance squared
  distance2 = 0.0;
  for (int idir = 0; idir < GLOBALDIM; idir++)
  {
    Real cur;
    cur = a_point[idir] - m_center[idir];

    distance2 += cur*cur;
  }

  // Return the difference between the squares (zero on the sphere)
  retval = distance2 - m_radius2;

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

Real HyperSphereIF::value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                          const IndexTM<Real,GLOBALDIM>& a_point) const
{
  Real retval= LARGEREALVAL;
  int maxDir = a_partialDerivative.maxDir(false);
  int derivativeOrder = a_partialDerivative.sum();

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
      if (a_partialDerivative[maxDir] == 2)
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

IndexTM<Real,GLOBALDIM> HyperSphereIF::normal(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  IndexTM<Real,GLOBALDIM> normal;
  Real norm = 0.0;
  for (int idir=0; idir<GLOBALDIM; idir++)
    {
      normal[idir] = a_point[idir] - m_center[idir];
      norm += normal[idir]*normal[idir];
    }
  norm = sqrt(norm);
  normal /= norm;
  if (!m_inside)
    {
      return -normal;
    }
  return normal ;
}

Vector<IndexTM<Real,GLOBALDIM> > HyperSphereIF::gradNormal(const IndexTM<Real,GLOBALDIM>& a_point)const
{
  Vector<IndexTM <Real,GLOBALDIM> >gradNorm;
  gradNorm.resize(GLOBALDIM);
  Real norm = 0.0;
  for (int kdir = 0 ; kdir < GLOBALDIM ; kdir++)
    {
      norm += (a_point[kdir] - m_center[kdir]) * (a_point[kdir] - m_center[kdir]);
    }
  norm = pow(norm,1.5);

  for (int idir=0; idir<GLOBALDIM; idir++)
    {
      IndexTM<Real,GLOBALDIM> vect = IndexTM<Real,GLOBALDIM>::Zero;
      for (int jdir=0; jdir<GLOBALDIM; jdir++)
        {
          if (idir==jdir)
            {
              for (int kdir = 0 ; kdir < GLOBALDIM ; kdir++)
                {
                  if (kdir != idir)
                    {
                      vect[jdir] += (a_point[kdir] - m_center[kdir]) * (a_point[kdir] - m_center[kdir]);
                    }
                }
            }
          else
            {
              vect[jdir] = -(a_point[idir] - m_center[idir]) * (a_point[jdir] - m_center[jdir]);
            }
        }
      vect /= norm;
      gradNorm[idir] = vect;
    }
  if (!m_inside)
    {
      for (int idir = 0 ; idir < GLOBALDIM ; idir++)
        {
          gradNorm[idir] = -gradNorm[idir];
        }
    }
  return gradNorm;

}

 BaseIF* HyperSphereIF::newImplicitFunction() const
{
  HyperSphereIF* spherePtr = new HyperSphereIF(m_radius,
                                               m_center,
                                               m_inside);

  return static_cast<BaseIF*>(spherePtr);
}


#include "NamespaceFooter.H"
