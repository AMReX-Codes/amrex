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

#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include <cmath>

#include "BaseIF.H"
#include "HyperPlaneIF.H"

#include "NamespaceHeader.H"

HyperPlaneIF::HyperPlaneIF(const IndexTM<Real,GLOBALDIM> & a_normal,
                           const IndexTM<Real,GLOBALDIM> & a_point,
                           const bool                    & a_normalIn)
 :m_normal(a_normal),
  m_point(a_point),
  m_normalIn(a_normalIn)
{
}

HyperPlaneIF::HyperPlaneIF(const HyperPlaneIF& a_inputIF)
{
  // Remember the parameters
  m_normal   = a_inputIF.m_normal;
  m_point    = a_inputIF.m_point;
  m_normalIn = a_inputIF.m_normalIn;
}

HyperPlaneIF::~HyperPlaneIF()
{
}

Real HyperPlaneIF::value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                         const IndexTM<Real,GLOBALDIM>& a_point) const
{
  int order = a_partialDerivative.sum();

  Real retval = LARGEREALVAL;

  if (order == 0)
    {
      retval = value(a_point);
    }
  else if (order == 1)
    {
      bool doAbs = true;
      retval = m_normal[a_partialDerivative.maxDir(doAbs)];

      if (m_normalIn)
        {
          retval = -retval;
        }
    }
  else
    {
      retval = 0.0;
    }

  return retval;
}

Real HyperPlaneIF::value(const RealVect & a_point) const
{
  IndexTM<Real,GLOBALDIM> pt;

  //check does SpaceDim = GLOBALDIM
  if (GLOBALDIM==3 && SpaceDim==2)
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

Real HyperPlaneIF::value(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  Real retval = 0.0;

  for (int idir = 0 ; idir < GLOBALDIM ; idir++)
    {
      retval += (a_point[idir] - m_point[idir]) * m_normal[idir];
    }

  if (m_normalIn)
    {
      retval = -retval;
    }

  return retval;
}

IndexTM<Real,GLOBALDIM> HyperPlaneIF::normal(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  IndexTM<Real,GLOBALDIM> normal = m_normal;

  Real norm = 0.0;
  for (int idir=0 ; idir < GLOBALDIM ; idir++)
    {
      norm += m_normal[idir] * m_normal[idir];
    }

  normal /= sqrt(norm);

  if (m_normalIn)
    {
      normal = -normal;
    }

  return normal;
}

Vector<IndexTM<Real,GLOBALDIM> > HyperPlaneIF::gradNormal(const IndexTM<Real,GLOBALDIM>& a_point)const
{
  Vector<IndexTM<Real,GLOBALDIM> > gradNorm;
  gradNorm.resize(GLOBALDIM);

  for (int idir = 0 ; idir < GLOBALDIM ; idir++)
    {
      IndexTM<Real,GLOBALDIM> vect;
      for (int jdir = 0 ; jdir < GLOBALDIM ; jdir++)
        {
          vect[jdir] = 0.0;
        }
      gradNorm[idir] = vect;
    }

  return gradNorm;
}

BaseIF* HyperPlaneIF::newImplicitFunction() const
{
  HyperPlaneIF* hyperPlanePtr = new HyperPlaneIF(m_normal,
                                                 m_point,
                                                 m_normalIn);

  return static_cast<BaseIF*>(hyperPlanePtr);
}

#include "NamespaceFooter.H"
