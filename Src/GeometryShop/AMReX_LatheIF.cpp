#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LatheIF.H"
#include "UnionIF.H"

#include "NamespaceHeader.H"

LatheIF::LatheIF(const BaseIF& a_impFunc1,
                 const bool&                 a_inside)
{
  // Make a copy of the implicit function and note only one was given
  m_impFunc1 = a_impFunc1.newImplicitFunction();
  m_impFunc2 = NULL;

  // Save inside flag
  m_inside = a_inside;
}

LatheIF::LatheIF(const BaseIF& a_impFunc1,
                 const BaseIF& a_impFunc2,
                 const RealVect&             a_point,
                 const bool&                 a_inside)
{
  // Make a copy of both implicit functions
  m_impFunc1 = a_impFunc1.newImplicitFunction();
  m_impFunc2 = a_impFunc2.newImplicitFunction();

  // Save the point of rotation
  m_point = a_point;

  // Save inside flag
  m_inside = a_inside;
}

LatheIF::LatheIF(const LatheIF& a_inputIF)
{
  // Make a copy of both implicit functions (the second may be NULL)
  m_impFunc1 = a_inputIF.m_impFunc1->newImplicitFunction();
  m_impFunc2 = a_inputIF.m_impFunc2->newImplicitFunction();

  // Save the point of rotation
  m_point = a_inputIF.m_point;

  // Save inside flag
  m_inside = a_inputIF.m_inside;
}

LatheIF::~LatheIF()
{
  delete m_impFunc1;
  if (m_impFunc2 != NULL)
  {
    delete m_impFunc2;
  }
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

#if CH_SPACEDIM == 2
  RealVect coord(r,0.0);

  retval =  m_impFunc1->value(coord);
#elif CH_SPACEDIM == 3
  Real z = a_point[2];
  Real r1,z1;

  if (m_impFunc2 == NULL)
  {
    r1 = r;
    z1 = z;
  }
  else
  {
    Real theta = atan2(y,x);
    RealVect coord(theta,0.0,0.0);

    Real angle = -m_impFunc2->value(coord);

    r -= m_point[0];
    z -= m_point[1];

    r1 = cos(angle)*r - sin(angle)*z;
    z1 = sin(angle)*r + cos(angle)*z;

    r1 += m_point[0];
    z1 += m_point[1];
  }

  RealVect coord2(r1,z1,0.0);

  retval = m_impFunc1->value(coord2);
#else
  MayDay::Abort("need higher dim in LatheIF\n");
#endif

  // Change the sign to change inside to outside
  if (!m_inside)
  {
    retval = -retval;
  }

  return retval;
}

GeometryService::InOut LatheIF::InsideOutside(const RealVect& lo, const RealVect& hi) const
{
  GeometryService::InOut rtn = GeometryService::Irregular;
#if CH_SPACEDIM == 3
  //project the 3D box defined by (lo,hi) onto box in x-z plane, (loPlane,hiPlane)
  Vector<Real> r;
  r.resize(4);
  r[0] = sqrt(lo[0]*lo[0]+lo[1]*lo[1]);
  r[1] = sqrt(lo[0]*lo[0]+hi[1]*hi[1]);
  r[2] = sqrt(hi[0]*hi[0]+lo[1]*lo[1]);
  r[3] = sqrt(hi[0]*hi[0]+hi[1]*hi[1]);
  RealVect loPlane,hiPlane;
//   loPlane[1] = lo[1];//in constant y-plane
//   hiPlane[1] = hi[1];//in constant y-plane
  loPlane[1] = 0.;//in constant y-plane
  hiPlane[1] = 0.;//in constant y-plane
  loPlane[2] = lo[2];
  hiPlane[2] = hi[2];
  if (m_impFunc2 == NULL)
    {
      loPlane[0] = r[0];
      hiPlane[0] = r[0];
      for (int i = 1; i < 4; i++)
        {
          if (r[i] < loPlane[0])
          {
            loPlane[0] = r[i];
          }
          if (r[i] > hiPlane[0])
          {
            hiPlane[0] = r[i];
          }
        }
    }
  else
    {
      MayDay::Abort("LatheIF::InsideOutside should not get here if m_impFunc2 != NULL\n");
    }
//   //now check intersection of projected box with surface (just as in PlaneIF)
//   RealVect len = hiPlane-loPlane;
//   Real firstValue = value(loPlane);
//   GeometryService::InOut rtn;
//   if ( firstValue < 0 )
//     {
//       rtn = GeometryService::Regular;
//     }
//   else
//     {
//       rtn = GeometryService::Covered;
//     }
//   Box unit(IntVect::Zero, IntVect::Unit);
//   for (BoxIterator b(unit); b.ok(); ++b)
//     {
//       RealVect corner = loPlane + RealVect(b())*len;
//       Real functionValue = value(corner);

//       if (functionValue * firstValue <= 0.0 )
//       {
//         return GeometryService::Irregular;
//       }
//     }
  rtn = m_impFunc1->InsideOutside(loPlane,hiPlane);
#else
  MayDay::Abort("LatheIF insideOutside for fastIntersection only works in 3D\n");
#endif
  return rtn;
}

BaseIF* LatheIF::newImplicitFunction() const
{
  LatheIF* lathePtr;

  if (m_impFunc2 == NULL)
  {
   lathePtr = new LatheIF(*m_impFunc1,m_inside);
  }
  else
  {
   lathePtr = new LatheIF(*m_impFunc1,*m_impFunc2,m_point,m_inside);
  }

  return static_cast<BaseIF*>(lathePtr);
}

#include "NamespaceFooter.H"
