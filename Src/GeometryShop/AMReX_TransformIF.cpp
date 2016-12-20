#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PolyGeom.H"
#include "TransformIF.H"

#include "NamespaceHeader.H"


TransformIF::TransformIF(const BaseIF& a_impFunc)
{
  m_impFunc = a_impFunc.newImplicitFunction();

  matrixIdentity(m_transform);
  matrixIdentity(m_invTransform);
}

TransformIF::TransformIF(const BaseIF& a_impFunc,
                         const Real    a_transform[SpaceDim+1][SpaceDim+1],
                         const Real    a_invTransform[SpaceDim+1][SpaceDim+1])
{
  m_impFunc = a_impFunc.newImplicitFunction();

  // Copy the matrices over
  for (int i = 0; i <= SpaceDim; i++)
  {
    for (int j = 0; j <= SpaceDim; j++)
    {
      m_transform   [i][j] = a_transform   [i][j];
      m_invTransform[i][j] = a_invTransform[i][j];
    }
  }
}

TransformIF::TransformIF(const TransformIF& a_inputIF)
{
  m_impFunc = a_inputIF.m_impFunc->newImplicitFunction();

  // Copy the matrices over
  for (int i = 0; i <= SpaceDim; i++)
  {
    for (int j = 0; j <= SpaceDim; j++)
    {
      m_transform   [i][j] = a_inputIF.m_transform   [i][j];
      m_invTransform[i][j] = a_inputIF.m_invTransform[i][j];
    }
  }
}

TransformIF::~TransformIF()
{
  delete m_impFunc;
}

void TransformIF::translate(const RealVect& a_trans)
{
  // Temporary transformation matrix
  Real temp[SpaceDim+1][SpaceDim+1];

  // Translation
  RealVect trans = a_trans;

  // Create the forward transformation matrix
  matrixTranslate(temp,trans);

  // Update the overall forward transform
  matrixMultiply(m_transform,temp,m_transform);

  // Inverse translation
  trans = RealVect::Zero;
  trans -= a_trans;

  // Create the inverse transformation matrix
  matrixTranslate(temp,trans);

  // Update the overall inverse transform
  matrixMultiply(m_invTransform,m_invTransform,temp);
}

void TransformIF::scale(const Real& a_scale)
{
  // Temporary transformation matrix
  Real temp[SpaceDim+1][SpaceDim+1];

  // Set up a scaling vector
  RealVect scale(D_DECL(a_scale,a_scale,a_scale));

  // Create the forward transformation matrix
  matrixScale(temp,scale);

  // Update the overall forward transform
  matrixMultiply(m_transform,temp,m_transform);

  // Set up an inverse scaling vector
  RealVect invScale(D_DECL(1.0/a_scale,1.0/a_scale,1.0/a_scale));

  // Create the inverse transformation matrix
  matrixScale(temp,invScale);

  // Update the overall inverse transform
  matrixMultiply(m_invTransform,m_invTransform,temp);
}

void TransformIF::scale(const RealVect& a_scale)
{
  // Temporary transformation matrix
  Real temp[SpaceDim+1][SpaceDim+1];

  // Set up a scaling vector
  RealVect scale = a_scale;

  // Create the forward transformation matrix
  matrixScale(temp,scale);

  // Update the overall forward transform
  matrixMultiply(m_transform,temp,m_transform);

  // Inverse scaling
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    scale[idir] = 1.0 / scale[idir];
  }

  // Create the inverse transformation matrix
  matrixScale(temp,scale);

  // Update the overall inverse transform
  matrixMultiply(m_invTransform,m_invTransform,temp);
}

void TransformIF::rotate(const Real&     a_angle,
                         const RealVect& a_point,
                         const RealVect& a_axis)
{
  // Temporary transformation matrix
  Real temp[SpaceDim+1][SpaceDim+1];

  // Transformation matrix for forward transformation
  Real forTrans[SpaceDim+1][SpaceDim+1];

  // Transformation matrix for inverse transformation
  Real invTrans[SpaceDim+1][SpaceDim+1];

  // Start by translating a_point to the origin
  RealVect invTranslate = RealVect::Zero;
  invTranslate -= a_point;

  // Update the forward and inverse transformations
  matrixTranslate(forTrans,invTranslate);
  matrixTranslate(invTrans,invTranslate);

#if CH_SPACEDIM == 3
  // Rotate everything so a_axis aligns with the z-axis
  Real theta = 0.0;
  Real sinTheta = 0.0;
  Real cosTheta = 1.0;

  Real phi = 0.0;
  Real sinPhi = 0.0;
  Real cosPhi = 1.0;

  if (a_axis[0] != 0)
  {
    // Rotate about the z-axis to make the x component of a_axis zero
    theta = atan2(a_axis[0],a_axis[1]);
    sinTheta = sin(theta);
    cosTheta = cos(theta);

    // Set up the rotation matrix
    matrixIdentity(temp);

    temp[0][0] =  cosTheta;
    temp[0][1] = -sinTheta;
    temp[1][0] =  sinTheta;
    temp[1][1] =  cosTheta;

    // Update the forward and inverse transformations
    matrixMultiply(forTrans,temp,forTrans);
    matrixMultiply(invTrans,temp,invTrans);
  }

  Real newY = sinTheta*a_axis[0] + cosTheta*a_axis[1];
  Real newZ = a_axis[2];

  // Rotate about the x-axis to make the y component of a_axis zero
  phi = atan2(newY,newZ);
  sinPhi = sin(phi);
  cosPhi = cos(phi);

  // Set up the rotation matrix
  matrixIdentity(temp);

  temp[1][1] =  cosPhi;
  temp[1][2] = -sinPhi;
  temp[2][1] =  sinPhi;
  temp[2][2] =  cosPhi;

  // Update the forward and inverse transformations
  matrixMultiply(forTrans,temp,forTrans);
  matrixMultiply(invTrans,temp,invTrans);
#endif

  // Now do a rotation in the xy plane(s)
  matrixIdentity(temp);

  temp[0][0] =  cos(a_angle);
  temp[0][1] = -sin(a_angle);
  temp[1][0] =  sin(a_angle);
  temp[1][1] =  cos(a_angle);

  // Update the forward transformation
  matrixMultiply(forTrans,temp,forTrans);

  // Now do a negative rotation in the xy plane(s)
  matrixIdentity(temp);

  temp[0][0] =  cos(-a_angle);
  temp[0][1] = -sin(-a_angle);
  temp[1][0] =  sin(-a_angle);
  temp[1][1] =  cos(-a_angle);

  // Update the inverse transformation
  matrixMultiply(invTrans,temp,invTrans);

#if CH_SPACEDIM == 3
  // Set up the rotation matrix
  matrixIdentity(temp);

  temp[1][1] =  cosPhi;
  temp[1][2] =  sinPhi;
  temp[2][1] = -sinPhi;
  temp[2][2] =  cosPhi;

  // Update the forward and inverse transformations
  matrixMultiply(forTrans,temp,forTrans);
  matrixMultiply(invTrans,temp,invTrans);

  if (a_axis[0] != 0)
  {
    // Set up the rotation matrix
    matrixIdentity(temp);

    temp[0][0] =  cosTheta;
    temp[0][1] =  sinTheta;
    temp[1][0] = -sinTheta;
    temp[1][1] =  cosTheta;

    // Update the forward and inverse transformations
    matrixMultiply(forTrans,temp,forTrans);
    matrixMultiply(invTrans,temp,invTrans);
  }
#endif

  // Now translate the origin back to a_point
  matrixTranslate(temp,a_point);

  // Update the forward and inverse transformations
  matrixMultiply(forTrans,temp,forTrans);
  matrixMultiply(invTrans,temp,invTrans);

  // Update the overall forward transform
  matrixMultiply(m_transform,forTrans,m_transform);

  // Update the overall inverse transform
  matrixMultiply(m_invTransform,m_invTransform,invTrans);
}

void TransformIF::rotate(const RealVect& a_axis1,
                         const RealVect& a_axis2,
                         const RealVect& a_point)
{
  RealVect axis3 = PolyGeom::cross(a_axis1,a_axis2);
  // Cross product seems to be backwards
  axis3 *= -1.0;

  Real axis1Length = sqrt(PolyGeom::dot(a_axis1,a_axis1));
  Real axis2Length = sqrt(PolyGeom::dot(a_axis2,a_axis2));
  Real axis3Length = sqrt(PolyGeom::dot(axis3,axis3));

  Real angle = acos(PolyGeom::dot(a_axis1,a_axis2)/(axis1Length*axis2Length));

#if CH_SPACEDIM == 2
  if (axis3[0] < 0.0)
  {
    angle *= -1.0;
  }
#endif

  if (axis3Length != 0.0)
  {
    rotate(angle,a_point,axis3);
  }
}

Real TransformIF::value(const RealVect& a_point) const
{
  // Function value at the inverse transform of a_point
  Real retval;

  // The inverse transform of a_point
  RealVect invPoint;

  // Inverse transform a_point to invPoint
  vectorMultiply(invPoint,m_invTransform,a_point);

  // Get the function value at invPoint
  retval = m_impFunc->value(invPoint);

  return retval;
}

Real TransformIF::value(const IndexTM<Real,GLOBALDIM>& a_point) const
{
  RealVect point;
  for (int idir = 0 ; idir < GLOBALDIM ; idir++)
    {
      point[idir] = a_point[idir];
    }

  return value(point);
}

BaseIF* TransformIF::newImplicitFunction() const
{
  TransformIF* transformPtr = new TransformIF(*m_impFunc,
                                              m_transform,
                                              m_invTransform);

  return static_cast<BaseIF*>(transformPtr);
}

bool TransformIF::fastIntersection(const RealVect& a_lo,
                                   const RealVect& a_hi) const
{
//   return m_impFunc->fastIntersection(a_lo,a_hi);
  return false;
}

GeometryService::InOut TransformIF::InsideOutside(const RealVect& a_lo,
                                                  const RealVect& a_hi) const
{
  // The size of the incoming bounding box
  RealVect len;
  len = a_hi - a_lo;

  // The bounding box corners in the pre-transform space
  RealVect newLo,newHi;

  // Get the lowest corner
  RealVect corner = a_lo;

  // The inverse transform of corner
  RealVect invCorner;

  // Inverse transform corner to invCorner
  vectorMultiply(invCorner,m_invTransform,corner);

  // Initialize the new corners
  newLo = invCorner;
  newHi = invCorner;

  // Build a unit box for iteration
  Box unit(IntVect::Zero, IntVect::Unit);

  // Iterate over all corners
  for (BoxIterator b(unit); b.ok(); ++b)
  {
    // Get the current corner
    corner = a_lo + RealVect(b())*len;

    // Inverse transform corner to invCorner
    vectorMultiply(invCorner,m_invTransform,corner);

    // Build the bounding box in the pre-transform space
    newLo.min(invCorner);
    newHi.max(invCorner);
  }

  // Report "InOut" using the pre-transform bounding box and the original
  // implicit function
  return m_impFunc->InsideOutside(newLo,newHi);
}

void TransformIF::vectorMultiply(IndexTM<Real,GLOBALDIM>& m_outPoint,
                                 const Real m_intrans[SpaceDim+1][SpaceDim+1],
                                 const IndexTM<Real,GLOBALDIM>& m_inPoint) const
{
  RealVect a_outPoint;
  RealVect a_inPoint;
  for (int idir = 0 ; idir < GLOBALDIM ; idir++)
    {
      a_outPoint[idir] = m_outPoint[idir];
      a_inPoint[idir] = m_inPoint[idir];
    }
  vectorMultiply(a_outPoint,m_intrans,a_inPoint);
  for (int idir = 0 ; idir < GLOBALDIM ; idir++)
    {
      m_outPoint[idir] = a_outPoint[idir];
    }
}

void TransformIF::vectorMultiply(RealVect&       m_outPoint,
                                 const Real      m_inTrans[SpaceDim+1][SpaceDim+1],
                                 const RealVect& m_inPoint) const
{
  // Temporary for input and output vector (in homogeneous coordinates)
  Real in[SpaceDim+1];
  Real out[SpaceDim+1];

  // Multiply the two inputs into the temporary matrix
  for (int i = 0; i < SpaceDim; i++)
  {
    in[i] = m_inPoint[i];
  }
  in[SpaceDim] = 1.0;

  // Multiply the vector by the matrix to get a vector
  for (int i = 0; i <= SpaceDim; i++)
  {
    Real total = 0.0;

    for (int k = 0; k <= SpaceDim; k++)
    {
      total += m_inTrans[i][k] * in[k];
    }

    out[i] = total;
  }

  Real scaling = 1.0 / out[SpaceDim];

  // Convert the transformed vector in homogeneous coordinates to a RaalVect
  for (int i = 0; i < SpaceDim; i++)
  {
      m_outPoint[i] = scaling * out[i];
  }
}

void TransformIF::matrixIdentity(Real m_trans[SpaceDim+1][SpaceDim+1]) const
{
  // Set up the identify
  for (int i = 0; i <= SpaceDim; i++)
  {
    for (int j = 0; j <= SpaceDim; j++)
    {
      m_trans[i][j] = 0.0;
    }

    m_trans[i][i] = 1.0;
  }
}

void TransformIF::matrixMultiply(Real       m_outTrans[SpaceDim+1][SpaceDim+1],
                                 const Real m_inTrans1[SpaceDim+1][SpaceDim+1],
                                 const Real m_inTrans2[SpaceDim+1][SpaceDim+1]) const
{
  // Temporary transformation matrix
  Real temp[SpaceDim+1][SpaceDim+1];

  // Multiply the two inputs into the temporary matrix
  for (int i = 0; i <= SpaceDim; i++)
  {
    for (int j = 0; j <= SpaceDim; j++)
    {
      Real total = 0.0;

      for (int k = 0; k <= SpaceDim; k++)
      {
        total += m_inTrans1[i][k] * m_inTrans2[k][j];
      }

      temp[i][j] = total;
    }
  }

  // Copy the temporary matrix into the returned matrix
  for (int i = 0; i <= SpaceDim; i++)
  {
    for (int j = 0; j <= SpaceDim; j++)
    {
      m_outTrans[i][j] = temp[i][j];
    }
  }
}

void TransformIF::matrixTranslate(Real            a_trans[SpaceDim+1][SpaceDim+1],
                                  const RealVect& a_translate) const
{
  // Start with the identity
  matrixIdentity(a_trans);

  // Set appropriate entries
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    a_trans[idir][SpaceDim] = a_translate[idir];
  }
}

void TransformIF::matrixScale(Real            a_trans[SpaceDim+1][SpaceDim+1],
                              const RealVect& a_scale) const
{
  // Start with the identity
  matrixIdentity(a_trans);

  // Set appropriate entries
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    a_trans[idir][idir] = a_scale[idir];
  }
}
#include "NamespaceFooter.H"
