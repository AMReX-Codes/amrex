#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TRANSFORMIF_H_
#define _TRANSFORMIF_H_

#include "MayDay.H"
#include "RealVect.H"
#include "Vector.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function is a tranformed version of the implicit function
    passed in the constructor.  Various transformations (in 2D and 3D) are
    allowed including translating, scaling, and rotations.
 */
class TransformIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying the implicit function
   */
  TransformIF(const BaseIF& a_impFunc);

  /// Copy constructor
  TransformIF(const TransformIF& a_inputIF);

  /// Destructor
  virtual ~TransformIF();

  ///
  /**
      Translate by a_trans
   */
  void translate(const RealVect& a_trans);

  ///
  /**
      Scale (uniformly) by a_scale
   */
  void scale(const Real& a_scale);

  ///
  /**
      Scale in each direction (idir) by a_scale[idir]
   */
  void scale(const RealVect& a_scale);

  ///
  /**
      In 2D - rotate by a_angle about a_point, a_axis is ignored.
      In 3D - rotate by a_angle about the line in space specified by
      a_point and and a_axis.
   */
  void rotate(const Real&     a_angle,
              const RealVect& a_point = RealVect::Zero,
              const RealVect& a_axis  = BASISREALV(0));

  ///
  /**
      About a_point, Rotate vector a_axis1 to align with vector a_axis2
   */
  void rotate(const RealVect& a_axis1,
              const RealVect& a_axis2,
              const RealVect& a_point = RealVect::Zero);

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  Real value(const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  virtual bool fastIntersection(const RealVect& a_lo,
                                const RealVect& a_hi) const ;

  virtual GeometryService::InOut InsideOutside(const RealVect& a_lo,
                                               const RealVect& a_hi) const ;

  ///
  /**
     Pass this call onto the IFs contained in this IF class.
  */
  virtual void boxLayoutChanged(const DisjointBoxLayout & a_newBoxLayout,
                                const RealVect          & a_dx)
  {
    m_impFunc->boxLayoutChanged(a_newBoxLayout,a_dx);
  }

protected:
  ///
  /**
      Constructor for the factory
   */
  TransformIF(const BaseIF& a_impFunc,
              const Real    a_transform[SpaceDim+1][SpaceDim+1],
              const Real    a_invTransform[SpaceDim+1][SpaceDim+1]);

  // transformation (in homogeneous coordinates)
  Real m_transform[SpaceDim+1][SpaceDim+1];

  // inverse transformation (in homogeneous coordinates)
  Real m_invTransform[SpaceDim+1][SpaceDim+1];

  BaseIF* m_impFunc; // implicit function to transform

  void vectorMultiply(RealVect&       m_outPoint,
                      const Real      m_intrans[SpaceDim+1][SpaceDim+1],
                      const RealVect& m_inPoint) const;

  void vectorMultiply(IndexTM<Real,GLOBALDIM>&       m_outPoint,
                      const Real                     m_intrans[SpaceDim+1][SpaceDim+1],
                      const IndexTM<Real,GLOBALDIM>& m_inPoint) const;

  void matrixIdentity(Real m_trans[SpaceDim+1][SpaceDim+1]) const;

  void matrixMultiply(Real       m_outTrans[SpaceDim+1][SpaceDim+1],
                      const Real m_intrans1[SpaceDim+1][SpaceDim+1],
                      const Real m_intrans2[SpaceDim+1][SpaceDim+1]) const;

  void matrixTranslate(Real            m_trans[SpaceDim+1][SpaceDim+1],
                       const RealVect& m_translate) const;

  void matrixScale(Real            m_trans[SpaceDim+1][SpaceDim+1],
                   const RealVect& m_scale) const;

private:
  TransformIF()
  {
    MayDay::Abort("TransformIF uses strong construction");
  }

  void operator=(const TransformIF& a_inputIF)
  {
    MayDay::Abort("TransformIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
