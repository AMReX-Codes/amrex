#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ELLIPSOIDIF_H_
#define _ELLIPSOIDIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies an ellipsoid.
 */
class EllipsoidIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying ellipsoid "radii" (a_radii), center (a_center),
      and whether the domain is on the inside (a_inside).
   */
  EllipsoidIF(const RealVect& a_radii,
              const RealVect& a_center,
              const bool&     a_inside);

  /// Copy constructor
  EllipsoidIF(const EllipsoidIF& a_inputIF);

  /// Destructor
  virtual ~EllipsoidIF();

  ///
  /**
      Return the parameter information
   */
  virtual void GetParams(RealVect& a_radii,
                         RealVect& a_center,
                         bool&     a_inside) const;

  ///
  /**
      Set the parameter information
   */
  virtual void SetParams(const RealVect& a_radii,
                         const RealVect& a_center,
                         const bool&     a_inside);

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  RealVect m_radii;     // radii
  RealVect m_center;    // center
  bool     m_inside;    // inside flag

  RealVect m_radii2;    // precomputed radii squared

private:
  EllipsoidIF()
  {
    MayDay::Abort("EllipsoidIF uses strong construction");
  }

  void operator=(const EllipsoidIF& a_inputIF)
  {
    MayDay::Abort("EllipsoidIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
