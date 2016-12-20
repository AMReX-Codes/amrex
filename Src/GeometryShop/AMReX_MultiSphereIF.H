#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTISPHEREIF_H_
#define _MULTISPHEREIF_H_

#include "MayDay.H"
#include "RealVect.H"

#include "BaseIF.H"
#include "ComplementIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a union of spheres.
 */
class MultiSphereIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying sphere radii (a_radii), centers (a_centers), and
      whether the domain is on the inside (a_inside).
   */
  MultiSphereIF(const Vector<Real>&     a_radii,
                const Vector<RealVect>& a_centers,
                const bool&             a_inside);

  /// Copy constructor
  MultiSphereIF(const MultiSphereIF& a_inputIF);

  /// Destructor
  virtual ~MultiSphereIF();

  ///
  /**
      Return the parameter information
   */
  virtual void GetParams(Vector<Real>&     a_radii,
                         Vector<RealVect>& a_centers,
                         bool&             a_inside) const;

  ///
  /**
      Set the parameter information
   */
  virtual void SetParams(const Vector<Real>&     a_radii,
                         const Vector<RealVect>& a_centers,
                         const bool&             a_inside);

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  virtual bool fastIntersection(const RealVect& a_low, const RealVect& a_high) const
  {return false;}

  virtual GeometryService::InOut InsideOutside(const RealVect& a_low, const RealVect& a_high) const;

  Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
             const IndexTM<Real,GLOBALDIM>& a_point) const;
protected:
  int              m_numSpheres;  // number of spheres
  Vector<Real>     m_radii;       // radi
  Vector<RealVect> m_centers;     // centers
  bool             m_inside;      // inside flag

  ComplementIF*    m_multiSphere; // the IF object for the spheres

private:
  MultiSphereIF()
  {
    MayDay::Abort("MultiSphereIF uses strong construction");
  }

  void operator=(const MultiSphereIF& a_inputIF)
  {
    MayDay::Abort("MultiSphereIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
