#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _HYPERSPHEREIF_H_
#define _HYPERSPHEREIF_H_

#include "MayDay.H"
#include "IndexTM.H"
#include "Notation.H"
#include "BaseIF.H"
#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a sphere.
 */
class HyperSphereIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying sphere radius (a_radius), center (a_center), and
      whether the domain is on the inside (a_inside).
   */
  HyperSphereIF(const Real                   & a_radius,
                const IndexTM<Real,GLOBALDIM>& a_center,
                const bool                   & a_inside);

  /// Copy constructor
  HyperSphereIF(const HyperSphereIF& a_inputIF);

  /// Destructor
  virtual ~HyperSphereIF();
  ///
  /**
      Return the parameter information
   */
  virtual void GetParams(Real                   & a_radius,
                         IndexTM<Real,GLOBALDIM>& a_center,
                         bool                   & a_inside) const;

  ///
  /**
      Set the parameter information
   */
  virtual void SetParams(const Real                   & a_radius,
                         const IndexTM<Real,GLOBALDIM>& a_center,
                         const bool                   & a_inside);

  ///
  /**
      Return the value of the function at a_point.
   */

  virtual Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                     const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual Real value(const RealVect& a_point) const;

  virtual Real value(const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual IndexTM<Real,GLOBALDIM> normal(const IndexTM<Real,GLOBALDIM>& a_point) const ;

  virtual Vector<IndexTM<Real,GLOBALDIM> > gradNormal(const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

protected:
  Real              m_radius;    // radius
  IndexTM<Real,GLOBALDIM> m_center;    // center
  bool              m_inside;    // inside flag
  Real              m_radius2;   // precomputed radius squared

private:
  HyperSphereIF()
  {
    MayDay::Abort("HyperSphereIF uses strong construction");
  }

  void operator=(const HyperSphereIF& a_inputIF)
  {
    MayDay::Abort("HyperSphereIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
