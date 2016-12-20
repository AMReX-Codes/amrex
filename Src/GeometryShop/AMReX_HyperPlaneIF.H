#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _HYPERPLANEIF_H_
#define _HYPERPLANEIF_H_

#include "MayDay.H"
#include "IndexTM.H"
#include "REAL.H"

#include "BaseIF.H"
#include "NamespaceHeader.H"

///
/**
    This implicit function specifies a plane templated in dimension.
 */
class HyperPlaneIF: public BaseIF
{
public:
  ///
  /**
     Constructor specifying the coefficients of the plane normal, a_normal,
     and a point on the plane, a_point.
    */
  HyperPlaneIF(const IndexTM<Real,GLOBALDIM>& a_normal,
               const IndexTM<Real,GLOBALDIM>& a_point,
               const bool                   & a_normalIn);

  /// Copy constructor
  HyperPlaneIF(const HyperPlaneIF& a_inputIF);

  /// Destructor
  virtual ~HyperPlaneIF();
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

  //protected:

  // normal to the plane
  IndexTM<Real,GLOBALDIM> m_normal;

  // a point on the plane
  IndexTM<Real,GLOBALDIM> m_point;

  // (pt - m_point) . m_normal > 0 is inside if m_normalIn = true
  bool m_normalIn;

private:
  HyperPlaneIF()
  {
    MayDay::Abort("HyperPlaneIF uses strong construction");
  }

  void operator=(const HyperPlaneIF& a_inputIF)
  {
    MayDay::Abort("HyperPlaneIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
