#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _COMPLEMENTIF_H_
#define _COMPLEMENTIF_H_

#include "MayDay.H"
#include "RealVect.H"
#include "Vector.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function is the negative of the implicit function with which
    it is constructed.  This results in an implicit function whose interior is
    the complement of the interior of the given implicit function (minus the
    the surface).
 */
class ComplementIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying the implicit function and whether to complement
      (default = true)
   */
  ComplementIF(const BaseIF& a_impFunc,
               const bool&   a_complement = true);

  /// Copy constructor with complement boolean (default = true)
  ComplementIF(const ComplementIF& a_inputIF,
               const bool&         a_complement = true);

  /// Destructor
  virtual ~ComplementIF();

  ///
  /**
      Return the parameter information
   */
  virtual void GetParams(bool& a_complement) const;

  ///
  /**
      Set the parameter information
   */
  virtual void SetParams(const bool& a_complement);

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  ///
  /**
      Return the value of the function at a_point (of type IndexTM).
   */
  virtual Real value(const IndexTM<Real,GLOBALDIM>& a_point) const;

 ///
  /**
     Return the derivative of the function at a_point (of type IndexTM).

  */
  virtual Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                     const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  virtual bool fastIntersection(const RealVect& a_low, const RealVect& a_high) const
  {
    return m_impFunc->fastIntersection(a_low, a_high);
  }

  virtual GeometryService::InOut InsideOutside(const RealVect& a_low, const RealVect& a_high) const;

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
  bool m_complement;               // complement the function (default = true)

  BaseIF* m_impFunc; // implicit function to complement

private:
  ComplementIF()
  {
    MayDay::Abort("ComplementIF uses strong construction");
  }

  void operator=(const ComplementIF& a_inputIF)
  {
    MayDay::Abort("ComplementIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
