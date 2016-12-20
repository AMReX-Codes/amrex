#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _UNIONIF_H_
#define _UNIONIF_H_

#include "MayDay.H"
#include "RealVect.H"
#include "Vector.H"

#include "BaseIF.H"

#include "NamespaceHeader.H"

///
/**
    This implicit function is the minimum of two (of more) implicit
    functions with which it is constructed.  This results in an implicit
    function whose interior is the union of the interiors of the given
    implicit functions.
 */
class UnionIF: public BaseIF
{
public:
  ///
  /**
      Constructor specifying two implicit functions as inputs.
   */
  UnionIF(const BaseIF& a_impFunc1,
          const BaseIF& a_impFunc2);

  ///
  /**
      Constructor specifying any number of implicit functions as inputs.
   */
  UnionIF(const Vector<BaseIF*>& a_impFuncs);

  /// Copy constructor
  UnionIF(const UnionIF& a_inputIF);

  /// Destructor
  virtual ~UnionIF();

  ///
  /**
      Return the value of the function at a_point.
   */
  virtual Real value(const RealVect& a_point) const;

  virtual Real value(const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual Real value(const IndexTM<int,GLOBALDIM> & a_partialDerivative,
                     const IndexTM<Real,GLOBALDIM>& a_point) const;

  virtual BaseIF* newImplicitFunction() const;

  void findClosest(const IndexTM<Real,GLOBALDIM>& a_point,int& closestIF) const;

  virtual bool fastIntersection(const RealVect& a_low,
                                const RealVect& a_high) const;

  virtual GeometryService::InOut InsideOutside(const RealVect& a_low,
                                               const RealVect& a_high) const ;

  ///
  /**
     Pass this call onto the IFs contained in this IF class.
  */
  virtual void boxLayoutChanged(const DisjointBoxLayout & a_newBoxLayout,
                                const RealVect          & a_dx)
  {
    for (int i = 0; i < m_numFuncs; i++)
    {
      m_impFuncs[i]->boxLayoutChanged(a_newBoxLayout,a_dx);
    }
  }

protected:
  int             m_numFuncs; // number of implicit functions
  Vector<BaseIF*> m_impFuncs; // implicit functions

private:
  UnionIF()
  {
    MayDay::Abort("UnionIF uses strong construction");
  }

  void operator=(const UnionIF& a_inputIF)
  {
    MayDay::Abort("UnionIF doesn't allow assignment");
  }
};

#include "NamespaceFooter.H"
#endif
