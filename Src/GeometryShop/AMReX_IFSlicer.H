
#ifndef AMREX_IFSLICER_H_
#define AMREX_IFSLICER_H_

#include <map>

#include "AMReX_Vector.H"
#include "AMReX_REAL.H"
#include "AMReX_IndexTM.H"
#include "AMReX_Notation.H"
#include "AMReX_BaseIF.H"

using amrex::BaseIF;
using amrex::IntVect;
using amrex::RealVect;

using std::endl;

template <int dim> class IFSlicer
{
public:

  /// null constructor
  /**
     null constructor
  */
  IFSlicer();

  /// copy constructor
  /**
     copy constructor
  */
  IFSlicer(const IFSlicer<dim> & a_IFSlicer);

  /// makes a reduced IFSlicer
  /**
     makes a reduced IFSlicer
  */
  IFSlicer(IFSlicer<dim+1> const * a_slicer,
           const int             & a_fixedComp,
           const Real            & a_fixedValue);

  void setParams(IFSlicer<dim+1> const * a_slicer,
                 const int             & a_fixedComp,
                 const Real            & a_fixedValue);

  /// destructor
  /**
     destructor
  */
  virtual ~IFSlicer();

  /// Return the partial derivative evaluated at a_point.
  /**
     Represent the partial derivative as a list of integers.
     Return the partial derivative evaluated at a_point.
  */
  virtual Real value(const IndexTM<int, dim> & a_partialDerivative,
                     const IndexTM<Real,dim> & a_point)const;


  /// output function
  /**
     output function
  */
  void print(ostream& out) const;

protected:
  // The higher dimensional IFSlicer being further sliced
  IFSlicer<dim + 1> const * m_IFSlicer;

  // The component direction being held fixed and its value
  int  m_fixedComp;
  Real m_fixedValue;

private:
  //
  // These functions are made private to prevent use of the default
  // functions provided by the C++ compiler.
  //
  IFSlicer<dim>& operator=(const IFSlicer<dim>&);

};

template<> class IFSlicer<BL_SPACEDIM>
{
public:
  /// null constructor
  /**
     null constructor
  */
  IFSlicer();

  /// copy constructor
  /**
     copy constructor
  */
  IFSlicer(const IFSlicer<BL_SPACEDIM> & a_IFSlicer);

  /// if dim == BL_SPACEDIM, construct from BaseIF
  /**
     if dim == BL_SPACEDIM, construct from BaseIF
  */
  IFSlicer(const BaseIF & a_implicitFunction);

  /// destructor
  /**
     destructor
  */
  virtual ~IFSlicer();

  /// Return the partial derivative evaluated at a_point.
  /**
     Represent the partial derivative as a list of integers.
     Return the partial derivative evaluated at a_point.
  */
  virtual Real value(const IndexTM<int, BL_SPACEDIM> & a_partialDerivative,
                     const IndexTM<Real,BL_SPACEDIM> & a_point)const
    {
      RealVect rvpt;
      IntVect  ivdr;
      for(int idir = 0; idir < BL_SPACEDIM; idir++)
      {
        rvpt[idir] = a_point[idir];
        ivdr[idir] = a_partialDerivative[idir];
      }
      return value(ivdr, rvpt);
    }

  virtual Real value(const IntVect  & a_partialDerivative,
                     const RealVect & a_point) const;


  void print(ostream& a_out) const;

protected:
  // The implicit function being sliced.  At BL_SPACEDIM this is just a
  // container class to start the slicing process.
  BaseIF* m_implicitFunction;

private:
  //
  // These functions are made private to prevent use of the default
  // functions provided by the C++ compiler.
  //
  IFSlicer<BL_SPACEDIM>& operator=(const IFSlicer<BL_SPACEDIM>&);
};


#include "AMReX_IFSlicerImplem.H"

#endif
