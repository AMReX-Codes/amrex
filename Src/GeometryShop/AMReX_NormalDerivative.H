
#ifndef AMREX_NORMALDERIVATIVE_H_
#define AMREX_NORMALDERIVATIVE_H_

#include "AMReX_Notation.H"
#include "AMReX_IFSlicer.H"

using std::map;
using std::pair;
/// This computes the derivatives of the normal of a sliced implicit function
/**
   This computes the derivatives of the normal of a sliced implicit function
*/
template <int dim> class NormalDerivative
{
public:
  typedef IndexTM<int,dim>  IvDim;
  typedef IndexTM<Real,dim> RvDim;

  // This represents a product of partial derivatives of a function.  Each
  // unique partial derivative is represented by a "IvDim" which specifies the
  // number of partial derivatives taken in each coordinate direction.  Each
  // entry in the map is one such term and the map's value for the entry is
  // the exponent of the entry.  The product is simply the product of all the
  // entries each with this exponent.
  typedef map<IvDim,int,LexLT<IvDim> > DerivativeProduct;

  // This represents the product of partial derivatives of a function (given
  // by the first term) times another function to a power (given by the second
  // term).  The second function is the reciprocal of the magnitude of the
  // gradient of the original function.  This second function (to various
  // powers) appears in the normal and all it's derivatives.  Further, the
  // normal and all its derivatives can be written as a sum of terms of this
  // form.
  typedef pair<DerivativeProduct,int> PartialDerivativeTerm;

  /// Null constructor
  /**
     Null constructor
  */
  NormalDerivative();

  /// Destructor
  /**
     Destructor
  */
  virtual ~NormalDerivative();

  /// Evaluate derivatives of the normal of an IFSlicer class
  /**
     Evaluate a derivative of the normal, a_multiIndex specifies how many
     derivatives to take in each coordinate direction, a_direction specifies
     which component of the normal use, a_point is the point in space to
     evaluate the derivative, and a_ifSlicer is a sliced function whose
     gradient is the normal.
  */
  virtual Real evaluate(const IvDim         & a_multiIndex,
                        const int           & a_direction,
                        const RvDim         & a_point,
                        const IFSlicer<dim> * a_ifSlicer);

  Real getMagnitudeOfGradient();

protected:
  // Expand and evaluate the multi-index partial derivative of a
  // PartialDerivativeTerm recursively.  If the multi-index is zero (i.e., no
  // further derivatives) then simply evaluate the PartialDerivateTerm at
  // "a_point" using the implicit function.  If the multi-index isn't zerom,
  // explicitly compute one partial derivative which is a sum of
  // PartialDerivativeTerm's and call "expand" with each of these terms and a
  // reduced multi-index (which will eventually be zero).  The sum the results
  // and return that sum.
  Real expand(const IvDim                 & a_multiIndex,
              const PartialDerivativeTerm & a_term,
              const RvDim                 & a_point,
              const IFSlicer<dim>         * a_ifSlicer) const;

  // The value of the magnitude of the gradient of the function at "a_point"
  Real m_magnitudeOfGradient;


private:
  NormalDerivative(const NormalDerivative& a_input)
    {
    }

  void operator=(const NormalDerivative& a_input)
    {
    }
};

/// This computes the derivatives of the normal of a sliced implicit function
/**
   This computes the derivatives of the normal of a sliced implicit function.
   Specialization so an additional method can be defined.
*/
template <> class NormalDerivative<BL_SPACEDIM>
{
public:
  typedef IndexTM<int,BL_SPACEDIM>  IvDim;
  typedef IndexTM<Real,BL_SPACEDIM> RvDim;

  // This represents a product of partial derivatives of a function.  Each
  // unique partial derivative is represented by a "IvDim" which specifies the
  // number of partial derivatives taken in each coordinate direction.  Each
  // entry in the map is one such term and the map's value for the entry is
  // the exponent of the entry.  The product is simply the product of all the
  // entries each with this exponent.
  typedef map<IvDim,int,LexLT<IvDim> > DerivativeProduct;

  // This represents the product of partial derivatives of a function (given
  // by the first term) times another function to a power (given by the second
  // term).  The second function is the reciprocal of the magnitude of the
  // gradient of the original function.  This second function (to various
  // powers) appears in the normal and all it's derivatives.  Further, the
  // normal and all its derivatives can be written as a sum of terms of this
  // form.
  typedef pair<DerivativeProduct,int> PartialDerivativeTerm;

  /// Null constructor
  /**
     Null constructor
  */
  NormalDerivative();

  /// Destructor
  /**
     Destructor
  */
  virtual ~NormalDerivative();

  /// Evaluate derivatives of the normal of a BaseIF subclass
  /**
     Evaluate a derivative of the normal, a_multiIndex specifies how many
     derivatives to take in each coordinate direction, a_direction specifies
     which component of the normal use, a_point is the point in space to
     evaluate the derivative, and a_impFunc is a function whose gradient
     is the normal.
  */
  virtual Real evaluate(const IvDim  & a_multiIndex,
                        const int    & a_direction,
                        const RvDim  & a_point,
                        const BaseIF & a_impFunc);

  /// Evaluate derivatives of the normal of an IFSlicer class
  /**
     Evaluate a derivative of the normal, a_multiIndex specifies how many
     derivatives to take in each coordinate direction, a_direction specifies
     which component of the normal use, a_point is the point in space to
     evaluate the derivative, and a_ifSlicer is a sliced function whose
     gradient is the normal.
  */
  virtual Real evaluate(const IvDim               & a_multiIndex,
                        const int                 & a_direction,
                        const RvDim               & a_point,
                        const IFSlicer<BL_SPACEDIM> * a_ifSlicer);

  Real getMagnitudeOfGradient();
protected:
  // Expand and evaluate the multi-index partial derivative of a
  // PartialDerivativeTerm recursively.  If the multi-index is zero (i.e., no
  // further derivatives) then simply evaluate the PartialDerivateTerm at
  // "a_point" using the implicit function.  If the multi-index isn't zerom,
  // explicitly compute one partial derivative which is a sum of
  // PartialDerivativeTerm's and call "expand" with each of these terms and a
  // reduced multi-index (which will eventually be zero).  The sum the results
  // and return that sum.
  Real expand(const IvDim                 & a_multiIndex,
              const PartialDerivativeTerm & a_term,
              const RvDim                 & a_point,
              const IFSlicer<BL_SPACEDIM>   * a_ifSlicer) const;

  // The value of the magnitude of the gradient of the function at "a_point"
  Real m_magnitudeOfGradient;


private:
  NormalDerivative(const NormalDerivative& a_input)
    {
    }

  void operator=(const NormalDerivative& a_input)
    {
    }
};


#endif
