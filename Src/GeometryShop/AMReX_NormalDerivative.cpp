
#include "AMReX_NormalDerivative.H"

using amrex::BASISV;
// Null constructor
NormalDerivative<BL_SPACEDIM>::NormalDerivative()
{
}

// Destructor
NormalDerivative<BL_SPACEDIM>::~NormalDerivative()
{
}

// Evaluate derivatives of the normal of a BaseIF subclass
Real NormalDerivative<BL_SPACEDIM>::evaluate(const IvDim  & a_multiIndex,
                                           const int    & a_direction,
                                           const RvDim  & a_point,
                                           const BaseIF & a_implicitFunction)
{
  // Create an IFSlicer from the BaseIF, call "evaluate" with this, and return
  // the result
  IFSlicer<BL_SPACEDIM> ifSlicer(a_implicitFunction);

  return evaluate(a_multiIndex,a_direction,a_point,&ifSlicer);
}

// Evaluate derivatives of the normal of an IFSlicer class
Real NormalDerivative<BL_SPACEDIM>::evaluate(const IvDim               & a_multiIndex,
                                           const int                 & a_direction,
                                           const RvDim               & a_point,
                                           const IFSlicer<BL_SPACEDIM> * a_implicitFunction)
{
  // This is the numerator of the "a_direction" component of the normal, i.e.
  // the "a_direction" component of the gradient of the function.
  DerivativeProduct gradientComponent;
  gradientComponent[BASISV_TM<int,BL_SPACEDIM>(a_direction)] = 1;

  // Represent the "a_direction" component of the normal as the "a_direction"
  // component of gradient of the function, "gradientComponent", divided by
  // the magnitude of the gradient.
  int magnitudeOfGradientPower = 1;
  PartialDerivativeTerm normalComponent(gradientComponent,
                                        magnitudeOfGradientPower);

  // Compute and store the magnitude of the gradient of the function.
  m_magnitudeOfGradient = 0.0;

  for (int idir = 0; idir < BL_SPACEDIM; idir++)
  {
    Real firstPartial = a_implicitFunction->value(BASISV(idir),RealVect(D_DECL(a_point[0], a_point[1], a_point[2])));
    m_magnitudeOfGradient += firstPartial*firstPartial;
  }

  m_magnitudeOfGradient = sqrt(m_magnitudeOfGradient);

  Real value = expand(a_multiIndex,
                      normalComponent,
                      a_point,
                      a_implicitFunction);

  return value;
}

//get m_magnitudeOfGradient
Real NormalDerivative<BL_SPACEDIM>::getMagnitudeOfGradient()
{
  return m_magnitudeOfGradient;
}

// Expand and evaluate the multi-index partial derivative of a
// PartialDerivativeTerm recursively.  If the multi-index is zero (i.e., no
// further derivatives) then simply evaluate the PartialDerivateTerm at
// "a_point" using the implicit function.  If the multi-index isn't zerom,
// explicitly compute one partial derivative which is a sum of
// PartialDerivativeTerm's and call "expand" with each of these terms and a
// reduced multi-index (which will eventually be zero).  The sum the results
// and return that sum.
Real NormalDerivative<BL_SPACEDIM>::expand(const IvDim                 & a_multiIndex,
                                         const PartialDerivativeTerm & a_term,
                                         const RvDim                 & a_point,
                                         const IFSlicer<BL_SPACEDIM>   * a_implicitFunction) const
{
  Real value;

  int firstNonZero;
  for (firstNonZero = 0; firstNonZero < BL_SPACEDIM; firstNonZero++)
  {
    if (a_multiIndex[firstNonZero] != 0)
    {
      break;
    }
  }

  // No more derivatives, evaluate the current term
  if (firstNonZero == BL_SPACEDIM)
  {
    value = 1.0;

    // Evalute the needed partial derivatives and take the product of the
    // specified powers
    const DerivativeProduct& curDerivativeProduct = a_term.first;

    for (DerivativeProduct::const_iterator it=curDerivativeProduct.begin();
         it != curDerivativeProduct.end();
         ++it)
    {
      const IvDim& curMultiIndex = it->first;
      const int&   curExponent   = it->second;

      IntVect    ivind(D_DECL(curMultiIndex[0],curMultiIndex[1],curMultiIndex[2]));
      RealVect rvpoint(D_DECL(      a_point[0],      a_point[1],      a_point[2]));
      // Evaluate a single partial derivative to its power
      Real curValue = pow(a_implicitFunction->value(ivind, rvpoint), curExponent);

      value *= curValue;
    }

    if (m_magnitudeOfGradient != 0.0)
    {
      // Divide by the magnitude of the gradient including the exponent
      int curExponent = a_term.second;
      value /= pow(m_magnitudeOfGradient,curExponent);
    }
  }
  else
  {
    value = 0.0;

    // This is the (current) partial derivative we are going to apply to the
    // current product of partial derivatives
    IvDim curPartialDerivative = BASISV_TM<int,BL_SPACEDIM>(firstNonZero);

    // This is the remaining multi-index that we haven't applied
    IvDim reducedMultiIndex = a_multiIndex;
    reducedMultiIndex[firstNonZero]--;

    // Distribute the current partial derivative over the product of
    // derivatives using the product rule
    const DerivativeProduct& curDerivativeProduct = a_term.first;
    int curExponentOfMagnitudeOfGradient = a_term.second;

    // Loop through each term in the product
    for (DerivativeProduct::const_iterator it=curDerivativeProduct.begin();
         it != curDerivativeProduct.end();
         ++it)
    {
      // Get the current derivative multi-index and exponent
      const IvDim& curMultiIndex = it->first;
      const int&   curExponent   = it->second;

      // Create the next term in the product rule by copying the current
      // product and take the current partial derivative of the current
      // product term (including the exponent).
      DerivativeProduct newDerivativeProduct = curDerivativeProduct;

      // Handle the exponent of the current product term
      Real multiplier = curExponent;

      if (curExponent == 1)
      {
        // Erase the current product term if its current exponent is one
        newDerivativeProduct.erase(curMultiIndex);
      }
      else
      {
        // Otherwise, decrement the exponent by one
        newDerivativeProduct[curMultiIndex] -= 1;
      }

      // Generate the new product term
      IvDim newMultiIndex = curMultiIndex;
      newMultiIndex += curPartialDerivative;

      // Put it into the product
      newDerivativeProduct[newMultiIndex] += 1;

      // Put the new product together with magnitude of the gradient term
      PartialDerivativeTerm newTerm(newDerivativeProduct,
                                    curExponentOfMagnitudeOfGradient);

      // Evaluate this term in the product rule (recursively)
      Real curValue = multiplier * expand(reducedMultiIndex,
                                          newTerm,
                                          a_point,
                                          a_implicitFunction);

      // Add the result into the overall product rule sum
      value += curValue;
    }

    // Now handle the last term in the overall product (and product rule)
    // which is the inverse of the magnitude of the gradient with an exponent.

    // The derivative of the magnitude of the gradient results in a sum which
    // has to be handled term by term
    for (int idir = 0; idir < BL_SPACEDIM; idir++)
    {
      // Copy the current overall product
      DerivativeProduct newDerivativeProduct = curDerivativeProduct;

      // Create the two new terms in the product
      IvDim firstPartial = BASISV_TM<int,BL_SPACEDIM>(idir);
      IvDim secondPartial = firstPartial + curPartialDerivative;

      // Put them in the new overall product
      newDerivativeProduct[firstPartial] += 1;
      newDerivativeProduct[secondPartial] += 1;

      // Generate the new exponent for the magnitude of the gradient
      int newExponentOfMagnitudeOfGradient = curExponentOfMagnitudeOfGradient + 2;

      // The multiplier due to the exponent
      Real multiplier = -curExponentOfMagnitudeOfGradient;

      // Put the new product together with magnitude of the gradient term
      PartialDerivativeTerm newTerm(newDerivativeProduct,
                                    newExponentOfMagnitudeOfGradient);

      // Evaluate this term in the product rule (recursively)
      Real curValue = multiplier * expand(reducedMultiIndex,
                                          newTerm,
                                          a_point,
                                          a_implicitFunction);

      // Add the result into the overall product rule sum
      value += curValue;
    }
  }

  return value;
}


