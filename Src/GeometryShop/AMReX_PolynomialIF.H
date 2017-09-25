/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */


#ifndef _POLYNOMIALIF_H_
#define _POLYNOMIALIF_H_

#include "AMReX_RealVect.H"
#include "AMReX_IntVect.H"
#include "AMReX_BaseIF.H"

namespace amrex
{


///
/**
   This represents one term in a general polynomial
*/
  class PolyTerm
  {
  public:

    ///
    /**
       Coefficient of this polynomial term
    */
    Real coef;

    ///
    /**
       Powers of this polynomial term
    */
    IntVect powers;
  };

///
/**
   This implicit function specifies a general polynomial function.
*/
  class PolynomialIF: public BaseIF
  {
  public:
    ///
    /**
       Constructor specifying the term composing the polynomial (a_polynomial)
       and whether the domain is on the inside (a_inside), i.e. where the
       function is negative.
    */
    PolynomialIF(const Array<PolyTerm>& a_polynomial,
                 const bool&             a_inside);

    /// Destructor
    virtual ~PolynomialIF()
      {
      }

    ///
    /**
       Return the value of the function at a_point using a_polynomial
    */
    virtual Real value(const RealVect         & a_point,
                       const Array<PolyTerm> & a_polynomial) const;

    ///
    /**
       Return the value of the function at a_point.
    */
    virtual Real value(const RealVect& a_point) const;


    virtual BaseIF* newImplicitFunction() const;

  protected:
    Array<PolyTerm> m_polynomial;  // polynomial (term by term)
    bool             m_inside;      // inside flag

  private:
    PolynomialIF();

  };

}
#endif
