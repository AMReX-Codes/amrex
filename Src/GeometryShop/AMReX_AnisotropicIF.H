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

#ifndef _AnisotropicIF_H_
#define _AnisotropicIF_H_

#include "AMReX_RealVect.H"
#include "AMReX_BaseIF.H"


namespace amrex
{
  ///
  /**
     This adapts an implicit function so it can work with anisotropic dx
  */
  class AnisotropicIF: public BaseIF
  {
  public:
    ///
    /**
    */
    AnisotropicIF(const shared_ptr<BaseIF>& a_baseif,
                  const RealVect          & a_dxVector)
      {
        m_baseif   =  a_baseif;
        m_dxVector = a_dxVector;
      }

    /// Destructor
    virtual ~AnisotropicIF()
    {
    }

    ///
    /**
       Return the value of the function at a_point.  When delineating a domain,
       the level set value=0 represents the boundary and value<0 is inside the
       fluid.
    */
    virtual Real value(const RealVect& a_point) const
      {
        RealVect scaledPt;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          Real ratio = m_dxVector[idir]/m_dxVector[0];
          scaledPt[idir] = ratio*a_point[idir];
        }
        Real retval = m_baseif->value(scaledPt);
        return retval;
      }


  BaseIF* 
  newImplicitFunction() const
  {
    AnisotropicIF* newPtr = new AnisotropicIF(m_baseif,
                                              m_dxVector);


    return static_cast<BaseIF*>(newPtr);
  }
  protected:
    RealVect           m_dxVector;
    shared_ptr<BaseIF> m_baseif;
  };
}

#endif
