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

#ifndef _AnisotropicPLANEIF_H_
#define _AnisotropicPLANEIF_H_

#include "AMReX_RealVect.H"
#include "AMReX_BaseIF.H"
#include "AMReX_PlaneIF.H"

namespace amrex
{
  ///
  /**
     This implicit function specifies a plane.
  */
  class AnisotropicDxPlaneIF: public PlaneIF
  {
  public:
    ///
    /**
       Constructor specifying plane normal (a_normal), a point on the plance
       (a_point), and whether the domain is the same side as the normal
       (a_inside).
    */
    AnisotropicDxPlaneIF(const RealVect& a_normal,
                       const RealVect& a_point,
                       const bool&     a_inside,
                       const RealVect& a_dxVector)
      :PlaneIF(a_normal, a_point, a_inside)
      {
        m_dxVector = a_dxVector;
      }

    /// Destructor
    virtual ~AnisotropicDxPlaneIF()
    {
    }

    ///
    /**
       Return the value of the function at a_point.  When delineating a domain,
       the level set value=0 represents the boundary and value<0 is inside the
       fluid.
    */
    virtual Real value(const RealVect& a_point) const override
      {
        RealVect scaledPt;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          Real ratio = m_dxVector[idir]/m_dxVector[0];
//          Real ratio = m_dxVector[0]/m_dxVector[idir];
          scaledPt[idir] = ratio*a_point[idir];
        }
        Real retval = PlaneIF::value(scaledPt);
      }


  protected:
    RealVect m_dxVector;
  };
}

#endif
