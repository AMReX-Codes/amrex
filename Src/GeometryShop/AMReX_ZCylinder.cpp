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

#include "AMReX_ZCylinder.H"

namespace amrex
{

  ZCylinder::
  ZCylinder(const Real&     a_radius,
           const RealVect& a_center,
           const bool&     a_inside)
    
  {
    m_radius  = a_radius;
    m_radius2 = m_radius*m_radius;
    m_inside  = a_inside;
    m_center  = a_center;
  }

  Real
  ZCylinder::
  value(const RealVect& a_point) const
  {
    RealVect dist = a_point - m_center;
    Real distance2 = D_TERM(dist[0]*dist[0], + dist[1]*dist[1], + dist[2]*dist[2]);

    Real retval = distance2 - m_radius2;
    // Change the sign to change inside to outside
    if (!m_inside)
      {
        retval = -retval;
      }

    return retval;
  }

  BaseIF* 
  ZCylinder::
  newImplicitFunction() const
  {
    ZCylinder* spherePtr = new ZCylinder(m_radius,
                                         m_center,
                                         m_inside);

    return static_cast<BaseIF*>(spherePtr);
  }
  
}

