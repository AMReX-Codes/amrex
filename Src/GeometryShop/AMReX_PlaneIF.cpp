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


#include "AMReX_PlaneIF.H"


namespace amrex
{
  PlaneIF::
  PlaneIF(const RealVect& a_normal,
          const RealVect& a_point,
          const bool&     a_inside)
      : m_normal(a_normal),
        m_point (a_point),
        m_inside(a_inside)
  {}

  Real 
  PlaneIF::
  value(const RealVect & a_point) const
  {
    Real retval = 0;
    for (int idir = 0 ; idir < SpaceDim; idir++)
      {
        retval += (a_point[idir] - m_point[idir]) * m_normal[idir];
      }
    if(m_inside)
      {
        retval = -retval;
      }
    return retval;
  }

  BaseIF* 
  PlaneIF::
  newImplicitFunction() const
  {
    PlaneIF* newPtr = new PlaneIF(m_normal,
                                     m_point,
                                     m_inside);


    return static_cast<BaseIF*>(newPtr);
  }

}
