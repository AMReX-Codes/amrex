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

#include "AMReX_EllipsoidIF.H"

namespace amrex
{

  EllipsoidIF::EllipsoidIF(const RealVect& a_radii,
                           const RealVect& a_center,
                           const bool&     a_inside)
  {
    // Remember the parameters
    m_radii  = a_radii;
    m_center = a_center;
    m_inside = a_inside;

    // Precompute the radii squared
    m_radii2  = m_radii;
    m_radii2 *= m_radii;
  }

  EllipsoidIF::EllipsoidIF(const EllipsoidIF& a_inputIF)
  {
    // Remember the parameters
    m_radii  = a_inputIF.m_radii;
    m_center = a_inputIF.m_center;
    m_inside = a_inputIF.m_inside;

    // Precompute the radii squared
    m_radii2  = m_radii;
    m_radii2 *= m_radii;
  }

  EllipsoidIF::~EllipsoidIF()
  {
  }

  Real EllipsoidIF::value(const RealVect& a_point) const
  {
    Real retval;

    // Compute the equation of the ellipsoid
    Real sum;

    sum = 0.0;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      Real cur;
      cur = a_point[idir] - m_center[idir];

      sum += cur*cur / m_radii2[idir];
    }

    // The sum should be 1.0 on the surface of the ellipsoid
    retval = sum - 1.0;

    // Change the sign to change inside to outside
    if (!m_inside)
    {
      retval = -retval;
    }

    return retval;
  }

  BaseIF* EllipsoidIF::newImplicitFunction() const
  {
    EllipsoidIF* ellipsoidPtr = new EllipsoidIF(m_radii,
                                                m_center,
                                                m_inside);

    return static_cast<BaseIF*>(ellipsoidPtr);
  }
  
}

