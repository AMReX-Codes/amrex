
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

#include "AMReX_IntersectionIF.H"

namespace amrex
{
  IntersectionIF::~IntersectionIF()
  {
    // Delete all the copies
    for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
    {
      if (m_impFuncs[ifunc] != NULL)
      {
        delete m_impFuncs[ifunc];
      }
    }
  }

  Real IntersectionIF::value(const RealVect& a_point) const
  {

    // Maximum of the implicit functions values
    Real retval;

    retval = -1.0;

    // Find the maximum value and return it
    if (m_numFuncs > 0)
    {
      retval = m_impFuncs[0]->value(a_point);

      for (int ifunc = 1; ifunc < m_numFuncs; ifunc++)
      {
        Real cur;

        cur = m_impFuncs[ifunc]->value(a_point);
        if (cur > retval)
        {
          retval = cur;
        }
      }
    }

    return retval;
  }


  BaseIF* IntersectionIF::newImplicitFunction() const
  {
    IntersectionIF* intersectionPtr = new IntersectionIF(m_impFuncs);

    return static_cast<BaseIF*>(intersectionPtr);
  }

  IntersectionIF::IntersectionIF(const Array<BaseIF *>& a_impFuncs)
  {
    // Number of implicit function in intersection
    m_numFuncs = a_impFuncs.size();

    // Vector of implicit function pointers
    m_impFuncs.resize(m_numFuncs);

    // Make copies of the implicit functions

    for (int ifunc = 0; ifunc < m_numFuncs; ifunc++)
    {
      if (a_impFuncs[ifunc] == NULL)
      {
        m_impFuncs[ifunc] = NULL;
      }
      else
      {
        m_impFuncs[ifunc] = a_impFuncs[ifunc]->newImplicitFunction();
      }
    }
  }
}
