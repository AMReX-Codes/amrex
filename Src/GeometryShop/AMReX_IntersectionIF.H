
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

#ifndef _INTERSECTIONIF_H_
#define _INTERSECTIONIF_H_

#include "AMReX_RealVect.H"
#include "AMReX_BaseIF.H"
#include "AMReX_Array.H"

namespace amrex
{

///
/**
   This implicit function is the maximum of two (of more) implicit
   functions with which it is constructed.  This results in an implicit
   function whose interior is the intersection of the interiors of the given
   implicit functions.
*/
  class IntersectionIF: public BaseIF
  {
  public:
    ///
    /**
       Constructor specifying any number of implicit functions as inputs.
    */
    IntersectionIF(const Array<BaseIF*>& a_impFuncs);

    /// Copy constructor


    /// Destructor
    virtual ~IntersectionIF();

    ///
    /**
       Return the value of the function at a_point.
    */
    virtual Real value(const RealVect& a_point) const;

    virtual BaseIF* newImplicitFunction() const;


  protected:
    int             m_numFuncs; // number of implicit functions
    Array<BaseIF*> m_impFuncs; // implicit functions

  private:
    IntersectionIF(const IntersectionIF& a_inputIF);
    IntersectionIF();
    void operator=(const IntersectionIF& a_inputIF);

  };

}

#endif
