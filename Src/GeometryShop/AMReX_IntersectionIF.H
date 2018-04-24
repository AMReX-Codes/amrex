#ifndef AMREX_INTERSECTIONIF_H_
#define AMREX_INTERSECTIONIF_H_

#include "AMReX_RealVect.H"
#include "AMReX_BaseIF.H"
#include <AMReX_Array.H>
#include <AMReX_Vector.H>

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
    IntersectionIF(const Vector<BaseIF*>& a_impFuncs);

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
    Vector<BaseIF*> m_impFuncs; // implicit functions

  private:
    IntersectionIF(const IntersectionIF& a_inputIF);
    IntersectionIF();
    void operator=(const IntersectionIF& a_inputIF);

  };

}

#endif
