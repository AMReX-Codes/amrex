#ifndef AMREX_UNIONIF_H_
#define AMREX_UNIONIF_H_

#include "AMReX_RealVect.H"
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include "AMReX_BaseIF.H"

namespace amrex
{

  ///
  /**
     This implicit function is the minimum of two (of more) implicit
     functions with which it is constructed.  This results in an implicit
     function whose interior is the union of the interiors of the given
     implicit functions.
  */
  class UnionIF: public BaseIF
  {
  public:

    ///
    UnionIF(const Vector<BaseIF*>& a_impFuncs);
   
   
    /// Destructor
    virtual ~UnionIF();
   
    ///
    /**
       Return the value of the function at a_point.
    */
    virtual Real value(const RealVect& a_point) const;
   
   
    virtual BaseIF* newImplicitFunction() const;

   
    virtual Real derivative(const  IntVect& a_deriv,
                            const RealVect& a_point) const;
   
  protected:
    int             m_numFuncs; // number of implicit functions
    Vector<BaseIF*> m_impFuncs; // implicit functions
   
  private:
    UnionIF();
    /// Copy constructor
    UnionIF(const UnionIF& a_inputIF);
   
   
    void operator=(const UnionIF& a_inputIF);
   
  };
}
#endif
