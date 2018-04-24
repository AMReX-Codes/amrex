
#ifndef AMREX_SMOOTHINTERSECTION_H_
#define AMREX_SMOOTHINTERSECTION_H_

#include "AMReX_RealVect.H"
#include "AMReX_Vector.H"
#include "AMReX_BaseIF.H"

namespace amrex
{
///
/**
   This implicit function is the maximum of two (of more) implicit
   functions with which it is constructed.  This results in an implicit
   function whose interior is the intersection of the interiors of the given
   implicit functions.
   For the maximum, it uses the smooth absolute value Ae(a,b) \approx |a-b|
   max(a,b) = 0.5*(a+b + | a-b |)
*/
  class SmoothIntersection: public BaseIF
  {
  public:

    ///
    /**
       Constructor specifying any number of implicit functions as inputs.
    */
    SmoothIntersection(const Vector<BaseIF*>& a_impFuncs,
                       const Real           & a_delta);

    /// Destructor
    virtual ~SmoothIntersection();

    ///
    /**
       Return the value of the function at a_point.
    */
    virtual Real value(const RealVect& a_point) const;

    ///
    virtual Real derivative(const  IntVect& a_deriv,
                            const RealVect& a_point) const;


    virtual BaseIF* newImplicitFunction() const;

    ///max(a,b) = 0.5*(a+b + | a-b |)
    /**
       loops through all the functions
    */
    virtual Real smoothMax(const  IntVect & a_deriv,
                           const RealVect & a_point,
                           const int      & a_closestIF, 
                           const int      & a_nextClosestIF
      ) const;

    void findClosest(int            & a_closestIF, 
                     int            & a_nextClosestIF,
                     int            & a_numWithinDelta,
                     const RealVect & a_point) const;

  protected:
    int             m_numFuncs; // number of implicit functions
    Vector<BaseIF*> m_impFuncs; // implicit functions
    Real            m_delta   ; //smoothing length scale


    //debug stuff
    IntVect  m_ivDebug;
    Real     m_dxDebug;
    RealVect  m_rvDebug;

  

  private:
    SmoothIntersection()
      {
      }

    SmoothIntersection(const SmoothIntersection& a_inputIF)
      {
      }

    void operator=(const SmoothIntersection& a_inputIF)
      {
      }
  };
}
#endif
