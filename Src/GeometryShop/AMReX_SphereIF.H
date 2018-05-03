#ifndef _SPHEREIF_H_
#define _SPHEREIF_H_

#include "AMReX_RealVect.H"
#include "AMReX_BaseIF.H"


namespace amrex
{
  ///
  /**
     This implicit function specifies a sphere.
  */
  class SphereIF: public BaseIF
  {
  public:
    ///
    /**
       Constructor specifying sphere radius (a_radius), center (a_center), and
       whether the domain is on the inside (a_inside).
    */
    SphereIF(const Real&     a_radius,
             const RealVect& a_center,
             const bool&     a_inside);

    /// Destructor
    virtual ~SphereIF()
    {;}


    ///
    /**
       Return the value of the function at a_point.  When delineating a domain,
       the level set value=0 represents the boundary and value<0 is inside the
       fluid.
    */
    virtual Real value(const RealVect& a_point) const;


    ///return the partial derivative at the point
    virtual Real derivative(const  IntVect& a_deriv,
                            const RealVect& a_point) const;

    ///
    /**
       Return a newly allocated derived class.  The responsibility
       for deleting the memory is left to the calling function.
    */
    virtual BaseIF* newImplicitFunction() const;


  protected:
    Real              m_radius;    // radius
    RealVect          m_center;    // center
    bool              m_inside;    // inside flag
    Real              m_radius2;   // precomputed radius squared
  private:
    SphereIF();

    void operator=(const SphereIF& a_inputIF)
    {
      //MayDay::Abort("SphereIF doesn't allow assignment");
    }
  };
}//namespace amrex

#endif
