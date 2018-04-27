#ifndef AMREX_EllipsoidIF_H_
#define AMREX_EllipsoidIF_H_

#include "AMReX_RealVect.H"
#include "AMReX_BaseIF.H"


namespace amrex
{
///
/**
   This implicit function specifies an ellipsoid.
*/
  class EllipsoidIF: public BaseIF
  {
  public:
    ///
    /**
       Constructor specifying ellipsoid "radii" (a_radii), center (a_center),
       and whether the domain is on the inside (a_inside).
    */
    EllipsoidIF(const RealVect& a_radii,
                const RealVect& a_center,
                const bool&     a_inside);

    /// Copy constructor
    EllipsoidIF(const EllipsoidIF& a_inputIF);

    /// Destructor
    virtual ~EllipsoidIF();

    ///
    /**
       Return the value of the function at a_point.
    */
    virtual Real value(const RealVect& a_point) const;

    virtual BaseIF* newImplicitFunction() const;

    ///return the partial derivative at the point
    virtual Real derivative(const  IntVect& a_deriv,
                            const RealVect& a_point) const;
  protected:
    RealVect m_radii;     // radii
    RealVect m_center;    // center
    bool     m_inside;    // inside flag

    RealVect m_radii2;    // precomputed radii squared

  private:
    EllipsoidIF()
      {
      }

    void operator=(const EllipsoidIF& a_inputIF)
      {
      }
  };
}//namespace amrex

#endif
