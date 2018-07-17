#ifndef AMREX_SMOOTHABSOLUTEVALUEIF_H_
#define AMREX_SMOOTHABSOLUTEVALUEIF_H_


#include "AMReX_RealVect.H"
#include "AMReX_IntVect.H"
#include "AMReX_BaseIF.H"

namespace amrex
{
///
/**
   Functions to take a smooth absolute value of the difference between two functions
   A_e(a,b) = |a-b|

   max(f, g) = (a + b + |a - b|)/2
   We need to make abs(a-b) smooth.

   We take the convolution of the absolute value with cos^4(pi w/(2 delta))
*/
  class SmoothAbsoluteValue
  {
  public:

    ///
    /**
       Set a known f and g (from short list I have programmed) to check answer against
       1 = double ramp with 
       ramp_normal_0 = 2. 1. 0.
       ramp_normal_1 = -2. 1. 0.
       ramp_point_0 = 0. 1. 0.
       ramp_point_1 = 1. 1. 0.

       2 = double sphere with
       sphere_center_0 = 0 0 0
       sphere_center_1 = 1 0 0
       sphere_radius_0 = 0.75
       sphere_radius_1 = 0.75
    */
    static void setKnownFunction(int a_whichfunc)
      {
        s_knownFunc = a_whichfunc;
      }
  
    ///
    SmoothAbsoluteValue(const BaseIF*   a_f,
                        const BaseIF*   a_g,
                        const Real  &   a_delta)
      {
        m_f  = a_f;
        m_g  = a_g;
        m_d  = a_delta;
        m_pi = 4.*atan(1.0);
      }

    ///
    virtual ~SmoothAbsoluteValue()
      {; }

    ///
    virtual Real smoothAbsFMinusG(const  IntVect& a_deriv,
                                  const RealVect& a_point) const;

    ///
    /**
       returns -1 if w < -delta, 1 if w > delta, 0 otherwise
       reduces to regular |f-g| unless case == 0
    */
    void getWCase(int            & a_case,
                  Real           & a_wval,
                  const RealVect & a_point)const;

  protected:

    ///
    /**
       Here is the logic of this stuff.   We have three cases that reduce 
       to two cases.    w = f(x) - g(x)
       case  1: (w > delta):  ---- whole integral is above zero
       answer = abs(w)
       case -1: (w < - delta): ---- whole integral is below zero
       answer = abs(w)
       case  0: (-delta <= w <= delta)  --- have to split integral into above and below
       answer = functionAem();
    */


    virtual Real valueAem(const RealVect& a_point) const;

    ///
    virtual Real firstDerivAem(const  IntVect& a_deriv,
                               const RealVect& a_point) const;

    ///
    virtual Real secondDerivAem(const  IntVect& a_deriv,
                                const RealVect& a_point) const;

    ///
    virtual Real thirdDerivAem(const  IntVect& a_deriv,
                               const RealVect& a_point) const;

    ///
    virtual Real fourthDerivAem(const  IntVect& a_deriv,
                                const RealVect& a_point) const;

    ///just checks nan
    bool isBogus(const Real& a_number) const;

    ///if s_knownFunc is set, check against known answer
    void checkAgainstKnown(const    Real  & a_myAnswer,
                           const  IntVect & a_deriv,
                           const RealVect & a_point) const;

                              
    //the two implicit functions are owned by others
    const BaseIF*   m_f;
    const BaseIF*   m_g;
    //delta = the smoothing length scale
    Real m_d;

    //pi, you know, pi = 4atan(1)
    Real m_pi;
    ///for debugging against known functions
    static int s_knownFunc;
  private:

    SmoothAbsoluteValue(const SmoothAbsoluteValue& a_inputIF)
      {
      }

    SmoothAbsoluteValue()
      {
      }

    void operator=(const SmoothAbsoluteValue& a_inputIF)
      {
      }
  };
///
/**
   offset sphere test case
   The mac daddy of test cases
   Test case where 
   f: (x-  0  )^2 + (y-  0  )^2 + (z-  0  )^2 -  (1/2)^2;
   g: (x-(1/2))^2 + (y-(1/2))^2 + (z-(1/2))^2 -  (1/4)^2;
   sphere_center_0 = 0 0 0
   sphere_center_1 = 0.5 0.5 0.5
   sphere_radius_0 = 0.5
   sphere_radius_1 = 0.25
**/
  class OffsetSphereExact
  {
  public:
    OffsetSphereExact(const Real & a_delta,
                      const Real & a_pi)
      {
        m_d  = a_delta;
        m_pi = a_pi;
      }

    ~OffsetSphereExact()
      {
      }

    ///this one calls the other functions
    Real value(const  IntVect& a_deriv,
               const RealVect& a_point) const;

  
    ///
    Real valueAem(const RealVect& a_point) const;

  
    ///
    Real firstDerivAem(const  IntVect& a_deriv,
                       const RealVect& a_point) const;

  
    ///
    Real secondDerivAem(const  IntVect& a_deriv,
                        const RealVect& a_point) const;

  
    ///
    Real thirdDerivAem(const  IntVect& a_deriv,
                       const RealVect& a_point) const;

  
    ///
    Real fourthDerivAem(const  IntVect& a_deriv,
                        const RealVect& a_point) const;
  private:
    void getXYIXIY(int & a_ix, int& a_iy, 
                   Real& a_x, Real& a_y,
                   const int &   a_xderivVal,
                   const  IntVect& a_deriv,
                   const RealVect& a_point) const;

    OffsetSphereExact()
      {
      }

    Real m_d;
    Real m_pi;

  };


///
/**
   Double sphere test case
   Test case where 
   f(x,y,z) = x^2 + y^2 + z^2 -  0.75^2
   g(x,y,z) = (x-1)^2 + y^2 + z^2 - 0.75^2
   g(x,y,z) = (-2*(x-1) + (y - 1))
   f(x,y.z) - g(x,y,z) = x^2 - (x-1)^2
   sphere_center_0 = 0 0 0
   sphere_center_1 = 1 0 0
   sphere_radius_0 = 0.75
   sphere_radius_1 = 0.75
**/
  class DoubleSphereExact
  {
  public:
    DoubleSphereExact(const Real & a_delta,
                      const Real & a_pi)
      {
        m_d  = a_delta;
        m_pi = a_pi;
      }

    ~DoubleSphereExact()
      {
      }

    ///this one calls the other functions
    Real value(const  IntVect& a_deriv,
               const RealVect& a_point) const;

  
    ///
    Real valueAem(const RealVect& a_point) const;

  
    ///
    Real firstDerivAem(const  IntVect& a_deriv,
                       const RealVect& a_point) const;

  
    ///
    Real secondDerivAem(const  IntVect& a_deriv,
                        const RealVect& a_point) const;

  
    ///
    Real thirdDerivAem(const  IntVect& a_deriv,
                       const RealVect& a_point) const;

  
    ///
    Real fourthDerivAem(const  IntVect& a_deriv,
                        const RealVect& a_point) const;
  private:
    DoubleSphereExact()
      {
      }

    Real m_d;
    Real m_pi;

  };

///
/**
   Double sphere test case
   Test case where 
   f(x,y) = 2*(x-0) + (y - 1)
   g(x,y) =-2*(x-1) + (y - 1)
   in the input file:
   ramp_normal_0 = 2. 1. 0.
   ramp_normal_1 = -2. 1. 0.
   ramp_point_0 = 0. 1. 0.
   ramp_point_1 = 1. 1. 0.
**/
  class DoubleRampExact
  {
  public:
    DoubleRampExact(const Real & a_delta,
                    const Real & a_pi)
      {
        m_d  = a_delta;
        m_pi = a_pi;
      }

    ~DoubleRampExact()
      {
      }

    ///this one calls the other functions
    Real value(const  IntVect& a_deriv,
               const RealVect& a_point) const;

  
    ///
    Real valueAem(const RealVect& a_point) const;

  
    ///
    Real firstDerivAem(const  IntVect& a_deriv,
                       const RealVect& a_point) const;

  
    ///
    Real secondDerivAem(const  IntVect& a_deriv,
                        const RealVect& a_point) const;

  
    ///
    Real thirdDerivAem(const  IntVect& a_deriv,
                       const RealVect& a_point) const;

  
    ///
    Real fourthDerivAem(const  IntVect& a_deriv,
                        const RealVect& a_point) const;
  private:

    DoubleRampExact()
      {
      }

    Real m_d;
    Real m_pi;

  };
}
#endif
