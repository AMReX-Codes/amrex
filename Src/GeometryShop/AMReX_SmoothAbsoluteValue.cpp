
#include "AMReX_SmoothAbsoluteValue.H"
#include "AMReX_Factorial.H"
#include "AMReX_EBArith.H"

///For these derivatives, let us be thankful for 
/// A) the existence of maxima
/// B) trig functions having uncomplicated derivatives
/// and most importantly
/// C) The equality of mixed derivs.

namespace amrex
{
  int SmoothAbsoluteValue::s_knownFunc = -1;

///if s_knownFunc is set, check against known answer
  void 
  SmoothAbsoluteValue::
  checkAgainstKnown(const    Real  & a_myAnswer,
                    const  IntVect & a_deriv,
                    const RealVect & a_point) const
  {
    if((s_knownFunc == 1) && (a_deriv.sum() <= 4))
    {
      Real tol = 1.0e-9;
      DoubleRampExact cramp(m_d, m_pi);
      Real rightAns = cramp.value(a_deriv, a_point);
      if(std::abs(rightAns - a_myAnswer) > tol)
      {
        pout() << "error detected at point = " << a_point << ", deriv = " << a_deriv <<  endl;
        amrex::Warning("double ramp fault");
      }
    }
    else if((s_knownFunc == 2) && (a_deriv.sum() <= 4))
    {
      Real tol = 1.0e-7;
      DoubleSphereExact cramp(m_d, m_pi);
      Real rightAns = cramp.value(a_deriv, a_point);
      if(std::abs(rightAns - a_myAnswer) > tol)
      {
        pout() << "error detected at point = " << a_point << ", deriv = " << a_deriv <<  endl;
        amrex::Warning("double sphere fault");
      }
    }
    else if((s_knownFunc == 3) && (a_deriv.sum() <= 4) && (SpaceDim == 2))
    {
      Real tol = 1.0e-7;
      OffsetSphereExact cramp(m_d, m_pi);
      Real rightAns = cramp.value(a_deriv, a_point);
      if(std::abs(rightAns - a_myAnswer) > tol)
      {
        pout() << "error detected at point = " << a_point << ", deriv = " << a_deriv <<  endl;
        amrex::Warning("offset sphere fault");
      }
    }
  }

///
  bool 
  SmoothAbsoluteValue::
  isBogus(const Real& a_number) const
  {
    bool retval = false;
    if(a_number != a_number)
      retval = true;

    return retval;
  }
  Real 
  SmoothAbsoluteValue::
  valueAem(const RealVect& a_point) const
  {

    Real wval; 
    int icase;
    getWCase(icase, wval, a_point);
    BL_ASSERT(icase == 0);
    //(%i2) Aem: integrate((4/(3*d))*y*cos(%pi*((w-y)/(2*d)))^4, y,0,w+d) - integrate((4/(3*d))*y*cos(%pi*((w-y)/(2*d)))^4, y,w-d,0);
    //
    //Is  w+d  positive, negative, or zero?
    //
    //positive;
    //(%o2) 4*((12*%pi^2*d*w+(6*%pi^2-15)*d^2)/(32*%pi^2)
    //        -(d^2*cos(2*%pi*w/d)+16*d^2*cos(%pi*w/d)-6*%pi^2*w^2)/(32*%pi^2))
    //       /(3*d)
    //       -4*((d^2*cos(2*%pi*w/d)+16*d^2*cos(%pi*w/d)-6*%pi^2*w^2)/(32*%pi^2)
    //          +(12*%pi^2*d*w+(15-6*%pi^2)*d^2)/(32*%pi^2))
    //        /(3*d)
    //(%i3) expand(%);
    //
    //(%o3) -d*cos(2*%pi*w/d)/(12*%pi^2)-4*d*cos(%pi*w/d)/(3*%pi^2)+w^2/(2*d)
    //                                  -5*d/(4*%pi^2)+d/2
    //(%i4) Aem: -d*cos(2*%pi*(f(x,y,z) - g(x,y,z))/d)/(12*%pi^2)-4*d*cos(%pi*(f(x,y,z) - g(x,y,z))/d)/(3*%pi^2)+(f(x,y,z) - g(x,y,z))^2/(2*d)
    //                                  -5*d/(4*%pi^2)+d/2;
    //
    //(%o4) -d*cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)/(12*%pi^2)
    //       -4*d*cos(%pi*(f(x,y,z)-g(x,y,z))/d)/(3*%pi^2)
    //       +(f(x,y,z)-g(x,y,z))^2/(2*d)-5*d/(4*%pi^2)+d/2
    //(%i5) expand(%);
    //
    //(%o5) -d*cos(2*%pi*g(x,y,z)/d-2*%pi*f(x,y,z)/d)/(12*%pi^2)
    //       -4*d*cos(%pi*g(x,y,z)/d-%pi*f(x,y,z)/d)/(3*%pi^2)+g(x,y,z)^2/(2*d)
    //       -f(x,y,z)*g(x,y,z)/d+f(x,y,z)^2/(2*d)-5*d/(4*%pi^2)+d/2
    //(%i47) 
    ///relevant line:
    //     -d*cos(2*%pi*w/d)/(12*%pi^2)
    //   -4*d*cos(  %pi*w/d)/( 3*%pi^2)
    //  +w^2/(2*d)-5*d/(4*%pi^2)+d/2
    //                                    
    Real retval =  
      -m_d*cos(2*m_pi*wval/m_d)/(12*POW(m_pi,2))
      -4*m_d*cos(  m_pi*wval/m_d)/(3* POW(m_pi,2)) 
      +POW(wval,2)/(2*m_d) - 5*m_d/(4*POW(m_pi,2))+m_d/2;

    if(isBogus(retval))
    {
      amrex::Error("bogosity detected aem0");
    }

    return retval;
  }
///
  Real 
  SmoothAbsoluteValue::
  firstDerivAem(const  IntVect& a_deriv,
                const RealVect& a_point) const
  {
    BL_ASSERT(a_deriv.sum() == 1);
    Real wval; 
    int icase;
    getWCase(icase, wval, a_point);
    BL_ASSERT(icase == 0);
    //values of the two functions
    Real fval = (*m_f).value(a_point);
    Real gval = (*m_g).value(a_point);

    //derivatives of the two functions
    Real dfx = (*m_f).value(a_deriv, a_point);
    Real dgx = (*m_g).value(a_deriv, a_point);

    ///(%i7) Aem_x: diff(Aem,x);                                    
    ///                                                             
    ///(%o7) sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1) -'diff(g(x,y,z),x,1))/(6*%pi)                                              
    ///   +4*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1) -'diff(g(x,y,z),x,1))/(3*%pi)                                             
    ///               +(f(x,y,z)-g(x,y,z)   )*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))/d

    Real retval = 
      sin(2*m_pi*(fval - gval)/m_d)*(dfx - dgx)/(6*m_pi)
      +4*sin(  m_pi*(fval - gval)/m_d)*(dfx - dgx)/(3*m_pi)
      +             (fval - gval     )*(dfx - dgx)/(m_d);

    if(isBogus(retval))
    {
      amrex::Error("bogosity detected aem1");
    }
    return retval;
  }
///
  Real 
  SmoothAbsoluteValue::
  secondDerivAem(const  IntVect& a_deriv,
                 const RealVect& a_point) const
  {
    BL_ASSERT(a_deriv.sum() == 2);
    Real wval; 
    int icase;
    getWCase(icase, wval, a_point);
    BL_ASSERT(icase == 0);

    //true for dxx dyy dzz, false for dxy dyz dxz
    //the two cases have different forms.
    bool doublex = (a_deriv.max() == 2);

    Real retval = 0;
    if(doublex)
    {
      int ix;
      bool found= false;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        if(a_deriv[idir] == 2)
        {
          ix = idir;
          found = true;
        }
      }
      if(!found)
      {
        amrex::Error("logic error aem2.0");
      }

      IntVect dx  = BASISV(ix);
      IntVect dxx = 2*dx;
      //derivatives
      Real  dfx  = (*m_f).value( dx    , a_point);
      Real  dfxx = (*m_f).value( dxx   , a_point);

      Real  dgx  = (*m_g).value( dx    , a_point);
      Real  dgxx = (*m_g).value( dxx   , a_point);


      Real fval = (*m_f).value(a_point);
      Real gval = (*m_g).value(a_point);

      ///(%i14) Aem_xx: diff(Aem_x,x);                                                  
      ///                                                                               
      ///(%o14) sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))/(6*%pi)                                                               
      ///    +4*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))/(3*%pi)                                                              
      ///      +cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2/(3*d)                                                                
      ///     4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2/(3*d)
      ///                +(f(x,y,z)-g(x,y,z)   )*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))/d       
      ///                                       +('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2/d                                       
      ///

      retval = 
        sin(   2*m_pi*(fval - gval)/m_d)*   (dfxx - dgxx  )/(6*m_pi)   
        +4*sin(  m_pi*(fval - gval)/m_d)*   (dfxx - dgxx  )/(3*m_pi)   
        +  cos(2*m_pi*(fval - gval)/m_d)*POW(dfx  - dgx, 2)/(3*m_d) 
        +4*cos(  m_pi*(fval - gval)/m_d)*POW(dfx  - dgx, 2)/(3*m_d) 
        +             (fval - gval)*        (dfxx - dgxx  )/m_d 
        +                                POW(dfx  - dgx, 2)/m_d;
      if(isBogus(retval))
      {
        amrex::Error("bogosity detected aem2.0");
      }
    }
    else
    {
      int ix, iy;
      bool foundx= false;
      bool foundy= false;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        if(a_deriv[idir] == 1)
        {
          if(!foundx)
          {
            ix = idir;
            foundx = true;
          }
          else
          {
            iy = idir;
            foundy = true;
          }
        }
      }
      if(!foundx || !foundy)
      {
        amrex::Error("logic error aem2.1");
      }
      IntVect dx  = BASISV(ix);
      IntVect dy  = BASISV(iy);
      IntVect dxy = dx + dy;

      Real dfx  = (*m_f).value( dx  , a_point);
      Real dfy  = (*m_f).value( dy  , a_point);
      Real dfxy = (*m_f).value( dxy , a_point);
 
      Real dgx  = (*m_g).value( dx  , a_point);
      Real dgy  = (*m_g).value( dy  , a_point);
      Real dgxy = (*m_g).value( dxy , a_point);

      Real fval = (*m_f).value(a_point);
      Real gval = (*m_g).value(a_point);

      //(%i21) Aem_xy: diff(Aem_x,y);
      //
      //(%o21)     cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d)
      //        +4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d)
      //        +                                   ('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/d
      //        +  sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))/(6*%pi)
      //        +4*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))/(3*%pi)
      //        +            (f(x,y,z)-g(x,y,z)   )*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))/d


      //(%i21) Aem_xy: diff(Aem_x,y);
      //
      retval = cos(2*m_pi*(fval - gval)/m_d)*(dfx  - dgx )*(dfy - dgy)/(3*m_d)
        +    4*cos(  m_pi*(fval - gval)/m_d)*(dfx  - dgx )*(dfy - dgy)/(3*m_d)
        +                                    (dfx  - dgx )*(dfy - dgy)/m_d
        +      sin(2*m_pi*(fval - gval)/m_d)*(dfxy - dgxy)/(6*m_pi)
        +    4*sin(  m_pi*(fval - gval)/m_d)*(dfxy - dgxy)/(3*m_pi)
        +                 (fval - gval)     *(dfxy - dgxy)/m_d;

      if(isBogus(retval))
      {
        amrex::Error("bogosity detected aem2.1");
      }
    }
    return retval;
  }
///
  Real
  SmoothAbsoluteValue::
  thirdDerivAem(const  IntVect& a_deriv,
                const RealVect& a_point) const
  {
    BL_ASSERT(a_deriv.sum() == 3);
    Real wval; 
    int icase;
    getWCase(icase, wval, a_point);
    BL_ASSERT(icase == 0);
    //these are the different forms that the third deriv can take
    bool xxx = (a_deriv.max() == 3);
    bool xxy = (a_deriv.max() == 2);
#if BL_SPACEDIM==3
    bool xyz = (a_deriv.max() == 1);
#endif
    Real retval;
    if(xxx)
    {
      int ix;
      bool found = false;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        if(a_deriv[idir] == 3)
        {
          ix = idir;
          found = true;
        }
      }
      if(!found)
      {
        amrex::Error("logic error aem3.0");
      }
      IntVect  dx   = BASISV(ix);
      IntVect  dxx  = 2*dx;
      IntVect  dxxx = 3*dx;        

      Real  dfx   = (*m_f).value( dx   , a_point);
      Real  dfxx  = (*m_f).value( dxx  , a_point);
      Real  dfxxx = (*m_f).value( dxxx , a_point);

      Real  dgx   = (*m_g).value( dx   , a_point);
      Real  dgxx  = (*m_g).value( dxx  , a_point);
      Real  dgxxx = (*m_g).value( dxxx , a_point);

      Real fval  = (*m_f).value(a_point);
      Real gval  = (*m_g).value(a_point);

      //(%i28) Aem_xxx: diff(Aem_xx,x);
      //
      //(%o28) 
      //             sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,3)-'diff(g(x,y,z),x,3))/(6*%pi)
      //        +  4*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,3)-'diff(g(x,y,z),x,3))/(3*%pi)
      //        +              (f(x,y,z)-g(x,y,z)   )*('diff(f(x,y,z),x,3)-'diff(g(x,y,z),x,3))/d
      //        +    cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))/d
      //        +  4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))/d
      //        +                                   3*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))/d
      //        -2*%pi*sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^3/(3*d^2)
      //        -4*%pi*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^3/(3*d^2)


      retval =   sin(2*m_pi*(fval - gval)/m_d)*   (dfxxx - dgxxx)/(6*m_pi)
        +      4*sin(  m_pi*(fval - gval)/m_d)*   (dfxxx - dgxxx)/(3*m_pi)
        +                   (fval - gval     )*   (dfxxx - dgxxx)/m_d
        +        cos(2*m_pi*(fval - gval)/m_d)*   (dfx   - dgx  )*(dfxx - dgxx)/m_d
        +      4*cos(  m_pi*(fval - gval)/m_d)*   (dfx   - dgx  )*(dfxx - dgxx)/m_d
        +                                    3*   (dfx   - dgx  )*(dfxx - dgxx)/m_d
        -2*m_pi*sin( 2*m_pi*(fval - gval)/m_d)*POW(dfx   - dgx,3)/(3*POW(m_d, 2))
        -4*m_pi*sin(   m_pi*(fval - gval)/m_d)*POW(dfx   - dgx,3)/(3*POW(m_d, 2));
      if(isBogus(retval))
      {
        amrex::Error("bogosity detected aemxxx");
      }
    }
    else if(xxy)
    {
      int ix, iy;
      bool foundx = false;
      bool foundy = false;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        if(a_deriv[idir] == 2)
        {
          ix = idir;
          foundx = true;
        }
        else if(a_deriv[idir] == 1)
        {
          iy = idir;
          foundy = true;
        }
      }
      if(!foundx || !foundy)
      {
        amrex::Error("logic error aem3.1");
      }

      IntVect  dx   = BASISV(ix);
      IntVect  dy   = BASISV(iy);
      IntVect  dxx  = 2*dx;
      IntVect  dxy  = dx  + dy;
      IntVect  dxxy = dxx + dy;


      Real  dfx   = (*m_f).value( dx   , a_point);
      Real  dfy   = (*m_f).value( dy   , a_point);
      Real  dfxx  = (*m_f).value( dxx  , a_point);
      Real  dfxy  = (*m_f).value( dxy  , a_point);
      Real  dfxxy = (*m_f).value( dxxy , a_point);

      Real  dgx   = (*m_g).value( dx   , a_point);
      Real  dgy   = (*m_g).value( dy   , a_point);
      Real  dgxx  = (*m_g).value( dxx  , a_point);
      Real  dgxy  = (*m_g).value( dxy  , a_point);
      Real  dgxxy = (*m_g).value( dxxy , a_point);

      Real fval  = (*m_f).value(a_point);
      Real gval  = (*m_g).value(a_point);

      //(%i38) Aem_xxy: diff(Aem_xx,y);
      //
      //(%o38) 
      //               cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d)
      //        +    4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d)
      //        +                                       ('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/d
      //        -2*%pi*sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d^2)
      //        -4*%pi*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d^2)
      //        +      sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2,y,1)-'diff(g(x,y,z),x,2,y,1))/(6*%pi)
      //        +    4*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2,y,1)-'diff(g(x,y,z),x,2,y,1))/(3*%pi)
      //        +                (f(x,y,z)-g(x,y,z)   )*('diff(f(x,y,z),x,2,y,1)-'diff(g(x,y,z),x,2,y,1))/d
      //        +    2*cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))/(3*d)
      //        +    8*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))/(3*d)
      //        +                                     2*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))/d


      retval =  cos(2*m_pi*(fval - gval)/m_d)*(   dfxx  - dgxx )*(dfy  - dgy )/(3*m_d)
        +     4*cos(  m_pi*(fval - gval)/m_d)*(   dfxx  - dgxx )*(dfy  - dgy )/(3*m_d)
        +                                     (   dfxx  - dgxx )*(dfy  - dgy )/m_d
        -2*m_pi*sin(2*m_pi*(fval - gval)/m_d)*POW(dfx   - dgx,2)*(dfy  - dgy )/(3*POW(m_d, 2))
        -4*m_pi*sin(  m_pi*(fval - gval)/m_d)*POW(dfx   - dgx,2)*(dfy  - dgy )/(3*POW(m_d, 2))
        +     2*cos(2*m_pi*(fval - gval)/m_d)*(   dfx   - dgx  )*(dfxy - dgxy)/(3*m_d)
        +     8*cos(  m_pi*(fval - gval)/m_d)*(   dfx   - dgx  )*(dfxy - dgxy)/(3*m_d)
        +                                   2*(   dfx   - dgx  )*(dfxy - dgxy)/m_d
        +       sin(2*m_pi*(fval - gval)/m_d)*(   dfxxy - dgxxy)/(6*m_pi)
        +     4*sin(  m_pi*(fval - gval)/m_d)*(   dfxxy - dgxxy)/(3*m_pi)
        +                  (fval - gval     )*(   dfxxy - dgxxy)/m_d;

      if(isBogus(retval))
      {
        amrex::Error("bogosity detected aemxxy");
      }
    }
#if BL_SPACEDIM == 3
    else if(xyz)
    {

      IntVect  dx   = BASISV(0);
      IntVect  dy   = BASISV(1);
      IntVect  dz   = BASISV(2);
      IntVect  dxy  = dx  + dy;
      IntVect  dxz  = dx  + dz;
      IntVect  dyz  = dy  + dz;
      IntVect  dxyz = dx  + dy + dz;

      Real fval  = (*m_f).value(a_point);
      Real gval  = (*m_g).value(a_point);

      Real  dfx   = (*m_f).value( dx   , a_point);
      Real  dfy   = (*m_f).value( dy   , a_point);
      Real  dfz   = (*m_f).value( dz   , a_point);
      Real  dfxy  = (*m_f).value( dxy  , a_point);
      Real  dfxz  = (*m_f).value( dxz  , a_point);
      Real  dfyz  = (*m_f).value( dyz  , a_point);
      Real  dfxyz = (*m_f).value( dxyz , a_point);

      Real  dgx   = (*m_g).value( dx   , a_point);
      Real  dgy   = (*m_g).value( dy   , a_point);
      Real  dgz   = (*m_g).value( dz   , a_point);
      Real  dgxy  = (*m_g).value( dxy  , a_point);
      Real  dgxz  = (*m_g).value( dxz  , a_point);
      Real  dgyz  = (*m_g).value( dyz  , a_point);
      Real  dgxyz = (*m_g).value( dxyz , a_point);


      //(%i35) Aem_xyz: diff(Aem_xy,z);
      //
      //(%o35) 
      //        -2*%pi*sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))*('diff(f(x,y,z),z,1)-'diff(g(x,y,z),z,1))/(3*d^2)
      //        -4*%pi*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))*('diff(f(x,y,z),z,1)-'diff(g(x,y,z),z,1))/(3*d^2)
      //        +      cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))*('diff(f(x,y,z),z,1)-'diff(g(x,y,z),z,1))/(3*d)
      //        +    4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))*('diff(f(x,y,z),z,1)-'diff(g(x,y,z),z,1))/(3*d)
      //        +                                       ('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))*('diff(f(x,y,z),z,1)-'diff(g(x,y,z),z,1))/d
      //        +      cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),y,1,z,1)-'diff(g(x,y,z),y,1,z,1))/(3*d)
      //        +    4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),y,1,z,1)-'diff(g(x,y,z),y,1,z,1))/(3*d)
      //        +                                       ('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),y,1,z,1)-'diff(g(x,y,z),y,1,z,1))/d
      //        +      cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1,z,1)-'diff(g(x,y,z),x,1,z,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d)
      //        +    4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1,z,1)-'diff(g(x,y,z),x,1,z,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d)
      //        +                                       ('diff(f(x,y,z),x,1,z,1)-'diff(g(x,y,z),x,1,z,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/d
      //        +      sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1,y,1,z,1)-'diff(g(x,y,z),x,1,y,1,z,1))/(6*%pi)
      //        +    4*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1,y,1,z,1)-'diff(g(x,y,z),x,1,y,1,z,1))/(3*%pi)
      //        +                (f(x,y,z)-g(x,y,z)   )*('diff(f(x,y,z),x,1,y,1,z,1)-'diff(g(x,y,z),x,1,y,1,z,1))/d



      retval = -2*m_pi*sin(2*m_pi*(fval - gval)/m_d)*(dfx   - dgx  )*(dfy  - dgy  )*(dfz  - dgz  )/(3*POW(m_d, 2))
        -       4*m_pi*sin(  m_pi*(fval - gval)/m_d)*(dfx   - dgx  )*(dfy  - dgy  )*(dfz  - dgz  )/(3*POW(m_d, 2))
        +              cos(2*m_pi*(fval - gval)/m_d)*(dfxy  - dgxy )*(dfz  - dgz  )/(3*m_d)
        +            4*cos(  m_pi*(fval - gval)/m_d)*(dfxy  - dgxy )*(dfz  - dgz  )/(3*m_d)
        +                                            (dfxy  - dgxy )*(dfz  - dgz  )/m_d
        +              cos(2*m_pi*(fval - gval)/m_d)*(dfx   - dgx  )*(dfyz - dgyz )/(3*m_d)
        +            4*cos(  m_pi*(fval - gval)/m_d)*(dfx   - dgx  )*(dfyz - dgyz )/(3*m_d)
        +                                            (dfx   - dgx  )*(dfyz - dgyz )/m_d
        +              cos(2*m_pi*(fval - gval)/m_d)*(dfxz  - dgxz )*(dfy  - dgy  )/(3*m_d)
        +            4*cos(  m_pi*(fval - gval)/m_d)*(dfxz  - dgxz )*(dfy  - dgy  )/(3*m_d)
        +                                            (dfxz  - dgxz )*(dfy  - dgy  )/m_d
        +              sin(2*m_pi*(fval - gval)/m_d)*(dfxyz - dgxyz)/(6*m_pi)
        +            4*sin(  m_pi*(fval - gval)/m_d)*(dfxyz - dgxyz)/(3*m_pi)
        +                         (fval - gval     )*(dfxyz - dgxyz)/m_d;
      if(isBogus(retval))
      {
        amrex::Error("bogosity detected aemxyz");
      }
    }
#endif
    else
    {
      amrex::Error("missed a case third deriv aem");
    }
    return retval;
  }
///
  Real
  SmoothAbsoluteValue::
  fourthDerivAem(const  IntVect& a_deriv,
                 const RealVect& a_point) const
  {
    BL_ASSERT(a_deriv.sum() == 4);

    //these are the different forms that the fourth deriv can take
    bool xxxx = (a_deriv.max() == 4);
    bool xxxy = (a_deriv.max() == 3);
#if BL_SPACEDIM==2
    bool xxyy = (a_deriv.max() == 2);
#else 
    //there are more cases in 3D
    bool xxyy = ((a_deriv.max() == 2) && (a_deriv.min() == 0));
    bool xxyz = ((a_deriv.max() == 2) && (a_deriv.min() == 1));
#endif

    Real retval;
    if(xxxx)
    {
      bool found = false;
      int ix;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        if(a_deriv[idir] == 4)
        {
          ix = idir;
          found = true;
          break;
        }
      }
      if(!found) amrex::Error("logic error gs4.0");

      IntVect  dx    =    BASISV(ix);
      IntVect  dxx   =  2*BASISV(ix);
      IntVect  dxxx  =  3*BASISV(ix);
      IntVect  dxxxx =  4*BASISV(ix);


      Real fval  = (*m_f).value(a_point);
      Real gval  = (*m_g).value(a_point);

      Real  dfx    =  (*m_f).value( dx    , a_point);
      Real  dfxx   =  (*m_f).value( dxx   , a_point);
      Real  dfxxx  =  (*m_f).value( dxxx  , a_point);
      Real  dfxxxx =  (*m_f).value( dxxxx , a_point);

      Real  dgx    =  (*m_g).value( dx    , a_point);
      Real  dgxx   =  (*m_g).value( dxx   , a_point);
      Real  dgxxx  =  (*m_g).value( dxxx  , a_point);
      Real  dgxxxx =  (*m_g).value( dxxxx , a_point);
      //(%i51) Aem_xxxx: diff(Aem_xxx,x);
      //
      //(%o51)        sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,4)-'diff(g(x,y,z),x,4))/(6*%pi)
      //    +       4*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,4)-'diff(g(x,y,z),x,4))/(3*%pi)
      //    +                   (f(x,y,z)-g(x,y,z)   )*('diff(f(x,y,z),x,4)-'diff(g(x,y,z),x,4))/d
      //    +       4*cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,3)-'diff(g(x,y,z),x,3))/(3*d)
      //    +      16*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,3)-'diff(g(x,y,z),x,3))/(3*d)
      //    +                                        4*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,3)-'diff(g(x,y,z),x,3))/d
      //    +         cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))^2/d
      //    +       4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))^2/d
      //    +                                        3*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))^2/d
      //    -   4*%pi*sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))/d^2
      //    -   8*%pi*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))/d^2
      //    - 4*%pi^2*cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^4/(3*d^3)
      //    - 4*%pi^2*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^4/(3*d^3)



      retval =  
        +              sin(2*m_pi*(fval - gval)/m_d)*   (dfxxxx - dgxxxx)/(6*m_pi)
        +            4*sin(  m_pi*(fval - gval)/m_d)*   (dfxxxx - dgxxxx)/(3*m_pi)
        +                         (fval - gval     )*   (dfxxxx - dgxxxx)/m_d
        +            4*cos(2*m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*(dfxxx  - dgxxx )/(3*m_d)
        +           16*cos(  m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*(dfxxx  - dgxxx )/(3*m_d)
        +                                          4*   (dfx    - dgx   )*(dfxxx  - dgxxx )/m_d
        +              cos(2*m_pi*(fval - gval)/m_d)*POW(dfxx   - dgxx,2)/m_d
        +            4*cos(  m_pi*(fval - gval)/m_d)*POW(dfxx   - dgxx,2)/m_d
        +                                          3*POW(dfxx   - dgxx,2)/m_d
        -       4*m_pi*sin(2*m_pi*(fval - gval)/m_d)*POW(dfx    - dgx ,2)*(dfxx   - dgxx  )/(POW(m_d,2))
        -       8*m_pi*sin(  m_pi*(fval - gval)/m_d)*POW(dfx    - dgx ,2)*(dfxx   - dgxx  )/(POW(m_d,2))
        -4*POW(m_pi,2)*cos(2*m_pi*(fval - gval)/m_d)*POW(dfx    - dgx ,4)/(3*POW(m_d,3))
        -4*POW(m_pi,2)*cos(  m_pi*(fval - gval)/m_d)*POW(dfx    - dgx ,4)/(3*POW(m_d,3));

      

      if(isBogus(retval))
      {
        amrex::Error("bogosity detected aemxxxx");
      }
    }
    else if(xxxy)
    {
      int ix,iy;
      bool foundx=false, foundy=false;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        if(a_deriv[idir] == 3)
        {
          ix = idir;
          foundx = true;
        }
        else if(a_deriv[idir] == 1)
        {
          iy = idir;
          foundy = true;
        }
      }
      if((!foundx)||(!foundy)) amrex::Error("logic error sab4.05");


      IntVect  dx    =   BASISV(ix);
      IntVect  dy    =   BASISV(iy);
      IntVect  dxx   = 2*BASISV(ix);
      IntVect  dxxx  = 3*BASISV(ix);
      IntVect  dxy   = dx   + dy;
      IntVect  dxxy  = dxx  + dy;
      IntVect  dxxxy = dxxx + dy;

      Real fval  = (*m_f).value(a_point);
      Real gval  = (*m_g).value(a_point);

      Real  dfx    = (*m_f).value(  dx    , a_point);
      Real  dfy    = (*m_f).value(  dy    , a_point);
      Real  dfxx   = (*m_f).value(  dxx   , a_point);
      Real  dfxxx  = (*m_f).value(  dxxx  , a_point);
      Real  dfxy   = (*m_f).value(  dxy   , a_point);
      Real  dfxxy  = (*m_f).value(  dxxy  , a_point);
      Real  dfxxxy = (*m_f).value(  dxxxy , a_point);

      Real  dgx    = (*m_g).value(  dx    , a_point);
      Real  dgy    = (*m_g).value(  dy    , a_point);
      Real  dgxx   = (*m_g).value(  dxx   , a_point);
      Real  dgxxx  = (*m_g).value(  dxxx  , a_point);
      Real  dgxy   = (*m_g).value(  dxy   , a_point);
      Real  dgxxy  = (*m_g).value(  dxxy  , a_point);
      Real  dgxxxy = (*m_g).value(  dxxxy , a_point);


      //(%i58) Aem_xxxy: diff(Aem_xxx,y);
      //
      //(%o58) 
      //                cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,3)-'diff(g(x,y,z),x,3))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d)
      //      +       4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,3)-'diff(g(x,y,z),x,3))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d)
      //      +                                          ('diff(f(x,y,z),x,3)-'diff(g(x,y,z),x,3))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/d
      //      -   2*%pi*sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/d^2
      //      -   4*%pi*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/d^2
      //      - 4*%pi^2*cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^3*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d^3)
      //      - 4*%pi^2*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^3*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d^3)
      //      +         cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,2,y,1)-'diff(g(x,y,z),x,2,y,1))/d
      //      +       4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,2,y,1)-'diff(g(x,y,z),x,2,y,1))/d
      //      +                                        3*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,2,y,1)-'diff(g(x,y,z),x,2,y,1))/d
      //      +         cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))/d
      //      +       4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))/d
      //      +                                        3*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))/d
      //      -   2*%pi*sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))/d^2
      //      -   4*%pi*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))/d^2
      //      +         sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,3,y,1)-'diff(g(x,y,z),x,3,y,1))/(6*%pi)
      //      +       4*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,3,y,1)-'diff(g(x,y,z),x,3,y,1))/(3*%pi)
      //      +                   (f(x,y,z)-g(x,y,z)   )*('diff(f(x,y,z),x,3,y,1)-'diff(g(x,y,z),x,3,y,1))/d


      retval =         cos(2*m_pi*(fval - gval)/m_d)*   (dfxxx  - dgxxx )*(dfy   - dgy   )/(3*m_d)
        +            4*cos(  m_pi*(fval - gval)/m_d)*   (dfxxx  - dgxxx )*(dfy   - dgy   )/(3*m_d)
        +                                               (dfxxx  - dgxxx )*(dfy   - dgy   )/m_d
        -       2*m_pi*sin(2*m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*(dfxx  - dgxx  )*(dfy   - dgy   )/POW(m_d,2)
        -       4*m_pi*sin(  m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*(dfxx  - dgxx  )*(dfy   - dgy   )/POW(m_d,2)
        -4*POW(m_pi,2)*cos(2*m_pi*(fval - gval)/m_d)*POW(dfx    - dgx, 3)*(dfy   - dgy   )/(3*POW(m_d,3))
        -4*POW(m_pi,2)*cos(  m_pi*(fval - gval)/m_d)*POW(dfx    - dgx, 3)*(dfy   - dgy   )/(3*POW(m_d,3))
        +              cos(2*m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*(dfxxy - dgxxy )/m_d
        +            4*cos(  m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*(dfxxy - dgxxy )/m_d
        +                                          3*   (dfx    - dgx   )*(dfxxy - dgxxy )/m_d
        +              cos(2*m_pi*(fval - gval)/m_d)*   (dfxy   - dgxy  )*(dfxx  - dgxx  )/m_d
        +            4*cos(  m_pi*(fval - gval)/m_d)*   (dfxy   - dgxy  )*(dfxx  - dgxx  )/m_d
        +                                          3*   (dfxy   - dgxy  )*(dfxx  - dgxx  )/m_d
        -       2*m_pi*sin(2*m_pi*(fval - gval)/m_d)*POW(dfx    - dgx ,2)*(dfxy  - dgxy  )/POW(m_d,2)
        -       4*m_pi*sin(  m_pi*(fval - gval)/m_d)*POW(dfx    - dgx ,2)*(dfxy  - dgxy  )/POW(m_d,2)
        +              sin(2*m_pi*(fval - gval)/m_d)*   (dfxxxy - dgxxxy)/(6*m_pi)
        +            4*sin(  m_pi*(fval - gval)/m_d)*   (dfxxxy - dgxxxy)/(3*m_pi)
        +                         (fval - gval     )*   (dfxxxy - dgxxxy)/m_d;

      if(isBogus(retval))
      {
        amrex::Error("bogosity detected aemxxxy");
      }
    }
    else if(xxyy)
    {
      int ix, iy;
      bool foundx, foundy;
      foundx = false;
      foundy = false;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        if(a_deriv[idir] == 2)
        {
          if(!foundx)
          {
            ix = idir;
            foundx = true;
          }
          else
          {
            iy = idir;
            foundy = true;
          }
        }
      }
      if(!foundx || !foundy) amrex::Error("logic error aemxxyy");

      IntVect  dx    =   BASISV(ix);
      IntVect  dy    =   BASISV(iy);
      IntVect  dxx   = 2*BASISV(ix);
      IntVect  dyy   = 2*BASISV(iy);
      IntVect  dxy   = dx   + dy ;
      IntVect  dxxy  = dxx  + dy ;
      IntVect  dxxyy = dxx  + dyy;
      IntVect  dxyy  = dx   + dyy;

      Real fval  = (*m_f).value(a_point);
      Real gval  = (*m_g).value(a_point);

      Real  dfx    =(*m_f).value( dx     , a_point);
      Real  dfy    =(*m_f).value( dy     , a_point);
      Real  dfxx   =(*m_f).value( dxx    , a_point);
      Real  dfyy   =(*m_f).value( dyy    , a_point);
      Real  dfxy   =(*m_f).value( dxy    , a_point);
      Real  dfxxy  =(*m_f).value( dxxy   , a_point);
      Real  dfxxyy =(*m_f).value( dxxyy  , a_point);
      Real  dfxyy  =(*m_f).value( dxyy   , a_point);

      Real  dgx    =(*m_g).value( dx     , a_point);
      Real  dgy    =(*m_g).value( dy     , a_point);
      Real  dgxx   =(*m_g).value( dxx    , a_point);
      Real  dgyy   =(*m_g).value( dyy    , a_point);
      Real  dgxy   =(*m_g).value( dxy    , a_point);
      Real  dgxxy  =(*m_g).value( dxxy   , a_point);
      Real  dgxxyy =(*m_g).value( dxxyy  , a_point);
      Real  dgxyy  =(*m_g).value( dxyy   , a_point);

                                                        
      //(%i71) Aem_xxyy: diff(Aem_xxy,y);
      //
      //(%o71) 
      //                  cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*('diff(f(x,y,z),y,2)-'diff(g(x,y,z),y,2))/(3*d)
      //        +       4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*('diff(f(x,y,z),y,2)-'diff(g(x,y,z),y,2))/(3*d)
      //        +                                          ('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*('diff(f(x,y,z),y,2)-'diff(g(x,y,z),y,2))/d
      //        -   2*%pi*sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2*('diff(f(x,y,z),y,2)-'diff(g(x,y,z),y,2))/(3*d^2)
      //        -   4*%pi*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2*('diff(f(x,y,z),y,2)-'diff(g(x,y,z),y,2))/(3*d^2)
      //        -   2*%pi*sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*  ('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))^2/(3*d^2)
      //        -   4*%pi*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*  ('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))^2/(3*d^2)
      //        - 4*%pi^2*cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))^2/(3*d^3)
      //        - 4*%pi^2*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))^2/(3*d^3)
      //        +       2*cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2,y,1)-'diff(g(x,y,z),x,2,y,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d)
      //        +       8*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2,y,1)-'diff(g(x,y,z),x,2,y,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d)
      //        +                                        2*('diff(f(x,y,z),x,2,y,1)-'diff(g(x,y,z),x,2,y,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/d
      //        -   8*%pi*sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d^2)
      //        -  16*%pi*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d^2)
      //        +       2*cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,y,2)-'diff(g(x,y,z),x,1,y,2))/(3*d)
      //        +       8*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,y,2)-'diff(g(x,y,z),x,1,y,2))/(3*d)
      //        +                                        2*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,y,2)-'diff(g(x,y,z),x,1,y,2))/d
      //        +       2*cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))^2/(3*d)
      //        +       8*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))^2/(3*d)
      //        +                                        2*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))^2/d
      //        +         sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2,y,2)-'diff(g(x,y,z),x,2,y,2))/(6*%pi)
      //        +       4*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2,y,2)-'diff(g(x,y,z),x,2,y,2))/(3*%pi)
      //        +                   (f(x,y,z)-g(x,y,z)   )*('diff(f(x,y,z),x,2,y,2)-'diff(g(x,y,z),x,2,y,2))/d


      retval      =     cos(2*m_pi*(fval - gval)/m_d)*   (dfxx   - dgxx  )*   (dfyy  - dgyy  )/(3*m_d)
        +             4*cos(  m_pi*(fval - gval)/m_d)*   (dfxx   - dgxx  )*   (dfyy  - dgyy  )/(3*m_d)
        +                                                (dfxx   - dgxx  )*   (dfyy  - dgyy  )/m_d
        -        2*m_pi*sin(2*m_pi*(fval - gval)/m_d)*POW(dfx    - dgx, 2)*   (dfyy  - dgyy  )/(3*POW(m_d,2))
        -        4*m_pi*sin(  m_pi*(fval - gval)/m_d)*POW(dfx    - dgx, 2)*   (dfyy  - dgyy  )/(3*POW(m_d,2))
        -        2*m_pi*sin(2*m_pi*(fval - gval)/m_d)*   (dfxx   - dgxx  )*POW(dfy   - dgy, 2)/(3*POW(m_d,2))
        -        4*m_pi*sin(  m_pi*(fval - gval)/m_d)*   (dfxx   - dgxx  )*POW(dfy   - dgy, 2)/(3*POW(m_d,2))
        - 4*POW(m_pi,2)*cos(2*m_pi*(fval - gval)/m_d)*POW(dfx    - dgx, 2)*POW(dfy   - dgy, 2)/(3*POW(m_d,3))
        - 4*POW(m_pi,2)*cos(  m_pi*(fval - gval)/m_d)*POW(dfx    - dgx, 2)*POW(dfy   - dgy, 2)/(3*POW(m_d,3))
        +             2*cos(2*m_pi*(fval - gval)/m_d)*   (dfxxy  - dgxxy )*   (dfy   - dgy   )/(3*m_d)
        +             8*cos(  m_pi*(fval - gval)/m_d)*   (dfxxy  - dgxxy )*   (dfy   - dgy   )/(3*m_d)
        +                                           2*   (dfxxy  - dgxxy )*   (dfy   - dgy   )/m_d
        -        8*m_pi*sin(2*m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*   (dfxy  - dgxy  )*(dfy - dgy)/(3*POW(m_d,2))
        -       16*m_pi*sin(  m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*   (dfxy  - dgxy  )*(dfy - dgy)/(3*POW(m_d,2))
        +             2*cos(2*m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*   (dfxyy - dgxyy )/(3*m_d)
        +             8*cos(  m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*   (dfxyy - dgxyy )/(3*m_d)
        +                                           2*   (dfx    - dgx   )*   (dfxyy - dgxyy )/m_d
        +             2*cos(2*m_pi*(fval - gval)/m_d)*POW(dfxy   - dgxy,2)/(3*m_d)
        +             8*cos(  m_pi*(fval - gval)/m_d)*POW(dfxy   - dgxy,2)/(3*m_d)
        +                                           2*POW(dfxy   - dgxy,2)/m_d
        +               sin(2*m_pi*(fval - gval)/m_d)*   (dfxxyy - dgxxyy)/(6*m_pi)
        +             4*sin(  m_pi*(fval - gval)/m_d)*   (dfxxyy - dgxxyy)/(3*m_pi)
        +                          (fval - gval     )*   (dfxxyy - dgxxyy)/m_d;

      if(isBogus(retval))
      {
        amrex::Error("bogosity detected aemxxyy");
      }
    }
#if BL_SPACEDIM==3
    else if(xxyz)
    {
      int  ix, iy, iz;
      bool foundx = false, foundy = false, foundz = false;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        if(a_deriv[idir] == 2)
        {
          ix = idir;
          foundx = true;
        }
        else if(a_deriv[idir] == 1)
        {
          if(!foundy)
          {
            iy = idir;
            foundy = true;
          }
          else
          {
            iz = idir;
            foundz = true;
          }
        }
        else
        {
          amrex::Error("logic error gs4.2");
        }
      }
      if(!foundx || !foundy || !foundz) amrex::Error("logic error xxyz");
      IntVect  dx     =   BASISV(ix);
      IntVect  dy     =   BASISV(iy);
      IntVect  dz     =   BASISV(iz);
      IntVect  dxx    = 2*BASISV(ix);
      IntVect  dxy    = dx   + dy ;
      IntVect  dxz    = dx   + dz ;
      IntVect  dyz    = dy   + dz ;
      IntVect  dxxy   = dxx  + dy ;
      IntVect  dxxz   = dxx  + dz ;
      IntVect  dxyz   = dx   + dy + dz;
      IntVect  dxxyz  = dxz  + dy + dz;

      Real fval  = (*m_f).value(a_point);
      Real gval  = (*m_g).value(a_point);

      Real   dfx     =(*m_f).value( dx     , a_point);
      Real   dfy     =(*m_f).value( dy     , a_point);
      Real   dfz     =(*m_f).value( dz     , a_point);
      Real   dfxx    =(*m_f).value( dxx    , a_point);
      Real   dfxy    =(*m_f).value( dxy    , a_point);
      Real   dfxz    =(*m_f).value( dxz    , a_point);
      Real   dfyz    =(*m_f).value( dyz    , a_point);
      Real   dfxxy   =(*m_f).value( dxxy   , a_point);
      Real   dfxxz   =(*m_f).value( dxxz   , a_point);
      Real   dfxyz   =(*m_f).value( dxyz   , a_point);
      Real   dfxxyz  =(*m_f).value( dxxyz  , a_point);

      Real   dgx     =(*m_g).value( dx     , a_point);
      Real   dgy     =(*m_g).value( dy     , a_point);
      Real   dgz     =(*m_g).value( dz     , a_point);
      Real   dgxx    =(*m_g).value( dxx    , a_point);
      Real   dgxy    =(*m_g).value( dxy    , a_point);
      Real   dgxz    =(*m_g).value( dxz    , a_point);
      Real   dgyz    =(*m_g).value( dyz    , a_point);
      Real   dgxxy   =(*m_g).value( dxxy   , a_point);
      Real   dgxxz   =(*m_g).value( dxxz   , a_point);
      Real   dgxyz   =(*m_g).value( dxyz   , a_point);
      Real   dgxxyz  =(*m_g).value( dxxyz  , a_point);



      //(%i78) Aem_xxyz: diff(Aem_xxy,z);
      //holy crap that is a big formula 
      //(%o78) 
      //        -    2*%pi*sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))*('diff(f(x,y,z),z,1)-'diff(g(x,y,z),z,1))/(3*d^2)
      //        -    4*%pi*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))*('diff(f(x,y,z),z,1)-'diff(g(x,y,z),z,1))/(3*d^2)
      //        -  4*%pi^2*cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))*('diff(f(x,y,z),z,1)-'diff(g(x,y,z),z,1))/(3*d^3)
      //        -  4*%pi^2*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))*('diff(f(x,y,z),z,1)-'diff(g(x,y,z),z,1))/(3*d^3)
      //        +          cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2,y,1)-'diff(g(x,y,z),x,2,y,1))*('diff(f(x,y,z),z,1)-'diff(g(x,y,z),z,1))/(3*d)
      //        +        4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2,y,1)-'diff(g(x,y,z),x,2,y,1))*('diff(f(x,y,z),z,1)-'diff(g(x,y,z),z,1))/(3*d)
      //        +                                           ('diff(f(x,y,z),x,2,y,1)-'diff(g(x,y,z),x,2,y,1))*('diff(f(x,y,z),z,1)-'diff(g(x,y,z),z,1))/d
      //        -    4*%pi*sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))*('diff(f(x,y,z),z,1)-'diff(g(x,y,z),z,1))/(3*d^2)
      //        -    8*%pi*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))*('diff(f(x,y,z),z,1)-'diff(g(x,y,z),z,1))/(3*d^2)
      //        +          cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*('diff(f(x,y,z),y,1,z,1)-'diff(g(x,y,z),y,1,z,1))/(3*d)
      //        +        4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*('diff(f(x,y,z),y,1,z,1)-'diff(g(x,y,z),y,1,z,1))/(3*d)
      //        +                                           ('diff(f(x,y,z),x,2)-'diff(g(x,y,z),x,2))*('diff(f(x,y,z),y,1,z,1)-'diff(g(x,y,z),y,1,z,1))/d
      //        -    2*%pi*sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2*('diff(f(x,y,z),y,1,z,1)-'diff(g(x,y,z),y,1,z,1))/(3*d^2)
      //        -    4*%pi*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))^2*('diff(f(x,y,z),y,1,z,1)-'diff(g(x,y,z),y,1,z,1))/(3*d^2)
      //        +          cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2,z,1)-'diff(g(x,y,z),x,2,z,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d)
      //        +        4*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2,z,1)-'diff(g(x,y,z),x,2,z,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d)
      //        +                                           ('diff(f(x,y,z),x,2,z,1)-'diff(g(x,y,z),x,2,z,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/d
      //        -    4*%pi*sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,z,1)-'diff(g(x,y,z),x,1,z,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d^2)
      //        -    8*%pi*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,z,1)-'diff(g(x,y,z),x,1,z,1))*('diff(f(x,y,z),y,1)-'diff(g(x,y,z),y,1))/(3*d^2)
      //        +        2*cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))*('diff(f(x,y,z),x,1,z,1)-'diff(g(x,y,z),x,1,z,1))/(3*d)
      //        +        8*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))*('diff(f(x,y,z),x,1,z,1)-'diff(g(x,y,z),x,1,z,1))/(3*d)
      //        +                                       2*('diff(f(x,y,z),x,1,y,1)-'diff(g(x,y,z),x,1,y,1))*('diff(f(x,y,z),x,1,z,1)-'diff(g(x,y,z),x,1,z,1))/d
      //        +        2*cos(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,y,1,z,1)-'diff(g(x,y,z),x,1,y,1,z,1))/(3*d)
      //        +        8*cos(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,y,1,z,1)-'diff(g(x,y,z),x,1,y,1,z,1))/(3*d)
      //        +                                         2*('diff(f(x,y,z),x,1)-'diff(g(x,y,z),x,1))*('diff(f(x,y,z),x,1,y,1,z,1)-'diff(g(x,y,z),x,1,y,1,z,1))/d
      //        +          sin(2*%pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2,y,1,z,1)-'diff(g(x,y,z),x,2,y,1,z,1))/(6*%pi)
      //        +        4*sin(  %pi*(f(x,y,z)-g(x,y,z))/d)*('diff(f(x,y,z),x,2,y,1,z,1)-'diff(g(x,y,z),x,2,y,1,z,1))/(3*%pi)
      //        +                    (f(x,y,z)-g(x,y,z)   )*('diff(f(x,y,z),x,2,y,1,z,1)-'diff(g(x,y,z),x,2,y,1,z,1))/d

      retval = 
        -         2*m_pi*sin(2*m_pi*(fval - gval)/m_d)*   (dfxx   - dgxx  )*(dfy    - dgy   )*(dfz - dgz)/(3*POW(m_d,2))
        -         4*m_pi*sin(  m_pi*(fval - gval)/m_d)*   (dfxx   - dgxx  )*(dfy    - dgy   )*(dfz - dgz)/(3*POW(m_d,2))
        -  4*POW(m_pi,2)*cos(2*m_pi*(fval - gval)/m_d)*POW(dfx    - dgx ,2)*(dfy    - dgy   )*(dfz - dgz)/(3*POW(m_d,3))
        -  4*POW(m_pi,2)*cos(  m_pi*(fval - gval)/m_d)*POW(dfx    - dgx ,2)*(dfy    - dgy   )*(dfz - dgz)/(3*POW(m_d,3))
        +                cos(2*m_pi*(fval - gval)/m_d)*   (dfxxy  - dgxxy )*(dfz    - dgz   )/(3*m_d)
        +              4*cos(  m_pi*(fval - gval)/m_d)*   (dfxxy  - dgxxy )*(dfz    - dgz   )/(3*m_d)
        +                                                 (dfxxy  - dgxxy )*(dfz    - dgz   )/m_d
        -         4*m_pi*sin(2*m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*(dfxy   - dgxy  )*(dfz - dgz)/(3*POW(m_d,2))
        -         8*m_pi*sin(  m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*(dfxy   - dgxy  )*(dfz - dgz)/(3*POW(m_d,2))
        +                cos(2*m_pi*(fval - gval)/m_d)*   (dfxx   - dgxx  )*(dfyz   - dgyz  )/(3*m_d)
        +              4*cos(  m_pi*(fval - gval)/m_d)*   (dfxx   - dgxx  )*(dfyz   - dgyz  )/(3*m_d)
        +                                                 (dfxx   - dgxx  )*(dfyz   - dgyz  )/m_d
        -         2*m_pi*sin(2*m_pi*(fval - gval)/m_d)*POW(dfx    - dgx, 2)*(dfyz   - dgyz  )/(3*POW(m_d,2))
        -         4*m_pi*sin(  m_pi*(fval - gval)/m_d)*POW(dfx    - dgx, 2)*(dfyz   - dgyz  )/(3*POW(m_d,2))
        +                cos(2*m_pi*(fval - gval)/m_d)*   (dfxxz  - dgxxz )*(dfy    - dgy   )/(3*m_d)
        +              4*cos(  m_pi*(fval - gval)/m_d)*   (dfxxz  - dgxxz )*(dfy    - dgy   )/(3*m_d)
        +                                                 (dfxxz  - dgxxz )*(dfy    - dgy   )/m_d
        -         4*m_pi*sin(2*m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*(dfxz   - dgxz  )*(dfy - dgy)/(3*POW(m_d,2))
        -         8*m_pi*sin(  m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*(dfxz   - dgxz  )*(dfy - dgy)/(3*POW(m_d,2))
        +              2*cos(2*m_pi*(fval - gval)/m_d)*   (dfxy   - dgxy  )*(dfxz   - dgxz  )/(3*m_d)
        +              8*cos(  m_pi*(fval - gval)/m_d)*   (dfxy   - dgxy  )*(dfxz   - dgxz  )/(3*m_d)
        +                                            2*   (dfxy   - dgxy  )*(dfxz   - dgxz  )/m_d
        +              2*cos(2*m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*(dfxyz  - dgxyz )/(3*m_d)
        +              8*cos(  m_pi*(fval - gval)/m_d)*   (dfx    - dgx   )*(dfxyz  - dgxyz )/(3*m_d)
        +                                            2*   (dfx    - dgx   )*(dfxyz  - dgxyz )/m_d
        +                sin(2*m_pi*(fval - gval)/m_d)*   (dfxxyz - dgxxyz)/(6*m_pi)
        +              4*sin(  m_pi*(fval - gval)/m_d)*   (dfxxyz - dgxxyz)/(3*m_pi)
        +                           (fval - gval     )*   (dfxxyz - dgxxyz)/m_d;

      if(isBogus(retval))
      {
        amrex::Error("bogosity detected aemxxyz");
      }
    }
#endif
    else
    {
      amrex::Error("missed a case in SmoothAbsoluteValue 4th deriv");
    }
  
    if(isBogus(retval))
    {
      amrex::Error("bogosity detected");
    }

    return retval;
  }
///w = f(x) - g(x)
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
  void
  SmoothAbsoluteValue::
  getWCase(int            & a_case,
           Real           & a_wval,
           const RealVect & a_point)const
  {
    Real fofx = (*m_f).value(a_point);
    Real gofx = (*m_g).value(a_point);
    a_wval  = fofx-gofx;
    if(a_wval >= m_d)
    {
      a_case = 1;
    }
    else if(a_wval <= -m_d)
    {
      a_case = -1;
    }
    else
    {
      a_case = 0;
    }
  }
///
  Real 
  SmoothAbsoluteValue::
  smoothAbsFMinusG(const  IntVect& a_deriv,
                   const RealVect& a_point) const
  {
    Real retval = 0;
    int order = a_deriv.sum();
    Real wval; 
    int icase;
    getWCase(icase, wval, a_point);
    if(icase != 0)
    {
      //when outside the limits of smoothing, reduces to regular absolute value
      Real fval = (*m_f).value(a_deriv, a_point);
      Real gval = (*m_g).value(a_deriv, a_point);
      retval = std::abs(fval-gval);
    }
    else
    {
      if(order == 0)
      {
        retval = valueAem(a_point);
      }
      else if(order == 1)
      {
        retval = firstDerivAem(a_deriv, a_point);
      }
      else if(order == 2)
      {
        retval = secondDerivAem(a_deriv, a_point);
      }
      else if(order == 3)
      {
        retval = thirdDerivAem(a_deriv, a_point);
      }
      else if(order == 4)
      {
        retval = fourthDerivAem(a_deriv, a_point);
      }
      else
      {
        Real fval = (*m_f).value(a_deriv, a_point);
        Real gval = (*m_g).value(a_deriv, a_point);
        retval = std::abs(fval-gval);
      }
      checkAgainstKnown(retval, a_deriv, a_point);
    }
    if(isBogus(retval))
    {
      amrex::Error("bogosity detected sabvalmain");
    }
    return retval;
  }
/////////////////////////////////////////////////////////////
/// double ramp debugging class
/////////////////////////////////////////////////////////////
  Real
  DoubleRampExact::
  value(const IntVect & a_deriv,
        const RealVect& a_point) const
  {
    Real retval = 0;
    int order = a_deriv.sum();
    if(order == 0)
    {
      retval = valueAem(a_point);
    }
    else if(order == 1)
    {
      retval = firstDerivAem(a_deriv, a_point);
    }
    else if(order == 2)
    {
      retval = secondDerivAem(a_deriv, a_point);
    }
    else if(order == 3)
    {
      retval = thirdDerivAem(a_deriv, a_point);
    }
    else if(order == 4)
    {
      retval = fourthDerivAem(a_deriv, a_point);
    }
    else
    {
      amrex::Error("derivs higher than 4 not programmed in this class");
    }
    return retval;
  }
////
  Real
  DoubleRampExact::
  valueAem(const RealVect& a_point) const
  {
    //(%o5) "now replace with known f and g"
    //(%i6) Aem: 4*((12*%pi^2*d*( ( 2*(x-0) + (y - 1))- (-2*(x-1) + (y - 1)))+(6*%pi^2-15)*d^2)/(32*%pi^2)
    //        -(d^2*cos(2*%pi*( ( 2*(x-0) + (y - 1))- (-2*(x-1) + (y - 1)))/d)+16*d^2*cos(%pi*( ( 2*(x-0) + (y - 1))- (-2*(x-1) + (y - 1)))/d)
    //                                          -6*%pi^2*( ( 2*(x-0) + (y - 1))- (-2*(x-1) + (y - 1)))^2)
    //         /(32*%pi^2))
    //       /(3*d)
    //       -4*((d^2*cos(2*%pi*( ( 2*(x-0) + (y - 1))- (-2*(x-1) + (y - 1)))/d)
    //          +16*d^2*cos(%pi*( ( 2*(x-0) + (y - 1))- (-2*(x-1) + (y - 1)))/d)-6*%pi^2*( ( 2*(x-0) + (y - 1))- (-2*(x-1) + (y - 1)))^2)
    //          /(32*%pi^2)
    //          +(12*%pi^2*d*( ( 2*(x-0) + (y - 1))- (-2*(x-1) + (y - 1)))+(15-6*%pi^2)*d^2)/(32*%pi^2))
    //        /(3*d);
    //
    //(%o6) 4*((12*%pi^2*d*(2*x+2*(x-1))+(6*%pi^2-15)*d^2)/(32*%pi^2)
    //        -(d^2*cos(2*%pi*(2*x+2*(x-1))/d)+16*d^2*cos(%pi*(2*x+2*(x-1))/d)
    //                                        -6*%pi^2*(2*x+2*(x-1))^2)
    //         /(32*%pi^2))
    //       /(3*d)
    //       -4*((d^2*cos(2*%pi*(2*x+2*(x-1))/d)+16*d^2*cos(%pi*(2*x+2*(x-1))/d)
    //                                          -6*%pi^2*(2*x+2*(x-1))^2)
    //          /(32*%pi^2)
    //          +(12*%pi^2*d*(2*x+2*(x-1))+(15-6*%pi^2)*d^2)/(32*%pi^2))
    //        /(3*d)
    //(%i7) expand(%);
    //
    //(%o7)  -d*  cos(8*%pi*x/d-4*%pi/d)/(12*%pi^2)
    //       -4*d*cos(4*%pi*x/d-2*%pi/d)/( 3*%pi^2)
    //       + 8*x^2/d-8*x/d
    //        -5*d/(4*%pi^2)
    //       +d/2 +2/d

    Real x  = a_point[0];
    Real retval = 
      -m_d*  cos(8*m_pi*x/m_d-4*m_pi/m_d)/(12*POW(m_pi,2))
      -4*m_d*cos(4*m_pi*x/m_d-2*m_pi/m_d)/( 3*POW(m_pi,2))
      + 8*POW(x,2)/m_d
      - 8*x/m_d
      - 5*m_d/(4*POW(m_pi,2))
      + m_d/2 +2/m_d;

    return retval;
  }
////
  Real
  DoubleRampExact::
  firstDerivAem(const IntVect & a_deriv,
                const RealVect& a_point)   const
  {
    Real retval = 0; //the initialization here is important
    if(a_deriv[0] == 1)
    {
      //only x derivs matter for this func
      //(%i8) Aem_x: diff(Aem,x);
      //
      //(%o8) 4*(3*d/2-(-8*%pi*d*sin(2*%pi*(2*x+2*(x-1))/d)
      //               -64*%pi*d*sin(%pi*(2*x+2*(x-1))/d)-48*%pi^2*(2*x+2*(x-1)))
      //               /(32*%pi^2))
      //       /(3*d)
      //       -4*((-8*%pi*d*sin(2*%pi*(2*x+2*(x-1))/d)
      //          -64*%pi*d*sin(%pi*(2*x+2*(x-1))/d)-48*%pi^2*(2*x+2*(x-1)))
      //          /(32*%pi^2)
      //          +3*d/2)
      //        /(3*d)
      //(%i3) Aem_x: 4*(3*d/2-(-8*%pi*d*sin(2*%pi*(2*x+2*(x-1))/d)
      //               -64*%pi*d*sin(%pi*(2*x+2*(x-1))/d)-48*%pi^2*(2*x+2*(x-1)))
      //               /(32*%pi^2))
      //       /(3*d)
      //       -4*((-8*%pi*d*sin(2*%pi*(2*x+2*(x-1))/d)
      //          -64*%pi*d*sin(%pi*(2*x+2*(x-1))/d)-48*%pi^2*(2*x+2*(x-1)))
      //          /(32*%pi^2)
      //          +3*d/2)
      //        /(3*d);
      //
      //(%o3) 4*(3*d/2-(-8*%pi*d*sin(2*%pi*(2*x+2*(x-1))/d)
      //               -64*%pi*d*sin(%pi*(2*x+2*(x-1))/d)-48*%pi^2*(2*x+2*(x-1)))
      //               /(32*%pi^2))
      // /(3*d)
      // -4*((-8*%pi*d*sin(2*%pi*(2*x+2*(x-1))/d)
      //    -64*%pi*d*sin(%pi*(2*x+2*(x-1))/d)-48*%pi^2*(2*x+2*(x-1)))
      //    /(32*%pi^2)
      //    +3*d/2)
      //  /(3*d)
      //(%i4) expand(%);
      //
      //(%o4) 
      //  2*sin(8*%pi*x/d-4*%pi/d)/(3*%pi)+
      //16 *sin(4*%pi*x/d-2*%pi/d)/(3*%pi)
      //+16*x/d-8/d

      Real x  = a_point[0];
      retval = 
        2*  sin(8*m_pi*x/m_d-4*m_pi/m_d)/(3*m_pi)+
        16 *sin(4*m_pi*x/m_d-2*m_pi/m_d)/(3*m_pi)
        +16*x/m_d-8/m_d;    
    }
    return retval;
  }
////
  Real
  DoubleRampExact::
  secondDerivAem(const IntVect & a_deriv,
                 const RealVect& a_point)  const
  {
    Real retval = 0; //the initialization here is important
    if(a_deriv[0] == 2)
    {
      //only x derivs matter for this func
      //(%i12) Aem_xx: diff(Aem_x,x);
      //
      // (%o12) -(
      //          - 64*%pi^2*cos(2*%pi*(2*x+2*(x-1))/d)
      //          -256*%pi^2*cos(  %pi*(2*x+2*(x-1))/d)
      //          -192*%pi^2)/(12*%pi^2*d)
      //(%i12) Aem_xx: diff(Aem_x,x);
      //
      Real x  = a_point[0];
      retval = -(- 64*POW(m_pi,2)*cos(2*m_pi*(2*x+2*(x-1))/m_d)
                 -256*POW(m_pi,2)*cos(  m_pi*(2*x+2*(x-1))/m_d)
                 -192*POW(m_pi,2))/(12*POW(m_pi,2)*m_d);
    }
    return retval;
  }
////
  Real
  DoubleRampExact::
  thirdDerivAem(const IntVect & a_deriv,
                const RealVect& a_point) const
  {
    Real retval = 0; //the initialization here is important
    if(a_deriv[0] == 3)
    {
      //only x derivs matter for this func
      //(%i20) Aem_xxx: diff(Aem_xx,x);
      //
      //(%o20) -( 512*%pi^3*sin(2*%pi*(2*x+2*(x-1))/d)/d
      //        +1024*%pi^3*sin(  %pi*(2*x+2*(x-1))/d)/d
      //          )/(12*%pi^2*d)

      Real x  = a_point[0];
      retval = -(512  *POW(m_pi,3)*sin(2*m_pi*(2*x+2*(x-1))/m_d)/m_d
                 +1024*POW(m_pi,3)*sin(  m_pi*(2*x+2*(x-1))/m_d)/m_d
        )/(12*POW(m_pi,2)*m_d);
    }
    return retval;
  }
////
  Real
  DoubleRampExact::
  fourthDerivAem(const IntVect & a_deriv,
                 const RealVect& a_point) const
  {
    Real retval = 0; //the initialization here is important
    if(a_deriv[0] == 4)
    {
      //only x derivs matter for this func
      //(%i33) Aem_xxxx: diff(Aem_xxx,x);
      //
      //(%o33) -(4096*%pi^4*cos(2*%pi*(2*x+2*(x-1))/d)/d^2
      //        +4096*%pi^4*cos(  %pi*(2*x+2*(x-1))/d)/d^2)
      //        /(12*%pi^2*d)
      //

      Real x  = a_point[0];
      retval =  -(4096 *POW(m_pi,4)*cos(2*m_pi*(2*x+2*(x-1))/m_d)/POW(m_d,2)
                  +4096*POW(m_pi,4)*cos(  m_pi*(2*x+2*(x-1))/m_d)/POW(m_d,2)
        )/(12*POW(m_pi,2)*m_d);
    }
    //because inside = false
    retval = -retval;
    return retval;
  }


/////////////////////////////////////////////////////////////
/// double sphere debugging class
/////////////////////////////////////////////////////////////
  Real
  DoubleSphereExact::
  value(const IntVect & a_deriv,
        const RealVect& a_point) const
  {
    Real retval = 0;
    int order = a_deriv.sum();
    if(order == 0)
    {
      retval = valueAem(a_point);
    }
    else if(order == 1)
    {
      retval = firstDerivAem(a_deriv, a_point);
    }
    else if(order == 2)
    {
      retval = secondDerivAem(a_deriv, a_point);
    }
    else if(order == 3)
    {
      retval = thirdDerivAem(a_deriv, a_point);
    }
    else if(order == 4)
    {
      retval = fourthDerivAem(a_deriv, a_point);
    }
    else
    {
      amrex::Error("derivs higher than 4 not programmed in this class");
    }
    return retval;
  }
////
  Real
  DoubleSphereExact::
  valueAem(const RealVect& a_point) const
  {

    //(%i3) Aem: integrate((4/(3*d))*y*cos(%pi*((w-y)/(2*d)))^4, y,0,w+d) - integrate((4/(3*d))*y*cos(%pi*((w-y)/(2*d)))^4, y,w-d,0);
    //  Is  w+d  positive, negative, or zero?
    //    positive;
    //(%o3) 4*((12*%pi^2*d*w+(6*%pi^2-15)*d^2)/(32*%pi^2)
    //        -(d^2*cos(2*%pi*w/d)+16*d^2*cos(%pi*w/d)-6*%pi^2*w^2)/(32*%pi^2))
    //       /(3*d)
    //       -4*((d^2*cos(2*%pi*w/d)+16*d^2*cos(%pi*w/d)-6*%pi^2*w^2)/(32*%pi^2)
    //          +(12*%pi^2*d*w+(15-6*%pi^2)*d^2)/(32*%pi^2))
    //        /(3*d)
    //(%i4) f: x^2 + y^2 + z^2 -  0.75^2;
    //
    //(%o4) z^2+y^2+x^2-0.5625
    //(%i5) g:  (x-1)^2 + y^2 + z^2 - 0.75^2;
    //
    //(%o5) z^2+y^2+(x-1)^2-0.5625
    //(%i6) w: f - g;
    //
    //(%o6) x^2-(x-1)^2
    //(%i7) Aem: 4*((12*%pi^2*d*(x^2-(x-1)^2)+(6*%pi^2-15)*d^2)/(32*%pi^2)
    //        -(d^2*cos(2*%pi*(x^2-(x-1)^2)/d)+16*d^2*cos(%pi*(x^2-(x-1)^2)/d)-6*%pi^2*(x^2-(x-1)^2)^2)/(32*%pi^2))
    //       /(3*d)
    //       -4*((d^2*cos(2*%pi*(x^2-(x-1)^2)/d)+16*d^2*cos(%pi*(x^2-(x-1)^2)/d)-6*%pi^2*(x^2-(x-1)^2)^2)/(32*%pi^2)
    //          +(12*%pi^2*d*(x^2-(x-1)^2)+(15-6*%pi^2)*d^2)/(32*%pi^2))
    //        /(3*d);
    //
    //(%o7) 4*((12*%pi^2*d*(x^2-(x-1)^2)+(6*%pi^2-15)*d^2)/(32*%pi^2)
    //        -(d^2*cos(2*%pi*(x^2-(x-1)^2)/d)+16*d^2*cos(%pi*(x^2-(x-1)^2)/d)
    //                                        -6*%pi^2*(x^2-(x-1)^2)^2)
    //         /(32*%pi^2))
    //       /(3*d)
    //       -4*((d^2*cos(2*%pi*(x^2-(x-1)^2)/d)+16*d^2*cos(%pi*(x^2-(x-1)^2)/d)
    //                                          -6*%pi^2*(x^2-(x-1)^2)^2)
    //          /(32*%pi^2)
    //          +(12*%pi^2*d*(x^2-(x-1)^2)+(15-6*%pi^2)*d^2)/(32*%pi^2))
    //        /(3*d)
    //expand(%);
    //
    //(%o4) -d*cos(4*%pi*x/d-2*%pi/d)/(12*%pi^2)
    //    -4*d*cos(2*%pi*x/d-  %pi/d)/(3*%pi^2)
    //     +2*x^2/d-2*x/d-5*d/(4*%pi^2)+d/2
    //     +1/(2*d)
    Real x  = a_point[0];
    Real retval = 
      -  m_d*cos(4*m_pi*x/m_d-2*m_pi/m_d)/(12*POW(m_pi,2))
      -4*m_d*cos(2*m_pi*x/m_d-  m_pi/m_d)/( 3*POW(m_pi,2))
      +2*POW(x,2)/m_d-2*x/m_d   -5*m_d/( 4*POW(m_pi,2))+m_d/2
      + 1/(2*m_d);

    return retval;
  }
////
  Real
  DoubleSphereExact::
  firstDerivAem(const IntVect & a_deriv,
                const RealVect& a_point)   const
  {
    Real retval = 0; //the initialization here is important
    if(a_deriv[0] == 1)
    {
      //only x derivs matter for this func
      //(%i8) Aem_x: diff(Aem,x);
      //
      //(%o8) 4*(3*d*(2*x-2*(x-1))/8-(-2*%pi*d*(2*x-2*(x-1))
      //                                *sin(2*%pi*(x^2-(x-1)^2)/d)
      //                             -16*%pi*d*(2*x-2*(x-1))*sin(%pi*(x^2-(x-1)^2)/d)
      //                             -12*%pi^2*(2*x-2*(x-1))*(x^2-(x-1)^2))
      //                             /(32*%pi^2))
      //       /(3*d)
      //       -4*((-2*%pi*d*(2*x-2*(x-1))*sin(2*%pi*(x^2-(x-1)^2)/d)
      //          -16*%pi*d*(2*x-2*(x-1))*sin(%pi*(x^2-(x-1)^2)/d)
      //          -12*%pi^2*(2*x-2*(x-1))*(x^2-(x-1)^2))
      //          /(32*%pi^2)
      //          +3*d*(2*x-2*(x-1))/8)
      //        /(3*d)
      //(%i9) expand(%);
      //
      //(%o9) sin(4*%pi*x/d-2*%pi/d)/(3*%pi)
      //   +8*sin(2*%pi*x/d-  %pi/d)/(3*%pi)
      //   +4*x/d-2/d
      //
      Real x  = a_point[0];
      retval = sin(4*m_pi*x/m_d-2*m_pi/m_d)/(3*m_pi)
        +    8*sin(2*m_pi*x/m_d-  m_pi/m_d)/(3*m_pi)
        + 4*x/m_d-2/m_d;
      

    }
    return retval;
  }
////
  Real
  DoubleSphereExact::
  secondDerivAem(const IntVect & a_deriv,
                 const RealVect& a_point)  const
  {
    Real retval = 0; //the initialization here is important
    if(a_deriv[0] == 2)
    {
      //only x derivs matter for this func
      //(%i16) Aem_xx: diff(Aem_x,x);
      //expand(%);
      //
      //(%o16) -(-4*%pi^2*(2*x-2*(x-1))^2*cos(2*%pi*(x^2-(x-1)^2)/d)
      //        -16*%pi^2*(2*x-2*(x-1))^2*cos(%pi*(x^2-(x-1)^2)/d)
      //        -12*%pi^2*(2*x-2*(x-1))^2)
      //        /(12*%pi^2*d)
      //(%i17) 
      //(%o17) 4*cos(4*%pi*x/d-2*%pi/d)/(3*d)
      //     +16*cos(2*%pi*x/d-  %pi/d)/(3*d)
      //     +4/d

      Real x  = a_point[0];
      retval = 4*cos(4*m_pi*x/m_d-2*m_pi/m_d)/(3*m_d)
        +     16*cos(2*m_pi*x/m_d-  m_pi/m_d)/(3*m_d)
        +4/m_d;

    }
    return retval;
  }
////
  Real
  DoubleSphereExact::
  thirdDerivAem(const IntVect & a_deriv,
                const RealVect& a_point) const
  {
    Real retval = 0; //the initialization here is important
    if(a_deriv[0] == 3)
    {
      //(%i32) Aem_xxx: diff(Aem_xx,x);
      //expand(%);
      //
      //(%o32) -(8*%pi^3*(2*x-2*(x-1))^3*sin(2*%pi*(x^2-(x-1)^2)/d)/d
      //        +16*%pi^3*(2*x-2*(x-1))^3*sin(%pi*(x^2-(x-1)^2)/d)/d)
      //        /(12*%pi^2*d)
      //(%i33) 
      //(%o33) -16*%pi*sin(4*%pi*x/d-2*%pi/d)/(3*d^2)
      //       -32*%pi*sin(2*%pi*x/d-  %pi/d)/(3*d^2)
      //only x derivs matter for this func
      Real x  = a_point[0];
      retval =
        -16*m_pi*sin(4*m_pi*x/m_d-2*m_pi/m_d)/(3*POW(m_d,2))
        -32*m_pi*sin(2*m_pi*x/m_d-  m_pi/m_d)/(3*POW(m_d,2));


    }
    return retval;
  }
////
  Real
  DoubleSphereExact::
  fourthDerivAem(const IntVect & a_deriv,
                 const RealVect& a_point) const
  {
    Real retval = 0; //the initialization here is important
    if(a_deriv[0] == 4)
    {
      //only x derivs matter for this func
      //(%i58) Aem_xxxx: diff(Aem_xxx,x);
      //expand(%);
      //
      //(%o58) -(16*%pi^4*(2*x-2*(x-1))^4*cos(2*%pi*(x^2-(x-1)^2)/d)/d^2
      //        +16*%pi^4*(2*x-2*(x-1))^4*cos(%pi*(x^2-(x-1)^2)/d)/d^2)
      //        /(12*%pi^2*d)
      //(%i59) 
      //(%o59) -64*%pi^2*cos(4*%pi*x/d-2*%pi/d)/(3*d^3)
      //       -64*%pi^2*cos(2*%pi*x/d-  %pi/d)/(3*d^3)

      Real x  = a_point[0];
      retval =  
        -64*POW(m_pi,2)*cos(4*m_pi*x/m_d-2*m_pi/m_d)/(3*POW(m_d,3))
        -64*POW(m_pi,2)*cos(2*m_pi*x/m_d-  m_pi/m_d)/(3*POW(m_d,3));

    }

    return retval;
  }
/////////////////////////////////////////////////////////////
/// offset double sphere debugging class
/// it gets ugly from here down
/////////////////////////////////////////////////////////////
  Real
  OffsetSphereExact::
  value(const IntVect & a_deriv,
        const RealVect& a_point) const
  {
    if(SpaceDim != 2)
    {
      amrex::Error("offset sphere now currently only handles 2d");
    }

    Real retval = 0;
    int order = a_deriv.sum();
    if(order == 0)
    {
      retval = valueAem(a_point);
    }
    else if(order == 1)
    {
      retval = firstDerivAem(a_deriv, a_point);
    }
    else if(order == 2)
    {
      retval = secondDerivAem(a_deriv, a_point);
    }
    else if(order == 3)
    {
      retval = thirdDerivAem(a_deriv, a_point);
    }
    else if(order == 4)
    {
      retval = fourthDerivAem(a_deriv, a_point);
    }
    else
    {
      amrex::Error("derivs higher than 4 not programmed in this class");
    }
    return retval;
  }
////
  Real
  OffsetSphereExact::
  valueAem(const RealVect& a_point) const
  {
    BL_ASSERT(SpaceDim == 2);

    //  (%i3) f: (x-  0  )^2 + (y-  0  )^2  -  (1/2)^2;
    //
    //  (%o3) y^2+x^2-1/4
    //  (%i4) g: (x-(1/2))^2 + (y-(1/2))^2  -  (1/4)^2;
    //
    //  (%o4) (y-1/2)^2+(x-1/2)^2-1/16
    //  (%i5) w: f- g;
    //
    //  (%o5) y^2-(y-1/2)^2+x^2-(x-1/2)^2-3/16
    //  (%i6) Aem: integrate((4/(3*d))*q*cos(%pi*((w-q)/(2*d)))^4, q,0,w+d) - integrate((4/(3*d))*q*cos(%pi*((w-q)/(2*d)))^4, q,w-d,0);
    //
    //  Is 16*y+16*x+16*d-11 positive, negative or zero?
    //
    //    positive;
    //  (%o6) 4*((48*%pi^2*d*y+48*%pi^2*d*x+(24*%pi^2-60)*d^2-33*%pi^2*d)/(128*%pi^2)
    //           -(128*d^2*cos((16*%pi*y+16*%pi*x-11*%pi)/(8*d))
    //             +2048*d^2*cos((16*%pi*y+16*%pi*x-11*%pi)/(16*d))-768*%pi^2*y^2
    //             +(1056*%pi^2-1536*%pi^2*x)*y-768*%pi^2*x^2+1056*%pi^2*x-363*%pi^2)
    //           /(4096*%pi^2))
    //    /(3*d)
    //    -4*((128*d^2*cos((16*%pi*y+16*%pi*x-11*%pi)/(8*d))
    //         +2048*d^2*cos((16*%pi*y+16*%pi*x-11*%pi)/(16*d))-768*%pi^2*y^2
    //         +(1056*%pi^2-1536*%pi^2*x)*y-768*%pi^2*x^2+1056*%pi^2*x-363*%pi^2)
    //        /(4096*%pi^2)
    //        +(48*%pi^2*d*y+48*%pi^2*d*x+(60-24*%pi^2)*d^2-33*%pi^2*d)/(128*%pi^2))
    //    /(3*d)
    //    (%i7) expand(%);
    ////relevant line:
    //  (%o7) 
    //      -d*cos(2*%pi*y/d+2*%pi*x/d-11*%pi/( 8*d))/(12*%pi^2)
    //    -4*d*cos(  %pi*y/d+  %pi*x/d-11*%pi/(16*d))/( 3*%pi^2)
    //    +y^2/(2*d)+x*y/d-11*y/(16*d)
    //    +x^2/(2*d)      -11*x/(16*d)
    //    -5*d/(4*%pi^2)+d/2+121/(512*d)

    Real x = a_point[0];
    Real y = a_point[1];
    Real retval =
      -  m_d*cos(2*m_pi*y/m_d+2*m_pi*x/m_d-11*m_pi/( 8*m_d))/(12*POW(m_pi,2))
      -4*m_d*cos(  m_pi*y/m_d+  m_pi*x/m_d-11*m_pi/(16*m_d))/( 3*POW(m_pi,2))
      +POW(y,2)/(2*m_d)+x*y/m_d-11*y/(16*m_d)
      +POW(x,2)/(2*m_d)        -11*x/(16*m_d)
      -5*m_d/(4*POW(m_pi,2))+m_d/2+121/(512*m_d);
    return retval;
  }
////
///I am sick of typing the many variations of this
  void
  OffsetSphereExact::
  getXYIXIY(int& a_ix, int& a_iy, Real& a_x, Real& a_y,
            const int &   a_xderivVal,
            const  IntVect& a_deriv,
            const RealVect& a_point) const
  {
    BL_ASSERT(SpaceDim == 2);

    int ix;
    bool found= false;
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      if(a_deriv[idir] == a_xderivVal)
      {
        ix = idir;
        found = true;
      }
    }
    if(!found)
    {
      amrex::Error("logic error offsph1.0");
    }
    //only works because spacedim==2
    int iy  = 1-ix;
    Real x = a_point[ix];
    Real y = a_point[iy];

    a_ix = ix;
    a_iy = iy;
    a_x  = x;
    a_y  = y;
  }

////
  Real
  OffsetSphereExact::
  firstDerivAem(const IntVect & a_deriv,
                const RealVect& a_point)   const
  {
    BL_ASSERT(SpaceDim == 2);
    BL_ASSERT(a_deriv.sum() == 1);
    //(%i8) Aem_x: diff(Aem,x);
    //
    //(%o8) 4*(3*d/8-(-256*%pi*d*sin((16*%pi*y+16*%pi*x-11*%pi)/(8*d))
    //               -2048*%pi*d*sin((16*%pi*y+16*%pi*x-11*%pi)/(16*d))-1536*%pi^2*y
    //               -1536*%pi^2*x+1056*%pi^2)
    //               /(4096*%pi^2))
    // /(3*d)
    // -4*((-256*%pi*d*sin((16*%pi*y+16*%pi*x-11*%pi)/(8*d))
    //    -2048*%pi*d*sin((16*%pi*y+16*%pi*x-11*%pi)/(16*d))-1536*%pi^2*y
    //    -1536*%pi^2*x+1056*%pi^2)
    //    /(4096*%pi^2)
    //    +3*d/8)
    //  /(3*d)
    //(%i9) expand(%);
    //relevant line:
    //(%o9) sin(2*%pi*y/d+2*%pi*x/d-11*%pi/( 8*d))/(6*%pi)
    // +  4*sin(  %pi*y/d+  %pi*x/d-11*%pi/(16*d))/(3*%pi)
    //   +y/d+x/d-11/(16*d)

    int ix, iy;
    Real x,y;
    int xderivVal = 1;
    getXYIXIY(ix,iy,x,y,xderivVal, 
              a_deriv, a_point);

    Real retval = 
      +    sin(2*m_pi*y/m_d+2*m_pi*x/m_d-11*m_pi/( 8*m_d))/(6*m_pi)
      +  4*sin(  m_pi*y/m_d+  m_pi*x/m_d-11*m_pi/(16*m_d))/(3*m_pi)
      +   y/m_d+x/m_d-11/(16*m_d);

    return retval;
  }
////
  Real
  OffsetSphereExact::
  secondDerivAem(const IntVect & a_deriv,
                 const RealVect& a_point)  const
  {
    BL_ASSERT(SpaceDim == 2);
    BL_ASSERT(a_deriv.sum() == 2);
    Real retval = 0; 
    if(a_deriv.max() == 2)
    {
      int ix, iy;
      Real x,y;
      int xderivVal = 2;
      getXYIXIY(ix,iy,x,y,xderivVal, 
                a_deriv, a_point);

      //(%i13) Aem_xx: diff(Aem_x,x);
      //
      //(%o13) -(-512*%pi^2*cos((16*%pi*y+16*%pi*x-11*%pi)/(8*d))
      // -2048*%pi^2*cos((16*%pi*y+16*%pi*x-11*%pi)/(16*d))-1536*%pi^2)
      // /(1536*%pi^2*d)
      //(%i14) expand(%);
      //
      //(%o14) 
      //    cos(2*%pi*y/d+2*%pi*x/d-11*%pi/( 8*d))/(3*d)
      // +4*cos(  %pi*y/d+  %pi*x/d-11*%pi/(16*d))/(3*d)+1/d
      retval = 
        +  cos(2*m_pi*y/m_d+2*m_pi*x/m_d-11*m_pi/( 8*m_d))/(3*m_d)
        +4*cos(  m_pi*y/m_d+  m_pi*x/m_d-11*m_pi/(16*m_d))/(3*m_d)+1/m_d;

    }
    else if(a_deriv.max() == 1)
    {
      Real x = a_point[0];
      Real y = a_point[1];
      //(%i18) Aem_xy: diff(Aem_x,y);
      //
      //(%o18) -(-512*%pi^2*cos((16*%pi*y+16*%pi*x-11*%pi)/(8*d))
      // -2048*%pi^2*cos((16*%pi*y+16*%pi*x-11*%pi)/(16*d))-1536*%pi^2)
      // /(1536*%pi^2*d)
      //(%i19) expand(%);
      //
      //(%o19) 
      //    cos(2*%pi*y/d+2*%pi*x/d-11*%pi/( 8*d))/(3*d)
      // +4*cos(  %pi*y/d+  %pi*x/d-11*%pi/(16*d))/(3*d)+1/d

      retval = 
        +  cos(2*m_pi*y/m_d+2*m_pi*x/m_d-11*m_pi/( 8*m_d))/(3*m_d)
        +4*cos(  m_pi*y/m_d+  m_pi*x/m_d-11*m_pi/(16*m_d))/(3*m_d)+1/m_d;

    }
    else
    {
      amrex::Error("logic error offsph2");
    }
    return retval;
  }
////
  Real
  OffsetSphereExact::
  thirdDerivAem(const IntVect & a_deriv,
                const RealVect& a_point) const
  {
    BL_ASSERT(SpaceDim == 2);
    BL_ASSERT(a_deriv.sum() == 3);
    Real retval = 0; 
    if(a_deriv.max() == 3)
    {
      //(%o20) "3 0 0"
      //(%i21) Aem_xxx: diff(Aem_xx,x);
      //
      //(%o21) -(1024*%pi^3*sin((16*%pi*y+16*%pi*x-11*%pi)/(8*d))/d
      // +2048*%pi^3*sin((16*%pi*y+16*%pi*x-11*%pi)/(16*d))/d)
      // /(1536*%pi^2*d)
      //(%i22) expand(%);
      //
      //(%o22) 
      // -2*%pi*sin(2*%pi*y/d+2*%pi*x/d-11*%pi/( 8*d))/(3*d^2)
      // -4*%pi*sin(  %pi*y/d+  %pi*x/d-11*%pi/(16*d))/(3*d^2)
      int ix, iy;
      Real x,y;
      int xderivVal = 3;
      getXYIXIY(ix,iy,x,y,xderivVal, 
                a_deriv, a_point);

      retval =
        -2*m_pi*sin(2*m_pi*y/m_d+2*m_pi*x/m_d-11*m_pi/( 8*m_d))/(3*POW(m_d,2))
        -4*m_pi*sin(  m_pi*y/m_d+  m_pi*x/m_d-11*m_pi/(16*m_d))/(3*POW(m_d,2));

    }
    else if(a_deriv.max() == 2)
    {
      //(%o25) "2 1 0"
      //(%i26) Aem_xxy: diff(Aem_xx,y);
      //
      //(%o26) -(1024*%pi^3*sin((16*%pi*y+16*%pi*x-11*%pi)/(8*d))/d
      // +2048*%pi^3*sin((16*%pi*y+16*%pi*x-11*%pi)/(16*d))/d)
      // /(1536*%pi^2*d)
      //(%i27) expand(%);
      //
      //(%o27) 
      // -2*%pi*sin(2*%pi*y/d+2*%pi*x/d-11*%pi/( 8*d))/(3*d^2)
      // -4*%pi*sin(  %pi*y/d+  %pi*x/d-11*%pi/(16*d))/(3*d^2)
      int ix, iy;
      Real x,y;
      int xderivVal = 2;
      getXYIXIY(ix,iy,x,y,xderivVal, 
                a_deriv, a_point);
      retval = 
        -2*m_pi*sin(2*m_pi*y/m_d+2*m_pi*x/m_d-11*m_pi/( 8*m_d))/(3*POW(m_d,2))
        -4*m_pi*sin(  m_pi*y/m_d+  m_pi*x/m_d-11*m_pi/(16*m_d))/(3*POW(m_d,2));

    }
    else
    {
      amrex::Error("logic error offsph3 (maybe running in 3d?)");
    }
    return retval;
  }
////
  Real
  OffsetSphereExact::
  fourthDerivAem(const IntVect & a_deriv,
                 const RealVect& a_point) const
  {
    BL_ASSERT(SpaceDim == 2);
    BL_ASSERT(a_deriv.sum() == 4);
    Real retval = 0; 
    if(a_deriv.max() == 4)
    {
      //(%i31) Aem_xxxx: diff(Aem_xxx,x);
      //
      //(%o31) -(2048*%pi^4*cos((16*%pi*y+16*%pi*x-11*%pi)/(8*d))/d^2
      // +2048*%pi^4*cos((16*%pi*y+16*%pi*x-11*%pi)/(16*d))/d^2)
      // /(1536*%pi^2*d)
      //(%i32) expand(%);
      //
      //(%o32)
      // -4*%pi^2*cos(2*%pi*y/d+2*%pi*x/d-11*%pi/( 8*d))/(3*d^3)
      // -4*%pi^2*cos(  %pi*y/d+  %pi*x/d-11*%pi/(16*d))/(3*d^3)
      int ix, iy;
      Real x,y;
      int xderivVal = 4;
      getXYIXIY(ix,iy,x,y,xderivVal, 
                a_deriv, a_point);
      retval = 
        -4*POW(m_pi,2)*cos(2*m_pi*y/m_d+2*m_pi*x/m_d-11*m_pi/( 8*m_d))/(3*POW(m_d,3))
        -4*POW(m_pi,2)*cos(  m_pi*y/m_d+  m_pi*x/m_d-11*m_pi/(16*m_d))/(3*POW(m_d,3));
    }
    else if(a_deriv.max() == 3)
    {
      //(%o35) "3 1 0"
      //(%i36) Aem_xxxy: diff(Aem_xxx,y);
      //
      //(%o36) -(2048*%pi^4*cos((16*%pi*y+16*%pi*x-11*%pi)/(8*d))/d^2
      // +2048*%pi^4*cos((16*%pi*y+16*%pi*x-11*%pi)/(16*d))/d^2)
      // /(1536*%pi^2*d)
      //(%i37) expand(%);
      //
      //(%o37) 
      // -4*%pi^2*cos(2*%pi*y/d+2*%pi*x/d-11*%pi/( 8*d))/(3*d^3)
      // -4*%pi^2*cos(  %pi*y/d+  %pi*x/d-11*%pi/(16*d))/(3*d^3)
      int ix, iy;
      Real x,y;
      int xderivVal = 3;
      getXYIXIY(ix,iy,x,y,xderivVal, 
                a_deriv, a_point);
      retval = 
        -4*POW(m_pi,2)*cos(2*m_pi*y/m_d+2*m_pi*x/m_d-11*m_pi/( 8*m_d))/(3*POW(m_d,3))
        -4*POW(m_pi,2)*cos(  m_pi*y/m_d+  m_pi*x/m_d-11*m_pi/(16*m_d))/(3*POW(m_d,3));
    }
    else if(a_deriv.max() == 2)
    {
      //(%o40) "2 2 0"
      //(%i41) Aem_xxyy: diff(Aem_xxy,y);
      //
      //(%o41) -(2048*%pi^4*cos((16*%pi*y+16*%pi*x-11*%pi)/(8*d))/d^2
      // +2048*%pi^4*cos((16*%pi*y+16*%pi*x-11*%pi)/(16*d))/d^2)
      // /(1536*%pi^2*d)
      //(%i42) expand(%);
      //
      //(%o42) 
      // -4*%pi^2*cos(2*%pi*y/d+2*%pi*x/d-11*%pi/( 8*d))/(3*d^3)
      // -4*%pi^2*cos(  %pi*y/d+  %pi*x/d-11*%pi/(16*d))/(3*d^3)

      Real x = a_point[0];
      Real y = a_point[1];
      retval = 
        -4*POW(m_pi,2)*cos(2*m_pi*y/m_d+2*m_pi*x/m_d-11*m_pi/( 8*m_d))/(3*POW(m_d,3))
        -4*POW(m_pi,2)*cos(  m_pi*y/m_d+  m_pi*x/m_d-11*m_pi/(16*m_d))/(3*POW(m_d,3));
    }
    else
    {
      amrex::Error("logic error offsph4");
    }
    return retval;
  }

}         
