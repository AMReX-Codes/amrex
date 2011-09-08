      double precision function dpowint(xlow,xhigh,ylow,yhigh,xval,pow)
c-----------------------------------------------------------------------
c
c Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
c                 for a value of x and a power pow.
c
c-----------------------------------------------------------------------

      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      parameter(EPSLON=1.0e-20)

      if((xhigh-xlow).lt.EPSLON) then
            dpowint = (yhigh+ylow)/2.0
      else
            dpowint = ylow + (yhigh-ylow)* 
     +                (((xval-xlow)/(xhigh-xlow))**pow)
      end if

      return
      end
