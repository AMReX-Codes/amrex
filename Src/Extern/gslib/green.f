      subroutine green(value,hexrep,gfrac)
c-----------------------------------------------------------------------
c
c Provided with a real value ``value'' this subroutine returns the green
c portion of the color specification.
c
c Note common block "color" and call to "hexa"
c
c-----------------------------------------------------------------------
      real*8            value
      character       hexrep*2,hexa*2
      common /color/  cmin,cmax,cint(4),cscl
      hexrep = '00'
      if(value.lt.cint(1))then
c
c Scale it between (0,0):
c
            integ  = 0
      else if((value.ge.cint(1)).and.(value.lt.cint(2)))then
c
c Scale it between (0,255):
c
            integ = int((value-cint(1))/(cint(2)-cint(1))*255.)
            if(integ.gt.255) integ = 255
            if(integ.lt.0)   integ = 0
      else if((value.ge.cint(2)).and.(value.lt.cint(4)))then
c
c Scale it between (255,255):
c
            integ  = 255
      else if(value.ge.cint(4))then
c
c Scale it between (255,0):
c
            integ=int((cmax-value)/(cmax-cint(4))*255.)
            if(integ.gt.255) integ = 255
            if(integ.lt.0)   integ = 0
      end if
c
c Establish coding and return:
c
      gfrac  = real(integ) / 255.
      hexrep = hexa(integ)
      return
      end
