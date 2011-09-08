      subroutine blue(value,hexrep,bfrac)
c-----------------------------------------------------------------------
c
c Provided with a real value ``value'' this subroutine returns the blue
c portion of the color specification.
c
c Note common block "color" and call to "hexa"
c
c-----------------------------------------------------------------------
      real*8            value
      character       hexrep*2,hexa*2
      common /color/  cmin,cmax,cint(4),cscl
      hexrep = '00'
      if(value.lt.cint(2))then
c
c Scale it between (255,255):
c
            integ  = 255
      else if((value.ge.cint(2)).and.(value.lt.cint(3)))then
c
c Scale it between (255,0):
c
            integ = int((cint(3)-value)/(cint(3)-cint(2))*255.)
            if(integ.gt.255) integ = 255
            if(integ.lt.0)   integ = 0
      else if(value.ge.cint(3))then
c
c Scale it between (0,0):
c
            integ  = 0
      end if
c
c Establish coding and return:
c
      bfrac  = real(integ) / 255.
      hexrep = hexa(integ)
      return
      end
