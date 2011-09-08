      subroutine red(value,hexrep,rfrac)
c-----------------------------------------------------------------------
c
c Provided with a real value ``value'' this subroutine returns the red
c portion of the color specification.
c
c Note common block "color" and call to "hexa"
c
c-----------------------------------------------------------------------
      real            value
      character       hexrep*2,hexa*2
      common /color/  cmin,cmax,cint(4),cscl
      hexrep='00'
      if(value.lt.cint(1))then
c
c Scale it between (y0,0):
c
            integ=int((cint(1)-value)/(cint(1)-cmin)*cscl)
            if(integ.gt.255) integ = 255
            if(integ.lt.0)   integ = 0
      else if((value.ge.cint(1)).and.(value.lt.cint(3)))then
c
c Scale it between (0,0):
c
            integ  = 0
      else if((value.ge.cint(3)).and.(value.lt.cint(4)))then
c
c Scale it between (0,255):
c
            integ = int((value-cint(3))/(cint(4)-cint(3))*255.)
            if(integ.gt.255) integ = 255
            if(integ.lt.0)   integ = 0
      else if(value.ge.cint(4))then
c
c Scale it between (255,255):
c
            integ  = 255
      end if
c
c Establish coding and return:
c
      rfrac  = real(integ) / 255.
      hexrep = hexa(integ)
      return
      end
