      subroutine numtext(value,str)
c-----------------------------------------------------------------------
c
c This subroutine will write a value into the string so that the label
c of the value can be written on the postscript file with the gramma
c of string.
c
c-----------------------------------------------------------------------
      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      character str*12
      real*8      test
c
c Write the number to a text string:
c
      test = abs(value)
      write(str,'(f12.0)') value
      if(test.le.999999.0) write(str,'(f12.2)') value
      if(test.le.   999.0) write(str,'(f12.3)') value
      if(test.le.     0.9) write(str,'(f12.4)') value
      if(test.le.    0.09) write(str,'(f12.5)') value
      if(test.le.   0.009) write(str,'(f12.6)') value
      if(test.le.  0.0009) write(str,'(f12.7)') value
      if(test.le. 0.00009) write(str,'(f12.8)') value
      if(test.le.0.000009) write(str,'(f12.9)') value
      if(test.eq.     0.0) write(str,'(f12.1)') value
c
c Return with the text string containing the number:
c
      return
      end
