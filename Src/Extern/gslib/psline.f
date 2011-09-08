      subroutine psline(np,x,y,lwidt,idsh)
c-----------------------------------------------------------------------
c
c              Write Postscript line commands to a file
c              ****************************************
c
c
c CALLING ARGUMENTS:
c
c  np           the number of points in the x and y array to join
c  x()          array of x values in the range xmin to xmax
c  y()          array of y values in the range ymin to ymax
c  lwidt        the width of the line (1.0 = dark, 0.5 = light)
c  idsh         Dashing Index
c
c NOTES:
c
c  1. The pxmin,pxmax,.. variables are in the standard 1/72 inch 
c     resolution of the postscript page. If a different scale is 
c     going to be used in the printing set pscl to the scale.
c
c  2. If "idsh" is zero then no dashing is perfomed
c
c
c-----------------------------------------------------------------------
      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      real*8      x(*),y(*),lwidt,lwold
      character dash(10)*24
c
c Common Block for Postscript Output Unit and Scaling:
c
      common /psdata/ pscl,pxmin,pxmax,pymin,pymax,xmin,
     +                xmax,ymin,ymax,lpsout
      save   lwold
c
c Dash Patterns:
c
      data dash/'[40 20] 0 setdash       ',
     +          '[13 14 13 20] 0 setdash ',
     +          '[12 21 4 21] 0 setdash  ',
     +          '[10 10] 0 setdash       ',
     +          '[20 20] 0 setdash       ',
     +          '[30 30] 0 setdash       ',
     +          '[40 40] 0 setdash       ',
     +          '[50 50] 0 setdash       ',
     +          '[50 50] 0 setdash       ',
     +          '[50 50] 0 setdash       '/
c
c Change the line width if necessary:
c
      if(pscl.lt.0.01) pscl = 1.0
      if(idsh.gt.10)   idsh = 10 
      if(lwidt.ne.lwold) then
            width = lwidt/pscl
            write(lpsout,100) width
 100        format(f6.3,' setlinewidth')
            lwold = lwidt
      endif
c
c Start a new path and loop through the points:
c
      if(idsh.gt.0) write(lpsout,'(a24)') dash(idsh)
      write(lpsout,101)
 101  format('n')
      do i=1,np
            ix = int(resc(xmin,xmax,pxmin,pxmax,x(i))/pscl)
            iy = int(resc(ymin,ymax,pymin,pymax,y(i))/pscl)
            if(i.eq.1) then
                  write(lpsout,102) ix,iy
 102              format(i5,1x,i5,' m')
            else
                  write(lpsout,103) ix,iy
 103              format(i5,1x,i5,' l')
            endif
      end do      
      write(lpsout,104)
 104  format('s')
      if(idsh.gt.0) write(lpsout,105)
 105  format('[] 0 setdash')
c
c Finished - Return to calling program:
c
      return
      end
