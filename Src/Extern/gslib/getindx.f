      subroutine getindx(n,min,siz,loc,index,inflag)
c-----------------------------------------------------------------------
c
c     Gets the coordinate index location of a point within a grid
c     ***********************************************************
c
c
c n       number of "nodes" or "cells" in this coordinate direction
c min     origin at the center of the first cell
c siz     size of the cells
c loc     location of the point being considered
c index   output index within [1,n]
c inflag  true if the location is actually in the grid (false otherwise
c         e.g., if the location is outside then index will be set to
c         nearest boundary
c
c
c
c-----------------------------------------------------------------------

      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      integer   n,index
      real*8      min,siz,loc
      logical   inflag
c
c Compute the index of "loc":
c
      index = int( (loc-min)/siz + 1.5 )
c
c Check to see if in or out:
c
      if(index.lt.1) then
            index  = 1
            inflag = .false.
      else if(index.gt.n) then
            index  = n
            inflag = .false.
      else
            inflag = .true.
      end if
c
c Return to calling program:
c
      return
      end
