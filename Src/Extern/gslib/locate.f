      subroutine locate(xx,n,is,ie,x,j)
c-----------------------------------------------------------------------
c
c Given an array "xx" of length "n", and given a value "x", this routine
c returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
c must be monotonic, either increasing or decreasing.  j=is-1 or j=ie is
c returned to indicate that x is out of range.
c
c Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
c-----------------------------------------------------------------------

      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      dimension xx(n)
c
c Initialize lower and upper methods:
c
      if(is.le.0) is = 1
      jl = is-1
      ju = ie
      if(xx(n).le.x) then
            j = ie
            return
      end if
c
c If we are not done then compute a midpoint:
c
 10   if(ju-jl.gt.1) then
            jm = (ju+jl)/2
c
c Replace the lower or upper limit with the midpoint:
c
            if((xx(ie).gt.xx(is)).eqv.(x.gt.xx(jm))) then
                  jl = jm
            else
                  ju = jm
            endif
            go to 10
      endif
c
c Return with the array index:
c
      j = jl
      return
      end
