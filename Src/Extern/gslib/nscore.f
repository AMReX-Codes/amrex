      subroutine nscore(nd,vr,tmin,tmax,iwt,wt,tmp,lout,vrg,ierror)
c-----------------------------------------------------------------------
c
c              Transform Univariate Data to Normal Scores
c              ******************************************
c
c This subroutibe takes "nd" data "vr(i),i=1,...,nd" possibly weighted
c by "wt(i),i=,...,nd" and returns the normal scores transform N(0,1)
c as "vrg(i),i=1,...,nd".  The extra storage array "tmp" is required
c so that the data can be returned in the same order (just in case
c there are associated arrays like the coordinate location).
c
c
c
c INPUT VARIABLES:
c
c   nd               Number of data (no missing values)
c   vr(nd)           Data values to be transformed
c   tmin,tmax        data trimming limits
c   iwt              =0, equal weighted; =1, then apply weight
c   wt(nd)           Weight for each data (don't have to sum to 1.0)
c   tmp(nd)          Temporary storage space for sorting
c   lout             if > 0 then transformation table will be written
c
c
c
c OUTPUT VARIABLES:
c
c   vrg(nd)          normal scores
c   ierror           error flag (0=error free,1=problem)
c
c
c
c EXTERNAL REFERENCES:
c
c   gauinv           Calculates the inverse of a Gaussian cdf
c   sortem           sorts a number of arrays according to a key array
c
c
c
c-----------------------------------------------------------------------

      implicit integer (i-n)
      implicit double precision (a-h,o-z)

      parameter(EPSLON=1.0e-20)
      real*8      vr(nd),wt(nd),vrg(nd),tmp(nd)
      real*8    pd
c
c Sort the data in ascending order and calculate total weight:
c
      ierror = 0
      twt    = 0.d0
      do i=1,nd
            tmp(i) = real(i)
            if(vr(i).ge.tmin.and.vr(i).lt.tmax) then
                  if(iwt.eq.0) then
                        twt = twt + 1.
                  else
                        twt = twt + wt(i)
                  end if
            end if
      end do
      if(nd.lt.1.or.twt.lt.EPSLON) then
            ierror = 1
            return
      end if
      call sortem(1,nd,vr,2,wt,tmp,d,e,f,g,h)
c
c Compute the cumulative probabilities:
c
      oldcp = 0.d0
      cp    = 0.d0
      do i=1,nd
            cp     =  cp + wt(i) / twt
            wt(i)  = (cp + oldcp)/ 2.0
            oldcp  =  cp
            call gauinv(dble(wt(i)),vrg(i),ierr)
            if(lout.gt.0) write(lout,'(f12.5,1x,f12.5)') vr(i),vrg(i)
      end do
c
c Get the arrays back in original order:
c
      call sortem(1,nd,tmp,3,wt,vr,vrg,e,f,g,h)
c
c Finished:
c
      return
      end
