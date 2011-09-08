      subroutine ordrel(ivtype,ncut,ccdf,ccdfo,nviol,aviol,xviol)
c-----------------------------------------------------------------------
c
c                 Correct Order Relation Problems
c                 *******************************
c
c This subroutine identifies and corrects order relation problems in a
c conditional distribution known at a specified number of cutoffs.
c
c
c
c INPUT VARIABLES:
c
c   ivtype           variable type (0=categorical, 1=continuous)
c   ncut             number of cutoffs
c   ccdf(i)          input ccdf values
c
c
c OUTPUT VARIABLES:
c
c   ccdfo            corrected ccdf values
c   nviol()          number of order relation violations
c   aviol()          average magnitude of the order relation violations
c   xviol()          maximum magnitude of the order relation violations
c
c
c
c PROGRAMMING NOTES:
c
c   1. the arrays ccdf1 and ccdf2 are used for temporary storage of the
c      ccdf corrected sequentially upwards and downwards.  The program
c      execution will be stopped if the memory allocation of these two
c      arrays is not sufficient.
c   
c
c
c-----------------------------------------------------------------------
      implicit integer (i-n)
      implicit double precision (a-h,o-z)


      parameter(MAXCUT=100)
      real*8      ccdf(*),ccdfo(*),aviol(*),xviol(*)
      real*8      ccdf1(MAXCUT),ccdf2(MAXCUT)
      integer   nviol(*)
c
c Make sure there is enough temporary storage: 
c
      if(ncut.gt.MAXCUT) then
            write(*,100) MAXCUT,ncut
 100        format('There is not enough temporary storage allocated'
     +          ,/,'in subroutine ordrel: increase and recompile'
     +          ,/,'      available = ',i3
     +          ,/,'      required  = ',i3)
            stop
      endif
c
c Make sure conditional cdf is within [0,1]:
c
      do i=1,ncut
            if(ccdf(i).lt.0.d0) then
                  ccdf1(i) = 0.d0
                  ccdf2(i) = 0.d0
            else if(ccdf(i).gt.1.0) then
                  ccdf1(i) = 1.0
                  ccdf2(i) = 1.0
            else
                  ccdf1(i) = ccdf(i)
                  ccdf2(i) = ccdf(i)
            endif
      end do
c
c Correct sequentially up, then down, and then average:
c
      if(ivtype.eq.0) then
            sumcdf = 0.d0
            do i=1,ncut
                  sumcdf = sumcdf + ccdf1(i)
            end do
            if(sumcdf.le.0.d0) sumcdf = 1.0
            do i=1,ncut
                  ccdfo(i) = ccdf1(i) / sumcdf
            end do
      else
            do i=2,ncut
                  if(ccdf1(i).lt.ccdf1(i-1)) ccdf1(i) = ccdf1(i-1)
            end do
            do i=ncut-1,1,-1
                  if(ccdf2(i).gt.ccdf2(i+1)) ccdf2(i) = ccdf2(i+1)
            end do
            do i=1,ncut
                  ccdfo(i) = 0.5*(ccdf1(i)+ccdf2(i))
            end do
      end if
c
c Accumulate error statistics:
c
      do i=1,ncut
            if(ccdf(i).ne.ccdfo(i)) then
                  viol = abs(ccdf(i)-ccdfo(i))
                  nviol(i) = nviol(i) + 1
                  aviol(i) = aviol(i) + viol
                  xviol(i) = max(xviol(i),viol)
            endif
      end do
c
c Return with corrected CDF:
c
      return
      end
