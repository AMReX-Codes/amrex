
module amrex_interp_module

  use amrex_fort_module
  use amrex_constants_module
  use amrex_bc_types_module
  use amrex_error_module

  implicit none

contains

#define IX_PROJ(A,B) (A+B*iabs(A))/B-iabs(A)

#undef  SLX
#undef  SLY
#undef  SLXY

! ::: 
! ::: --------------------------------------------------------------
! ::: cbinterp:  cell centered bilinear interpolation
! ::: 
! ::: NOTE: it is assumed that the coarse grid array is
! ::: large enough to define interpolated values
! ::: in the region fblo:fbhi on the fine grid
! ::: 
! ::: Inputs/Outputs
! ::: fine        <=>  (modify) fine grid array
! ::: fine_l1,fine_h1   =>  (const)  index limits of fine grid
! ::: fb_l1,fb_h1     =>  (const)  subregion of fine grid to get values
! ::: 
! ::: crse         =>  (const)  coarse grid data 
! ::: crse_l1,crse_h1   =>  (const)  index limits of coarse grid
! ::: 
! ::: lratio       =>  (const)  refinement ratio between levels
! ::: nvar         =>  (const)  number of components in array
! ::: 
! ::: TEMPORARY ARRAYS
! ::: slx,sly,slxy =>  1-D slope arrays
! ::: strip        =>  1-D temp array
! ::: --------------------------------------------------------------
! ::: 
    subroutine AMREX_CBINTERP (crse, crse_l1,crse_h1, cb_l1,cb_h1, &
                              fine, fine_l1,fine_h1, fb_l1,fb_h1, &
                              lratio, nvar, &
                              sl, num_slp, strip, strip_lo, strip_hi, &
                              actual_comp,actual_state) bind(c,name='amrex_cbinterp')

      implicit none

      integer crse_l1,crse_h1
      integer cb_l1,cb_h1
      integer fine_l1,fine_h1
      integer fb_l1,fb_h1
      integer lratio, nvar
      integer num_slp
      integer actual_comp,actual_state
      integer strip_lo, strip_hi
      real(amrex_real)  fine(fine_l1:fine_h1, nvar)
      real(amrex_real)  crse(crse_l1:crse_h1, nvar)
      real(amrex_real)    sl(cb_l1:cb_h1,num_slp)
      real(amrex_real) strip(strip_lo:strip_hi)

#define SLX 1
#define SLY 2
#define SLXY 3

      ! local var
      integer lx, ly
      integer hrat, ic, jc, jfn, jfc, i, j, n
      real(amrex_real) x, y
      real(amrex_real) denom

      denom = one/dble(2*lratio)
      hrat = lratio/2
      do 200 n = 1, nvar 
      ! first fill a strip that will fit
          do ic = cb_l1,cb_h1-1
            sl(ic,SLX) = crse(ic+1,n)-crse(ic,n)
	  enddo

            do lx = 0, lratio-1
              do ic = cb_l1, cb_h1-1
                i = ic*lratio + lx
                x = denom*(two*lx + one)
                strip(i) = crse(ic,n) + x*sl(ic,SLX) 
	      enddo
	    enddo

            ! stuff into output array
            do i = fb_l1, fb_h1 
              fine(i,n) = strip(i-hrat)
            enddo
230       continue
200   continue

    end subroutine AMREX_CBINTERP

#undef  SLX
#undef  SLY
#undef  SLXY

! ::: 
! ::: --------------------------------------------------------------
! ::: quartinterp: quartic conservative interpolation from coarse grid to
! ::: subregion of fine grid defined by (fblo,fbhi)
! ::: 
! ::: Inputs/Outputs
! ::: fine        <=>  (modify) fine grid array
! ::: flo,fhi      =>  (const)  index limits of fine grid
! ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
! ::: nvar         =>  (const)  number of variables in state vector
! ::: lratiox      =>  (const)  refinement ratio between levels
! ::: 
! ::: crse         =>  (const)  coarse grid data
! ::: clo,chi      =>  (const)  index limits of crse grid
! ::: cblo,cbhi    =>  (const)  coarse grid region containing fblo,fbhi and widen by 2 or 4 cells
! :::
! ::: cb2lo,cb2hi  =>  (const)  coarse grid region containing fblo,fbhi
! ::: fb2lo,fb2hi  =>  (const)  fine version of cb2. It could be wider than fb
! ::: 
! ::: TEMPORARY ARRAYS
! ::: ftmp         =>  1-D temp array
! ::: --------------------------------------------------------------
! ::: 
     subroutine AMREX_QUARTINTERP (fine, fine_l1,fine_h1, &
                                  fblo, fbhi, fb2lo, fb2hi, &
                                  crse, crse_l1,crse_h1, &
                                  cblo, cbhi, cb2lo, cb2hi, &
                                  nvar, &
                                  lratiox, &
                                  ftmp, &
                                  bc,actual_comp,actual_state) bind(c,name='amrex_quartinterp')

       implicit none

       integer fine_l1,fine_h1
       integer crse_l1,crse_h1
       integer fblo(1), fbhi(1), fb2lo(1), fb2hi(1)
       integer cblo(1), cbhi(1), cb2lo(1), cb2hi(1)
       integer lratiox, nvar
       integer bc(1,2,nvar)
       integer actual_comp,actual_state
       real(amrex_real) fine(fine_l1:fine_h1,nvar)
       real(amrex_real) crse(crse_l1:crse_h1,nvar)
       real(amrex_real) ftmp(fb2lo(1):fb2hi(1))

!      Local variables
       integer i, ii, n
       real(amrex_real) cL(-2:2)
!       real(amrex_real) cR(-2:2)
       data cL/ -0.01171875D0,  0.0859375D0, 0.5d0, -0.0859375D0, &
                 0.01171875D0 /
!$$$       data cR/  0.01171875D0, -0.0859375D0, 0.5d0,  0.0859375D0, &
!$$$                -0.01171875D0 /
       
       if (lratiox .eq. 2) then
          do n = 1, nvar
             do i = cb2lo(1), cb2hi(1)
                ii = 2*i
                ftmp(ii  ) = 2.d0*(cL(-2)*crse(i-2,n) &
                     +             cL(-1)*crse(i-1,n) &
                     +             cL( 0)*crse(i  ,n) &
                     +             cL( 1)*crse(i+1,n) &
                     +             cL( 2)*crse(i+2,n))
                ftmp(ii+1) = 2.d0*crse(i,n)-ftmp(ii)
!$$$                ftmp(ii+1) = 2.d0*(cR(-2)*crse(i-2,n) &
!$$$                     +             cR(-1)*crse(i-1,n) &
!$$$                     +             cR( 0)*crse(i  ,n) &
!$$$                     +             cR( 1)*crse(i+1,n) &
!$$$                     +             cR( 2)*crse(i+2,n))
             enddo
             do ii = fblo(1), fbhi(1)
                fine(ii,n) = ftmp(ii)
             enddo
          enddo
       else if (lratiox .eq. 4) then
          call amrex_error('AMREX_QUARTINTERP: refinement ratio = 4 TODO')
       else
          call amrex_error('AMREX_QUARTINTERP: unsupported refinement ratio')
       endif

     end subroutine AMREX_QUARTINTERP

end module amrex_interp_module
