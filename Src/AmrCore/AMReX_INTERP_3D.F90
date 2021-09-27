
module amrex_interp_module

  use amrex_fort_module
  use amrex_constants_module
  use amrex_bc_types_module
  use amrex_error_module

  implicit none

contains

    subroutine AMREX_CQINTERP (fine, fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3, &
                              fb_l1,fb_l2,fb_l3,fb_h1,fb_h2,fb_h3, &
                              nvar, lratiox, lratioy, lratioz, crse, &
                              clo, chi, cb_l1,cb_l2,cb_l3,cb_h1,cb_h2,cb_h3, &
                              fslo, fshi, cslope, clen, fslope, fdat, &
                              flen, voff, bc, limslope, &
                              fvcx, fvcy, fvcz, cvcx, cvcy, cvcz, &
                              actual_comp, actual_state) bind(c,name='amrex_cqinterp')

      implicit none

      integer fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3
      integer fb_l1,fb_l2,fb_l3,fb_h1,fb_h2,fb_h3
      integer cb_l1,cb_l2,cb_l3,cb_h1,cb_h2,cb_h3
      integer fslo(3), fshi(3)
      integer nvar, lratiox, lratioy, lratioz
      integer bc(3,2,nvar)
      integer clen, flen, clo, chi, limslope
      integer actual_comp,actual_state
      real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,fine_l3:fine_h3,nvar)
      real(amrex_real) crse(clo:chi, nvar)
      real(amrex_real) cslope(clo:chi, 3)
      real(amrex_real) fslope(flen, 3)
      real(amrex_real) fdat(flen)
      real(amrex_real) voff(flen)
      real(amrex_real) fvcx(fb_l1:fb_h1+1)
      real(amrex_real) fvcy(fb_l2:fb_h2+1)
      real(amrex_real) fvcz(fb_l3:fb_h3+1)
      real(amrex_real) cvcx(cb_l1:cb_h1+1)
      real(amrex_real) cvcy(cb_l2:cb_h2+1)
      real(amrex_real) cvcz(cb_l3:cb_h3+1)

      call bl_abort('QUADRATIC INTERP NOT IMPLEMENTED IN 3-D')

    end subroutine AMREX_CQINTERP

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
! ::: lratio[xyz]  =>  (const)  refinement ratio between levels
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
! ::: ctmp         =>  2-D temp array
! ::: ctmp2        =>  2-D temp array
! ::: --------------------------------------------------------------
! :::
     subroutine AMREX_QUARTINTERP (fine, fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3, &
                                  fblo, fbhi, fb2lo, fb2hi, &
                                  crse, crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3, &
                                  cblo, cbhi, cb2lo, cb2hi, &
                                  nvar, &
                                  lratiox, lratioy, lratioz, &
                                  ftmp, ctmp, ctmp2, &
                                  bc,actual_comp,actual_state) bind(c,name='amrex_quartinterp')

       implicit none

       integer fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3
       integer crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3
       integer fblo(3), fbhi(3), fb2lo(3), fb2hi(3)
       integer cblo(3), cbhi(3), cb2lo(3), cb2hi(3)
       integer nvar,lratiox,lratioy,lratioz
       integer bc(3,2,nvar)
       integer actual_comp,actual_state
       real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,fine_l3:fine_h3,nvar)
       real(amrex_real) crse(crse_l1:crse_h1,crse_l2:crse_h2,crse_l3:crse_h3,nvar)
       real(amrex_real) ftmp(fb2lo(1):fb2hi(1))
       real(amrex_real) ctmp(cblo(1):cbhi(1),0:lratioy-1)
       real(amrex_real) ctmp2(cblo(1):cbhi(1),cblo(2):cbhi(2),0:lratioz-1)

!      Local variables
       integer i,j,k,ii,jj,kk,n,iry,irz
       real(amrex_real) cL(-2:2)
!       real(amrex_real) cR(-2:2)
       data cL/ -0.01171875D0,  0.0859375D0, 0.5d0, -0.0859375D0, &
                 0.01171875D0 /
!       data cR/  0.01171875D0, -0.0859375D0, 0.5d0,  0.0859375D0, &
!                -0.01171875D0 /

       if (lratiox.eq.2 .and. lratioy.eq.2 .and. lratioz.eq.2) then

          do n = 1, nvar
          do k = cb2lo(3), cb2hi(3)

             do j = cblo(2), cbhi(2)
             do i = cblo(1), cbhi(1)
                ctmp2(i,j,0) = 2.d0*(cL(-2)*crse(i,j,k-2,n) &
                     +               cL(-1)*crse(i,j,k-1,n) &
                     +               cL( 0)*crse(i,j,k  ,n) &
                     +               cL( 1)*crse(i,j,k+1,n) &
                     +               cL( 2)*crse(i,j,k+2,n))
                ctmp2(i,j,1) = 2.d0*crse(i,j,k,n) - ctmp2(i,j,0)
!$$$                ctmp2(i,j,1) = 2.d0*(cR(-2)*crse(i,j,k-2,n) &
!$$$                     +               cR(-1)*crse(i,j,k-1,n) &
!$$$                     +               cR( 0)*crse(i,j,k  ,n) &
!$$$                     +               cR( 1)*crse(i,j,k+1,n) &
!$$$                     +               cR( 2)*crse(i,j,k+2,n))
             enddo
             enddo

             do irz = 0, 1
                kk = k*2+irz
                if (kk.ge.fblo(3) .and. kk.le.fbhi(3)) then

                   do j = cb2lo(2), cb2hi(2)

                      do i = cblo(1), cbhi(1)
                         ctmp(i,0) = 2.d0*(cL(-2)*ctmp2(i,j-2,irz) &
                              +            cL(-1)*ctmp2(i,j-1,irz) &
                              +            cL( 0)*ctmp2(i,j  ,irz) &
                              +            cL( 1)*ctmp2(i,j+1,irz) &
                              +            cL( 2)*ctmp2(i,j+2,irz))
                         ctmp(i,1) = 2.d0*ctmp2(i,j,irz) - ctmp(i,0)
!$$$                         ctmp(i,1) = 2.d0*(cR(-2)*ctmp2(i,j-2,irz) &
!$$$                              +            cR(-1)*ctmp2(i,j-1,irz) &
!$$$                              +            cR( 0)*ctmp2(i,j  ,irz) &
!$$$                              +            cR( 1)*ctmp2(i,j+1,irz) &
!$$$                              +            cR( 2)*ctmp2(i,j+2,irz))
                      enddo

                      do iry = 0, 1
                         jj = j*2+iry

                         if (jj.ge.fblo(2).and.jj.le.fbhi(2)) then
                            do i = cb2lo(1), cb2hi(1)
                               ii = 2*i
                               ftmp(ii  ) = 2.d0*(cL(-2)*ctmp(i-2,iry) &
                                    +             cL(-1)*ctmp(i-1,iry) &
                                    +             cL( 0)*ctmp(i  ,iry) &
                                    +             cL( 1)*ctmp(i+1,iry) &
                                    +             cL( 2)*ctmp(i+2,iry))
                               ftmp(ii+1) = 2.d0*ctmp(i,iry) - ftmp(ii)
!$$$                               ftmp(ii+1) = 2.d0*(cR(-2)*ctmp(i-2,iry) &
!$$$                                    +             cR(-1)*ctmp(i-1,iry) &
!$$$                                    +             cR( 0)*ctmp(i  ,iry) &
!$$$                                    +             cR( 1)*ctmp(i+1,iry) &
!$$$                                    +             cR( 2)*ctmp(i+2,iry))
                            enddo
                            do ii = fblo(1), fbhi(1)
                               fine(ii,jj,kk,n) = ftmp(ii)
                            enddo
                         endif  ! if (jj.ge.......
                      enddo  ! do iry

                   enddo  ! do j

                endif  ! if (kk.ge.......
             enddo  ! do irz

          enddo  ! do k
          enddo  ! do n

       else if (lratiox.eq.4 .and. lratioy.eq.4 .and. lratioz.eq.4) then
          call amrex_error('AMREX_QUARTINTERP: refinement ratio = 4 TODO')
       else
          call amrex_error('AMREX_QUARTINTERP: unsupported refinement ratio')
       endif

     end subroutine AMREX_QUARTINTERP

end module amrex_interp_module
