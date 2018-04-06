
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_INTERP_F.H"
#include <AMReX_ArrayLim.H>

#define IX_PROJ(A,B) (A+B*iabs(A))/B-iabs(A)
#define SDIM 2


! ::: --------------------------------------------------------------
! ::: nbinterp:  node based bilinear interpolation
! :::
! ::: INPUTS/OUTPUTS
! ::: fine        <=>  (modify) fine grid array
! ::: fine_l1,fine_l2,fine_h1,fine_h2   =>  (const)  index limits of fine grid
! ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
! :::
! ::: crse         =>  (const)  coarse grid data widened by 1 zone
! ::: crse_l1,crse_l2,crse_h1,crse_h2   =>  (const)  index limits of coarse grid
! :::
! ::: lratio(3)    =>  (const)  refinement ratio between levels
! ::: nvar         =>  (const)  number of components in array
! ::: num_slp      =>  (const)  number of types of slopes
! :::
! ::: TEMPORARY ARRAYS
! ::: sl           =>  num_slp 1-D slope arrays
! ::: --------------------------------------------------------------
! :::
    subroutine FORT_NBINTERP (crse, crse_l1,crse_l2,crse_h1,crse_h2, cb_l1,cb_l2,cb_h1,cb_h2, &
                              fine, fine_l1,fine_l2,fine_h1,fine_h2, fb_l1,fb_l2,fb_h1,fb_h2, &
                              lratiox, lratioy, nvar, &
                              sl, num_slp, &
                              actual_comp,actual_state)

      implicit none

      integer crse_l1,crse_l2,crse_h1,crse_h2
      integer cb_l1,cb_l2,cb_h1,cb_h2
      integer fine_l1,fine_l2,fine_h1,fine_h2
      integer fb_l1,fb_l2,fb_h1,fb_h2
      integer lratiox, lratioy, nvar
      integer num_slp
      integer actual_comp,actual_state
      REAL_T  fine(fine_l1:fine_h1,fine_l2:fine_h2,nvar)
      REAL_T  crse(crse_l1:crse_h1,crse_l2:crse_h2,nvar)
      REAL_T  sl(DIM1(cb),num_slp)

#define  SLX 1
#define  SLY 2
#define  SLXY 3

      integer lx, ly
      integer i, j, ifn, jfn, n
      integer ilo, ihi, jlo, jhi
      integer jstrtFine, jstopFine, istrtFine, istopFine

      REAL_T fx, fy
      REAL_T RX, RY, RXY
      REAL_T dx0, d0x, dx1
      REAL_T slope

      slope(i,j,n,fx,fy) = crse(i,j,n) + &
                           fx*sl(i,SLX) + fy*sl(i,SLY) + fx*fy*sl(i,SLXY)

      RX = one/dble(lratiox)
      RY = one/dble(lratioy)
      RXY = RX*RY

!     NOTES:
!         1) (i, j) loop over the coarse cells
!         2) ?strtFine and ?stopFine are the beginning and ending fine cell
!            indices corresponding to the current coarse cell.  ?stopFine
!            is restricted for the last coarse cell in each direction since
!            for this cell we only need to do the face and not the fine nodes
!            inside this cell.
!         3) (lx, ly) as well as ?lo and ?hi refer to the fine node indices
!            as an offset from ?strtFine.

      do 100 n = 1, nvar
        do 120 j = ARG_L2(cb), ARG_H2(cb)
          jstrtFine = j * lratioy
          jstopFine = jstrtFine + lratioy - 1
          if (j .eq. ARG_H2(cb)) jstopFine = jstrtFine

          jlo = max(ARG_L2(fb),jstrtFine) - jstrtFine
          jhi = min(ARG_H2(fb),jstopFine) - jstrtFine

!         ::::: compute slopes :::::
!
!         NOTE: The IF logic in the calculation of the slopes is to
!               prevent stepping out of bounds on the coarse data when
!               computing the slopes on the ARG_H?(cb) cells.  These
!               slopes actually are not used since they are multiplied by
!               zero.

          do i = ARG_L1(cb), ARG_H1(cb)
            dx0 = zero
            if (i .NE. ARG_H1(cb)) dx0 = crse(i+1,j,n) - crse(i,j,n)

            d0x = zero
            if (j .NE. ARG_H2(cb)) d0x = crse(i,j+1,n) - crse(i,j,n)

            dx1 = zero
            if (i .NE. ARG_H1(cb) .and. j .NE. ARG_H2(cb)) &
              dx1 = crse(i+1,j+1,n) - crse(i,j+1,n)

            sl(i,SLX) = RX*dx0
            sl(i,SLY) = RY*d0x
            sl(i,SLXY) = RXY*(dx1 - dx0)
          end do

!         ::::: compute fine strip of interpolated data

          do ly = jlo, jhi
            jfn = lratioy * j + ly
            fy = dble(ly)

            do i = ARG_L1(cb), ARG_H1(cb)
              istrtFine = i * lratiox
              istopFine = istrtFine + lratiox - 1
              if (i .eq. ARG_H1(cb)) istopFine = istrtFine

              ilo = max(ARG_L1(fb),istrtFine) - istrtFine
              ihi = min(ARG_H1(fb),istopFine) - istrtFine

              do lx = ilo, ihi
                ifn = lratiox * i + lx
                fx = dble(lx)

                fine(ifn,jfn,n) = slope(i,j,n,fx,fy)
              end do
            end do
          end do
120     continue
100   continue

#undef  SLX
#undef  SLY
#undef  SLXY

    end subroutine FORT_NBINTERP
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
! ::: fine_l1,fine_l2,fine_h1,fine_h2   =>  (const)  index limits of fine grid
! ::: fb_l1,fb_l2,fb_h1,fb_h2     =>  (const)  subregion of fine grid to get values
! ::: 
! ::: crse         =>  (const)  coarse grid data 
! ::: crse_l1,crse_l2,crse_h1,crse_h2   =>  (const)  index limits of coarse grid
! ::: 
! ::: lratio(2)    =>  (const)  refinement ratio between levels
! ::: nvar         =>  (const)  number of components in array
! ::: 
! ::: TEMPORARY ARRAYS
! ::: slx,sly,slxy =>  1-D slope arrays
! ::: strip        =>  1-D temp array
! ::: --------------------------------------------------------------
! ::: 
    subroutine FORT_CBINTERP (crse, crse_l1,crse_l2,crse_h1,crse_h2, cb_l1,cb_l2,cb_h1,cb_h2, &
                              fine, fine_l1,fine_l2,fine_h1,fine_h2, fb_l1,fb_l2,fb_h1,fb_h2, &
                              lratiox, lratioy, nvar, &
                              sl, num_slp, strip, strip_lo, strip_hi, &
                              actual_comp,actual_state)

      implicit none

      integer crse_l1,crse_l2,crse_h1,crse_h2
      integer cb_l1,cb_l2,cb_h1,cb_h2
      integer fine_l1,fine_l2,fine_h1,fine_h2
      integer fb_l1,fb_l2,fb_h1,fb_h2
      integer lratiox, lratioy, nvar
      integer num_slp
      integer actual_comp,actual_state
      integer strip_lo, strip_hi
      REAL_T  fine(fine_l1:fine_h1,fine_l2:fine_h2, nvar)
      REAL_T  crse(crse_l1:crse_h1,crse_l2:crse_h2, nvar)
      REAL_T  sl(DIM1(cb),num_slp)
      REAL_T  strip(strip_lo:strip_hi)

#define SLX 1
#define SLY 2
#define SLXY 3

      integer lx, ly, hratx, hraty, ic, jc, jfn, jfc, i, n
      REAL_T x, y, denomx, denomy

      denomx = one/dble(2*lratiox)
      denomy = one/dble(2*lratioy)

      hratx = lratiox/2
      hraty = lratioy/2

      do n = 1, nvar 
         do jc = ARG_L2(cb), ARG_H2(cb)-1 
            do ic = ARG_L1(cb), ARG_H1(cb)-1
               sl(ic,SLX) = crse(ic+1,jc,n)-crse(ic,jc,n)
               sl(ic,SLY) = crse(ic,jc+1,n)-crse(ic,jc,n)
               sl(ic,SLXY) = crse(ic+1,jc+1,n)-crse(ic+1,jc,n) &
                    - crse(ic  ,jc+1,n)+crse(ic  ,jc,n)
            end do
            do ly = 0, lratioy-1 
               jfn = jc*lratioy + ly
               jfc = jfn + hraty
               if (jfc .ge. ARG_L2(fb)  .and.  jfc .le. ARG_H2(fb)) then
                  y = denomy*(two*ly + one)
                  do lx = 0, lratiox-1
                     do ic = ARG_L1(cb), ARG_H1(cb)-1
                        i = ic*lratiox + lx
                        x = denomx*(two*lx + one)
                        strip(i) = crse(ic,jc,n) + x*sl(ic,SLX) + &
                                   y*sl(ic,SLY) + x*y*sl(ic,SLXY)
                     end do
                  end do
                  do i = ARG_L1(fb), ARG_H1(fb) 
                     fine(i,jfc,n) = strip(i-hratx)
                  end do
               end if
            end do
         end do
      end do

    end subroutine FORT_CBINTERP

#undef  SLX
#undef  SLY
#undef  SLXY

! ::: 
! ::: --------------------------------------------------------------
! ::: linccinterp:   linear conservative interpolation from coarse grid to
! ::: subregion of fine grid defined by (fblo,fbhi)
! ::: 
! ::: The interpolation is linear in that it uses a
! ::: a limiting scheme that preserves the value of 
! ::: any linear combination of the
! ::: coarse grid data components--e.g.,
! ::: if sum_ivar a(ic,jc,ivar)*fab(ic,jc,ivar) = 0, then
! ::: sum_ivar a(ic,jc,ivar)*fab(if,jf,ivar) = 0 is satisfied
! ::: in all fine cells if,jf covering coarse cell ic,jc.
! ::: 
! ::: If lin_limit = 0, the interpolation scheme is identical to
! ::: the used in ccinterp for limslope=1; the results should
! ::: be exactly the same -- difference = hard 0.
! ::: 
! ::: Inputs/Outputs
! ::: fine        <=>  (modify) fine grid array
! ::: flo,fhi      =>  (const)  index limits of fine grid
! ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
! ::: nvar         =>  (const)  number of variables in state vector
! ::: lratio(2)    =>  (const)  refinement ratio between levels
! ::: 
! ::: crse         =>  (const)  coarse grid data widended by 1 zone
! ::: clo,chi      =>  (const)  index limits of crse grid
! ::: cslo,cshi    =>  (const)  coarse grid index limits where
! :::				slopes are to be defined. This is
! :::				the projection of (fblo,fbhi) down
! :::				to the coarse level 
! ::: ucslope      =>  (modify) temp array of unlimited coarse grid slopes
! ::: lcslope      =>  (modify) temp array of limited coarse grid slopes
! ::: slope_factor =>  (modify) temp array of slope limiting factors
! ::: lin_limit    =>  (const)  != 0 => do linear slope limiting scheme
! :::
! ::: --------------------------------------------------------------
! ::: 
     subroutine FORT_LINCCINTERP (fine, fine_l1,fine_l2,fine_h1,fine_h2, fblo, fbhi, &
                                  fvcb_l1,fvcb_l2,fvcb_h1,fvcb_h2, &
                                  crse, crse_l1,crse_l2,crse_h1,crse_h2, cvcb_l1,cvcb_l2,cvcb_h1,cvcb_h2, &
                                  uc_xslope, lc_xslope, xslope_factor, &
                                  uc_yslope, lc_yslope, yslope_factor, &
                                  cslope_l1,cslope_l2,cslope_h1,cslope_h2, &
                                  cslopelo, cslopehi, &
                                  nvar, lratiox, lratioy, &
                                  bc, lim_slope, lin_limit, &
                                  fvcx, fvcy, cvcx, cvcy, &
                                  voffx, voffy, alpha, cmax, cmin, &
                                  actual_comp,actual_state)

       implicit none

       integer fine_l1,fine_l2,fine_h1,fine_h2
       integer crse_l1,crse_l2,crse_h1,crse_h2
       integer fvcb_l1,fvcb_l2,fvcb_h1,fvcb_h2
       integer cvcb_l1,cvcb_l2,cvcb_h1,cvcb_h2
       integer cslope_l1,cslope_l2,cslope_h1,cslope_h2
       integer fblo(2), fbhi(2)
       integer cslopelo(2), cslopehi(2)
       integer lratiox, lratioy, nvar
       integer lim_slope, lin_limit
       integer bc(2,2,nvar)
       integer actual_comp,actual_state
       REAL_T fine(fine_l1:fine_h1,fine_l2:fine_h2,nvar)
       REAL_T crse(crse_l1:crse_h1,crse_l2:crse_h2, nvar)
       REAL_T uc_xslope(cslope_l1:cslope_h1,cslope_l2:cslope_h2,nvar)
       REAL_T lc_xslope(cslope_l1:cslope_h1,cslope_l2:cslope_h2,nvar)
       REAL_T xslope_factor(cslope_l1:cslope_h1,cslope_l2:cslope_h2)
       REAL_T uc_yslope(cslope_l1:cslope_h1,cslope_l2:cslope_h2,nvar)
       REAL_T lc_yslope(cslope_l1:cslope_h1,cslope_l2:cslope_h2,nvar)
       REAL_T yslope_factor(cslope_l1:cslope_h1,cslope_l2:cslope_h2)
       REAL_T alpha(cslope_l1:cslope_h1,cslope_l2:cslope_h2,nvar)
       REAL_T cmax(cslope_l1:cslope_h1,cslope_l2:cslope_h2,nvar)
       REAL_T cmin(cslope_l1:cslope_h1,cslope_l2:cslope_h2,nvar)
       REAL_T fvcx(DIM1(fvcb))
       REAL_T fvcy(DIM2(fvcb))
       REAL_T voffx(DIM1(fvcb))
       REAL_T voffy(DIM2(fvcb))
       REAL_T cvcx(DIM1(cvcb))
       REAL_T cvcy(DIM2(cvcb))

#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)

       integer n
       integer i, ic
       integer j, jc
       REAL_T cen, forw, back, slp
       REAL_T factorn, denom
       REAL_T fxcen, cxcen, fycen, cycen
       REAL_T orig_corr_fact,corr_fact
       REAL_T dummy_fine
       logical xok, yok
       integer ncbx, ncby
       integer ioff,joff
       integer voff_lo(2),voff_hi(2)

       ncbx = cslopehi(1)-cslopelo(1)+1
       ncby = cslopehi(2)-cslopelo(2)+1

       voff_lo(1) = cslopelo(1) * lratiox
       voff_lo(2) = cslopelo(2) * lratioy
       voff_hi(1) = (cslopehi(1)+1) * lratiox - 1
       voff_hi(2) = (cslopehi(2)+1) * lratioy - 1

       xok = (ncbx .ge. 2)
       yok = (ncby .ge. 2)

       do j = voff_lo(2),voff_hi(2)
         jc = IX_PROJ(j,lratioy)
         fycen = half*(fvcy(j)+fvcy(j+1))
         cycen = half*(cvcy(jc)+cvcy(jc+1))
         voffy(j) = (fycen-cycen)/(cvcy(jc+1)-cvcy(jc))
       end do
       do i = voff_lo(1),voff_hi(1)
          ic = IX_PROJ(i,lratiox)
          fxcen = half*(fvcx(i)+fvcx(i+1))
          cxcen = half*(cvcx(ic)+cvcx(ic+1))
          voffx(i) = (fxcen-cxcen)/(cvcx(ic+1)-cvcx(ic))
       end do

       do n = 1, nvar

! ...     Initialize alpha = 1 and define cmax and cmin as neighborhood max/mins.

          do j = cslopelo(2),cslopehi(2)
             do i = cslopelo(1), cslopehi(1)
                alpha(i,j,n) = 1.d0
                cmax(i,j,n) = crse(i,j,n)
                cmin(i,j,n) = crse(i,j,n)
                do joff = -1,1
                do ioff = -1,1
                  cmax(i,j,n) = max(cmax(i,j,n),crse(i+ioff,j+joff,n))
                  cmin(i,j,n) = min(cmin(i,j,n),crse(i+ioff,j+joff,n))
                end do
                end do
             end do
          end do

       end do

! ...  Compute unlimited and limited slopes

       do n = 1, nvar

          do j=cslopelo(2), cslopehi(2)
             do i=cslopelo(1), cslopehi(1)
                uc_xslope(i,j,n) = half*(crse(i+1,j,n)-crse(i-1,j,n))
                cen  = uc_xslope(i,j,n)
                forw = two*(crse(i+1,j,n)-crse(i,j,n))
                back = two*(crse(i,j,n)-crse(i-1,j,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,zero,forw*back>=zero)
                lc_xslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))
             end do
          end do

          if (bclo(1,n) .eq. EXT_DIR .or. bclo(1,n).eq.HOEXTRAP) then
            i = cslopelo(1)
            if (xok) then
                do j=cslopelo(2), cslopehi(2)
                   uc_xslope(i,j,n)  = -sixteen/fifteen*crse(i-1,j,n) &
                       + half*crse(i,j,n) &
                        + two3rd*crse(i+1,j,n) - tenth*crse(i+2,j,n)
                end do
            else
                do j=cslopelo(2), cslopehi(2)
                   uc_xslope(i,j,n)  = fourth * ( &
                     crse(i+1,j,n) + five*crse(i,j,n) - six*crse(i-1,j,n) )
                end do
            endif
            do j=cslopelo(2), cslopehi(2)
               cen  = uc_xslope(i,j,n)
               forw = two*(crse(i+1,j,n)-crse(i,j,n))
               back = two*(crse(i,j,n)-crse(i-1,j,n))
               slp  = min(abs(forw),abs(back))
               slp  = merge(slp,zero,forw*back>=zero)
               lc_xslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))
            end do
          end if

          if (bchi(1,n) .eq. EXT_DIR .or. bchi(1,n).eq.HOEXTRAP) then
            i = cslopehi(1)
            if (xok) then
                do j=cslopelo(2), cslopehi(2)
                   uc_xslope(i,j,n) = sixteen/fifteen*crse(i+1,j,n) &
                        - half*crse(i,j,n) &
                        - two3rd*crse(i-1,j,n) + tenth*crse(i-2,j,n)
                end do
            else
                do j=cslopelo(2), cslopehi(2)
                   uc_xslope(i,j,n) = -fourth * ( &
                     crse(i-1,j,n) + five*crse(i,j,n) - six*crse(i+1,j,n) )
                end do
            endif
            do j=cslopelo(2), cslopehi(2)
               cen  = uc_xslope(i,j,n)
               forw = two*(crse(i+1,j,n)-crse(i,j,n))
               back = two*(crse(i,j,n)-crse(i-1,j,n))
               slp  = min(abs(forw),abs(back))
               slp  = merge(slp,zero,forw*back>=zero)
               lc_xslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))
            end do
          end if

          do j=cslopelo(2), cslopehi(2)
             do i=cslopelo(1), cslopehi(1)
                uc_yslope(i,j,n) = half*(crse(i,j+1,n)-crse(i,j-1,n))
                cen  = uc_yslope(i,j,n)
                forw = two*(crse(i,j+1,n)-crse(i,j,n))
                back = two*(crse(i,j,n)-crse(i,j-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,zero,forw*back>=zero)
                lc_yslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))
             end do
          end do

          if (bclo(2,n) .eq. EXT_DIR .or. bclo(2,n).eq.HOEXTRAP) then
             j = cslopelo(2)
             if (yok) then
                do i=cslopelo(1), cslopehi(1)
                   uc_yslope(i,j,n)  = -sixteen/fifteen*crse(i,j-1,n) &
                        + half*crse(i,j,n) &
                        + two3rd*crse(i,j+1,n) - tenth*crse(i,j+2,n)
                end do
             else
                do i=cslopelo(1), cslopehi(1)
                   uc_yslope(i,j,n)  = fourth * ( &
                     crse(i,j+1,n) + five*crse(i,j,n) - six*crse(i,j-1,n) )
                end do
             endif
             do i=cslopelo(1), cslopehi(1)
                cen  = uc_yslope(i,j,n)
                forw = two*(crse(i,j+1,n)-crse(i,j,n))
                back = two*(crse(i,j,n)-crse(i,j-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,zero,forw*back>=zero)
                lc_yslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))
             end do
          end if

          if (bchi(2,n) .eq. EXT_DIR .or. bchi(2,n).eq.HOEXTRAP) then
             j = cslopehi(2)
             if (yok) then
                do i=cslopelo(1), cslopehi(1)
                   uc_yslope(i,j,n) = sixteen/fifteen*crse(i,j+1,n) &
                        - half*crse(i,j,n) &
                        - two3rd*crse(i,j-1,n) + tenth*crse(i,j-2,n)
                end do
             else
                do i=cslopelo(1), cslopehi(1)
                   uc_yslope(i,j,n) = -fourth * ( &
                     crse(i,j-1,n) + five*crse(i,j,n) - six*crse(i,j+1,n) )
                end do
             endif
             do i=cslopelo(1), cslopehi(1)
                cen  = uc_yslope(i,j,n)
                forw = two*(crse(i,j+1,n)-crse(i,j,n))
                back = two*(crse(i,j,n)-crse(i,j-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,zero,forw*back>=zero)
                lc_yslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))
             end do
          end if

       end do

       if (lim_slope.eq.0) then

! ...    Do the interpolation using unlimited slopes.

          do n = 1, nvar
             do j = fblo(2), fbhi(2)
                jc = IX_PROJ(j,lratioy)
                do i = fblo(1), fbhi(1)
                   ic = IX_PROJ(i,lratiox)
                   fine(i,j,n) = crse(ic,jc,n) &
                        + voffx(i)*uc_xslope(ic,jc,n) &
                        + voffy(j)*uc_yslope(ic,jc,n)
                end do
             end do
          end do

       else 

         if (lin_limit.eq.1) then

! ...      compute linear limited slopes
!          Note that the limited and the unlimited slopes
!          have the same sign, and it is assumed that they do.
!
! ... --> compute slope factors

           do j=cslopelo(2), cslopehi(2)
             do i=cslopelo(1), cslopehi(1)
                xslope_factor(i,j) = one
                yslope_factor(i,j) = one
             end do
           end do

           do n = 1, nvar
             do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                   denom = uc_xslope(i,j,n)
                   denom = merge(denom,one,denom.ne.zero)
                   factorn = lc_xslope(i,j,n)/denom
                   factorn = merge(one,factorn,denom.eq.zero)
                   xslope_factor(i,j) = min(xslope_factor(i,j),factorn)
                   denom = uc_yslope(i,j,n)
                   denom = merge(denom,one,denom.ne.zero)
                   factorn = lc_yslope(i,j,n)/denom
                   factorn = merge(one,factorn,denom.eq.zero)
                   yslope_factor(i,j) = min(yslope_factor(i,j),factorn)
                end do
             end do
           end do

! ... -->  compute linear limited slopes

           do n = 1, nvar
             do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                   lc_xslope(i,j,n) = xslope_factor(i,j)*uc_xslope(i,j,n)
                   lc_yslope(i,j,n) = yslope_factor(i,j)*uc_yslope(i,j,n)
                end do
             end do
           end do

         else

!          Limit slopes so as to not introduce new maxs or mins.


            do n = 1, nvar
               do j = voff_lo(2),voff_hi(2)
                  jc = IX_PROJ(j,lratioy)

                  do i = voff_lo(1),voff_hi(1)
                     ic = IX_PROJ(i,lratiox)

                     orig_corr_fact = voffx(i)*lc_xslope(ic,jc,n) &
                          + voffy(j)*lc_yslope(ic,jc,n) 
                     dummy_fine = crse(ic,jc,n) + orig_corr_fact
                     if ( (dummy_fine .gt. cmax(ic,jc,n)) .and. &
                          (abs(orig_corr_fact) .gt. 1.e-10*abs(crse(ic,jc,n)))) then
                        corr_fact = (cmax(ic,jc,n) - crse(ic,jc,n)) / orig_corr_fact
                        alpha(ic,jc,n) = min(alpha(ic,jc,n),corr_fact)
                     endif
                     if ( (dummy_fine .lt. cmin(ic,jc,n)) .and. &
                          (abs(orig_corr_fact) .gt. 1.e-10*abs(crse(ic,jc,n)))) then
                        corr_fact = (cmin(ic,jc,n) - crse(ic,jc,n)) / orig_corr_fact
                        alpha(ic,jc,n) = min(alpha(ic,jc,n),corr_fact)
                     endif

#ifndef NDEBUG
                     if (alpha(ic,jc,n) .lt. 0.d0) then
                        print *,'OOPS - ALPHA SHOULD BE POSITIVE IN CCINTERP '
                        print *,'ALPHA = ',alpha(ic,jc,n)
                        print *,'AT (I,J,N) = ',ic,jc,n
                        print *,'ORIG_CORR_FACT = ',orig_corr_fact
                        call bl_abort(" ")
                     endif
                     if (alpha(ic,jc,n) .gt. 1.d0) then
                        print *,'OOPS - ALPHA SHOULD BE <= 1.0 IN CCINTERP '
                        print *,'ALPHA = ',alpha(ic,jc,n)
                        print *,'AT (I,J,N) = ',ic,jc,n
                        print *,'ORIG_CORR_FACT = ',orig_corr_fact
                        call bl_abort(" ")
                     endif
#endif
                  end do
               end do
            end do

         end if

! ...    Do the interpolation with limited slopes.

          do n = 1, nvar
            do j = fblo(2), fbhi(2)
               jc = IX_PROJ(j,lratioy)
               do i = fblo(1), fbhi(1)
                  ic = IX_PROJ(i,lratiox)
                  fine(i,j,n) = crse(ic,jc,n) + alpha(ic,jc,n)* &
                     ( voffx(i)*lc_xslope(ic,jc,n) &
                      +voffy(j)*lc_yslope(ic,jc,n) )
               end do
            end do
          end do

       end if

     end subroutine FORT_LINCCINTERP

     subroutine FORT_CQINTERP (fine, fine_l1,fine_l2,fine_h1,fine_h2, &
                               fb_l1, fb_l2, fb_h1, fb_h2, &
                               nvar, lratiox, lratioy, crse, clo, chi, &
                               cb_l1, cb_l2, cb_h1, cb_h2, &
                               fslo, fshi, cslope, clen, fslope, fdat, &
                               flen, voff, bc, limslope, &
                               fvcx, fvcy, cvcx, cvcy, &
                               actual_comp,actual_state)

      implicit none

      integer fine_l1,fine_l2,fine_h1,fine_h2
      integer fslo(2), fshi(2)
      integer fb_l1, fb_l2, fb_h1, fb_h2
      integer cb_l1, cb_l2, cb_h1, cb_h2
      integer clo, chi
      integer lratiox, lratioy, nvar, clen, flen, limslope
      integer bc(2,2,nvar)
      integer actual_comp,actual_state
      REAL_T fine(fine_l1:fine_h1,fine_l2:fine_h2,nvar)
      REAL_T crse(clo:chi, nvar)
      REAL_T cslope(clo:chi, 5)
      REAL_T fslope(flen, 5)
      REAL_T fdat(flen)
      REAL_T voff(flen)
      REAL_T fvcx(fb_l1:fb_h1+1)
      REAL_T fvcy(fb_l2:fb_h2+1)
      REAL_T cvcx(cb_l1:cb_h1+1)
      REAL_T cvcy(cb_l2:cb_h2+1)

#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)

      integer n, fn
      integer i, ic, ioff
      integer j, jc, joff
      integer ist, jst
      REAL_T cen
      REAL_T fcen, ccen
      REAL_T diffxy,diffxx,diffyy
      REAL_T yoff
      integer ncbx, ncby
      integer ncsx
      integer jslo
      integer icc, istart, iend
      logical xok, yok

      ncbx = cb_h1-cb_l1+1
      ncby = cb_h2-cb_l2+1
      xok = (ncbx .ge. 2)
      yok = (ncby .ge. 2)
      ncsx = ncbx+2
      ist = 1
      jst = ncsx
      jslo = cb_l2-1

      do i = fb_l1, fb_h1 
         fn = i-fslo(1)+1
         ic = IX_PROJ(i,lratiox)
         fcen = half*(fvcx(i)+fvcx(i+1))
         ccen = half*(cvcx(ic)+cvcx(ic+1))
         voff(fn) = (fcen-ccen)/(cvcx(ic+1)-cvcx(ic))
      end do   

      do n = 1, nvar 
        do i = clo, chi 
          crse(i,n) = merge(crse(i,n),zero,abs(crse(i,n)).gt.1.0d-50)
        end do
      end do

      do 290 n = 1, nvar 

            do i = 1, clen 
               cen = half*(crse(i+ist,n)-crse(i-ist,n))
               diffxy = fourth*(crse(i+ist+jst,n)+crse(i-ist-jst,n) &
                               -crse(i-ist+jst,n)-crse(i+ist-jst,n))
               diffxx = crse(i+ist,n)-two*crse(i,n)+crse(i-ist,n)
               cslope(i,1)=cen
               cslope(i,3)=diffxx
               cslope(i,5)=diffxy
            end do
            if (xok) then
               if (bclo(1,n) .eq. EXT_DIR .or. bclo(1,n).eq.HOEXTRAP) then
                  do i = 1, clen, jst 
                     cen  = -sixteen/fifteen*crse(i-ist,n) + half*crse(i,n) &
                          + two3rd*crse(i+ist,n) - tenth*crse(i+2*ist,n)
                     cslope(i,1)=cen
                     cslope(i,3)=zero
                     cslope(i,5)=zero
                  end do
               end if
               if (bchi(1,n) .eq. EXT_DIR .or. bchi(1,n).eq.HOEXTRAP) then
                  do i = ncbx, clen, jst 
                     cen = sixteen/fifteen*crse(i+ist,n) - half*crse(i,n) &
                          - two3rd*crse(i-ist,n) + tenth*crse(i-2*ist,n)
                     cslope(i,1)=cen
                     cslope(i,3)=zero
                     cslope(i,5)=zero
                  end do
               end if
            end if

            do i = 1, clen 
               cen  = half*(crse(i+jst,n)-crse(i-jst,n))
               diffyy = crse(i+jst,n)-two*crse(i,n)+crse(i-jst,n)
               cslope(i,2)=cen
               cslope(i,4)=diffyy
            end do
            if (yok) then
               if (bclo(2,n) .eq. EXT_DIR .or. bclo(2,n).eq.HOEXTRAP) then
                  do i = 1, ncbx 
                     cen  = -sixteen/fifteen*crse(i-jst,n) + half*crse(i,n) &
                          + two3rd*crse(i+jst,n) - tenth*crse(i+2*jst,n)
                     cslope(i,2)=cen
                     cslope(i,4)=zero
                     cslope(i,5)=zero
                  end do
               end if
               if (bchi(2,n) .eq. EXT_DIR .or. bchi(2,n).eq.HOEXTRAP) then
                  do i = clen-ncbx,clen 
                     cen = sixteen/fifteen*crse(i+jst,n) - half*crse(i,n) &
                          - two3rd*crse(i-jst,n) + tenth*crse(i-2*jst,n)
                     cslope(i,2)=cen
                     cslope(i,4)=zero
                     cslope(i,5)=zero
                  end do
               end if
            end if

            do 360 jc = cb_l2, cb_h2 
               do 370 ioff = 1, lratiox 
                  icc = clo + ist + jst*(jc-jslo)
                  istart = ioff
                  iend = ioff + (ncbx-1)*lratiox
                  do 380 fn = istart, iend, lratiox 
                     fslope(fn,1) = cslope(icc,1)
                     fslope(fn,2) = cslope(icc,2)
                     fslope(fn,3) = cslope(icc,3)
                     fslope(fn,4) = cslope(icc,4)
                     fslope(fn,5) = cslope(icc,5)
                     fdat(fn) = crse(icc,n)
                     icc = icc + ist
380               continue
370            continue

               do 390 joff = 0, lratioy-1 
                  j = lratioy*jc + joff
                  if ((j.lt.fb_l2).or.(j.gt.fb_h2)) goto 390
                  fcen = half*(fvcy(j)+fvcy(j+1))
                  ccen = half*(cvcy(jc)+cvcy(jc+1))
                  yoff = (fcen-ccen)/(cvcy(jc+1)-cvcy(jc))

                  do 400 i = fb_l1, fb_h1 
                     fn = i-fslo(1)+1
                     fine(i,j,n) = fdat(fn) &
                          + voff(fn)              *fslope(fn,1) &
                          + yoff                  *fslope(fn,2) &
                          + half*voff(fn)*voff(fn)*fslope(fn,3) &
                          + half*yoff    *yoff    *fslope(fn,4) &
                          + voff(fn)*yoff         *fslope(fn,5)
400               continue
390            continue
360         continue

290   continue

    end subroutine FORT_CQINTERP
! ::: 
! ::: --------------------------------------------------------------
! ::: pcinterp:  cell centered piecewise constant interpolation
! ::: 
! ::: Inputs/Outputs
! ::: fine        <=>  (modify) fine grid array
! ::: flo,fhi      =>  (const)  index limits of fine grid
! ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
! ::: 
! ::: crse         =>  (const)  coarse grid data 
! ::: clo,chi      =>  (const)  index limits of coarse grid
! ::: cblo,cbhi    =>  (const) coarse grid region containing fblo,fbhi
! ::: 
! ::: longdir      =>  (const)  which index direction is longest (1 or 2)
! ::: lratio(2)    =>  (const)  refinement ratio between levels
! ::: nvar         =>  (const)  number of components in array
! ::: 
! ::: TEMPORARY ARRAYS
! ::: ftmp         =>  1-D temp array
! ::: --------------------------------------------------------------
! ::: 
    subroutine FORT_PCINTERP (crse,crse_l1,crse_l2,crse_h1,crse_h2,cblo,cbhi, &
                              fine,fine_l1,fine_l2,fine_h1,fine_h2,fblo,fbhi, &
                              longdir,lratiox,lratioy,nvar, &
                              ftmp,ftmp_lo,ftmp_hi, &
                              actual_comp,actual_state)

      implicit none

      integer crse_l1,crse_l2,crse_h1,crse_h2
      integer cblo(2), cbhi(2)
      integer fine_l1,fine_l2,fine_h1,fine_h2
      integer fblo(2), fbhi(2)
      integer ftmp_lo, ftmp_hi
      integer nvar, lratiox, lratioy, longdir
      integer actual_comp,actual_state
      REAL_T  crse(crse_l1:crse_h1,crse_l2:crse_h2, nvar)
      REAL_T  fine(fine_l1:fine_h1,fine_l2:fine_h2, nvar)
      REAL_T  ftmp(ftmp_lo:ftmp_hi)

      integer i, j, ic, jc, ioff, joff, n

      if (longdir .eq. 1) then
         do n = 1, nvar
         do jc = cblo(2), cbhi(2)
	    j = jc*lratioy
	    do ioff = 0, lratiox-1
	       do ic = cblo(1), cbhi(1)
	          i = lratiox*ic + ioff
	          ftmp(i) = crse(ic,jc,n)
               end do
	    end do
	    do joff = 0, lratioy-1
	       j = lratioy*jc + joff
	       if (j.ge.fblo(2).and.j.le.fbhi(2)) then
	          do i = fblo(1), fbhi(1)
		     fine(i,j,n) = ftmp(i)
		  end do
	       end if
	    end do
	 end do
	 end do
      else
         do n = 1, nvar
         do ic = cblo(1), cbhi(1)
	    i = ic*lratiox
	    do joff = 0, lratioy-1
	       do jc = cblo(2), cbhi(2)
	          j = lratioy*jc + joff
	          ftmp(j) = crse(ic,jc,n)
               end do
	    end do
	    do ioff = 0, lratiox-1
	       i = lratiox*ic + ioff
	       if (i.ge.fblo(1).and.i.le.fbhi(1)) then
	          do j = fblo(2), fbhi(2)
		     fine(i,j,n) = ftmp(j)
		  end do
	       end if
	    end do
	 end do
	 end do
      end if

    end subroutine FORT_PCINTERP

! ::: 
! ::: --------------------------------------------------------------
! ::: protect_interp:   redo interpolation if the result of linccinterp
! ::: generates under- or overshoots.
! ::: 
! ::: 
! ::: Inputs/Outputs
! ::: fine        <=>  (modify) fine grid array
! ::: flo,fhi      =>  (const)  index limits of fine grid
! ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
! ::: cblo,cbhi    =>  (const)  coarse equivalent of fblo,fbhi
! ::: nvar         =>  (const)  number of variables in state vector
! ::: lratio(3)    =>  (const)  refinement ratio between levels
! ::: 
! ::: crse         =>  (const)  coarse grid data widended by 1 zone
! ::: clo,chi      =>  (const)  index limits of crse grid
! :::
! ::: --------------------------------------------------------------
! ::: 

    subroutine FORT_PROTECT_INTERP (fine, fine_l1,fine_l2,fine_h1,fine_h2, fblo, fbhi, &
                                    crse, crse_l1,crse_l2,crse_h1,crse_h2, cblo, cbhi, &
                                    fvcx, fvcy, &
                                    fb_l1, fb_l2, fb_h1, fb_h2, &
                                    cvcx, cvcy, &
                                    cb_l1, cb_l2, cb_h1, cb_h2, &
                                    fine_state, state_l1,state_l2,state_h1,state_h2, &
                                    nvar, lratiox, lratioy, bc)

      implicit none

      integer fine_l1,fine_l2,fine_h1,fine_h2
      integer crse_l1,crse_l2,crse_h1,crse_h2
      integer state_l1,state_l2,state_h1,state_h2
      integer fblo(2), fbhi(2)
      integer cblo(2), cbhi(2)
      integer fb_l1, fb_l2, fb_h1, fb_h2
      integer cb_l1, cb_l2, cb_h1, cb_h2
      integer lratiox, lratioy, nvar
      integer bc(2,2,nvar)
      REAL_T fine(fine_l1:fine_h1,fine_l2:fine_h2,nvar)
      REAL_T crse(crse_l1:crse_h1,crse_l2:crse_h2, nvar)
      REAL_T fine_state(state_l1:state_h1,state_l2:state_h2, nvar)
      REAL_T fvcx(fb_l1:fb_h1)
      REAL_T fvcy(fb_l2:fb_h2)
      REAL_T cvcx(cb_l1:cb_h1)
      REAL_T cvcy(cb_l2:cb_h2)

      integer rMAX
      parameter (rMAX = 32)
      REAL_T alpha, sumN, sumP, negVal, posVal
      REAL_T crseTot, crseTotnew
      REAL_T orig_fine(0:rMAX-1,0:rMAX-1)
      REAL_T fvol,cvol
      integer redo_me
      integer ilo,ihi,jlo,jhi
      integer i,j,ic,jc,n
      integer icase

      if (MAX(lratiox,lratioy).gt.rMAX) then
         print *,'rMAX in INTERP_2D::FORT_PROTECT_INTERP must be >= ',MAX(lratiox,lratioy)
         call bl_abort(" ")
      endif

      do jc = cblo(2), cbhi(2)
      do ic = cblo(1), cbhi(1)

         ilo = max(lratiox*ic            ,fine_l1)
         ihi = min(lratiox*ic+(lratiox-1),fine_h1)
         jlo = max(lratioy*jc            ,fine_l2)
         jhi = min(lratioy*jc+(lratioy-1),fine_h2)

         do n = 2, nvar-1

            redo_me = 0
            do j = jlo,jhi
            do i = ilo,ihi
               if ((fine_state(i,j,n)+fine(i,j,n)) .lt. 0.d0) redo_me = 1
            enddo
            enddo

! ****************************************************************************************
!
!           If all the fine values are non-negative after the original interpolated 
!            correction, then we do nothing here.
!
!           If any of the fine values are negative after the original interpolated
!            correction, then we do our best.
!
!           Special cases:
!
!             1) Coarse correction > 0, and fine_state has some cells with 
!                negative values which will be filled before adding to the other cells.
!                Use the correction to bring negative cells to zero, then
!                distribute the remaining positive proportionally.
!
!             2) Coarse correction > 0, and correction can not make them all
!                positive.  Add correction only to the negative cells, in proportion
!                to their magnitude.
!
!             3) Coarse correction < 0, and fine_state DOES NOT have enough
!                  have enough positive state to absorb it.  Here we bring
!                  all the positive fine cells to zero then distribute the remaining
!                  negative amount in such a way as to make them all as close to the
!                  same negative value as possible.
!
!             4) Coarse correction < 0, fine_state has enough
!                  positive state to absorb it without making any fine 
!                  cells negative, BUT fine_state+fine is currently negative
!                  in at least one fine cell.  Here just take a constant percentage
!                  away from each positive and don't touch the negatives.
!
!             crseTot = volume-weighted sum of all interpolated values of the correction,
!                       which is equivalent to the total volume-weighted coarse correction
!             SumN = volume-weighted sum of all negative values of fine_state
!             SumP = volume-weighted sum of all positive values of fine_state
!
! ****************************************************************************************

            if (redo_me .eq. 1) then

               icase = 0

               do j = jlo,jhi
               do i = ilo,ihi
                  orig_fine(i-ilo,j-jlo) = fine(i,j,n)
               enddo
               enddo

               crseTot = 0.d0
               do j = jlo,jhi
               do i = ilo,ihi
                  fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                  crseTot = crseTot + fvol * fine(i,j,n)
               enddo
               enddo

               cvol = (cvcx(ic+1)-cvcx(ic)) * (cvcy(jc+1)-cvcy(jc))

               sumN = zero
               sumP = zero
               do j = jlo,jhi
               do i = ilo,ihi
                  fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                  if (fine_state(i,j,n) .le. 0.d0) then
                    sumN = SumN + fvol * fine_state(i,j,n)
                  else
                    sumP = sumP + fvol * fine_state(i,j,n)
                  endif
               enddo
               enddo

               if (crseTot .gt. 0.d0 .and. crseTot .ge. abs(sumN)) then
!              Here we want to fill in the negative values first, then add
!                the remaining positive proportionally.

                   icase = 1
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,n) .le. 0.d0) then
                        fine(i,j,n) = -fine_state(i,j,n)
                      endif
                   enddo
                   enddo

                   if (sumP > 0.d0) then

                     alpha = (crseTot - abs(sumN)) / sumP

                     do j = jlo,jhi
                     do i = ilo,ihi
                       if (fine_state(i,j,n) .ge. 0.d0) then
                         fine(i,j,n) = alpha * fine_state(i,j,n)
                       endif
                     enddo
                     enddo

                   else

                     posVal = (crseTot - abs(sumN)) / cvol

                     do j = jlo,jhi
                     do i = ilo,ihi
                       fine(i,j,n) = fine(i,j,n) + posVal
                     enddo
                     enddo

                   endif
            
                 endif

               if (crseTot .gt. 0.d0 .and. crseTot .lt. abs(sumN)) then
!              Here we don't have enough positive correction to fill all the
!                negative values of state, so we just try to fill them proportionally
!                and don't add any correction to the states already positive.

                   icase = 2
                   alpha = crseTot / abs(sumN)

                   do j = jlo,jhi
                   do i = ilo,ihi
                     if (fine_state(i,j,n) .lt. 0.d0) then
                       fine(i,j,n) = alpha * abs(fine_state(i,j,n))
                     else 
                       fine(i,j,n) = 0.d0
                     endif
                   enddo
                   enddo

               endif

               if (crseTot .lt. 0.d0 .and. abs(crseTot) .gt. sumP) then
!              Here we don't have enough positive states to absorb all the
!                negative correction, so we want to end up with all the fine
!                cells having the same negative value.

                   icase = 3
                   negVal = (sumP + sumN + crseTot)/cvol

                   do j = jlo,jhi
                   do i = ilo,ihi
                      fine(i,j,n) = negVal - fine_state(i,j,n)
                   enddo
                   enddo

               endif

               if (crseTot .lt. 0.d0 .and. abs(crseTot) .lt. sumP &
                                     .and. (sumP+sumN+crseTot) .gt. 0.d0) then
!              Here we have enough positive states to absorb all the
!                negative correction *and* redistribute to make negative cells
!                positive. 

                   icase = 4
                   alpha = (crseTot + sumN) / sumP

                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,n) .lt. 0.d0) then
                        fine(i,j,n) = -fine_state(i,j,n)
                      else
                        fine(i,j,n) = alpha * fine_state(i,j,n)
                      endif  
                   enddo
                   enddo

               endif

               if (crseTot .lt. 0.d0 .and. abs(crseTot) .lt. sumP &
                                     .and. (sumP+sumN+crseTot) .le. 0.d0) then
!              Here we have enough positive states to absorb all the
!                negative correction, but not to fix the states already negative. 
!                We bring all the positive states to zero, and use whatever 
!                remaining positiveness from the states to help the negative states.

                   icase = 5
                   alpha = (crseTot + sumP) / sumN

                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,n) .gt. 0.d0) then
                        fine(i,j,n) = -fine_state(i,j,n)
                      else 
                        fine(i,j,n) = alpha * fine_state(i,j,n)
                      endif
                   enddo
                   enddo

               endif

               crseTotnew   = 0.d0
               do j = jlo,jhi
               do i = ilo,ihi
                  fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                  crseTotnew   = crseTotnew   + fvol * fine(i,j,n)
               enddo
               enddo

               if (abs(crseTotnew - crseTot)/cvol .gt. 1.e-8) then
                  print *,' '
                  print *,'BLEW CONSERVATION with ICASE = ',icase
                  print *,'AT COARSE CELL ',ic,jc,' AND COMPONENT ',n
                  print *,'CRSETOT NEW OLD ',crseTotnew, crseTot
                  print *,'CVOL ',cvol
                  print *,'SUMP SUMN ',sumP,sumN
                  do j = jlo,jhi
                  do i = ilo,ihi
                     fvol = (fvcx(i+1)-fvcx(i)) * (fvcy(j+1)-fvcy(j))
                     print *,'FINE OLD NEW ',i,j,orig_fine(i-ilo,j-jlo), &
                                             fine(i,j,n), fine_state(i,j,n), &
                                             fvol
                     if (abs(fvol) .lt. 1.d-50) then
                       print *,'MAKING FVOL ',fvcx(i+1),fvcx(i),fvcy(j+1),fvcy(j)
                     endif
                  enddo
                  enddo
               endif

!              do j = jlo,jhi
!              do i = ilo,ihi
!                 if ((fine_state(i,j,n) + fine(i,j,n)) .lt. 0.d0) then
!                    print *,'STILL NEGATIVE AT ',i,j,n
!                    print *,'AT COARSE CELL ',ic,jc
!                    print *,'FINE STATE ',fine_state(i,j,n)
!                    print *,'FINE CORRECTION ',fine(i,j,n)
!                    print *,'CRSETOT ',crseTot
!                    print *,'SUMN / SUMP ',sumN, sumP
!                    print *,' '
!                 endif
!              enddo
!              enddo
!              enddo

!           End (if redo .eq. 1)
            endif

         enddo

!     Set sync for density (n=1) to sum of spec sync (2:nvar-1)
         do j = jlo,jhi
         do i = ilo,ihi
            fine(i,j,1) = 0.d0
            do n = 2,nvar-1
               fine(i,j,1) = fine(i,j,1) + fine(i,j,n)
            enddo
         enddo
         enddo

!     End of coarse index loops
      enddo
      enddo

    end subroutine FORT_PROTECT_INTERP

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
! ::: lratio[xy]   =>  (const)  refinement ratio between levels
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
! ::: --------------------------------------------------------------
! ::: 
     subroutine FORT_QUARTINTERP (fine, fine_l1,fine_l2,fine_h1,fine_h2, &
                                  fblo, fbhi, fb2lo, fb2hi, &
                                  crse, crse_l1,crse_l2,crse_h1,crse_h2, &
                                  cblo, cbhi, cb2lo, cb2hi, &
                                  nvar, &
                                  lratiox, lratioy, &
                                  ftmp, ctmp, &
                                  bc,actual_comp,actual_state)

       implicit none

       integer fine_l1,fine_l2,fine_h1,fine_h2
       integer crse_l1,crse_l2,crse_h1,crse_h2
       integer fblo(2), fbhi(2), fb2lo(2), fb2hi(2)
       integer cblo(2), cbhi(2), cb2lo(2), cb2hi(2)
       integer nvar,lratiox,lratioy
       integer bc(2,2,nvar)
       integer actual_comp,actual_state
       REAL_T fine(fine_l1:fine_h1,fine_l2:fine_h2,nvar)
       REAL_T crse(crse_l1:crse_h1,crse_l2:crse_h2,nvar)
       REAL_T ftmp(fb2lo(1):fb2hi(1))
       REAL_T ctmp(cblo(1):cbhi(1),0:lratioy-1)

!      Local variables
       integer i,j,ii,jj,n,iry
       REAL_T cL(-2:2)
!       REAL_T cR(-2:2)
       data cL/ -0.01171875D0,  0.0859375D0, 0.5d0, -0.0859375D0, &
                 0.01171875D0 /
!       data cR/  0.01171875D0, -0.0859375D0, 0.5d0,  0.0859375D0, &
!                -0.01171875D0 /
       
       if (lratiox.eq.2 .and. lratioy.eq.2) then
          do n = 1, nvar
             do j = cb2lo(2), cb2hi(2)
                do i = cblo(1), cbhi(1)
                   ctmp(i,0) = 2.d0*(cL(-2)*crse(i,j-2,n) &
                        +            cL(-1)*crse(i,j-1,n) &
                        +            cL( 0)*crse(i,j  ,n) &
                        +            cL( 1)*crse(i,j+1,n) &
                        +            cL( 2)*crse(i,j+2,n))
                   ctmp(i,1) = 2.d0*crse(i,j,n) - ctmp(i,0)
!$$$                   ctmp(i,1) = 2.d0*(cR(-2)*crse(i,j-2,n) & 
!$$$                        +            cR(-1)*crse(i,j-1,n) &
!$$$                        +            cR( 0)*crse(i,j  ,n) &
!$$$                        +            cR( 1)*crse(i,j+1,n) &
!$$$                        +            cR( 2)*crse(i,j+2,n))
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
!$$$                         ftmp(ii+1) = 2.d0*(cR(-2)*ctmp(i-2,iry) & 
!$$$                              +             cR(-1)*ctmp(i-1,iry) &
!$$$                              +             cR( 0)*ctmp(i  ,iry) &
!$$$                              +             cR( 1)*ctmp(i+1,iry) &
!$$$                              +             cR( 2)*ctmp(i+2,iry))
                      enddo
                      do ii = fblo(1), fbhi(1)
                         fine(ii,jj,n) = ftmp(ii)
                      enddo
                   endif
                enddo
             enddo
          enddo
       else if (lratiox.eq.4 .and. lratioy.eq.4) then
!      todo
          write(6,*) 'FORT_QUARTINTERP: refinement ratio = 4 TODO'
          stop
       else
          write(6,*) 'FORT_QUARTINTERP: unsupported refinement ratio'
          stop
       endif

     end subroutine FORT_QUARTINTERP
