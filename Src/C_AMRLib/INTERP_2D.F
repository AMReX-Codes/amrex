c
c $Id: INTERP_2D.F,v 1.11 2002-08-29 22:14:34 car Exp $
c

#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "INTERP_F.H"
#include <ArrayLim.H>

#define IX_PROJ(A,B) (A+B*iabs(A))/B-iabs(A)
#define SDIM 2


c ::: --------------------------------------------------------------
c ::: nbinterp:  node based bilinear interpolation
c :::
c ::: INPUTS/OUTPUTS
c ::: fine        <=>  (modify) fine grid array
c ::: DIMS(fine)   =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c :::
c ::: crse         =>  (const)  coarse grid data widened by 1 zone
c ::: DIMS(crse)   =>  (const)  index limits of coarse grid
c :::
c ::: lratio(3)    =>  (const)  refinement ratio between levels
c ::: nvar         =>  (const)  number of components in array
c ::: num_slp      =>  (const)  number of types of slopes
c :::
c ::: TEMPORARY ARRAYS
c ::: sl           =>  num_slp 1-D slope arrays
c ::: --------------------------------------------------------------
c :::
      subroutine FORT_NBINTERP (crse, DIMS(crse), DIMS(cb),
     $                          fine, DIMS(fine), DIMS(fb),
     $                          lratiox, lratioy, nvar,
     $                          sl, num_slp)

      integer DIMDEC(crse)
      integer DIMDEC(cb)
      integer DIMDEC(fine)
      integer DIMDEC(fb)
      integer lratiox, lratioy, nvar
      integer num_slp
      REAL_T  fine(DIMV(fine),nvar)
      REAL_T  crse(DIMV(crse),nvar)
      REAL_T  sl(DIM1(cb),num_slp)

#define  SLX 1
#define  SLY 2
#define  SLXY 3

c ::: local var
      integer lx, ly
      integer i, j, ifn, jfn, n
      integer ilo, ihi, jlo, jhi
      integer jstrtFine, jstopFine, istrtFine, istopFine

      REAL_T fx, fy
      REAL_T RX, RY, RXY
      REAL_T dx0, d0x, dx1
      REAL_T slope

      slope(i,j,n,fx,fy) = crse(i,j,n) +
     &                     fx*sl(i,SLX) + fy*sl(i,SLY) + fx*fy*sl(i,SLXY)

      RX = one/dfloat(lratiox)
      RY = one/dfloat(lratioy)
      RXY = RX*RY

c
c     NOTES:
c         1) (i, j) loop over the coarse cells
c         2) ?strtFine and ?stopFine are the beginning and ending fine cell
c            indices corresponding to the current coarse cell.  ?stopFine
c            is restricted for the last coarse cell in each direction since
c            for this cell we only need to do the face and not the fine nodes
c            inside this cell.
c         3) (lx, ly) as well as ?lo and ?hi refer to the fine node indices
c            as an offset from ?strtFine.
c
      do 100 n = 1, nvar
        do 120 j = ARG_L2(cb), ARG_H2(cb)
          jstrtFine = j * lratioy
          jstopFine = jstrtFine + lratioy - 1
          if (j .eq. ARG_H2(cb)) jstopFine = jstrtFine

          jlo = max(ARG_L2(fb),jstrtFine) - jstrtFine
          jhi = min(ARG_H2(fb),jstopFine) - jstrtFine


c
c         ::::: compute slopes :::::
c
c         NOTE: The IF logic in the calculation of the slopes is to
c               prevent stepping out of bounds on the coarse data when
c               computing the slopes on the ARG_H?(cb) cells.  These
c               slopes actually are not used since they are multiplied by
c               zero.
c
          do i = ARG_L1(cb), ARG_H1(cb)
            dx0 = zero
            if (i .NE. ARG_H1(cb)) dx0 = crse(i+1,j,n) - crse(i,j,n)

            d0x = zero
            if (j .NE. ARG_H2(cb)) d0x = crse(i,j+1,n) - crse(i,j,n)

            dx1 = zero
            if (i .NE. ARG_H1(cb) .and. j .NE. ARG_H2(cb))
     $        dx1 = crse(i+1,j+1,n) - crse(i,j+1,n)

            sl(i,SLX) = RX*dx0
            sl(i,SLY) = RY*d0x
            sl(i,SLXY) = RXY*(dx1 - dx0)
          end do

c
c         ::::: compute fine strip of interpolated data
c
          do ly = jlo, jhi
            jfn = lratioy * j + ly
            fy = dfloat(ly)

            do i = ARG_L1(cb), ARG_H1(cb)
              istrtFine = i * lratiox
              istopFine = istrtFine + lratiox - 1
              if (i .eq. ARG_H1(cb)) istopFine = istrtFine

              ilo = max(ARG_L1(fb),istrtFine) - istrtFine
              ihi = min(ARG_H1(fb),istopFine) - istrtFine

              do lx = ilo, ihi
                ifn = lratiox * i + lx
                fx = dfloat(lx)

                fine(ifn,jfn,n) = slope(i,j,n,fx,fy)
              end do
            end do
          end do

c       ::::: end of do j
120     continue
c     ::::: end of do n
100   continue

#undef  SLX
#undef  SLY
#undef  SLXY

      end


c ::: 
c ::: --------------------------------------------------------------
c ::: cbinterp:  cell centered bilinear interpolation
c ::: 
c ::: NOTE: it is assumed that the coarse grid array is
c ::: large enough to define interpolated values
c ::: in the region fblo:fbhi on the fine grid
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: DIMS(fine)   =>  (const)  index limits of fine grid
c ::: DIMS(fb)     =>  (const)  subregion of fine grid to get values
c ::: 
c ::: crse         =>  (const)  coarse grid data 
c ::: DIMS(crse)   =>  (const)  index limits of coarse grid
c ::: 
c ::: lratio(2)    =>  (const)  refinement ratio between levels
c ::: nvar         =>  (const)  number of components in array
c ::: 
c ::: TEMPORARY ARRAYS
c ::: slx,sly,slxy =>  1-D slope arrays
c ::: strip        =>  1-D temp array
c ::: --------------------------------------------------------------
c ::: 
      subroutine FORT_CBINTERP (crse, DIMS(crse), DIMS(cb),
     $                          fine, DIMS(fine), DIMS(fb),
     $		                lratiox, lratioy, nvar,
     $                          sl, num_slp, strip, strip_lo, strip_hi)

      integer DIMDEC(crse)
      integer DIMDEC(cb)
      integer DIMDEC(fine)
      integer DIMDEC(fb)
      integer lratiox, lratioy, nvar
      integer num_slp
      integer strip_lo, strip_hi
      REAL_T  fine(DIMV(fine), nvar)
      REAL_T  crse(DIMV(crse), nvar)
      REAL_T  sl(DIM1(cb),num_slp)
      REAL_T  strip(strip_lo:strip_hi)

#define SLX 1
#define SLY 2
#define SLXY 3

c ::: local var
      integer lx, ly, hratx, hraty, ic, jc, jfn, jfc, i, j, n
      REAL_T x, y, denomx, denomy

      denomx = one/dfloat(2*lratiox)
      denomy = one/dfloat(2*lratioy)

      hratx = lratiox/2
      hraty = lratioy/2

      do n = 1, nvar 
         do jc = ARG_L2(cb), ARG_H2(cb)-1 

c ::: first fill a strip that will fit

            do ic = ARG_L1(cb), ARG_H1(cb)-1
               sl(ic,SLX) = crse(ic+1,jc,n)-crse(ic,jc,n)
               sl(ic,SLY) = crse(ic,jc+1,n)-crse(ic,jc,n)
               sl(ic,SLXY) = crse(ic+1,jc+1,n)-crse(ic+1,jc,n)
     $              - crse(ic  ,jc+1,n)+crse(ic  ,jc,n)
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
                        strip(i) = crse(ic,jc,n) + x*sl(ic,SLX) +
     $                             y*sl(ic,SLY) + x*y*sl(ic,SLXY)
                     end do
                  end do

c ::: stuff into output array

                  do i = ARG_L1(fb), ARG_H1(fb) 
                     fine(i,jfc,n) = strip(i-hratx)
                  end do
               end if
            end do
         end do
      end do

      end

#undef  SLX
#undef  SLY
#undef  SLXY

c ::: 
c ::: --------------------------------------------------------------
c ::: ccinterp:   conservative interpolation from coarse grid to
c ::: subregion of fine grid defined by (fblo,fbhi)
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: flo,fhi      =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: nvar         =>  (const)  number of variables in state vector
c ::: lratio(2)    =>  (const)  refinement ratio between levels
c ::: 
c ::: crse         =>  (const)  coarse grid data widended by 1 zone
c ::: and unrolled
c ::: clo,chi      =>  (const)  one dimensional limits of crse grid
c ::: cslo,cshi    =>  (const)  coarse grid index limits where
c :::				slopes are to be defined. This is
c :::				the projection of (fblo,fbhi) down
c :::				to the coarse level 
c ::: fslo,fshi    =>  (const)  fine grid index limits where
c :::				slopes are needed.  This is the
c :::				refinement of (cslo,cshi) and
c :::				contains but may not be identical
c :::				to (fblo,fbhi).
c ::: cslope       =>  (modify) temp array coarse grid slopes
c ::: clen         =>  (const)  length of coarse gtid slopes
c ::: fslope       =>  (modify) temp array for fine grid slope
c ::: flen         =>  (const)  length of fine grid slope array
c ::: fdat         =>  (const)  temp array for fine grid data
c ::: limslope     =>  (const)  != 0 => limit slopes
c :::
c ::: NOTE: data must be sent in so that 
c :::	    cslope(1,*) and crse(1,*) are associated with
c :::	    the same cell
c :::
c ::: EXAMPLE:
c ::: Suppose the patch called "fine" has index extent:
c ::: 
c ::: floi1 = 3, fhii1 = 12
c ::: floi2 = 8, fhii2 = 20
c ::: 
c ::: suppose the subergion of this patch that is to be filled 
c ::: by interpolation has index extent:
c ::: 
c ::: fblo(1) = 5, fbhi(1) = 10
c ::: fblo(2) = 13, fbhi(2) = 20
c ::: 
c ::: suppose the refinement ratio is 2
c ::: 
c ::: Then the coarsening of this subregion (to level 0) is
c ::: 
c ::: cb_l1 = 2  cb_h1 = 5         (ncbx = 4)
c ::: cb_l2 = 6  cb_h2 = 10        (ncby = 5)
c ::: 
c ::: In order to compute slopes, we need one extra row of
c ::: coarse grid zones:
c ::: 
c ::: cslo(1) = 1  cshi(1) = 6         (ncsx = 6)
c ::: cslo(2) = 5  cshi(2) = 11        (ncsy = 7)
c ::: 
c ::: This is the size of the coarse grid array of data that filpatch 
c ::: has filled at level 0.
c ::: The "cslope" and "crse" arrays are this size.
c ::: 
c ::: In order to unroll the slope calculation we make these arrays look
c ::: like 1-D arrays.  The mapping from 2-D to 1-D is as fillows:
c ::: 
c ::: The point (cb_l(1),cb_l(2)) -> 1
c ::: The point (cslo(1),cslo(2)) -> clo = 1 - 1 - ncsx = -6
c ::: 
c ::: The point (cb_h1,cb_h2) -> clen = ncby*ncsx - 2 = 5*6-2 = 28
c ::: The point (cshi(1),cshi(2)) -> chi = clo + ncsx*ncsy - 1 
c :::                                    =  -6 +    6*7    - 1 = 35
c ::: 
c :::      -------------------------------------------------
c :::      |       |       |       |       |       |  chi  |  
c :::  11  |   30  |   31  |   32  |   33  |   34  |   35  |   cshi(2)
c :::      |       |       |       |       |       |       |
c :::      -------------------------------------------------
c :::      |       |       |       |       |  clen |       |  
c :::  10  |   24  |   25  |   26  |   27  |   28  |   29  |   cb_h(2)
c :::      |       |       |       |       |       |       |
c :::      -------------------------------------------------
c :::      |       |       |       |       |       |       |  
c :::   9  |   18  |   19  |   20  |   21  |   22  |   23  |  
c :::      |       |       |       |       |       |       |
c :::      -------------------------------------------------
c :::      |       |       |       |       |       |       |  
c :::   8  |   12  |   13  |   14  |   15  |   16  |   17  |  
c :::      |       |       |       |       |       |       |
c :::      -------------------------------------------------
c :::      |       |       |       |       |       |       |  
c :::   7  |    6  |    7  |    8  |    9  |   10  |   11  |  
c :::      |       |       |       |       |       |       |
c :::      -------------------------------------------------
c :::      |       |       |       |       |       |       |  
c :::   6  |    0  |    1  |    2  |    3  |    4  |    5  |   cb_l(2)
c :::      |       |       |       |       |       |       |
c :::      -------------------------------------------------
c :::      |  clo  |       |       |       |       |       |  
c :::   5  |   -6  |   -5  |   -4  |   -3  |   -2  |   -1  |   cslo(2)
c :::      |       |       |       |       |       |       |
c :::      -------------------------------------------------
c :::          1       2       3       4       5       6
c :::               cb_l1                   cb_h1
c :::       cslo(1)                                 cshi(1)
c ::: 
c ::: 
c ::: In the 1-D coordinates:
c :::    ist = 1    = stride in I direction
c :::    jst = 6    = stride in J direction  (ncsx)
c ::: 
c ::: --------------------------------------------------------------
c ::: 
      subroutine FORT_CCINTERP (fine, DIMS(fine), 
     $                          fb_l1, fb_l2, fb_h1, fb_h2,
     $                          nvar, lratiox, lratioy, crse, clo, chi, 
     $                          cb_l1, cb_l2, cb_h1, cb_h2,
     $		                fslo, fshi, cslope, clen, fslope, fdat,
     $                          flen, voff, bc, limslope,
     $                          fvcx, fvcy, cvcx, cvcy)

      integer DIMDEC(fine)
      integer fslo(2), fshi(2)
      integer fb_l1, fb_l2, fb_h1, fb_h2
      integer cb_l1, cb_l2, cb_h1, cb_h2
      integer clo, chi
      integer lratiox, lratioy, nvar, clen, flen, limslope
      integer bc(2,2,nvar)
      REAL_T fine(DIMV(fine),nvar)
      REAL_T crse(clo:chi, nvar)
      REAL_T cslope(clo:chi, 2)
      REAL_T fslope(flen, 2)
      REAL_T fdat(flen)
      REAL_T voff(flen)
      REAL_T fvcx(fb_l1:fb_h1+1)
      REAL_T fvcy(fb_l2:fb_h2+1)
      REAL_T cvcx(cb_l1:cb_h1+1)
      REAL_T cvcy(cb_l2:cb_h2+1)

#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)

c ::: local var
      integer n, fn
      integer i, ic, ioff
      integer j, jc, joff
      integer ist, jst
      REAL_T cen, forw, back, slp, sgn
      REAL_T fcen, ccen
      REAL_T xoff, yoff
      integer ncbx, ncby
      integer ncsx, ncsy
      integer islo, jslo
      integer icc, istart, iend
      logical xok, yok

      ncbx = cb_h1-cb_l1+1
      ncby = cb_h2-cb_l2+1
      xok = (ncbx .ge. 2)
      yok = (ncby .ge. 2)
      ncsx = ncbx+2
      ncsy = ncby+2
      ist = 1
      jst = ncsx
      islo = cb_l1-1
      jslo = cb_l2-1

      do i = fb_l1, fb_h1 
         fn = i-fslo(1)+1
         ic = IX_PROJ(i,lratiox)
         fcen = half*(fvcx(i)+fvcx(i+1))
         ccen = half*(cvcx(ic)+cvcx(ic+1))
         voff(fn) = (fcen-ccen)/(cvcx(ic+1)-cvcx(ic))
      end do

c ::: ::::: added to prevent underflow for small crse values
      do n = 1, nvar 
        do i = clo, chi 
          crse(i,n) = cvmgt(crse(i,n),zero,abs(crse(i,n)).gt.1.0D-20)
        end do
      end do

      do 290 n = 1, nvar 

c ::: ::::: compute slopes in x direction
         if (limslope .ne. 0) then
            do i = 1, clen 
               cen = half*(crse(i+ist,n)-crse(i-ist,n))
               forw = two*(crse(i+ist,n)-crse(i,n))
               back = two*(crse(i,n)-crse(i-ist,n))
               slp  = min(abs(forw),abs(back))
               slp  = cvmgp(slp,zero,forw*back)
               cslope(i,1)=sign(one,cen)*min(slp,abs(cen))
            end do
            if (xok) then
               if (bclo(1,n) .eq. EXT_DIR .or. bclo(1,n).eq.HOEXTRAP) then
                  do i = 1, clen, jst 
                     cen  = -sixteen/fifteen*crse(i-ist,n) + half*crse(i,n)
     $                    + two3rd*crse(i+ist,n) - tenth*crse(i+2*ist,n)
                     sgn  = sign(one,cen)
                     forw = two*(crse(i+ist,n)-crse(i,n))
                     back = two*(crse(i,n)-crse(i-ist,n))
                     slp  = min(abs(forw),abs(back))
                     slp  = cvmgp(slp,zero,forw*back)
                     cslope(i,1)=sgn*min(slp,abs(cen))
                  end do
               end if
               if (bchi(1,n) .eq. EXT_DIR .or. bchi(1,n).eq.HOEXTRAP) then
                  do i = ncbx, clen, jst 
                     cen = sixteen/fifteen*crse(i+ist,n) - half*crse(i,n)
     $                    - two3rd*crse(i-ist,n) + tenth*crse(i-2*ist,n)
                     sgn  = sign(one,cen)
                     forw = two*(crse(i+ist,n)-crse(i,n))
                     back = two*(crse(i,n)-crse(i-ist,n))
                     slp  = min(abs(forw),abs(back))
                     slp  = cvmgp(slp,zero,forw*back)
                     cslope(i,1)=sgn*min(slp,abs(cen))
                  end do
               end if
            end if
         else
            do i = 1, clen 
               cen = half*(crse(i+ist,n)-crse(i-ist,n))
               cslope(i,1)=cen
            end do
            if (xok) then
               if (bclo(1,n) .eq. EXT_DIR .or. bclo(1,n).eq.HOEXTRAP) then
                  do i = 1, clen, jst 
                     cen  = -sixteen/fifteen*crse(i-ist,n) + half*crse(i,n)
     $                    + two3rd*crse(i+ist,n) - tenth*crse(i+2*ist,n)
                     cslope(i,1)=cen
                  end do
               end if
               if (bchi(1,n) .eq. EXT_DIR .or. bchi(1,n).eq.HOEXTRAP) then
                  do i = ncbx, clen, jst 
                     cen = sixteen/fifteen*crse(i+ist,n) - half*crse(i,n)
     $                    - two3rd*crse(i-ist,n) + tenth*crse(i-2*ist,n)
                     cslope(i,1)=cen
                  end do
               end if
            end if
         end if

c ::: ::::: compute slopes in y direction
         if (limslope .ne. 0) then
            do i = 1, clen 
               cen  = half*(crse(i+jst,n)-crse(i-jst,n))
               forw = two*(crse(i+jst,n)-crse(i,n))
               back = two*(crse(i,n)-crse(i-jst,n))
               slp  = min(abs(forw),abs(back))
               slp  = cvmgp(slp,zero,forw*back)
               cslope(i,2)=sign(one,cen)*min(slp,abs(cen))
            end do
            if (yok) then
               if (bclo(2,n) .eq. EXT_DIR .or. bclo(2,n).eq.HOEXTRAP) then
                  do i = 1, ncbx 
                     cen  = -sixteen/fifteen*crse(i-jst,n) + half*crse(i,n)
     $                    + two3rd*crse(i+jst,n) - tenth*crse(i+2*jst,n)
                     sgn  = sign(one,cen)
                     forw = two*(crse(i+jst,n)-crse(i,n))
                     back = two*(crse(i,n)-crse(i-jst,n))
                     slp  = min(abs(forw),abs(back))
                     slp  = cvmgp(slp,zero,forw*back)
                     cslope(i,2)=sgn*min(slp,abs(cen))
                  end do
               end if
               if (bchi(2,n) .eq. EXT_DIR .or. bchi(2,n).eq.HOEXTRAP) then
                  do i = clen-ncbx,clen 
                     cen = sixteen/fifteen*crse(i+jst,n) - half*crse(i,n)
     $                    - two3rd*crse(i-jst,n) + tenth*crse(i-2*jst,n)
                     sgn  = sign(one,cen)
                     forw = two*(crse(i+jst,n)-crse(i,n))
                     back = two*(crse(i,n)-crse(i-jst,n))
                     slp  = min(abs(forw),abs(back))
                     slp  = cvmgp(slp,zero,forw*back)
                     cslope(i,2)=sgn*min(slp,abs(cen))
                  end do
               end if
            end if
         else
            do i = 1, clen 
               cen  = half*(crse(i+jst,n)-crse(i-jst,n))
               cslope(i,2)=cen
            end do
            if (yok) then
               if (bclo(2,n) .eq. EXT_DIR .or. bclo(2,n).eq.HOEXTRAP) then
                  do i = 1, ncbx 
                     cen  = -sixteen/fifteen*crse(i-jst,n) + half*crse(i,n)
     $                    + two3rd*crse(i+jst,n) - tenth*crse(i+2*jst,n)
                     cslope(i,2)=cen
                  end do
               end if
               if (bchi(2,n) .eq. EXT_DIR .or. bchi(2,n).eq.HOEXTRAP) then
                  do i = clen-ncbx,clen 
                     cen = sixteen/fifteen*crse(i+jst,n) - half*crse(i,n)
     $                    - two3rd*crse(i-jst,n) + tenth*crse(i-2*jst,n)
                     cslope(i,2)=cen
                  end do
               end if
            end if
         end if

         do 360 jc = cb_l2, cb_h2 
c           ::::: strip out a fine grid slope vector
            do 370 ioff = 1, lratiox 
               icc = clo + ist + jst*(jc-jslo)
               istart = ioff
               iend = ioff + (ncbx-1)*lratiox
               do 380 fn = istart, iend, lratiox 
                  fslope(fn,1) = cslope(icc,1)
                  fslope(fn,2) = cslope(icc,2)
                  fdat(fn) = crse(icc,n)
                  icc = icc + ist
380            continue
370         continue

            do 390 joff = 0, lratioy-1 
               j = lratioy*jc + joff
               if ((j.lt.fb_l2).or.(j.gt.fb_h2)) goto 390
               fcen = half*(fvcy(j)+fvcy(j+1))
               ccen = half*(cvcy(jc)+cvcy(jc+1))
               yoff = (fcen-ccen)/(cvcy(jc+1)-cvcy(jc))

               do 400 i = fb_l1, fb_h1 
                  fn = i-fslo(1)+1
                  fine(i,j,n) = fdat(fn) + voff(fn)*fslope(fn,1)
     &                                   +     yoff*fslope(fn,2)
400            continue
390         continue
360      continue

290   continue

      end
c ::: 
c ::: --------------------------------------------------------------
c ::: linccinterp:   linear conservative interpolation from coarse grid to
c ::: subregion of fine grid defined by (fblo,fbhi)
c ::: 
c ::: The interpolation is linear in that it uses a
c ::: a limiting scheme that preserves the value of 
c ::: any linear combination of the
c ::: coarse grid data components--e.g.,
c ::: if sum_ivar a(ic,jc,ivar)*fab(ic,jc,ivar) = 0, then
c ::: sum_ivar a(ic,jc,ivar)*fab(if,jf,ivar) = 0 is satisfied
c ::: in all fine cells if,jf covering coarse cell ic,jc.
c ::: 
c ::: If lin_limit = 0, the interpolation scheme is identical to
c ::: the used in ccinterp for limslope=1; the results should
c ::: be exactly the same -- difference = hard 0.
c :::
c ::: Unlike FORT_CCINTERP, this routine does not do any clever unrolling
c ::: and it does not use any 1-d strip--all calculations are done
c ::: on full 2-d arrays.  The onlu concession to vectorization
c ::: is that the innermost loops are longest.
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: flo,fhi      =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: nvar         =>  (const)  number of variables in state vector
c ::: lratio(2)    =>  (const)  refinement ratio between levels
c ::: 
c ::: crse         =>  (const)  coarse grid data widended by 1 zone
c ::: clo,chi      =>  (const)  index limits of crse grid
c ::: cslo,cshi    =>  (const)  coarse grid index limits where
c :::				slopes are to be defined. This is
c :::				the projection of (fblo,fbhi) down
c :::				to the coarse level 
c ::: ucslope      =>  (modify) temp array of unlimited coarse grid slopes
c ::: lcslope      =>  (modify) temp array of limited coarse grid slopes
c ::: slope_factor =>  (modify) temp array of slope limiting factors
c ::: lin_limit    =>  (const)  != 0 => do linear slope limiting scheme
c :::
c ::: --------------------------------------------------------------
c ::: 
      subroutine FORT_LINCCINTERP (fine, DIMS(fine), fblo, fbhi, 
     &                          DIMS(fvcb), 
     &                          crse, DIMS(crse), DIMS(cvcb),
     &                          uc_xslope, lc_xslope, xslope_factor,
     &                          uc_yslope, lc_yslope, yslope_factor,
     &                          DIMS(cslope),
     &                          cslopelo, cslopehi,
     $                          nvar, lratiox, lratioy, 
     $                          bc, lin_limit,
     $                          fvcx, fvcy, cvcx, cvcy,
     &                          voffx, voffy)

      implicit none

      integer DIMDEC(fine)
      integer DIMDEC(crse)
      integer DIMDEC(fvcb)
      integer DIMDEC(cvcb)
      integer DIMDEC(cslope)
      integer fblo(2), fbhi(2)
      integer cslopelo(2), cslopehi(2)
      integer lratiox, lratioy, nvar, lin_limit
      integer bc(2,2,nvar)
      REAL_T fine(DIMV(fine),nvar)
      REAL_T crse(DIMV(crse), nvar)
      REAL_T uc_xslope(DIMV(cslope),nvar)
      REAL_T lc_xslope(DIMV(cslope),nvar)
      REAL_T xslope_factor(DIMV(cslope))
      REAL_T uc_yslope(DIMV(cslope),nvar)
      REAL_T lc_yslope(DIMV(cslope),nvar)
      REAL_T yslope_factor(DIMV(cslope))
      REAL_T fvcx(DIM1(fvcb))
      REAL_T fvcy(DIM2(fvcb))
      REAL_T voffx(DIM1(fvcb))
      REAL_T voffy(DIM2(fvcb))
      REAL_T cvcx(DIM1(cvcb))
      REAL_T cvcy(DIM2(cvcb))

#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)

c ::: local var
      integer n 
      integer i, ic
      integer j, jc
      REAL_T cen, forw, back, slp, sgn
      REAL_T factorn, denom, slope_factor, dummy
      REAL_T fxcen, cxcen, fycen, cycen
      logical xok, yok
      integer ncbx, ncby

      ncbx = cslopehi(1)-cslopelo(1)+1
      ncby = cslopehi(2)-cslopelo(2)+1

      xok = (ncbx .ge. 2)
      yok = (ncby .ge. 2)

      do j = fblo(2), fbhi(2)
        jc = IX_PROJ(j,lratioy)
        fycen = half*(fvcy(j)+fvcy(j+1))
        cycen = half*(cvcy(jc)+cvcy(jc+1))
        voffy(j) = (fycen-cycen)/(cvcy(jc+1)-cvcy(jc))
      end do
      do i = fblo(1), fbhi(1)
         ic = IX_PROJ(i,lratiox)
         fxcen = half*(fvcx(i)+fvcx(i+1))
         cxcen = half*(cvcx(ic)+cvcx(ic+1))
         voffx(i) = (fxcen-cxcen)/(cvcx(ic+1)-cvcx(ic))
      end do

      if(ncbx.gt.ncby)then

c=============== CASE 1: x direction is long direction ===================

c ... added to prevent underflow for small crse values

        do n = 1, nvar 
          do j = cslopelo(2)-1,cslopehi(2)+1
            do i = cslopelo(1)-1, cslopehi(1)+1 
              crse(i,j,n) = cvmgt(crse(i,j,n),zero,abs(crse(i,j,n)).gt.1.0e-20)
            end do
          end do
        end do

c ... computed unlimited and limited slopes

        do n = 1, nvar 

c ... --> in x direction

          do j=cslopelo(2), cslopehi(2)
            do i=cslopelo(1), cslopehi(1)
              uc_xslope(i,j,n) = half*(crse(i+1,j,n)-crse(i-1,j,n))

c ... note: the following 6 lines of code is repeated in two other places.
c           A similar code snippet appears three times in the y slope
c           calculation.  Although it looks wasteful, writing the code
c           this way sped up the routine by ~10% (on DEC-alpha). So leave 
c           it alone unless you can make it faster -- rbp

              cen  = uc_xslope(i,j,n)
              forw = two*(crse(i+1,j,n)-crse(i,j,n))
              back = two*(crse(i,j,n)-crse(i-1,j,n))
              slp  = min(abs(forw),abs(back))
              slp  = cvmgp(slp,zero,forw*back)
              lc_xslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))

             end do
          end do
          if (xok) then
            if (bclo(1,n) .eq. EXT_DIR .or. bclo(1,n).eq.HOEXTRAP) then
              i = cslopelo(1)
              do j=cslopelo(2), cslopehi(2)
                uc_xslope(i,j,n)  = -sixteen/fifteen*crse(i-1,j,n) 
     &                      + half*crse(i,j,n)
     $                      + two3rd*crse(i+1,j,n) - tenth*crse(i+2,j,n)

                cen  = uc_xslope(i,j,n)
                forw = two*(crse(i+1,j,n)-crse(i,j,n))
                back = two*(crse(i,j,n)-crse(i-1,j,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_xslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))

               end do
            end if
            if (bchi(1,n) .eq. EXT_DIR .or. bchi(1,n).eq.HOEXTRAP) then
              i = cslopehi(1)
              do j=cslopelo(2), cslopehi(2)
                uc_xslope(i,j,n) = sixteen/fifteen*crse(i+1,j,n) 
     &                      - half*crse(i,j,n)
     $                      - two3rd*crse(i-1,j,n) + tenth*crse(i-2,j,n)

                cen  = uc_xslope(i,j,n)
                forw = two*(crse(i+1,j,n)-crse(i,j,n))
                back = two*(crse(i,j,n)-crse(i-1,j,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_xslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))

               end do
            end if
          end if

c ... --> in y direction

          do j=cslopelo(2), cslopehi(2)
            do i=cslopelo(1), cslopehi(1)
              uc_yslope(i,j,n) = half*(crse(i,j+1,n)-crse(i,j-1,n))

              cen  = uc_yslope(i,j,n)
              forw = two*(crse(i,j+1,n)-crse(i,j,n))
              back = two*(crse(i,j,n)-crse(i,j-1,n))
              slp  = min(abs(forw),abs(back))
              slp  = cvmgp(slp,zero,forw*back)
              lc_yslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))

             end do
          end do
          if (yok) then
            if (bclo(2,n) .eq. EXT_DIR .or. bclo(2,n).eq.HOEXTRAP) then
              j = cslopelo(2)
              do i=cslopelo(1), cslopehi(1)
                uc_yslope(i,j,n)  = -sixteen/fifteen*crse(i,j-1,n) 
     &                      + half*crse(i,j,n)
     $                      + two3rd*crse(i,j+1,n) - tenth*crse(i,j+2,n)

                cen  = uc_yslope(i,j,n)
                forw = two*(crse(i,j+1,n)-crse(i,j,n))
                back = two*(crse(i,j,n)-crse(i,j-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_yslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))

               end do
            end if
            if (bchi(2,n) .eq. EXT_DIR .or. bchi(2,n).eq.HOEXTRAP) then
              j = cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                uc_yslope(i,j,n) = sixteen/fifteen*crse(i,j+1,n) 
     &                    - half*crse(i,j,n)
     $                    - two3rd*crse(i,j-1,n) + tenth*crse(i,j-2,n)

                cen  = uc_yslope(i,j,n)
                forw = two*(crse(i,j+1,n)-crse(i,j,n))
                back = two*(crse(i,j,n)-crse(i,j-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_yslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))

               end do
            end if
          end if
        end do

        if (lin_limit.eq.1)then

c ... compute linear limited slopes
c     Note that the limited and the unlimited slopes
c     have the same sign, and it is assumed that they do.

c ... --> compute slope factors

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
                denom = cvmgt(denom,one,denom.ne.zero)
                factorn = lc_xslope(i,j,n)/denom
                factorn = cvmgt(one,factorn,denom.eq.zero)
                xslope_factor(i,j) = min(xslope_factor(i,j),factorn)

                denom = uc_yslope(i,j,n)
                denom = cvmgt(denom,one,denom.ne.zero)
                factorn = lc_yslope(i,j,n)/denom
                factorn = cvmgt(one,factorn,denom.eq.zero)
                yslope_factor(i,j) = min(yslope_factor(i,j),factorn)
              end do
            end do
          end do

c ... --> compute linear limited slopes

          do n = 1, nvar 
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                lc_xslope(i,j,n) = xslope_factor(i,j)*uc_xslope(i,j,n)
                lc_yslope(i,j,n) = yslope_factor(i,j)*uc_yslope(i,j,n)
              end do
            end do
          end do
        end if

c ... do the interpolation

        do n = 1, nvar
          do j = fblo(2), fbhi(2)
            jc = IX_PROJ(j,lratioy)
            do i = fblo(1), fbhi(1)
              ic = IX_PROJ(i,lratiox)
              fine(i,j,n) = crse(ic,jc,n) + voffx(i)*lc_xslope(ic,jc,n)
     &                                    + voffy(j)*lc_yslope(ic,jc,n)
            end do
          end do
        end do

      else

c=============== CASE 2: y direction is long direction ===================

c ... added to prevent underflow for small crse values

        do n = 1, nvar 
          do i = cslopelo(1)-1, cslopehi(1)+1 
            do j = cslopelo(2)-1,cslopehi(2)+1
              crse(i,j,n) = cvmgt(crse(i,j,n),zero,abs(crse(i,j,n)).gt.1.0e-20)
            end do
          end do
        end do

c ... computed unlimited and limited slopes

        do n = 1, nvar 

c ... --> in x direction

          do i=cslopelo(1), cslopehi(1)
            do j=cslopelo(2), cslopehi(2)
              uc_xslope(i,j,n) = half*(crse(i+1,j,n)-crse(i-1,j,n))

              cen  = uc_xslope(i,j,n)
              forw = two*(crse(i+1,j,n)-crse(i,j,n))
              back = two*(crse(i,j,n)-crse(i-1,j,n))
              slp  = min(abs(forw),abs(back))
              slp  = cvmgp(slp,zero,forw*back)
              lc_xslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))

             end do
          end do
          if (xok) then
            if (bclo(1,n) .eq. EXT_DIR .or. bclo(1,n).eq.HOEXTRAP) then
              i = cslopelo(1)
              do j=cslopelo(2), cslopehi(2)
                uc_xslope(i,j,n)  = -sixteen/fifteen*crse(i-1,j,n) 
     &                      + half*crse(i,j,n)
     $                      + two3rd*crse(i+1,j,n) - tenth*crse(i+2,j,n)

                cen  = uc_xslope(i,j,n)
                forw = two*(crse(i+1,j,n)-crse(i,j,n))
                back = two*(crse(i,j,n)-crse(i-1,j,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_xslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))

               end do
            end if
            if (bchi(1,n) .eq. EXT_DIR .or. bchi(1,n).eq.HOEXTRAP) then
              i = cslopehi(1)
              do j=cslopelo(2), cslopehi(2)
                uc_xslope(i,j,n) = sixteen/fifteen*crse(i+1,j,n) 
     &                      - half*crse(i,j,n)
     $                      - two3rd*crse(i-1,j,n) + tenth*crse(i-2,j,n)

                cen  = uc_xslope(i,j,n)
                forw = two*(crse(i+1,j,n)-crse(i,j,n))
                back = two*(crse(i,j,n)-crse(i-1,j,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_xslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))

               end do
            end if
          end if

c ... --> in y direction

          do i=cslopelo(1), cslopehi(1)
            do j=cslopelo(2), cslopehi(2)
              uc_yslope(i,j,n) = half*(crse(i,j+1,n)-crse(i,j-1,n))

              cen  = uc_yslope(i,j,n)
              forw = two*(crse(i,j+1,n)-crse(i,j,n))
              back = two*(crse(i,j,n)-crse(i,j-1,n))
              slp  = min(abs(forw),abs(back))
              slp  = cvmgp(slp,zero,forw*back)
              lc_yslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))

            end do
          end do
          if (yok) then
            if (bclo(2,n) .eq. EXT_DIR .or. bclo(2,n).eq.HOEXTRAP) then
              j = cslopelo(2)
              do i=cslopelo(1), cslopehi(1)
                uc_yslope(i,j,n)  = -sixteen/fifteen*crse(i,j-1,n) 
     &                      + half*crse(i,j,n)
     $                      + two3rd*crse(i,j+1,n) - tenth*crse(i,j+2,n)

                cen  = uc_yslope(i,j,n)
                forw = two*(crse(i,j+1,n)-crse(i,j,n))
                back = two*(crse(i,j,n)-crse(i,j-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_yslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))

              end do
            end if
            if (bchi(2,n) .eq. EXT_DIR .or. bchi(2,n).eq.HOEXTRAP) then
              j = cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                uc_yslope(i,j,n) = sixteen/fifteen*crse(i,j+1,n) 
     &                    - half*crse(i,j,n)
     $                    - two3rd*crse(i,j-1,n) + tenth*crse(i,j-2,n)

                cen  = uc_yslope(i,j,n)
                forw = two*(crse(i,j+1,n)-crse(i,j,n))
                back = two*(crse(i,j,n)-crse(i,j-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_yslope(i,j,n)=sign(one,cen)*min(slp,abs(cen))

              end do
            end if
          end if
        end do

        if (lin_limit.eq.1)then

c ... compute linear limited slopes
c     Note that the limited and the unlimited slopes
c     have the same sign, and it is assumed that they do.

c ... --> compute slope factors

          do i=cslopelo(1), cslopehi(1)
            do j=cslopelo(2), cslopehi(2)
              xslope_factor(i,j) = 1.0
              yslope_factor(i,j) = 1.0
            end do
          end do

          do n = 1, nvar 
            do i=cslopelo(1), cslopehi(1)
              do j=cslopelo(2), cslopehi(2)
                denom = uc_xslope(i,j,n)
                denom = cvmgt(denom,one,denom.ne.zero)
                factorn = lc_xslope(i,j,n)/denom
                factorn = cvmgt(one,factorn,denom.eq.zero)
                xslope_factor(i,j) = min(xslope_factor(i,j),factorn)

                denom = uc_yslope(i,j,n)
                denom = cvmgt(denom,one,denom.ne.zero)
                factorn = lc_yslope(i,j,n)/denom
                factorn = cvmgt(one,factorn,denom.eq.zero)
                yslope_factor(i,j) = min(yslope_factor(i,j),factorn)
              end do
            end do
          end do

c ... --> compute linear limited slopes

          do n = 1, nvar 
            do i=cslopelo(1), cslopehi(1)
              do j=cslopelo(2), cslopehi(2)
                lc_xslope(i,j,n) = xslope_factor(i,j)*uc_xslope(i,j,n)
                lc_yslope(i,j,n) = yslope_factor(i,j)*uc_yslope(i,j,n)
              end do
            end do
          end do
        end if

c ... do the interpolation

        do n = 1, nvar
          do i = fblo(1), fbhi(1)
            ic = IX_PROJ(i,lratiox)
            do j = fblo(2), fbhi(2)
              jc = IX_PROJ(j,lratioy)
              fine(i,j,n) = crse(ic,jc,n) + voffx(i)*lc_xslope(ic,jc,n)
     &                                    + voffy(j)*lc_yslope(ic,jc,n)
            end do
          end do
        end do

      end if
      end
c ::: 
c ::: --------------------------------------------------------------
c ::: 

      subroutine FORT_CQINTERP (fine, DIMS(fine), 
     $                          fb_l1, fb_l2, fb_h1, fb_h2,
     $                          nvar, lratiox, lratioy, crse, clo, chi, 
     $                          cb_l1, cb_l2, cb_h1, cb_h2,
     $		                fslo, fshi, cslope, clen, fslope, fdat,
     $                          flen, voff, bc, limslope,
     $                          fvcx, fvcy, cvcx, cvcy)

      integer DIMDEC(fine)
      integer fslo(2), fshi(2)
      integer fb_l1, fb_l2, fb_h1, fb_h2
      integer cb_l1, cb_l2, cb_h1, cb_h2
      integer clo, chi
      integer lratiox, lratioy, nvar, clen, flen, limslope
      integer bc(2,2,nvar)
      REAL_T fine(DIMV(fine),nvar)
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

c ::: local var
      integer n, fn
      integer i, ic, ioff
      integer j, jc, joff
      integer ist, jst
      REAL_T cen, forw, back, slp, sgn
      REAL_T fcen, ccen
      REAL_T diffxy,diffxx,diffyy
      REAL_T xoff, yoff
      integer ncbx, ncby
      integer ncsx, ncsy
      integer islo, jslo
      integer icc, istart, iend
      logical xok, yok

      ncbx = cb_h1-cb_l1+1
      ncby = cb_h2-cb_l2+1
      xok = (ncbx .ge. 2)
      yok = (ncby .ge. 2)
      ncsx = ncbx+2
      ncsy = ncby+2
      ist = 1
      jst = ncsx
      islo = cb_l1-1
      jslo = cb_l2-1

      do i = fb_l1, fb_h1 
         fn = i-fslo(1)+1
         ic = IX_PROJ(i,lratiox)
         fcen = half*(fvcx(i)+fvcx(i+1))
         ccen = half*(cvcx(ic)+cvcx(ic+1))
         voff(fn) = (fcen-ccen)/(cvcx(ic+1)-cvcx(ic))
      end do   

c ::: ::::: added to prevent underflow for small crse values
      do n = 1, nvar 
        do i = clo, chi 
          crse(i,n) = cvmgt(crse(i,n),zero,abs(crse(i,n)).gt.1.0e-20)
        end do
      end do

      do 290 n = 1, nvar 

c ::: ::::: compute slopes in x direction
            do i = 1, clen 
               cen = half*(crse(i+ist,n)-crse(i-ist,n))
               diffxy = fourth*(crse(i+ist+jst,n)+crse(i-ist-jst,n)
     $                         -crse(i-ist+jst,n)-crse(i+ist-jst,n))
               diffxx = crse(i+ist,n)-two*crse(i,n)+crse(i-ist,n)
               cslope(i,1)=cen
               cslope(i,3)=diffxx
               cslope(i,5)=diffxy
            end do
            if (xok) then
               if (bclo(1,n) .eq. EXT_DIR .or. bclo(1,n).eq.HOEXTRAP) then
                  do i = 1, clen, jst 
                     cen  = -sixteen/fifteen*crse(i-ist,n) + half*crse(i,n)
     $                    + two3rd*crse(i+ist,n) - tenth*crse(i+2*ist,n)
                     cslope(i,1)=cen
                     cslope(i,3)=zero
                     cslope(i,5)=zero
                  end do
               end if
               if (bchi(1,n) .eq. EXT_DIR .or. bchi(1,n).eq.HOEXTRAP) then
                  do i = ncbx, clen, jst 
                     cen = sixteen/fifteen*crse(i+ist,n) - half*crse(i,n)
     $                    - two3rd*crse(i-ist,n) + tenth*crse(i-2*ist,n)
                     cslope(i,1)=cen
                     cslope(i,3)=zero
                     cslope(i,5)=zero
                  end do
               end if
            end if

c ::: ::::: compute slopes in y direction
            do i = 1, clen 
               cen  = half*(crse(i+jst,n)-crse(i-jst,n))
               diffyy = crse(i+jst,n)-two*crse(i,n)+crse(i-jst,n)
               cslope(i,2)=cen
               cslope(i,4)=diffyy
            end do
            if (yok) then
               if (bclo(2,n) .eq. EXT_DIR .or. bclo(2,n).eq.HOEXTRAP) then
                  do i = 1, ncbx 
                     cen  = -sixteen/fifteen*crse(i-jst,n) + half*crse(i,n)
     $                    + two3rd*crse(i+jst,n) - tenth*crse(i+2*jst,n)
                     cslope(i,2)=cen
                     cslope(i,4)=zero
                     cslope(i,5)=zero
                  end do
               end if
               if (bchi(2,n) .eq. EXT_DIR .or. bchi(2,n).eq.HOEXTRAP) then
                  do i = clen-ncbx,clen 
                     cen = sixteen/fifteen*crse(i+jst,n) - half*crse(i,n)
     $                    - two3rd*crse(i-jst,n) + tenth*crse(i-2*jst,n)
                     cslope(i,2)=cen
                     cslope(i,4)=zero
                     cslope(i,5)=zero
                  end do
               end if
            end if

            do 360 jc = cb_l2, cb_h2 
c           ::::: strip out a fine grid slope vector
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
                     fine(i,j,n) = fdat(fn) + voff(fn)*fslope(fn,1)
     &                    + voff(fn)*voff(fn)*fslope(fn,3)
     &                    + yoff*fslope(fn,2)
     &                    + yoff*yoff*fslope(fn,4)
     &                    + voff(fn)*yoff*fslope(fn,5)
400               continue
390            continue
360         continue

290   continue

      end

c ::: 
c ::: --------------------------------------------------------------
c ::: pcinterp:  cell centered piecewise constant interpolation
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: flo,fhi      =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: 
c ::: crse         =>  (const)  coarse grid data 
c ::: clo,chi      =>  (const)  index limits of coarse grid
c ::: cblo,cbhi    =>  (const) coarse grid region containing fblo,fbhi
c ::: 
c ::: longdir      =>  (const)  which index direction is longest (1 or 2)
c ::: lratio(2)    =>  (const)  refinement ratio between levels
c ::: nvar         =>  (const)  number of components in array
c ::: 
c ::: TEMPORARY ARRAYS
c ::: ftmp         =>  1-D temp array
c ::: --------------------------------------------------------------
c ::: 
      subroutine FORT_PCINTERP (crse,DIMS(crse),cblo,cbhi,
     &                          fine,DIMS(fine),fblo,fbhi,
     &                          longdir,lratiox,lratioy,nvar,
     &                          ftmp,ftmp_lo,ftmp_hi)
      integer DIMDEC(crse)
      integer cblo(2), cbhi(2)
      integer DIMDEC(fine)
      integer fblo(2), fbhi(2)
      integer ftmp_lo, ftmp_hi
      integer nvar, lratiox, lratioy, longdir
      REAL_T  crse(DIMV(crse), nvar)
      REAL_T  fine(DIMV(fine), nvar)
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

      end





