c
c $Id: INTERP_3D.F,v 1.11 2002-08-29 22:14:34 car Exp $
c

#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "INTERP_F.H"
#include "ArrayLim.H"

#define IX_PROJ(A,B) (A+B*iabs(A))/B-iabs(A)
#define SDIM 3

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
c ::: strip        =>  1-D temp array
c ::: --------------------------------------------------------------
c :::
      subroutine FORT_NBINTERP (crse, DIMS(crse), DIMS(cb),
     $                          fine, DIMS(fine), DIMS(fb),
     $                          lratiox, lratioy, lratioz, nvar,
     $                          sl, num_slp)

      integer DIMDEC(crse)
      integer DIMDEC(cb)
      integer DIMDEC(fine)
      integer DIMDEC(fb)
      integer lratiox, lratioy, lratioz, nvar
      integer num_slp
      REAL_T  fine(DIMV(fine),nvar)
      REAL_T  crse(DIMV(crse),nvar)
      REAL_T  sl(DIM1(cb),num_slp)

#define  SLX 1
#define  SLY 2
#define  SLZ 3
#define  SLXY 4
#define  SLXZ 5
#define  SLYZ 6
#define  SLXYZ 7

c ::: local var
      integer lx, ly, lz
      integer i, j, k, ifn, jfn, kfn, n
      integer ilo, ihi, jlo, jhi, klo, khi
      integer kstrtFine, kstopFine, jstrtFine, jstopFine, istrtFine, istopFine

      REAL_T fx, fy,fz
      REAL_T RX, RY, RZ, RXY, RXZ, RYZ, RXYZ
      REAL_T dx00, d0x0, d00x, dx10, dx01, d0x1, dx11
      REAL_T slope

      slope(i,j,k,n,fx,fy,fz) = crse(i,j,k,n) +
     &     fx*sl(i,SLX) + fy*sl(i,SLY) + fz*sl(i,SLZ) +
     &     fx*fy*sl(i,SLXY) + fx*fz*sl(i,SLXZ) + fy*fz*sl(i,SLYZ) +
     &     fx*fy*fz*sl(i,SLXYZ)

      RX = one/dfloat(lratiox)
      RY = one/dfloat(lratioy)
      RZ = one/dfloat(lratioz)
      RXY = RX*RY
      RXZ = RX*RZ
      RYZ = RY*RZ
      RXYZ = RX*RY*RZ

c
c     NOTES:
c         1) (i, j, k) loop over the coarse cells
c         2) ?strtFine and ?stopFine are the beginning and ending fine cell
c            indices corresponding to the current coarse cell.  ?stopFine
c            is restricted for the last coarse cell in each direction since
c            for this cell we only need to do the face and not the fine nodes
c            inside this cell.
c         3) (lx, ly, lz) as well as ?lo and ?hi refer to the fine node indices
c            as an offset from ?strtFine.
c
      do 100 n = 1, nvar
        do 110 k = ARG_L3(cb), ARG_H3(cb)
          kstrtFine = k * lratioz
          kstopFine = kstrtFine + lratioz - 1
          if (k .eq. ARG_H3(cb)) kstopFine = kstrtFine

          klo = max(ARG_L3(fb),kstrtFine) - kstrtFine
          khi = min(ARG_H3(fb),kstopFine) - kstrtFine

          do 120 j = ARG_L2(cb), ARG_H2(cb)
            jstrtFine = j * lratioy
            jstopFine = jstrtFine + lratioy - 1
            if (j .eq. ARG_H2(cb)) jstopFine = jstrtFine

            jlo = max(ARG_L2(fb),jstrtFine) - jstrtFine
            jhi = min(ARG_H2(fb),jstopFine) - jstrtFine

c
c           ::::: compute slopes :::::
c
c           NOTE: The IF logic in the calculation of the slopes is to
c                 prevent stepping out of bounds on the coarse data when
c                 computing the slopes on the ARG_H?(cb) cells.  These
c                 slopes actually are not used since they are multiplied by
c                 zero.
c
            do i = ARG_L1(cb), ARG_H1(cb)
              dx00 = zero
              if (i .NE. ARG_H1(cb)) dx00 = crse(i+1,j,k,n) - crse(i,j,k,n)

              d0x0 = zero
              if (j .NE. ARG_H2(cb)) d0x0 = crse(i,j+1,k,n) - crse(i,j,k,n)

              d00x = zero
              if (k .NE. ARG_H3(cb)) d00x = crse(i,j,k+1,n) - crse(i,j,k,n)

              dx10 = zero
              if (i .NE. ARG_H1(cb) .and. j .NE. ARG_H2(cb))
     $          dx10 = crse(i+1,j+1,k,n) - crse(i,j+1,k,n)

              dx01 = zero
              if (i .NE. ARG_H1(cb) .and. k .NE. ARG_H3(cb))
     $          dx01 = crse(i+1,j,k+1,n) - crse(i,j,k+1,n)

              d0x1 = zero
              if (j .NE. ARG_H2(cb) .and. k .NE. ARG_H3(cb))
     $          d0x1 = crse(i,j+1,k+1,n) - crse(i,j,k+1,n)

              dx11 = zero
              if (i .NE. ARG_H1(cb) .and. j .NE. ARG_H2(cb)
     $                              .and. k .NE. ARG_H3(cb)) 
     $          dx11 = crse(i+1,j+1,k+1,n) - crse(i,j+1,k+1,n)

              sl(i,SLX) = RX*dx00
              sl(i,SLY) = RY*d0x0
              sl(i,SLZ) = RZ*d00x

              sl(i,SLXY) = RXY*(dx10 - dx00)
              sl(i,SLXZ) = RXZ*(dx01 - dx00)
              sl(i,SLYZ) = RYZ*(d0x1 - d0x0)

              sl(i,SLXYZ) = RXYZ*(dx11 - dx01 - dx10 + dx00)
            end do

c
c           ::::: compute fine strip of interpolated data
c
            do lz = klo, khi
              kfn = lratioz * k + lz
              fz = dfloat(lz)

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

                    fine(ifn,jfn,kfn,n) = slope(i,j,k,n,fx,fy,fz)
                  end do
                end do
              end do
            end do

c         ::::: end of do j
120       continue
c       ::::: end of do k
110     continue
c     ::::: end of do n
100   continue

#undef  SLX
#undef  SLY
#undef  SLZ
#undef  SLXY
#undef  SLXZ
#undef  SLYZ
#undef  SLXYZ

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
c ::: lratio(3)    =>  (const)  refinement ratio between levels
c ::: nvar         =>  (const)  number of components in array
c ::: 
c ::: TEMPORARY ARRAYS
c ::: slx,sly,slxy =>  1-D slope arrays
c ::: strip        =>  1-D temp array
c ::: --------------------------------------------------------------
c ::: 
      subroutine FORT_CBINTERP (crse, DIMS(crse), DIMS(cb),
     $                          fine, DIMS(fine), DIMS(fb),
     $                          lratiox, lratioy, lratioz, nvar,
     $                          sl, num_slp, strip, strip_lo, strip_hi)

      integer DIMDEC(crse)
      integer DIMDEC(cb)
      integer DIMDEC(fine)
      integer DIMDEC(fb)
      integer lratiox, lratioy, lratioz, nvar
      integer num_slp
      integer strip_lo, strip_hi
      REAL_T fine(DIMV(fine),nvar)
      REAL_T crse(DIMV(crse),nvar)
      REAL_T strip(strip_lo:strip_hi)
      REAL_T sl(DIM1(cb), num_slp)

c ::: local var
#if 0
      integer lx, ly
      integer hrat, ic, jc, jfn, jfc, i, j, n
      REAL_T x, y
      REAL_T denom
#endif
      write(6,*) "FORT_CBINTERP not implemented"
      stop

      end

c ::: --------------------------------------------------------------
c ::: ccinterp:   conservative interpolation from coarse grid to
c ::: subregion of fine grid defined by (fblo,fbhi)
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: DIMS(fine)   =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: nvar         =>  (const)  number of variables in state vector
c ::: lratio(3)    =>  (const)  refinement ratio between levels
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
c ::: 2-D EXAMPLE:
c ::: Suppose the patch called "fine" has index extent:
c ::: 
c ::: floi1 = 3, fhii1 = 12
c ::: floi2 = 8, fhii2 = 20
c ::: 
c ::: suppose the subergion of this patch that is to be filled 
c ::: by interpolation has index extent:
c ::: 
c ::: fb_l1 = 5,  fb_h1 = 10
c ::: fb_l2 = 13, fb_h2 = 20
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
c ::: The point (cb_l1,cb_l2) -> 1
c ::: The point (cslo(1),cslo(2)) -> clo = 1 - 1 - ncsx = -6
c ::: 
c ::: The point (cbhi(1),cbhi(2)) -> clen = ncby*ncsx - 2 = 5*6-2 = 28
c ::: The point (cshi(1),cshi(2)) -> chi = clo + ncsx*ncsy - 1 
c :::                                    =  -6 +    6*7    - 1 = 35
c ::: 
c :::      -------------------------------------------------
c :::      |       |       |       |       |       |  chi  |  
c :::  11  |   30  |   31  |   32  |   33  |   34  |   35  |   cshi(2)
c :::      |       |       |       |       |       |       |
c :::      -------------------------------------------------
c :::      |       |       |       |       |  clen |       |  
c :::  10  |   24  |   25  |   26  |   27  |   28  |   29  |   cb_h2
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
c :::   6  |    0  |    1  |    2  |    3  |    4  |    5  |   cb_l2
c :::      |       |       |       |       |       |       |
c :::      -------------------------------------------------
c :::      |  clo  |       |       |       |       |       |  
c :::   5  |   -6  |   -5  |   -4  |   -3  |   -2  |   -1  |   cslo(2)
c :::      |       |       |       |       |       |       |
c :::      -------------------------------------------------
c :::          1       2       3       4       5       6
c :::               cb_l1                    cb_h1
c :::       cslo(1)                                 cshi(1)
c ::: 
c ::: 
c ::: In the 1-D coordinates:
c :::    ist = 1    = stride in I direction
c :::    jst = 6    = stride in J direction  (ncsx)
c ::: 
c ::: --------------------------------------------------------------
      subroutine FORT_CCINTERP (fine, DIMS(fine), 
     $                          DIMS(fb),
     $                          nvar,lratiox,lratioy,lratioz,crse,clo,
     $                          chi, DIMS(cb),
     $		                fslo, fshi, cslope, clen, fslope, fdat,
     $                          flen, voff, bc, limslope,
     $                          fvcx, fvcy, fvcz, cvcx, cvcy, cvcz)

      integer DIMDEC(fine)
      integer DIMDEC(fb)
      integer DIMDEC(cb)
      integer fslo(3), fshi(3)
      integer nvar, lratiox, lratioy, lratioz
      integer bc(3,2,nvar)
      integer clen, flen, clo, chi, limslope
      REAL_T fine(DIMV(fine),nvar)
      REAL_T crse(clo:chi, nvar)
      REAL_T cslope(clo:chi, 3)
      REAL_T fslope(flen, 3)
      REAL_T fdat(flen)
      REAL_T voff(flen)
      REAL_T fvcx(fb_l1:fb_h1+1)
      REAL_T fvcy(fb_l2:fb_h2+1)
      REAL_T fvcz(fb_l3:fb_h3+1)
      REAL_T cvcx(cb_l1:cb_h1+1)
      REAL_T cvcy(cb_l2:cb_h2+1)
      REAL_T cvcz(cb_l3:cb_h3+1)

#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)

c ::: local var
      integer n, fn
      integer i, ii, ic, ioff
      integer j, jj, jc, joff
      integer k, kk, kc, koff
      integer ist, jst, kst
      integer cslo(3), cshi(3)
      REAL_T cen, forw, back, slp, sgn
      REAL_T fcen, ccen
      REAL_T xoff, yoff, zoff
      integer ncbx, ncby, ncbz
      integer ncsx, ncsy, ncsz
      integer islo, jslo, kslo
      integer icc, istart, iend
      integer ilo, ihi, jlo, jhi, klo, khi
      logical xok, yok, zok

c     :::::: helpful statement functions
      integer sloc
      integer strd
      REAL_T  slplft, slprgt
      sloc(i,j,k) = clo+i-cslo(1)+ncsx*(j-cslo(2)+ncsy*(k-cslo(3)))

      slplft(i,strd,n) =  -sixteen/fifteen*crse(i-strd,n)
     $		         + half*crse(i,n)
     $                   + two3rd*crse(i+strd,n)
     $			 - tenth*crse(i+2*strd,n)
      slprgt(i,strd,n) =   sixteen/fifteen*crse(i+strd,n)
     $		         - half*crse(i,n)
     $                   - two3rd*crse(i-strd,n)
     $			 + tenth*crse(i-2*strd,n)

      ncbx = cb_h1-cb_l1+1
      ncby = cb_h2-cb_l2+1
      ncbz = cb_h3-cb_l3+1
      cslo(1) = cb_l1-1
      cshi(1) = cb_h1+1
      cslo(2) = cb_l2-1
      cshi(2) = cb_h2+1
      cslo(3) = cb_l3-1
      cshi(3) = cb_h3+1
      xok = (cb_h1-cb_l1+1 .ge. 2)
      yok = (cb_h2-cb_l2+1 .ge. 2)
      zok = (cb_h3-cb_l3+1 .ge. 2)
      ncsx = ncbx+2
      ncsy = ncby+2
      ncsz = ncbz+2
      ist = 1
      jst = ncsx
      kst = ncsx*ncsy
      islo = cb_l1-1
      jslo = cb_l2-1
      kslo = cb_l3-1

      do i = fb_l1, fb_h1 
         fn = i-fslo(1)+1
         ic = IX_PROJ(i,lratiox)
         fcen = half*(fvcx(i)+fvcx(i+1))
         ccen = half*(cvcx(ic)+cvcx(ic+1))
         voff(fn) = (fcen-ccen)/(cvcx(ic+1)-cvcx(ic))
      end do

      do 100 n = 1, nvar 

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
                      cen  = slplft(i,ist,n)
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
                      cen = slprgt(i,ist,n)
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
                cslope(i,1) = half*(crse(i+ist,n)-crse(i-ist,n))
             end do
             if (xok) then
                if (bclo(1,n) .eq. EXT_DIR .or. bclo(1,n).eq.HOEXTRAP) then
                   do i = 1, clen, jst 
                      cslope(i,1) = slplft(i,ist,n)
                   end do
                end if
                if (bchi(1,n) .eq. EXT_DIR .or. bchi(1,n).eq.HOEXTRAP) then
                   do i = ncbx, clen, jst 
                      cslope(i,1) = slprgt(i,ist,n)
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
                   do k = cb_l3, cb_h3
                      ilo = sloc(cb_l1,cb_l2,k)
                      ihi = sloc(cb_h1,cb_l2,k)
                      do i = ilo, ihi
                         cen  = slplft(i,jst,n)
                         sgn  = sign(one,cen)
                         forw = two*(crse(i+jst,n)-crse(i,n))
                         back = two*(crse(i,n)-crse(i-jst,n))
                         slp  = min(abs(forw),abs(back))
                         slp  = cvmgp(slp,zero,forw*back)
                         cslope(i,2)=sgn*min(slp,abs(cen))
                      end do
                   end do
                end if
                if (bchi(2,n) .eq. EXT_DIR .or. bchi(2,n).eq.HOEXTRAP) then
                   do k = cb_l3, cb_h3
                      ilo = sloc(cb_l1,cb_h2,k)
                      ihi = sloc(cb_h1,cb_h2,k)
                      do i = ilo, ihi
                         cen  = slprgt(i,jst,n)
                         sgn  = sign(one,cen)
                         forw = two*(crse(i+jst,n)-crse(i,n))
                         back = two*(crse(i,n)-crse(i-jst,n))
                         slp  = min(abs(forw),abs(back))
                         slp  = cvmgp(slp,zero,forw*back)
                         cslope(i,2)=sgn*min(slp,abs(cen))
                      end do
                   end do
                end if
             end if
	  else
            do i = 1, clen 
              cslope(i,2) = half*(crse(i+jst,n)-crse(i-jst,n))
           end do
           if (yok) then
              if (bclo(2,n) .eq. EXT_DIR .or. bclo(2,n).eq.HOEXTRAP) then
                 do k = cb_l2, cb_h3
                    ilo = sloc(cb_l1,cb_l2,k)
                    ihi = sloc(cb_h1,cb_l2,k)
                    do i = ilo, ihi
                       cslope(i,2) = slplft(i,jst,n)
                    end do
                 end do
              end if
              if (bchi(2,n) .eq. EXT_DIR .or. bchi(2,n).eq.HOEXTRAP) then
                 do k = cb_l3, cb_h3
                    ilo = sloc(cb_l1,cb_h2,k)
                    ihi = sloc(cb_h1,cb_h2,k)
                    do i = ilo, ihi
                       cslope(i,2) = slprgt(i,jst,n)
                    end do
                 end do
              end if
           end if
         end if
c ::: ::::: compute slopes in z direction
          if (limslope .ne. 0) then
             do i = 1, clen 
                cen  = half*(crse(i+kst,n)-crse(i-kst,n))
                forw = two*(crse(i+kst,n)-crse(i,n))
                back = two*(crse(i,n)-crse(i-kst,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                cslope(i,3)=sign(one,cen)*min(slp,abs(cen))
             end do
             if (zok) then
                if (bclo(3,n) .eq. EXT_DIR .or. bclo(3,n).eq.HOEXTRAP) then
                   ilo = sloc(cb_l1,cb_l2,cb_l3)
                   ihi = sloc(cb_h1,cb_h2,cb_l3)
                   do i = ilo, ihi
                      cen  = slplft(i,kst,n)
                      sgn  = sign(one,cen)
                      forw = two*(crse(i+kst,n)-crse(i,n))
                      back = two*(crse(i,n)-crse(i-kst,n))
                      slp  = min(abs(forw),abs(back))
                      slp  = cvmgp(slp,zero,forw*back)
                      cslope(i,3)=sgn*min(slp,abs(cen))
                   end do
                end if
                if (bchi(3,n) .eq. EXT_DIR .or. bchi(3,n).eq.HOEXTRAP) then
                   ilo = sloc(cb_l1,cb_l2,cb_h3)
                   ihi = sloc(cb_h1,cb_h2,cb_h3)
                   do i = ilo, ihi
                      cen  = slprgt(i,kst,n)
                      sgn  = sign(one,cen)
                      forw = two*(crse(i+kst,n)-crse(i,n))
                      back = two*(crse(i,n)-crse(i-kst,n))
                      slp  = min(abs(forw),abs(back))
                      slp  = cvmgp(slp,zero,forw*back)
                      cslope(i,3)=sgn*min(slp,abs(cen))
                   end do
                end if
             end if
	  else
             do i = 1, clen 
                cslope(i,3) = half*(crse(i+kst,n)-crse(i-kst,n))
             end do
             if (zok) then
                if (bclo(3,n) .eq. EXT_DIR .or. bclo(3,n).eq.HOEXTRAP) then
                   ilo = sloc(cb_l1,cb_l2,cb_l3)
                   ihi = sloc(cb_h1,cb_h2,cb_l3)
                   do i = ilo, ihi
                      cslope(i,3) = slplft(i,kst,n)
                   end do
                end if
                if (bchi(3,n) .eq. EXT_DIR .or. bchi(3,n).eq.HOEXTRAP) then
                   ilo = sloc(cb_l1,cb_l2,cb_h3)
                   ihi = sloc(cb_h1,cb_h2,cb_h3)
                   do i = ilo, ihi
                      cslope(i,3) = slprgt(i,kst,n)
                   end do
                end if
             end if
	  end if

          do 200 kc = cb_l3, cb_h3 
             do 250 jc = cb_l2, cb_h2 

c :::           ::::: strip out fine grid slope vectors
                do ioff = 1, lratiox
                   icc = sloc(cb_l1,jc,kc)
                   istart = ioff
                   iend = ioff + (ncbx-1)*lratiox
                   do fn = istart, iend, lratiox 
                      fslope(fn,1) = cslope(icc,1)
                      fslope(fn,2) = cslope(icc,2)
                      fslope(fn,3) = cslope(icc,3)
                      fdat(fn) = crse(icc,n)
                      icc = icc + ist
                   end do
                end do

c :::           ::::: interpolate to fine grid
                jj = lratioy*jc
                jlo = max(jj,fb_l2) - jj
                jhi = min(jj+lratioy-1,fb_h2) - jj
                kk = lratioz*kc
                klo = max(kk,fb_l3) - kk
                khi = min(kk+lratioz-1,fb_h3) - kk
                do koff = klo, khi
                   k = lratioz*kc + koff
                   fcen = half*(fvcz(k) +fvcz(k+1))
                   ccen = half*(cvcz(kc)+cvcz(kc+1))
                   zoff = (fcen-ccen)/(cvcz(kc+1)-cvcz(kc))
                   do joff = jlo, jhi
                      j = lratioy*jc + joff
                      fcen = half*(fvcy(j) +fvcy(j+1))
                      ccen = half*(cvcy(jc)+cvcy(jc+1))
                      yoff = (fcen-ccen)/(cvcy(jc+1)-cvcy(jc))
                      do i = fb_l1, fb_h1
                         fn = i-fslo(1)+1
                         fine(i,j,k,n) = fdat(fn) + 
     $                        voff(fn)*fslope(fn,1)+
     $                        yoff*fslope(fn,2)+
     $                        zoff*fslope(fn,3)
                      end do
                   end do
                end do
250          continue
200       continue

100   continue

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
     &                          uc_zslope, lc_zslope, zslope_factor,
     &                          DIMS(cslope),
     &                          cslopelo, cslopehi,
     $                          nvar, lratiox, lratioy, lratioz, 
     $                          bc, lin_limit,
     $                          fvcx, fvcy, fvcz, cvcx, cvcy, cvcz,
     $                          voffx,voffy,voffz)

      implicit none

      integer DIMDEC(fine)
      integer DIMDEC(crse)
      integer DIMDEC(fvcb)
      integer DIMDEC(cvcb)
      integer DIMDEC(cslope)
      integer fblo(3), fbhi(3)
      integer cslopelo(3), cslopehi(3)
      integer lratiox, lratioy, lratioz, nvar, lin_limit
      integer bc(3,2,nvar)
      REAL_T fine(DIMV(fine),nvar)
      REAL_T crse(DIMV(crse), nvar)
      REAL_T uc_xslope(DIMV(cslope),nvar)
      REAL_T lc_xslope(DIMV(cslope),nvar)
      REAL_T xslope_factor(DIMV(cslope))
      REAL_T uc_yslope(DIMV(cslope),nvar)
      REAL_T lc_yslope(DIMV(cslope),nvar)
      REAL_T yslope_factor(DIMV(cslope))
      REAL_T uc_zslope(DIMV(cslope),nvar)
      REAL_T lc_zslope(DIMV(cslope),nvar)
      REAL_T zslope_factor(DIMV(cslope))
      REAL_T fvcx(DIM1(fvcb))
      REAL_T fvcy(DIM2(fvcb))
      REAL_T fvcz(DIM3(fvcb))
      REAL_T voffx(DIM1(fvcb))
      REAL_T voffy(DIM2(fvcb))
      REAL_T voffz(DIM3(fvcb))       
      REAL_T cvcx(DIM1(cvcb))
      REAL_T cvcy(DIM2(cvcb))
      REAL_T cvcz(DIM3(cvcb))

#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)

c ::: local var
      integer n 
      integer i, ic
      integer j, jc
      integer k, kc
      REAL_T cen, forw, back, slp
      REAL_T factorn, denom
      REAL_T fxcen, cxcen, fycen, cycen, fzcen, czcen
      logical xok, yok, zok
      integer ncbx, ncby, ncbz

      ncbx = cslopehi(1)-cslopelo(1)+1
      ncby = cslopehi(2)-cslopelo(2)+1
      ncbz = cslopehi(3)-cslopelo(3)+1
      xok = (ncbx .ge. 2)
      yok = (ncby .ge. 2)
      zok = (ncbz .ge. 2)

      do k = fblo(3), fbhi(3)
        kc = IX_PROJ(k,lratioz)
        fzcen = half*(fvcz(k)+fvcz(k+1))
        czcen = half*(cvcz(kc)+cvcz(kc+1))
        voffz(k) = (fzcen-czcen)/(cvcz(kc+1)-cvcz(kc))
      end do
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

      if(ncbx.ge.ncby.and.ncbx.ge.ncbz)then

c=============== CASE 1: x direction is long direction ======================

c ... added to prevent underflow for small crse values

        do n = 1, nvar 
          do k = cslopelo(3)-1,cslopehi(3)+1
            do j = cslopelo(2)-1,cslopehi(2)+1
              do i = cslopelo(1)-1, cslopehi(1)+1 
                crse(i,j,k,n) = cvmgt(crse(i,j,k,n),zero,abs(crse(i,j,k,n)).gt.1.0D-20)
              end do
            end do
          end do
        end do

c ... computed unlimited and limited slopes

        do n = 1, nvar 

c ... --> in x direction

          do k=cslopelo(3), cslopehi(3)
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                uc_xslope(i,j,k,n) = half*(crse(i+1,j,k,n)-crse(i-1,j,k,n))


c ... note: the following 6 lines of code is repeated in two other places.
c           Similar code snippets appear three times in both the y and z
c           slope calculations.  Although it looks wasteful, writing the code
c           this way sped up the routine by ~10% (on DEC-alpha). So leave 
c           it alone unless you can make it faster -- rbp

                cen  = uc_xslope(i,j,k,n)
                forw = two*(crse(i+1,j,k,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i-1,j,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_xslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

             end do
            end do
          end do
          if (xok) then
            if (bclo(1,n) .eq. EXT_DIR .or. bclo(1,n).eq.HOEXTRAP) then
              i = cslopelo(1)
              do k=cslopelo(3), cslopehi(3)
                do j=cslopelo(2), cslopehi(2)
                  uc_xslope(i,j,k,n)  = -sixteen/fifteen*crse(i-1,j,k,n) 
     &                        + half*crse(i,j,k,n)
     $                        + two3rd*crse(i+1,j,k,n) - tenth*crse(i+2,j,k,n)

                  cen  = uc_xslope(i,j,k,n)
                  forw = two*(crse(i+1,j,k,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i-1,j,k,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_xslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                end do
              end do
            end if
            if (bchi(1,n) .eq. EXT_DIR .or. bchi(1,n).eq.HOEXTRAP) then
              i = cslopehi(1)
              do k=cslopelo(3), cslopehi(3)
                do j=cslopelo(2), cslopehi(2)
                  uc_xslope(i,j,k,n) = sixteen/fifteen*crse(i+1,j,k,n) 
     &                      - half*crse(i,j,k,n)
     $                      - two3rd*crse(i-1,j,k,n) + tenth*crse(i-2,j,k,n)

                  cen  = uc_xslope(i,j,k,n)
                  forw = two*(crse(i+1,j,k,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i-1,j,k,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_xslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                end do
              end do
            end if
          end if

c ... --> in y direction

          do k=cslopelo(3), cslopehi(3)
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                uc_yslope(i,j,k,n) = half*(crse(i,j+1,k,n)-crse(i,j-1,k,n))

                cen  = uc_yslope(i,j,k,n)
                forw = two*(crse(i,j+1,k,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i,j-1,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_yslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

               end do
            end do
          end do
          if (yok) then
            if (bclo(2,n) .eq. EXT_DIR .or. bclo(2,n).eq.HOEXTRAP) then
              j = cslopelo(2)
              do k=cslopelo(3), cslopehi(3)
                do i=cslopelo(1), cslopehi(1)
                  uc_yslope(i,j,k,n)  = -sixteen/fifteen*crse(i,j-1,k,n) 
     &                        + half*crse(i,j,k,n)
     $                        + two3rd*crse(i,j+1,k,n) - tenth*crse(i,j+2,k,n)

                  cen  = uc_yslope(i,j,k,n)
                  forw = two*(crse(i,j+1,k,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i,j-1,k,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_yslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                end do
              end do
            end if
            if (bchi(2,n) .eq. EXT_DIR .or. bchi(2,n).eq.HOEXTRAP) then
              j = cslopehi(2)
              do k=cslopelo(3), cslopehi(3)
                do i=cslopelo(1), cslopehi(1)
                  uc_yslope(i,j,k,n) = sixteen/fifteen*crse(i,j+1,k,n) 
     &                      - half*crse(i,j,k,n)
     $                      - two3rd*crse(i,j-1,k,n) + tenth*crse(i,j-2,k,n)

                  cen  = uc_yslope(i,j,k,n)
                  forw = two*(crse(i,j+1,k,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i,j-1,k,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_yslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                end do
              end do
            end if
          end if

c ... --> in z direction

          do k=cslopelo(3), cslopehi(3)
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                uc_zslope(i,j,k,n) = half*(crse(i,j,k+1,n)-crse(i,j,k-1,n))

                cen  = uc_zslope(i,j,k,n)
                forw = two*(crse(i,j,k+1,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i,j,k-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_zslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

              end do
            end do
          end do
          if (zok) then
            if (bclo(3,n) .eq. EXT_DIR .or. bclo(3,n).eq.HOEXTRAP) then
              k = cslopelo(3)
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  uc_zslope(i,j,k,n)  = -sixteen/fifteen*crse(i,j,k-1,n) 
     &                        + half*crse(i,j,k,n)
     $                        + two3rd*crse(i,j,k+1,n) - tenth*crse(i,j,k+2,n)
 
                  cen  = uc_zslope(i,j,k,n)
                  forw = two*(crse(i,j,k+1,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i,j,k-1,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_zslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                end do
              end do
            end if
            if (bchi(3,n) .eq. EXT_DIR .or. bchi(3,n).eq.HOEXTRAP) then
              k = cslopehi(3)
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  uc_zslope(i,j,k,n) = sixteen/fifteen*crse(i,j,k+1,n) 
     &                      - half*crse(i,j,k,n)
     $                      - two3rd*crse(i,j,k-1,n) + tenth*crse(i,j,k-2,n)
  
                  cen  = uc_zslope(i,j,k,n)
                  forw = two*(crse(i,j,k+1,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i,j,k-1,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_zslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

               end do
              end do
            end if
          end if
        end do

        if (lin_limit.eq.1)then

c ... compute linear limited slopes
c     Note that the limited and the unlimited slopes
c     have the same sign, and it is assumed that they do.

c ... --> compute slope factors

          do k=cslopelo(3), cslopehi(3)
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                xslope_factor(i,j,k) = one
                yslope_factor(i,j,k) = one
                zslope_factor(i,j,k) = one
              end do
            end do
          end do

          do n = 1, nvar 
            do k=cslopelo(3), cslopehi(3)
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  denom = uc_xslope(i,j,k,n)
                  denom = cvmgt(denom,one,denom.ne.zero)
                  factorn = lc_xslope(i,j,k,n)/denom
                  factorn = cvmgt(one,factorn,denom.eq.zero)
                  xslope_factor(i,j,k) = min(xslope_factor(i,j,k),factorn)

                  denom = uc_yslope(i,j,k,n)
                  denom = cvmgt(denom,one,denom.ne.zero)
                  factorn = lc_yslope(i,j,k,n)/denom
                  factorn = cvmgt(one,factorn,denom.eq.zero)
                  yslope_factor(i,j,k) = min(yslope_factor(i,j,k),factorn)

                  denom = uc_zslope(i,j,k,n)
                  denom = cvmgt(denom,one,denom.ne.zero)
                  factorn = lc_zslope(i,j,k,n)/denom
                  factorn = cvmgt(one,factorn,denom.eq.zero)
                  zslope_factor(i,j,k) = min(zslope_factor(i,j,k),factorn)
                end do
              end do
            end do
          end do

c ... --> compute linear limited slopes

          do n = 1, nvar 
            do k=cslopelo(3), cslopehi(3)
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  lc_xslope(i,j,k,n) = xslope_factor(i,j,k)*uc_xslope(i,j,k,n)
                  lc_yslope(i,j,k,n) = yslope_factor(i,j,k)*uc_yslope(i,j,k,n)
                  lc_zslope(i,j,k,n) = zslope_factor(i,j,k)*uc_zslope(i,j,k,n)
                end do
              end do
            end do
          end do
        end if

c ... do the interpolation

        do n = 1, nvar
          do k = fblo(3), fbhi(3)
            kc = IX_PROJ(k,lratioz)
            do j = fblo(2), fbhi(2)
              jc = IX_PROJ(j,lratioy)
              do i = fblo(1), fbhi(1)
                ic = IX_PROJ(i,lratiox)
                fine(i,j,k,n) = crse(ic,jc,kc,n) + voffx(i)*lc_xslope(ic,jc,kc,n)
     &                                      + voffy(j)*lc_yslope(ic,jc,kc,n)
     &                                      + voffz(k)*lc_zslope(ic,jc,kc,n)
              end do
            end do
          end do
        end do

      else if(ncby.ge.ncbx.and.ncby.ge.ncbz)then

c=============== CASE 2: y direction is long direction ======================

c ... added to prevent underflow for small crse values

        do n = 1, nvar 
          do k = cslopelo(3)-1,cslopehi(3)+1
            do i = cslopelo(1)-1, cslopehi(1)+1 
              do j = cslopelo(2)-1,cslopehi(2)+1
                crse(i,j,k,n) = cvmgt(crse(i,j,k,n),zero,abs(crse(i,j,k,n)).gt.1.0D-20)
              end do
            end do
          end do
        end do

c ... computed unlimited and limited slopes

        do n = 1, nvar 

c ... --> in x direction

          do k=cslopelo(3), cslopehi(3)
            do i=cslopelo(1), cslopehi(1)
              do j=cslopelo(2), cslopehi(2)
                uc_xslope(i,j,k,n) = half*(crse(i+1,j,k,n)-crse(i-1,j,k,n))

                cen  = uc_xslope(i,j,k,n)
                forw = two*(crse(i+1,j,k,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i-1,j,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_xslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

              end do
            end do
          end do
          if (xok) then
            if (bclo(1,n) .eq. EXT_DIR .or. bclo(1,n).eq.HOEXTRAP) then
              i = cslopelo(1)
              do k=cslopelo(3), cslopehi(3)
                do j=cslopelo(2), cslopehi(2)
                  uc_xslope(i,j,k,n)  = -sixteen/fifteen*crse(i-1,j,k,n) 
     &                        + half*crse(i,j,k,n)
     $                        + two3rd*crse(i+1,j,k,n) - tenth*crse(i+2,j,k,n)

                  cen  = uc_xslope(i,j,k,n)
                  forw = two*(crse(i+1,j,k,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i-1,j,k,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_xslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                end do
              end do
            end if
            if (bchi(1,n) .eq. EXT_DIR .or. bchi(1,n).eq.HOEXTRAP) then
              i = cslopehi(1)
              do k=cslopelo(3), cslopehi(3)
                do j=cslopelo(2), cslopehi(2)
                  uc_xslope(i,j,k,n) = sixteen/fifteen*crse(i+1,j,k,n) 
     &                      - half*crse(i,j,k,n)
     $                      - two3rd*crse(i-1,j,k,n) + tenth*crse(i-2,j,k,n)

                  cen  = uc_xslope(i,j,k,n)
                  forw = two*(crse(i+1,j,k,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i-1,j,k,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_xslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                end do
              end do
            end if
          end if

c ... --> in y direction

          do k=cslopelo(3), cslopehi(3)
            do i=cslopelo(1), cslopehi(1)
              do j=cslopelo(2), cslopehi(2)
                uc_yslope(i,j,k,n) = half*(crse(i,j+1,k,n)-crse(i,j-1,k,n))

                cen  = uc_yslope(i,j,k,n)
                forw = two*(crse(i,j+1,k,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i,j-1,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_yslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

              end do
            end do
          end do
          if (yok) then
            if (bclo(2,n) .eq. EXT_DIR .or. bclo(2,n).eq.HOEXTRAP) then
              j = cslopelo(2)
              do k=cslopelo(3), cslopehi(3)
                do i=cslopelo(1), cslopehi(1)
                  uc_yslope(i,j,k,n)  = -sixteen/fifteen*crse(i,j-1,k,n) 
     &                        + half*crse(i,j,k,n)
     $                        + two3rd*crse(i,j+1,k,n) - tenth*crse(i,j+2,k,n)

                  cen  = uc_yslope(i,j,k,n)
                  forw = two*(crse(i,j+1,k,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i,j-1,k,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_yslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                end do
              end do
            end if
            if (bchi(2,n) .eq. EXT_DIR .or. bchi(2,n).eq.HOEXTRAP) then
              j = cslopehi(2)
              do k=cslopelo(3), cslopehi(3)
                do i=cslopelo(1), cslopehi(1)
                  uc_yslope(i,j,k,n) = sixteen/fifteen*crse(i,j+1,k,n) 
     &                      - half*crse(i,j,k,n)
     $                      - two3rd*crse(i,j-1,k,n) + tenth*crse(i,j-2,k,n)

                  cen  = uc_yslope(i,j,k,n)
                  forw = two*(crse(i,j+1,k,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i,j-1,k,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_yslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                end do
              end do
            end if
          end if

c ... --> in z direction

          do k=cslopelo(3), cslopehi(3)
            do i=cslopelo(1), cslopehi(1)
              do j=cslopelo(2), cslopehi(2)
                uc_zslope(i,j,k,n) = half*(crse(i,j,k+1,n)-crse(i,j,k-1,n))

                cen  = uc_zslope(i,j,k,n)
                forw = two*(crse(i,j,k+1,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i,j,k-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_zslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

               end do
            end do
          end do
          if (zok) then
            if (bclo(3,n) .eq. EXT_DIR .or. bclo(3,n).eq.HOEXTRAP) then
              k = cslopelo(3)
              do i=cslopelo(1), cslopehi(1)
                do j=cslopelo(2), cslopehi(2)
                  uc_zslope(i,j,k,n)  = -sixteen/fifteen*crse(i,j,k-1,n) 
     &                        + half*crse(i,j,k,n)
     $                        + two3rd*crse(i,j,k+1,n) - tenth*crse(i,j,k+2,n)

                  cen  = uc_zslope(i,j,k,n)
                  forw = two*(crse(i,j,k+1,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i,j,k-1,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_zslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                 end do
              end do
            end if
            if (bchi(3,n) .eq. EXT_DIR .or. bchi(3,n).eq.HOEXTRAP) then
              k = cslopehi(3)
              do i=cslopelo(1), cslopehi(1)
                do j=cslopelo(2), cslopehi(2)
                  uc_zslope(i,j,k,n) = sixteen/fifteen*crse(i,j,k+1,n) 
     &                      - half*crse(i,j,k,n)
     $                      - two3rd*crse(i,j,k-1,n) + tenth*crse(i,j,k-2,n)

                  cen  = uc_zslope(i,j,k,n)
                  forw = two*(crse(i,j,k+1,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i,j,k-1,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_zslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                end do
              end do
            end if
          end if
        end do

        if (lin_limit.eq.1)then

c ... compute linear limited slopes
c     Note that the limited and the unlimited slopes
c     have the same sign, and it is assumed that they do.

c ... --> compute slope factors

          do k=cslopelo(3), cslopehi(3)
            do i=cslopelo(1), cslopehi(1)
              do j=cslopelo(2), cslopehi(2)
                xslope_factor(i,j,k) = one
                yslope_factor(i,j,k) = one
                zslope_factor(i,j,k) = one
              end do
            end do
          end do

          do n = 1, nvar 
            do k=cslopelo(3), cslopehi(3)
              do i=cslopelo(1), cslopehi(1)
                do j=cslopelo(2), cslopehi(2)
                  denom = uc_xslope(i,j,k,n)
                  denom = cvmgt(denom,one,denom.ne.zero)
                  factorn = lc_xslope(i,j,k,n)/denom
                  factorn = cvmgt(one,factorn,denom.eq.zero)
                  xslope_factor(i,j,k) = min(xslope_factor(i,j,k),factorn)

                  denom = uc_yslope(i,j,k,n)
                  denom = cvmgt(denom,one,denom.ne.zero)
                  factorn = lc_yslope(i,j,k,n)/denom
                  factorn = cvmgt(one,factorn,denom.eq.zero)
                  yslope_factor(i,j,k) = min(yslope_factor(i,j,k),factorn)

                  denom = uc_zslope(i,j,k,n)
                  denom = cvmgt(denom,one,denom.ne.zero)
                  factorn = lc_zslope(i,j,k,n)/denom
                  factorn = cvmgt(one,factorn,denom.eq.zero)
                  zslope_factor(i,j,k) = min(zslope_factor(i,j,k),factorn)
                end do
              end do
            end do
          end do

c ... --> compute linear limited slopes

          do n = 1, nvar 
            do k=cslopelo(3), cslopehi(3)
              do i=cslopelo(1), cslopehi(1)
                do j=cslopelo(2), cslopehi(2)
                  lc_xslope(i,j,k,n) = xslope_factor(i,j,k)*uc_xslope(i,j,k,n)
                  lc_yslope(i,j,k,n) = yslope_factor(i,j,k)*uc_yslope(i,j,k,n)
                  lc_zslope(i,j,k,n) = zslope_factor(i,j,k)*uc_zslope(i,j,k,n)
                end do
              end do
            end do
          end do
        end if

c ... do the interpolation

        do n = 1, nvar
          do k = fblo(3), fbhi(3)
            kc = IX_PROJ(k,lratioz)
            do i = fblo(1), fbhi(1)
              ic = IX_PROJ(i,lratiox)
              do j = fblo(2), fbhi(2)
                jc = IX_PROJ(j,lratioy)
                fine(i,j,k,n) = crse(ic,jc,kc,n) + voffx(i)*lc_xslope(ic,jc,kc,n)
     &                                      + voffy(j)*lc_yslope(ic,jc,kc,n)
     &                                      + voffz(k)*lc_zslope(ic,jc,kc,n)
              end do
            end do
          end do
        end do

      else if(ncbz.ge.ncbx.and.ncbz.ge.ncby)then

c=============== CASE 3: z direction is long direction ======================

c ... added to prevent underflow for small crse values

        do n = 1, nvar 
          do j = cslopelo(2)-1,cslopehi(2)+1
            do i = cslopelo(1)-1, cslopehi(1)+1 
              do k = cslopelo(3)-1,cslopehi(3)+1
                crse(i,j,k,n) = cvmgt(crse(i,j,k,n),zero,abs(crse(i,j,k,n)).gt.1.0D-20)
              end do
            end do
          end do
        end do

c ... computed unlimited and limited slopes

        do n = 1, nvar 

c ... --> in x direction

          do j=cslopelo(2), cslopehi(2)
            do i=cslopelo(1), cslopehi(1)
              do k=cslopelo(3), cslopehi(3)
                uc_xslope(i,j,k,n) = half*(crse(i+1,j,k,n)-crse(i-1,j,k,n))

                cen  = uc_xslope(i,j,k,n)
                forw = two*(crse(i+1,j,k,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i-1,j,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_xslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

               end do
            end do
          end do
          if (xok) then
            if (bclo(1,n) .eq. EXT_DIR .or. bclo(1,n).eq.HOEXTRAP) then
              i = cslopelo(1)
              do j=cslopelo(2), cslopehi(2)
                do k=cslopelo(3), cslopehi(3)
                  uc_xslope(i,j,k,n)  = -sixteen/fifteen*crse(i-1,j,k,n) 
     &                        + half*crse(i,j,k,n)
     $                        + two3rd*crse(i+1,j,k,n) - tenth*crse(i+2,j,k,n)

                  cen  = uc_xslope(i,j,k,n)
                  forw = two*(crse(i+1,j,k,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i-1,j,k,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_xslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                 end do
              end do
            end if
            if (bchi(1,n) .eq. EXT_DIR .or. bchi(1,n).eq.HOEXTRAP) then
              i = cslopehi(1)
              do j=cslopelo(2), cslopehi(2)
                do k=cslopelo(3), cslopehi(3)
                  uc_xslope(i,j,k,n) = sixteen/fifteen*crse(i+1,j,k,n) 
     &                      - half*crse(i,j,k,n)
     $                      - two3rd*crse(i-1,j,k,n) + tenth*crse(i-2,j,k,n)

                  cen  = uc_xslope(i,j,k,n)
                  forw = two*(crse(i+1,j,k,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i-1,j,k,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_xslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                end do
              end do
            end if
          end if

c ... --> in y direction

          do j=cslopelo(2), cslopehi(2)
            do i=cslopelo(1), cslopehi(1)
              do k=cslopelo(3), cslopehi(3)
                uc_yslope(i,j,k,n) = half*(crse(i,j+1,k,n)-crse(i,j-1,k,n))

                cen  = uc_yslope(i,j,k,n)
                forw = two*(crse(i,j+1,k,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i,j-1,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_yslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

               end do
            end do
          end do
          if (yok) then
            if (bclo(2,n) .eq. EXT_DIR .or. bclo(2,n).eq.HOEXTRAP) then
              j = cslopelo(2)
              do i=cslopelo(1), cslopehi(1)
                do k=cslopelo(3), cslopehi(3)
                  uc_yslope(i,j,k,n)  = -sixteen/fifteen*crse(i,j-1,k,n) 
     &                        + half*crse(i,j,k,n)
     $                        + two3rd*crse(i,j+1,k,n) - tenth*crse(i,j+2,k,n)

                  cen  = uc_yslope(i,j,k,n)
                  forw = two*(crse(i,j+1,k,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i,j-1,k,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_yslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                 end do
              end do
            end if
            if (bchi(2,n) .eq. EXT_DIR .or. bchi(2,n).eq.HOEXTRAP) then
              j = cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                do k=cslopelo(3), cslopehi(3)
                  uc_yslope(i,j,k,n) = sixteen/fifteen*crse(i,j+1,k,n) 
     &                      - half*crse(i,j,k,n)
     $                      - two3rd*crse(i,j-1,k,n) + tenth*crse(i,j-2,k,n)
 
                  cen  = uc_yslope(i,j,k,n)
                  forw = two*(crse(i,j+1,k,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i,j-1,k,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_yslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

               end do
              end do
            end if
          end if

c ... --> in z direction

          do j=cslopelo(2), cslopehi(2)
            do i=cslopelo(1), cslopehi(1)
              do k=cslopelo(3), cslopehi(3)
                uc_zslope(i,j,k,n) = half*(crse(i,j,k+1,n)-crse(i,j,k-1,n))

                cen  = uc_zslope(i,j,k,n)
                forw = two*(crse(i,j,k+1,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i,j,k-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = cvmgp(slp,zero,forw*back)
                lc_zslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

               end do
            end do
          end do
          if (zok) then
            if (bclo(3,n) .eq. EXT_DIR .or. bclo(3,n).eq.HOEXTRAP) then
              k = cslopelo(3)
              do i=cslopelo(1), cslopehi(1)
                do j=cslopelo(2), cslopehi(2)
                  uc_zslope(i,j,k,n)  = -sixteen/fifteen*crse(i,j,k-1,n) 
     &                        + half*crse(i,j,k,n)
     $                        + two3rd*crse(i,j,k+1,n) - tenth*crse(i,j,k+2,n)

                  cen  = uc_zslope(i,j,k,n)
                  forw = two*(crse(i,j,k+1,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i,j,k-1,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_zslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                 end do
              end do
            end if
            if (bchi(3,n) .eq. EXT_DIR .or. bchi(3,n).eq.HOEXTRAP) then
              k = cslopehi(3)
              do i=cslopelo(1), cslopehi(1)
                do j=cslopelo(2), cslopehi(2)
                  uc_zslope(i,j,k,n) = sixteen/fifteen*crse(i,j,k+1,n) 
     &                      - half*crse(i,j,k,n)
     $                      - two3rd*crse(i,j,k-1,n) + tenth*crse(i,j,k-2,n)

                  cen  = uc_zslope(i,j,k,n)
                  forw = two*(crse(i,j,k+1,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i,j,k-1,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = cvmgp(slp,zero,forw*back)
                  lc_zslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))

                end do
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
              do k=cslopelo(3), cslopehi(3)
                xslope_factor(i,j,k) = one
                yslope_factor(i,j,k) = one
                zslope_factor(i,j,k) = one
              end do
            end do
          end do

          do n = 1, nvar 
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                do k=cslopelo(3), cslopehi(3)
                  denom = uc_xslope(i,j,k,n)
                  denom = cvmgt(denom,one,denom.ne.zero)
                  factorn = lc_xslope(i,j,k,n)/denom
                  factorn = cvmgt(one,factorn,denom.eq.zero)
                  xslope_factor(i,j,k) = min(xslope_factor(i,j,k),factorn)

                  denom = uc_yslope(i,j,k,n)
                  denom = cvmgt(denom,one,denom.ne.zero)
                  factorn = lc_yslope(i,j,k,n)/denom
                  factorn = cvmgt(one,factorn,denom.eq.zero)
                  yslope_factor(i,j,k) = min(yslope_factor(i,j,k),factorn)

                  denom = uc_zslope(i,j,k,n)
                  denom = cvmgt(denom,one,denom.ne.zero)
                  factorn = lc_zslope(i,j,k,n)/denom
                  factorn = cvmgt(one,factorn,denom.eq.zero)
                  zslope_factor(i,j,k) = min(zslope_factor(i,j,k),factorn)
                end do
              end do
            end do
          end do

c ... --> compute linear limited slopes

          do n = 1, nvar 
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                do k=cslopelo(3), cslopehi(3)
                  lc_xslope(i,j,k,n) = xslope_factor(i,j,k)*uc_xslope(i,j,k,n)
                  lc_yslope(i,j,k,n) = yslope_factor(i,j,k)*uc_yslope(i,j,k,n)
                  lc_zslope(i,j,k,n) = zslope_factor(i,j,k)*uc_zslope(i,j,k,n)
                end do
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
              do k = fblo(3), fbhi(3)
                kc = IX_PROJ(k,lratioz)
                fine(i,j,k,n) = crse(ic,jc,kc,n) + voffx(i)*lc_xslope(ic,jc,kc,n)
     &                                      + voffy(j)*lc_yslope(ic,jc,kc,n)
     &                                      + voffz(k)*lc_zslope(ic,jc,kc,n)
              end do
            end do
          end do
        end do
      end if

      end
c ::: 
c ::: --------------------------------------------------------------
c ::: 

      subroutine FORT_CQINTERP (fine, DIMS(fine), 
     $                          DIMS(fb),
     $                          nvar, lratiox, lratioy, lratioz, crse,
     $                          clo, chi, DIMS(cb),
     $		                fslo, fshi, cslope, clen, fslope, fdat,
     $                          flen, voff, bc, limslope,
     $                          fvcx, fvcy, fvcz, cvcx, cvcy, cvcz)

      integer DIMDEC(fine)
      integer DIMDEC(fb)
      integer DIMDEC(cb)
      integer fslo(3), fshi(3)
      integer nvar, lratiox, lratioy, lratioz
      integer bc(3,2,nvar)
      integer clen, flen, clo, chi, limslope
      REAL_T fine(DIMV(fine),nvar)
      REAL_T crse(clo:chi, nvar)
      REAL_T cslope(clo:chi, 3)
      REAL_T fslope(flen, 3)
      REAL_T fdat(flen)
      REAL_T voff(flen)
      REAL_T fvcx(fb_l1:fb_h1+1)
      REAL_T fvcy(fb_l2:fb_h2+1)
      REAL_T fvcz(fb_l3:fb_h3+1)
      REAL_T cvcx(cb_l1:cb_h1+1)
      REAL_T cvcy(cb_l2:cb_h2+1)
      REAL_T cvcz(cb_l3:cb_h3+1)

      print *,'QUADRATIC INTERP NOT IMPLEMEMNTED IN 3-D'
      stop
      end

# if 0
c 
c -----------------------------------------------------------------
c THIS IS A SCALAR VERSION OF THE ABOVE CODE
c -----------------------------------------------------------------
c
      subroutine FORT_CCINTERP2 (fine, DIMS(fine),
     $                           fb_l1, fb_l2, fb_l3,
     $                           fb_h1, fb_h2, fb_h3,
     $                           nvar, lratiox, lratioy, lratioz, crse,
     $                           clo, chi, cb_l1, cb_l2, cb_l3,
     $                           cb_h1, cb_h2, cb_h3,
     $		                 fslo, fshi, cslope, clen, fslope, fdat,
     $                           flen, voff, bc, limslope)

      integer DIMDEC(fine)
      integer fb_l1, fb_l2, fb_l3
      integer fb_h1, fb_h2, fb_h3
      integer cb_l1, cb_l2, cb_l3
      integer cb_h1, cb_h2, cb_h3
      integer fslo(3), fshi(3)
      integer bc(3,2,nvar)
      integer lratiox, lratioy, lratioz, nvar, clen, flen, clo, chi
      integer limslope
      REAL_T fine(DIMV(fine),nvar)
      REAL_T crse(cb_l1-1:cb_h1+1, cb_l2-1:cb_h2+1, 
     &            cb_l3-1:cb_h3+1, nvar )
      REAL_T cslope(cb_l1-1:cb_h1+1, cb_l2-1:cb_h2+1, 
     &              cb_l3-1:cb_h3+1, 3 )
      REAL_T fslope(flen, 3)
      REAL_T fdat(flen)
      REAL_T voff(flen)

#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)

c ::: local var
      integer n, fn
      integer i, ii, ic, ioff
      integer j, jj, jc, joff
      integer k, kk, kc, koff
      REAL_T hafratx, hafraty, hafratz, volratiox, volratioy, volratioz
      REAL_T cen, forw, back, slp
      REAL_T xoff, yoff, zoff
      REAL_T cx, cy, cz
      REAL_T sgn
      integer ilo, ihi, jlo, jhi, klo, khi

c     :::::: helpful statement functions
      REAL_T  slplox, slphix, slploy, slphiy, slploz, slphiz
      REAL_T  c1, c2, c3, c4

      slplox(i,j,k,n) =  - c1*crse(i-1,j,k,n)
     $		         + c2*crse(i  ,j,k,n)
     $                   + c3*crse(i+1,j,k,n)
     $			 - c4*crse(i+2,j,k,n)
      slphix(i,j,k,n) =    c1*crse(i+1,j,k,n)
     $		         - c2*crse(i  ,j,k,n)
     $                   - c3*crse(i-1,j,k,n)
     $			 + c4*crse(i-2,j,k,n)
      slploy(i,j,k,n) =  - c1*crse(i,j-1,k,n)
     $		         + c2*crse(i,j  ,k,n)
     $                   + c3*crse(i,j+1,k,n)
     $			 - c4*crse(i,j+2,k,n)
      slphiy(i,j,k,n) =    c1*crse(i,j+1,k,n)
     $		         - c2*crse(i,j  ,k,n)
     $                   - c3*crse(i,j-1,k,n)
     $			 + c4*crse(i,j-2,k,n)
      slploz(i,j,k,n) =  - c1*crse(i,j,k-1,n)
     $		         + c2*crse(i,j,k  ,n)
     $                   + c3*crse(i,j,k+1,n)
     $			 - c4*crse(i,j,k+2,n)
      slphiz(i,j,k,n) =    c1*crse(i,j,k+1,n)
     $		         - c2*crse(i,j,k  ,n)
     $                   - c3*crse(i,j,k-1,n)
     $			 + c4*crse(i,j,k-2,n)
      
      c1 = sixteen/fifteen
      c2 = half
      c3 = two3rd
      c4 = tenth

      hafratx = half*dfloat(lratiox-1)
      hafraty = half*dfloat(lratioy-1)
      hafratz = half*dfloat(lratioz-1)

      volratiox = one/dfloat(lratiox)
      volratioy = one/dfloat(lratioy)
      volratioz = one/dfloat(lratioz)

      do n = 1, nvar
      do kc = cb_l3, cb_h3
      do jc = cb_l2, cb_h2
      do ic = cb_l1, cb_h1

c        ::::: compute x slopes
	 if (limslope .ne. 0) then
            cen  = half*(crse(ic+1,jc,kc,n)-crse(ic-1,jc,kc,n))
            forw = two*(crse(ic+1,jc,kc,n)-crse(ic,jc,kc,n))
            back = two*(crse(ic,jc,kc,n) - crse(ic-1,jc,kc,n))
	    slp  = min(abs(forw),abs(back))
	    slp  = cvmgp(slp,zero,forw*back)
	    sgn  = sign(one,cen)
            cx   = sgn*min(slp,abs(cen))
            if (ic.eq.cb_l1 .and. (bclo(1,n) .eq. EXT_DIR 
     &	        .or. bclo(1,n).eq.HOEXTRAP)) then
	        cen  = slplox(ic,jc,kc,n)
                cx   = sgn*min(slp,abs(cen))
            end if
            if (ic.eq.cb_h1 .and. (bchi(1,n) .eq. EXT_DIR 
     &          .or. bchi(1,n).eq.HOEXTRAP)) then
                cen  = slphix(ic,jc,kc,n)
                cx   = sgn*min(slp,abs(cen))
            end if
	 else
	    cx = half*(crse(ic+1,jc,kc,n)-crse(ic-1,jc,kc,n))
            if (ic.eq.cb_l1 .and. (bclo(1,n) .eq. EXT_DIR 
     &	        .or. bclo(1,n).eq.HOEXTRAP)) then
	        cx  = slplox(ic,jc,kc,n)
            end if
            if (ic.eq.cb_h1 .and. (bchi(1,n) .eq. EXT_DIR 
     &          .or. bchi(1,n).eq.HOEXTRAP)) then
                cx  = slphix(ic,jc,kc,n)
            end if
	 end if

c	 ::::: slopes in the Y direction
	 if (limslope .ne. 0) then
            cen  = half*(crse(ic,jc+1,kc,n)-crse(ic,jc-1,kc,n))
            forw = two*(crse(ic,jc+1,kc,n)-crse(ic,jc,kc,n))
            back = two*(crse(ic,jc,kc,n) - crse(ic,jc-1,kc,n))
	    slp  = min(abs(forw),abs(back))
	    slp  = cvmgp(slp,zero,forw*back)
	    sgn  = sign(one,cen)
            cy   = sgn*min(slp,abs(cen))
            if (jc.eq.cb_l2 .and. (bclo(2,n) .eq. EXT_DIR 
     &	        .or. bclo(2,n).eq.HOEXTRAP)) then
	        cen  = slploy(ic,jc,kc,n)
                cy   = sgn*min(slp,abs(cen))
            end if
            if (jc.eq.cb_h2 .and. (bchi(2,n) .eq. EXT_DIR 
     &          .or. bchi(2,n).eq.HOEXTRAP)) then
                cen  = slphiy(ic,jc,kc,n)
                cy   = sgn*min(slp,abs(cen))
            end if
	 else
	    cy = half*(crse(ic,jc+1,kc,n)-crse(ic,jc-1,kc,n))
            if (jc.eq.cb_l2 .and. (bclo(2,n) .eq. EXT_DIR 
     &	        .or. bclo(2,n).eq.HOEXTRAP)) then
	        cy   = slploy(ic,jc,kc,n)
            end if
            if (ic.eq.cb_h2 .and. (bchi(2,n) .eq. EXT_DIR 
     &          .or. bchi(2,n).eq.HOEXTRAP)) then
                cy   = slphiy(ic,jc,kc,n)
            end if
	 end if

c	 ::::: slopes in the Z direction
	 if (limslope .ne. 0) then
            cen  = half*(crse(ic,jc,kc+1,n)-crse(ic,jc,kc-1,n))
            forw = two*(crse(ic,jc,kc+1,n)-crse(ic,jc,kc,n))
            back = two*(crse(ic,jc,kc,n) - crse(ic,jc,kc-1,n))
	    slp  = min(abs(forw),abs(back))
	    slp  = cvmgp(slp,zero,forw*back)
	    sgn  = sign(one,cen)
            cz   = sgn*min(slp,abs(cen))
            if (kc.eq.cb_l3 .and. (bclo(3,n) .eq. EXT_DIR 
     &	        .or. bclo(3,n).eq.HOEXTRAP)) then
	        cen  = slploz(ic,jc,kc,n)
                cz   = sgn*min(slp,abs(cen))
            end if
            if (kc.eq.cb_h3 .and. (bchi(3,n) .eq. EXT_DIR 
     &          .or. bchi(3,n).eq.HOEXTRAP)) then
                cen  = slphiz(ic,jc,kc,n)
                cz   = sgn*min(slp,abs(cen))
            end if
	 else
	    cz = half*(crse(ic,jc,kc+1,n)-crse(ic,jc,kc-1,n))
            if (kc.eq.cb_l3 .and. (bclo(3,n) .eq. EXT_DIR 
     &	        .or. bclo(3,n).eq.HOEXTRAP)) then
	        cz   = slploz(ic,jc,kc,n)
            end if
            if (kc.eq.cb_h3 .and. (bchi(3,n) .eq. EXT_DIR 
     &          .or. bchi(3,n).eq.HOEXTRAP)) then
                cz   = slphiz(ic,jc,kc,n)
            end if
	 end if

c	 ::::: now interpolate to fine grid
	 ii  = lratiox*ic
	 ilo = max(ii,fb_l1) - ii
	 ihi = min(ii+lratiox-1,fb_h1) - ii
         jj  = lratioy*jc
	 jlo = max(jj,fb_l2) - jj
	 jhi = min(jj+lratioy-1,fb_h2) - jj
	 kk  = lratioz*kc
	 klo = max(kk,fb_l3) - kk
	 khi = min(kk+lratioz-1,fb_h2) - kk

	 do koff = klo, khi
	    k = lratioz*kc + koff
	    zoff = dfloat(koff)-hafratz
	    do joff = jlo, jhi
	       j = lratioy*jc + joff
	       yoff = dfloat(joff)-hafraty
	       do ioff = ilo, ihi
	          i = lratiox*ic + ioff
		  xoff = dfloat(ioff)-hafratx
		  fine(i,j,k,n) = crse(ic,jc,kc,n) + 
     &				( volratiox*xoff*cx + volratioy*yoff*cy 
     &                          + volratioz*zoff*cz )
	       end do
	    end do
	 end do
 
      end do
      end do
      end do
      end do

      end
#endif


c ::: 
c ::: --------------------------------------------------------------
c ::: pcinterp:  cell centered piecewise constant interpolation
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fine grid array
c ::: DIMS(fine)   =>  (const)  index limits of fine grid
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: 
c ::: crse         =>  (const)  coarse grid data 
c ::: DIMS(crse)   =>  (const)  index limits of coarse grid
c ::: cblo,cbhi    =>  (const) coarse grid region containing fblo,fbhi
c ::: 
c ::: longdir      =>  (const)  which index direction is longest (1 or 2)
c ::: lratio(3)    =>  (const)  refinement ratio between levels
c ::: nvar         =>  (const)  number of components in array
c ::: 
c ::: TEMPORARY ARRAYS
c ::: ftmp         =>  1-D temp array
c ::: --------------------------------------------------------------
c ::: 
      subroutine FORT_PCINTERP (crse,DIMS(crse),cblo,cbhi,
     &                          fine,DIMS(fine),fblo,fbhi,
     &                          longdir,lratiox,lratioy,lratioz,nvar,
     &                          ftmp, ftmp_lo, ftmp_hi)
      integer DIMDEC(crse)
      integer cblo(3), cbhi(3)
      integer DIMDEC(fine)
      integer fblo(3), fbhi(3)
      integer nvar, lratiox, lratioy, lratioz, longdir
      integer ftmp_lo, ftmp_hi
      REAL_T  crse(DIMV(crse), nvar)
      REAL_T  fine(DIMV(fine), nvar)
      REAL_T  ftmp(ftmp_lo:ftmp_hi)

      integer i, j, k, ic, jc, kc, ioff, joff, koff, n
      integer istrt, istop, jstrt, jstop, kstrt, kstop
      integer ilo, ihi, jlo, jhi, klo, khi

      if (longdir .eq. 1) then
         do n = 1, nvar
	 do kc = cblo(3), cbhi(3)
	    kstrt = kc*lratioz
	    kstop = (kc+1)*lratioz - 1
	    klo = max(fblo(3),kstrt)
	    khi = min(fbhi(3),kstop)
            do jc = cblo(2), cbhi(2)

c	       ::::: fill strip in i direction
	       do ioff = 0, lratiox-1
	          do ic = cblo(1), cbhi(1)
	             i = lratiox*ic + ioff
	             ftmp(i) = crse(ic,jc,kc,n)
                  end do
	       end do

c	       ::::: stuff into fine array
	       jstrt = jc*lratioy
	       jstop = (jc+1)*lratioy - 1
	       jlo = max(fblo(2),jstrt)
	       jhi = min(fbhi(2),jstop)
	       do k = klo, khi
	       do j = jlo, jhi
	       do i = fblo(1), fbhi(1)
	          fine(i,j,k,n) = ftmp(i)
	       end do
	       end do
	       end do
	    end do
	 end do
	 end do
      else if (longdir.eq.2) then
         do n = 1, nvar
	 do kc = cblo(3), cbhi(3)
	    kstrt = kc*lratioz
	    kstop = (kc+1)*lratioz - 1
	    klo = max(fblo(3),kstrt)
	    khi = min(fbhi(3),kstop)
            do ic = cblo(1), cbhi(1)

c	       ::::: fill strip in j direction
	       do joff = 0, lratioy-1
	          do jc = cblo(2), cbhi(2)
	             j = lratioy*jc + joff
	             ftmp(j) = crse(ic,jc,kc,n)
                  end do
	       end do

c	       ::::: stuff into fine array
	       istrt = ic*lratiox
	       istop = (ic+1)*lratiox - 1
	       ilo = max(fblo(1),istrt)
	       ihi = min(fbhi(1),istop)
	       do k = klo, khi
	       do i = ilo, ihi
	       do j = fblo(2), fbhi(2)
	          fine(i,j,k,n) = ftmp(j)
	       end do
	       end do
	       end do
	    end do
	 end do
	 end do
      else
         do n = 1, nvar
	 do ic = cblo(1), cbhi(1)
	    istrt = ic*lratiox
	    istop = (ic+1)*lratiox - 1
	    ilo = max(fblo(1),istrt)
	    ihi = min(fbhi(1),istop)
            do jc = cblo(2), cbhi(2)

c	       ::::: fill strip in k direction
	       do koff = 0, lratioz-1
	          do kc = cblo(3), cbhi(3)
	             k = lratioz*kc + koff
	             ftmp(k) = crse(ic,jc,kc,n)
                  end do
	       end do

c	       ::::: stuff into fine array
	       jstrt = jc*lratioy
	       jstop = (jc+1)*lratioy - 1
	       jlo = max(fblo(2),jstrt)
	       jhi = min(fbhi(2),jstop)
	       do i = ilo, ihi
	       do j = jlo, jhi
	       do k = fblo(3), fbhi(3)
	          fine(i,j,k,n) = ftmp(k)
	       end do
	       end do
	       end do
	    end do
	 end do
	 end do
      end if

      end





