#undef  BL_LANG_CC
#define BL_LANG_FORT

//#include "REAL.H"
#include "CONSTANTS.H"
#include "BC_TYPES.H"
#include "INTERP_F.H"
#include "ArrayLim.H"

#include "ONE_D.H"


#define IX_PROJ(A,B) (A+B*iabs(A))/B-iabs(A)
#define SDIM BL_SPACEDIM

c ::: --------------------------------------------------------------
c ::: nbinterp:  node based bilinear interpolation
c ::: 
c ::: INPUTS/OUTPUTS
c ::: fine        <=>  (modify) fine grid array
c ::: DIMS(fine)   =>  (const)  index limits of fine grid
c ::: DIMS(fb)     =>  (const)  subregion of fine grid to get values
c ::: 
c ::: crse         =>  (const)  coarse grid data widened by 1 zone
c ::: DIMS(crse)   =>  (const)  index limits of coarse grid
c ::: 
c ::: lratio       =>  (const)  refinement ratio between levels
c ::: nvar         =>  (const)  number of components in array
c ::: num_slp      =>  (const)  number of types of slopes
c ::: strip_lo,strip_hi  =>  (const)  index limits of strip array
c ::: 
c :::  ::: TEMPORARY ARRAYS
c ::: sl           =>  num_slp 1-D slope arrays
c ::: strip        =>  1-D temp array
c ::: --------------------------------------------------------------
c ::: 
      subroutine FORT_NBINTERP (crse, DIMS(crse), DIMS(cb),
     $                          fine, DIMS(fine), DIMS(fb),
     $		                lratio, nvar,
     $                          sl, num_slp, strip, strip_lo, strip_hi)

      integer DIMDEC(crse)
      integer DIMDEC(cb)
      integer DIMDEC(fine)
      integer DIMDEC(fb)
      integer lratio, nvar
      integer num_slp
      integer strip_lo, strip_hi
      REAL_T  fine(DIMV(fine), nvar)
      REAL_T  crse(DIMV(crse), nvar)
      REAL_T    sl(DIM1(cb),num_slp)
      REAL_T strip(strip_lo:strip_hi,0:lratio)

#define  SLX 1
#define  SLY 2
#define  SLXY 3

c ::: local var
      integer lx, ly, lz
      integer i,j,k,ii,jj,kk,n
      integer ibeg,iend,jstrs,jends,jbeg,jend
      integer lys,lye
      REAL_T invratio

      invratio = one/dfloat(lratio)
      ibeg = max( ARG_L1(cb)*lratio, ARG_L1(fine) )
      iend = min( ARG_H1(cb)*lratio, ARG_H1(fine) )

      do 100 n = 1, nvar 
c         ::: first fill a strip that will fit
          do i = ARG_L1(cb), ARG_H1(cb)-1 
            sl(i,SLX) = invratio*(crse(i+1,n)-crse(i,n))
          enddo
          i = ARG_H1(cb)

            do lx = 0, lratio-1
              do i = ARG_L1(cb), ARG_H1(cb)-1 
                ii = i*lratio + lx
                strip(ii,0) = crse(i,n)
     $		             + dfloat(lx)*sl(i,SLX) 
	      enddo
	    enddo
            i = ARG_H1(cb)
            ii = i*lratio
            strip(ii,0) = crse(i,n) 
c ::: copy on intersection
            do i = ibeg,iend 
              fine(i,n) = strip(i,0)
            enddo
100   continue
      return
      end

#undef  SLX
#undef  SLY
#undef  SLXY

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
c ::: lratio       =>  (const)  refinement ratio between levels
c ::: nvar         =>  (const)  number of components in array
c ::: 
c ::: TEMPORARY ARRAYS
c ::: slx,sly,slxy =>  1-D slope arrays
c ::: strip        =>  1-D temp array
c ::: --------------------------------------------------------------
c ::: 
      subroutine FORT_CBINTERP (crse, DIMS(crse), DIMS(cb),
     $                          fine, DIMS(fine), DIMS(fb),
     $		                lratio, nvar,
     $                          sl, num_slp, strip, strip_lo, strip_hi)

      integer DIMDEC(crse)
      integer DIMDEC(cb)
      integer DIMDEC(fine)
      integer DIMDEC(fb)
      integer lratio, nvar
      integer num_slp
      integer strip_lo, strip_hi
      REAL_T  fine(DIMV(fine), nvar)
      REAL_T  crse(DIMV(crse), nvar)
      REAL_T    sl(DIM1(cb),num_slp)
      REAL_T strip(strip_lo:strip_hi)

#define SLX 1
#define SLY 2
#define SLXY 3

c ::: local var
      integer lx, ly
      integer hrat, ic, jc, jfn, jfc, i, j, n
      REAL_T x, y
      REAL_T denom

      denom = one/dfloat(2*lratio)
      hrat = lratio/2
      do 200 n = 1, nvar 
c ::: first fill a strip that will fit
          do ic = ARG_L1(cb),ARG_H1(cb)-1
            sl(ic,SLX) = crse(ic+1,n)-crse(ic,n)
	  enddo

            do lx = 0, lratio-1
              do ic = ARG_L1(cb), ARG_H1(cb)-1
                i = ic*lratio + lx
                x = denom*(two*lratio + one)
                strip(i) = crse(ic,n) + x*sl(ic,SLX) 
	      enddo
	    enddo

c ::: stuff into output array
            do i = ARG_L1(fb), ARG_H1(fb) 
              fine(i,n) = strip(i-hrat)
            enddo
230       continue
200   continue

      return
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
c ::: lratio       =>  (const)  refinement ratio between levels
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
#if 0
      subroutine FORT_CCINTERP (fine, DIMS(fine), 
     $                          fb_l1, fb_h1, 
     $                          nvar, lratio, crse, clo, chi, 
     $                          cb_l1, cb_h1,
     $		                fslo, fshi, cslope, clen, fslope, fdat,
     $                          flen, voff, bc, limslope,
     $                          fvcx, cvcx)


      integer DIMDEC(fine)
      integer fslo(1), fshi(1)
      integer fb_l1, fb_h1
      integer cb_l1, cb_h1
      integer clo, chi
      integer lratio, nvar, clen, flen, limslope
      integer bc(1,2,nvar)
      REAL_T fine(DIMV(fine),nvar)
      REAL_T crse(clo:chi, nvar)
      REAL_T cslope(clo:chi, 2)
      REAL_T fslope(flen, 2)
      REAL_T fdat(flen)
      REAL_T voff(flen)
      REAL_T fvcx(fb_l1:fb_h1+1)
      REAL_T cvcx(cb_l1:cb_h1+1)

#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)

c ::: local var
c ::: local var
      integer n, fn
      integer i, ic, ioff
      integer j, jc, joff
      integer ist, jst
      REAL_T hafrat, volratio
      REAL_T cen, forw, back, slp, sgn
      REAL_T fcen, ccen
      REAL_T xoff, yoff
      integer ncbx, ncby
      integer ncsx, ncsy
      integer islo, jslo
      integer icc, istart, iend
      integer lenx, leny, maxlen
      logical xok, yok

      hafrat = half*dfloat(lratio-1)
      volratio = one/dfloat(lratio)

      ncbx = cb_h1-cb_l1+1
      ncsx = ncbx+2
      ist = 1
      do 200 i = fb_l1, fb_h1
          fn = i-fslo(1)+1
          ic = IX_PROJ(i,lratio)
          fcen = half*(fvcx(i)+fvcx(i+1))
          ccen = half*(cvcx(ic)+cvcx(ic+1))
          voff(fn) = (fcen-ccen)/(cvcx(ic+1)-cvcx(ic))
200   continue
      do 210 n = 1, nvar

c ::: ::::: compute slopes in x direction
          do 220 i = 1, clen
              cen = half*(crse(i+ist,n)-crse(i-ist,n))
              forw = crse(i+ist,n)-crse(i,n)
              back = crse(i,n)-crse(i-ist,n)
              slp = sign(one,cen)*min(abs(cen),abs(forw),abs(back))
              cslope(i,1)=cvmgp(slp,zero,forw*back)
220       continue

c ::: .....::::: strip out a fine grid slope vector
          do 230 ioff = 1, lratio
              icc = clo + ist
              istart = ioff
              iend = ioff + (ncbx-1)*lratio
              do 240 fn = istart, iend, lratio
                  fslope(fn,1) = cslope(icc,1)
                  fdat(fn) = crse(icc,n)
                  icc = icc + ist
240           continue
230       continue
          do 250 i = fb_l1, fb_h1
              fn = i-fslo(1)+1
              fine(i,n) = fdat(fn) + voff(fn)*fslope(fn,1)
250       continue
210   continue

      return
      end

#endif




# if 1

      subroutine FORT_CCINTERP (fine, DIMS(fine), 
     $                          fb_l1, fb_h1, 
     $                          nvar, lratio, crse, clo, chi, 
     $                          cb_l1, cb_h1,
     $		                fslo, fshi, cslope, clen, fslope, fdat,
     $                          flen, voff, bc, limslope,
     $                          fvcx, cvcx)

      integer DIMDEC(fine)
      integer fslo(1), fshi(1)
      integer fb_l1, fb_h1
      integer cb_l1, cb_h1
      integer clo, chi
      integer lratio, nvar, clen, flen, limslope
      integer bc(1,2,nvar)
      REAL_T fine(DIMV(fine),nvar)
      REAL_T crse(clo:chi, nvar)
      REAL_T cslope(clo:chi, 2)
      REAL_T fslope(flen, 2)
      REAL_T fdat(flen)
      REAL_T voff(flen)
      REAL_T fvcx(fb_l1:fb_h1+1)
      REAL_T cvcx(cb_l1:cb_h1+1)
#if 0
#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)
#endif

c ::: local var
      integer n, fn
      integer i, ic, ioff
      integer j, jc, joff
      integer ist, jst
      REAL_T hafrat, volratio
      REAL_T cen, forw, back, slp, sgn
      REAL_T fcen, ccen
      REAL_T xoff, yoff
      integer ncbx, ncby
      integer ncsx, ncsy
      integer islo, jslo
      integer icc, istart, iend
      integer lenx, leny, maxlen
      logical xok, yok

      hafrat = half*dfloat(lratio-1)
      volratio = one/dfloat(lratio)

      ncbx = cb_h1-cb_l1+1
      xok = (ncbx .ge. 2)
      ncsx = ncbx+2
      ist = 1
      islo = cb_l1-1
      jst = ncsx
      lenx = fb_h1-fb_l1+1
         do i = fb_l1, fb_h1 
          fn = i-fslo(1)+1
          ic = IX_PROJ(i,lratio)
          fcen = half*(fvcx(i)+fvcx(i+1))
          ccen = half*(cvcx(ic)+cvcx(ic+1))
          voff(fn) = (fcen-ccen)/(cvcx(ic+1)-cvcx(ic))
        enddo   


c ::: ::::: added to prevent underflow for small crse values
      do n = 1, nvar 
        do i = clo, chi 
          crse(i,n) = cvmgt(crse(i,n),zero,abs(crse(i,n)).gt.1.0e-20)
        enddo
      enddo

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
            enddo
            if (xok) then
c               if (bclo(1,n) .eq. EXT_DIR .or. bclo(1,n).eq.HOEXTRAP) then
               if (bc(1,1,n) .eq. EXT_DIR .or. bc(1,1,n).eq.HOEXTRAP) then
                  do i = 1, clen, jst 
                     cen  = -sixteen/fifteen*crse(i-ist,n) + half*crse(i,n)
     $                    + two3rd*crse(i+ist,n) - tenth*crse(i+2*ist,n)
                     sgn  = sign(one,crse(i+ist,n)-crse(i-ist,n))
                     forw = two*(crse(i+ist,n)-crse(i,n))
                     back = two*(crse(i,n)-crse(i-ist,n))
                     slp  = min(abs(forw),abs(back))
                     slp  = cvmgp(slp,zero,forw*back)
                     cslope(i,1)=sgn*min(slp,abs(cen))
                  enddo
               endif
               if (bc(1,2,n) .eq. EXT_DIR .or. bc(1,2,n).eq.HOEXTRAP) then
                  do i = ncbx, clen, jst 
                     cen = sixteen/fifteen*crse(i+ist,n) - half*crse(i,n)
     $                    - two3rd*crse(i-ist,n) + tenth*crse(i-2*ist,n)
                     sgn  = sign(one,crse(i+ist,n)-crse(i-ist,n))
                     forw = two*(crse(i+ist,n)-crse(i,n))
                     back = two*(crse(i,n)-crse(i-ist,n))
                     slp  = min(abs(forw),abs(back))
                     slp  = cvmgp(slp,zero,forw*back)
                     cslope(i,1)=sgn*min(slp,abs(cen))
                  enddo
               endif
            endif
         else
            do i = 1, clen 
               cen = half*(crse(i+ist,n)-crse(i-ist,n))
               cslope(i,1)=cen
            enddo
            if (xok) then
c               if (bclo(1,n) .eq. EXT_DIR .or. bclo(1,n).eq.HOEXTRAP) then
               if (bc(1,1,n) .eq. EXT_DIR .or. bc(1,1,n).eq.HOEXTRAP) then
                  do i = 1, clen, jst 
                     cen  = -sixteen/fifteen*crse(i-ist,n) + half*crse(i,n)
     $                    + two3rd*crse(i+ist,n) - tenth*crse(i+2*ist,n)
                     cslope(i,1)=cen
                  enddo
               endif
               if (bc(1,2,n) .eq. EXT_DIR .or. bc(1,2,n).eq.HOEXTRAP) then
                  do i = ncbx, clen, jst 
                     cen = sixteen/fifteen*crse(i+ist,n) - half*crse(i,n)
     $                    - two3rd*crse(i-ist,n) + tenth*crse(i-2*ist,n)
                     cslope(i,1)=cen
                  enddo
               endif
            endif
         endif

c           ::::: strip out a fine grid slope vector
            do 370 ioff = 1, lratio 
              icc = clo + ist 
              istart = ioff
              iend = ioff + (ncbx-1)*lratio
              do 380 fn = istart, iend, lratio 
                fslope(fn,1) = cslope(icc,1)
c                fslope(fn,2) = cslope(icc,2)
                fdat(fn) = crse(icc,n)
                icc = icc + ist
380           continue
370         continue

              do 400 i = fb_l1, fb_h1 
                fn = i-fslo(1)+1
                fine(i,n) = fdat(fn) + voff(fn)*fslope(fn,1)
400           continue
391         continue

290   continue

	

      return
      end
#endif

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
c ::: ratio        =>  (const)  refinement ratio between levels
c ::: nvar         =>  (const)  number of components in array
c ::: 
c ::: TEMPORARY ARRAYS
c ::: ftmp         =>  1-D temp array
c ::: --------------------------------------------------------------
c ::: 
      subroutine FORT_PCINTERP (crse,DIMS(crse),cblo,cbhi,
     &                          fine,DIMS(fine),fblo,fbhi,
     &                          longdir,lratio,nvar,ftmp,ftmp_lo,ftmp_hi)
      integer DIMDEC(crse)
      integer cblo(2), cbhi(2)
      integer DIMDEC(fine)
      integer fblo(2), fbhi(2)
      integer ftmp_lo, ftmp_hi
      integer nvar, lratio, longdir
      REAL_T  crse(DIMV(crse), nvar)
      REAL_T  fine(DIMV(fine), nvar)
      REAL_T  ftmp(ftmp_lo:ftmp_hi)
	REAL_T sumrho

      integer i, j, ic, jc, ioff, joff, n

         do n = 1, nvar
	    do ioff = 0, lratio-1
	       do ic = cblo(1), cbhi(1)
	          i = lratio*ic + ioff
	          ftmp(i) = crse(ic,n)
               enddo
	    enddo
	    do i = fblo(1), fbhi(1)
	       fine(i,n) = ftmp(i)
	    enddo
	 enddo

#if 0
	do i = fblo(1), fbhi(1)
	   sumrho = fine(i,5)+fine(i,8)
	   if(abs(sumrho-fine(i,1)) .gt. 1.d-15) then
		write(6,*)'  sum of rhos .ne. total '
   	   endif
	enddo
#endif

      return
      end
