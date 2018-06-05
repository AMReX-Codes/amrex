
module amrex_interp_module

  use amrex_fort_module
  use amrex_constants_module
  use amrex_bc_types_module

  implicit none

contains

#define IX_PROJ(A,B) (A+B*iabs(A))/B-iabs(A)

! ::: --------------------------------------------------------------
! ::: nbinterp:  node based bilinear interpolation
! ::: 
! ::: INPUTS/OUTPUTS
! ::: fine        <=>  (modify) fine grid array
! ::: fine_l1,fine_h1   =>  (const)  index limits of fine grid
! ::: fb_l1,fb_h1     =>  (const)  subregion of fine grid to get values
! ::: 
! ::: crse         =>  (const)  coarse grid data widened by 1 zone
! ::: crse_l1,crse_h1   =>  (const)  index limits of coarse grid
! ::: 
! ::: lratio       =>  (const)  refinement ratio between levels
! ::: nvar         =>  (const)  number of components in array
! ::: num_slp      =>  (const)  number of types of slopes
! ::: 
! :::  ::: TEMPORARY ARRAYS
! ::: sl           =>  num_slp 1-D slope arrays
! ::: --------------------------------------------------------------
! ::: 
    subroutine AMREX_NBINTERP (crse, crse_l1,crse_h1, cb_l1,cb_h1, &
                              fine, fine_l1,fine_h1, fb_l1,fb_h1, &
                              lratio, nvar, &
                              sl, num_slp, &
                              actual_comp,actual_state) bind(c,name='amrex_nbinterp')

      implicit none

      integer crse_l1,crse_h1
      integer cb_l1,cb_h1
      integer fine_l1,fine_h1
      integer fb_l1,fb_h1
      integer lratio, nvar
      integer num_slp
      integer actual_comp,actual_state
      real(amrex_real)  fine(fine_l1:fine_h1, nvar)
      real(amrex_real)  crse(crse_l1:crse_h1, nvar)
      real(amrex_real)    sl(cb_l1:cb_h1,num_slp)
      real(amrex_real) strip(cb_l1*lratio:cb_h1*lratio)

#define  SLX 1
#define  SLY 2
#define  SLXY 3

      ! local var
      integer lx, ly, lz
      integer i,j,k,ii,jj,kk,n
      integer ibeg,iend,jstrs,jends,jbeg,jend
      integer lys,lye
      real(amrex_real) invratio

      invratio = one/dble(lratio)
      ibeg = max( cb_l1*lratio, fine_l1 )
      iend = min( cb_h1*lratio, fine_h1 )

      do 100 n = 1, nvar 
          ! first fill a strip that will fit
          do i = cb_l1, cb_h1-1 
            sl(i,SLX) = invratio*(crse(i+1,n)-crse(i,n))
          enddo
          i = cb_h1

            do lx = 0, lratio-1
              do i = cb_l1, cb_h1-1 
                ii = i*lratio + lx
                strip(ii) = crse(i,n) + dble(lx)*sl(i,SLX) 
	      enddo
	    enddo
            i = cb_h1
            ii = i*lratio
            strip(ii) = crse(i,n) 
            ! copy on intersection
            do i = ibeg,iend 
              fine(i,n) = strip(i)
            enddo
100   continue

    end subroutine AMREX_NBINTERP

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
                x = denom*(two*lratio + one)
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
! ::: ccinterp:   conservative interpolation from coarse grid to
! ::: subregion of fine grid defined by (fblo,fbhi)
! ::: 
! ::: Inputs/Outputs
! ::: fine        <=>  (modify) fine grid array
! ::: flo,fhi      =>  (const)  index limits of fine grid
! ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
! ::: nvar         =>  (const)  number of variables in state vector
! ::: lratio       =>  (const)  refinement ratio between levels
! ::: 
! ::: crse         =>  (const)  coarse grid data widended by 1 zone
! ::: and unrolled
! ::: clo,chi      =>  (const)  one dimensional limits of crse grid
! ::: cslo,cshi    =>  (const)  coarse grid index limits where
! :::				slopes are to be defined. This is
! :::				the projection of (fblo,fbhi) down
! :::				to the coarse level 
! ::: fslo,fshi    =>  (const)  fine grid index limits where
! :::				slopes are needed.  This is the
! :::				refinement of (cslo,cshi) and
! :::				contains but may not be identical
! :::				to (fblo,fbhi).
! ::: cslope       =>  (modify) temp array coarse grid slopes
! ::: clen         =>  (const)  length of coarse gtid slopes
! ::: fslope       =>  (modify) temp array for fine grid slope
! ::: flen         =>  (const)  length of fine grid slope array
! ::: fdat         =>  (const)  temp array for fine grid data
! ::: limslope     =>  (const)  != 0 => limit slopes
! :::
! ::: NOTE: data must be sent in so that 
! :::	    cslope(1,*) and crse(1,*) are associated with
! :::	    the same cell
! :::
! ::: EXAMPLE:
! ::: Suppose the patch called "fine" has index extent:
! ::: 
! ::: floi1 = 3, fhii1 = 12
! ::: floi2 = 8, fhii2 = 20
! ::: 
! ::: suppose the subergion of this patch that is to be filled 
! ::: by interpolation has index extent:
! ::: 
! ::: fblo(1) = 5, fbhi(1) = 10
! ::: fblo(2) = 13, fbhi(2) = 20
! ::: 
! ::: suppose the refinement ratio is 2
! ::: 
! ::: Then the coarsening of this subregion (to level 0) is
! ::: 
! ::: cb_l1 = 2  cb_h1 = 5         (ncbx = 4)
! ::: cb_l2 = 6  cb_h2 = 10        (ncby = 5)
! ::: 
! ::: In order to compute slopes, we need one extra row of
! ::: coarse grid zones:
! ::: 
! ::: cslo(1) = 1  cshi(1) = 6         (ncsx = 6)
! ::: cslo(2) = 5  cshi(2) = 11        (ncsy = 7)
! ::: 
! ::: This is the size of the coarse grid array of data that filpatch 
! ::: has filled at level 0.
! ::: The "cslope" and "crse" arrays are this size.
! ::: 
! ::: In order to unroll the slope calculation we make these arrays look
! ::: like 1-D arrays.  The mapping from 2-D to 1-D is as follows:
! ::: 
! ::: The point (cb_l(1),cb_l(2)) -> 1
! ::: The point (cslo(1),cslo(2)) -> clo = 1 - 1 - ncsx = -6
! ::: 
! ::: The point (cb_h1,cb_h2) -> clen = ncby*ncsx - 2 = 5*6-2 = 28
! ::: The point (cshi(1),cshi(2)) -> chi = clo + ncsx*ncsy - 1 
! :::                                    =  -6 +    6*7    - 1 = 35
! ::: 
! :::      -------------------------------------------------
! :::      |       |       |       |       |       |  chi  |  
! :::  11  |   30  |   31  |   32  |   33  |   34  |   35  |   cshi(2)
! :::      |       |       |       |       |       |       |
! :::      -------------------------------------------------
! :::      |       |       |       |       |  clen |       |  
! :::  10  |   24  |   25  |   26  |   27  |   28  |   29  |   cb_h(2)
! :::      |       |       |       |       |       |       |
! :::      -------------------------------------------------
! :::      |       |       |       |       |       |       |  
! :::   9  |   18  |   19  |   20  |   21  |   22  |   23  |  
! :::      |       |       |       |       |       |       |
! :::      -------------------------------------------------
! :::      |       |       |       |       |       |       |  
! :::   8  |   12  |   13  |   14  |   15  |   16  |   17  |  
! :::      |       |       |       |       |       |       |
! :::      -------------------------------------------------
! :::      |       |       |       |       |       |       |  
! :::   7  |    6  |    7  |    8  |    9  |   10  |   11  |  
! :::      |       |       |       |       |       |       |
! :::      -------------------------------------------------
! :::      |       |       |       |       |       |       |  
! :::   6  |    0  |    1  |    2  |    3  |    4  |    5  |   cb_l(2)
! :::      |       |       |       |       |       |       |
! :::      -------------------------------------------------
! :::      |  clo  |       |       |       |       |       |  
! :::   5  |   -6  |   -5  |   -4  |   -3  |   -2  |   -1  |   cslo(2)
! :::      |       |       |       |       |       |       |
! :::      -------------------------------------------------
! :::          1       2       3       4       5       6
! :::               cb_l1                   cb_h1
! :::       cslo(1)                                 cshi(1)
! ::: 
! ::: 
! ::: In the 1-D coordinates:
! :::    ist = 1    = stride in I direction
! :::    jst = 6    = stride in J direction  (ncsx)
! ::: 
! ::: --------------------------------------------------------------
! ::: 
#if 0
    subroutine AMREX_CCINTERP (fine, fine_l1,fine_h1, &
                              fb_l1, fb_h1, &
                              nvar, lratio, crse, clo, chi, &
                              cb_l1, cb_h1, &
                              fslo, fshi, cslope, clen, fslope, fdat, &
                              flen, voff, bc, limslope, &
                              fvcx, cvcx, &
                              actual_comp,actual_state) bind(c,name='amrex_ccinterp')

      implicit none

      integer fine_l1,fine_h1
      integer fslo(1), fshi(1)
      integer fb_l1, fb_h1
      integer cb_l1, cb_h1
      integer clo, chi
      integer lratio, nvar, clen, flen, limslope
      integer bc(1,2,nvar)
      integer actual_comp,actual_state
      real(amrex_real) fine(fine_l1:fine_h1,nvar)
      real(amrex_real) crse(clo:chi, nvar)
      real(amrex_real) cslope(clo:chi, 2)
      real(amrex_real) fslope(flen, 2)
      real(amrex_real) fdat(flen)
      real(amrex_real) voff(flen)
      real(amrex_real) fvcx(fb_l1:fb_h1+1)
      real(amrex_real) cvcx(cb_l1:cb_h1+1)

#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)

      ! local var
      integer n, fn
      integer i, ic, ioff
      integer j, jc, joff
      integer ist, jst
      real(amrex_real) hafrat, volratio
      real(amrex_real) cen, forw, back, slp, sgn
      real(amrex_real) fcen, ccen
      real(amrex_real) xoff, yoff
      integer ncbx, ncby
      integer ncsx, ncsy
      integer islo, jslo
      integer icc, istart, iend
      integer lenx, leny, maxlen
      logical xok, yok

      hafrat = half*dble(lratio-1)
      volratio = one/dble(lratio)

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

          ! compute slopes in x direction
          do 220 i = 1, clen
              cen = half*(crse(i+ist,n)-crse(i-ist,n))
              forw = crse(i+ist,n)-crse(i,n)
              back = crse(i,n)-crse(i-ist,n)
              slp = sign(one,cen)*min(abs(cen),abs(forw),abs(back))
              cslope(i,1)=merge(slp,zero,forw*back>=zero)
220       continue

          ! strip out a fine grid slope vector
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

    end subroutine AMREX_CCINTERP

#endif




# if 1

    subroutine AMREX_CCINTERP (fine, fine_l1,fine_h1, &
                              fb_l1, fb_h1, &
                              nvar, lratio, crse, clo, chi, &
                              cb_l1, cb_h1, &
                              fslo, fshi, cslope, clen, fslope, fdat, &
                              flen, voff, bc, limslope, &
                              fvcx, cvcx, &
                              cmax, cmin, alpha, &
                              actual_comp,actual_state) bind(c,name='amrex_ccinterp')

      implicit none

      integer fine_l1,fine_h1
      integer fslo(1), fshi(1)
      integer fb_l1, fb_h1
      integer cb_l1, cb_h1
      integer clo, chi
      integer lratio, nvar, clen, flen, limslope
      integer actual_comp,actual_state
      integer bc(1,2,nvar)
      real(amrex_real) fine(fine_l1:fine_h1,nvar)
      real(amrex_real) crse(clo:chi, nvar)
      real(amrex_real) cslope(clo:chi, 2)
      real(amrex_real) fslope(flen, 2)
      real(amrex_real) fdat(flen)
      real(amrex_real) voff(flen)
      real(amrex_real) fvcx(fb_l1:fb_h1+1)
      real(amrex_real) cvcx(cb_l1:cb_h1+1)
      real(amrex_real) cmax, cmin, alpha

#if 0
#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)
#endif

      ! local var
      integer n, fn
      integer i, ic, ioff
      integer j, jc, joff
      integer ist, jst
      real(amrex_real) hafrat, volratio
      real(amrex_real) cen, forw, back, slp, sgn
      real(amrex_real) fcen, ccen
      real(amrex_real) xoff, yoff
      integer ncbx, ncby
      integer ncsx, ncsy
      integer islo, jslo
      integer icc, istart, iend
      integer lenx, leny, maxlen
      logical xok, yok

      hafrat = half*dble(lratio-1)
      volratio = one/dble(lratio)

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


      ! added to prevent underflow for small crse values
      do n = 1, nvar 
        do i = clo, chi 
          crse(i,n) = merge(crse(i,n),zero,abs(crse(i,n)).gt.1.0d-50)
        enddo
      enddo

      do 290 n = 1, nvar 

         ! compute slopes in x direction
         if (limslope .ne. 0) then
            do i = 1, clen 
               cen = half*(crse(i+ist,n)-crse(i-ist,n))
               forw = two*(crse(i+ist,n)-crse(i,n))
               back = two*(crse(i,n)-crse(i-ist,n))
               slp  = min(abs(forw),abs(back))
               slp  = merge(slp,zero,forw*back>=zero)
               cslope(i,1)=sign(one,cen)*min(slp,abs(cen))
            enddo
            if (xok) then
!               if (bclo(1,n) .eq. amrex_bc_ext_dir .or. bclo(1,n).eq.amrex_bc_hoextrap) then
               if (bc(1,1,n) .eq. amrex_bc_ext_dir .or. bc(1,1,n).eq.amrex_bc_hoextrap) then
                  do i = 1, clen, jst 
                     cen  = -sixteen/fifteen*crse(i-ist,n) + half*crse(i,n) &
                          + two3rd*crse(i+ist,n) - tenth*crse(i+2*ist,n)
                     sgn  = sign(one,crse(i+ist,n)-crse(i-ist,n))
                     forw = two*(crse(i+ist,n)-crse(i,n))
                     back = two*(crse(i,n)-crse(i-ist,n))
                     slp  = min(abs(forw),abs(back))
                     slp  = merge(slp,zero,forw*back>=zero)
                     cslope(i,1)=sgn*min(slp,abs(cen))
                  enddo
               endif
               if (bc(1,2,n) .eq. amrex_bc_ext_dir .or. bc(1,2,n).eq.amrex_bc_hoextrap) then
                  do i = ncbx, clen, jst 
                     cen = sixteen/fifteen*crse(i+ist,n) - half*crse(i,n) &
                          - two3rd*crse(i-ist,n) + tenth*crse(i-2*ist,n)
                     sgn  = sign(one,crse(i+ist,n)-crse(i-ist,n))
                     forw = two*(crse(i+ist,n)-crse(i,n))
                     back = two*(crse(i,n)-crse(i-ist,n))
                     slp  = min(abs(forw),abs(back))
                     slp  = merge(slp,zero,forw*back>=zero)
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
!               if (bclo(1,n) .eq. amrex_bc_ext_dir .or. bclo(1,n).eq.amrex_bc_hoextrap) then
               if (bc(1,1,n) .eq. amrex_bc_ext_dir .or. bc(1,1,n).eq.amrex_bc_hoextrap) then
                  do i = 1, clen, jst 
                     cen  = -sixteen/fifteen*crse(i-ist,n) + half*crse(i,n) &
                          + two3rd*crse(i+ist,n) - tenth*crse(i+2*ist,n)
                     cslope(i,1)=cen
                  enddo
               endif
               if (bc(1,2,n) .eq. amrex_bc_ext_dir .or. bc(1,2,n).eq.amrex_bc_hoextrap) then
                  do i = ncbx, clen, jst 
                     cen = sixteen/fifteen*crse(i+ist,n) - half*crse(i,n) &
                          - two3rd*crse(i-ist,n) + tenth*crse(i-2*ist,n)
                     cslope(i,1)=cen
                  enddo
               endif
            endif
         endif

            ! strip out a fine grid slope vector
            do 370 ioff = 1, lratio 
              icc = clo + ist 
              istart = ioff
              iend = ioff + (ncbx-1)*lratio
              do 380 fn = istart, iend, lratio 
                fslope(fn,1) = cslope(icc,1)
!                fslope(fn,2) = cslope(icc,2)
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

    end subroutine AMREX_CCINTERP

#endif

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
! ::: ratio        =>  (const)  refinement ratio between levels
! ::: nvar         =>  (const)  number of components in array
! ::: 
! ::: TEMPORARY ARRAYS
! ::: ftmp         =>  1-D temp array
! ::: --------------------------------------------------------------
! ::: 
    subroutine AMREX_PCINTERP (crse,crse_l1,crse_h1,cblo,cbhi, &
                              fine,fine_l1,fine_h1,fblo,fbhi, &
                              longdir,lratio,nvar,ftmp,ftmp_lo,ftmp_hi, &
                              actual_comp,actual_state) bind(c,name='amrex_pcinterp')

      implicit none

      integer crse_l1,crse_h1
      integer cblo(1), cbhi(1)
      integer fine_l1,fine_h1
      integer fblo(1), fbhi(1)
      integer ftmp_lo, ftmp_hi
      integer nvar, lratio, longdir
      integer actual_comp,actual_state
      real(amrex_real)  crse(crse_l1:crse_h1, nvar)
      real(amrex_real)  fine(fine_l1:fine_h1, nvar)
      real(amrex_real)  ftmp(ftmp_lo:ftmp_hi)

      ! Local variables    
      real(amrex_real) sumrho
      integer i, ic, ioff, n

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

    end subroutine AMREX_PCINTERP

! ::: 
! ::: --------------------------------------------------------------
! ::: linccinterp:   linear conservative interpolation from coarse grid to
! ::: subregion of fine grid defined by (fblo,fbhi)
! ::: 
! ::: The interpolation is linear in that it uses a
! ::: a limiting scheme that preserves the value of 
! ::: any linear combination of the
! ::: coarse grid data components--e.g.,
! ::: if sum_ivar a(ic,ivar)*fab(ic,ivar) = 0, then
! ::: sum_ivar a(ic,ivar)*fab(if,ivar) = 0 is satisfied
! ::: in all fine cells if covering coarse cell ic.
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
! ::: lratio(1)    =>  (const)  refinement ratio between levels
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
     subroutine AMREX_LINCCINTERP (fine, fine_l1,fine_h1, fblo, fbhi, &
                                  fvcb_l1,fvcb_h1, &
                                  crse, crse_l1,crse_h1, cvcb_l1,cvcb_h1, &
                                  uc_xslope, lc_xslope, xslope_factor, &
                                  cslope_l1,cslope_h1, &
                                  cslopelo, cslopehi, &
                                  nvar, lratiox, &
                                  bc, lim_slope, lin_limit, &
                                  fvcx, cvcx, &
                                  voffx, alpha, cmax, cmin, &
                                  actual_comp,actual_state) bind(c,name='amrex_linccinterp')

       implicit none

       integer fine_l1,fine_h1
       integer crse_l1,crse_h1
       integer fvcb_l1,fvcb_h1
       integer cvcb_l1,cvcb_h1
       integer cslope_l1,cslope_h1
       integer fblo(1), fbhi(1)
       integer cslopelo(1), cslopehi(1)
       integer lratiox, nvar
       integer lim_slope, lin_limit
       integer bc(1,2,nvar)
       integer actual_comp,actual_state
       real(amrex_real) fine(fine_l1:fine_h1,nvar)
       real(amrex_real) crse(crse_l1:crse_h1, nvar)
       real(amrex_real) uc_xslope(cslope_l1:cslope_h1,nvar)
       real(amrex_real) lc_xslope(cslope_l1:cslope_h1,nvar)
       real(amrex_real) xslope_factor(cslope_l1:cslope_h1)
       real(amrex_real) alpha(cslope_l1:cslope_h1,nvar)
       real(amrex_real) cmax(cslope_l1:cslope_h1,nvar)
       real(amrex_real) cmin(cslope_l1:cslope_h1,nvar)
       real(amrex_real) fvcx(fvcb_l1:fvcb_h1)
       real(amrex_real) voffx(fvcb_l1:fvcb_h1)
       real(amrex_real) cvcx(cvcb_l1:cvcb_h1)

#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)

       integer n
       integer i, ic
       real(amrex_real) cen, forw, back, slp
       real(amrex_real) factorn, denom
       real(amrex_real) fxcen, cxcen, fycen, cycen
       real(amrex_real) orig_corr_fact,corr_fact
       real(amrex_real) dummy_fine
       logical xok, yok
       integer ncbx, ncby
       integer ioff
       integer voff_lo(1),voff_hi(1)

       ncbx = cslopehi(1)-cslopelo(1)+1

       voff_lo(1) = cslopelo(1) * lratiox
       voff_hi(1) = (cslopehi(1)+1) * lratiox - 1

       xok = (ncbx .ge. 2)

       do i = voff_lo(1),voff_hi(1)
          ic = IX_PROJ(i,lratiox)
          fxcen = half*(fvcx(i)+fvcx(i+1))
          cxcen = half*(cvcx(ic)+cvcx(ic+1))
          voffx(i) = (fxcen-cxcen)/(cvcx(ic+1)-cvcx(ic))
       end do

       do n = 1, nvar

! ...     Prevent underflow for small crse values.

          do i = cslopelo(1)-1, cslopehi(1)+1
             crse(i,n) = merge(crse(i,n),zero,abs(crse(i,n)).gt.1.0d-50)
          end do

! ...     Initialize alpha = 1 and define cmax and cmin as neighborhood max/mins.

          do i = cslopelo(1), cslopehi(1)
             alpha(i,n) = 1.d0
             cmax(i,n) = crse(i,n)
             cmin(i,n) = crse(i,n)
             do ioff = -1,1
               cmax(i,n) = max(cmax(i,n),crse(i+ioff,n))
               cmin(i,n) = min(cmin(i,n),crse(i+ioff,n))
             end do
          end do

       end do

! ...  Compute unlimited and limited slopes

       do n = 1, nvar

          do i=cslopelo(1), cslopehi(1)
             uc_xslope(i,n) = half*(crse(i+1,n)-crse(i-1,n))
             cen  = uc_xslope(i,n)
             forw = two*(crse(i+1,n)-crse(i,n))
             back = two*(crse(i,n)-crse(i-1,n))
             slp  = min(abs(forw),abs(back))
             slp  = merge(slp,zero,forw*back>=zero)
             lc_xslope(i,n)=sign(one,cen)*min(slp,abs(cen))
          end do

          if (bclo(1,n) .eq. amrex_bc_ext_dir .or. bclo(1,n).eq.amrex_bc_hoextrap) then
            i = cslopelo(1)
            if (xok) then
                uc_xslope(i,n)  = -sixteen/fifteen*crse(i-1,n) &
                    + half*crse(i,n) &
                     + two3rd*crse(i+1,n) - tenth*crse(i+2,n)
            else
                uc_xslope(i,n)  = fourth * ( &
                  crse(i+1,n) + five*crse(i,n) - six*crse(i-1,n) )
            endif
            cen  = uc_xslope(i,n)
            forw = two*(crse(i+1,n)-crse(i,n))
            back = two*(crse(i,n)-crse(i-1,n))
            slp  = min(abs(forw),abs(back))
            slp  = merge(slp,zero,forw*back>=zero)
            lc_xslope(i,n)=sign(one,cen)*min(slp,abs(cen))
          end if

          if (bchi(1,n) .eq. amrex_bc_ext_dir .or. bchi(1,n).eq.amrex_bc_hoextrap) then
            i = cslopehi(1)
            if (xok) then
                uc_xslope(i,n) = sixteen/fifteen*crse(i+1,n) &
                     - half*crse(i,n) &
                     - two3rd*crse(i-1,n) + tenth*crse(i-2,n)
            else
                uc_xslope(i,n) = -fourth * ( &
                  crse(i-1,n) + five*crse(i,n) - six*crse(i+1,n) )
            endif
            cen  = uc_xslope(i,n)
            forw = two*(crse(i+1,n)-crse(i,n))
            back = two*(crse(i,n)-crse(i-1,n))
            slp  = min(abs(forw),abs(back))
            slp  = merge(slp,zero,forw*back>=zero)
            lc_xslope(i,n)=sign(one,cen)*min(slp,abs(cen))
          end if

       end do

       if (lim_slope.eq.0) then

! ...    Do the interpolation using unlimited slopes.

          do n = 1, nvar
             do i = fblo(1), fbhi(1)
                ic = IX_PROJ(i,lratiox)
                fine(i,n) = crse(ic,n) &
                     + voffx(i)*uc_xslope(ic,n)
             end do
          end do

       else 

         if (lin_limit.eq.1) then

! ...      compute linear limited slopes
!          Note that the limited and the unlimited slopes
!          have the same sign, and it is assumed that they do.

! ... --> compute slope factors

           do i=cslopelo(1), cslopehi(1)
             xslope_factor(i) = one
           end do

           do n = 1, nvar
             do i=cslopelo(1), cslopehi(1)
                denom = uc_xslope(i,n)
                denom = merge(denom,one,denom.ne.zero)
                factorn = lc_xslope(i,n)/denom
                factorn = merge(one,factorn,denom.eq.zero)
                xslope_factor(i) = min(xslope_factor(i),factorn)
             end do
           end do

! ... -->  compute linear limited slopes

           do n = 1, nvar
             do i=cslopelo(1), cslopehi(1)
                lc_xslope(i,n) = xslope_factor(i)*uc_xslope(i,n)
             end do
           end do

         else

!          Limit slopes so as to not introduce new maxs or mins.

            do n = 1, nvar
               do i = voff_lo(1),voff_hi(1)
                  ic = IX_PROJ(i,lratiox)

                     orig_corr_fact = voffx(i)*lc_xslope(ic,n)
                     dummy_fine = crse(ic,n) + orig_corr_fact
                     if ( (dummy_fine .gt. cmax(ic,n)) .and. &
                          (abs(orig_corr_fact) .gt. 1.e-10*abs(crse(ic,n)))) then
                        corr_fact = (cmax(ic,n) - crse(ic,n)) / orig_corr_fact
                        alpha(ic,n) = min(alpha(ic,n),corr_fact)
                     endif
                     if ( (dummy_fine .lt. cmin(ic,n)) .and. &
                          (abs(orig_corr_fact) .gt. 1.e-10*abs(crse(ic,n)))) then
                        corr_fact = (cmin(ic,n) - crse(ic,n)) / orig_corr_fact
                        alpha(ic,n) = min(alpha(ic,n),corr_fact)
                     endif

#ifndef NDEBUG
                     if (alpha(ic,n) .lt. 0.d0) then
                        print *,'OOPS - ALPHA SHOULD BE POSITIVE IN CCINTERP '
                        print *,'ALPHA = ',alpha(ic,n)
                        print *,'AT (I,N) = ',ic,n
                        print *,'ORIG_CORR_FACT = ',orig_corr_fact
                        call bl_abort(" ")
                     endif
                     if (alpha(ic,n) .gt. 1.d0) then
                        print *,'OOPS - ALPHA SHOULD BE <= 1.0 IN CCINTERP '
                        print *,'ALPHA = ',alpha(ic,n)
                        print *,'AT (I,N) = ',ic,n
                        print *,'ORIG_CORR_FACT = ',orig_corr_fact
                        call bl_abort(" ")
                     endif
#endif
               end do
            end do

         end if

! ...    Do the interpolation with limited slopes.

          do n = 1, nvar
            do i = fblo(1), fbhi(1)
               ic = IX_PROJ(i,lratiox)
               fine(i,n) = crse(ic,n) + alpha(ic,n)* &
                    voffx(i)*lc_xslope(ic,n)
            end do
          end do

       end if

     end subroutine AMREX_LINCCINTERP

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
!      todo
          write(6,*) 'AMREX_QUARTINTERP: refinement ratio = 4 TODO'
          stop
       else
          write(6,*) 'AMREX_QUARTINTERP: unsupported refinement ratio'
          stop
       endif

     end subroutine AMREX_QUARTINTERP

end module amrex_interp_module
