c ::: SCCS stuff "@(#)FABUTIL_3D.F	3.1\t6/25/93"

#define   FORT_CINTERP   cinterp3d
#define   FORT_PCINTERP  pcinterp3d
#define   FORT_VCINTERP  vcinterp3d
#define   FORT_CARTGRIDMINMAX  cartgridminmax3d

#undef AMREX_XSDK

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"

#define IX_PROJ(A,B) (A+B*iabs(A))/B-iabs(A)


c ::: --------------------------------------------------------------
      subroutine FORT_CARTGRIDMINMAX (data, dlo1, dlo2, dlo3, dhi1,
     $                                dhi2, dhi3, vfracdata, vfeps,
     $                                dmin, dmax)
      implicit none

      integer dlo1, dlo2, dlo3, dhi1, dhi2, dhi3
      REAL_T data(dlo1:dhi1 ,dlo2:dhi2, dlo3:dhi3)
      REAL_T vfracdata(dlo1:dhi1 ,dlo2:dhi2, dlo3:dhi3)
      REAL_T vfeps, dmin, dmax

      integer i, j, k

      dmax = -bigreal
      dmin = bigreal
      do 410 k = dlo3, dhi3
        do 420 j = dlo2, dhi2
          do 430 i = dlo1, dhi1
              if ( (vfracdata(i,j,k).ge.vfeps)) then
                dmax = max(dmax,data(i,j,k))
                dmin = min(dmin,data(i,j,k))
              endif
430       continue
420     continue
410   continue

      return
      end


c ::: --------------------------------------------------------------
      subroutine FORT_CINTERP (fine,floi1,floi2,floi3,fhii1,fhii2,fhii3,
     $  fblo,fbhi,nvar,lratio, crse,clo,chi,cblo,cbhi,fslo,fshi, cslope,
     $ clen,fslope,fdat,flen,voff)

      integer floi1,floi2,floi3, fhii1,fhii2,fhii3
      integer fblo(3), fbhi(3)
      integer cblo(3), cbhi(3)
      integer fslo(3), fshi(3)
      integer lratio, nvar, clen, flen, clo, chi
      REAL_T fine(floi1 :fhii1 ,floi2 :fhii2 ,floi3 :fhii3, nvar)
      REAL_T crse(clo:chi, nvar)
      REAL_T cslope(clo:chi, 3)
      REAL_T fslope(flen, 3)
      REAL_T fdat(flen)
      REAL_T voff(flen)

c ::: NOTE: data must be sent in so that 
c ::: cslope(1,*) and crse(1,*) are associated with
c ::: the same cell

c ::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ::: cinterp:   conservative interpolation from coarse grid to
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
c ::: slopes are to be defined. This is
c ::: the projection of (fblo,fbhi) down
c ::: to the coarse level 
c ::: fslo,fshi    =>  (const)  fine grid index limits where
c ::: slopes are needed.  This is the
c ::: refinement of (cslo,cshi) and
c ::: contains but may not be identical
c ::: to (fblo,fbhi).
c ::: cslope       =>  (modify) temp array coarse grid slopes
c ::: clen         =>  (const)  length of coarse gtid slopes
c ::: fslope       =>  (modify) temp array for fine grid slope
c ::: flen         =>  (const)  length of fine grid slope array
c ::: fdat         =>  (const)  temp array for fine grid data
c ::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c ::: local var
      integer n, fn
      integer i, ic, ioff
      integer j, jc, joff
      integer k, kc, koff
      integer ist, jst, kst
      REAL_T hafrat, volratio
      REAL_T cen, forw, back, slp
      REAL_T xoff, yoff, zoff
      integer ncbx, ncby, ncbz
      integer ncsx, ncsy, ncsz
      integer islo, jslo, kslo
      integer icc, istart, iend
      integer lenx, leny, lenz, maxlen

      hafrat = half*float(lratio-1)
      volratio = one/float(lratio)

      ncbx = cbhi(1)-cblo(1)+1
      ncby = cbhi(2)-cblo(2)+1
      ncbz = cbhi(3)-cblo(3)+1
      ncsx = ncbx+2
      ncsy = ncby+2
      ncsz = ncbz+2
      ist = 1
      jst = ncsx
      kst = ncsx*ncsy
      islo = cblo(1)-1
      jslo = cblo(2)-1
      kslo = cblo(3)-1
      lenx = fbhi(1)-fblo(1)+1
      leny = fbhi(2)-fblo(2)+1
      lenz = fbhi(3)-fblo(3)+1
      maxlen = max(lenx,leny,lenz)
      if (maxlen .eq. lenx) then
          do 100 i = fblo(1), fbhi(1) 
              fn = i-fslo(1)+1
              ioff = mod(fn-1,lratio)
              voff(fn) = float(ioff)-hafrat
100       continue
      else if (maxlen .eq. leny) then
          do 110 j = fblo(2), fbhi(2) 
              fn = j-fslo(2)+1
              joff = mod(fn-1,lratio)
              voff(fn) = float(joff)-hafrat
110       continue
      else
          do 120 k = fblo(3), fbhi(3) 
              fn = k-fslo(3)+1
              koff = mod(fn-1,lratio)
              voff(fn) = float(koff)-hafrat
120       continue
      end if
      do 130 n = 1, nvar 

c ::: ::::: compute slopes in x direction
          do 140 i = 1, clen 
              cen = half*(crse(i+ist,n)-crse(i-ist,n))
              forw = crse(i+ist,n)-crse(i,n)
              back = crse(i,n)-crse(i-ist,n)
              slp = sign(one,cen)*min(abs(cen),abs(forw),abs(back))
              cslope(i,1)=merge(slp,zero,forw*back>=0.0d0)
140       continue
c ::: ::::: compute slopes in y direction
          do 150 i = 1, clen 
              cen = half*(crse(i+jst,n)-crse(i-jst,n))
              forw = crse(i+jst,n)-crse(i,n)
              back = crse(i,n)-crse(i-jst,n)
              slp = sign(one,cen)*min(abs(cen),abs(forw),abs(back))
              cslope(i,2)=merge(slp,zero,forw*back>=0.0d0)
150       continue
c ::: ::::: compute slopes in z direction
          do 160 i = 1, clen 
              cen = half*(crse(i+kst,n)-crse(i-kst,n))
              forw = crse(i+kst,n)-crse(i,n)
              back = crse(i,n)-crse(i-kst,n)
              slp = sign(one,cen)*min(abs(cen),abs(forw),abs(back))
              cslope(i,3)=merge(slp,zero,forw*back>=0.0d0)
160       continue

          if (maxlen .eq. lenx) then
              do 170 kc = cblo(3), cbhi(3) 
                  do 180 jc = cblo(2), cbhi(2) 

c ::: ..,.......::::: strip out a fine grid slope vector
                      do 190 ioff = 1, lratio 
                          icc = clo + ist + jst*(jc-jslo) + kst*(kc-
     $                    kslo)
                          istart = ioff
                          iend = ioff + (ncbx-1)*lratio
                          do 200 fn = istart, iend, lratio 
                              fslope(fn,1) = cslope(icc,1)
                              fslope(fn,2) = cslope(icc,2)
                              fslope(fn,3) = cslope(icc,3)
                              fdat(fn) = crse(icc,n)
                              icc = icc + ist
200                       continue
190                   continue

                      do 210 koff = 0, lratio-1 
                          k = lratio*kc + koff
                          if (k .lt. fblo(3)) then
                              goto 210
c                             --- next ---
                          end if
                          if (k .gt. fbhi(3)) then
                              goto 211
c                             --- break ---
                          end if
                          zoff = float(koff)-hafrat
                          do 220 joff = 0, lratio-1 
                              j = lratio*jc + joff
                              if (j .lt. fblo(2)) then
                                  goto 220
c                                 --- next ---
                              end if
                              if (j .gt. fbhi(2)) then
                                  goto 221
c                                 --- break ---
                              end if
                              yoff = float(joff)-hafrat

                              do 230 i = fblo(1), fbhi(1) 
                                  fn = i-fslo(1)+1
                                  fine(i,j,k,n) = fdat(fn) + volratio* 
     $                            (voff(fn)*fslope(fn,1)+yoff*fslope(fn,
     $                            2)+ zoff*fslope(fn,3))
230                           continue
220                       continue
221                       continue
210                   continue
211                   continue
180               continue
170           continue
          else if (maxlen .eq. leny) then
              do 240 kc = cblo(3), cbhi(3) 
                  do 250 ic = cblo(1), cbhi(1) 

c ::: ..,.......::::: strip out a fine grid slope vector
                      do 260 joff = 1, lratio 
                          icc = clo + ist*(ic-islo) + jst + kst*(kc-
     $                    kslo)
                          istart = joff
                          iend = joff + (ncby-1)*lratio
                          do 270 fn = istart, iend, lratio 
                              fslope(fn,1) = cslope(icc,1)
                              fslope(fn,2) = cslope(icc,2)
                              fslope(fn,3) = cslope(icc,3)
                              fdat(fn) = crse(icc,n)
                              icc = icc + jst
270                       continue
260                   continue

                      do 280 koff = 0, lratio-1 
                          k = lratio*kc + koff
                          if (k .lt. fblo(3)) then
                              goto 280
c                             --- next ---
                          end if
                          if (k .gt. fbhi(3)) then
                              goto 281
c                             --- break ---
                          end if
                          zoff = float(koff)-hafrat
                          do 290 ioff = 0, lratio-1 
                              i = lratio*ic + ioff
                              if (i .lt. fblo(1)) then
                                  goto 290
c                                 --- next ---
                              end if
                              if (i .gt. fbhi(1)) then
                                  goto 291
c                                 --- break ---
                              end if
                              xoff = float(ioff)-hafrat

                              do 300 j = fblo(2), fbhi(2) 
                                  fn = j-fslo(2)+1
                                  fine(i,j,k,n) = fdat(fn) + volratio* 
     $                            (xoff*fslope(fn,1)+voff(fn)*fslope(fn,
     $                            2)+ zoff*fslope(fn,3))
300                           continue
290                       continue
291                       continue
280                   continue
281                   continue
250               continue
240           continue
          else
              do 310 jc = cblo(2), cbhi(2) 
                  do 320 ic = cblo(1), cbhi(1) 

c ::: ..,.......::::: strip out a fine grid slope vector
                      do 330 koff = 1, lratio 
                          icc = clo + ist*(ic-islo) + jst*(jc-jslo) + 
     $                    kst
                          istart = koff
                          iend = koff + (ncbz-1)*lratio
                          do 340 fn = istart, iend, lratio 
                              fslope(fn,1) = cslope(icc,1)
                              fslope(fn,2) = cslope(icc,2)
                              fslope(fn,3) = cslope(icc,3)
                              fdat(fn) = crse(icc,n)
                              icc = icc + kst
340                       continue
330                   continue

                      do 350 joff = 0, lratio-1 
                          j = lratio*jc + joff
                          if (j .lt. fblo(2)) then
                              goto 350
c                             --- next ---
                          end if
                          if (j .gt. fbhi(2)) then
                              goto 351
c                             --- break ---
                          end if
                          yoff = float(joff)-hafrat
                          do 360 ioff = 0, lratio-1 
                              i = lratio*ic + ioff
                              if (i .lt. fblo(1)) then
                                  goto 360
c                                 --- next ---
                              end if
                              if (i .gt. fbhi(1)) then
                                  goto 361
c                                 --- break ---
                              end if
                              xoff = float(ioff)-hafrat

                              do 370 k = fblo(3), fbhi(3) 
                                  fn = k-fslo(3)+1
                                  fine(i,j,k,n) = fdat(fn) + volratio* 
     $                            (xoff*fslope(fn,1)+yoff*fslope(fn,2)+ 
     $                            voff(fn)*fslope(fn,3))
370                           continue
360                       continue
361                       continue
350                   continue
351                   continue
320               continue
310           continue
          end if
130   continue

      return
      end


c ::: --------------------------------------------------------------
      subroutine FORT_PCINTERP (fine,floi1,floi2,floi3,fhii1,fhii2,
     $ fhii3,fblo,fbhi,lrat,nvar,crse,cloi1,cloi2,cloi3,chii1,chii2,
     $ chii3,cblo,cbhi,temp,tloi,thii)

      integer floi1,floi2,floi3
      integer fhii1,fhii2,fhii3
      integer cloi1,cloi2,cloi3
      integer chii1,chii2,chii3

      integer fblo(3), fbhi(3)
      integer cblo(3), cbhi(3)
      integer lrat, nvar, tloi, thii
      REAL_T fine(floi1 :fhii1 ,floi2 :fhii2 ,floi3 :fhii3, nvar)
      REAL_T crse(cloi1 :chii1 ,cloi2 :chii2 ,cloi3 :chii3, nvar)
      REAL_T temp(tloi:thii + 1)
c ::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ::: pcinterp:  use piecewise constant interpolation to define
c ::: values on the subregion of the fine FAB defined
c ::: by (fblo,fbhi).
c ::: 
c ::: Inputs/Outputs
c ::: fine        <=>  (modify) fab to get interpolated values
c ::: flo,fhi      =>  (const)  index limits of fine
c ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
c ::: crse         =>  (const)  fab holding coarse grid values
c ::: clo,chi      =>  (const)  index limits of src
c ::: cblo,cbhi    =>  (const)  subregion of coarse grid holding values
c ::: temp         =>  (modify) temporary space for vectorization
c ::: tlo,thi      =>  (const)  index limits of temp space
c ::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ::: local var
      integer i,j,k,ic,jc,kc,ioff,n
      integer ixproj,ii,ll
      ixproj(ii,ll) = (ii + ll*iabs(ii))/ll - iabs(ii)

      do 130 k = fblo(3), fbhi(3) 
          kc = ixproj(k,lrat)
          do 140 j = fblo(2), fbhi(2) 
              jc = ixproj(j,lrat)
              do 150 n = 1, nvar 
                  do 160 ioff = 0, lrat-1 
                      do 170 ic = cblo(1),cbhi(1) 
                          i = lrat*ic + ioff
                          temp(i) = crse(ic,jc,kc,n)
170                   continue
160               continue
                  do 180 i = fblo(1), fbhi(1) 
                      fine(i,j,k,n) = temp(i)
180               continue
150           continue
140       continue
130   continue

      return
      end


c ::: --------------------------------------------------------------
      subroutine FORT_VCINTERP (fine,floi1,floi2,floi3,fhii1,fhii2,
     $ fhii3,fbloi1,fbloi2,fbloi3, fbhii1,fbhii2,fbhii3,nvar,lratio, 
     $ crse,clo,chi,cblo,cbhi, fslo,fshi,cvloi1,cvloi2,cvloi3,cvhii1,
     $ cvhii2,cvhii3, cslope,clen,fslope,fdat,flen,voff, cvc1,fvc1,cvc2,
     $ fvc2,cvc3,fvc3)

      integer floi1,floi2,floi3
      integer fhii1,fhii2,fhii3
      integer fbloi1,fbloi2,fbloi3
      integer fbhii1,fbhii2,fbhii3
      integer cvloi1,cvloi2,cvloi3
      integer cvhii1,cvhii2,cvhii3

      integer cblo(3), cbhi(3)
      integer fslo(3), fshi(3)
      integer lratio, nvar, clen, flen, clo, chi
      REAL_T fine(floi1 :fhii1 ,floi2 :fhii2 ,floi3 :fhii3, nvar)
      REAL_T crse(clo:chi, nvar)
      REAL_T cslope(clo:chi, 3)
      REAL_T fslope(flen, 3)
      REAL_T fdat(flen)
      REAL_T voff(flen)

      REAL_T cvc1(cvloi1 :cvhii1+1)
      REAL_T cvc2(cvloi2 :cvhii2+1)
      REAL_T cvc3(cvloi3 :cvhii3+1)
      REAL_T fvc1(fbloi1 :fbhii1+1)
      REAL_T fvc2(fbloi2 :fbhii2+1)
      REAL_T fvc3(fbloi3 :fbhii3+1)

c ::: NOTE: data must be sent in so that 
c ::: cslope(1,*) and crse(1,*) are associated with
c ::: the same cell

c ::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ::: vcinterp:   conservative interpolation in volume coordinates
c ::: from coarse grid to
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
c ::: slopes are to be defined. This is
c ::: the projection of (fblo,fbhi) down
c ::: to the coarse level 
c ::: fslo,fshi    =>  (const)  fine grid index limits where
c ::: slopes are needed.  This is the
c ::: refinement of (cslo,cshi) and
c ::: contains but may not be identical
c ::: to (fblo,fbhi).
c ::: cslope       =>  (modify) temp array coarse grid slopes
c ::: clen         =>  (const)  length of coarse gtid slopes
c ::: fslope       =>  (modify) temp array for fine grid slope
c ::: flen         =>  (const)  length of fine grid slope array
c ::: fdat         =>  (const)  temp array for fine grid data
c ::: cvlo,cvhi    =>  (const)  coarse grid index limits where
c ::: volume coordinates are defined
c ::: cvc1         =>  (const)  coarse grid volume coords, x-dir
c ::: cvc2         =>  (const)  coarse grid volume coords, y_dir
c ::: cvc3         =>  (const)  coarse grid volume coords, z-dir
c ::: fvc1         =>  (const)  fine grid volume coords, x-dir
c ::: fvc2         =>  (const)  fine grid volume coords, y-dir
c ::: fvc3         =>  (const)  fine grid volume coords, z-dir
c ::: ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c ::: local var
      integer n, fn
      integer i, ic, ioff
      integer j, jc, joff
      integer k, kc, koff
      integer ist, jst, kst
      REAL_T cen, forw, back, slp
      REAL_T xoff, yoff, zoff
      REAL_T fcen, ccen
      integer ncbx, ncby, ncbz
      integer ncsx, ncsy, ncsz
      integer islo, jslo, kslo
      integer icc, istart, iend
      integer lenx, leny, lenz, maxlen

      ncbx = cbhi(1)-cblo(1)+1
      ncby = cbhi(2)-cblo(2)+1
      ncbz = cbhi(3)-cblo(3)+1
      ncsx = ncbx+2
      ncsy = ncby+2
      ncsz = ncbz+2
      ist = 1
      jst = ncsx
      kst = ncsx*ncsy
      islo = cblo(1)-1
      jslo = cblo(2)-1
      kslo = cblo(3)-1
      lenx = fbhii1-fbloi1+1
      leny = fbhii2-fbloi2+1
      lenz = fbhii3-fbloi3+1
      maxlen = max(lenx,leny,lenz)
      if (maxlen .eq. lenx) then
          do 380 i = fbloi1, fbhii1 
              fn = i-fslo(1)+1
              ic = IX_PROJ(i,lratio)
              fcen = half*(fvc1(i)+fvc1(i+1))
              ccen = half*(cvc1(ic)+cvc1(ic+1))
              voff(fn) = (fcen-ccen)/(cvc1(ic+1)-cvc1(ic))
380       continue
      else if (maxlen .eq. leny) then
          do 390 j = fbloi2, fbhii2 
              fn = j-fslo(2)+1
              jc = IX_PROJ(j,lratio)
              fcen = half*(fvc2(j)+fvc2(j+1))
              ccen = half*(cvc2(jc)+cvc2(jc+1))
              voff(fn) = (fcen-ccen)/(cvc2(jc+1)-cvc2(jc))
390       continue
      else
          do 400 k = fbloi3, fbhii3 
              fn = k-fslo(3)+1
              kc = IX_PROJ(k,lratio)
              fcen = half*(fvc3(k)+fvc3(k+1))
              ccen = half*(cvc3(kc)+cvc3(kc+1))
              voff(fn) = (fcen-ccen)/(cvc3(kc+1)-cvc3(kc))
400       continue
      end if
      do 410 n = 1, nvar 

c ::: ::::: compute slopes in x direction
          do 420 i = 1, clen 
              cen = half*(crse(i+ist,n)-crse(i-ist,n))
              forw = crse(i+ist,n)-crse(i,n)
              back = crse(i,n)-crse(i-ist,n)
              slp = sign(one,cen)*min(abs(cen),abs(forw),abs(back))
              cslope(i,1)=merge(slp,zero,forw*back>=0.0d0)
420       continue
c ::: ::::: compute slopes in y direction
          do 430 i = 1, clen 
              cen = half*(crse(i+jst,n)-crse(i-jst,n))
              forw = crse(i+jst,n)-crse(i,n)
              back = crse(i,n)-crse(i-jst,n)
              slp = sign(one,cen)*min(abs(cen),abs(forw),abs(back))
              cslope(i,2)=merge(slp,zero,forw*back>=0.0d0)
430       continue
c ::: ::::: compute slopes in z direction
          do 440 i = 1, clen 
              cen = half*(crse(i+kst,n)-crse(i-kst,n))
              forw = crse(i+kst,n)-crse(i,n)
              back = crse(i,n)-crse(i-kst,n)
              slp = sign(one,cen)*min(abs(cen),abs(forw),abs(back))
              cslope(i,3)=merge(slp,zero,forw*back>=0.0d0)
440       continue

          if (maxlen .eq. lenx) then
              do 450 kc = cblo(3), cbhi(3) 
                  do 460 jc = cblo(2), cbhi(2) 

c ::: ..,.......::::: strip out a fine grid slope vector
                      do 470 ioff = 1, lratio 
                          icc = clo + ist + jst*(jc-jslo) + kst*(kc-
     $                    kslo)
                          istart = ioff
                          iend = ioff + (ncbx-1)*lratio
                          do 480 fn = istart, iend, lratio 
                              fslope(fn,1) = cslope(icc,1)
                              fslope(fn,2) = cslope(icc,2)
                              fslope(fn,3) = cslope(icc,3)
                              fdat(fn) = crse(icc,n)
                              icc = icc + ist
480                       continue
470                   continue

                      do 490 koff = 0, lratio-1 
                          k = lratio*kc + koff
                          if (k .lt. fbloi3) then
                              goto 490
c                             --- next ---
                          end if
                          if (k .gt. fbhii3) then
                              goto 491
c                             --- break ---
                          end if
                          fcen = half*(fvc3(k)+fvc3(k+1))
                          ccen = half*(cvc3(kc)+cvc3(kc+1))
                          zoff = (fcen-ccen)/(cvc3(kc+1)-cvc3(kc))
                          do 500 joff = 0, lratio-1 
                              j = lratio*jc + joff
                              if (j .lt. fbloi2) then
                                  goto 500
c                                 --- next ---
                              end if
                              if (j .gt. fbhii2) then
                                  goto 501
c                                 --- break ---
                              end if
                              fcen = half*(fvc2(j)+fvc2(j+1))
                              ccen = half*(cvc2(jc)+cvc2(jc+1))
                              yoff = (fcen-ccen)/(cvc2(jc+1)-cvc2(jc))

                              do 510 i = fbloi1, fbhii1 
                                  fn = i-fslo(1)+1
                                  fine(i,j,k,n) = fdat(fn) + voff(fn)*
     $                            fslope(fn,1)+yoff*fslope(fn,2)+ zoff*
     $                            fslope(fn,3)
510                           continue
500                       continue
501                       continue
490                   continue
491                   continue
460               continue
450           continue
          else if (maxlen .eq. leny) then
              do 520 kc = cblo(3), cbhi(3) 
                  do 530 ic = cblo(1), cbhi(1) 

c ::: ..,.......::::: strip out a fine grid slope vector
                      do 540 joff = 1, lratio 
                          icc = clo + ist*(ic-islo) + jst + kst*(kc-
     $                    kslo)
                          istart = joff
                          iend = joff + (ncby-1)*lratio
                          do 550 fn = istart, iend, lratio 
                              fslope(fn,1) = cslope(icc,1)
                              fslope(fn,2) = cslope(icc,2)
                              fslope(fn,3) = cslope(icc,3)
                              fdat(fn) = crse(icc,n)
                              icc = icc + jst
550                       continue
540                   continue

                      do 560 koff = 0, lratio-1 
                          k = lratio*kc + koff
                          if (k .lt. fbloi3) then
                              goto 560
c                             --- next ---
                          end if
                          if (k .gt. fbhii3) then
                              goto 561
c                             --- break ---
                          end if
                          fcen = half*(fvc3(k)+fvc3(k+1))
                          ccen = half*(cvc3(kc)+cvc3(kc+1))
                          zoff = (fcen-ccen)/(cvc3(kc+1)-cvc3(kc))
                          do 570 ioff = 0, lratio-1 
                              i = lratio*ic + ioff
                              if (i .lt. fbloi1) then
                                  goto 570
c                                 --- next ---
                              end if
                              if (i .gt. fbhii1) then
                                  goto 571
c                                 --- break ---
                              end if
                              fcen = half*(fvc1(i)+fvc1(i+1))
                              ccen = half*(cvc1(ic)+cvc1(ic+1))
                              xoff = (fcen-ccen)/(cvc1(ic+1)-cvc1(ic))

                              do 580 j = fbloi2, fbhii2 
                                  fn = j-fslo(2)+1
                                  fine(i,j,k,n) = fdat(fn) + xoff*
     $                            fslope(fn,1)+voff(fn)*fslope(fn,2)+ 
     $                            zoff*fslope(fn,3)
580                           continue
570                       continue
571                       continue
560                   continue
561                   continue
530               continue
520           continue
          else
              do 590 jc = cblo(2), cbhi(2) 
                  do 600 ic = cblo(1), cbhi(1) 

c ::: ..,.......::::: strip out a fine grid slope vector
                      do 610 koff = 1, lratio 
                          icc = clo + ist*(ic-islo) + jst*(jc-jslo) + 
     $                    kst
                          istart = koff
                          iend = koff + (ncbz-1)*lratio
                          do 620 fn = istart, iend, lratio 
                              fslope(fn,1) = cslope(icc,1)
                              fslope(fn,2) = cslope(icc,2)
                              fslope(fn,3) = cslope(icc,3)
                              fdat(fn) = crse(icc,n)
                              icc = icc + kst
620                       continue
610                   continue

                      do 630 joff = 0, lratio-1 
                          j = lratio*jc + joff
                          if (j .lt. fbloi2) then
                              goto 630
c                             --- next ---
                          end if
                          if (j .gt. fbhii2) then
                              goto 631
c                             --- break ---
                          end if
                          fcen = half*(fvc2(j)+fvc2(j+1))
                          ccen = half*(cvc2(jc)+cvc2(jc+1))
                          yoff = (fcen-ccen)/(cvc2(jc+1)-cvc2(jc))
                          do 640 ioff = 0, lratio-1 
                              i = lratio*ic + ioff
                              if (i .lt. fbloi1) then
                                  goto 640
c                                 --- next ---
                              end if
                              if (i .gt. fbhii1) then
                                  goto 641
c                                 --- break ---
                              end if
                              fcen = half*(fvc1(i)+fvc1(i+1))
                              ccen = half*(cvc1(ic)+cvc1(ic+1))
                              xoff = (fcen-ccen)/(cvc1(ic+1)-cvc1(ic))

                              do 650 k = fbloi3, fbhii3 
                                  fn = k-fslo(3)+1
                                  fine(i,j,k,n) = fdat(fn) + xoff*
     $                            fslope(fn,1)+yoff*fslope(fn,2)+ 
     $                            voff(fn)*fslope(fn,3)
650                           continue
640                       continue
641                       continue
630                   continue
631                   continue
600               continue
590           continue
          end if
410   continue

      return
      end


