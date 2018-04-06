
module amrex_interp_module

  use amrex_fort_module
  use amrex_constants_module

  implicit none

  include 'AMReX_bc_types.fi'

contains

#define IX_PROJ(A,B) (A+B*iabs(A))/B-iabs(A)

! ::: --------------------------------------------------------------
! ::: nbinterp:  node based bilinear interpolation
! :::
! ::: INPUTS/OUTPUTS
! ::: fine        <=>  (modify) fine grid array
! ::: fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3   =>  (const)  index limits of fine grid
! ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
! :::
! ::: crse         =>  (const)  coarse grid data widened by 1 zone
! ::: crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3   =>  (const)  index limits of coarse grid
! :::
! ::: lratio(3)    =>  (const)  refinement ratio between levels
! ::: nvar         =>  (const)  number of components in array
! ::: num_slp      =>  (const)  number of types of slopes
! :::
! ::: TEMPORARY ARRAYS
! ::: sl           =>  num_slp 1-D slope arrays
! ::: --------------------------------------------------------------
! :::
    subroutine FORT_NBINTERP (crse, crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3, cb_l1,cb_l2,cb_l3,cb_h1,cb_h2,cb_h3, &
                              fine, fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3, fb_l1,fb_l2,fb_l3,fb_h1,fb_h2,fb_h3, &
                              lratiox, lratioy, lratioz, nvar, &
                              sl, num_slp, &
                              actual_comp, actual_state) bind(c,name='amrex_nbinterp')

      implicit none

      integer crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3
      integer cb_l1,cb_l2,cb_l3,cb_h1,cb_h2,cb_h3
      integer fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3
      integer fb_l1,fb_l2,fb_l3,fb_h1,fb_h2,fb_h3
      integer lratiox, lratioy, lratioz, nvar
      integer num_slp
      integer actual_comp,actual_state
      real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,fine_l3:fine_h3,nvar)
      real(amrex_real) crse(crse_l1:crse_h1,crse_l2:crse_h2,crse_l3:crse_h3,nvar)
      real(amrex_real) sl(cb_l1:cb_h1,num_slp)

      integer ioff, joff, koff
      integer i, j, k, ic, jc, kc, n
      integer ilo, ihi, jlo, jhi, klo, khi
      integer iratio, jratio, kratio
      integer kstrt, kstop, jstrt, jstop, istrt, istop

      real(amrex_real) fx, fy,fz
      real(amrex_real) sx, sy, sz, sxy, sxz, syz, sxyz
      real(amrex_real) RX, RY, RZ, RXY, RXZ, RYZ, RXYZ
      real(amrex_real) dx00, d0x0, d00x, dx10, dx01, d0x1, dx11

      RX   = one/lratiox
      RY   = one/lratioy
      RZ   = one/lratioz
      RXY  = RX*RY
      RXZ  = RX*RZ
      RYZ  = RY*RZ
      RXYZ = RX*RY*RZ

      do n = 1, nvar

         do kc = cb_l3, cb_h3-1

            kratio = lratioz-1
            if (kc .eq. cb_h3-1) kratio = lratioz
            kstrt = kc*lratioz
            kstop = kstrt + kratio
            klo = max(fine_l3,kstrt) - kstrt
            khi = min(fine_h3,kstop) - kstrt
            do jc = cb_l2, cb_h2-1

               jratio = lratioy-1
               if (jc .eq. cb_h2-1) jratio = lratioy
               jstrt = jc*lratioy
               jstop = jstrt + jratio
               jlo = max(fine_l2,jstrt) - jstrt
               jhi = min(fine_h2,jstop) - jstrt

               do ic = cb_l1, cb_h1-1

                  iratio = lratiox-1
                  if (ic .eq. cb_h1-1) iratio = lratiox
                  istrt = ic*lratiox
                  istop = istrt + iratio
                  ilo = max(fine_l1,istrt) - istrt
                  ihi = min(fine_h1,istop) - istrt
                  !
                  ! ::::: compute slopes
                  !
                  dx00 = crse(ic+1,jc,kc,n) - crse(ic,jc,kc,n)
                  d0x0 = crse(ic,jc+1,kc,n) - crse(ic,jc,kc,n)
                  d00x = crse(ic,jc,kc+1,n) - crse(ic,jc,kc,n)

                  dx10 = crse(ic+1,jc+1,kc,n) - crse(ic,jc+1,kc,n)
                  dx01 = crse(ic+1,jc,kc+1,n) - crse(ic,jc,kc+1,n)
                  d0x1 = crse(ic,jc+1,kc+1,n) - crse(ic,jc,kc+1,n)

                  dx11 = crse(ic+1,jc+1,kc+1,n) - crse(ic,jc+1,kc+1,n)

                  sx   = RX*dx00
                  sy   = RY*d0x0
                  sz   = RZ*d00x
                  sxy  = RXY*(dx10 - dx00)
                  sxz  = RXZ*(dx01 - dx00)
                  syz  = RYZ*(d0x1 - d0x0)
                  sxyz = RXYZ*(dx11 - dx01 - dx10 + dx00)
                  !
                  ! ::::: interpolate to fine grid
                  !
                  do koff = klo, khi
                     k = lratioz*kc + koff
                     fz = koff
                     do joff = jlo, jhi
                        j = lratioy*jc + joff
                        fy = joff
                        do ioff = ilo, ihi
                           i = lratiox*ic + ioff
                           fx = ioff
                           fine(i,j,k,n) = crse(ic,jc,kc,n) + &
                                fx*sx + fy*sy + fz*sz + &
                                fx*fy*sxy + fx*fz*sxz + fy*fz*syz + &
                                fx*fy*fz*sxyz
                        end do
                     end do
                  end do           

               end do
            end do
         end do
      end do

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
! ::: fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3   =>  (const)  index limits of fine grid
! ::: fb_l1,fb_l2,fb_l3,fb_h1,fb_h2,fb_h3     =>  (const)  subregion of fine grid to get values
! ::: 
! ::: crse         =>  (const)  coarse grid data 
! ::: crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3   =>  (const)  index limits of coarse grid
! ::: 
! ::: lratio(3)    =>  (const)  refinement ratio between levels
! ::: nvar         =>  (const)  number of components in array
! ::: 
! ::: TEMPORARY ARRAYS
! ::: slx,sly,slxy =>  1-D slope arrays
! ::: strip        =>  1-D temp array
! ::: --------------------------------------------------------------
! ::: 
    subroutine FORT_CBINTERP (crse, crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3, cb_l1,cb_l2,cb_l3,cb_h1,cb_h2,cb_h3, &
                              fine, fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3, fb_l1,fb_l2,fb_l3,fb_h1,fb_h2,fb_h3, &
                              lratiox, lratioy, lratioz, nvar, &
                              sl, num_slp, strip, strip_lo, strip_hi, &
                              actual_comp, actual_state) bind(c,name='amrex_cbinterp')

      implicit none

      integer crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3
      integer cb_l1,cb_l2,cb_l3,cb_h1,cb_h2,cb_h3
      integer fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3
      integer fb_l1,fb_l2,fb_l3,fb_h1,fb_h2,fb_h3
      integer lratiox, lratioy, lratioz, nvar
      integer num_slp
      integer strip_lo, strip_hi
      integer actual_comp,actual_state
      real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,fine_l3:fine_h3,nvar)
      real(amrex_real) crse(crse_l1:crse_h1,crse_l2:crse_h2,crse_l3:crse_h3,nvar)
      real(amrex_real) strip(strip_lo:strip_hi)
      real(amrex_real) sl(cb_l1:cb_h1, num_slp)

      call bl_abort("FORT_CBINTERP not implemented")

    end subroutine FORT_CBINTERP
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
! ::: that used in ccinterp for limslope=1; the results should
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
    subroutine FORT_LINCCINTERP (fine, fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3, fblo, fbhi, &
                                 fvcb_l1,fvcb_l2,fvcb_l3,fvcb_h1,fvcb_h2,fvcb_h3, &
                                 crse, crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3, cvcb_l1,cvcb_l2,cvcb_l3,cvcb_h1,cvcb_h2,cvcb_h3, &
                                 uc_xslope, lc_xslope, xslope_factor, &
                                 uc_yslope, lc_yslope, yslope_factor, &
                                 uc_zslope, lc_zslope, zslope_factor, &
                                 cslope_l1,cslope_l2,cslope_l3,cslope_h1,cslope_h2,cslope_h3, &
                                 cslopelo, cslopehi, &
                                 nvar, lratiox, lratioy, lratioz, &
                                 bc, lim_slope, lin_limit, &
                                 fvcx, fvcy, fvcz, cvcx, cvcy, cvcz, &
                                 voffx,voffy,voffz, alpha, cmax, cmin, &
                                 actual_comp, actual_state) bind(c,name='amrex_linccinterp')

      implicit none

      integer fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3
      integer crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3
      integer fvcb_l1,fvcb_l2,fvcb_l3,fvcb_h1,fvcb_h2,fvcb_h3
      integer cvcb_l1,cvcb_l2,cvcb_l3,cvcb_h1,cvcb_h2,cvcb_h3
      integer cslope_l1,cslope_l2,cslope_l3,cslope_h1,cslope_h2,cslope_h3
      integer fblo(3), fbhi(3)
      integer cslopelo(3), cslopehi(3)
      integer lratiox, lratioy, lratioz, nvar
      integer lim_slope, lin_limit
      integer bc(3,2,nvar)
      integer actual_comp,actual_state
      real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,fine_l3:fine_h3,nvar)
      real(amrex_real) crse(crse_l1:crse_h1,crse_l2:crse_h2,crse_l3:crse_h3, nvar)
      real(amrex_real) uc_xslope(cslope_l1:cslope_h1,cslope_l2:cslope_h2,cslope_l3:cslope_h3,nvar)
      real(amrex_real) lc_xslope(cslope_l1:cslope_h1,cslope_l2:cslope_h2,cslope_l3:cslope_h3,nvar)
      real(amrex_real) xslope_factor(cslope_l1:cslope_h1,cslope_l2:cslope_h2,cslope_l3:cslope_h3)
      real(amrex_real) uc_yslope(cslope_l1:cslope_h1,cslope_l2:cslope_h2,cslope_l3:cslope_h3,nvar)
      real(amrex_real) lc_yslope(cslope_l1:cslope_h1,cslope_l2:cslope_h2,cslope_l3:cslope_h3,nvar)
      real(amrex_real) yslope_factor(cslope_l1:cslope_h1,cslope_l2:cslope_h2,cslope_l3:cslope_h3)
      real(amrex_real) uc_zslope(cslope_l1:cslope_h1,cslope_l2:cslope_h2,cslope_l3:cslope_h3,nvar)
      real(amrex_real) lc_zslope(cslope_l1:cslope_h1,cslope_l2:cslope_h2,cslope_l3:cslope_h3,nvar)
      real(amrex_real) zslope_factor(cslope_l1:cslope_h1,cslope_l2:cslope_h2,cslope_l3:cslope_h3)
      real(amrex_real) alpha(cslope_l1:cslope_h1,cslope_l2:cslope_h2,cslope_l3:cslope_h3,nvar)
      real(amrex_real)  cmax(cslope_l1:cslope_h1,cslope_l2:cslope_h2,cslope_l3:cslope_h3,nvar)
      real(amrex_real)  cmin(cslope_l1:cslope_h1,cslope_l2:cslope_h2,cslope_l3:cslope_h3,nvar)
      real(amrex_real) fvcx(fvcb_l1:fvcb_h1)
      real(amrex_real) fvcy(fvcb_l2:fvcb_h2)
      real(amrex_real) fvcz(fvcb_l3:fvcb_h3)
      real(amrex_real) voffx(fvcb_l1:fvcb_h1)
      real(amrex_real) voffy(fvcb_l2:fvcb_h2)
      real(amrex_real) voffz(fvcb_l3:fvcb_h3)       
      real(amrex_real) cvcx(cvcb_l1:cvcb_h1)
      real(amrex_real) cvcy(cvcb_l2:cvcb_h2)
      real(amrex_real) cvcz(cvcb_l3:cvcb_h3)

#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)

      integer n 
      integer i, ic
      integer j, jc
      integer k, kc
      real(amrex_real) cen, forw, back, slp
      real(amrex_real) factorn, denom
      real(amrex_real) fxcen, cxcen, fycen, cycen, fzcen, czcen
      real(amrex_real) corr_fact,orig_corr_fact
      real(amrex_real) dummy_fine
      logical xok, yok, zok
      integer ncbx, ncby, ncbz
      integer ioff,joff,koff

      integer voff_lo(3), voff_hi(3)

      ncbx = cslopehi(1)-cslopelo(1)+1
      ncby = cslopehi(2)-cslopelo(2)+1
      ncbz = cslopehi(3)-cslopelo(3)+1
      xok = (ncbx .ge. 2)
      yok = (ncby .ge. 2)
      zok = (ncbz .ge. 2)

      voff_lo(1) = cslopelo(1) * lratiox
      voff_lo(2) = cslopelo(2) * lratioy
      voff_lo(3) = cslopelo(3) * lratioz
      voff_hi(1) = (cslopehi(1)+1) * lratiox - 1
      voff_hi(2) = (cslopehi(2)+1) * lratioy - 1
      voff_hi(3) = (cslopehi(3)+1) * lratioz - 1

      do k = voff_lo(3),voff_hi(3)
        kc = IX_PROJ(k,lratioz)
        fzcen = half*(fvcz(k)+fvcz(k+1))
        czcen = half*(cvcz(kc)+cvcz(kc+1))
        voffz(k) = (fzcen-czcen)/(cvcz(kc+1)-cvcz(kc))
      end do
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

!     Initialize alpha = 1 and define cmax and cmin as neighborhood max/mins.

          do k = cslopelo(3),cslopehi(3)
            do j = cslopelo(2),cslopehi(2)
              do i = cslopelo(1), cslopehi(1)
                alpha(i,j,k,n) = 1.d0
                cmax(i,j,k,n) = crse(i,j,k,n)
                cmin(i,j,k,n) = crse(i,j,k,n)
                do koff = -1,1
                do joff = -1,1
                do ioff = -1,1
                  cmax(i,j,k,n) = max(cmax(i,j,k,n),crse(i+ioff,j+joff,k+koff,n))
                  cmin(i,j,k,n) = min(cmin(i,j,k,n),crse(i+ioff,j+joff,k+koff,n))
                end do
                end do
                end do
              end do
            end do
          end do
      end do

!     Computed unlimited and limited slopes

      do n = 1, nvar 

          do k=cslopelo(3), cslopehi(3)
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                uc_xslope(i,j,k,n) = half*(crse(i+1,j,k,n)-crse(i-1,j,k,n))
                cen  = uc_xslope(i,j,k,n)
                forw = two*(crse(i+1,j,k,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i-1,j,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,zero,forw*back>=zero)
                lc_xslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))
             end do
            end do
          end do

          if (bclo(1,n) .eq. EXT_DIR .or. bclo(1,n).eq.HOEXTRAP) then
            i = cslopelo(1)
            if (xok) then
              do k=cslopelo(3), cslopehi(3)
                do j=cslopelo(2), cslopehi(2)
                  uc_xslope(i,j,k,n)  = -sixteen/fifteen*crse(i-1,j,k,n) &
                              + half*crse(i,j,k,n) &
                              + two3rd*crse(i+1,j,k,n) - tenth*crse(i+2,j,k,n)
                end do
              end do
            else
              do k=cslopelo(3), cslopehi(3)
                do j=cslopelo(2), cslopehi(2)
                  uc_xslope(i,j,k,n)  = fourth * ( &
                     crse(i+1,j,k,n) + five*crse(i,j,k,n) - six*crse(i-1,j,k,n) )
                end do
              end do
            endif
            do k=cslopelo(3), cslopehi(3)
              do j=cslopelo(2), cslopehi(2)
                cen  = uc_xslope(i,j,k,n)
                forw = two*(crse(i+1,j,k,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i-1,j,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,zero,forw*back>=zero)
                lc_xslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))
              end do
            end do
          end if

          if (bchi(1,n) .eq. EXT_DIR .or. bchi(1,n).eq.HOEXTRAP) then
            i = cslopehi(1)
            if (xok) then
              do k=cslopelo(3), cslopehi(3)
                do j=cslopelo(2), cslopehi(2)
                  uc_xslope(i,j,k,n) = sixteen/fifteen*crse(i+1,j,k,n) &
                            - half*crse(i,j,k,n) &
                            - two3rd*crse(i-1,j,k,n) + tenth*crse(i-2,j,k,n)
                end do
              end do
            else 
              do k=cslopelo(3), cslopehi(3)
                do j=cslopelo(2), cslopehi(2)
                  uc_xslope(i,j,k,n)  = -fourth * ( &
                     crse(i-1,j,k,n) + five*crse(i,j,k,n) - six*crse(i+1,j,k,n) )
                end do
              end do
            endif
            do k=cslopelo(3), cslopehi(3)
              do j=cslopelo(2), cslopehi(2)
                cen  = uc_xslope(i,j,k,n)
                forw = two*(crse(i+1,j,k,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i-1,j,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,zero,forw*back>=zero)
                lc_xslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))
              end do
            end do
          end if

          do k=cslopelo(3), cslopehi(3)
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                uc_yslope(i,j,k,n) = half*(crse(i,j+1,k,n)-crse(i,j-1,k,n))
                cen  = uc_yslope(i,j,k,n)
                forw = two*(crse(i,j+1,k,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i,j-1,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,zero,forw*back>=zero)
                lc_yslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))
               end do
            end do
          end do

          if (bclo(2,n) .eq. EXT_DIR .or. bclo(2,n).eq.HOEXTRAP) then
            j = cslopelo(2)
            if (yok) then
              do k=cslopelo(3), cslopehi(3)
                do i=cslopelo(1), cslopehi(1)
                  uc_yslope(i,j,k,n)  = -sixteen/fifteen*crse(i,j-1,k,n) &
                              + half*crse(i,j,k,n) &
                              + two3rd*crse(i,j+1,k,n) - tenth*crse(i,j+2,k,n)
                end do
              end do
            else
              do k=cslopelo(3), cslopehi(3)
                do i=cslopelo(1), cslopehi(1)
                  uc_yslope(i,j,k,n)  = fourth * ( &
                     crse(i,j+1,k,n) + five*crse(i,j,k,n) - six*crse(i,j-1,k,n) )
                end do
              end do
            endif
            do k=cslopelo(3), cslopehi(3)
              do i=cslopelo(1), cslopehi(1)
                cen  = uc_yslope(i,j,k,n)
                forw = two*(crse(i,j+1,k,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i,j-1,k,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,zero,forw*back>=zero)
                lc_yslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))
              end do
            end do
          end if

          if (bchi(2,n) .eq. EXT_DIR .or. bchi(2,n).eq.HOEXTRAP) then
            j = cslopehi(2)
            if (yok) then
              do k=cslopelo(3), cslopehi(3)
                do i=cslopelo(1), cslopehi(1)
                  uc_yslope(i,j,k,n) = sixteen/fifteen*crse(i,j+1,k,n) &
                            - half*crse(i,j,k,n) &
                            - two3rd*crse(i,j-1,k,n) + tenth*crse(i,j-2,k,n)
                end do
              end do
            else
              do k=cslopelo(3), cslopehi(3)
                do i=cslopelo(1), cslopehi(1)
                  uc_yslope(i,j,k,n)  = -fourth * ( &
                     crse(i,j-1,k,n) + five*crse(i,j,k,n) - six*crse(i,j+1,k,n) )
                end do
              end do
            endif
              do k=cslopelo(3), cslopehi(3)
                do i=cslopelo(1), cslopehi(1)
                  cen  = uc_yslope(i,j,k,n)
                  forw = two*(crse(i,j+1,k,n)-crse(i,j,k,n))
                  back = two*(crse(i,j,k,n)-crse(i,j-1,k,n))
                  slp  = min(abs(forw),abs(back))
                  slp  = merge(slp,zero,forw*back>=zero)
                  lc_yslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))
                end do
              end do
          end if

          do k=cslopelo(3), cslopehi(3)
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                uc_zslope(i,j,k,n) = half*(crse(i,j,k+1,n)-crse(i,j,k-1,n))
                cen  = uc_zslope(i,j,k,n)
                forw = two*(crse(i,j,k+1,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i,j,k-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,zero,forw*back>=zero)
                lc_zslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))
              end do
            end do
          end do

          if (bclo(3,n) .eq. EXT_DIR .or. bclo(3,n).eq.HOEXTRAP) then
            k = cslopelo(3)
            if (zok) then
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  uc_zslope(i,j,k,n)  = -sixteen/fifteen*crse(i,j,k-1,n) &
                              + half*crse(i,j,k,n) &
                              + two3rd*crse(i,j,k+1,n) - tenth*crse(i,j,k+2,n)
                end do
              end do
            else
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  uc_zslope(i,j,k,n)  = fourth * ( &
                     crse(i,j,k+1,n) + five*crse(i,j,k,n) - six*crse(i,j,k-1,n) )
                end do
              end do
            endif
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                cen  = uc_zslope(i,j,k,n)
                forw = two*(crse(i,j,k+1,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i,j,k-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,zero,forw*back>=zero)
                lc_zslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))
              end do
            end do
          end if

          if (bchi(3,n) .eq. EXT_DIR .or. bchi(3,n).eq.HOEXTRAP) then
            k = cslopehi(3)
            if (zok) then
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  uc_zslope(i,j,k,n) = sixteen/fifteen*crse(i,j,k+1,n) &
                            - half*crse(i,j,k,n) &
                            - two3rd*crse(i,j,k-1,n) + tenth*crse(i,j,k-2,n)
               end do
              end do
            else
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  uc_zslope(i,j,k,n)  = -fourth * ( &
                     crse(i,j,k-1,n) + five*crse(i,j,k,n) - six*crse(i,j,k+1,n) )
               end do
              end do
            endif
            do j=cslopelo(2), cslopehi(2)
              do i=cslopelo(1), cslopehi(1)
                cen  = uc_zslope(i,j,k,n)
                forw = two*(crse(i,j,k+1,n)-crse(i,j,k,n))
                back = two*(crse(i,j,k,n)-crse(i,j,k-1,n))
                slp  = min(abs(forw),abs(back))
                slp  = merge(slp,zero,forw*back>=zero)
                lc_zslope(i,j,k,n)=sign(one,cen)*min(slp,abs(cen))
             end do
            end do
          end if

       end do

       if (lim_slope.eq.0) then

!        Do the interpolation using unlimited slopes.

          do n = 1, nvar
             do k = fblo(3), fbhi(3)
               kc = IX_PROJ(k,lratioz)
               do j = fblo(2), fbhi(2)
                  jc = IX_PROJ(j,lratioy)
                  do i = fblo(1), fbhi(1)
                     ic = IX_PROJ(i,lratiox)
                     fine(i,j,k,n) = crse(ic,jc,kc,n) &
                          + voffx(i)*uc_xslope(ic,jc,kc,n) &
                          + voffy(j)*uc_yslope(ic,jc,kc,n) &
                          + voffz(k)*uc_zslope(ic,jc,kc,n)
                  end do
               end do
             end do
          end do

       else

         if (lin_limit.eq.1)then

!     compute linear limited slopes
!     Note that the limited and the unlimited slopes
!     have the same sign, and it is assumed that they do.
!
!     compute slope factors

          xslope_factor = one
          yslope_factor = one
          zslope_factor = one

          do n = 1, nvar 
            do k=cslopelo(3), cslopehi(3)
              do j=cslopelo(2), cslopehi(2)
                do i=cslopelo(1), cslopehi(1)
                  denom = uc_xslope(i,j,k,n)
                  denom = merge(denom,one,denom.ne.zero)
                  factorn = lc_xslope(i,j,k,n)/denom
                  factorn = merge(one,factorn,denom.eq.zero)
                  xslope_factor(i,j,k) = min(xslope_factor(i,j,k),factorn)

                  denom = uc_yslope(i,j,k,n)
                  denom = merge(denom,one,denom.ne.zero)
                  factorn = lc_yslope(i,j,k,n)/denom
                  factorn = merge(one,factorn,denom.eq.zero)
                  yslope_factor(i,j,k) = min(yslope_factor(i,j,k),factorn)

                  denom = uc_zslope(i,j,k,n)
                  denom = merge(denom,one,denom.ne.zero)
                  factorn = lc_zslope(i,j,k,n)/denom
                  factorn = merge(one,factorn,denom.eq.zero)
                  zslope_factor(i,j,k) = min(zslope_factor(i,j,k),factorn)
                end do
              end do
            end do
          end do

!         Compute linear limited slopes

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

         else

!         Limit slopes so as to not introduce new maxs or mins.

            do n = 1,nvar
               do k = voff_lo(3),voff_hi(3)
                  kc = IX_PROJ(k,lratioz)
                  do j = voff_lo(2),voff_hi(2)
                     jc = IX_PROJ(j,lratioy)
                     do i = voff_lo(1),voff_hi(1)
                        ic = IX_PROJ(i,lratiox)
                        orig_corr_fact = voffx(i)*lc_xslope(ic,jc,kc,n) &
                             + voffy(j)*lc_yslope(ic,jc,kc,n) &
                             + voffz(k)*lc_zslope(ic,jc,kc,n)
                        dummy_fine = crse(ic,jc,kc,n) + orig_corr_fact
                        if ((dummy_fine .gt. cmax(ic,jc,kc,n)) .and. &
                             (abs(orig_corr_fact) .gt. 1.e-10*abs(crse(ic,jc,kc,n)))) then
                           corr_fact = (cmax(ic,jc,kc,n) - crse(ic,jc,kc,n)) / orig_corr_fact
                           alpha(ic,jc,kc,n) = min(alpha(ic,jc,kc,n),corr_fact)
                        endif

                        if ((dummy_fine .lt. cmin(ic,jc,kc,n)) .and. &
                             (abs(orig_corr_fact) .gt. 1.e-10*abs(crse(ic,jc,kc,n)))) then
                           corr_fact = (cmin(ic,jc,kc,n) - crse(ic,jc,kc,n)) / orig_corr_fact
                           alpha(ic,jc,kc,n) = min(alpha(ic,jc,kc,n),corr_fact)
                        endif

#if 0
                        if (alpha(ic,jc,kc,n) .lt. 0.d0) then
                           print *,'OOPS - ALPHA SHOULD BE POSITIVE IN CCINTERP '
                           print *,'ALPHA = ',alpha(ic,jc,kc,n)
                           print *,'AT (I,J,N) = ',ic,jc,kc,n
                           print *,'ORIG_CORR_FACT = ',orig_corr_fact
                           call bl_abort(" ")
                        endif
                        if (alpha(ic,jc,kc,n) .gt. 1.d0) then
                           print *,'OOPS - ALPHA SHOULD BE <= 1.0 IN CCINTERP '
                           print *,'ALPHA = ',alpha(ic,jc,kc,n)
                           print *,'AT (I,J,N) = ',ic,jc,kc,n
                           print *,'ORIG_CORR_FACT = ',orig_corr_fact
                           call bl_abort(" ")
                        endif
#endif
                     end do
                  end do
               end do
            end do
         end if

!       Do the interpolation with limited slopes.

        do n = 1, nvar
          do k = fblo(3), fbhi(3)
            kc = IX_PROJ(k,lratioz)
            do j = fblo(2), fbhi(2)
              jc = IX_PROJ(j,lratioy)
              do i = fblo(1), fbhi(1)
                ic = IX_PROJ(i,lratiox)
                fine(i,j,k,n) = crse(ic,jc,kc,n) + alpha(ic,jc,kc,n) * &
                     ( voffx(i)*lc_xslope(ic,jc,kc,n) &
                      +voffy(j)*lc_yslope(ic,jc,kc,n) &
                      +voffz(k)*lc_zslope(ic,jc,kc,n) )
              end do
            end do
          end do
        end do

      end if

    end subroutine FORT_LINCCINTERP

    subroutine FORT_CQINTERP (fine, fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3, &
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

      call bl_abort('QUADRATIC INTERP NOT IMPLEMEMNTED IN 3-D')

    end subroutine FORT_CQINTERP

# if 0

! -----------------------------------------------------------------
! THIS IS A SCALAR VERSION OF THE ABOVE CODE
! -----------------------------------------------------------------

    subroutine FORT_CCINTERP2 (fine, fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3, &
                               fb_l1, fb_l2, fb_l3, &
                               fb_h1, fb_h2, fb_h3, &
                               nvar, lratiox, lratioy, lratioz, crse, &
                               clo, chi, cb_l1, cb_l2, cb_l3, &
                               cb_h1, cb_h2, cb_h3, &
                               fslo, fshi, cslope, clen, fslope, fdat, &
                               flen, voff, bc, limslope, &
                               actual_comp, actual_state) bind(c,name='amrex_ccinterp2')

      implicit none

      integer fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3
      integer fb_l1, fb_l2, fb_l3
      integer fb_h1, fb_h2, fb_h3
      integer cb_l1, cb_l2, cb_l3
      integer cb_h1, cb_h2, cb_h3
      integer fslo(3), fshi(3)
      integer bc(3,2,nvar)
      integer lratiox, lratioy, lratioz, nvar, clen, flen, clo, chi
      integer limslope
      integer actual_comp,actual_state
      real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,fine_l3:fine_h3,nvar)
      real(amrex_real) crse(cb_l1-1:cb_h1+1, cb_l2-1:cb_h2+1, cb_l3-1:cb_h3+1, nvar )
      real(amrex_real) cslope(cb_l1-1:cb_h1+1, cb_l2-1:cb_h2+1, cb_l3-1:cb_h3+1, 3 )
      real(amrex_real) fslope(flen, 3)
      real(amrex_real) fdat(flen)
      real(amrex_real) voff(flen)

#define bclo(i,n) bc(i,1,n)
#define bchi(i,n) bc(i,2,n)

! ::: local var
      integer n, fn
      integer i, ii, ic, ioff
      integer j, jj, jc, joff
      integer k, kk, kc, koff
      real(amrex_real) hafratx, hafraty, hafratz, volratiox, volratioy, volratioz
      real(amrex_real) cen, forw, back, slp
      real(amrex_real) xoff, yoff, zoff
      real(amrex_real) cx, cy, cz
      real(amrex_real) sgn
      integer ilo, ihi, jlo, jhi, klo, khi

!     :::::: helpful statement functions
      real(amrex_real)  slplox, slphix, slploy, slphiy, slploz, slphiz
      real(amrex_real)  c1, c2, c3, c4

      slplox(i,j,k,n) =  - c1*crse(i-1,j,k,n) &
      		         + c2*crse(i  ,j,k,n) &
                         + c3*crse(i+1,j,k,n) &
      			 - c4*crse(i+2,j,k,n)
      slphix(i,j,k,n) =    c1*crse(i+1,j,k,n) &
      		         - c2*crse(i  ,j,k,n) &
                         - c3*crse(i-1,j,k,n) &
      			 + c4*crse(i-2,j,k,n)
      slploy(i,j,k,n) =  - c1*crse(i,j-1,k,n) &
      		         + c2*crse(i,j  ,k,n) &
                         + c3*crse(i,j+1,k,n) &
      			 - c4*crse(i,j+2,k,n)
      slphiy(i,j,k,n) =    c1*crse(i,j+1,k,n) &
      		         - c2*crse(i,j  ,k,n) & 
                         - c3*crse(i,j-1,k,n) &
      			 + c4*crse(i,j-2,k,n)
      slploz(i,j,k,n) =  - c1*crse(i,j,k-1,n) &
      		         + c2*crse(i,j,k  ,n) &
                         + c3*crse(i,j,k+1,n) &
      			 - c4*crse(i,j,k+2,n)
      slphiz(i,j,k,n) =    c1*crse(i,j,k+1,n) &
      		         - c2*crse(i,j,k  ,n) &
                         - c3*crse(i,j,k-1,n) &
      			 + c4*crse(i,j,k-2,n)

      c1 = sixteen/fifteen
      c2 = half
      c3 = two3rd
      c4 = tenth

      hafratx = half*dble(lratiox-1)
      hafraty = half*dble(lratioy-1)
      hafratz = half*dble(lratioz-1)

      volratiox = one/dble(lratiox)
      volratioy = one/dble(lratioy)
      volratioz = one/dble(lratioz)

      do n = 1, nvar
      do kc = cb_l3, cb_h3
      do jc = cb_l2, cb_h2
      do ic = cb_l1, cb_h1

!        ::::: compute x slopes
	 if (limslope .ne. 0) then
            cen  = half*(crse(ic+1,jc,kc,n)-crse(ic-1,jc,kc,n))
            forw = two*(crse(ic+1,jc,kc,n)-crse(ic,jc,kc,n))
            back = two*(crse(ic,jc,kc,n) - crse(ic-1,jc,kc,n))
	    slp  = min(abs(forw),abs(back))
	    slp  = merge(slp,zero,forw*back>=zero)
	    sgn  = sign(one,cen)
            cx   = sgn*min(slp,abs(cen))
            if (ic.eq.cb_l1 .and. (bclo(1,n) .eq. EXT_DIR &
      	        .or. bclo(1,n).eq.HOEXTRAP)) then
	        cen  = slplox(ic,jc,kc,n)
                cx   = sgn*min(slp,abs(cen))
            end if
            if (ic.eq.cb_h1 .and. (bchi(1,n) .eq. EXT_DIR &
                .or. bchi(1,n).eq.HOEXTRAP)) then
                cen  = slphix(ic,jc,kc,n)
                cx   = sgn*min(slp,abs(cen))
            end if
	 else
	    cx = half*(crse(ic+1,jc,kc,n)-crse(ic-1,jc,kc,n))
            if (ic.eq.cb_l1 .and. (bclo(1,n) .eq. EXT_DIR &
     	        .or. bclo(1,n).eq.HOEXTRAP)) then
	        cx  = slplox(ic,jc,kc,n)
            end if
            if (ic.eq.cb_h1 .and. (bchi(1,n) .eq. EXT_DIR &
                .or. bchi(1,n).eq.HOEXTRAP)) then
                cx  = slphix(ic,jc,kc,n)
            end if
	 end if

!	 ::::: slopes in the Y direction
	 if (limslope .ne. 0) then
            cen  = half*(crse(ic,jc+1,kc,n)-crse(ic,jc-1,kc,n))
            forw = two*(crse(ic,jc+1,kc,n)-crse(ic,jc,kc,n))
            back = two*(crse(ic,jc,kc,n) - crse(ic,jc-1,kc,n))
	    slp  = min(abs(forw),abs(back))
	    slp  = merge(slp,zero,forw*back>=zero)
	    sgn  = sign(one,cen)
            cy   = sgn*min(slp,abs(cen))
            if (jc.eq.cb_l2 .and. (bclo(2,n) .eq. EXT_DIR &
      	        .or. bclo(2,n).eq.HOEXTRAP)) then
	        cen  = slploy(ic,jc,kc,n)
                cy   = sgn*min(slp,abs(cen))
            end if
            if (jc.eq.cb_h2 .and. (bchi(2,n) .eq. EXT_DIR &
                .or. bchi(2,n).eq.HOEXTRAP)) then
                cen  = slphiy(ic,jc,kc,n)
                cy   = sgn*min(slp,abs(cen))
            end if
	 else
	    cy = half*(crse(ic,jc+1,kc,n)-crse(ic,jc-1,kc,n))
            if (jc.eq.cb_l2 .and. (bclo(2,n) .eq. EXT_DIR &
      	        .or. bclo(2,n).eq.HOEXTRAP)) then
	        cy   = slploy(ic,jc,kc,n)
            end if
            if (ic.eq.cb_h2 .and. (bchi(2,n) .eq. EXT_DIR &
                .or. bchi(2,n).eq.HOEXTRAP)) then
                cy   = slphiy(ic,jc,kc,n)
            end if
	 end if

!	 ::::: slopes in the Z direction
	 if (limslope .ne. 0) then
            cen  = half*(crse(ic,jc,kc+1,n)-crse(ic,jc,kc-1,n))
            forw = two*(crse(ic,jc,kc+1,n)-crse(ic,jc,kc,n))
            back = two*(crse(ic,jc,kc,n) - crse(ic,jc,kc-1,n))
	    slp  = min(abs(forw),abs(back))
	    slp  = merge(slp,zero,forw*back>=zero)
	    sgn  = sign(one,cen)
            cz   = sgn*min(slp,abs(cen))
            if (kc.eq.cb_l3 .and. (bclo(3,n) .eq. EXT_DIR &
      	        .or. bclo(3,n).eq.HOEXTRAP)) then
	        cen  = slploz(ic,jc,kc,n)
                cz   = sgn*min(slp,abs(cen))
            end if
            if (kc.eq.cb_h3 .and. (bchi(3,n) .eq. EXT_DIR &
                .or. bchi(3,n).eq.HOEXTRAP)) then
                cen  = slphiz(ic,jc,kc,n)
                cz   = sgn*min(slp,abs(cen))
            end if
	 else
	    cz = half*(crse(ic,jc,kc+1,n)-crse(ic,jc,kc-1,n))
            if (kc.eq.cb_l3 .and. (bclo(3,n) .eq. EXT_DIR &
      	        .or. bclo(3,n).eq.HOEXTRAP)) then
	        cz   = slploz(ic,jc,kc,n)
            end if
            if (kc.eq.cb_h3 .and. (bchi(3,n) .eq. EXT_DIR &
                .or. bchi(3,n).eq.HOEXTRAP)) then
                cz   = slphiz(ic,jc,kc,n)
            end if
	 end if

!	 ::::: now interpolate to fine grid
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
	    zoff = dble(koff)-hafratz
	    do joff = jlo, jhi
	       j = lratioy*jc + joff
	       yoff = dble(joff)-hafraty
	       do ioff = ilo, ihi
	          i = lratiox*ic + ioff
		  xoff = dble(ioff)-hafratx
		  fine(i,j,k,n) = crse(ic,jc,kc,n) + &
      				( volratiox*xoff*cx + volratioy*yoff*cy &
                                + volratioz*zoff*cz )
	       end do
	    end do
	 end do
 
      end do
      end do
      end do
      end do

    end subroutine FORT_CCINTERP2
#endif

! ::: 
! ::: --------------------------------------------------------------
! ::: pcinterp:  cell centered piecewise constant interpolation
! ::: 
! ::: Inputs/Outputs
! ::: fine        <=>  (modify) fine grid array
! ::: fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3   =>  (const)  index limits of fine grid
! ::: fblo,fbhi    =>  (const)  subregion of fine grid to get values
! ::: 
! ::: crse         =>  (const)  coarse grid data 
! ::: crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3   =>  (const)  index limits of coarse grid
! ::: cblo,cbhi    =>  (const) coarse grid region containing fblo,fbhi
! ::: 
! ::: longdir      =>  (const)  which index direction is longest (1 or 2)
! ::: lratio(3)    =>  (const)  refinement ratio between levels
! ::: nvar         =>  (const)  number of components in array
! ::: 
! ::: TEMPORARY ARRAYS
! ::: ftmp         =>  1-D temp array
! ::: --------------------------------------------------------------
! ::: 
    subroutine FORT_PCINTERP (crse,crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3,cblo,cbhi, &
                              fine,fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3,fblo,fbhi, &
                              longdir,lratiox,lratioy,lratioz,nvar, &
                              ftmp, ftmp_lo, ftmp_hi, &
                              actual_comp, actual_state) bind(c,name='amrex_pcinterp')

      implicit none

      integer crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3
      integer cblo(3), cbhi(3)
      integer fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3
      integer fblo(3), fbhi(3)
      integer nvar, lratiox, lratioy, lratioz, longdir
      integer ftmp_lo, ftmp_hi
      integer actual_comp,actual_state
      real(amrex_real)  crse(crse_l1:crse_h1,crse_l2:crse_h2,crse_l3:crse_h3, nvar)
      real(amrex_real)  fine(fine_l1:fine_h1,fine_l2:fine_h2,fine_l3:fine_h3, nvar)
      real(amrex_real)  ftmp(ftmp_lo:ftmp_hi)

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

!	       ::::: fill strip in i direction
	       do ioff = 0, lratiox-1
	          do ic = cblo(1), cbhi(1)
	             i = lratiox*ic + ioff
	             ftmp(i) = crse(ic,jc,kc,n)
                  end do
	       end do

!	       ::::: stuff into fine array
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

!	       ::::: fill strip in j direction
	       do joff = 0, lratioy-1
	          do jc = cblo(2), cbhi(2)
	             j = lratioy*jc + joff
	             ftmp(j) = crse(ic,jc,kc,n)
                  end do
	       end do

!	       ::::: stuff into fine array
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

!	       ::::: fill strip in k direction
	       do koff = 0, lratioz-1
	          do kc = cblo(3), cbhi(3)
	             k = lratioz*kc + koff
	             ftmp(k) = crse(ic,jc,kc,n)
                  end do
	       end do

!	       ::::: stuff into fine array
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
    subroutine FORT_PROTECT_INTERP (fine, fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3, fblo, fbhi, &
                                    crse, crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3, cblo, cbhi, &
                                    fine_state, state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                                    nvar, lratiox,lratioy,lratioz, bc) bind(c,name='amrex_protect_interp')

      implicit none

      integer fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3
      integer crse_l1,crse_l2,crse_l3,crse_h1,crse_h2,crse_h3
      integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
      integer fblo(3), fbhi(3)
      integer cblo(3), cbhi(3)
      integer lratiox, lratioy, lratioz, nvar
      integer bc(3,2,nvar)
      real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,fine_l3:fine_h3,nvar)
      real(amrex_real) crse(crse_l1:crse_h1,crse_l2:crse_h2,crse_l3:crse_h3, nvar)
      real(amrex_real) fine_state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3, nvar)

      integer rMAX
      parameter (rMAX = 16)
      real(amrex_real) alpha, sumN, sumP, crseTot, negVal, posVal
      real(amrex_real) sum_fine_new,sum_fine_old
      real(amrex_real) orig_fine(0:rMAX-1,0:rMAX-1,0:rMAX-1)
      integer redo_me
      integer ilo,ihi,jlo,jhi,klo,khi
      integer i,j,k,ic,jc,kc,n
      integer numFineCells
      integer icase

      if (MAX(lratiox,lratioy,lratioz).gt.rMAX) then
         print *,'rMAX in INTERP_3D::FORT_PROTECT_INTERP must be >= ', &
              MAX(lratiox,lratioy,lratioz)
         call bl_abort(" ")
      endif

      do kc = cblo(3), cbhi(3)
      do jc = cblo(2), cbhi(2)
      do ic = cblo(1), cbhi(1)

         ilo = max(lratiox*ic            ,fine_l1)
         ihi = min(lratiox*ic+(lratiox-1),fine_h1)
         jlo = max(lratioy*jc            ,fine_l2)
         jhi = min(lratioy*jc+(lratioy-1),fine_h2)
         klo = max(lratioz*kc            ,fine_l3)
         khi = min(lratioz*kc+(lratioz-1),fine_h3)

         do n = 2, nvar-1

            redo_me = 0
            do k = klo,khi
            do j = jlo,jhi
            do i = ilo,ihi
               if ((fine_state(i,j,k,n)+fine(i,j,k,n)) .lt. 0.d0) redo_me = 1
            enddo
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
!             crseTot = sum of all interpolated values of the correction,
!                       which is equivalent to the coarse correction * ratio**3
!             SumN = sum of all negative values of fine_state
!             SumP = sum of all positive values of fine_state
!
! ****************************************************************************************

            if (redo_me .eq. 1) then

               icase = 0
               sum_fine_old = 0.d0
               do k = klo,khi
               do j = jlo,jhi
               do i = ilo,ihi
                  sum_fine_old = sum_fine_old + fine(i,j,k,n)
                  orig_fine(i-ilo,j-jlo,k-klo) = fine(i,j,k,n)
               enddo
               enddo
               enddo

               crseTot = sum_fine_old
               numFineCells = (ihi-ilo+1) * (jhi-jlo+1) * (khi-klo+1)

               sumN = zero
               sumP = zero
               do k = klo,khi
               do j = jlo,jhi
               do i = ilo,ihi
                  if (fine_state(i,j,k,n) .le. 0.d0) then
                    sumN = SumN + fine_state(i,j,k,n)
                  else
                    sumP = sumP + fine_state(i,j,k,n)
                  endif
               enddo
               enddo
               enddo

               if (crseTot .gt. 0.d0 .and. crseTot .ge. abs(sumN)) then
!              Here we want to fill in the negative values first, then add
!                the remaining positive proportionally.

                   icase = 1
                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,k,n) .le. 0.d0) then
                        fine(i,j,k,n) = -fine_state(i,j,k,n)
                      endif
                   enddo
                   enddo
                   enddo

                   if (sumP .gt. 0.d0) then

                    alpha = (crseTot - abs(sumN)) / sumP

                    do k = klo,khi
                    do j = jlo,jhi
                    do i = ilo,ihi
                       if (fine_state(i,j,k,n) .ge. 0.d0) then
                         fine(i,j,k,n) = alpha * fine_state(i,j,k,n)
                       endif
                    enddo
                    enddo
                    enddo

                  else

                    posVal = (crseTot - abs(sumN)) / float(numFineCells)

                    do k = klo,khi
                    do j = jlo,jhi
                    do i = ilo,ihi
                       fine(i,j,k,n) = fine(i,j,k,n) + posVal
                    enddo
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

                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,k,n) .lt. 0.d0) then
                        fine(i,j,k,n) = alpha * abs(fine_state(i,j,k,n))
                      else 
                        fine(i,j,k,n) = 0.d0
                      endif
                   enddo
                   enddo
                   enddo

               endif

               if (crseTot .lt. 0.d0 .and. abs(crseTot) .gt. sumP) then
!              Here we don't have enough positive states to absorb all the
!                negative correction, so we want to end up with all the fine
!                cells having the same negative value.

                   icase = 3
                   negVal = (sumP + sumN + crseTot)/float(numFineCells)

                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      fine(i,j,k,n) = negVal - fine_state(i,j,k,n)
                   enddo
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

                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,k,n) .lt. 0.d0) then
                        fine(i,j,k,n) = -fine_state(i,j,k,n)
                      else
                        fine(i,j,k,n) = alpha * fine_state(i,j,k,n)
                      endif  
                   enddo
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

                   do k = klo,khi
                   do j = jlo,jhi
                   do i = ilo,ihi
                      if (fine_state(i,j,k,n) .gt. 0.d0) then
                        fine(i,j,k,n) = -fine_state(i,j,k,n)
                      else 
                        fine(i,j,k,n) = alpha * fine_state(i,j,k,n)
                      endif
                   enddo
                   enddo
                   enddo

               endif

               sum_fine_new = 0.d0
               do k = klo,khi
               do j = jlo,jhi
               do i = ilo,ihi
                  sum_fine_new = sum_fine_new + fine(i,j,k,n)
               enddo
               enddo
               enddo

               if (abs(sum_fine_new - sum_fine_old) .gt. 1.e-8) then
                  print *,' '
                  print *, &
                       'PROTECT_INTERP: BLEW CONSERVATION with ICASE = ' &
                       ,icase
                  print *,'AT COARSE CELL ',ic,jc,kc,' AND COMPONENT ',n
                  print *,'NEW SUM / OLD SUM ', &
                       sum_fine_new, sum_fine_old
                  print *,'CRSETOT ',crseTot
                  print *,'SUMP SUMN ',sumP,sumN
                  do k = klo,khi
                  do j = jlo,jhi
                  do i = ilo,ihi
                     print *,'FINE OLD NEW ', &
                          i,j,k,orig_fine(i-ilo,j-jlo,k-klo), &
                          fine(i,j,k,n), fine_state(i,j,k,n)
                  enddo
                  enddo
                  enddo
               endif

!              do k = klo,khi
!              do j = jlo,jhi
!              do i = ilo,ihi
!                 if ((fine_state(i,j,k,n) + fine(i,j,k,n)) .lt. 0.d0) then
!                    print *,'STILL NEGATIVE AT ',i,j,k,n
!                    print *,'AT COARSE CELL ',ic,jc,kc
!                    print *,'FINE STATE ',fine_state(i,j,k,n)
!                    print *,'FINE CORRECTION ',fine(i,j,k,n)
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
         do k = klo,khi
         do j = jlo,jhi
         do i = ilo,ihi
            fine(i,j,k,1) = 0.d0
            do n = 2,nvar-1
               fine(i,j,k,1) = fine(i,j,k,1) + fine(i,j,k,n)
            enddo
         enddo
         enddo
         enddo

!     End of coarse index loops
      enddo
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
     subroutine FORT_QUARTINTERP (fine, fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3, &
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
!      todo
          write(6,*) 'FORT_QUARTINTERP: refinement ratio = 4 TODO'
          stop
       else
          write(6,*) 'FORT_QUARTINTERP: unsupported refinement ratio'
          stop
       endif

     end subroutine FORT_QUARTINTERP

end module amrex_interp_module
