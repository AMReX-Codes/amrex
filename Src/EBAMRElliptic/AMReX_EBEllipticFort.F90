#include "AMReX_CONSTANTS.H"

module ebfnd_average_module

  !     since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
  !     for e.g., #if (BL_SPACEDIM == 1) statements.

  implicit none

  public

contains

  subroutine ebfnd_average( &
       coar, coar_lo, coar_hi,  coar_nco, &
       fine, fine_lo, fine_hi,  fine_nco, &
       coarboxlo,coarboxhi, &
       refboxlo,  refboxhi, refrat, &
       isrc, idst, ncomp) &
       bind(C, name="ebfnd_average")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: iif,jjf,kkf, coar_nco, fine_nco
    integer      :: iic,jjc,kkc, refrat, ncomp, isrc, idst
    integer      :: ib, jb, kb, ivar, ivarc, ivarf, idir
    integer      :: coar_lo(0:2),coar_hi(0:2)
    integer      :: fine_lo(0:2),fine_hi(0:2)
    integer      :: coarboxlo(0:2), coarboxhi(0:2)
    integer      :: refboxlo(0:2),   refboxhi(0:2)
    real(c_real) :: numfinepercoar, fineval, coarval
    real(c_real) :: fine(fine_lo(0):fine_hi(0),fine_lo(1):fine_hi(1),fine_lo(2):fine_hi(2), 0:fine_nco-1)
    real(c_real) :: coar(coar_lo(0):coar_hi(0),coar_lo(1):coar_hi(1),coar_lo(2):coar_hi(2), 0:coar_nco-1)

    numfinepercoar = one
    do idir = 1, BL_SPACEDIM
       numfinepercoar = numfinepercoar*refrat
    enddo

    do ivar = 0, ncomp-1
       ivarf = isrc + ivar
       ivarc = idst + ivar
       do kkc = coarboxlo(2), coarboxhi(2)
          do jjc = coarboxlo(1), coarboxhi(1)
             do iic = coarboxlo(0), coarboxhi(0)

                coar(iic,jjc,kkc, ivarc) = zero

                do kb = refboxlo(2), refboxhi(2)
                   do jb = refboxlo(1), refboxhi(1)
                      do ib = refboxlo(0), refboxhi(0)

                         iif = refrat*iic + ib
                         jjf = refrat*jjc + jb
                         kkf = refrat*kkc + kb
                         coar(iic,jjc,kkc, ivarc) = coar(iic, jjc, kkc,ivarc) + fine(iif, jjf, kkf, ivarf)
                         fineval = fine(iif, jjf, kkf, ivarf)
                         coarval = coar(iic,jjc,kkc, ivarc)

                      enddo
                   enddo
                enddo

                coar(iic,jjc,kkc, ivarc) = coar(iic, jjc, kkc,ivarc)/numfinepercoar

             enddo
          enddo
       enddo
    enddo


  end subroutine ebfnd_average

end module ebfnd_average_module



module ebfnd_pwcinterp_module

  !     since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
  !     for e.g., #if (BL_SPACEDIM == 1) statements.

  implicit none

  public

contains

  subroutine ebfnd_pwcinterp( &
       fine, fine_lo, fine_hi, fine_nco,  &
       coar, coar_lo, coar_hi, coar_nco,  &
       fineboxlo,fineboxhi, &
       refrat, isrc, idst, ncomp) &
       bind(C, name="ebfnd_pwcinterp")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: iif,jjf,kkf,  ncomp, ivar, ivarf, ivarc
    integer      :: iic,jjc,kkc, refrat, coar_nco, fine_nco, isrc, idst
    integer      :: coar_lo(0:2),coar_hi(0:2)
    integer      :: fine_lo(0:2),fine_hi(0:2)
    integer      :: fineboxlo(0:2), fineboxhi(0:2)

    real(c_real) :: fine(fine_lo(0):fine_hi(0),fine_lo(1):fine_hi(1),fine_lo(2):fine_hi(2), 0:fine_nco-1)
    real(c_real) :: coar(coar_lo(0):coar_hi(0),coar_lo(1):coar_hi(1),coar_lo(2):coar_hi(2), 0:coar_nco-1)

    do ivar = 0, ncomp-1
       ivarc = isrc + ivar
       ivarf = idst + ivar

       do kkf = fineboxlo(2), fineboxhi(2)
          do jjf = fineboxlo(1), fineboxhi(1)
             do iif = fineboxlo(0), fineboxhi(0)

                iic = iif/refrat
                jjc = jjf/refrat
                kkc = kkf/refrat

                fine(iif, jjf, kkf, ivarf) = coar(iic,jjc,kkc, ivarc)

             enddo
          enddo
       enddo
    enddo


  end subroutine ebfnd_pwcinterp

end module ebfnd_pwcinterp_module



module ebfnd_pwlinterp_nobound_module

  !     since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
  !     for e.g., #if (BL_SPACEDIM == 1) statements.

  implicit none

  public

contains

  ! this is piecewise linear interp all in one pass - away from boundaries
  subroutine ebfnd_pwlinterp_nobound( &
       fine, fine_lo, fine_hi, fine_nco,  &
       coar, coar_lo, coar_hi, coar_nco,  &
       fineboxlo,fineboxhi, & 
       refrat, isrc, idst, ncomp) &
       bind(C, name="ebfnd_pwlinterp_nobound")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: iif,jjf,kkf,  ncomp, ivar, ivarf, ivarc
    integer      :: iic,jjc,kkc, refrat, coar_nco, fine_nco, isrc, idst
    integer      :: coar_lo(0:2),coar_hi(0:2)
    integer      :: fine_lo(0:2),fine_hi(0:2)
    integer      :: fineboxlo(0:2), fineboxhi(0:2)

    real(c_real) :: fine(fine_lo(0):fine_hi(0),fine_lo(1):fine_hi(1),fine_lo(2):fine_hi(2), 0:fine_nco-1)
    real(c_real) :: coar(coar_lo(0):coar_hi(0),coar_lo(1):coar_hi(1),coar_lo(2):coar_hi(2), 0:coar_nco-1)
    real(c_real) :: xdist, ydist, zdist, xslope, yslope, zslope, dxf, dxc
    real(c_real) :: xcoar, ycoar, zcoar, xfine, yfine, zfine

    dxf = one
    dxc = refrat

    do ivar = 0, ncomp-1
       ivarc = isrc + ivar
       ivarf = idst + ivar

       ! first do the easy bit--not near any boundaries
       do kkf = fineboxlo(2), fineboxhi(2)
          do jjf = fineboxlo(1), fineboxhi(1)
             do iif = fineboxlo(0), fineboxhi(0)

                iic = iif/refrat
                jjc = jjf/refrat
                kkc = kkf/refrat

                xcoar = dxc*(iic + half)
                ycoar = dxc*(jjc + half)
                zcoar = zero

                xfine = dxf*(iif + half)
                yfine = dxf*(jjf + half)
                zfine = zero

                xslope = half*(coar(iic+1,jjc  ,kkc  , ivarc)-  coar(iic-1,jjc  ,kkc  , ivarc))/dxc 
                yslope = half*(coar(iic  ,jjc+1,kkc  , ivarc)-  coar(iic  ,jjc-1,kkc  , ivarc))/dxc
                zslope = zero

                xdist = xfine - xcoar
                ydist = yfine - ycoar
                zdist = zero

#if BL_SPACEDIM == 3
                zcoar = dxc*(kkc + half)
                zfine = dxf*(kkf + half)
                zdist = zfine - zcoar
                zslope = half*(coar(iic  ,jjc  ,kkc+1, ivarc)-  coar(iic  ,jjc  ,kkc-1, ivarc))/dxc
#endif
                fine(iif, jjf, kkf, ivarf) = coar(iic,jjc,kkc, ivarc) &
                     + xdist * xslope  &
                     + ydist * yslope  &
                     + zdist * zslope 

             enddo
          enddo
       enddo


    enddo


  end subroutine ebfnd_pwlinterp_nobound

end module ebfnd_pwlinterp_nobound_module



module ebfnd_pwlincr_at_bound_module

  !     since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
  !     for e.g., #if (BL_SPACEDIM == 1) statements.

  implicit none

  public

contains

  ! this is piecewise linear incrementing by the slope*dist in one direction
  subroutine ebfnd_pwl_incr_at_bound( &
       fine, fine_lo, fine_hi, fine_nco,  &
       coar, coar_lo, coar_hi, coar_nco,  &
       has_lo, has_hi, idir,  &
       coar_lo_box_lo,  coar_lo_box_hi, &
       coar_hi_box_lo,  coar_hi_box_hi, &
       coar_cen_box_lo, coar_cen_box_hi, &
       refboxlo, refboxhi, &
       refrat, isrc, idst, ncomp) &
         bind(C, name="ebfnd_pwl_incr_at_bound")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none
    integer      :: has_lo,  has_hi, idir, ioff, joff, koff
    integer      :: coar_lo_box_lo(0:2),  coar_lo_box_hi(0:2)
    integer      :: coar_hi_box_lo(0:2),  coar_hi_box_hi(0:2)
    integer      :: coar_cen_box_lo(0:2),  coar_cen_box_hi(0:2)
    integer      :: refboxlo(0:2),   refboxhi(0:2), ib, jb, kb

    integer      :: iif,jjf,kkf,  ncomp, ivar, ivarf, ivarc
    integer      :: iic,jjc,kkc, refrat, coar_nco, fine_nco, isrc, idst
    integer      :: coar_lo(0:2),coar_hi(0:2)
    integer      :: fine_lo(0:2),fine_hi(0:2)


    real(c_real) :: fine(fine_lo(0):fine_hi(0),fine_lo(1):fine_hi(1),fine_lo(2):fine_hi(2), 0:fine_nco-1)
    real(c_real) :: coar(coar_lo(0):coar_hi(0),coar_lo(1):coar_hi(1),coar_lo(2):coar_hi(2), 0:coar_nco-1)
    real(c_real) :: xdist, xslope,  dxf, dxc,  indf, indc
    real(c_real) :: xcoar, xfine, finevalold, finevalnew, coarhi, coarlo

    dxf = one
    dxc = refrat
    ioff = 0
    joff = 0
    koff = 0
    if(idir .eq. 0) then
       ioff = 1
    else if (idir .eq. 1) then
       joff = 1
    else if (idir .eq. 2) then
       koff = 1
    endif



    do ivar = 0, ncomp-1
       ivarc = isrc + ivar
       ivarf = idst + ivar

       !first the center bit
       do kkc = coar_cen_box_lo(2),  coar_cen_box_hi(2)
          do jjc = coar_cen_box_lo(1),  coar_cen_box_hi(1)
             do iic = coar_cen_box_lo(0),  coar_cen_box_hi(0)

                do kb = refboxlo(2), refboxhi(2)
                   do jb = refboxlo(1), refboxhi(1)
                      do ib = refboxlo(0), refboxhi(0)

                         iif = refrat*iic + ib
                         jjf = refrat*jjc + jb
                         kkf = refrat*kkc + kb
                         if(idir .eq. 0) then
                            indf = iif
                            indc = iic
                         else if (idir .eq. 1) then
                            indf = jjf
                            indc = jjc

                         else if (idir .eq. 2) then
                            indf = kkf
                            indc = kkc
                         endif

                         xcoar = dxc*(indc + half)
                         xfine = dxf*(indf + half)
                         xdist = xfine - xcoar

                         coarhi = coar(iic+ioff,jjc+joff,kkc+koff, ivarc)
                         coarlo = coar(iic-ioff,jjc-joff,kkc-koff, ivarc)
                         xslope = (coarhi - coarlo)/(two*dxc)
                         finevalold = fine(iif, jjf, kkf, ivarf)


                         finevalnew = finevalold  + xdist * xslope
                         fine(iif, jjf, kkf, ivarf) = finevalnew

                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo

       !now use one-sided diff for lo and hi
       if(has_hi .eq. 1) then

          do kkc = coar_hi_box_lo(2),  coar_hi_box_hi(2)
             do jjc = coar_hi_box_lo(1),  coar_hi_box_hi(1)
                do iic = coar_hi_box_lo(0),  coar_hi_box_hi(0)

                   do kb = refboxlo(2), refboxhi(2)
                      do jb = refboxlo(1), refboxhi(1)
                         do ib = refboxlo(0), refboxhi(0)

                            iif = refrat*iic + ib
                            jjf = refrat*jjc + jb
                            kkf = refrat*kkc + kb
                            if(idir .eq. 0) then
                               indf = iif
                               indc = iic
                            else if (idir .eq. 1) then
                               indf = jjf
                               indc = jjc

                            else if (idir .eq. 2) then
                               indf = kkf
                               indc = kkc
                            endif

                            xcoar = dxc*(indc + half)
                            xfine = dxf*(indf + half)
                            xdist = xfine - xcoar

                            coarhi = coar(iic     ,jjc     ,kkc     , ivarc)
                            coarlo = coar(iic-ioff,jjc-joff,kkc-koff, ivarc)

                            xslope = (coarhi - coarlo)/dxc
                            fine(iif, jjf, kkf, ivarf) = fine(iif, jjf, kkf, ivarf)  + xdist * xslope

                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       endif

       if(has_lo .eq. 1) then

          do kkc = coar_lo_box_lo(2),  coar_lo_box_hi(2)
             do jjc = coar_lo_box_lo(1),  coar_lo_box_hi(1)
                do iic = coar_lo_box_lo(0),  coar_lo_box_hi(0)

                   do kb = refboxlo(2), refboxhi(2)
                      do jb = refboxlo(1), refboxhi(1)
                         do ib = refboxlo(0), refboxhi(0)

                            iif = refrat*iic + ib
                            jjf = refrat*jjc + jb
                            kkf = refrat*kkc + kb
                            if(idir .eq. 0) then
                               indf = iif
                               indc = iic
                            else if (idir .eq. 1) then
                               indf = jjf
                               indc = jjc

                            else if (idir .eq. 2) then
                               indf = kkf
                               indc = kkc
                            endif

                            xcoar = dxc*(indc + half)
                            xfine = dxf*(indf + half)
                            xdist = xfine - xcoar
                            coarhi = coar(iic+ioff,jjc+joff,kkc+koff, ivarc)
                            coarlo = coar(iic     ,jjc     ,kkc     , ivarc)

                            xslope = (coarhi - coarlo)/dxc

                            fine(iif, jjf, kkf, ivarf) = fine(iif, jjf, kkf, ivarf)  + xdist * xslope

                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       endif

    enddo


  end subroutine ebfnd_pwl_incr_at_bound

end module ebfnd_pwlincr_at_bound_module


module ebfnd_pwqinterp_nobound_module

  !     since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
  !     for e.g., #if (BL_SPACEDIM == 1) statements.

  implicit none

  public

contains

  ! this is piecewise linear interp all in one pass - away from boundaries
  subroutine ebfnd_pwqinterp_nobound( &
       fine, fine_lo, fine_hi, fine_nco,  &
       coar, coar_lo, coar_hi, coar_nco,  &
       fineboxlo,fineboxhi, & 
       refrat, isrc, idst, ncomp) &
       bind(C, name="ebfnd_pwqinterp_nobound")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: iif,jjf,kkf,  ncomp, ivar, ivarf, ivarc
    integer      :: iic,jjc,kkc, refrat, coar_nco, fine_nco, isrc, idst
    integer      :: coar_lo(0:2),coar_hi(0:2)
    integer      :: fine_lo(0:2),fine_hi(0:2)
    integer      :: fineboxlo(0:2), fineboxhi(0:2)

    real(c_real) :: fine(fine_lo(0):fine_hi(0),fine_lo(1):fine_hi(1),fine_lo(2):fine_hi(2), 0:fine_nco-1)
    real(c_real) :: coar(coar_lo(0):coar_hi(0),coar_lo(1):coar_hi(1),coar_lo(2):coar_hi(2), 0:coar_nco-1)
    real(c_real) :: xdist, ydist, zdist, xslope, yslope, zslope, dxf, dxc
    real(c_real) :: xcoar, ycoar, zcoar, xfine, yfine, zfine, coarval, finevalnew
    real(c_real) :: dxx, dyy, dzz, dxy, dxz, dyz, coarhix, coarhiy, coarhiz,coarlox, coarloy, coarloz

    dxf = one
    dxc = refrat

    do ivar = 0, ncomp-1
       ivarc = isrc + ivar
       ivarf = idst + ivar

       ! first do the easy bit--not near any boundaries
       do kkf = fineboxlo(2), fineboxhi(2)
          do jjf = fineboxlo(1), fineboxhi(1)
             do iif = fineboxlo(0), fineboxhi(0)

                iic = iif/refrat
                jjc = jjf/refrat
                kkc = kkf/refrat

                xcoar = dxc*(iic + half)
                ycoar = dxc*(jjc + half)
                zcoar = zero

                xfine = dxf*(iif + half)
                yfine = dxf*(jjf + half)
                zfine = zero
                
                coarhix = coar(iic+1,jjc  ,kkc  , ivarc)
                coarlox = coar(iic-1,jjc  ,kkc  , ivarc)
                coarhiy = coar(iic  ,jjc+1,kkc  , ivarc)
                coarloy = coar(iic  ,jjc-1,kkc  , ivarc)
                coarhiz = 0
                coarloz = 0
                xslope = (coarhix - coarlox)/(two*dxc) 
                yslope = (coarhiy - coarloy)/(two*dxc) 
                zslope = zero

                dxx = (  coar(iic+1,jjc  ,kkc  ,ivarc) &
                     +   coar(iic-1,jjc  ,kkc  ,ivarc)-  &
                     two*coar(iic  ,jjc  ,kkc  ,ivarc))/(dxc*dxc) 
                dyy = (  coar(iic  ,jjc+1,kkc  ,ivarc) &
                     +   coar(iic  ,jjc-1,kkc  ,ivarc)-  &
                     two*coar(iic  ,jjc  ,kkc  ,ivarc))/(dxc*dxc) 

                dxy = (coar(iic+1,jjc+1,kkc ,ivarc) &
                     + coar(iic-1,jjc-1,kkc, ivarc) &
                     - coar(iic+1,jjc-1,kkc ,ivarc) &
                     - coar(iic-1,jjc+1,kkc, ivarc))/(four*dxc*dxc) 
                dxz = zero
                dyz = zero
                dzz = zero
                xdist = xfine - xcoar
                ydist = yfine - ycoar
                zdist = zero

#if BL_SPACEDIM == 3
                coarhiz = coar(iic  ,jjc  ,kkc+1, ivarc)
                coarloz = coar(iic  ,jjc  ,kkc-1, ivarc)
                zslope = (coarhiz - coarloz)/(two*dxc) 
                dxz = (coar(iic+1,jjc  ,kkc+1,ivarc) &
                     + coar(iic-1,jjc  ,kkc-1,ivarc) &
                     - coar(iic+1,jjc  ,kkc-1,ivarc) &
                     - coar(iic-1,jjc  ,kkc+1,ivarc))/(four*dxc*dxc) 
                dyz = (coar(iic  ,jjc+1,kkc+1,ivarc)&
                     + coar(iic  ,jjc-1,kkc-1,ivarc)&
                     - coar(iic  ,jjc-1,kkc+1,ivarc) &
                     - coar(iic  ,jjc+1,kkc-1,ivarc))/(four*dxc*dxc) 

                zcoar = dxc*(kkc + half)
                zfine = dxf*(kkf + half)
                zdist = zfine - zcoar
                dzz = (coar(iic  ,jjc  ,kkc+1, ivarc) + coar(iic  ,jjc  ,kkc-1, ivarc)-  two*coar(iic,jjc,kkc,ivarc))/(dxc*dxc) 
#endif
!               begin debug
                dxx = zero
                dyy = zero
                dxy = zero
!               end debug                
                coarval = coar(iic,jjc,kkc, ivarc) 
                finevalnew = coarval &
                     + xdist * xslope  &
                     + ydist * yslope  &
                     + zdist * zslope  &
                     + half*xdist*xdist*dxx &
                     + half*ydist*ydist*dyy &
                     + half*zdist*zdist*dzz &
                     + xdist*ydist*dxy &
                     + xdist*zdist*dxz &
                     + ydist*zdist*dyz 
                fine(iif, jjf, kkf, ivarf) = finevalnew

             enddo
          enddo
       enddo


    enddo


  end subroutine ebfnd_pwqinterp_nobound

end module ebfnd_pwqinterp_nobound_module
