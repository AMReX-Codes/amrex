
module ebfnd_ebamrtools_module

  !     since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
  !     for e.g., #if (BL_SPACEDIM == 1) statements.

  implicit none

  public

contains
  integer function imatebamrt(i, j)
    implicit none
    integer i, j, retval
    retval = 0
    if (i.eq.j) then
       retval = 1
    endif

    imatebamrt = retval
  end function imatebamrt

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

    numfinepercoar = 1.0d0
    do idir = 1, BL_SPACEDIM
       numfinepercoar = numfinepercoar*refrat
    enddo

    do ivar = 0, ncomp-1
       ivarf = isrc + ivar
       ivarc = idst + ivar
       do kkc = coarboxlo(2), coarboxhi(2)
          do jjc = coarboxlo(1), coarboxhi(1)
             do iic = coarboxlo(0), coarboxhi(0)

                coar(iic,jjc,kkc, ivarc) = 0.0d0

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


  subroutine ebfnd_average_face( &
       coar, coar_lo, coar_hi,  coar_nco, &
       fine, fine_lo, fine_hi,  fine_nco, &
       coarboxlo,coarboxhi, &
       refboxlo,refboxhi, &
       facedir, refrat, isrc, idst, ncomp) &
       bind(C, name="ebfnd_average_face")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: iif,jjf,kkf, coar_nco, fine_nco
    integer      :: iic,jjc,kkc, refrat, ncomp, isrc, idst
    integer      :: ivar, ivarc, ivarf, facedir
    integer      :: coar_lo(0:2),coar_hi(0:2)
    integer      :: fine_lo(0:2),fine_hi(0:2), ii, jj, kk
    integer      :: coarboxlo(0:2), coarboxhi(0:2)
    integer      :: refboxlo(0:2), refboxhi(0:2)
    real(c_real) :: fineval, coarval, numfinepercoar
    real(c_real) :: fine(fine_lo(0):fine_hi(0),fine_lo(1):fine_hi(1),fine_lo(2):fine_hi(2), 0:fine_nco-1)
    real(c_real) :: coar(coar_lo(0):coar_hi(0),coar_lo(1):coar_hi(1),coar_lo(2):coar_hi(2), 0:coar_nco-1)


    !number of fine faces per coarse face
    numfinepercoar = refrat
#if BL_SPACEDIM==3
    numfinepercoar = refrat*refrat
#endif

    do ivar = 0, ncomp-1
       ivarf = isrc + ivar
       ivarc = idst + ivar
       do kkc = coarboxlo(2), coarboxhi(2)
          do jjc = coarboxlo(1), coarboxhi(1)
             do iic = coarboxlo(0), coarboxhi(0)

                coar(iic,jjc,kkc, ivarc) = 0.0d0
                do kk = refboxlo(2), refboxhi(2)
                   do jj = refboxlo(1), refboxhi(1)
                      do ii = refboxlo(0), refboxhi(0)

                         iif = iic*refrat + ii
                         jjf = jjc*refrat + jj
                         kkf = kkc*refrat + kk


                         fineval = fine(iif, jjf, kkf, ivarf)
                         coarval = coar(iic, jjc, kkc, ivarc)
                         coar(iic,jjc,kkc, ivarc) = coarval + fineval

                      enddo
                   enddo
                enddo

                coar(iic,jjc,kkc, ivarc) = coar(iic, jjc, kkc,ivarc)/numfinepercoar

             enddo
          enddo
       enddo
    enddo


  end subroutine ebfnd_average_face


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

    dxf = 1.0d0
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

                xcoar = dxc*(iic + 0.5d0)
                ycoar = dxc*(jjc + 0.5d0)
                zcoar = 0.0d0

                xfine = dxf*(iif + 0.5d0)
                yfine = dxf*(jjf + 0.5d0)
                zfine = 0.0d0

                xslope = 0.5d0*(coar(iic+1,jjc  ,kkc  , ivarc)-  coar(iic-1,jjc  ,kkc  , ivarc))/dxc 
                yslope = 0.5d0*(coar(iic  ,jjc+1,kkc  , ivarc)-  coar(iic  ,jjc-1,kkc  , ivarc))/dxc
                zslope = 0.0d0

                xdist = xfine - xcoar
                ydist = yfine - ycoar
                zdist = 0.0d0

#if BL_SPACEDIM == 3
                zcoar = dxc*(kkc + 0.5d0)
                zfine = dxf*(kkf + 0.5d0)
                zdist = zfine - zcoar
                zslope = 0.5d0*(coar(iic  ,jjc  ,kkc+1, ivarc)-  coar(iic  ,jjc  ,kkc-1, ivarc))/dxc
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

    dxf = 1.0d0
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

                         xcoar = dxc*(indc + 0.5d0)
                         xfine = dxf*(indf + 0.5d0)
                         xdist = xfine - xcoar

                         coarhi = coar(iic+ioff,jjc+joff,kkc+koff, ivarc)
                         coarlo = coar(iic-ioff,jjc-joff,kkc-koff, ivarc)
                         xslope = (coarhi - coarlo)/(2.0d0*dxc)
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

                            xcoar = dxc*(indc + 0.5d0)
                            xfine = dxf*(indf + 0.5d0)
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

                            xcoar = dxc*(indc + 0.5d0)
                            xfine = dxf*(indf + 0.5d0)
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

    dxf = 1.0d0
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

                xcoar = dxc*(iic + 0.5d0)
                ycoar = dxc*(jjc + 0.5d0)
                zcoar = 0.0d0

                xfine = dxf*(iif + 0.5d0)
                yfine = dxf*(jjf + 0.5d0)
                zfine = 0.0d0
                
                coarhix = coar(iic+1,jjc  ,kkc  , ivarc)
                coarlox = coar(iic-1,jjc  ,kkc  , ivarc)
                coarhiy = coar(iic  ,jjc+1,kkc  , ivarc)
                coarloy = coar(iic  ,jjc-1,kkc  , ivarc)
                coarhiz = 0
                coarloz = 0
                xslope = (coarhix - coarlox)/(2.0d0*dxc) 
                yslope = (coarhiy - coarloy)/(2.0d0*dxc) 
                zslope = 0.0d0

                dxx = (  coar(iic+1,jjc  ,kkc  ,ivarc) &
                     +   coar(iic-1,jjc  ,kkc  ,ivarc)-  &
                     2.0d0*coar(iic  ,jjc  ,kkc  ,ivarc))/(dxc*dxc) 
                dyy = (  coar(iic  ,jjc+1,kkc  ,ivarc) &
                     +   coar(iic  ,jjc-1,kkc  ,ivarc)-  &
                     2.0d0*coar(iic  ,jjc  ,kkc  ,ivarc))/(dxc*dxc) 

                dxy = (coar(iic+1,jjc+1,kkc ,ivarc) &
                     + coar(iic-1,jjc-1,kkc, ivarc) &
                     - coar(iic+1,jjc-1,kkc ,ivarc) &
                     - coar(iic-1,jjc+1,kkc, ivarc))/(4.0d0*dxc*dxc) 
                dxz = 0.0d0
                dyz = 0.0d0
                dzz = 0.0d0
                xdist = xfine - xcoar
                ydist = yfine - ycoar
                zdist = 0.0d0

#if BL_SPACEDIM == 3
                coarhiz = coar(iic  ,jjc  ,kkc+1, ivarc)
                coarloz = coar(iic  ,jjc  ,kkc-1, ivarc)
                zslope = (coarhiz - coarloz)/(2.0d0*dxc) 
                dxz = (coar(iic+1,jjc  ,kkc+1,ivarc) &
                     + coar(iic-1,jjc  ,kkc-1,ivarc) &
                     - coar(iic+1,jjc  ,kkc-1,ivarc) &
                     - coar(iic-1,jjc  ,kkc+1,ivarc))/(4.0d0*dxc*dxc) 
                dyz = (coar(iic  ,jjc+1,kkc+1,ivarc)&
                     + coar(iic  ,jjc-1,kkc-1,ivarc)&
                     - coar(iic  ,jjc-1,kkc+1,ivarc) &
                     - coar(iic  ,jjc+1,kkc-1,ivarc))/(4.0d0*dxc*dxc) 

                zcoar = dxc*(kkc + 0.5d0)
                zfine = dxf*(kkf + 0.5d0)
                zdist = zfine - zcoar
                dzz = (coar(iic  ,jjc  ,kkc+1, ivarc) + coar(iic  ,jjc  ,kkc-1, ivarc)-  2.0d0*coar(iic,jjc,kkc,ivarc))/(dxc*dxc) 
#endif

                coarval = coar(iic,jjc,kkc, ivarc) 
                finevalnew = coarval &
                     + xdist * xslope  &
                     + ydist * yslope  &
                     + zdist * zslope  &
                     + 0.5d0*xdist*xdist*dxx &
                     + 0.5d0*ydist*ydist*dyy &
                     + 0.5d0*zdist*zdist*dzz &
                     + xdist*ydist*dxy &
                     + xdist*zdist*dxz &
                     + ydist*zdist*dyz 
                fine(iif, jjf, kkf, ivarf) = finevalnew

             enddo
          enddo
       enddo


    enddo


  end subroutine ebfnd_pwqinterp_nobound


  subroutine ebfnd_divflux( &
       divflux, divflux_lo, divflux_hi, divflux_nco,  &
       fluxfa0, fluxfa0_lo, fluxfa0_hi, fluxfa0_nco,  &
       fluxfa1, fluxfa1_lo, fluxfa1_hi, fluxfa1_nco,  &
       fluxfa2, fluxfa2_lo, fluxfa2_hi, fluxfa2_nco,  &
       gridlo,gridhi, &
       multbyarea, dx, isrc, idst, ncomp) &
       bind(C, name="ebfnd_divflux")

    use amrex_fort_module, only : dim=>amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: iif,jjf,kkf,  ncomp, ivar, ivardivf, ivarflux, multbyarea
    integer      :: fluxfa0_nco
    integer      :: fluxfa1_nco
    integer      :: fluxfa2_nco
    integer      :: divflux_nco, isrc, idst
    integer      :: fluxfa0_lo(0:2),fluxfa0_hi(0:2)
    integer      :: fluxfa1_lo(0:2),fluxfa1_hi(0:2)
    integer      :: fluxfa2_lo(0:2),fluxfa2_hi(0:2)
    integer      :: divflux_lo(0:2),divflux_hi(0:2)
    integer      :: gridlo(0:2), gridhi(0:2)

    real(c_real) :: dx, xterm, yterm, zterm
    real(c_real) :: denom, inv_denom
    integer      :: d
    real(c_real) :: divflux(divflux_lo(0):divflux_hi(0),divflux_lo(1):divflux_hi(1),divflux_lo(2):divflux_hi(2), 0:divflux_nco-1)
    real(c_real) :: fluxfa0(fluxfa0_lo(0):fluxfa0_hi(0),fluxfa0_lo(1):fluxfa0_hi(1),fluxfa0_lo(2):fluxfa0_hi(2), 0:fluxfa0_nco-1)
    real(c_real) :: fluxfa1(fluxfa1_lo(0):fluxfa1_hi(0),fluxfa1_lo(1):fluxfa1_hi(1),fluxfa1_lo(2):fluxfa1_hi(2), 0:fluxfa1_nco-1)
    real(c_real) :: fluxfa2(fluxfa2_lo(0):fluxfa2_hi(0),fluxfa2_lo(1):fluxfa2_hi(1),fluxfa2_lo(2):fluxfa2_hi(2), 0:fluxfa2_nco-1)

    if(multbyarea .eq. 1) then
       denom = dx
    else
       denom = 1.0d0
       do d=1,dim
          denom = denom * dx
       enddo
    endif
    inv_denom = 1.0d0/denom

    do ivar = 0, ncomp-1
       ivarflux = isrc + ivar
       ivardivf = idst + ivar

       do kkf = gridlo(2), gridhi(2)
          do jjf = gridlo(1), gridhi(1)
             do iif = gridlo(0), gridhi(0)

                xterm = fluxfa0(iif+1, jjf  , kkf  , ivarflux) - fluxfa0(iif, jjf, kkf, ivarflux)
                yterm = fluxfa1(iif  , jjf+1, kkf  , ivarflux) - fluxfa1(iif, jjf, kkf, ivarflux)
                zterm = 0.0d0
#if BL_SPACEDIM==3
                zterm = fluxfa2(iif  , jjf  , kkf+1, ivarflux) - fluxfa2(iif, jjf, kkf, ivarflux) 
#endif
                divflux(iif, jjf, kkf, ivardivf) = (xterm + yterm + zterm) * inv_denom

             enddo
          enddo
       enddo
    enddo


  end subroutine ebfnd_divflux


  subroutine ebfnd_gradlim( &
       gph, gph_lo, gph_hi, gph_nco,  &
       phi, phi_lo, phi_hi, phi_nco,  &
       gridlo,gridhi, &
       ncomp, dx, dolimiting) &
       bind(C, name="ebfnd_gradlim")

    use amrex_fort_module, only : amrex_spacedim, c_real=>amrex_real

    implicit none

    integer      :: i,j,k,  ncomp, ivar
    integer      :: phi_nco, gradcomp, vecdir,ii,jj,kk
    integer      :: gph_nco,  dolimiting
    integer      :: phi_lo(0:2), phi_hi(0:2)
    integer      :: gph_lo(0:2), gph_hi(0:2)
    integer      :: gridlo(0:2), gridhi(0:2)

    real(c_real) :: dx, dplo, dphi, dpce, dplim
    real(c_real) :: gph(gph_lo(0):gph_hi(0),gph_lo(1):gph_hi(1),gph_lo(2):gph_hi(2), 0:gph_nco-1)
    real(c_real) :: phi(phi_lo(0):phi_hi(0),phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2), 0:phi_nco-1)

    if(dolimiting .eq. 1) then
       do ivar = 0, ncomp-1
          do vecdir = 0, BL_SPACEDIM-1

             gradcomp = BL_SPACEDIM*ivar + vecdir
             ii = imatebamrt(0, vecdir)
             jj = imatebamrt(1, vecdir)
             kk = imatebamrt(2, vecdir)

             do k = gridlo(2), gridhi(2)
                do j = gridlo(1), gridhi(1)
                   do i = gridlo(0), gridhi(0)

                      dphi = phi(i+ii, j+jj, k+kk, ivar) - phi(i   ,j   ,k   , ivar) 
                      dplo = phi(i   , j   , k   , ivar) - phi(i-ii,j-jj,k-kk, ivar) 
                      dpce = 0.5d0*(dplo + dphi)
                      !van Leer limiting
                      if(dplo*dphi .lt. 0.0d0) then
                         dplim = 0.0d0
                      else
                         dplim = min(2.0d0*abs(dplo), 2.0d0*abs(dphi))
                         dplim = min(dplim, abs(dpce))
                         dplim = dplim*sign(1.0d0, dplo)
                      endif

                      gph(i, j, k, gradcomp) = dplim/dx

                   enddo
                enddo
             enddo
          enddo
       enddo
    else
       do ivar = 0, ncomp-1
          do vecdir = 0, BL_SPACEDIM-1

             gradcomp = BL_SPACEDIM*ivar + vecdir
             ii = imatebamrt(0, vecdir)
             jj = imatebamrt(1, vecdir)
             kk = imatebamrt(2, vecdir)

             do k = gridlo(2), gridhi(2)
                do j = gridlo(1), gridhi(1)
                   do i = gridlo(0), gridhi(0)
                      dphi = phi(i+ii, j+jj, k+kk, ivar) - phi(i   ,j   ,k   , ivar) 
                      dplo = phi(i   , j   , k   , ivar) - phi(i-ii,j-jj,k-kk, ivar) 
                      dpce = 0.5d0*(dplo + dphi)
                      ! no limiting here
                      gph(i, j , k , gradcomp) = dpce/dx

                   enddo
                enddo
             enddo
          enddo
       enddo
    endif


  end subroutine ebfnd_gradlim

end module ebfnd_ebamrtools_module
