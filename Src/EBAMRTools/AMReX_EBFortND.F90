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
    integer      :: coar_lo(3),coar_hi(3)
    integer      :: fine_lo(3),fine_hi(3)
    integer      :: coarboxlo(3), coarboxhi(3)
    integer      :: refboxlo(3),   refboxhi(3)
    real(c_real) :: numfinepercoar
    real(c_real) :: fine(fine_lo(1):fine_hi(1),fine_lo(2):fine_hi(2),fine_lo(3):fine_hi(3), 0:fine_nco-1)
    real(c_real) :: coar(coar_lo(1):coar_hi(1),coar_lo(2):coar_hi(2),coar_lo(3):coar_hi(3), 0:coar_nco-1)

    numfinepercoar = 1.0
    do idir = 1, BL_SPACEDIM
       numfinepercoar = numfinepercoar*refrat
    enddo

    do ivar = 0, ncomp-1
       ivarf = isrc + ivar
       ivarc = idst + ivar
       do kkc = coarboxlo(3), coarboxhi(3)
          do jjc = coarboxlo(2), coarboxhi(2)
             do iic = coarboxlo(1), coarboxhi(1)

                coar(iic,jjc,kkc, ivarc) = 0.0

                do kb = refboxlo(3), refboxhi(3)
                   do jb = refboxlo(2), refboxhi(2)
                      do ib = refboxlo(1), refboxhi(1)

                         iif = refrat*iic + ib
                         jjf = refrat*jjc + jb
                         kkf = refrat*kkc + kb
                         coar(iic,jjc,kkc, ivarc) = coar(iic, jjc, kkc,ivarc) + fine(iif, jjf, kkf, ivarf)

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
    integer      :: coar_lo(3),coar_hi(3)
    integer      :: fine_lo(3),fine_hi(3)
    integer      :: fineboxlo(3), fineboxhi(3)

    real(c_real) :: fine(fine_lo(1):fine_hi(1),fine_lo(2):fine_hi(2),fine_lo(3):fine_hi(3), 0:fine_nco-1)
    real(c_real) :: coar(coar_lo(1):coar_hi(1),coar_lo(2):coar_hi(2),coar_lo(3):coar_hi(3), 0:coar_nco-1)

    do ivar = 0, ncomp-1
       ivarc = isrc + ivar
       ivarf = idst + ivar

       do kkf = fineboxlo(3), fineboxhi(3)
          do jjf = fineboxlo(2), fineboxhi(2)
             do iif = fineboxlo(1), fineboxhi(1)

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
