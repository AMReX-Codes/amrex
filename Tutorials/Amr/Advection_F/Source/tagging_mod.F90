module tagging_module

  use iso_c_binding

  use amrex_fort_module

  implicit none

  private

  public :: tag_phi_error

contains

  subroutine tag_phi_error (level, time, lo, hi, phi, philo, phihi, tag, taglo, taghi, phierr, &
       settag, cleartag)
    integer, intent(in) :: level, lo(3), hi(3), philo(4), phihi(4), taglo(4), taghi(4)
    real(amrex_real) , intent(in   ) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
    character(kind=c_char), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(amrex_real), intent(in) :: time, phierr
    character(kind=c_char), intent(in) :: settag, cleartag

    integer :: i,j,k

    real(amrex_real) :: phierr_this

    phierr_this = phierr

    if (level >= 1) then
       !
       !  This is here for testing only! 
       !  Remove this if you use this as a template.
       !
       if (time > 0.75_amrex_real .and. time < 1.0_amrex_real) then
          phierr_this = 2.0_amrex_real * phierr
       end if
    end if

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (phi(i,j,k) .ge. phierr_this) then
                tag(i,j,k) = settag
             endif
          enddo
       enddo
    enddo
  end subroutine tag_phi_error

end module tagging_module

