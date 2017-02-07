module error_estimate_module

  use iso_c_binding

  use amrex_amr_module

  use my_amr_module, only : phi_new

  implicit none

  private

  public :: error_estimate

contains

  subroutine error_estimate (lev, cp, t, settag, cleartag) bind(c)
    integer, intent(in), value :: lev
    type(c_ptr), intent(in), value :: cp
    real(amrex_real), intent(in), value :: t
    character(c_char), intent(in), value :: settag, cleartag

    real(amrex_real), allocatable, save :: phierr(:)
    type(amrex_parmparse) :: pp
    type(amrex_tagboxarray) :: tag
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phiarr(:,:,:,:)
    character(c_char), contiguous, pointer :: tagarr(:,:,:,:)

    if (.not.allocated(phierr)) then
       call amrex_parmparse_build(pp, "myamr")
       call pp%getarr("phierr", phierr)
    end if

    tag%p = cp

    !$omp parallel private(mfi, bx, phiarr, tagarr)
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    do while(mfi%next())
       bx = mfi%tilebox()
       phiarr => phi_new(lev)%dataptr(mfi)
       tagarr => tag%dataptr(mfi)
       call tag_phi_error(bx%lo, bx%hi, &
            phiarr, lbound(phiarr), ubound(phiarr), &
            tagarr, lbound(tagarr), ubound(tagarr), &
            phierr(lev+1), settag, cleartag)  ! +1 because level starts with 0, but phierr starts with 1
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel

  end subroutine error_estimate

  subroutine tag_phi_error (lo, hi, phi, philo, phihi, tag, taglo, taghi, phierr, settag, cleartag)
    integer, intent(in) :: lo(3), hi(3), philo(4), phihi(4), taglo(4), taghi(4)
    real(amrex_real) , intent(in   ) :: phi(philo(1):phihi(1),philo(2):phihi(2),philo(3):phihi(3))
    character(c_char), intent(inout) :: tag(taglo(1):taghi(1),taglo(2):taghi(2),taglo(3):taghi(3))
    real(amrex_real), intent(in) :: phierr
    character(c_char), intent(in) :: settag, cleartag

    integer :: i,j,k

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (phi(i,j,k) .ge. phierr) then
                tag(i,j,k) = settag
             endif
          enddo
       enddo
    enddo
  end subroutine tag_phi_error

end module error_estimate_module

