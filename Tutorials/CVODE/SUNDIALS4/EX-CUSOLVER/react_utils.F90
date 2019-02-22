module react_utils_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module

  implicit none

contains

  subroutine init_state(lo, hi, &
                        state, s_lo, s_hi, ncomp, npts) bind(C, name="init_state")

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), ncomp)
    integer, intent(in) :: npts, ncomp

    integer :: i, j, k, n

    n = 0
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             state(i, j, k, 1) = ONE - real(i, kind=rt)/real(npts, kind=rt)
             state(i, j, k, 2) = real(i, kind=rt)/real(npts, kind=rt)/2.0e0
             state(i, j, k, 3) = real(i, kind=rt)/real(npts, kind=rt)/2.0e0

             n = n + 1
          enddo
       enddo
    enddo

  end subroutine init_state


  subroutine get_state(state, s_lo, s_hi, ncomp, i, j, k, c, fval) bind(C, name="get_state")
    integer, intent(in) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), ncomp)
    integer, intent(in) :: ncomp
    integer, intent(in) :: i, j, k, c
    real(rt), intent(inout) :: fval
    fval = state(i, j, k, c)
  end subroutine get_state

  
  subroutine set_state(state, s_lo, s_hi, ncomp, i, j, k, c, fval) bind(C, name="set_state")
    integer, intent(in) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), ncomp)
    integer, intent(in) :: ncomp
    integer, intent(in) :: i, j, k, c
    real(rt), intent(in) :: fval
    state(i, j, k, c) = fval
  end subroutine set_state

end module react_utils_module
