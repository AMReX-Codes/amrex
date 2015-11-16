
module flux_module
  implicit none

  private

  public :: cmpflx

contains

  subroutine cmpflx(lo, hi, qm, qp, qlo, qhi, v, vlo, vhi, f, flo, fhi)
    integer, intent(in) :: lo(2), hi(2), qlo(2), qhi(2), vlo(2), vhi(2), flo(2), fhi(2)
    double precision, intent(in   ) :: qm(qlo(1):qhi(1),qlo(2):qhi(2))
    double precision, intent(in   ) :: qp(qlo(1):qhi(1),qlo(2):qhi(2))
    double precision, intent(in   ) :: v (vlo(1):vhi(1),vlo(2):vhi(2))
    double precision, intent(inout) :: f (flo(1):fhi(1),flo(2):fhi(2))

    integer :: i, j

    ! Upwind 
    do    j = lo(2),hi(2)
       do i = lo(1),hi(1)
          if (v(i,j) .gt. 0.d0) then
             f(i,j) = qm(i,j) * v(i,j)
          else if (v(i,j) .lt. 0.d0) then
             f(i,j) = qp(i,j) * v(i,j)
          else 
             f(i,j) = 0.5d0 * (qm(i,j) + qp(i,j)) * v(i,j)
          endif
       end do
    end do

  end subroutine cmpflx

end module flux_module
