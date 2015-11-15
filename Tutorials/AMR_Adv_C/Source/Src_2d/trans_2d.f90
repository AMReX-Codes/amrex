
module transverse_module

  implicit none

  private

  public :: transy, transx

contains
  
  subroutine transy(lo, hi, hdtdy, &
                    qm, qp, qlo, qhi, &
                    fy, flo, fhi)
    integer, intent(in) :: lo(2), hi(2), qlo(2), qhi(2), flo(2), fhi(2)
    double precision, intent(in) :: hdtdy
    double precision, intent(inout) :: qm(qlo(1):qhi(1),qlo(2):qhi(2))
    double precision, intent(inout) :: qp(qlo(1):qhi(1),qlo(2):qhi(2))
    double precision, intent(in   ) :: fy(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i, j

    do j = lo(2), hi(2)
       i = lo(1)-1
       qm(i+1,j) = qm(i+1,j) + hdtdy * (fy(i,j)-fy(i,j+1)) 

       do i = lo(1), hi(1)-1
          qp(i  ,j) = qp(i  ,j) + hdtdy * (fy(i,j)-fy(i,j+1))
          qm(i+1,j) = qm(i+1,j) + hdtdy * (fy(i,j)-fy(i,j+1))
       end do

       i = hi(1)
       qp(i,j) = qp(i,j) + hdtdy * (fy(i,j)-fy(i,j+1))
    end do
  end subroutine transy

  subroutine transx(lo, hi, hdtdx, &
                    qm, qp, qlo, qhi, &
                    fx, flo, fhi)
    integer, intent(in) :: lo(2), hi(2), qlo(2), qhi(2), flo(2), fhi(2)
    double precision, intent(in) :: hdtdx
    double precision, intent(inout) :: qm(qlo(1):qhi(1),qlo(2):qhi(2))
    double precision, intent(inout) :: qp(qlo(1):qhi(1),qlo(2):qhi(2))
    double precision, intent(in   ) :: fx(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i, j

    j = lo(2)-1
    do i = lo(1), hi(1)
       qm(i,j+1) = qm(i,j+1) + hdtdx * (fx(i,j)-fx(i+1,j))
    end do

    do j = lo(2), hi(2)-1
       do i = lo(1), hi(1)
          qp(i,j  ) = qp(i,j  ) + hdtdx * (fx(i,j)-fx(i+1,j))
          qm(i,j+1) = qm(i,j+1) + hdtdx * (fx(i,j)-fx(i+1,j))
       end do
    end do

    j = hi(2)
    do i = lo(1), hi(1)
       qp(i,j  ) = qp(i,j  ) + hdtdx * (fx(i,j)-fx(i+1,j))
    end do

  end subroutine transx

end module transverse_module
