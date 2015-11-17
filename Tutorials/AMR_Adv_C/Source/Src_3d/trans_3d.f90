
module transverse_module

  implicit none

  private

  public :: transy, transx, transx_term, transy_term

contains
  
  subroutine transy(lo, hi, hdtdy, &
                    vx, vxlo, vxhi, &
                    vdqdy, dlo, dhi, &
                    qx, qlo, qhi)
    integer, intent(in) :: lo(2), hi(2), vxlo(2), vxhi(2), dlo(2), dhi(2), qlo(2), qhi(2)
    double precision, intent(in) :: hdtdy
    double precision, intent(in   ) :: vx   (vxlo(1):vxhi(1),vxlo(2):vxhi(2))
    double precision, intent(  out) :: vdqdy( dlo(1): dhi(1), dlo(2): dhi(2))
    double precision, intent(inout) :: qx   ( qlo(1): qhi(1), qlo(2): qhi(2))

    integer :: i, j

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (vx(i,j) .ge. 0.d0) then
             qx(i,j) = qx(i,j) - hdtdy * vdqdy(i-1,j)
          else
             qx(i,j) = qx(i,j) - hdtdy * vdqdy(i  ,j)
          end if
       end do
    end do

  end subroutine transy

  subroutine transx(lo, hi, hdtdx, &
                    vy, vylo, vyhi, &
                    udqdx, dlo, dhi, &
                    qy, qlo, qhi)
    integer, intent(in) :: lo(2), hi(2), vylo(2), vyhi(2), dlo(2), dhi(2), qlo(2), qhi(2)
    double precision, intent(in) :: hdtdx
    double precision, intent(in   ) :: vy   (vylo(1):vyhi(1),vylo(2):vyhi(2))
    double precision, intent(  out) :: udqdx( dlo(1): dhi(1), dlo(2): dhi(2))
    double precision, intent(inout) :: qy   ( qlo(1): qhi(1), qlo(2): qhi(2))

    integer :: i, j

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (vy(i,j) .ge. 0.d0) then
             qy(i,j) = qy(i,j) - hdtdx * udqdx(i,j-1)
          else
             qy(i,j) = qy(i,j) - hdtdx * udqdx(i,j  )
          end if
       end do
    end do

  end subroutine transx


  subroutine transx_term(lo, hi, &
                         vx, vxlo, vxhi, &
                         qx, qxlo, qxhi, &
                         udqdx, dlo, dhi)
    integer, intent(in) :: lo(2), hi(2), vxlo(2), vxhi(2), qxlo(2), qxhi(2), dlo(2), dhi(2)
    double precision, intent(in ) :: vx   (vxlo(1):vxhi(1),vxlo(2):vxhi(2))
    double precision, intent(in ) :: qx   (qxlo(1):qxhi(1),qxlo(2):qxhi(2))
    double precision, intent(out) :: udqdx( dlo(1): dhi(1), dlo(2): dhi(2))

    integer :: i, j

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          udqdx(i,j) = 0.5d0*(vx(i,j)+vx(i+1,j))*(qx(i+1,j)-qx(i,j))
       end do
    end do

  end subroutine transx_term


  subroutine transy_term(lo, hi, &
                         vy, vylo, vyhi, &
                         qy, qylo, qyhi, &
                         vdqdy, dlo, dhi)
    integer, intent(in) :: lo(2), hi(2), vylo(2), vyhi(2), qylo(2), qyhi(2), dlo(2), dhi(2)
    double precision, intent(in ) :: vy   (vylo(1):vyhi(1),vylo(2):vyhi(2))
    double precision, intent(in ) :: qy   (qylo(1):qyhi(1),qylo(2):qyhi(2))
    double precision, intent(out) :: vdqdy( dlo(1): dhi(1), dlo(2): dhi(2))

    integer :: i, j

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          vdqdy(i,j) = 0.5d0*(vy(i,j)+vy(i,j+1))*(qy(i,j+1)-qy(i,j))
       end do
    end do

  end subroutine transy_term

end module transverse_module
