
module trace_module
  implicit none

  private
  public :: tracex, tracey

contains

  subroutine tracex(lo, hi, hdtdx, &
                    dq, dqlo, dqhi, &
                    q ,  qlo,  qhi, &
                    vx, vxlo, vxhi, &
                    qm, qp, mplo, mphi)

    use slope_module, only : slopex

    integer, intent(in) :: lo(2), hi(2), dqlo(2), dqhi(2), qlo(2), qhi(2),&
         vxlo(2), vxhi(2), mplo(2), mphi(2)
    double precision, intent(in)  :: hdtdx
    double precision              :: dq(dqlo(1):dqhi(1),dqlo(2):dqhi(2))
    double precision, intent(in ) ::  q( qlo(1): qhi(1), qlo(2): qhi(2))
    double precision, intent(in ) :: vx(vxlo(1):vxhi(1),vxlo(2):vxhi(2))
    double precision, intent(out) :: qm(mplo(1):mphi(1),mplo(2):mphi(2))
    double precision, intent(out) :: qp(mplo(1):mphi(1),mplo(2):mphi(2))
    
    integer :: i, j

    call slopex((/lo(1)-1,lo(2)/), hi, &
                q, qlo, qhi, &
                dq, dqlo, dqhi)
 
    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (vx(i,j) .ge. 0.d0) then
             qm(i,j) = q(i-1,j) + (0.5d0 - vx(i,j)*hdtdx) * dq(i-1,j)
             qp(i,j) = q(i,j)
          else
             qm(i,j) = q(i-1,j)
             qp(i,j) = q(i,j) - (0.5d0 + vx(i,j)*hdtdx) * dq(i,j)
          end if
       end do
    end do

  end subroutine tracex

  subroutine tracey(lo, hi, hdtdy, &
                    dq, dqlo, dqhi, &
                    q ,  qlo,  qhi, &
                    vy, vylo, vyhi, &
                    qm, qp, mplo, mphi)

    use slope_module, only : slopey

    integer, intent(in) :: lo(2), hi(2), dqlo(2), dqhi(2), qlo(2), qhi(2),&
         vylo(2), vyhi(2), mplo(2), mphi(2)
    double precision, intent(in)  :: hdtdy
    double precision              :: dq(dqlo(1):dqhi(1),dqlo(2):dqhi(2))
    double precision, intent(in ) ::  q( qlo(1): qhi(1), qlo(2): qhi(2))
    double precision, intent(in ) :: vy(vylo(1):vyhi(1),vylo(2):vyhi(2))
    double precision, intent(out) :: qm(mplo(1):mphi(1),mplo(2):mphi(2))
    double precision, intent(out) :: qp(mplo(1):mphi(1),mplo(2):mphi(2))

    integer :: i, j

    call slopey((/lo(1),lo(2)-1/), hi, &
                q, qlo, qhi, &
                dq, dqlo, dqhi)
 
    do j    = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (vy(i,j) .ge. 0.d0) then
             qm(i,j) = q(i,j-1) + (0.5d0 - vy(i,j)*hdtdy) * dq(i,j-1)
             qp(i,j) = q(i,j)
          else
             qm(i,j) = q(i,j-1)
             qp(i,j) = q(i,j) - (0.5d0 + vy(i,j)*hdtdy) * dq(i,j)
          end if
       end do
    end do
 
  end subroutine tracey

end module trace_module
