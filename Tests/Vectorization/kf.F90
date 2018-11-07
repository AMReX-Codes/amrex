module kf_module

  implicit none

  double precision, parameter :: gamma = 1.4d0
  double precision, parameter :: smallr = 1.d-19
  double precision, parameter :: smallp = 1.d-10

  integer, parameter :: URHO  = 0;
  integer, parameter :: UMX   = 1;
  integer, parameter :: UMY   = 2;
  integer, parameter :: UMZ   = 3;
  integer, parameter :: UEDEN = 4;
  integer, parameter :: UEINT = 5;
  integer, parameter :: UTEMP = 6;
  integer, parameter :: ncons = 7
 
  integer, parameter :: QRHO   = 0;
  integer, parameter :: QU     = 1;
  integer, parameter :: QV     = 2;
  integer, parameter :: QW     = 3;
  integer, parameter :: QPRES  = 4;
  integer, parameter :: QCS    = 5;
  integer, parameter :: QEINT  = 6;
  integer, parameter :: QTEMP  = 7;
  integer, parameter :: nprim  = 8;

contains

  subroutine ctoprim_f (lo,hi,u,ulo,uhi,q,qlo,qhi) bind(c,name='ctoprim_f')
    integer, intent(in) :: lo(3),hi(3),ulo(3),uhi(3),qlo(3),qhi(3)
    double precision, intent(in   ) :: u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),0:ncons-1)
    double precision, intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),0:nprim-1)
    
    integer :: i,j,k
    double precision :: rhoinv, kineng
    
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             q(i,j,k,qrho) = max(smallr,u(i,j,k,urho))
             rhoinv = 1.d0/q(i,j,k,qrho)
             q(i,j,k,qu) = u(i,j,k,umx)*rhoinv
             q(i,j,k,qv) = u(i,j,k,umy)*rhoinv
             q(i,j,k,qw) = u(i,j,k,umz)*rhoinv
             kineng = 0.5d0*q(i,j,k,qrho)*(q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2)
             q(i,j,k,qeint) = (u(i,j,k,ueden)-kineng) * rhoinv
             if (q(i,j,k,qeint) .le. 0.d0) then
                q(i,j,k,qeint) = u(i,j,k,ueint) * rhoinv
             end if
             q(i,j,k,qpres) = max(smallp,(gamma-1.d0)*q(i,j,k,qeint)*q(i,j,k,qrho))
             q(i,j,k,qcs) = sqrt(gamma*q(i,j,k,qpres)*rhoinv)
             q(i,j,k,qtemp) = 0.d0
          end do
       end do
    end do
  end subroutine ctoprim_f
  
  subroutine flux_to_dudt_f (lo, hi, dudt, ulo, uhi, fx, xlo, xhi, &
       fy, ylo, yhi, fz, zlo, zhi, dxinv, nc) bind(c,name="flux_to_dudt_f")
    integer, dimension(3), intent(in) :: lo, hi, ulo, uhi, xlo, xhi, ylo, yhi, zlo, zhi
    integer, intent(in) :: nc
    double precision, intent(inout) :: dudt(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),nc)
    double precision, intent(inout) :: fx  (xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3),nc)
    double precision, intent(inout) :: fy  (ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3),nc)
    double precision, intent(inout) :: fz  (zlo(1):zhi(1),zlo(2):zhi(2),zlo(3):zhi(3),nc)
    double precision, intent(in) :: dxinv(3)
    
    integer :: i,j,k,n
    
    do n = 1, nc
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dudt(i,j,k,n) = dxinv(1) * (fx(i,j,k,n) - fx(i+1,j,k,n)) &
                     +          dxinv(2) * (fy(i,j,k,n) - fy(i,j+1,k,n)) &
                     +          dxinv(3) * (fz(i,j,k,n) - fz(i,j,k+1,n))
             end do
          end do
       end do
    end do
    
  end subroutine flux_to_dudt_f

end module kf_module
