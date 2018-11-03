
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
