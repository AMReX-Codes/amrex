subroutine lst_comp_norm ( &
     lo, hi, &
     soln, s_l1, s_h1, &
     exac, e_l1, e_h1, &
     mask, m_l1, m_h1, &
     volb, v_l1, v_h1, &
     norm2, norm0, volume, nsoln)

  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: nsoln
  integer, intent(in) :: s_l1, s_h1
  integer, intent(in) :: e_l1, e_h1
  integer, intent(in) :: m_l1, m_h1
  integer, intent(in) :: v_l1, v_h1
  double precision, intent(inout) :: norm2(0:nsoln-1), norm0(0:nsoln-1), volume
  double precision, intent(in) :: soln(s_l1:s_h1, 0:nsoln-1)
  double precision, intent(in) :: exac(e_l1:e_h1)
  double precision, intent(in) :: mask(m_l1:m_h1)
  double precision, intent(in) :: volb(v_l1:v_h1)

  integer :: i, isoln
  double precision :: error

  do i = lo(1), hi(1)
     if (mask(i) .gt. 0.d0) then
        volume = volume + volb(i)
        
        do isoln = 0, nsoln-1
           error = soln(i,isoln)-exac(i)
           norm2(isoln) = norm2(isoln) + error**2 * volb(i)
           norm0(isoln) = max(norm0(isoln), abs(error))
        end do
     end if
  end do

end subroutine lst_comp_norm
