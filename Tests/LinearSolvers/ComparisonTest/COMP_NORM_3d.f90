subroutine lst_comp_norm ( &
     lo, hi, &
     soln, s_l1, s_l2, s_l3, s_h1, s_h2, s_h3,  &
     exac, e_l1, e_l2, e_l3, e_h1, e_h2, e_h3,  &
     mask, m_l1, m_l2, m_l3, m_h1, m_h2, m_h3,  &
     volb, v_l1, v_l2, v_l3, v_h1, v_h2, v_h3,  &
     norm2, norm0, volume, nsoln)

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: nsoln
  integer, intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
  integer, intent(in) :: e_l1, e_l2, e_l3, e_h1, e_h2, e_h3
  integer, intent(in) :: m_l1, m_l2, m_l3, m_h1, m_h2, m_h3
  integer, intent(in) :: v_l1, v_l2, v_l3, v_h1, v_h2, v_h3
  double precision, intent(inout) :: norm2(0:nsoln-1), norm0(0:nsoln-1), volume
  double precision, intent(in) :: soln(s_l1:s_h1, s_l2:s_h2, s_l3:s_h3, 0:nsoln-1)
  double precision, intent(in) :: exac(e_l1:e_h1, e_l2:e_h2, e_l3:e_h3)
  double precision, intent(in) :: mask(m_l1:m_h1, m_l2:m_h2, m_l3:m_h3)
  double precision, intent(in) :: volb(v_l1:v_h1, v_l2:v_h2, v_l3:v_h3)

  integer :: i, j, k, isoln
  double precision :: error

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           if (mask(i,j,k) .gt. 0.d0) then
              volume = volume + volb(i,j,k)

              do isoln = 0, nsoln-1
                 error = soln(i,j,k,isoln)-exac(i,j,k)
                 norm2(isoln) = norm2(isoln) + error**2 * volb(i,j,k)
                 norm0(isoln) = max(norm0(isoln), abs(error))
              end do
           end if
        end do
     end do
  end do

end subroutine lst_comp_norm
