subroutine lst_comp_norm ( &
     lo, hi, &
     soln, s_l1, s_h1, &
     exac, e_l1, e_h1, &
     mask, m_l1, m_h1, &
     volb, v_l1, v_h1, &
     norm2, norm0, volume, nsoln, iCpp, iF90)

  implicit none

  integer, intent(in) :: lo(1), hi(1)
  integer, intent(in) :: nsoln, iCpp, iF90
  integer, intent(in) :: s_l1, s_h1
  integer, intent(in) :: e_l1, e_h1
  integer, intent(in) :: m_l1, m_h1
  integer, intent(in) :: v_l1, v_h1
  double precision, intent(inout) :: norm2(0:nsoln-1), norm0(0:nsoln-1), volume
  double precision, intent(in) :: soln(s_l1:s_h1, 0:nsoln-1)
  double precision, intent(in) :: exac(e_l1:e_h1)
  double precision, intent(in) :: mask(m_l1:m_h1)
  double precision, intent(in) :: volb(v_l1:v_h1)

  integer :: i
  double precision :: error

  do i = lo(1), hi(1)
     if (mask(i) .gt. 0.d0) then
        volume = volume + volb(i)
        
        if (iCpp .ge. 0) then
           error = soln(i,iCpp)-exac(i)
           norm2(iCpp) = norm2(iCpp) + error**2 * volb(i)
           norm0(iCpp) = max(norm0(iCpp), abs(error))
        end if
        
        if (iF90 .ge. 0) then
           error = soln(i,iF90)-exac(i)
           norm2(iF90) = norm2(iF90) + error**2 * volb(i)
           norm0(iF90) = max(norm0(iF90), abs(error))
        end if
     end if
  end do

end subroutine lst_comp_norm
