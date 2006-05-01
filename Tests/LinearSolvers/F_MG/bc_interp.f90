module intrp_module

  use bl_types

  implicit none

  integer, parameter, private :: d = dp_t

  real(d), parameter :: dirichlet(5,7) = &
       reshape((/ &
                2.0_d,  -1.0_d,        0.0_d,        0.0_d,       0.0_d,        0.0_d, 0.0_d, &
          8.0_d/3.0_d,  -2.0_d,  1.0_d/3.0_d,        0.0_d,       0.0_d,        0.0_d, 0.0_d, &
         16.0_d/5.0_d,  -3.0_d,        1.0_d, -1.0_d/5.0_d,       0.0_d,        0.0_d, 0.0_d, &
       128.0_d/35.0_d,  -4.0_d,        2.0_d, -4.0_d/5.0_d, 1.0_d/7.0_d,        0.0_d, 0.0_d, &
       256.0_d/63.0_d,  -5.0_d, 10.0_d/3.0_d,       -2.0_d, 5.0_d/7.0_d, -1.0_d/9.0_d, 0.0_d  &
       /), (/5,7/))

  real(d), parameter :: neumann(5,7) = &
       reshape((/ &
                 -1.0_d,           1.0_d,             0.0_d,            0.0_d,           0.0_d,            0.0_d, 0.0_d, &
                 -1.0_d,           1.0_d,             0.0_d,            0.0_d,           0.0_d,            0.0_d, 0.0_d, &
         -24.0_d/23.0_d,   21.0_d/23.0_d,      3.0_d/23.0_d,    -1.0_d/23.0_d,           0.0_d,            0.0_d, 0.0_d, & 
         -12.0_d/11.0_d,   17.0_d/22.0_d,      9.0_d/22.0_d,    -5.0_d/22.0_d,    1.0_d/22.0_d,            0.0_d, 0.0_d, &
       -640.0_d/563.0_d, 335.0_d/563.0_d, 1430.0_d/1689.0_d, -370.0_d/563.0_d, 145.0_d/563.0_d, -71.0_d/1689.0_d, 0.0_d  &
       /), (/5,7/))

  real(d), parameter :: corn(5,6) = &
       reshape((/ &
       1.0_d,   0.0_d,  0.0_d,  0.0_d, 0.0_d, 0.0_d, &
       2.0_d,  -1.0_d,  0.0_d,  0.0_d, 0.0_d, 0.0_d, &
       3.0_d,  -3.0_d,  1.0_d,  0.0_d, 0.0_d, 0.0_d, &
       4.0_d,  -6.0_d,  4.0_d, -1.0_d, 0.0_d, 0.0_d, &
       5.0_d, -10.0_d, 10.0_d, -5.0_d, 1.0_d, 0.0_d  &
       /), (/5,6/))

  type bc_rec
     integer :: dim = 0
     integer, pointer :: ord(:) => Null()
     real(dp_t), pointer :: cf(:,:) => Null()
     logical, pointer :: m_xl(:) => Null()
     logical, pointer :: m_xh(:) => Null()
     logical, pointer :: m_yl(:) => Null()
     logical, pointer :: m_yh(:) => Null()
  end type bc_rec

contains

  subroutine bc_interp_2d(ff, lo, bcr, &
       v_xl, v_xh, &
       v_yl, v_yh  &
       )
       
    integer, intent(in) :: lo(:)
    real(dp_t), intent(inout) :: ff(lo(1)-1:,lo(2)-1:)
    real(dp_t), intent(in) :: v_xl(lo(2):)
    real(dp_t), intent(in) :: v_xh(lo(2):)
    real(dp_t), intent(in) :: v_yl(lo(2):)
    real(dp_t), intent(in) :: v_yh(lo(2):)
    type(bc_rec), intent(in) :: bcr
    integer :: n, i, j
    integer :: hi(size(lo))

    hi(1) = ubound(ff,1)-1; hi(2) = ubound(ff,2)-1

    ! X Face
    do j = lo(2), hi(2)
       if ( bcr%m_xl(j) ) then
          ff(lo(1)-1, j) = v_xl(j)*bcr%cf(1,0) + sum(ff(lo(1):lo(1)+bcr%ord(1):+1,j)*bcr%cf(1,1:bcr%ord(1)))
       end if
       if ( bcr%m_xh(j) ) then
          ff(hi(1)+1, j) = v_xh(j)*bcr%cf(1,0) + sum(ff(hi(1):hi(1)-bcr%ord(1)+1:-1,j)*bcr%cf(1,1:bcr%ord(1)))
       end if
    end do

    ! Y Face
    do i = lo(1), hi(1)
       if ( bcr%m_yl(i) ) then
          ff(i,lo(2)-1) = v_yl(i)*bcr%cf(2,0) + sum(ff(i,lo(2):lo(2)+bcr%ord(2):+1)*bcr%cf(2,1:bcr%ord(2)))
       end if
       if ( bcr%m_yh(i) ) then
          ff(i,hi(2)+1) = v_yh(i)*bcr%cf(2,0) + sum(ff(i,hi(2):hi(2)-bcr%ord(2)+1:-1)*bcr%cf(2,1:bcr%ord(2)))
       end if
    end do

  end subroutine bc_interp_2d

end module intrp_module
