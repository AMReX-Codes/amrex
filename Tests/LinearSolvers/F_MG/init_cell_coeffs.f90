module init_cell_coeffs_module 

    use BoxLib
    use ml_layout_module
    use multifab_module
    use bl_error_module

    implicit none

contains

  subroutine init_cell_coeffs(mla, cell_coeffs, pd, coeffs_type)

    type(ml_layout), intent(inout) :: mla
    type(multifab ), intent(inout) :: cell_coeffs
    type(box)      , intent(in   ) :: pd
    integer        , intent(in   ) :: coeffs_type

    real(kind=dp_t) :: c_norm

    if (coeffs_type .eq. 1 .or. coeffs_type .eq. 2) then
       call coeffs_init_1(cell_coeffs,pd,coeffs_type)
    else
       call bl_error("Dont know this coeffs type")
    end if

    call multifab_fill_boundary(cell_coeffs)

    c_norm = norm_inf(cell_coeffs)
    if ( parallel_ioprocessor() ) &
       print *,'Max norm of coefficient array ',c_norm

  end subroutine init_cell_coeffs

  subroutine coeffs_init_1(cell_coeffs,pd,coeffs_type)

    type( multifab), intent(inout) :: cell_coeffs
    type(box)      , intent(in   ) :: pd
    integer       , intent(in   )  :: coeffs_type

    type(box)                :: bx
    real(kind=dp_t)          :: dx
    real(kind=dp_t), pointer :: cp(:,:,:,:)

    integer         :: i,dm,nx,ny
    real(kind=dp_t) :: fac_2d, fac_3d

    nx = pd%hi(1) - pd%lo(1) + 1
    ny = pd%hi(2) - pd%lo(2) + 1

    dm = cell_coeffs%dim
    dx = 1.d0 / dble(nx)

    if (coeffs_type .eq. 1) then
        fac_2d = 160.d0
        fac_3d = 640.d0
    else if (coeffs_type .eq. 2) then
        fac_2d = 1600.d0
        fac_3d = 6400.d0
    end if

    do i = 1, nfabs(cell_coeffs)
       bx = get_ibox(cell_coeffs, i)
       cp => dataptr(cell_coeffs,i)
       if (dm.eq.2) then
          call init_coeffs_2d(cp(:,:,1,1),bx,dx,fac_2d)
       else if (dm.eq.3) then
          call init_coeffs_3d(cp(:,:,:,1),bx,dx,fac_3d)
       end if
    end do

  end subroutine coeffs_init_1

  subroutine init_coeffs_2d(cc_coeffs,bx,dx,fac)

      type(box)       , intent(in   ) :: bx
      double precision, intent(inout) :: cc_coeffs(bx%lo(1)-1:,bx%lo(2)-1:)
      double precision, intent(in   ) :: dx,fac

      integer :: i,j
      double precision :: x, y

      do j = bx%lo(2)-1,bx%hi(2)+1
      do i = bx%lo(1)-1,bx%hi(1)+1
         x = (dble(i)+0.5d0) * dx 
         y = (dble(j)+0.5d0) * dx 
         cc_coeffs(i,j) = 1.d0 + fac * (x**4 - x**2) * (y**4 - y**2)
      end do
      end do

  end subroutine init_coeffs_2d

  subroutine init_coeffs_3d(cell_coeffs,bx,dx,fac)

      type(box)       , intent(in   ) :: bx
      double precision, intent(inout) :: cell_coeffs(bx%lo(1)-1:,&
                                                     bx%lo(2)-1:,bx%lo(3)-1:)
      double precision, intent(in   ) :: dx,fac

      integer :: i,j,k
      double precision :: x, y, z

      do k = bx%lo(3)-1,bx%hi(3)+1
      do j = bx%lo(2)-1,bx%hi(2)+1
      do i = bx%lo(1)-1,bx%hi(1)+1
         x = (dble(i)+0.5d0) * dx
         y = (dble(j)+0.5d0) * dx
         z = (dble(k)+0.5d0) * dx
         cell_coeffs(i,j,k) = 1.d0 - fac * (x**4 - x**2) &
                            * (y**4 - y**2) * (z**4 - z**2)
      end do
      end do
      end do

  end subroutine init_coeffs_3d

end module init_cell_coeffs_module
