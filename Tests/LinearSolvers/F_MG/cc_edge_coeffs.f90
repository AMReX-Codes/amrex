module cc_edge_coeffs_module 

    use BoxLib
    use fabio_module
    use ml_layout_module
    use multifab_module
    use bl_error_module

    implicit none

contains

  subroutine make_edge_coeffs(mla, n, edge_coeffs, pd, coeffs_type, fabio)

    type(ml_layout), intent(inout) :: mla
    type(multifab ), intent(inout) :: edge_coeffs(:)
    type(box)     , intent(in   )  :: pd
    integer       , intent(in   )  :: n, coeffs_type
    logical       , intent(in   )  :: fabio

    real(kind=dp_t) :: c_norm
    type(multifab)  :: cell_coeffs

    call multifab_build(cell_coeffs,mla%la(n),1,1)

    if (coeffs_type .eq. 1 .or. coeffs_type .eq. 2) then
       call coeffs_init_1(cell_coeffs,edge_coeffs(:),pd, coeffs_type)
    else
       call bl_error("Dont know this coeffs type")
    end if

    c_norm = norm_inf(cell_coeffs)
    if ( parallel_ioprocessor() ) &
       print *,'Max norm of coefficient array ',c_norm

    if ( fabio ) then
      call fabio_write(cell_coeffs, "coeffs","Header")
    end if

    call multifab_destroy(cell_coeffs)

  end subroutine make_edge_coeffs

  subroutine coeffs_init_1(cell_coeffs,edge_coeffs,pd,coeffs_type)

    type( multifab), intent(inout) :: cell_coeffs
    type( multifab), intent(inout) :: edge_coeffs(:)
    type(box)      , intent(in   ) :: pd
    integer       , intent(in   )  :: coeffs_type

    type(box)                :: bx
    real(kind=dp_t)          :: dx
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: epx(:,:,:,:)
    real(kind=dp_t), pointer :: epy(:,:,:,:)
    real(kind=dp_t), pointer :: epz(:,:,:,:)

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

    call multifab_fill_boundary(cell_coeffs)

    do i = 1, nfabs(cell_coeffs)
       bx = get_ibox(cell_coeffs, i)
       cp  => dataptr(cell_coeffs,i)
       epx => dataptr(edge_coeffs(1),i)
       epy => dataptr(edge_coeffs(2),i)
       if (dm.eq.2) then
          call cc_to_edges_2d(cp(:,:,1,1),epx(:,:,1,1),epy(:,:,1,1),bx)
       else if (dm.eq.3) then
          epz => dataptr(edge_coeffs(3),i)
          call cc_to_edges_3d(cp(:,:,:,1),epx(:,:,:,1),epy(:,:,:,1),epz(:,:,:,1),bx)
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

  subroutine cc_to_edges_2d(cc_coeffs,x_edges,y_edges,bx)

      type(box)       , intent(in   ) :: bx
      double precision, intent(in   ) :: cc_coeffs(bx%lo(1)-1:,bx%lo(2)-1:)
      double precision, intent(inout) :: x_edges(bx%lo(1):,bx%lo(2):)
      double precision, intent(inout) :: y_edges(bx%lo(1):,bx%lo(2):)

      integer :: i,j

      do j = bx%lo(2),bx%hi(2)
      do i = bx%lo(1),bx%hi(1)+1
            x_edges(i,j) = 0.5d0 * (cc_coeffs(i,j) + cc_coeffs(i-1,j))
      end do
      end do

      do i = bx%lo(1),bx%hi(1)
      do j = bx%lo(2),bx%hi(2)+1
            y_edges(i,j) = 0.5d0 * (cc_coeffs(i,j) + cc_coeffs(i,j-1))
      end do
      end do

  end subroutine cc_to_edges_2d

  subroutine cc_to_edges_3d(cc_coeffs,x_edges,y_edges,z_edges,bx)

      type(box)       , intent(in   ) :: bx
      double precision, intent(in   ) :: cc_coeffs(bx%lo(1)-1:,bx%lo(2)-1:,bx%lo(3)-1:)
      double precision, intent(inout) :: x_edges(bx%lo(1):,bx%lo(2):,bx%lo(3):)
      double precision, intent(inout) :: y_edges(bx%lo(1):,bx%lo(2):,bx%lo(3):)
      double precision, intent(inout) :: z_edges(bx%lo(1):,bx%lo(2):,bx%lo(3):)

      integer :: i,j,k

      do k = bx%lo(3),bx%hi(3)
      do j = bx%lo(2),bx%hi(2)
      do i = bx%lo(1),bx%hi(1)+1
            x_edges(i,j,k) = 0.5d0 * (cc_coeffs(i,j,k) + cc_coeffs(i-1,j,k))
      end do
      end do
      end do

      do k = bx%lo(3),bx%hi(3)
      do j = bx%lo(2),bx%hi(2)+1
      do i = bx%lo(1),bx%hi(1)
            y_edges(i,j,k) = 0.5d0 * (cc_coeffs(i,j,k) + cc_coeffs(i,j-1,k))
      end do
      end do
      end do

      do k = bx%lo(3),bx%hi(3)+1
      do j = bx%lo(2),bx%hi(2)
      do i = bx%lo(1),bx%hi(1)
            z_edges(i,j,k) = 0.5d0 * (cc_coeffs(i,j,k) + cc_coeffs(i,j,k-1))
      end do
      end do
      end do

  end subroutine cc_to_edges_3d


end module cc_edge_coeffs_module
