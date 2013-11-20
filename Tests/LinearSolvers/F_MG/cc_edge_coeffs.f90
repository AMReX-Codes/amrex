module cc_edge_coeffs_module 

    use BoxLib
    use ml_layout_module
    use multifab_module
    use bl_error_module

    implicit none

contains

  subroutine cell_to_edge_coeffs(cell_coeffs,edge_coeffs)

    type(multifab ), intent(in   ) :: cell_coeffs
    type(multifab ), intent(inout) :: edge_coeffs(:)

    type(box)                :: bx
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: epx(:,:,:,:)
    real(kind=dp_t), pointer :: epy(:,:,:,:)
    real(kind=dp_t), pointer :: epz(:,:,:,:)

    integer :: i, dm

    dm = cell_coeffs%dim

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

  end subroutine cell_to_edge_coeffs

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
