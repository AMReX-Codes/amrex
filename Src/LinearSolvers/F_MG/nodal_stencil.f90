module nodal_stencil_module

  use bl_constants_module
  use bl_types
  use multifab_module
  use bc_functions_module
  use stencil_types_module

  implicit none

  public :: set_faces_edges_corners_2d,set_faces_edges_corners_3d

contains

  subroutine s_1d_nodal(ss, sg, mm, dh)
    real (kind = dp_t), intent(inout) :: ss(0:)
    real (kind = dp_t), intent(inout) :: sg(0:)
    integer           , intent(in   ) :: mm(:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    integer            :: nx
    real (kind = dp_t) :: f1

    nx = size(sg,dim=1) - 2

    if (bc_neumann(mm( 1),1,-1)) sg( 0) = sg(   1)
    if (bc_neumann(mm(nx),1,+1)) sg(nx) = sg(nx-1)

    f1 = ONE/dh(1)**2

    ss(:) = f1 * sg(:)

  end subroutine s_1d_nodal

  subroutine s_2d_nodal(ss, sg, mm, dh, lo, hi)

    integer           , intent(in   ) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: ss(lo(1)-1:,lo(2)-1:)
    real (kind = dp_t), intent(inout) :: sg(lo(1)-1:,lo(2)-1:)
    integer           , intent(in   ) :: mm(lo(1)  :,lo(2)  :)
    real (kind = dp_t), intent(in   ) :: dh(:)

    real (kind = dp_t) :: fac

    call set_faces_edges_corners_2d(sg, mm, lo, hi)

    ! This is not the full scaling for either stencil
    !   -- the cross stencil will need a factor of (1/2)
    !   -- the dense stencil will need a factor of (1/3)
    fac = ONE / (dh(1))**2

    ss(:,:) = fac * sg(:,:)

  end subroutine s_2d_nodal

  subroutine s_3d_nodal(ss, sg, mm, dh, lo, hi)

    integer           , intent(in   ) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: ss(0:,0:,0:)
    real (kind = dp_t), intent(inout) :: sg(0:,0:,0:)
    integer           , intent(inout) :: mm(:,:,:)
    real (kind = dp_t), intent(in   ) :: dh(:)

    real (kind = dp_t) :: fac

    call set_faces_edges_corners_3d(sg, mm, lo, hi)

    ! This is the right scaling for the cross stencil
    ! We end up needing to multiply by four later when we
    !    use this in the dense stencil.
    fac = (FOURTH / (dh(1))**2)

    ss(:,:,:) = fac * sg(:,:,:)

  end subroutine s_3d_nodal

  subroutine set_faces_edges_corners_2d(sg, mm, lo, hi)
 
    integer           , intent(in   ) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: sg(lo(1)-1:,lo(2)-1:)
    integer           , intent(in   ) :: mm(lo(1)  :,lo(2)  :)

    integer :: i, j
    !
    ! Set sg on edges at a Neumann boundary.
    !
    do i = lo(1),hi(1)
       if (bc_neumann(mm(i,lo(2)  ),2,-1)) sg(i,lo(2)-1) = sg(i,lo(2))
       if (bc_neumann(mm(i,hi(2)+1),2,+1)) sg(i,hi(2)+1) = sg(i,hi(2))
    end do

    do j = lo(2),hi(2)
       if (bc_neumann(mm(lo(1)  ,j),1,-1)) sg(lo(1)-1,j) = sg(lo(1),j)
       if (bc_neumann(mm(hi(1)+1,j),1,+1)) sg(hi(1)+1,j) = sg(hi(1),j)
    end do

    !
    ! Note: we do the corners *after* each of the edges has been done.
    !
    if (bc_neumann(mm(lo(1)  ,lo(2)  ),1,-1)) sg(lo(1)-1,lo(2)-1) = sg(lo(1),lo(2)-1)
    if (bc_neumann(mm(lo(1)  ,hi(2)+1),1,-1)) sg(lo(1)-1,hi(2)+1) = sg(lo(1),hi(2)+1)

    if (bc_neumann(mm(hi(1)+1,lo(2)  ),1,+1)) sg(hi(1)+1,lo(2)-1) = sg(hi(1),lo(2)-1)
    if (bc_neumann(mm(hi(1)+1,hi(2)+1),1,+1)) sg(hi(1)+1,hi(2)+1) = sg(hi(1),hi(2)+1)

    if (bc_neumann(mm(lo(1)  ,lo(2)  ),2,-1)) sg(lo(1)-1,lo(2)-1) = sg(lo(1)-1,lo(2))
    if (bc_neumann(mm(hi(1)+1,lo(2)  ),2,-1)) sg(hi(1)+1,lo(2)-1) = sg(hi(1)+1,lo(2))

    if (bc_neumann(mm(lo(1)  ,hi(2)+1),2,+1)) sg(lo(1)-1,hi(2)+1) = sg(lo(1)-1,hi(2))
    if (bc_neumann(mm(hi(1)+1,hi(2)+1),2,+1)) sg(hi(1)+1,hi(2)+1) = sg(hi(1)+1,hi(2))

  end subroutine set_faces_edges_corners_2d

  subroutine set_faces_edges_corners_3d(sg, mm, lo, hi)

    integer           , intent(in   ) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) :: sg(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    integer           , intent(in   ) :: mm(lo(1)  :,lo(2)  :,lo(3)  :)

    integer :: i, j, k
    !
    ! Set sg on faces at a Neumann boundary.
    !
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)
       if (bc_neumann(mm(i,j,lo(3)  ),3,-1)) sg(i,j,lo(3)-1) = sg(i,j,lo(3))
       if (bc_neumann(mm(i,j,hi(3)+1),3,+1)) sg(i,j,hi(3)+1) = sg(i,j,hi(3))
    end do
    end do

    do k = lo(3),hi(3)
    do i = lo(1),hi(1)
       if (bc_neumann(mm(i,lo(2)  ,k),2,-1)) sg(i,lo(2)-1,k) = sg(i,lo(2),k)
       if (bc_neumann(mm(i,hi(2)+1,k),2,+1)) sg(i,hi(2)+1,k) = sg(i,hi(2),k)
    end do
    end do

    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
       if (bc_neumann(mm(lo(1)  ,j,k),1,-1)) sg(lo(1)-1,j,k) = sg(lo(1),j,k)
       if (bc_neumann(mm(hi(1)+1,j,k),1,+1)) sg(hi(1)+1,j,k) = sg(hi(1),j,k)
    end do
    end do

    !
    ! Set sg on edges at a Neumann boundary.
    !
    do i = lo(1),hi(1)
       if (bc_neumann(mm(lo(1),lo(2)  ,lo(3)  ),2,-1)) sg(i,lo(2)-1,lo(3)-1) = sg(i,lo(2),lo(3)-1) 
       if (bc_neumann(mm(lo(1),lo(2)  ,hi(3)+1),2,-1)) sg(i,lo(2)-1,hi(3)+1) = sg(i,lo(2),hi(3)+1) 

       if (bc_neumann(mm(lo(1),hi(2)+1,lo(3)  ),2,+1)) sg(i,hi(2)+1,lo(3)-1) = sg(i,hi(2),lo(3)-1)
       if (bc_neumann(mm(lo(1),hi(2)+1,hi(3)+1),2,+1)) sg(i,hi(2)+1,hi(3)+1) = sg(i,hi(2),hi(3)+1)

       if (bc_neumann(mm(lo(1),lo(2)  ,lo(3)  ),3,-1)) sg(i,lo(2)-1,lo(3)-1) = sg(i,lo(2)-1,lo(3))
       if (bc_neumann(mm(lo(1),hi(2)+1,lo(3)  ),3,-1)) sg(i,hi(2)+1,lo(3)-1) = sg(i,hi(2)+1,lo(3))

       if (bc_neumann(mm(lo(1),lo(2)  ,hi(3)+1),3,+1)) sg(i,lo(2)-1,hi(3)+1) = sg(i,lo(2)-1,hi(3))
       if (bc_neumann(mm(lo(1),hi(2)+1,hi(3)+1),3,+1)) sg(i,hi(2)+1,hi(3)+1) = sg(i,hi(2)+1,hi(3))
    end do

    do j = lo(2),hi(2)
       if (bc_neumann(mm(lo(1)  ,j,lo(3)  ),1,-1)) sg(lo(1)-1,j,lo(3)-1) = sg(lo(1),j,lo(3)-1)
       if (bc_neumann(mm(lo(1)  ,j,hi(3)+1),1,-1)) sg(lo(1)-1,j,hi(3)+1) = sg(lo(1),j,hi(3)+1)

       if (bc_neumann(mm(hi(1)+1,j,lo(3)  ),1,+1)) sg(hi(1)+1,j,lo(3)-1) = sg(hi(1),j,lo(3)-1)
       if (bc_neumann(mm(hi(1)+1,j,hi(3)+1),1,+1)) sg(hi(1)+1,j,hi(3)+1) = sg(hi(1),j,hi(3)+1)

       if (bc_neumann(mm(lo(1)  ,j,lo(3)  ),3,-1)) sg(lo(1)-1,j,lo(3)-1) = sg(lo(1)-1,j,lo(3))  
       if (bc_neumann(mm(hi(1)+1,j,lo(3)  ),3,-1)) sg(hi(1)+1,j,lo(3)-1) = sg(hi(1)+1,j,lo(3))

       if (bc_neumann(mm(lo(1)  ,j,hi(3)+1),3,+1)) sg(lo(1)-1,j,hi(3)+1) = sg(lo(1)-1,j,hi(3))
       if (bc_neumann(mm(hi(1)+1,j,hi(3)+1),3,+1)) sg(hi(1)+1,j,hi(3)+1) = sg(hi(1)+1,j,hi(3))
    end do

    do k = lo(3),hi(3)
       if (bc_neumann(mm(lo(1)  ,lo(2)  ,k),1,-1)) sg(lo(1)-1,lo(2)-1,k) = sg(lo(1),lo(2)-1,k)
       if (bc_neumann(mm(lo(1)  ,hi(2)+1,k),1,-1)) sg(lo(1)-1,hi(2)+1,k) = sg(lo(1),hi(2)+1,k)

       if (bc_neumann(mm(hi(1)+1,lo(2)  ,k),1,+1)) sg(hi(1)+1,lo(2)-1,k) = sg(hi(1),lo(2)-1,k)
       if (bc_neumann(mm(hi(1)+1,hi(2)+1,k),1,+1)) sg(hi(1)+1,hi(2)+1,k) = sg(hi(1),hi(2)+1,k)

       if (bc_neumann(mm(lo(1)  ,lo(2)  ,k),2,-1)) sg(lo(1)-1,lo(2)-1,k) = sg(lo(1)-1,lo(2),k)
       if (bc_neumann(mm(hi(1)+1,lo(2)  ,k),2,-1)) sg(hi(1)+1,lo(2)-1,k) = sg(hi(1)+1,lo(2),k)

       if (bc_neumann(mm(lo(1)  ,hi(2)+1,k),2,+1)) sg(lo(1)-1,hi(2)+1,k) = sg(lo(1)-1,hi(2),k)
       if (bc_neumann(mm(hi(1)+1,hi(2)+1,k),2,+1)) sg(hi(1)+1,hi(2)+1,k) = sg(hi(1)+1,hi(2),k)
    end do

    if (bc_neumann(mm(lo(1),lo(2),lo(3)),1,-1)) sg(lo(1)-1,lo(2)-1,lo(3)-1) = sg(lo(1)  ,lo(2)-1,lo(3)-1) 
    if (bc_neumann(mm(lo(1),lo(2),lo(3)),2,-1)) sg(lo(1)-1,lo(2)-1,lo(3)-1) = sg(lo(1)-1,lo(2)  ,lo(3)-1) 
    if (bc_neumann(mm(lo(1),lo(2),lo(3)),3,-1)) sg(lo(1)-1,lo(2)-1,lo(3)-1) = sg(lo(1)-1,lo(2)-1,lo(3)  ) 

    if (bc_neumann(mm(hi(1)+1,lo(2),lo(3)),1,+1)) sg(hi(1)+1,lo(2)-1,lo(3)-1) = sg(hi(1)  ,lo(2)-1,lo(3)-1) 
    if (bc_neumann(mm(hi(1)+1,lo(2),lo(3)),2,-1)) sg(hi(1)+1,lo(2)-1,lo(3)-1) = sg(hi(1)+1,lo(2)  ,lo(3)-1) 
    if (bc_neumann(mm(hi(1)+1,lo(2),lo(3)),3,-1)) sg(hi(1)+1,lo(2)-1,lo(3)-1) = sg(hi(1)+1,lo(2)-1,lo(3))   

    if (bc_neumann(mm(lo(1),hi(2)+1,lo(3)),1,-1)) sg(lo(1)-1,hi(2)+1,lo(3)-1) = sg(lo(1)  ,hi(2)+1, lo(3)-1) 
    if (bc_neumann(mm(lo(1),hi(2)+1,lo(3)),2,+1)) sg(lo(1)-1,hi(2)+1,lo(3)-1) = sg(lo(1)-1,hi(2)  , lo(3)-1) 
    if (bc_neumann(mm(lo(1),hi(2)+1,lo(3)),3,-1)) sg(lo(1)-1,hi(2)+1,lo(3)-1) = sg(lo(1)-1,hi(2)+1, lo(3)  )   

    if (bc_neumann(mm(lo(1),lo(2),hi(3)+1),1,-1)) sg(lo(1)-1,lo(2)-1,hi(3)+1) = sg(lo(1)  ,lo(2)-1,hi(3)+1) 
    if (bc_neumann(mm(lo(1),lo(2),hi(3)+1),2,-1)) sg(lo(1)-1,lo(2)-1,hi(3)+1) = sg(lo(1)-1,lo(2)  ,hi(3)+1) 
    if (bc_neumann(mm(lo(1),lo(2),hi(3)+1),3,+1)) sg(lo(1)-1,lo(2)-1,hi(3)+1) = sg(lo(1)-1,lo(2)-1,hi(3)  ) 

    if (bc_neumann(mm(hi(1)+1,hi(2)+1,lo(3)),1,+1)) sg(hi(1)+1,hi(2)+1,lo(3)-1) = sg(hi(1)  ,hi(2)+1, lo(3)-1) 
    if (bc_neumann(mm(hi(1)+1,hi(2)+1,lo(3)),2,+1)) sg(hi(1)+1,hi(2)+1,lo(3)-1) = sg(hi(1)+1,hi(2)  , lo(3)-1) 
    if (bc_neumann(mm(hi(1)+1,hi(2)+1,lo(3)),3,-1)) sg(hi(1)+1,hi(2)+1,lo(3)-1) = sg(hi(1)+1,hi(2)+1, lo(3)  ) 

    if (bc_neumann(mm(hi(1)+1,lo(2),hi(3)+1),1,+1)) sg(hi(1)+1,lo(2)-1,hi(3)+1) = sg(hi(1)  ,lo(2)-1,hi(3)+1) 
    if (bc_neumann(mm(hi(1)+1,lo(2),hi(3)+1),2,-1)) sg(hi(1)+1,lo(2)-1,hi(3)+1) = sg(hi(1)+1,lo(2)  ,hi(3)+1) 
    if (bc_neumann(mm(hi(1)+1,lo(2),hi(3)+1),3,+1)) sg(hi(1)+1,lo(2)-1,hi(3)+1) = sg(hi(1)+1,lo(2)-1,hi(3)  ) 

    if (bc_neumann(mm(lo(1),hi(2)+1,hi(3)+1),1,-1)) sg(lo(1)-1,hi(2)+1,hi(3)+1) = sg(lo(1)  ,hi(2)+1,hi(3)+1) 
    if (bc_neumann(mm(lo(1),hi(2)+1,hi(3)+1),2,+1)) sg(lo(1)-1,hi(2)+1,hi(3)+1) = sg(lo(1)-1,hi(2)  ,hi(3)+1) 
    if (bc_neumann(mm(lo(1),hi(2)+1,hi(3)+1),3,+1)) sg(lo(1)-1,hi(2)+1,hi(3)+1) = sg(lo(1)-1,hi(2)+1,hi(3)  ) 

    if (bc_neumann(mm(hi(1)+1,hi(2)+1,hi(3)+1),1,+1)) sg(hi(1)+1,hi(2)+1,hi(3)+1) = sg(hi(1)  ,hi(2)+1,hi(3)+1) 
    if (bc_neumann(mm(hi(1)+1,hi(2)+1,hi(3)+1),2,+1)) sg(hi(1)+1,hi(2)+1,hi(3)+1) = sg(hi(1)+1,hi(2)  ,hi(3)+1) 
    if (bc_neumann(mm(hi(1)+1,hi(2)+1,hi(3)+1),3,+1)) sg(hi(1)+1,hi(2)+1,hi(3)+1) = sg(hi(1)+1,hi(2)+1,hi(3)) 

  end subroutine set_faces_edges_corners_3d

end module nodal_stencil_module
