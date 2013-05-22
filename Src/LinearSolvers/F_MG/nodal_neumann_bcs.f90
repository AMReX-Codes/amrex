module impose_neumann_bcs_module

  use bl_types
  use bc_functions_module

  implicit none

contains

  subroutine impose_neumann_bcs_1d(uu,mm,lo,ng)

    integer, intent(in) :: ng,lo(:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:)
    integer           , intent(in   ) :: mm(lo(1)   :)
    integer :: hi(1)

    hi(1) = lo(1) + size(mm,dim=1)-1

    if (bc_neumann(mm(lo(1)),1,-1)) uu(lo(1)-1) = uu(lo(1)+1)
    if (bc_neumann(mm(hi(1)),1,+1)) uu(hi(1)+1) = uu(hi(1)-1)

  end subroutine impose_neumann_bcs_1d

  subroutine impose_neumann_bcs_2d(uu,mm,lo,ng)

    integer, intent(in) :: ng,lo(:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:)
    integer           , intent(in   ) :: mm(lo(1)   :,lo(2)   :)
    integer :: i, j, hi(2)

    hi(1) = lo(1) + size(mm,dim=1)-1
    hi(2) = lo(2) + size(mm,dim=2)-1

    do i = lo(1),hi(1)
       if (bc_neumann(mm(i,lo(2)),2,-1)) uu(i,lo(2)-1) = uu(i,lo(2)+1)
       if (bc_neumann(mm(i,hi(2)),2,+1)) uu(i,hi(2)+1) = uu(i,hi(2)-1)
    end do

    do j = lo(2),hi(2)
       if (bc_neumann(mm(lo(1),j),1,-1)) uu(lo(1)-1,j) = uu(lo(1)+1,j)
       if (bc_neumann(mm(hi(1),j),1,+1)) uu(hi(1)+1,j) = uu(hi(1)-1,j)
    end do

    if (bc_neumann(mm(lo(1),lo(2)),1,-1)) uu(lo(1)-1,lo(2)-1) = uu(lo(1)+1,lo(2)-1) 
    if (bc_neumann(mm(lo(1),lo(2)),2,-1)) uu(lo(1)-1,lo(2)-1) = uu(lo(1)-1,lo(2)+1) 

    if (bc_neumann(mm(hi(1),lo(2)),1,+1)) uu(hi(1)+1,lo(2)-1) = uu(hi(1)-1,lo(2)-1) 
    if (bc_neumann(mm(hi(1),lo(2)),2,-1)) uu(hi(1)+1,lo(2)-1) = uu(hi(1)+1,lo(2)+1) 

    if (bc_neumann(mm(lo(1),hi(2)),1,-1)) uu(lo(1)-1,hi(2)+1) = uu(lo(1)+1,hi(2)+1) 
    if (bc_neumann(mm(lo(1),hi(2)),2,+1)) uu(lo(1)-1,hi(2)+1) = uu(lo(1)-1,hi(2)-1) 

    if (bc_neumann(mm(hi(1),hi(2)),1,+1)) uu(hi(1)+1,hi(2)+1) = uu(hi(1)-1,hi(2)+1) 
    if (bc_neumann(mm(hi(1),hi(2)),2,+1)) uu(hi(1)+1,hi(2)+1) = uu(hi(1)+1,hi(2)-1) 

  end subroutine impose_neumann_bcs_2d

  subroutine impose_neumann_bcs_3d(uu,mm,lo,ng)

    integer, intent(in) :: ng,lo(:)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer           , intent(in   ) :: mm(lo(1)   :,lo(2)   :,lo(3)   :)
    integer :: i, j, k, hi(3)

    hi(1) = lo(1) + size(mm,dim=1)-1
    hi(2) = lo(2) + size(mm,dim=2)-1
    hi(3) = lo(3) + size(mm,dim=3)-1
    !
    ! Faces
    !
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          i = lo(1)
          if (bc_neumann(mm(i,j,k),1,-1)) uu(i-1,j,k) = uu(i+1,j,k)
          i = hi(1)
          if (bc_neumann(mm(i,j,k),1,+1)) uu(i+1,j,k) = uu(i-1,j,k)
       end do
    end do

    do k = lo(3),hi(3)
       do i = lo(1),hi(1)
          j = lo(2)
          if (bc_neumann(mm(i,j,k),2,-1)) uu(i,j-1,k) = uu(i,j+1,k)
          j = hi(2)
          if (bc_neumann(mm(i,j,k),2,+1)) uu(i,j+1,k) = uu(i,j-1,k)
       end do
    end do

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          k = lo(3)
          if (bc_neumann(mm(i,j,k),3,-1)) uu(i,j,k-1) = uu(i,j,k+1)
          k = hi(3)
          if (bc_neumann(mm(i,j,k),3,+1)) uu(i,j,k+1) = uu(i,j,k-1)
       end do
    end do
    !
    ! Edges
    !
    do i = lo(1),hi(1)
      if (bc_neumann(mm(i,lo(2),lo(3)),2,-1)) uu(i,lo(2)-1,lo(3)-1) = uu(i,lo(2)+1,lo(3)-1)
      if (bc_neumann(mm(i,lo(2),hi(3)),2,-1)) uu(i,lo(2)-1,hi(3)+1) = uu(i,lo(2)+1,hi(3)+1)

      if (bc_neumann(mm(i,hi(2),lo(3)),2,+1)) uu(i,hi(2)+1,lo(3)-1) = uu(i,hi(2)-1,lo(3)-1)
      if (bc_neumann(mm(i,hi(2),hi(3)),2,+1)) uu(i,hi(2)+1,hi(3)+1) = uu(i,hi(2)-1,hi(3)+1)

      if (bc_neumann(mm(i,lo(2),lo(3)),3,-1)) uu(i,lo(2)-1,lo(3)-1) = uu(i,lo(2)-1,lo(3)+1)
      if (bc_neumann(mm(i,hi(2),lo(3)),3,-1)) uu(i,hi(2)+1,lo(3)-1) = uu(i,hi(2)+1,lo(3)+1)

      if (bc_neumann(mm(i,lo(2),hi(3)),3,+1)) uu(i,lo(2)-1,hi(3)+1) = uu(i,lo(2)-1,hi(3)-1)
      if (bc_neumann(mm(i,hi(2),hi(3)),3,+1)) uu(i,hi(2)+1,hi(3)+1) = uu(i,hi(2)+1,hi(3)-1)
    end do

    do j = lo(2),hi(2)
      if (bc_neumann(mm(lo(1),j,lo(3)),1,-1)) uu(lo(1)-1,j,lo(3)-1) = uu(lo(1)+1,j,lo(3)-1)
      if (bc_neumann(mm(lo(1),j,hi(3)),1,-1)) uu(lo(1)-1,j,hi(3)+1) = uu(lo(1)+1,j,hi(3)+1)

      if (bc_neumann(mm(hi(1),j,lo(3)),1,+1)) uu(hi(1)+1,j,lo(3)-1) = uu(hi(1)-1,j,lo(3)-1)
      if (bc_neumann(mm(hi(1),j,hi(3)),1,+1)) uu(hi(1)+1,j,hi(3)+1) = uu(hi(1)-1,j,hi(3)+1)

      if (bc_neumann(mm(lo(1),j,lo(3)),3,-1)) uu(lo(1)-1,j,lo(3)-1) = uu(lo(1)-1,j,lo(3)+1)
      if (bc_neumann(mm(hi(1),j,lo(3)),3,-1)) uu(hi(1)+1,j,lo(3)-1) = uu(hi(1)+1,j,lo(3)+1)

      if (bc_neumann(mm(lo(1),j,hi(3)),3,+1)) uu(lo(1)-1,j,hi(3)+1) = uu(lo(1)-1,j,hi(3)-1)
      if (bc_neumann(mm(hi(1),j,hi(3)),3,+1)) uu(hi(1)+1,j,hi(3)+1) = uu(hi(1)+1,j,hi(3)-1)
    end do

    do k = lo(3),hi(3)
      if (bc_neumann(mm(lo(1),lo(2),k),1,-1)) uu(lo(1)-1,lo(2)-1,k) = uu(lo(1)+1,lo(2)-1,k)
      if (bc_neumann(mm(lo(1),hi(2),k),1,-1)) uu(lo(1)-1,hi(2)+1,k) = uu(lo(1)+1,hi(2)+1,k)

      if (bc_neumann(mm(hi(1),lo(2),k),1,+1)) uu(hi(1)+1,lo(2)-1,k) = uu(hi(1)-1,lo(2)-1,k)
      if (bc_neumann(mm(hi(1),hi(2),k),1,+1)) uu(hi(1)+1,hi(2)+1,k) = uu(hi(1)-1,hi(2)+1,k)

      if (bc_neumann(mm(lo(1),lo(2),k),2,-1)) uu(lo(1)-1,lo(2)-1,k) = uu(lo(1)-1,lo(2)+1,k)
      if (bc_neumann(mm(hi(1),lo(2),k),2,-1)) uu(hi(1)+1,lo(2)-1,k) = uu(hi(1)+1,lo(2)+1,k)

      if (bc_neumann(mm(lo(1),hi(2),k),2,+1)) uu(lo(1)-1,hi(2)+1,k) = uu(lo(1)-1,hi(2)-1,k)
      if (bc_neumann(mm(hi(1),hi(2),k),2,+1)) uu(hi(1)+1,hi(2)+1,k) = uu(hi(1)+1,hi(2)-1,k)

    end do
    !
    ! Corners
    !
    i = lo(1)
    j = lo(2)
    k = lo(3)
    if (bc_neumann(mm(i,j,k),1,-1))  uu(i-1,j-1,k-1) = uu(i+1,j-1,k-1)
    if (bc_neumann(mm(i,j,k),2,-1))  uu(i-1,j-1,k-1) = uu(i-1,j+1,k-1)
    if (bc_neumann(mm(i,j,k),3,-1))  uu(i-1,j-1,k-1) = uu(i-1,j-1,k+1)

    i = hi(1)
    j = lo(2)
    k = lo(3)
    if (bc_neumann(mm(i,j,k),1,+1))  uu(i+1,j-1,k-1) = uu(i-1,j-1,k-1)
    if (bc_neumann(mm(i,j,k),2,-1))  uu(i+1,j-1,k-1) = uu(i+1,j+1,k-1)
    if (bc_neumann(mm(i,j,k),3,-1))  uu(i+1,j-1,k-1) = uu(i+1,j-1,k+1)

    i = lo(1)
    j = hi(2)
    k = lo(3)
    if (bc_neumann(mm(i,j,k),1,-1))  uu(i-1,j+1,k-1) = uu(i+1,j+1,k-1)
    if (bc_neumann(mm(i,j,k),2,+1))  uu(i-1,j+1,k-1) = uu(i-1,j-1,k-1)
    if (bc_neumann(mm(i,j,k),3,-1))  uu(i-1,j+1,k-1) = uu(i-1,j+1,k+1)

    i = lo(1)
    j = lo(2)
    k = hi(3)
    if (bc_neumann(mm(i,j,k),1,-1))  uu(i-1,j-1,k+1) = uu(i+1,j-1,k+1)
    if (bc_neumann(mm(i,j,k),2,-1))  uu(i-1,j-1,k+1) = uu(i-1,j+1,k+1)
    if (bc_neumann(mm(i,j,k),3,+1))  uu(i-1,j-1,k+1) = uu(i-1,j-1,k-1)

    i = hi(1)
    j = hi(2)
    k = lo(3)
    if (bc_neumann(mm(i,j,k),1,+1))  uu(i+1,j+1,k-1) = uu(i-1,j+1,k-1)
    if (bc_neumann(mm(i,j,k),2,+1))  uu(i+1,j+1,k-1) = uu(i+1,j-1,k-1)
    if (bc_neumann(mm(i,j,k),3,-1))  uu(i+1,j+1,k-1) = uu(i+1,j+1,k+1)

    i = hi(1)
    j = lo(2)
    k = hi(3)
    if (bc_neumann(mm(i,j,k),1,+1))  uu(i+1,j-1,k+1) = uu(i-1,j-1,k+1)
    if (bc_neumann(mm(i,j,k),2,-1))  uu(i+1,j-1,k+1) = uu(i+1,j+1,k+1)
    if (bc_neumann(mm(i,j,k),3,+1))  uu(i+1,j-1,k+1) = uu(i+1,j-1,k-1)

    i = lo(1)
    j = hi(2)
    k = hi(3)
    if (bc_neumann(mm(i,j,k),1,-1))  uu(i-1,j+1,k+1) = uu(i+1,j+1,k+1)
    if (bc_neumann(mm(i,j,k),2,+1))  uu(i-1,j+1,k+1) = uu(i-1,j-1,k+1)
    if (bc_neumann(mm(i,j,k),3,+1))  uu(i-1,j+1,k+1) = uu(i-1,j+1,k-1)

    i = hi(1)
    j = hi(2)
    k = hi(3)
    if (bc_neumann(mm(i,j,k),1,+1))  uu(i+1,j+1,k+1) = uu(i-1,j+1,k+1)
    if (bc_neumann(mm(i,j,k),2,+1))  uu(i+1,j+1,k+1) = uu(i+1,j-1,k+1)
    if (bc_neumann(mm(i,j,k),3,+1))  uu(i+1,j+1,k+1) = uu(i+1,j+1,k-1)

  end subroutine impose_neumann_bcs_3d

end module impose_neumann_bcs_module
