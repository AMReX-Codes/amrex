module nodal_enforce_module

  use bl_constants_module
  use bc_functions_module
  use mg_tower_module

  implicit none

contains

  subroutine enforce_dirichlet_rhs(divu_rhs,mgt)

    type(multifab) , intent(inout) :: divu_rhs(:)
    type(mg_tower) , intent(in   ) :: mgt(:)

    type(imultifab):: mask
    integer        :: i,n,ng_d
    integer        :: lo(get_dim(divu_rhs(1))),hi(get_dim(divu_rhs(1)))
    integer        , pointer ::   mp(:,:,:,:) 
    real(kind=dp_t), pointer :: divp(:,:,:,:) 

    ng_d = nghost(divu_rhs(1))

    do n = 1, size(mgt)
       mask = mgt(n)%mm(mgt(n)%nlevels)
       do i = 1, nfabs(divu_rhs(n))
          divp => dataptr(divu_rhs(n) ,i)
            mp => dataptr(mask        ,i)
          lo = lwb(get_box(divu_rhs(n),i))
          hi = upb(get_box(divu_rhs(n),i))
          select case (mgt(n)%dim)
          case (1)
             call enforce_dirichlet_1d(divp(:,1,1,1), mp(:,1,1,1), &
                  lo, hi, ng_d)
          case (2)
             call enforce_dirichlet_2d(divp(:,:,1,1), mp(:,:,1,1), &
                  lo, hi, ng_d)
          case (3)
             call enforce_dirichlet_3d(divp(:,:,:,1), mp(:,:,:,1), &
                  lo, hi, ng_d)
          end select
       end do
    end do

  end subroutine enforce_dirichlet_rhs

  ! ******************************************************************************** !

  subroutine enforce_dirichlet_1d(divu_rhs,mm,lo,hi,ng_d)

    integer        , intent(in   ) :: ng_d,lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: divu_rhs(lo(1)-ng_d:)
    integer        , intent(in   ) ::       mm(lo(1)     :)

    integer        :: i

    ! This takes care of any fine nodes on the coarse-fine interface
    !      OR                           on the outflow boundary
    do i = lo(1), hi(1)+1
        if (bc_dirichlet(mm(i),1,0)) divu_rhs(i) = ZERO
    end do

  end subroutine enforce_dirichlet_1d


  ! ******************************************************************************** !

  subroutine enforce_dirichlet_2d(divu_rhs,mm,lo,hi,ng_d)

    integer        , intent(in   ) :: ng_d,lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: divu_rhs(lo(1)-ng_d:,lo(2)-ng_d:)
    integer        , intent(in   ) ::       mm(lo(1)     :,lo(2)     :)

    integer        :: i,j

    ! This takes care of any fine nodes on the coarse-fine interface
    !      OR                           on the outflow boundary
    do j = lo(2), hi(2)+1
    do i = lo(1), hi(1)+1
        if (bc_dirichlet(mm(i,j),1,0)) divu_rhs(i,j) = ZERO
    end do
    end do

  end subroutine enforce_dirichlet_2d

  ! ******************************************************************************** !

  subroutine enforce_dirichlet_3d(divu_rhs,mm,lo,hi,ng_d)

    integer        , intent(in   ) :: ng_d,lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: divu_rhs(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:)
    integer        , intent(in   ) ::       mm(lo(1)     :,lo(2)     :,lo(3)     :)

    integer        :: i,j,k

    ! This takes care of any fine nodes on the coarse-fine interface
    !      OR                           on the outflow boundary
    do k = lo(3), hi(3)+1
    do j = lo(2), hi(2)+1
    do i = lo(1), hi(1)+1
        if (bc_dirichlet(mm(i,j,k),1,0)) divu_rhs(i,j,k) = ZERO
    end do
    end do
    end do

  end subroutine enforce_dirichlet_3d

end module nodal_enforce_module


