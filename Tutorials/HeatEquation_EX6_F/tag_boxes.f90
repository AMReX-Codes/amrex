module tag_boxes_module

  use multifab_module
  use bl_error_module

  implicit none 

contains

  subroutine tag_boxes(mf,tagboxes,dx,lev,aux_tag_mf)

    type( multifab)         , intent(in   ) :: mf
    type(lmultifab)         , intent(inout) :: tagboxes
    real(dp_t)              , intent(in   ) :: dx
    integer                 , intent(in   ) :: lev
    type(multifab), optional, intent(in   ) :: aux_tag_mf
    ! aux_tag_mf allows user to pass in additional multifabs for tagging logic
    
    ! local variables
    real(kind = dp_t), pointer :: mfp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(get_dim(mf)), hi(get_dim(mf)), ng

    if (present(aux_tag_mf)) then
       call bl_error("tag_boxes.f90: aux_tag_mf passed to tag_boxes without implementation")
    end if

    ng = nghost(mf)

    do i = 1, nfabs(mf)
       mfp => dataptr(mf, i)
       tp  => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))
       hi =  upb(get_box(tagboxes, i))
       select case (get_dim(mf))
       case (2)
          call tag_boxes_2d(tp(:,:,1,1),mfp(:,:,1,1),lo,hi,ng,dx,lev)
       case  (3)
          call tag_boxes_3d(tp(:,:,:,1),mfp(:,:,:,1),lo,hi,ng,dx,lev)
       end select
    end do

  end subroutine tag_boxes

  subroutine tag_boxes_2d(tagbox,mf,lo,hi,ng,dx,lev)

    integer          , intent(in   ) :: lo(:),hi(:),ng
    logical          , intent(  out) :: tagbox(lo(1)   :,lo(2)   :)
    real(kind = dp_t), intent(in   ) ::     mf(lo(1)-ng:,lo(2)-ng:)
    real(dp_t)       , intent(in   ) :: dx
    integer          , intent(in   ) :: lev

    ! local variables
    integer :: i,j

    ! initially say that we do not want to tag any cells for refinement
    tagbox = .false.

    select case(lev)
    case (1)
       ! tag all boxes where the first component of mf >= 1.01
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (mf(i,j) .gt. 1.01d0) then
                tagbox(i,j) = .true.
             end if
          end do
       enddo
    case (2)
       ! for level 2 tag all boxes where the first component of mf >= 1.1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (mf(i,j) .gt. 1.1d0) then
                tagbox(i,j) = .true.
             end if
          end do
       end do
    case default
       ! for level 3 or greater tag all boxes where the first component of mf >= 1.5
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (mf(i,j) .gt. 1.5d0) then
                tagbox(i,j) = .true.
             end if
          end do
       end do
    end select

  end subroutine tag_boxes_2d

  subroutine tag_boxes_3d(tagbox,mf,lo,hi,ng,dx,lev)

    integer          , intent(in   ) :: lo(:),hi(:),ng
    logical          , intent(  out) :: tagbox(lo(1)   :,lo(2)   :,lo(3)   :)
    real(kind = dp_t), intent(in   ) ::     mf(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(dp_t)       , intent(in   ) :: dx
    integer          , intent(in   ) :: lev

    ! local variables
    integer :: i,j,k

    ! initially say that we do not want to tag any cells for refinement
    tagbox = .false.

    select case(lev)
    case (1)
       ! tag all boxes where the first component of mf >= 1.01
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j,k) .gt. 1.01d0) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          enddo
       end do
    case (2)
       ! for level 2 tag all boxes where the first component of mf >= 1.1
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j,k) .gt. 1.1d0) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          end do
       end do
    case default
       ! for level 3 or greater tag all boxes where the first component of mf >= 1.5
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (mf(i,j,k) .gt. 1.5d0) then
                   tagbox(i,j,k) = .true.
                end if
             end do
          end do
       end do
    end select

  end subroutine tag_boxes_3d

end module tag_boxes_module
