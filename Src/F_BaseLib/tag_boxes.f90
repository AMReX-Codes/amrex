module tag_boxes_module

  use BoxLib
  use omp_module
  use f2kcli
  use list_box_module
  use boxarray_module
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use box_util_module
  use bl_IO_module
  use cluster_module
  use ml_layout_module

  implicit none 

contains

  subroutine tag_boxes(mf,tagboxes,dx,lev)

    type( multifab), intent(in   ) :: mf
    type(lmultifab), intent(inout) :: tagboxes
    real(dp_t)     , intent(in   ) :: dx
    integer        , intent(in   ) :: lev

    real(kind = dp_t), pointer :: sp(:,:,:,:)
    logical          , pointer :: tp(:,:,:,:)
    integer           :: i, lo(mf%dim), hi(mf%dim), ng

    ng = mf%ng

    do i = 1, mf%nboxes
       if ( multifab_remote(mf, i) ) cycle
       sp => dataptr(mf, i)
       tp => dataptr(tagboxes, i)
       lo =  lwb(get_box(tagboxes, i))
       hi =  upb(get_box(tagboxes, i))
       select case (mf%dim)
       case (2)
          call tag_boxes_2d(tp(:,:,1,1),sp(:,:,1,1),lo,hi,ng,dx,lev)
       case  (3)
          call tag_boxes_3d(tp(:,:,:,1),sp(:,:,:,1),lo,hi,ng,dx,lev)
       end select
    end do

  end subroutine tag_boxes

  subroutine tag_boxes_2d(tagbox,mf,lo,hi,ng,dx,lev)

    integer          , intent(in   ) :: lo(:),hi(:),ng
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):)
    real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:)
    real(dp_t)       , intent(in   ) :: dx
    integer, optional, intent(in   ) :: lev
    integer :: i,j,llev

    llev = 1; if (present(lev)) llev = lev

    tagbox = .false.

    select case(llev)
    case (1)
       ! tag all boxes with a density >= 1.01
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (mf(i,j) .gt. 1.01d0) then
                tagbox(i,j) = .true.
             end if
          end do
       enddo
    case (2)
       ! for level 2 tag all boxes with a density >= 1.1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             if (mf(i,j) .gt. 1.1d0) then
                tagbox(i,j) = .true.
             end if
          end do
       end do
    case default
       ! for level 3 or greater tag all boxes with a density >= 1.5
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
    logical          , intent(  out) :: tagbox(lo(1):,lo(2):,lo(3):)
    real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(dp_t)       , intent(in   ) :: dx
    integer, optional, intent(in   ) :: lev

    integer :: i,j,k,llev

    llev = 1; if (present(lev)) llev = lev

    tagbox = .false.

    select case(llev)
    case (1)
       ! tag all boxes with a density >= 1.01
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
       ! for level 2 tag all boxes with a density >= 1.1
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
       ! for level 3 or greater tag all boxes with a density >= 1.5
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
