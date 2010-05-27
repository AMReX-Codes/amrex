module tag_boxes_module

  use BoxLib
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

  subroutine tagboxes_coarsen(tagboxes,ctagboxes,ratio)

    type(lmultifab), intent(in   ) :: tagboxes
    type(lmultifab), intent(inout) :: ctagboxes
    integer,         intent(in   ) :: ratio

    integer          :: ii, i, j, k, ic, jc, kc
    integer          :: flo(tagboxes%dim), fhi(tagboxes%dim)
    type(layout)     :: cla
    type(boxarray)   :: cba
    logical, pointer :: fp(:,:,:,:), cp(:,:,:,:)

    call bl_assert(ncomp(tagboxes) == 1, 'tagboxes should only have one component')
    !
    ! ctagboxes should be an lmultifab on which build() has not been called.
    ! ctagboxes will be built on the appropriately grown & coarsen'd boxarray
    ! and will have no grow cells itself.  We want all the coarsen'd values
    ! to be in ctagboxes valid region.  Callers of this routine need to
    ! destroy both ctagboxes and its layout.
    !
    call boxarray_build_copy(cba, get_boxarray(tagboxes))

    call boxarray_grow(cba, nghost(tagboxes))

    call boxarray_coarsen(cba, ratio)
    !
    ! I'm playing a little fast & loose here.
    ! I'm assuming we don't really care about the domain or the periodicity.
    ! All we really need to get right is the mapping.
    !
    call build(cla, cba, get_pd(tagboxes%la), get_pmask(tagboxes%la), explicit_mapping = get_proc(tagboxes%la))

    call destroy(cba)

    call build(ctagboxes, cla, 1, 0)

    call setval(ctagboxes, .false., all = .true.)

    do ii = 1, tagboxes%nboxes
       if ( remote(tagboxes, ii) ) cycle

       fp  => dataptr(tagboxes,  ii)
       cp  => dataptr(ctagboxes, ii)

       flo = lwb(get_pbox(tagboxes, ii))
       fhi = upb(get_pbox(tagboxes, ii))

       select case (tagboxes%dim)
       case (2)
          do j = flo(2), fhi(2)
             jc = int_coarsen(j,ratio)
             do i = flo(1), fhi(1)
                ic = int_coarsen(i,ratio)
                if ( fp(i,j,1,1) ) cp(ic,jc,1,1) = .true.
             end do
          end do
       case  (3)
          do k = flo(3), fhi(3)
             kc = int_coarsen(k,ratio)
             do j = flo(2), fhi(2)
                jc = int_coarsen(j,ratio)
                do i = flo(1), fhi(1)
                   ic = int_coarsen(i,ratio)
                   if ( fp(i,j,k,1) ) cp(ic,jc,kc,1) = .true.
                end do
             end do
          end do
       end select
    end do

  end subroutine tagboxes_coarsen

end module tag_boxes_module
