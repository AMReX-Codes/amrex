module regrid_module

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

  contains

    subroutine make_new_grids(mla,mla_new,mf,dx_crse,buf_wid)

       type(ml_layout), intent(in   ) :: mla
       type(ml_layout), intent(  out) :: mla_new
       type(multifab) , intent(in   ) :: mf
       real(dp_t)     , intent(in   ) :: dx_crse
       integer        , intent(in   ) :: buf_wid

       type(ml_boxarray) :: mba_new
       type(boxarray) :: ba_new
       type(boxarray) :: ba_loc
       type(box     ) :: bx_loc

       integer :: dim
       logical :: pmask(mf%dim)

       pmask = layout_get_pmask(mf%la)

       call make_boxes(mf,ba_new,dx_crse,buf_wid)
       call boxarray_refine(ba_new,mla%mba%rr(1,:))

       call copy(mba_new,mla%mba)
       mba_new%bas(2) = ba_new

       dim = mf%dim
       call ml_layout_build(mla_new,mba_new,pmask)

    end subroutine make_new_grids

    subroutine make_boxes(mf,ba_new,dx_crse,buf_wid)

      type(multifab), intent(in   ) :: mf
      type(boxarray), intent(  out) :: ba_new
      real(dp_t)    , intent(in   ) :: dx_crse
      integer       , intent(in   ) :: buf_wid

      type(list_box_node), pointer :: bn
      type(lmultifab) :: tagboxes

      real(kind = dp_t), pointer :: sp(:,:,:,:)
      logical          , pointer :: tp(:,:,:,:)

      real(kind = dp_t) :: min_eff
      integer           :: i, j, k, dm
      integer, allocatable  :: lo(:)
      integer           :: minwidth

      dm = mf%dim
      allocate(lo(dm))

      call lmultifab_build(tagboxes,mf%la,1,0)

      do i = 1, mf%nboxes
        if ( multifab_remote(mf, i) ) cycle
        sp => dataptr(mf, i)
        tp => dataptr(tagboxes, i)
        lo =  lwb(get_box(tagboxes, i))
        select case (dm)
          case (2)
             call tag_boxes_2d(tp(:,:,1,1),sp(:,:,1,1),lo,ng_cell,dx_crse)
          case (3)
             call tag_boxes_3d(tp(:,:,:,1),sp(:,:,:,1),lo,ng_cell,dx_crse)
          end select
      end do

      minwidth = 2
      min_eff = .7

      call cluster(ba_new, tagboxes, minwidth, buf_wid, min_eff)

      call destroy(tagboxes)

    end subroutine make_boxes

    subroutine tag_boxes_2d(tagbox,mf,lo,ng,dx_crse)

      integer          , intent(in   ) :: lo(:),ng
      logical          , intent(  out) :: tagbox(lo(1):,lo(2):)
      real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:)
      real(kind = dp_t), intent(in   ) :: dx_crse

      integer :: i,j,nx,ny

      nx = size(tagbox,dim=1)
      ny = size(tagbox,dim=2)

      tagbox = .false.

      ! NOTE: THIS IS JUST A PLACEHOLDER LOOP -- CRITERION SHOULD REALLY BE USER-SET
      do j = lo(2),lo(2)+ny-1
         do i = lo(1),lo(1)+nx-1
            if (mf(i,j) .gt. 0.0) then
               tagbox(i,j) = .true.
            end if
         end do
      end do

    end subroutine tag_boxes_2d

    subroutine tag_boxes_3d(tagbox,mf,lo,ng,dx_crse)

      integer          , intent(in   ) :: lo(:),ng
      logical          , intent(  out) :: tagbox(lo(1):,lo(2):,lo(3):)
      real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
      real(kind = dp_t), intent(in   ) :: dx_crse

      integer :: i,j,k,nx,ny,nz

      nx = size(tagbox,dim=1)
      ny = size(tagbox,dim=2)
      nz = size(tagbox,dim=3)

      tagbox = .false.

      ! NOTE: THIS IS JUST A PLACEHOLDER ARRAY -- CRITERION SHOULD REALLY BE USER-SET
      do k = lo(3),lo(3)+nz-1
      do j = lo(2),lo(2)+ny-1
      do i = lo(1),lo(1)+nx-1
         if (mf(i,j,k) .gt. 0.0) then
            tagbox(i,j,k) = .true.
         end if
      end do
      end do
      end do

    end subroutine tag_boxes_3d

end module regrid_module
