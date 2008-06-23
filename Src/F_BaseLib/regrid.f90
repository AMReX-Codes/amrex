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


  implicit none 

  contains

    subroutine make_new_grids(mla,mla_new,mf,dx_crse,buf_wid,rr,lev,max_grid_size,new_grid)

       type(ml_layout), intent(in   ) :: mla
       type(ml_layout), intent(inout) :: mla_new
       type(multifab) , intent(in   ) :: mf
       real(dp_t)     , intent(in   ) :: dx_crse
       integer        , intent(in   ) :: buf_wid
       integer        , intent(in   ) :: max_grid_size

       ! ref_ratio now fixed across levels and dimensions -- fix me?
       integer        , intent(in   ) :: rr(:,:)

       integer, optional, intent(in   ) :: lev
       logical, optional, intent(  out) :: new_grid

       type(ml_boxarray) :: mba_new
       type(boxarray)    :: ba_new
       integer           :: llev,dim,n,i

       logical :: pmask(mf%dim)


       dim = mf%dim
       pmask = mla%pmask
       llev = 1; if (present(lev)) llev = lev

       call make_boxes(mf,ba_new,dx_crse,buf_wid,llev)
       call boxarray_maxsize(ba_new,max_grid_size)

       if (empty(ba_new)) then 
          new_grid = .false.
          call build(mba_new,llev,dim)
          do n = 1, llev-1
             mba_new%pd(n) = mla%mba%pd(n)
             mba_new%rr(n,:) = rr(n,:)
             call copy(mba_new%bas(n),mla%mba%bas(n))
          enddo
          mba_new%pd(llev) = mla%mba%pd(llev)
          call copy(mba_new%bas(llev),mla%mba%bas(llev))
          call  ml_layout_build(mla_new,mba_new,pmask)
          call destroy(mba_new)
          return
       endif

       new_grid = .true.
       
       call build(mba_new,llev+1,dim)
       do n = 1, llev
          mba_new%pd(n) = mla%mba%pd(n)
          mba_new%rr(n,:) = rr(n,:)
          call copy(mba_new%bas(n),mla%mba%bas(n))
       enddo
       call boxarray_refine(ba_new,rr(llev,:))
       mba_new%bas(llev+1) = ba_new
       mba_new%pd(llev+1) = refine(mba_new%pd(llev),rr(llev,:))

       call ml_layout_build(mla_new,mba_new,pmask)

       call destroy(mba_new)

    end subroutine make_new_grids

    subroutine make_boxes(mf,ba_new,dx_crse,buf_wid,lev)

      type(multifab), intent(in   ) :: mf
      type(boxarray), intent(  out) :: ba_new
      real(dp_t)    , intent(in   ) :: dx_crse
      integer       , intent(in   ) :: buf_wid
      integer, optional, intent(in   ) :: lev

      type(list_box_node), pointer :: bn
      type(lmultifab) :: tagboxes

      real(kind = dp_t), pointer :: sp(:,:,:,:)
      logical          , pointer :: tp(:,:,:,:)
      real(kind = dp_t) :: min_eff
      integer           :: i, j, k, dm
      integer, allocatable  :: lo(:)
      integer           :: minwidth
      integer           :: llev, ng_cell


      llev = 1; if (present(lev)) llev = lev
      dm = mf%dim
      ng_cell = mf%ng
      allocate(lo(dm))

      call lmultifab_build(tagboxes,mf%la,1,0)

      do i = 1, mf%nboxes
        if ( multifab_remote(mf, i) ) cycle
        sp => dataptr(mf, i)
        tp => dataptr(tagboxes, i)
        lo =  lwb(get_box(tagboxes, i))

        select case (dm)
          case (2)
             call tag_boxes_2d(tp(:,:,1,1),sp(:,:,1,1),lo,ng_cell,dx_crse,llev)
          case  (3)
             call tag_boxes_3d(tp(:,:,:,1),sp(:,:,:,1),lo,ng_cell,dx_crse)
          end select
      end do

      minwidth = 2
      min_eff = .7

      call cluster(ba_new, tagboxes, minwidth, buf_wid, min_eff)

      call destroy(tagboxes)
  
    end subroutine make_boxes

    subroutine tag_boxes_2d(tagbox,mf,lo,ng,dx_crse,lev)

      integer          , intent(in   ) :: lo(:),ng
      logical          , intent(  out) :: tagbox(lo(1):,lo(2):)
      real(kind = dp_t), intent(in   ) :: mf(lo(1)-ng:,lo(2)-ng:)
      real(kind = dp_t), intent(in   ) :: dx_crse
      integer, optional, intent(in   ) :: lev
      integer :: i,j,nx,ny, llev

      llev = 1; if (present(lev)) llev = lev
      nx = size(tagbox,dim=1)
      ny = size(tagbox,dim=2)

      tagbox = .false.

      select case(llev)
         case (1)
            ! tag all boxes with a density >= 1.1
            do j = lo(2),lo(2)+ny-1
               do i = lo(1),lo(1)+nx-1
                  if (mf(i,j) .gt. 1.1d0) then
                     tagbox(i,j) = .true.
                  end if
               end do
            enddo
         case (2)
            ! for level 2 tag all boxes with a density >= 2.5
            do j = lo(2),lo(2)+ny-1
               do i = lo(1),lo(1)+nx-1
                  if (mf(i,j) .gt. 1.5d0) then
                     tagbox(i,j) = .true.
                  end if
               end do
            end do
         case (3)
            ! for level 3 tag all boxes with a density >= 2.5
            do j = lo(2),lo(2)+ny-1
               do i = lo(1),lo(1)+nx-1
                  if (mf(i,j) .gt. 1.5d0) then
                     tagbox(i,j) = .true.
                  end if
               end do
            end do
         end select

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

    subroutine buffer(lev,mla,buff)

      integer,         intent(in   ) :: lev
      type(ml_layout), intent(inout) :: mla
      integer,         intent(in   ) :: buff

      type(lmultifab) :: boxes, boxes2
      type(boxarray)  :: ba_new, ba_lev
      integer         :: minwidth, n
      real(kind = dp_t) :: min_eff
      type(ml_boxarray) :: mba_new
      type(ml_layout)   :: mla_new
      type(layout)      :: la
      
      call lmultifab_build(boxes,mla%la(lev),1,0)
      call setval(boxes,.true.,all=.true.)
     
      minwidth = 2
      min_eff = .7

      call cluster(ba_new,boxes,minwidth,buff,min_eff)

      call destroy(boxes)
      call build(mba_new,lev,mla%dim)
      do n = 1, lev-1
         mba_new%pd(n) = mla%mba%pd(n)
         mba_new%rr(n,:) = mla%mba%rr(n,:)
         call copy(mba_new%bas(n),mla%mba%bas(n))
      enddo
      mba_new%pd(lev) = mla%mba%pd(lev)
      mba_new%bas(lev) = ba_new

      call ml_layout_build(mla_new,mba_new,mla%pmask)
      call destroy(mla)
      call ml_layout_build(mla,mba_new,mla_new%pmask)
      call destroy(mla_new)

      call destroy(mba_new)

    end subroutine buffer

end module regrid_module
