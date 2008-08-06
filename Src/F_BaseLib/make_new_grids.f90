module make_new_grids_module

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

  implicit none 

  integer        , parameter, private :: minwidth = 4
  real(kind=dp_t), parameter, private :: min_eff  = .7

  contains

    subroutine make_new_grids(la_crse,la_fine,mf,dx_crse,buf_wid,ref_ratio,lev_in,max_grid_size,new_grid)

       type(layout)     , intent(in   ) :: la_crse
       type(layout)     , intent(inout) :: la_fine
       type(multifab)   , intent(in   ) :: mf
       real(dp_t)       , intent(in   ) :: dx_crse
       integer          , intent(in   ) :: buf_wid
       integer          , intent(in   ) :: max_grid_size
       integer          , intent(in   ) :: ref_ratio
       integer, optional, intent(in   ) :: lev_in
       logical, optional, intent(  out) :: new_grid

       type(box)         :: pd_fine
       type(boxarray)    :: ba_new
       integer           :: llev

       llev = 1; if (present(lev_in)) llev = lev_in

       call make_boxes(mf,ba_new,dx_crse,buf_wid,llev)

       if (empty(ba_new)) then 

          new_grid = .false.

       else

          new_grid = .true.

          ! Need to divide max_grid_size by ref_ratio since we are creating
          !  the grids at the coarser resolution but want max_grid_size to
          !  apply at the fine resolution.
          call boxarray_maxsize(ba_new,max_grid_size/ref_ratio)
       
          call boxarray_refine(ba_new,ref_ratio)
 
          pd_fine = refine(layout_get_pd(la_crse),ref_ratio)
   
          call layout_build_ba(la_fine,ba_new,pd_fine,layout_get_pmask(la_crse))
   
       endif

       call destroy(ba_new)

    end subroutine make_new_grids

    subroutine make_boxes(mf,ba_new,dx_crse,buf_wid,lev)

      use tag_boxes_module

      type(multifab), intent(in   ) :: mf
      type(boxarray), intent(  out) :: ba_new
      real(dp_t)    , intent(in   ) :: dx_crse
      integer       , intent(in   ) :: buf_wid
      integer, optional, intent(in   ) :: lev

      integer         :: llev
      type(lmultifab) :: tagboxes

      llev = 1; if (present(lev)) llev = lev

      call lmultifab_build(tagboxes,mf%la,1,0) 
      call setval(tagboxes, .false.)

      call tag_boxes(mf,tagboxes,llev)

      if (lmultifab_count(tagboxes) == 0) &
          call bl_warn('No points tagged at level ',lev)

      call cluster(ba_new, tagboxes, minwidth, buf_wid, min_eff)

      call destroy(tagboxes)
  
    end subroutine make_boxes

    subroutine buffer(lev,la_np1,la_n,la_nm1,ba_new,ref_ratio,buff)

      integer          , intent(in   ) :: lev
      type(layout)     , intent(in   ) :: la_np1
      type(layout)     , intent(in   ) :: la_n
      type(layout)     , intent(in   ) :: la_nm1
      type(boxarray)   , intent(inout) :: ba_new
      integer          , intent(in   ) :: buff
      integer          , intent(in   ) :: ref_ratio

      type(boxarray)  :: boxes
      type(lmultifab) :: tagboxes
      integer         :: buff_c

      ! Coarsen the (n+1) boxes twice so that we define tagboxes on the (n-1) level.
      call copy(boxes,get_boxarray(la_np1))
      call boxarray_coarsen(boxes,ref_ratio*ref_ratio)

      call lmultifab_build(tagboxes,la_nm1,1,0)
      call setval(tagboxes, .false., all=.true.)
      call setval(tagboxes, .true., boxes)

      if (lmultifab_count(tagboxes) == 0) &
          call bl_warn('No points tagged at level ',lev)

      buff_c = (buff+1) / (ref_ratio*ref_ratio)
      call cluster(ba_new,tagboxes,minwidth,buff_c,min_eff)

      ! Now refine so we're back to the right level
      call boxarray_refine(ba_new,ref_ratio)

      call destroy(tagboxes)

      call destroy(boxes)

    end subroutine buffer

    subroutine enforce_proper_nesting(mba,la_array)

     use probin_module, only : nlevs, pmask, ref_ratio, max_grid_size

     type(ml_boxarray), intent(inout) :: mba
     type(layout)     , intent(inout) :: la_array(:)

     integer :: nl,i,ii,jj,ng_buffer
     type(list_box)                 :: bl
     type(box_intersector), pointer :: bi(:)
     type(boxarray)                 :: ba_new
     type(boxarray)                 :: ba_newest
     type(layout)                   :: la_old_comp
     type(boxarray)                 :: ba_old_comp
        
      nl = nlevs - 1

      ng_buffer = 8

      do while ( (nl .ge. 2) )

            if (.not. ml_boxarray_properly_nested(mba, ng_buffer, pmask, nl+1, nl+1)) then
               
                ! Buffer returns a boxarray "ba_new" that contains everything at level nl 
                !  that the level nl+1 level will need for proper nesting
                call buffer(nl,la_array(nl+1),la_array(nl),la_array(nl-1), &
                            ba_new,ref_ratio,ng_buffer)

                ! Merge the new boxarray "ba_new" with the existing box_array 
                ! mba%bas(nl) so that we get the union of points.
                call boxarray_complementIn(ba_old_comp,mba%pd(nl),mba%bas(nl))
                call build(la_old_comp,ba_old_comp,mba%pd(nl))

                ! Start to load bl with the boxes we had before in ba_old (aka mba%bas(nl)).
                do i = 1, mba%bas(nl)%nboxes
                   call push_back(bl,  mba%bas(nl)%bxs(i))
                end do

                ! Now load with the new boxes that are the intersection of 
                !  ba_new with the complement of ba_old (aka mba%bas(nl))
                do jj = 1, ba_new%nboxes
                   bi => layout_get_box_intersector(la_old_comp, ba_new%bxs(jj))
                   do ii = 1, size(bi)
                      call push_back(bl, bi(ii)%bx)
                   end do
                   deallocate(bi)
                end do

                call build(ba_newest,bl)
                call boxarray_simplify(ba_newest)
                call boxarray_maxsize(ba_newest,max_grid_size)

                ! Do some cleanup.
                call destroy(bl)
                call destroy(ba_new)
                call destroy(ba_old_comp)
                call destroy(la_old_comp)

                ! Replace mba%bas(nl) by ba_newest
                call destroy(mba%bas(nl))
                call copy(mba%bas(nl),ba_newest)
                call destroy(ba_newest)

                ! Destroy the old layout and build a new one.
                call destroy(la_array(nl))
                call layout_build_ba(la_array(nl),mba%bas(nl),mba%pd(nl),pmask)
  
                ! Double check we got the proper nesting right
                if (.not. ml_boxarray_properly_nested(mba, ng_buffer, pmask, nl+1, nl+1)) &
                  call bl_error('Still not properly nested, darn it')

            endif  !if not properly nested

            nl = nl - 1

      enddo ! do while

    end subroutine enforce_proper_nesting

end module make_new_grids_module
