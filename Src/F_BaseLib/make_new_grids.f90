module make_new_grids_module

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

  implicit none 

  contains

    subroutine make_new_grids(new_grid_flag,la_crse,la_fine,mf,dx_crse,buf_wid,ref_ratio, &
                              lev,max_grid_size,aux_tag_mf)

       logical                 , intent(  out) :: new_grid_flag
       type(layout)            , intent(in   ) :: la_crse
       type(layout)            , intent(inout) :: la_fine
       type(multifab)          , intent(in   ) :: mf
       real(dp_t)              , intent(in   ) :: dx_crse
       integer                 , intent(in   ) :: buf_wid
       integer                 , intent(in   ) :: max_grid_size
       integer                 , intent(in   ) :: ref_ratio
       integer                 , intent(in   ) :: lev
       type(multifab), optional, intent(in   ) :: aux_tag_mf

       type(box)         :: pd
       type(boxarray)    :: ba_new

       if (present(aux_tag_mf)) then
          call make_boxes(mf,ba_new,dx_crse,buf_wid,lev,aux_tag_mf)
       else
          call make_boxes(mf,ba_new,dx_crse,buf_wid,lev)
       endif

       if (empty(ba_new)) then 

          new_grid_flag = .false.

       else

          new_grid_flag = .true.

          ! Need to divide max_grid_size by ref_ratio since we are creating
          !  the grids at the coarser resolution but want max_grid_size to
          !  apply at the fine resolution.
          call boxarray_maxsize(ba_new,max_grid_size/ref_ratio)
       
          call boxarray_refine(ba_new,ref_ratio)
 
          pd = refine(layout_get_pd(la_crse),ref_ratio)
   
          call layout_build_ba(la_fine,ba_new,pd,layout_get_pmask(la_crse))
   
       endif

       call destroy(ba_new)

    end subroutine make_new_grids

    subroutine make_boxes(mf,ba_new,dx_crse,buf_wid,lev,aux_tag_mf)

      use tag_boxes_module

      type(multifab)          , intent(in   ) :: mf
      type(boxarray)          , intent(  out) :: ba_new
      real(dp_t)              , intent(in   ) :: dx_crse
      integer                 , intent(in   ) :: buf_wid
      integer       , optional, intent(in   ) :: lev
      type(multifab), optional, intent(in   ) :: aux_tag_mf

      integer         :: llev
      type(lmultifab) :: tagboxes

      llev = 1; if (present(lev)) llev = lev

      call lmultifab_build(tagboxes,get_layout(mf),1,0) 
      call setval(tagboxes, .false.)

      if (present(aux_tag_mf)) then
         call tag_boxes(mf,tagboxes,dx_crse,llev,aux_tag_mf)
      else
         call tag_boxes(mf,tagboxes,dx_crse,llev)
      endif

      if (lmultifab_count(tagboxes) == 0) then

          call bl_warn('No points tagged at level ',lev)

      else 

         call cluster(ba_new, tagboxes, buf_wid)

      endif

      call destroy(tagboxes)
  
    end subroutine make_boxes

    subroutine buffer(lev,la_np1,la_nm1,ba_new,ref_ratio,buff)

      integer          , intent(in   ) :: lev
      type(layout)     , intent(in   ) :: la_np1
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
      call cluster(ba_new,tagboxes,buff_c)

      ! Now refine so we're back to the right level
      call boxarray_refine(ba_new,ref_ratio)

      call destroy(tagboxes)

      call destroy(boxes)

    end subroutine buffer

    subroutine enforce_proper_nesting(mba,la_array,max_grid_size_2,max_grid_size_3_in)

      implicit none

      type(ml_boxarray), intent(inout) :: mba
      type(layout)     , intent(inout) :: la_array(:)
      integer          , intent(in   ) :: max_grid_size_2
      integer          , intent(in   ), optional  :: max_grid_size_3_in

      integer                        :: nl,n,i,j,ng_buffer,nlevs
      integer                        :: counter
      integer                        :: ref_ratio(mba%dim)
      logical                        :: pmask(mba%dim)
      logical                        :: all_properly_nested
      type(box)                      :: pd, bx
      type(list_box)                 :: bl
      type(box_intersector), pointer :: bi(:)
      type(boxarray)                 :: ba_new
      type(boxarray)                 :: ba_newest
      type(boxarray)                 :: ba_crse_fine
      type(layout)                   :: la_old_comp
      type(boxarray)                 :: ba_old_comp
      type(lmultifab)                :: tagboxes
      type(layout)                   :: latmp
         
      integer :: max_grid_size_3

      if (present (max_grid_size_3_in)) then
         max_grid_size_3 = max_grid_size_3_in
      else
         max_grid_size_3 = max_grid_size_2
      end if

      nlevs = mba%nlevel

      pmask = get_pmask(la_array(1))
      ng_buffer = 4

      all_properly_nested = .false.
      counter = 0

      do while ( .not. all_properly_nested )

        all_properly_nested = .true.
        counter = counter + 1

        if (counter .gt. nlevs + 1) call bl_error('Still not properly nested after nlevs + 1 tries')

        nl = nlevs - 1
        do while ( (nl .ge. 2) )

            if (.not. empty(mba%bas(nl+1))) then

              ! Test whether level nl+1 boxes are properly nested inside level nl boxes.
              if (.not. ml_boxarray_properly_nested(mba, ng_buffer, pmask, nl+1, nl+1)) then

                ref_ratio = mba%rr(nl,:)
               
                ! Buffer returns a boxarray "ba_new" that contains everything at level nl 
                !  that the level nl+1 level will need for proper nesting
                call buffer(nl,la_array(nl+1),la_array(nl-1), &
                            ba_new,ref_ratio(1),ng_buffer)

                ! Merge the new boxarray "ba_new" with the existing box_array 
                ! mba%bas(nl) so that we get the union of points.
                call build(latmp,mba%bas(nl),mba%pd(nl))
                call layout_boxarray_diff(ba_old_comp,mba%pd(nl),latmp)
                call destroy(latmp)

                call build(la_old_comp,ba_old_comp,mba%pd(nl),mapping = LA_LOCAL)
                ! LA_LOCAL ==> bypass processor distribution calculation.

                ! Start to load bl with the boxes we had before in ba_old (aka mba%bas(nl)).
                do i = 1, nboxes(mba%bas(nl))
                   call push_back(bl, get_box(mba%bas(nl),i))
                end do

                ! split up ba_new so the number of intersections per
                ! box isn't to big. i.e., there is more than one box in ba_new
                if (nl .eq. 2) then
                   call boxarray_maxsize(ba_new,max_grid_size_2)
                else
                   call boxarray_maxsize(ba_new,max_grid_size_3)
                end if

                ! Now load with the new boxes that are the intersection of 
                !  ba_new with the complement of ba_old (aka mba%bas(nl))
                do j = 1, nboxes(ba_new)
                   bi => layout_get_box_intersector(la_old_comp, get_box(ba_new,j))
                   do i = 1, size(bi)
                      call push_back(bl, bi(i)%bx)
                   end do
                   deallocate(bi)
                end do

                call build(ba_newest,bl)
                call boxarray_simplify(ba_newest)

                if (nl .eq. 2) then
                   call boxarray_maxsize(ba_newest,max_grid_size_2)
                else
                   call boxarray_maxsize(ba_newest,max_grid_size_3)
                end if

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

                ! Redo the boxes using the cluster algorithm
                call lmultifab_build(tagboxes,la_array(nl-1),1,0) 
                call setval(tagboxes, .false.)
                call copy(ba_crse_fine,mba%bas(nl))
                call boxarray_coarsen(ba_crse_fine,ref_ratio)
                call setval(tagboxes, .true., ba_crse_fine)
                call destroy(ba_crse_fine)
                call cluster(ba_new, tagboxes, 0)
                call destroy(tagboxes)
                if (nl .eq. 2) then
                   call boxarray_maxsize(ba_new,max_grid_size_2/ref_ratio)
                else
                   call boxarray_maxsize(ba_new,max_grid_size_3/ref_ratio)
                end if
                call boxarray_refine(ba_new,ref_ratio)

                ! Destroy the old boxarray level and put ba_new there.
                call destroy(mba%bas(nl))
                call copy(mba%bas(nl),ba_new)

                ! Destroy the old layout and build a new one from ba_new.
                call destroy(la_array(nl))
                call layout_build_ba(la_array(nl),ba_new,mba%pd(nl),pmask)
                call destroy(ba_new)
  
                ! Double check we got the proper nesting right

                if (.not. ml_boxarray_properly_nested(mba, ng_buffer, pmask, nl+1, nl+1)) &
                  all_properly_nested = .false.

              endif  !if not properly nested
            endif  ! if fine level has any boxes

            nl = nl - 1

         enddo ! do while

      enddo ! .not. properly nested

      ! Now test on whether any grids at level nl+1 come within one level nl cell of a physical boundary,
      !  and if they do, extend those boxes to the boundary.
      do n = 2,nlevs

         pd = mba%pd(n)

         if (.not. empty(mba%bas(n))) then

            do j = 1,mba%dim
               if (.not. pmask(j)) then
                  do i = 1, nboxes(mba,n)

                     if ( (  lwb(get_box(mba%bas(n),i),j) - lwb(pd,j)) .le. 2) then
                        bx = get_box(mba%bas(n),i)
                        call set_lwb(bx,j,lwb(pd,j))
                        call set_box(mba%bas(n),i,bx)
                     end if

                     if ( (upb(pd,j) - upb(get_box(mba%bas(n),i),j)) .le. 2) then
                        bx = get_box(mba%bas(n),i)
                        call set_upb(bx,j,upb(pd,j))
                        call set_box(mba%bas(n),i,bx)
                     end if

                  end do
               end if
            end do

         end if

      end do

    end subroutine enforce_proper_nesting

end module make_new_grids_module
