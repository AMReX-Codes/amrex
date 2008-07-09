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
  use ml_layout_module

  implicit none 

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
       integer           :: llev,dm,n,i

       logical :: pmask(mf%dim)

       dm = mf%dim
       pmask = layout_get_pmask(la_crse)
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
   
          call layout_build_ba(la_fine,ba_new,pd_fine,pmask)
   
          call destroy(ba_new)

       endif

    end subroutine make_new_grids

    subroutine make_boxes(mf,ba_new,dx_crse,buf_wid,lev)

      use tag_boxes_module

      type(multifab), intent(in   ) :: mf
      type(boxarray), intent(  out) :: ba_new
      real(dp_t)    , intent(in   ) :: dx_crse
      integer       , intent(in   ) :: buf_wid
      integer, optional, intent(in   ) :: lev

      type(list_box_node), pointer :: bn
      type(lmultifab) :: tagboxes

      real(kind = dp_t) :: min_eff
      integer           :: i, j, k, dm
      integer           :: minwidth
      integer           :: llev, ng_cell

      llev = 1; if (present(lev)) llev = lev
      dm = mf%dim
      ng_cell = mf%ng

      call lmultifab_build(tagboxes,mf%la,1,0)

      call setval(tagboxes,.false.)

      call tag_boxes(mf,tagboxes,llev)

      minwidth = 2
      min_eff = .7

      call cluster(ba_new, tagboxes, minwidth, buf_wid, min_eff)

      call destroy(tagboxes)
  
    end subroutine make_boxes

    subroutine buffer(lev,la_fine,ba_crse,ref_ratio,buff)

      integer          , intent(in   ) :: lev
      type(layout)     , intent(in   ) :: la_fine
      type(boxarray)   , intent(inout) :: ba_crse
      integer          , intent(in   ) :: buff
      integer          , intent(in   ) :: ref_ratio

      type(lmultifab) :: boxes
      integer         :: minwidth, n
      real(kind = dp_t) :: min_eff
      
      call lmultifab_build(boxes,la_fine,1,0)

      call setval(boxes,.false.)

      min_eff  = .7     
      minwidth = 2

      call cluster(ba_crse,boxes,minwidth,buff,min_eff)

      call boxarray_coarsen(ba_crse, ref_ratio)

    end subroutine buffer

end module make_new_grids_module
