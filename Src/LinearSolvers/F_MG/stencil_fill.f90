module stencil_fill_module

  use bl_types
  use bc_module
  use multifab_module
  use stencil_module

  implicit none

  real (kind = dp_t), private, parameter :: ZERO  = 0.0_dp_t
  real (kind = dp_t), private, parameter :: ONE   = 1.0_dp_t
  real (kind = dp_t), private, parameter :: TWO   = 2.0_dp_t
  real (kind = dp_t), private, parameter :: THREE = 3.0_dp_t
  real (kind = dp_t), private, parameter :: FOUR  = 4.0_dp_t
  real (kind = dp_t), private, parameter :: FIVE  = 5.0_dp_t
  real (kind = dp_t), private, parameter :: SIX   = 6.0_dp_t
  real (kind = dp_t), private, parameter :: SEVEN = 7.0_dp_t
  real (kind = dp_t), private, parameter :: EIGHT = 8.0_dp_t
  real (kind = dp_t), private, parameter :: TEN   = 10.0_dp_t
  real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t
  real (kind = dp_t), private, parameter :: THIRD = 1.0_dp_t/3.0_dp_t
  real (kind = dp_t), private, parameter :: FOUR_THIRD = 4.0_dp_t/3.0_dp_t

contains

  recursive subroutine stencil_fill_nodal_all_mglevels(mgt, sg, stencil_type)

    use coarsen_coeffs_module
    use mg_module

    type(mg_tower ), intent(inout) :: mgt
    type(multifab ), intent(inout) :: sg(:)
    integer        , intent(in   ) :: stencil_type

    type(multifab)                 :: stored_coeffs, stored_coeffs_grown
    type(multifab)                 :: new_coeffs_grown
    type(multifab), allocatable    :: coarse_coeffs(:)
    type(boxarray)                 :: ba_cc
    type(  layout)                 :: old_la_grown, new_la_grown
    integer                        :: i, maxlev, maxlev_bottom
    real(dp_t), pointer            :: sc_orig(:,:,:,:), sc_grown(:,:,:,:)

    maxlev = mgt%nlevels

    ! NOTE: sg(maxlev) comes in built and filled, but the other levels
    !       are not even built yet
    do i = maxlev-1, 1, -1
       call multifab_build(sg(i), mgt%ss(i)%la, 1, 1)
       call setval(sg(i), ZERO, 1, 1, all=.true.)
       call coarsen_cell_coeffs(sg(i+1),sg(i))
       call multifab_fill_boundary(sg(i))
    end do

    do i = maxlev, 1, -1
       call stencil_fill_nodal(mgt%ss(i), sg(i), mgt%dh(:,i), mgt%mm(i), &
                               mgt%face_type, stencil_type)
    end do

    if (associated(mgt%bottom_mgt)) then

       call multifab_build(stored_coeffs, mgt%ss(1)%la, 1, 1)
       call multifab_copy_c(stored_coeffs,1,sg(1),1,1,ng=sg(1)%ng)
       call multifab_fill_boundary(stored_coeffs)

       maxlev_bottom = mgt%bottom_mgt%nlevels
       allocate(coarse_coeffs(maxlev_bottom))
       call multifab_build(coarse_coeffs(maxlev_bottom),mgt%bottom_mgt%cc(maxlev_bottom)%la,1,1)
       call setval(coarse_coeffs(maxlev_bottom),ZERO,all=.true.)

       ! Grow the stored coefficients
       call boxarray_build_copy(ba_cc,get_boxarray(stored_coeffs))
       call boxarray_grow(ba_cc,1)
       call layout_build_ba(old_la_grown,ba_cc,pmask=mgt%ss(1)%la%lap%pmask, &
                            explicit_mapping=get_proc(mgt%ss(1)%la))
       call destroy(ba_cc)
       call multifab_build(stored_coeffs_grown,old_la_grown,1,ng=0)

       do i = 1, stored_coeffs_grown%nboxes
          if (remote(stored_coeffs_grown,i)) cycle
          sc_orig  => dataptr(stored_coeffs      ,i,get_pbox(stored_coeffs_grown,i),1,1)
          sc_grown => dataptr(stored_coeffs_grown,i,get_pbox(stored_coeffs_grown,i),1,1)
          sc_grown = sc_orig
       end do

       call boxarray_build_copy(ba_cc,get_boxarray(mgt%bottom_mgt%ss(maxlev_bottom)))
       call boxarray_grow(ba_cc,1)
       call layout_build_ba(new_la_grown,ba_cc,pmask=mgt%ss(1)%la%lap%pmask, &
                            explicit_mapping=get_proc(mgt%bottom_mgt%ss(maxlev_bottom)%la))
       call destroy(ba_cc)
       call multifab_build(new_coeffs_grown,new_la_grown,1,ng=0)
       call multifab_copy_c(new_coeffs_grown,1,stored_coeffs_grown,1,1)

       do i = 1, new_coeffs_grown%nboxes
          if (remote(new_coeffs_grown,i)) cycle
          sc_orig  => dataptr(coarse_coeffs(maxlev_bottom),i,get_pbox(new_coeffs_grown,i),1,1)
          sc_grown => dataptr(new_coeffs_grown    ,i,get_pbox(new_coeffs_grown,i),1,1)
          sc_orig = sc_grown
       end do

       call destroy(new_coeffs_grown)

       call stencil_fill_nodal_all_mglevels(mgt%bottom_mgt, coarse_coeffs, stencil_type)

       call destroy(coarse_coeffs(maxlev_bottom))
       deallocate(coarse_coeffs)

       call destroy(stored_coeffs)
       call destroy(stored_coeffs_grown)
       call destroy(old_la_grown)
       call destroy(new_la_grown)

    end if

    do i = maxlev-1, 1, -1
       call destroy(sg(i))
    end do

  end subroutine stencil_fill_nodal_all_mglevels

  subroutine stencil_fill_nodal(ss, sg, dh, mask, face_type, stencil_type)

    use stencil_nodal_module 

    type(multifab ), intent(inout) :: ss
    type(multifab ), intent(inout) :: sg
    real(kind=dp_t), intent(in   ) :: dh(:)
    type(imultifab), intent(inout) :: mask
    integer        , intent(in   ) :: face_type(:,:,:)
    integer        , intent(in   ) :: stencil_type

    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)

    type(box)                :: pd_periodic, bx, nbx, bx1
    type(boxarray)           :: bxa_periodic, bxa_temp
    integer                  :: i, ib, jb, kb, ib_lo, jb_lo, kb_lo
    integer                  :: shift_vect(ss%dim)
    type(list_box)           :: lb,nbxs
    type(box), allocatable   :: bxs(:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "stencil_fill_nodal")
    !
    ! Do this just to set everything in the mask to zero.
    !
    call setval(mask,BC_INT)
    !
    ! Construct a new boxarray that has periodically translated boxes as well
    ! as the original boxes.
    !
    pd_periodic = ss%la%lap%pd

    call boxarray_build_copy(bxa_periodic,ss%la%lap%bxa)

    if ( any(ss%la%lap%pmask) ) then
       !
       ! First trim out all boxes that can't effect periodicity.
       !
       do i = 1,nboxes(bxa_periodic)
          if ( .not. contains(pd_periodic, get_box(bxa_periodic,i), strict = .true.) ) then
             call push_back(lb, get_box(bxa_periodic,i))
          end if
       end do

       do i = 1,ss%dim
          if (ss%la%lap%pmask(i)) then
             pd_periodic = grow(grow(pd_periodic,1,i,-1),1,i,1)
          end if
       end do

       ib_lo = 1
       if ( ss%la%lap%pmask(1) )    ib_lo = -1

       jb_lo = 1
       if ( ss%dim .ge. 2) then
          if ( ss%la%lap%pmask(2) ) jb_lo = -1
       end if

       kb_lo = 1
       if ( ss%dim .ge. 3) then
          if ( ss%la%lap%pmask(3) ) kb_lo = -1
       end if

       do kb = kb_lo, 1
          do jb = jb_lo, 1
             do ib = ib_lo, 1
                call copy(bxa_temp,lb)

                shift_vect = 0

                if ( ss%la%lap%pmask(1) )    shift_vect(1) = ib * extent(ss%la%lap%pd,1)

                if ( ss%dim > 1 ) then
                   if ( ss%la%lap%pmask(2) ) shift_vect(2) = jb * extent(ss%la%lap%pd,2)
                end if
                if ( ss%dim > 2 ) then
                   if ( ss%la%lap%pmask(3) ) shift_vect(3) = kb * extent(ss%la%lap%pd,3)
                end if

                call boxarray_shift(bxa_temp,shift_vect)

                do i = 1, bxa_temp%nboxes
                   bx1 = intersection(bxa_temp%bxs(i),pd_periodic)
                   if ( .not. empty(bx1) ) then
                      call push_back(nbxs, bx1)
                   end if
                end do

                call destroy(bxa_temp)
             end do
          end do
       end do

       call destroy(lb)

       allocate(bxs(size(nbxs)))

       do i = 1, size(bxs)
          bxs(i) = front(nbxs)
          call pop_front(nbxs)
       end do

       call boxarray_add_clean_boxes(bxa_periodic,bxs,simplify = .false.)

       call destroy(nbxs)

    end if

    do i = 1, ss%nboxes
       if ( multifab_remote(ss,i) ) cycle

       sp => dataptr(ss,   i)
       cp => dataptr(sg,   i)
       mp => dataptr(mask, i)

       bx  = get_box(ss,i)
       nbx = get_ibox(ss, i)
       call stencil_set_bc_nodal(ss%dim, bx, nbx, i, mask, face_type, pd_periodic, bxa_periodic)

       select case (ss%dim)
       case (1)
          call s_simple_1d_nodal(sp(:,1,1,:), cp(:,1,1,1), mp(:,1,1,1), dh)
       case (2)
          if (stencil_type == ST_DENSE) then
            call s_dense_2d_nodal(sp(:,:,1,:), cp(:,:,1,1), mp(:,:,1,1), &
                                  face_type(i,:,:), dh)
          else if (stencil_type == ST_CROSS) then
            call s_cross_2d_nodal(sp(:,:,1,:), cp(:,:,1,1), mp(:,:,1,1), &
                                   face_type(i,:,:), dh)
          else 
            print *,'DONT KNOW THIS NODAL STENCIL TYPE ',stencil_type
            call bl_error('stencil_fill_nodal')
          end if
       case (3)
          if (stencil_type == ST_DENSE) then
            call s_dense_3d_nodal(sp(:,:,:,:), cp(:,:,:,1), mp(:,:,:,1), &
                                  face_type(i,:,:), dh)
          else if (stencil_type == ST_CROSS) then
            call s_cross_3d_nodal(sp(:,:,:,:), cp(:,:,:,1), mp(:,:,:,1), &
                                   face_type(i,:,:), dh)
          else 
            print*,'DONT KNOW THIS NODAL STENCIL TYPE ',stencil_type
            call bl_error('stencil_fill_nodal')
          end if
       end select
    end do

    call destroy(bxa_periodic)
    call destroy(bpt)

  end subroutine stencil_fill_nodal

  subroutine stencil_fill_one_sided(ss, sg, dh, mask, face_type)

    use stencil_nodal_module 

    type(multifab ), intent(inout) :: ss
    type(multifab ), intent(inout) :: sg
    real(kind=dp_t), intent(in   ) :: dh(:)
    type(imultifab), intent(inout) :: mask
    integer        , intent(in   ) :: face_type(:,:,:)

    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer                  :: i

    do i = 1, ss%nboxes
       if ( multifab_remote(ss,i) ) cycle

       sp => dataptr(ss,   i)
       cp => dataptr(sg,   i)
       mp => dataptr(mask, i)

       select case (ss%dim)
       case (1)
!         call s_simple_1d_one_sided(sp(:,1,1,:), cp(:,1,1,1), mp(:,1,1,1), face_type(i,1,:), dh)
       case (2)
          call s_simple_2d_one_sided(sp(:,:,1,:), cp(:,:,1,1), mp(:,:,1,1), &
                                     face_type(i,:,:), dh)
       case (3)
          call s_simple_3d_one_sided(sp(:,:,:,:), cp(:,:,:,1), mp(:,:,:,1), &
                                     face_type(i,:,:), dh)
       end select
    end do

  end subroutine stencil_fill_one_sided

  recursive subroutine stencil_fill_cc_all_mglevels(mgt, cell_coeffs, edge_coeffs, &
                                                    xa, xb, pxa, pxb, &
                                                    stencil_order, bc_face)

    use coarsen_coeffs_module
    use mg_module

    type(mg_tower) , intent(inout) :: mgt
    type(multifab) , intent(inout) :: cell_coeffs(:)
    type(multifab) , intent(inout) :: edge_coeffs(:,:)
    real(kind=dp_t), intent(in   ) :: xa(:), xb(:), pxa(:), pxb(:)
    integer        , intent(in   ) :: stencil_order
    integer        , intent(in   ) :: bc_face(:,:)

    ! Local variables
    integer                     :: i, d, dm, maxlev, maxlev_bottom
    real(dp_t)                  :: coarse_xa(mgt%dim),  coarse_xb(mgt%dim)
    real(dp_t)                  :: coarse_pxa(mgt%dim), coarse_pxb(mgt%dim)
    type(layout)                :: old_la_grown, new_la_grown
    type(boxarray)              :: ba_cc
    type(multifab)              :: old_edge_coeffs_grown
    type(multifab)              :: new_edge_coeffs_grown
    type(multifab), allocatable :: coarse_cell_coeffs(:)
    type(multifab), allocatable :: coarse_edge_coeffs(:,:)
    real(dp_t), pointer :: sc_orig(:,:,:,:), sc_grown(:,:,:,:)

        dm = mgt%dim
    maxlev = mgt%nlevels

    ! NOTE: coeffs(maxlev) comes in built and filled, but the other levels
    !       are not even built yet
    do i = maxlev-1, 1, -1
       call multifab_build(cell_coeffs(i), mgt%ss(i)%la, cell_coeffs(maxlev)%nc, cell_coeffs(maxlev)%ng)
       call setval(cell_coeffs(i), ZERO, 1, cell_coeffs(maxlev)%nc, all=.true.)
       call coarsen_cell_coeffs ( cell_coeffs(i+1), cell_coeffs(i))
       call multifab_fill_boundary(cell_coeffs(i))
    end do

    do i = maxlev-1, 1, -1
       do d = 1,dm
          call multifab_build_edge(edge_coeffs(i,d), mgt%ss(i)%la, edge_coeffs(maxlev,d)%nc, edge_coeffs(maxlev,d)%ng, d)
          call setval(edge_coeffs(i,d), ZERO, 1, edge_coeffs(maxlev,d)%nc, all=.true.)
       end do
       call coarsen_edge_coeffs(edge_coeffs(i+1,:),edge_coeffs(i,:))
       do d = 1,dm
          call multifab_fill_boundary(edge_coeffs(i,d))
       end do
    end do

    do i = maxlev, 1, -1
       call stencil_fill_cc(mgt%ss(i), cell_coeffs(i), edge_coeffs(i,:), &
                            mgt%dh(:,i), mgt%mm(i), xa, xb, pxa, pxb, stencil_order, bc_face) 
    end do

    if (associated(mgt%bottom_mgt)) then

       maxlev_bottom = mgt%bottom_mgt%nlevels

       ! First we copy just the cell-centered component -- which does not need ghost cells copied (I think)
       allocate(coarse_cell_coeffs(maxlev_bottom))
       call multifab_build(coarse_cell_coeffs(maxlev_bottom), mgt%bottom_mgt%cc(maxlev_bottom)%la, &
                           cell_coeffs(1)%nc,cell_coeffs(1)%ng)
       call setval(coarse_cell_coeffs(maxlev_bottom),ZERO,all=.true.)
       call multifab_copy_c(coarse_cell_coeffs(maxlev_bottom),1,cell_coeffs(1),1,cell_coeffs(1)%nc,ng=0)
       call multifab_fill_boundary(coarse_cell_coeffs(maxlev_bottom))

       ! Make space for the coarsened edge coefficients but don't copy directly
       allocate(coarse_edge_coeffs(maxlev_bottom,dm))
       do d = 1, dm 
          call multifab_build_edge(coarse_edge_coeffs(maxlev_bottom,d), mgt%bottom_mgt%cc(maxlev_bottom)%la, &
                                   edge_coeffs(1,d)%nc,edge_coeffs(1,d)%ng,d)
          call setval(coarse_edge_coeffs(maxlev_bottom,d),ZERO,all=.true.)
       end do

       do d = 1, dm
          call boxarray_build_copy(ba_cc,get_boxarray(edge_coeffs(1,d)))
          call boxarray_grow(ba_cc,edge_coeffs(1,d)%ng)
          call layout_build_ba(old_la_grown,ba_cc,pmask=mgt%ss(1)%la%lap%pmask, &
                               explicit_mapping=get_proc(mgt%ss(1)%la))
          call destroy(ba_cc)
          call multifab_build_edge(old_edge_coeffs_grown,old_la_grown,edge_coeffs(1,d)%nc,0,d)

          do i = 1, old_edge_coeffs_grown%nboxes
             if (remote(old_edge_coeffs_grown,i)) cycle
             sc_orig  => dataptr(edge_coeffs(1,d)     ,i,get_pbox(old_edge_coeffs_grown,i),1,edge_coeffs(1,d)%nc)
             sc_grown => dataptr(old_edge_coeffs_grown,i,get_pbox(old_edge_coeffs_grown,i),1,edge_coeffs(1,d)%nc)
             sc_grown = sc_orig
          end do

          call boxarray_build_copy(ba_cc,get_boxarray(mgt%bottom_mgt%ss(maxlev_bottom)))
          call boxarray_grow(ba_cc,edge_coeffs(1,d)%ng)
          call layout_build_ba(new_la_grown,ba_cc,pmask=mgt%ss(1)%la%lap%pmask, &
                               explicit_mapping=get_proc(mgt%bottom_mgt%ss(maxlev_bottom)%la))
          call destroy(ba_cc)
          call multifab_build_edge(new_edge_coeffs_grown,new_la_grown,edge_coeffs(1,d)%nc,0,d)
          call multifab_copy_c(new_edge_coeffs_grown,1,old_edge_coeffs_grown,1,nc=edge_coeffs(1,d)%nc)

          call destroy(old_edge_coeffs_grown)
          call destroy(old_la_grown)

          do i = 1, new_edge_coeffs_grown%nboxes
             if (remote(new_edge_coeffs_grown,i)) cycle
             sc_orig  => dataptr(coarse_edge_coeffs(maxlev_bottom,d),i,get_pbox(new_edge_coeffs_grown,i),1,edge_coeffs(1,d)%nc)
             sc_grown => dataptr(new_edge_coeffs_grown              ,i,get_pbox(new_edge_coeffs_grown,i),1,edge_coeffs(1,d)%nc)
             sc_orig = sc_grown
          end do

          call destroy(new_edge_coeffs_grown)
          call destroy(new_la_grown)

       end do

       do i = maxlev_bottom-1, 1, -1
          call multifab_build(coarse_cell_coeffs(i),mgt%bottom_mgt%ss(i)%la,cell_coeffs(1)%nc,cell_coeffs(1)%ng)
          call setval(coarse_cell_coeffs(i),ZERO,1,cell_coeffs(1)%nc,all=.true.)
          call coarsen_cell_coeffs(coarse_cell_coeffs(i+1),coarse_cell_coeffs(i))
          call multifab_fill_boundary(coarse_cell_coeffs(i))
       end do

       do i = maxlev_bottom-1, 1, -1
          do d = 1,dm
             call multifab_build(coarse_edge_coeffs(i,d),mgt%bottom_mgt%ss(i)%la,edge_coeffs(1,d)%nc,edge_coeffs(1,d)%ng)
             call setval(coarse_edge_coeffs(i,d),ZERO,1,edge_coeffs(1,d)%nc,all=.true.)
          end do
             call coarsen_edge_coeffs(coarse_edge_coeffs(i+1,:),coarse_edge_coeffs(i,:))
          do d = 1,dm
             call multifab_fill_boundary(coarse_edge_coeffs(i,d))
          end do
       end do

       coarse_xa = ZERO
       coarse_xb = ZERO
       coarse_pxa = ZERO
       coarse_pxb = ZERO

       call stencil_fill_cc_all_mglevels(mgt%bottom_mgt, coarse_cell_coeffs, coarse_edge_coeffs, &
                                         coarse_xa, coarse_xb, coarse_pxa, coarse_pxb, stencil_order, bc_face)

       call destroy(coarse_cell_coeffs(maxlev_bottom))
       deallocate(coarse_cell_coeffs)

       do d = 1,dm
          call destroy(coarse_edge_coeffs(maxlev_bottom,d))
       end do
       deallocate(coarse_edge_coeffs)

    end if

    do i = maxlev-1, 1, -1
       call destroy(cell_coeffs(i))
       do d = 1,dm
          call destroy(edge_coeffs(i,d))
       end do
    end do

  end subroutine stencil_fill_cc_all_mglevels

  subroutine stencil_fill_cc(ss, cell_coeffs, edge_coeffs, &
                             dh, mask, xa, xb, pxa, pxb, order, bc_face)

    use bl_prof_module

    type(multifab) , intent(inout) :: ss
    type(multifab) , intent(in   ) :: cell_coeffs
    type(multifab) , intent(in   ) :: edge_coeffs(:)
    real(kind=dp_t), intent(in   ) :: dh(:)
    type(imultifab), intent(inout) :: mask
    integer        , intent(in   ) :: order
    integer        , intent(in   ) :: bc_face(:,:)
    real(kind=dp_t), intent(in   ) :: xa(:), xb(:), pxa(:), pxb(:)

    type(box)                 :: bx, pd
    real(kind=dp_t)           :: lxa(ss%dim), lxb(ss%dim)

    real(kind=dp_t), pointer  ::  sp(:,:,:,:)
    real(kind=dp_t), pointer  :: ccp(:,:,:,:)
    real(kind=dp_t), pointer  :: xcp(:,:,:,:)
    real(kind=dp_t), pointer  :: ycp(:,:,:,:)
    real(kind=dp_t), pointer  :: zcp(:,:,:,:)
    integer        , pointer  ::  mp(:,:,:,:)
    integer                   :: i,ns,ng_b,ng_c,id,ncomp_coeffs
    logical                   :: minion_stencil
    type(bl_prof_timer), save :: bpt

    call build(bpt, "stencil_fill_cc")

    pd = layout_get_pd(ss%la)

    minion_stencil = .false.

    if (cell_coeffs%ng .eq. 2) then

       minion_stencil = .true.
       print *,'MINION STENCIL '

    endif 

    do i = 1, ss%nboxes
       if ( multifab_remote(ss,i) ) cycle
       bx = get_box(ss,i)
       call stencil_set_bc(ss, i, mask%fbs(i), bc_face)
       lxa = xa
       lxb = xb
       do id = 1,pd%dim
          if ( .not. ss%la%lap%pmask(id) ) then
             if ( bx%lo(id) == pd%lo(id) ) then
                lxa(id) = pxa(id)
             end if
             if ( bx%hi(id) == pd%hi(id) ) then
                lxb(id) = pxb(id)
             end if
          end if
       end do

       sp => dataptr(ss, i)
       ccp => dataptr(cell_coeffs, i)
       xcp => dataptr(edge_coeffs(1), i)
       mp => dataptr(mask, i)

       ng_c = cell_coeffs%ng
       ng_b = edge_coeffs(1)%ng

       if (minion_stencil) then

          ns   = ss%nc

          ycp => dataptr(edge_coeffs(2), i)

          select case (ss%dim)
          case (2)
             if (ns .eq. 7) then
                call s_minion_second_fill_2d(sp(:,:,1,:), ccp(:,:,1,1), ng_c, &
                                             xcp(:,:,1,1), ycp(:,:,1,1), ng_b, & 
                                             dh, mp(:,:,1,1), &
                                             bx%lo, bx%hi, lxa, lxb)
             else if (ns .eq. 9) then
                call s_minion_cross_fill_2d(sp(:,:,1,:), ccp(:,:,1,1), ng_c, &
                                            xcp(:,:,1,1), ycp(:,:,1,1), ng_b, & 
                                            dh, mp(:,:,1,1), &
                                            bx%lo, bx%hi, lxa, lxb)
             else if (ns .eq. 25) then
                call s_minion_full_fill_2d(sp(:,:,1,:), ccp(:,:,1,1), ng_c, &
                                           xcp(:,:,1,1), ycp(:,:,1,1), ng_b, & 
                                           dh, mp(:,:,1,1), &
                                           bx%lo, bx%hi, lxa, lxb)
             end if
          case (3)
             zcp => dataptr(edge_coeffs(3), i)
             if (ns .eq. 10) then
                call s_minion_second_fill_3d(sp(:,:,:,:), ccp(:,:,:,1), ng_c, &
                                             xcp(:,:,:,1), ycp(:,:,:,1), zcp(:,:,:,1), ng_b, & 
                                             dh, mp(:,:,:,1), &
                                             bx%lo, bx%hi, lxa, lxb)
             else if (ns .eq. 13) then
                call s_minion_cross_fill_3d(sp(:,:,:,:), ccp(:,:,:,1), ng_c, &
                                            xcp(:,:,:,1), ycp(:,:,:,1), zcp(:,:,:,1), ng_b, & 
                                            dh, mp(:,:,:,1), &
                                            bx%lo, bx%hi, lxa, lxb)
!            else if (ns .eq. 125) then
!               call s_minion_full_fill_3d(sp(:,:,:,:), ccp(:,:,:,1), ng_c, &
!                                          xcp(:,:,:,1), ycp(:,:,:,1), zcp(:,:,:,1), ng_b, & 
!                                          dh, mp(:,:,:,1), &
!                                          bx%lo, bx%hi, lxa, lxb)
             end if
          end select

       else

          ncomp_coeffs = multifab_ncomp(edge_coeffs(1)) 

          select case (ss%dim)
          case (1)
             call s_simple_1d_cc(sp(:,1,1,:), ccp(:,1,1,1), ng_c, xcp(:,1,1,1), ng_b, dh, &
                                 mp(:,1,1,1), bx%lo, bx%hi, lxa, lxb, order)
          case (2)
             ycp => dataptr(edge_coeffs(2), i)
             if (ncomp_coeffs > 1) then
                call s_simplen_2d_cc(sp(:,:,1,:), ccp(:,:,1,:), ng_c, xcp(:,:,1,:), ycp(:,:,1,:), ng_b, dh, &
                                     mp(:,:,1,1), bx%lo, bx%hi, lxa, lxb, order)
             else
                call s_simple_2d_cc(sp(:,:,1,:), ccp(:,:,1,1), ng_c, xcp(:,:,1,1), ycp(:,:,1,1), ng_b, dh, &
                                    mp(:,:,1,1), bx%lo, bx%hi, lxa, lxb, order)
             endif
          case (3)
             ycp => dataptr(edge_coeffs(2), i)
             zcp => dataptr(edge_coeffs(3), i)
             call s_simple_3d_cc(sp(:,:,:,:), ccp(:,:,:,1), ng_c, xcp(:,:,:,1), ycp(:,:,:,1), zcp(:,:,:,1), ng_b, &
                                 dh, mp(:,:,:,1), bx%lo, bx%hi, lxa, lxb, order)
          end select

       end if

    end do
    
    call destroy(bpt)

  end subroutine stencil_fill_cc

end module stencil_fill_module
