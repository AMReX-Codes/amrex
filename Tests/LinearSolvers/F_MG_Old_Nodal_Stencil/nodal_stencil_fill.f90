module nodal_stencil_fill_module

  use bl_types
  use bc_module
  use multifab_module
  use nodal_stencil_module

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
    use mg_tower_module

    type(mg_tower ), intent(inout) :: mgt
    type(multifab ), intent(inout) :: sg(:)
    integer        , intent(in   ) :: stencil_type

    type(multifab)                 :: stored_coeffs, stored_coeffs_grown
    type(multifab)                 :: new_coeffs_grown
    type(multifab), allocatable    :: coarse_coeffs(:)
    type(boxarray)                 :: ba_cc
    type(  layout)                 :: old_la_grown, new_la_grown
    integer                        :: i, maxlev, maxlev_bottom
    integer                        :: sg_ncomp
    real(dp_t), pointer            :: sc_orig(:,:,:,:), sc_grown(:,:,:,:)

    maxlev = mgt%nlevels
    sg_ncomp = ncomp(sg(maxlev))

    ! NOTE: sg(maxlev) comes in built and filled, but the other levels
    !       are not even built yet
    do i = maxlev-1, 1, -1
       call multifab_build(sg(i), get_layout(mgt%ss(i)), sg_ncomp, 1)
       call setval(sg(i), ZERO, all=.true.)
       call coarsen_cell_coeffs(sg(i+1),sg(i))
       call multifab_fill_boundary(sg(i))
    end do

    do i = maxlev, 1, -1
       call stencil_fill_nodal(mgt%ss(i), sg(i), mgt%dh(:,i), mgt%mm(i), &
                               mgt%face_type, stencil_type)
    end do

    if (associated(mgt%bottom_mgt)) then

       call multifab_build(stored_coeffs, get_layout(mgt%ss(1)), 1, 1)
       call multifab_copy_c(stored_coeffs,1,sg(1),1,1,ng = nghost(sg(1)))
       call multifab_fill_boundary(stored_coeffs)

       maxlev_bottom = mgt%bottom_mgt%nlevels
       allocate(coarse_coeffs(maxlev_bottom))
       call multifab_build(coarse_coeffs(maxlev_bottom),get_layout(mgt%bottom_mgt%cc(maxlev_bottom)),1,1)
       call setval(coarse_coeffs(maxlev_bottom),ZERO,all=.true.)

       ! Grow the stored coefficients
       call boxarray_build_copy(ba_cc,get_boxarray(stored_coeffs))
       call boxarray_grow(ba_cc,1)
       call layout_build_ba(old_la_grown,ba_cc,boxarray_bbox(ba_cc),pmask = get_pmask(get_layout(mgt%ss(1))), &
                            explicit_mapping=get_proc(get_layout(mgt%ss(1))))
       call destroy(ba_cc)
       call multifab_build(stored_coeffs_grown,old_la_grown,1,ng=0)

       do i = 1, nfabs(stored_coeffs_grown)
          sc_orig  => dataptr(stored_coeffs      ,i,get_pbox(stored_coeffs_grown,i),1,1)
          sc_grown => dataptr(stored_coeffs_grown,i,get_pbox(stored_coeffs_grown,i),1,1)
          sc_grown = sc_orig
       end do

       call boxarray_build_copy(ba_cc,get_boxarray(mgt%bottom_mgt%ss(maxlev_bottom)))
       call boxarray_grow(ba_cc,1)
       call layout_build_ba(new_la_grown,ba_cc,boxarray_bbox(ba_cc),pmask = get_pmask(get_layout(mgt%ss(1))), &
            explicit_mapping = get_proc(get_layout(mgt%bottom_mgt%ss(maxlev_bottom))))
       call destroy(ba_cc)
       call multifab_build(new_coeffs_grown,new_la_grown,1,ng=0)
       call multifab_copy_c(new_coeffs_grown,1,stored_coeffs_grown,1,1)

       do i = 1, nfabs(new_coeffs_grown)
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

    type(multifab ), intent(inout) :: ss
    type(multifab ), intent(inout) :: sg
    real(kind=dp_t), intent(in   ) :: dh(:)
    type(imultifab), intent(inout) :: mask
    integer        , intent(in   ) :: face_type(:,:,:)
    integer        , intent(in   ) :: stencil_type

    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)

    type(layout)             :: la
    type(box)                :: pd_periodic, bx, nbx, bx1, pd
    type(boxarray)           :: bxa_periodic, bxa_temp
    integer                  :: i, ib, jb, kb, ib_lo, jb_lo, kb_lo, dm
    integer                  :: shift_vect(get_dim(ss))
    integer                  :: lo(get_dim(ss)), hi(get_dim(ss))
    integer                  :: ng_sg
    type(list_box)           :: lb,nbxs
    type(box), allocatable   :: bxs(:)
    logical                  :: pmask(get_dim(ss))

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
    pd_periodic = get_pd(get_layout(ss))

    call boxarray_build_copy(bxa_periodic, get_boxarray(ss))

    pmask = get_pmask(get_layout(ss))

    dm = get_dim(ss)

    if ( any(pmask) ) then
       !
       ! First trim out all boxes that can't effect periodicity.
       !
       do i = 1,nboxes(bxa_periodic)
          if ( .not. contains(pd_periodic, get_box(bxa_periodic,i), strict = .true.) ) then
             call push_back(lb, get_box(bxa_periodic,i))
          end if
       end do

       do i = 1,dm
          if ( pmask(i) ) then
             pd_periodic = grow(grow(pd_periodic,1,i,-1),1,i,1)
          end if
       end do

       ib_lo = 1
       if ( pmask(1) )    ib_lo = -1

       jb_lo = 1
       if ( dm > 1 ) then
          if ( pmask(2) ) jb_lo = -1
       end if

       kb_lo = 1
       if ( dm > 2 ) then
          if ( pmask(3) ) kb_lo = -1
       end if

       pd = get_pd(get_layout(ss))

       do kb = kb_lo, 1
          do jb = jb_lo, 1
             do ib = ib_lo, 1
                call copy(bxa_temp,lb)

                shift_vect = 0

                if (    pmask(1) ) shift_vect(1) = ib * extent(pd,1)

                if ( dm > 1 ) then
                   if ( pmask(2) ) shift_vect(2) = jb * extent(pd,2)
                end if

                if ( dm > 2 ) then
                   if ( pmask(3) ) shift_vect(3) = kb * extent(pd,3)
                end if

                call boxarray_shift(bxa_temp,shift_vect)

                do i = 1, nboxes(bxa_temp)
                   bx1 = intersection(get_box(bxa_temp,i),pd_periodic)
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
    !
    ! Build layout on bxa_periodic.  We use a layout to make the stencil_set_bc_nodal() more efficient.
    !
    call build(la,bxa_periodic,get_pd(get_layout(ss)))

    ng_sg = nghost(sg)

    do i = 1, nfabs(ss)

       sp  => dataptr(ss,   i)
       cp  => dataptr(sg,   i)
       mp  => dataptr(mask, i)

       bx  = get_box(ss,i)
       nbx = get_ibox(ss, i)

       call stencil_set_bc_nodal(dm, bx, nbx, i, mask, face_type, pd_periodic, la)

       lo = lwb(get_box(sg,i))
       hi = upb(get_box(sg,i))

       select case (dm)
       case (1)
          call s_simple_1d_nodal(sp(:,:,1,1), cp(:,1,1,1), mp(:,1,1,1), dh)
       case (2)
          if (stencil_type == ND_DENSE_STENCIL) then
            call s_dense_2d_nodal(sp(:,:,:,1), cp(:,:,1,:), ng_sg, mp(:,:,1,1), &
                                  face_type(i,:,:), dh, lo, hi)
          else if (stencil_type == ND_CROSS_STENCIL) then
            call s_cross_2d_nodal(sp(:,:,:,1), cp(:,:,1,1), mp(:,:,1,1), &
                                   face_type(i,:,:), dh)
          else 
            print *,'DONT KNOW THIS NODAL STENCIL TYPE ',stencil_type
            call bl_error('stencil_fill_nodal')
          end if
       case (3)
          if (stencil_type == ND_DENSE_STENCIL) then
            call s_dense_3d_nodal(sp(:,:,:,:), cp(:,:,:,1), mp(:,:,:,1), dh)
          else if (stencil_type == ND_CROSS_STENCIL) then
            call s_cross_3d_nodal(sp(:,:,:,:), cp(:,:,:,1), mp(:,:,:,1), dh)
          else 
            print*,'DONT KNOW THIS NODAL STENCIL TYPE ',stencil_type
            call bl_error('stencil_fill_nodal')
          end if
       end select
    end do

    call destroy(bxa_periodic)
    call destroy(la)
    call destroy(bpt)

  end subroutine stencil_fill_nodal

  subroutine stencil_fill_one_sided(ss, sg, dh, mask, face_type)

    type(multifab ), intent(inout) :: ss
    type(multifab ), intent(inout) :: sg
    real(kind=dp_t), intent(in   ) :: dh(:)
    type(imultifab), intent(inout) :: mask
    integer        , intent(in   ) :: face_type(:,:,:)

    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    integer        , pointer :: mp(:,:,:,:)
    integer                  :: i, dm

    dm = get_dim(ss)

    do i = 1, nfabs(ss)

       sp  => dataptr(ss,   i)
       cp  => dataptr(sg,   i)
       mp  => dataptr(mask, i)

       select case (dm)
       case (1)
         call bl_error('s_simple_1d_one_sided() not implemented')
       case (2)
          call s_simple_2d_one_sided(sp(:,:,:,1), cp(:,:,1,1), mp(:,:,1,1), &
                                     face_type(i,:,:), dh)
       case (3)
          call s_simple_3d_one_sided(sp(:,:,:,:), cp(:,:,:,1), mp(:,:,:,1), dh)
       end select
    end do

  end subroutine stencil_fill_one_sided

end module nodal_stencil_fill_module
