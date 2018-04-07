module cc_stencil_fill_module

  use bl_types
  use bc_module
  use multifab_module
  use cc_stencil_module

  implicit none

  private

  public :: stencil_fill_cc_all_mglevels, stencil_fill_const_all_mglevels, stencil_fill_cc

contains

  recursive subroutine stencil_fill_cc_all_mglevels(mgt, cell_coeffs, edge_coeffs, &
                                                    xa, xb, pxa, pxb, &
                                                    stencil_order, bc_face, nc_opt)

    use coarsen_coeffs_module
    use mg_tower_module

    type(mg_tower) , intent(inout) :: mgt
    type(multifab) , intent(inout) :: cell_coeffs(:)
    type(multifab) , intent(inout) :: edge_coeffs(:,:)
    real(kind=dp_t), intent(in   ) :: xa(:), xb(:), pxa(:), pxb(:)
    integer        , intent(in   ) :: stencil_order
    integer        , intent(in   ) :: bc_face(:,:)
    integer        , intent(in   ), optional :: nc_opt

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
       call multifab_build(cell_coeffs(i), get_layout(mgt%ss(i)), ncomp(cell_coeffs(maxlev)), nghost(cell_coeffs(maxlev)))
       call setval(cell_coeffs(i), ZERO, 1, ncomp(cell_coeffs(maxlev)), all=.true.)
       call coarsen_cell_coeffs ( cell_coeffs(i+1), cell_coeffs(i))
       if (nghost(cell_coeffs(i)) .gt. 0) &
            call multifab_fill_boundary(cell_coeffs(i))
    end do

    do i = maxlev-1, 1, -1
       do d = 1,dm
          call multifab_build_edge(edge_coeffs(i,d), get_layout(mgt%ss(i)), ncomp(edge_coeffs(maxlev,d)), &
               nghost(edge_coeffs(maxlev,d)), d)
          call setval(edge_coeffs(i,d), ZERO, 1, ncomp(edge_coeffs(maxlev,d)), all=.true.)
       end do
       call coarsen_edge_coeffs(edge_coeffs(i+1,:),edge_coeffs(i,:))
       do d = 1,dm
          if (nghost(edge_coeffs(i,d)) .gt. 0) &
               call multifab_fill_boundary(edge_coeffs(i,d))
       end do
    end do

    do i = maxlev, 1, -1
       call stencil_fill_cc(mgt%ss(i), cell_coeffs(i), edge_coeffs(i,:), &
                            mgt%dh(:,i), mgt%mm(i), xa, xb, pxa, pxb, stencil_order, bc_face, nc_opt) 
    end do

    if (associated(mgt%bottom_mgt)) then

       maxlev_bottom = mgt%bottom_mgt%nlevels

       ! First we copy just the cell-centered component -- which does not need ghost cells copied (I think)
       allocate(coarse_cell_coeffs(maxlev_bottom))
       call multifab_build(coarse_cell_coeffs(maxlev_bottom), get_layout(mgt%bottom_mgt%cc(maxlev_bottom)), &
            ncomp(cell_coeffs(1)), nghost(cell_coeffs(1)))
       call setval(coarse_cell_coeffs(maxlev_bottom),ZERO,all=.true.)
       call multifab_copy_c(coarse_cell_coeffs(maxlev_bottom),1,cell_coeffs(1),1,ncomp(cell_coeffs(1)),ng=0)
       if (nghost(coarse_cell_coeffs(maxlev_bottom)) .gt. 0) &
            call multifab_fill_boundary(coarse_cell_coeffs(maxlev_bottom))

       ! Make space for the coarsened edge coefficients but don't copy directly
       allocate(coarse_edge_coeffs(maxlev_bottom,dm))
       do d = 1, dm 
          call multifab_build_edge(coarse_edge_coeffs(maxlev_bottom,d), get_layout(mgt%bottom_mgt%cc(maxlev_bottom)), &
               ncomp(edge_coeffs(1,d)), nghost(edge_coeffs(1,d)), d)
          call setval(coarse_edge_coeffs(maxlev_bottom,d),ZERO,all=.true.)
       end do

       do d = 1, dm
          call boxarray_build_copy(ba_cc,get_boxarray(edge_coeffs(1,d)))
          call boxarray_grow(ba_cc,nghost(edge_coeffs(1,d)))
          call layout_build_ba(old_la_grown,ba_cc,boxarray_bbox(ba_cc),pmask = get_pmask(get_layout(mgt%ss(1))), &
               explicit_mapping = get_proc(get_layout(mgt%ss(1))))
          call boxarray_destroy(ba_cc)
          call multifab_build_edge(old_edge_coeffs_grown,old_la_grown,ncomp(edge_coeffs(1,d)),0,d)

          do i = 1, nfabs(old_edge_coeffs_grown)
             sc_orig  => dataptr(edge_coeffs(1,d)     ,i,get_pbox(old_edge_coeffs_grown,i),1,ncomp(edge_coeffs(1,d)))
             sc_grown => dataptr(old_edge_coeffs_grown,i,get_pbox(old_edge_coeffs_grown,i),1,ncomp(edge_coeffs(1,d)))
             sc_grown = sc_orig
          end do

          call boxarray_build_copy(ba_cc,get_boxarray(mgt%bottom_mgt%ss(maxlev_bottom)))
          call boxarray_grow(ba_cc,nghost(edge_coeffs(1,d)))
          call layout_build_ba(new_la_grown,ba_cc,boxarray_bbox(ba_cc),pmask = get_pmask(get_layout(mgt%ss(1))), &
               explicit_mapping = get_proc(get_layout(mgt%bottom_mgt%ss(maxlev_bottom))))
          call boxarray_destroy(ba_cc)
          call multifab_build_edge(new_edge_coeffs_grown,new_la_grown,ncomp(edge_coeffs(1,d)),0,d)
          call multifab_copy_c(new_edge_coeffs_grown,1,old_edge_coeffs_grown,1,nc=ncomp(edge_coeffs(1,d)))

          call multifab_destroy(old_edge_coeffs_grown)
          call layout_destroy(old_la_grown)

          do i = 1, nfabs(new_edge_coeffs_grown)
             sc_orig  => dataptr(coarse_edge_coeffs(maxlev_bottom,d),i, &
                  get_pbox(new_edge_coeffs_grown,i),1,ncomp(edge_coeffs(1,d)))
             sc_grown => dataptr(new_edge_coeffs_grown              ,i, &
                  get_pbox(new_edge_coeffs_grown,i),1,ncomp(edge_coeffs(1,d)))
             sc_orig = sc_grown
          end do

          call multifab_destroy(new_edge_coeffs_grown)
          call layout_destroy(new_la_grown)

       end do

       coarse_xa = ZERO
       coarse_xb = ZERO
       coarse_pxa = ZERO
       coarse_pxb = ZERO

       call stencil_fill_cc_all_mglevels(mgt%bottom_mgt, coarse_cell_coeffs, coarse_edge_coeffs, &
                                         coarse_xa, coarse_xb, coarse_pxa, coarse_pxb, stencil_order, bc_face, nc_opt)

       call multifab_destroy(coarse_cell_coeffs(maxlev_bottom))
       deallocate(coarse_cell_coeffs)

       do d = 1,dm
          call multifab_destroy(coarse_edge_coeffs(maxlev_bottom,d))
       end do
       deallocate(coarse_edge_coeffs)

    end if

    do i = maxlev-1, 1, -1
       call multifab_destroy(cell_coeffs(i))
       do d = 1,dm
          call multifab_destroy(edge_coeffs(i,d))
       end do
    end do

  end subroutine stencil_fill_cc_all_mglevels

! ******************************************************************************

  subroutine stencil_fill_cc(ss, cell_coeffs, edge_coeffs, &
                             dh, mask, xa, xb, pxa, pxb, order, bc_face, nc_opt)

    use bl_prof_module

    type(multifab) , intent(inout) :: ss
    type(multifab) , intent(in   ) :: cell_coeffs
    type(multifab) , intent(in   ) :: edge_coeffs(:)
    real(kind=dp_t), intent(in   ) :: dh(:)
    type(imultifab), intent(inout) :: mask
    integer        , intent(in   ) :: order
    integer        , intent(in   ) :: bc_face(:,:)
    real(kind=dp_t), intent(in   ) :: xa(:), xb(:), pxa(:), pxb(:)
    integer        , intent(in   ), optional :: nc_opt

    type(box)                 :: bx, pd
    real(kind=dp_t)           :: lxa(get_dim(ss)), lxb(get_dim(ss))

    real(kind=dp_t), pointer  ::  sp(:,:,:,:)
    real(kind=dp_t), pointer  :: ccp(:,:,:,:)
    real(kind=dp_t), pointer  :: xcp(:,:,:,:)
    real(kind=dp_t), pointer  :: ycp(:,:,:,:)
    real(kind=dp_t), pointer  :: zcp(:,:,:,:)
    integer        , pointer  ::  mp(:,:,:,:)
    integer                   :: i,ns,ng_b,ng_c,id,ncomp_coeffs,dm
    integer                   :: lnc_opt
    logical                   :: minion_stencil, pmask(get_dim(ss))

    type(bl_prof_timer), save :: bpt

    call build(bpt, "stencil_fill_cc")

    lnc_opt = 0
    if (present (nc_opt)) lnc_opt = nc_opt

    pd = get_pd(get_layout(ss))

    minion_stencil = .false.

    if ( nghost(cell_coeffs) .eq. 2 ) then

       minion_stencil = .true.

    endif 

    pmask = get_pmask(get_layout(ss))

    dm = get_dim(ss)

    do i = 1, nfabs(ss)
       bx = get_box(ss,i)
       call stencil_set_bc(ss, i, mask, bc_face)
       lxa = xa
       lxb = xb
       do id = 1,pd%dim
          if ( .not. pmask(id) ) then
             if ( bx%lo(id) == pd%lo(id) ) then
                lxa(id) = pxa(id)
             end if
             if ( bx%hi(id) == pd%hi(id) ) then
                lxb(id) = pxb(id)
             end if
          end if
       end do

       sp  => dataptr(ss, i)
       ccp => dataptr(cell_coeffs, i)
       xcp => dataptr(edge_coeffs(1), i)
       mp  => dataptr(mask, i)

       ng_c = nghost(cell_coeffs)
       ng_b = nghost(edge_coeffs(1))

       if (minion_stencil) then

          ns   = ncomp(ss)

          ycp => dataptr(edge_coeffs(2), i)

          select case (dm)
          case (2)
             if (ns .eq. 7) then
                call s_simple_2d_cc(sp(:,:,:,1), ccp(:,:,1,1), ng_c, &
                                    xcp(:,:,1,1), ycp(:,:,1,1), ng_b, dh, &
                                    mp(:,:,1,1), bx%lo, bx%hi, lxa, lxb, order)
             else if (ns .eq. 9) then
                call s_minion_cross_fill_2d(sp(:,:,:,1), ccp(:,:,1,1), ng_c, &
                                            xcp(:,:,1,1), ycp(:,:,1,1), ng_b, & 
                                            dh, mp(:,:,1,1), &
                                            bx%lo, bx%hi)
             else if (ns .eq. 25) then
                call s_minion_full_fill_2d(sp(:,:,:,1), ccp(:,:,1,1), ng_c, &
                                           xcp(:,:,1,1), ycp(:,:,1,1), ng_b, & 
                                           dh, mp(:,:,1,1), &
                                           bx%lo, bx%hi)
             end if
          case (3)
             zcp => dataptr(edge_coeffs(3), i)
             if (ns .eq. 10) then
                call s_simple_3d_cc(sp(:,:,:,:), ccp(:,:,:,1), ng_c, &
                                    xcp(:,:,:,1), ycp(:,:,:,1), zcp(:,:,:,1), ng_b, &
                                    dh, mp(:,:,:,1), bx%lo, bx%hi, lxa, lxb, order)
             else if (ns .eq. 13) then
                call s_minion_cross_fill_3d(sp(:,:,:,:), ccp(:,:,:,1), ng_c, &
                                            xcp(:,:,:,1), ycp(:,:,:,1), zcp(:,:,:,1), ng_b, & 
                                            dh, mp(:,:,:,1), &
                                            bx%lo, bx%hi)
             else if (ns .eq. 61) then
                call s_minion_full_fill_3d(sp(:,:,:,:), ccp(:,:,:,1), ng_c, &
                                           xcp(:,:,:,1), ycp(:,:,:,1), zcp(:,:,:,1), ng_b, & 
                                           dh, mp(:,:,:,1), &
                                           bx%lo, bx%hi)
             end if
          end select

       else

          ncomp_coeffs = multifab_ncomp(edge_coeffs(1)) 

          select case (dm)
          case (1)
             call s_simple_1d_cc(sp(:,:,1,1), ccp(:,1,1,1), ng_c, xcp(:,1,1,1), ng_b, dh, &
                                 mp(:,1,1,1), bx%lo, bx%hi, lxa, lxb, order)
          case (2)
             ycp => dataptr(edge_coeffs(2), i)
             if (ncomp_coeffs > 1) then
                if (lnc_opt .eq. 0) then
                   call s_simplen_2d_cc(sp(:,:,:,1), ccp(:,:,1,:), ng_c, xcp(:,:,1,:), ycp(:,:,1,:), ng_b, dh, &
                        mp(:,:,1,1), bx%lo, bx%hi, lxa, lxb, order)
                elseif (lnc_opt .eq. 1) then
                   call s_simplem_2d_cc(sp(:,:,:,1), ccp(:,:,1,:), ng_c, xcp(:,:,1,:), ycp(:,:,1,:), ng_b, dh, &
                        mp(:,:,1,1), bx%lo, bx%hi, lxa, lxb, order)     
                elseif (lnc_opt .eq. 2) then
                   call s_simpleg_2d_cc(sp(:,:,:,1), ccp(:,:,1,:), ng_c, xcp(:,:,1,:), ycp(:,:,1,:), ng_b, dh, &
                        mp(:,:,1,1), bx%lo, bx%hi)
                end if
             else
                call s_simple_2d_cc(sp(:,:,:,1), ccp(:,:,1,1), ng_c,&
                                    xcp(:,:,1,1), ycp(:,:,1,1), ng_b, dh, &
                                    mp(:,:,1,1), bx%lo, bx%hi, lxa, lxb, order)
             endif
          case (3)
             ycp => dataptr(edge_coeffs(2), i)
             zcp => dataptr(edge_coeffs(3), i)
             if (ncomp_coeffs > 1) then
                 if (lnc_opt .eq. 2) then 
                    call s_simpleg_3d_cc(sp(:,:,:,:), ccp(:,:,:,:), ng_c, xcp(:,:,:,:), ycp(:,:,:,:), zcp(:,:,:,:), ng_b, &
                      dh, mp(:,:,:,1), bx%lo, bx%hi)
                 end if
              else
                 call s_simple_3d_cc(sp(:,:,:,:), ccp(:,:,:,1), ng_c, &
                                     xcp(:,:,:,1), ycp(:,:,:,1), zcp(:,:,:,1), ng_b, &
                                     dh, mp(:,:,:,1), bx%lo, bx%hi, lxa, lxb, order)
              endif
          end select

       end if

    end do
    
    call destroy(bpt)

  end subroutine stencil_fill_cc

! ******************************************************************************

  recursive subroutine stencil_fill_const_all_mglevels(mgt, alpha_const, beta_const, &
                                                       xa, xb, pxa, pxb, &
                                                       stencil_order, bc_face)

    use coarsen_coeffs_module
    use mg_tower_module

    type(mg_tower) , intent(inout) :: mgt
    real(kind=dp_t), intent(in   ) :: xa(:), xb(:), pxa(:), pxb(:)
    integer        , intent(in   ) :: stencil_order
    integer        , intent(in   ) :: bc_face(:,:)
    real(dp_t)     , intent(in   ) :: alpha_const, beta_const

    ! Local variables
    integer                     :: i, maxlev
    real(dp_t)                  :: coarse_xa(mgt%dim),  coarse_xb(mgt%dim)
    real(dp_t)                  :: coarse_pxa(mgt%dim), coarse_pxb(mgt%dim)

    maxlev = mgt%nlevels

    do i = maxlev, 1, -1
       call stencil_fill_const(mgt%ss(i), alpha_const, beta_const, &
                               mgt%dh(:,i), mgt%mm(i), &
                               xa, xb, pxa, pxb, stencil_order, bc_face) 
    end do

    if (associated(mgt%bottom_mgt)) then

       coarse_xa = ZERO
       coarse_xb = ZERO
       coarse_pxa = ZERO
       coarse_pxb = ZERO

       call stencil_fill_const_all_mglevels(mgt%bottom_mgt, alpha_const, beta_const, &
                                            coarse_xa,  coarse_xb, &
                                            coarse_pxa, coarse_pxb, stencil_order, bc_face)

    end if


  end subroutine stencil_fill_const_all_mglevels

! ******************************************************************************

  subroutine stencil_fill_const(ss, alpha_const, beta_const, &
                                dh, mask, xa, xb, pxa, pxb, order, bc_face)

    use bl_prof_module
    use stencil_util_module, only : make_ibc_stencil_fab, simple_ib_const

    type(multifab) , intent(inout) :: ss
    real(kind=dp_t), intent(in   ) :: dh(:)
    type(imultifab), intent(inout) :: mask
    integer        , intent(in   ) :: order
    integer        , intent(in   ) :: bc_face(:,:)
    real(kind=dp_t), intent(in   ) :: alpha_const, beta_const
    real(kind=dp_t), intent(in   ) :: xa(:), xb(:), pxa(:), pxb(:)

    type(box)                 :: bx, pd
    real(kind=dp_t)           :: lxa(get_dim(ss)), lxb(get_dim(ss))

    real(kind=dp_t), pointer  ::  sp(:,:,:,:)
    integer        , pointer  ::  mp(:,:,:,:)
    integer                   :: i,id,dm
    logical                   :: pmask(get_dim(ss)), intbox

    type(bl_prof_timer), save :: bpt

    call build(bpt, "stencil_fill_const")

    pd = get_pd(get_layout(ss))

    pmask = get_pmask(get_layout(ss))

    dm = get_dim(ss)

    do i = 1, nfabs(ss)
       bx = get_box(ss,i)
       call stencil_set_bc(ss, i, mask, bc_face, intbox=intbox)

       if (intbox) then

          call make_ibc_stencil_fab(ss, i, dm)
          sp  => dataptr(ss, i)

          select case (dm)
          case (1)
             call bl_error("simple_1d_ib_const not yet implemented")
          case default 
             call simple_ib_const(sp(:,1,1,1), alpha_const, beta_const, dh, dm)
          end select

       else

          lxa = xa
          lxb = xb
          do id = 1,pd%dim
             if ( .not. pmask(id) ) then
                if ( bx%lo(id) == pd%lo(id) ) then
                   lxa(id) = pxa(id)
                end if
                if ( bx%hi(id) == pd%hi(id) ) then
                   lxb(id) = pxb(id)
                end if
             end if
          end do
          
          sp  => dataptr(ss, i)
          mp  => dataptr(mask, i)
          
          select case (dm)
          case (1)
             call bl_error("simple_1d_const not yet implemented")
          case (2)
             call simple_2d_const(sp(:,:,:,1),alpha_const,beta_const,&
                                  dh,mp(:,:,1,1), &
                                  bx%lo, bx%hi, lxa, lxb, order)
          case (3)
             call simple_3d_const(sp(:,:,:,:),alpha_const,beta_const,&
                                  dh,mp(:,:,:,1), &
                                  bx%lo, bx%hi, lxa, lxb, order)
          end select

       end if

    end do
    
    call destroy(bpt)

  end subroutine stencil_fill_const

! ******************************************************************************

end module cc_stencil_fill_module
