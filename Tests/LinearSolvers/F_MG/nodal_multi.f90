subroutine t_nodal_ml_multigrid(mla, mgt, domain_bc, &
                                bottom_solver,do_diagnostics,eps,test, fabio, stencil_type)
  use BoxLib
  use stencil_module
  use coeffs_module
  use mg_module
  use list_box_module
  use ml_boxarray_module
  use ml_layout_module
  use itsol_module
  use bl_mem_stat_module
  use bl_timer_module
  use box_util_module
  use bl_IO_module
  use fabio_module
  use ml_nd_module
  use nodal_mask_module

  use ml_restriction_module
  use ml_prolongation_module
  use ml_interface_stencil_module
  use ml_util_module

  use bndry_reg_module

  implicit none

  type(ml_layout), intent(inout) :: mla
  type(mg_tower) , intent(inout) :: mgt(:)
  integer        , intent(in   ) :: domain_bc(:,:)
  integer        , intent(in   ) :: bottom_solver 
  integer        , intent(in   ) :: do_diagnostics 
  real(dp_t)     , intent(in   ) :: eps
  integer        , intent(in   ) :: test
  logical        , intent(in   ) :: fabio
  integer        , intent(in   ) :: stencil_type

  type(lmultifab), allocatable :: fine_mask(:)
  type( multifab), allocatable :: coeffs(:)
  type( multifab), allocatable :: full_soln(:)
  type( multifab), allocatable ::        rh(:)
  type( multifab), allocatable :: one_sided_ss(:)

  type(box)    :: pd
  type(layout) :: la
  integer      :: nlevs, n, i, dm, ns
  integer      :: stat
  real(dp_t)   :: snrm(2)

  logical, allocatable :: nodal(:)
  integer, allocatable :: ref_ratio(:,:)
  integer, allocatable :: lo(:), hi(:)

  logical :: make_solvable 

  make_solvable = .true.
  if ( any(domain_bc == BC_DIR) ) make_solvable = .false.

  dm = mla%dim

  allocate(lo(dm), hi(dm))
  allocate(nodal(dm))
  nodal = .True.

  nlevs = mla%nlevel

  allocate(full_soln(nlevs),rh(nlevs),fine_mask(nlevs))

  allocate(ref_ratio(nlevs-1,dm))
  do n = 1,nlevs-1
    ref_ratio(n,:) = mla%mba%rr(n,:)
  end do

  ! NOTE: we need to allocate for both stencils even for the dense stencil
  !       we don't actually create or use the multifabs
  allocate(one_sided_ss(2:nlevs))

  if (stencil_type .eq. ST_DENSE) then
     if (dm .eq. 3) then
       i = mgt(nlevs)%nlevels
       if ( (mgt(nlevs)%dh(1,i) .eq. mgt(nlevs)%dh(2,i)) .and. &
            (mgt(nlevs)%dh(1,i) .eq. mgt(nlevs)%dh(3,i)) ) then 
         ns = 21
       else
         ns = 27
       end if
     else if (dm .eq. 2) then
       ns = 9
     end if
     if ( parallel_ioprocessor() ) print *,'SETTING UP DENSE STENCIL WITH NS = ',ns
  else 
    ns = 2*dm+1
    do n = nlevs, 2, -1
      la = mla%la(n)
      call multifab_build(one_sided_ss(n), la, ns, 0, nodal)
      call setval(one_sided_ss(n), ZERO,all=.true.)
    end do
  end if

  do n = nlevs, 1, -1

     pd = mla%mba%pd(n)
     if ( parallel_ioprocessor() ) then
        print *,'LEVEL ',n
        print *,'PD ',extent(pd)
     end if
     la = mla%la(n)
     call multifab_build(full_soln(n), la, 1, 1, nodal)
     call multifab_build(       rh(n), la, 1, 1, nodal)
     call lmultifab_build(fine_mask(n), la, 1, 0, nodal)

     call setval(full_soln(n), ZERO,all=.true.)
     call setval(       rh(n), ZERO,all=.true.)
 
     if (n == nlevs) call mf_init(rh(n))

  end do

  !! Fill coefficient array

  do n = nlevs, 1, -1
   
     allocate(coeffs(mgt(n)%nlevels))

     la = get_layout(full_soln(n))
     pd = mla%mba%pd(n)

     if ( parallel_ioprocessor() ) print *,'N MG LEVS ',mgt(n)%nlevels

     call multifab_build(coeffs(mgt(n)%nlevels), la, 1, 1)
     call setval(coeffs(mgt(n)%nlevels), 0.0_dp_t, 1, all=.true.)
     do i = 1,layout_nboxes(la)
        call multifab_setval_bx(coeffs(mgt(n)%nlevels), 1.0_dp_t, get_box(la,i), all=.true.)
!       BEGIN HACK FOR 3-D MULTILEVEL PROBLEMS
        if (stencil_type .eq. ST_DENSE) then
         if (dm .eq. 3) then 
          if (nlevs .eq. 2) then 
            call multifab_setval_bx(coeffs(mgt(n)%nlevels), 0.5_dp_t, get_box(la,i), all=.true.)
          else if (nlevs > 2) then
!           FOR FAC 4
!           call multifab_setval_bx(coeffs(mgt(n)%nlevels), 0.0625_dp_t, get_box(la,i), all=.true.)
            call multifab_setval_bx(coeffs(mgt(n)%nlevels), 0.25_dp_t, get_box(la,i), all=.true.)
          end if
         end if
        end if
!       END HACK FOR 3-D MULTILEVEL PROBLEMS
     end do 
     call multifab_fill_boundary(coeffs(mgt(n)%nlevels))

     do i = mgt(n)%nlevels-1, 1, -1
        call multifab_build(coeffs(i), mgt(n)%ss(i)%la, 1, 1)
        call setval(coeffs(i), 0.0_dp_t, 1, all=.true.)
        call coarsen_coeffs(coeffs(i+1),coeffs(i))
        call multifab_fill_boundary(coeffs(i))
      end do

!    NOTE: we define the stencils with the finest dx.
     do i = mgt(n)%nlevels, 1, -1
        call stencil_fill_nodal(mgt(n)%ss(i), coeffs(i), mgt(n)%dh(:,i), &
                                mgt(n)%mm(i), mgt(n)%face_type, stencil_type)
 
     end do
     if (stencil_type .eq. ST_CROSS .and. n .gt. 1) then
        i = mgt(n)%nlevels
        call stencil_fill_one_sided(one_sided_ss(n), coeffs(i), mgt(n)%dh(:,i), &
                                    mgt(n)%mm(i), mgt(n)%face_type)
     end if

     call setval(fine_mask(n), val = .TRUE., all = .true.)

     if ( n < nlevs ) then
       call create_nodal_mask(n,fine_mask(n), &
                              mgt(n  )%mm(mgt(n  )%nlevels), &
                              mgt(n+1)%mm(mgt(n+1)%nlevels), &
                              mla)
     endif

     if ( test == 0 .and. n == 1 .and. bottom_solver == 3 ) then
        if ( parallel_ioprocessor() ) print *,'SPARSE BOTTOM SOLVER '
        call copy(mgt(n)%ss1, mgt(n)%ss(1))
        call copy(mgt(n)%mm1, mgt(n)%mm(1))
        if ( parallel_ioprocessor() ) then
           call sparse_nodal_build(mgt(n)%sparse_object, mgt(n)%ss1, &
                mgt(n)%mm1, mgt(n)%ss1%la, &
                mgt(n)%face_type, mgt(nlevs)%verbose)
        end if
     end if
     do i = mgt(n)%nlevels, 1, -1
        call multifab_destroy(coeffs(i))
     end do
     deallocate(coeffs)

  end do

  !!
  !! Solver Starts Here: We allow for MG or Sparse solver
  !!
  if (test == 0) then
     if ( parallel_ioprocessor() ) print *,'DOING MG SOLVER'
     continue
  else if (test == 3) then
     if ( parallel_ioprocessor() ) print *,'DOING SPARSE SOLVER AT TOP LEVEL '
     n = nlevs
     i = mgt(n)%nlevels
     if ( parallel_nprocs() > 1 ) then
        call bl_error("NODAL_MULTI: can't do parallel at top level")
     end if
     call sparse_nodal_build(mgt(n)%sparse_object, mgt(n)%ss(i), &
                             mgt(n)%mm(i), mgt(n)%ss(i)%la, &
                             mgt(n)%face_type,mgt(n)%verbose)
     call sparse_nodal_solve(mgt(n)%sparse_object, full_soln(n), rh(n), &
                             eps, mgt(n)%max_iter, mgt(n)%verbose, stat)
     if ( stat /= 0 ) then
        call bl_warn("SPARSE SOLVE FAILED DUE TO BREAKDOWN")
     end if
     if ( fabio ) then
       call fabio_ml_write(full_soln, ref_ratio(:,1), "sparse_soln_nodal")
     end if
     stop
  else 
     print *,'WE DO NOT SUPPORT THAT TEST TYPE FOR NODAL: ',test
     stop
  end if

  if ( fabio ) then
     call fabio_ml_write(rh, ref_ratio(:,1), "rh-init_nodal")
  end if

  snrm(2) = ml_norm_inf(rh,fine_mask)
  if ( parallel_IOProcessor() ) then
     print *, 'RHS MAX NORM ', snrm(2)
  end if

! ****************************************************************************

  call ml_nd(mla,mgt,rh,full_soln,fine_mask,one_sided_ss,ref_ratio,do_diagnostics,eps)

! ****************************************************************************

  if ( fabio ) then
     call fabio_ml_write(full_soln, ref_ratio(:,1), "soln_nodal")
  end if

  snrm(1) = ml_norm_l2(full_soln, ref_ratio, fine_mask)
  snrm(2) = ml_norm_inf(full_soln, fine_mask)
  if ( parallel_IOProcessor() ) then
     print *, 'SOLUTION MAX NORM ', snrm(2)
     print *, 'SOLUTION L2 NORM ', snrm(1)
  end if

  if ( parallel_IOProcessor() ) print *, 'MEMORY STATS'
  call print(multifab_mem_stats(),  " multifab before")
  call print(imultifab_mem_stats(), "imultifab before")
  call print(fab_mem_stats(),       "      fab before")
  call print(ifab_mem_stats(),      "     ifab before")
  call print(boxarray_mem_stats(),  " boxarray before")
  call print(boxassoc_mem_stats(),  " boxassoc before")
  call print(layout_mem_stats(),    "   layout before")

  do n = 1,nlevs
     call multifab_destroy(rh(n))
     call multifab_destroy(full_soln(n))
     call lmultifab_destroy(fine_mask(n))
  end do
  if (stencil_type == ST_CROSS) then
    do n = 2,nlevs
      call multifab_destroy(one_sided_ss(n))
    end do
  end if

contains

  subroutine mf_init(mf)
    type(multifab), intent(inout) :: mf
    integer i
    type(box) bx
    do i = 1, mf%nboxes; if ( remote(mf,i) ) cycle

!      Single point of non-zero RHS
       bx = get_ibox(mf,i)
       bx%lo(1:bx%dim) = (bx%hi(1:bx%dim) + bx%lo(1:bx%dim))/2
       bx%hi(1:bx%dim) = bx%lo(1:bx%dim)
       call setval(mf%fbs(i), ONE, bx)
       print *,'SETTING RHS TO  1 IN BOX ',i,' : ', bx%lo(1:bx%dim)

!      Single point of non-zero RHS: use this to make system solvable
       bx = get_ibox(mf,i)
       bx%lo(1       ) = (bx%hi(1       ) + bx%lo(1       ))/2 + 1
       bx%lo(2:bx%dim) = (bx%hi(2:bx%dim) + bx%lo(2:bx%dim))/2
       bx%hi(1:bx%dim) = bx%lo(1:bx%dim)
       call setval(mf%fbs(i), -ONE, bx)
       print *,'SETTING RHS TO -1 IN BOX ',i,' : ', bx%lo(1:bx%dim)

!      1-d strip of non-zero RHS in vertical
!      bx%lo(1) = (bx%hi(1) + bx%lo(1))/2
!      bx%hi(1) = bx%lo(1)

!      1-d strip of non-zero RHS in horizontal
!      bx%lo(2) = (bx%hi(2) + bx%lo(2))/2
!      bx%hi(2) = bx%lo(2)

    end do
  end subroutine mf_init

  subroutine mf_init1(mf)
    type(multifab), intent(inout) :: mf
    integer i
    type(box) bx
    type(box) rhs_box, rhs_intersect_box

    rhs_box%dim = mf%dim
    rhs_box%lo(1:rhs_box%dim) = 8
    rhs_box%hi(1:rhs_box%dim) = 8

    do i = 1, mf%nboxes

       if ( remote(mf,i) ) cycle
       bx = get_ibox(mf,i)
       rhs_intersect_box = box_intersection(bx,rhs_box)
       if (.not. empty(rhs_intersect_box)) then
         bx%lo(1:bx%dim) = lwb(rhs_intersect_box)
         bx%hi(1:bx%dim) = upb(rhs_intersect_box)
         call setval(mf%fbs(i), ONE, bx)
       end if
    end do

  end subroutine mf_init1

end subroutine t_nodal_ml_multigrid
