
subroutine t_cc_ml_multigrid(mla, mgt, rh, coeffs_type, domain_bc, do_diagnostics, stencil_order, fabio)

  use BoxLib
  use cc_stencil_module
  use cc_stencil_fill_module
  use ml_norm_module
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
  use ml_cc_restriction_module
  use ml_prolongation_module
  use ml_cc_module

  use bndry_reg_module

  use cc_rhs_module
  use cc_edge_coeffs_module
  use init_cell_coeffs_module

  implicit none

  type(ml_layout), intent(inout) :: mla
  type(mg_tower) , intent(inout) :: mgt(:)
  type( multifab), intent(inout) :: rh(:)

  integer        , intent(in   ) :: coeffs_type
  integer        , intent(in   ) :: domain_bc(:,:)
  integer        , intent(in   ) :: do_diagnostics 
  integer        , intent(in   ) :: stencil_order
  logical        , intent(in   ) :: fabio

  type(box      )                :: pd

  type(multifab), allocatable :: alpha(:)
  type(multifab), allocatable :: edge_coeffs(:,:)

  type( multifab), allocatable   :: full_soln(:)

  type(multifab)                 :: cell_coeffs

  type(layout)                   :: la
  real(dp_t)     , allocatable   :: xa(:), xb(:), pxa(:), pxb(:)

  integer        , allocatable   :: ref_ratio(:,:)
  integer                        :: d, n, dm, nlevs

  real(dp_t)                     :: snrm(2)

  real(dp_t)                     :: mac_beta

  dm = mla%dim

  nlevs = mla%nlevel

  allocate(full_soln(nlevs))
  allocate(xa(dm), xb(dm), pxa(dm), pxb(dm))

  allocate(ref_ratio(nlevs-1,dm))
  do n = 1,nlevs-1
    ref_ratio(n,:) = mla%mba%rr(n,:)
  end do

  ! NOTE THIS CHANGE: we now have stencil values which reach outside the
  !  grid in the case of Dirichlet bc's and skewed stencils

  do n = nlevs, 1, -1

     call multifab_build(full_soln(n), mla%la(n), 1, 1)
     call setval(full_soln(n), val = ZERO, all=.true.)

  end do

  !! Fill coefficient arrays

  mac_beta = 1.0_dp_t

  do n = nlevs, 1, -1

     la = mla%la(n)

     allocate(alpha(mgt(n)%nlevels))
     allocate(edge_coeffs(mgt(n)%nlevels,dm))

     call multifab_build(alpha(mgt(n)%nlevels),la,nc=1,ng=1)
     call setval(alpha(mgt(n)%nlevels),0.d0,all=.true.)

     do d = 1,dm
        call multifab_build_edge(edge_coeffs(mgt(n)%nlevels,d), la, nc=1, ng=0, dir=d)
     end do

     if (coeffs_type .eq. 0) then
        do d = 1,dm
           call setval(edge_coeffs(mgt(n)%nlevels,d),mac_beta,all=.true.)
        end do
     else 
        pd = mla%mba%pd(n)
        call multifab_build(cell_coeffs,mla%la(n),nc=1,ng=1)
        call init_cell_coeffs(mla,cell_coeffs,pd,coeffs_type)
        call cell_to_edge_coeffs(cell_coeffs,edge_coeffs(mgt(n)%nlevels,:))
        call multifab_destroy(cell_coeffs)
        if (n.gt.1) then
           do d = 1,dm
              call ml_edge_restriction(edge_coeffs(n-1,d),edge_coeffs(n,d),mla%mba%rr(n-1,:),d)
           end do
        end if
     end if

     pxa = ZERO
     pxb = ZERO
     if (n > 1) then
        xa = HALF*ref_ratio(n-1,:)*mgt(n)%dh(:,mgt(n)%nlevels)
        xb = HALF*ref_ratio(n-1,:)*mgt(n)%dh(:,mgt(n)%nlevels)
     else
        xa = ZERO
        xb = ZERO
     end if

     call stencil_fill_cc_all_mglevels(mgt(n), alpha, edge_coeffs, xa, xb, pxa, pxb, &
                                       stencil_order, domain_bc)

     call destroy(alpha(mgt(n)%nlevels))
     deallocate(alpha)

     do d = 1, dm
        call destroy(edge_coeffs(mgt(n)%nlevels,d))
     end do
     deallocate(edge_coeffs)
  end do

  if ( fabio ) then
    call fabio_ml_write(rh, ref_ratio(:,1), "rh-init_cc")
  end if

  snrm(2) = ml_norm_inf(rh,mla%mask)
! if ( parallel_IOProcessor() ) then
!    print *, 'RHS MAX NORM ', snrm(2)
! end if

! ****************************************************************************

  call ml_cc(mla, mgt, rh, full_soln, mla%mask, do_diagnostics)

! ****************************************************************************

  if ( fabio ) then
     call fabio_ml_write(full_soln, ref_ratio(:,1), "soln_cc")
  end if

! snrm(1) = ml_norm_l2(full_soln,ref_ratio,mla%mask)
! snrm(2) = ml_norm_inf(full_soln,mla%mask)
! if ( parallel_IOProcessor() ) then
!    print *, 'SOLUTION MAX NORM ', snrm(2)
!    print *, 'SOLUTION L2 NORM ', snrm(1)
! end if

! if ( parallel_IOProcessor() ) print *, 'MEMORY STATS'
! call print(multifab_mem_stats(),  " multifab before")
! call print(imultifab_mem_stats(), "imultifab before")
! call print(fab_mem_stats(),       "      fab before")
! call print(ifab_mem_stats(),      "     ifab before")
! call print(boxarray_mem_stats(),  " boxarray before")
! call print(boxassoc_mem_stats(),  " boxassoc before")
! call print(layout_mem_stats(),    "   layout before")

  do n = 1,nlevs
     call multifab_destroy(full_soln(n))
  end do

end subroutine t_cc_ml_multigrid
