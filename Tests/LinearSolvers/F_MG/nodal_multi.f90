subroutine t_nodal_ml_multigrid(mla, mgt, rh, coeffs_type, domain_bc, &
                                do_diagnostics,test, fabio, stencil_type)
  use BoxLib
  use nodal_stencil_module
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
  use ml_nodal_restriction_module
  use ml_prolongation_module
  use bndry_reg_module
  use ml_norm_module
  use coarsen_coeffs_module
  use nodal_stencil_fill_module
  use stencil_types_module
  use init_cell_coeffs_module

  use nodal_rhs_module

  implicit none

  type(ml_layout), intent(inout) :: mla
  type(mg_tower) , intent(inout) :: mgt(:)
  type(multifab) , intent(inout) :: rh(:)
  integer        , intent(in   ) :: coeffs_type
  integer        , intent(in   ) :: domain_bc(:,:)
  integer        , intent(in   ) :: do_diagnostics 
  integer        , intent(in   ) :: test
  logical        , intent(in   ) :: fabio
  integer        , intent(in   ) :: stencil_type

  type(lmultifab), allocatable :: fine_mask(:)
  type( multifab), allocatable :: cell_coeffs(:)
  type( multifab), allocatable :: full_soln(:)

  type(box)    :: pd
  type(layout) :: la
  integer      :: nlevs, n, i, dm
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

  allocate(full_soln(nlevs),fine_mask(nlevs))

  allocate(ref_ratio(nlevs-1,dm))
  do n = 1,nlevs-1
    ref_ratio(n,:) = mla%mba%rr(n,:)
  end do

  do n = nlevs, 1, -1

     pd = mla%mba%pd(n)
     if ( parallel_ioprocessor() ) &
        print *,'PD AT LEVEL ',n,' IS ',extent(pd)

     la = mla%la(n)
     call multifab_build(full_soln(n), la, 1, 1, nodal)
     call lmultifab_build(fine_mask(n), la, 1, 0, nodal)

     call setval(full_soln(n), ZERO,all=.true.)

  end do

  !! Fill coefficient array

  do n = nlevs, 1, -1
   
     allocate(cell_coeffs(mgt(n)%nlevels))

     la = get_layout(full_soln(n))

     pd = mla%mba%pd(n)

     call multifab_build(cell_coeffs(mgt(n)%nlevels), la, 1, 1)

     if (coeffs_type .eq. 1) then
         call setval(cell_coeffs(mgt(n)%nlevels), 1.0_dp_t, 1, all=.true.)
     else
         call init_cell_coeffs(mla,cell_coeffs(mgt(n)%nlevels),pd,coeffs_type)
     end if

     call stencil_fill_nodal_all_mglevels(mgt(n), cell_coeffs)

     call setval(fine_mask(n), val = .TRUE., all = .true.)

     if ( n < nlevs ) then
       call create_nodal_mask(n,fine_mask(n), &
                              mgt(n  )%mm(mgt(n  )%nlevels), &
                              mgt(n+1)%mm(mgt(n+1)%nlevels), &
                              mla)
     endif

     call multifab_destroy(cell_coeffs(mgt(n)%nlevels))
     deallocate(cell_coeffs)

  end do

  if ( fabio ) &
     call fabio_ml_write(rh, ref_ratio(:,1), "rh-init_nodal")

! ****************************************************************************

  call ml_nd(mla,mgt,rh,full_soln,fine_mask,do_diagnostics)

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
     call lmultifab_destroy(fine_mask(n))
  end do

contains

end subroutine t_nodal_ml_multigrid
