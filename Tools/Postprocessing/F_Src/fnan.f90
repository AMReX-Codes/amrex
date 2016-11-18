! report for each variable in a plotfile whether there is a NaN

program fnan

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use multifab_module
  use util_module

  implicit none

  type(plotfile) :: pf_a
  character (len=256) :: plotfile_a
  character (len=256) :: var
  integer :: unit_a
  logical :: zone_info

  real(kind=dp_t), pointer :: p_a(:,:,:,:)

  integer :: lo_a(MAX_SPACEDIM), hi_a(MAX_SPACEDIM)

  real(kind=dp_t) :: dx_a(MAX_SPACEDIM)

  integer :: nboxes_a

  integer :: n_a

  integer :: narg, farg
  character (len=256) :: fname

  integer :: i, j
  integer :: ii, jj, kk

  integer ir_a, ng

  integer :: itest

  real(kind=dp_t) :: pa

  integer :: dm
  type(box) :: bx_a

  logical :: has_nan
  

  ! process the command line arguments

  narg = command_argument_count()

  ! defaults
  plotfile_a = ""

  farg = 1
  call get_command_argument(farg, value = plotfile_a)
  farg = farg + 1

  
  !---------------------------------------------------------------------------
  ! build the plotfiles and do initial comparisons
  !---------------------------------------------------------------------------

  unit_a = unit_new()
  call build(pf_a, plotfile_a, unit_a)


  !---------------------------------------------------------------------------
  ! go variable-by-variable, then level-by-level and patch-by-patch and
  ! look for NaNs
  !---------------------------------------------------------------------------

  do n_a = 1, pf_a%nvars
     has_nan = .false.
     
     do i = 1, pf_a%flevel
        do j = 1, nboxes(pf_a, i)

           call fab_bind_comp_vec(pf_a, i, j, [n_a])
           p_a => dataptr(pf_a, i, j)

           call fab_contains_nan(p_a, volume(get_pbox(pf_a, i, j)), ir_a)

           if (ir_a == 1) then
              has_nan = .true.
           endif

           call fab_unbind(pf_a, i, j)
        enddo  ! boxes loop
     enddo

1000 format(1x, a24, ": ", a20)
     
     if (has_nan) then
        write(*,1000) trim(pf_a%names(n_a)), "has NaNs"
     else
        write(*,1000) trim(pf_a%names(n_a)), "clean"
     endif

  enddo
  
end program fnan
