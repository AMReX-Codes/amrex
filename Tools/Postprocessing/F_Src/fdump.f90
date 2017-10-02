! simply dump the data and zone indices for a specific variable out to
! the screen

program fdump

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use multifab_module
  use box_util_module

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

  integer :: n_a, n

  integer :: narg, farg
  character (len=256) :: fname

  integer :: i, j
  integer :: ii, jj, kk

  integer ir_a, ng

  integer :: itest

  real(kind=dp_t) :: pa

  character(len=256) :: minval_str, maxval_str
  real(kind=dp_t) :: minval, maxval

  character(len=256) :: idx_str
  integer :: idx(3)

  character(len=256) :: lev_str
  integer :: level

  integer :: dm
  type(box) :: bx_a

  minval = -1.d200
  maxval = 1.d200

  idx(:) = -1

  level = -1

  !---------------------------------------------------------------------------
  ! process the command line arguments

  narg = command_argument_count()

  ! defaults
  plotfile_a = ""

  var = ""

  farg = 1
  do while (farg <= narg)
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('--var')
        farg = farg + 1
        call get_command_argument(farg, value = var)

     ! Only print out data larger than this value.
     case ('--min')
        farg = farg + 1
        call get_command_argument(farg, value = minval_str)
        read(minval_str, *) minval

     ! Only print out data smaller than this value.
     case ('--max')
        farg = farg + 1
        call get_command_argument(farg, value = maxval_str)
        read(maxval_str, *) maxval

     ! Only print out data with this zone index.
     case ('--i')
        farg = farg + 1
        call get_command_argument(farg, value = idx_str)
        read(idx_str, *) idx(1)

     case ('--j')
        farg = farg + 1
        call get_command_argument(farg, value = idx_str)
        read(idx_str, *) idx(2)

     case ('--k')
        farg = farg + 1
        call get_command_argument(farg, value = idx_str)
        read(idx_str, *) idx(3)

     case ('--level')
        farg = farg + 1
        call get_command_argument(farg, value = lev_str)
        read(lev_str, *) level

     case default
        exit

     end select
     farg = farg + 1
  enddo

  call get_command_argument(farg, value = plotfile_a)
  farg = farg + 1

  
  !---------------------------------------------------------------------------
  ! build the plotfiles and do initial comparisons
  !---------------------------------------------------------------------------

  unit_a = unit_new()
  call build(pf_a, plotfile_a, unit_a)

  n_a = 0
  do n = 1, pf_a%nvars
     if (pf_a%names(n) == trim(var)) then
        n_a = n
        exit
     endif
  enddo
  
  dm = pf_a%dim

  !---------------------------------------------------------------------------
  ! go level-by-level and patch-by-patch and output the data
  !---------------------------------------------------------------------------

  do i = 1, pf_a%flevel

     if (level >= 0 .and. level /= i) cycle

     nboxes_a = nboxes(pf_a, i)
     do j = 1, nboxes_a

        bx_a = get_box(pf_a, i, j)                                                         
                                                                                           
        lo_a = 1                                                                           
        hi_a = 1                                                                           
        lo_a(1:dm) = lwb(bx_a)                                                             
        hi_a(1:dm) = upb(bx_a)        
        
        ng = nghost(pf_a, i, j)
        
        call fab_bind_comp_vec(pf_a, i, j, [n_a])

        p_a => dataptr(pf_a, i, j)

        call fab_contains_nan(p_a, volume(get_pbox(pf_a, i, j)), ir_a)
        if (ir_a == 1) then
           print *, "contains NaN"
        endif
        
        do kk = lo_a(3)-ng, hi_a(3)+ng
           if (idx(3) >= 0 .and. kk /= idx(3)) cycle
           do jj = lo_a(2)-ng, hi_a(2)+ng
              if (idx(2) >= 0 .and. jj /= idx(2)) cycle
              do ii = lo_a(1)-ng, hi_a(1)+ng
                 if (idx(1) >= 0 .and. ii /= idx(1)) cycle

                 if (p_a(ii,jj,kk,1) >= minval .and. p_a(ii,jj,kk,1) <= maxval) then
                    print *, ii, jj, kk, p_a(ii,jj,kk,1)
                 end if
              enddo
           enddo
        enddo

        call fab_unbind(pf_a, i, j)
     enddo  ! boxes loop
  enddo

end program fdump
