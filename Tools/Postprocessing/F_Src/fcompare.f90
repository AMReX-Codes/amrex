! Take 2 plotfiles as input and compare them zone by zone for differences.
!
! For the comparison to take place, the grids must be identical.
!

module util_module
  use bl_constants_module, only: dp_t
  implicit none
  type zone_t
     integer :: level
     integer :: box
     integer :: i
     integer :: j
     integer :: k
     real(kind=dp_t) :: max_abs_err
  end type zone_t
end module util_module

module reduction_mod
  implicit none

  integer, public :: norm
  integer, public :: n_a, n_b
  integer, public :: zone_info_var_a
  contains

    subroutine fort_error_reduce(omp_in, omp_out)
      use bl_constants_module, only: dp_t
      implicit none
      real(dp_t), dimension(:), intent(in) :: omp_in
      real(dp_t), dimension(:), intent(inout) :: omp_out
      if (norm == 0) then
        omp_out(n_a) = max(omp_out(n_a), omp_in(n_a))
      else
        omp_out(n_a) = omp_out(n_a) + omp_in(n_a)
      end if
    end subroutine fort_error_reduce

    subroutine fort_err_zone_reduce(omp_in, omp_out)
      use bl_constants_module, only: dp_t
      use util_module, only: zone_t
      implicit none
      type(zone_t), intent(in) :: omp_in
      type(zone_t), intent(inout) :: omp_out

      if (omp_in % max_abs_err > omp_out % max_abs_err) then
        omp_out % max_abs_err = omp_in % max_abs_err
        omp_out % level       = omp_in % level
        omp_out % box         = omp_in % box
        omp_out % i           = omp_in % i
        omp_out % j           = omp_in % j
        omp_out % k           = omp_in % k
      end if

    end subroutine fort_err_zone_reduce

    subroutine fort_init_err_zone(omp_priv)
      use util_module, only: zone_t
      implicit none
      type(zone_t), intent(out) :: omp_priv

      omp_priv % level = 0
      omp_priv % box = 0
      omp_priv % i = 0
      omp_priv % j = 0
      omp_priv % k = 0
      omp_priv % max_abs_err = 0.0
    end subroutine fort_init_err_zone

end module reduction_mod

program fcompare

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use multifab_module
  use util_module
  use reduction_mod, only: fort_error_reduce, &
                           fort_err_zone_reduce, &
                           fort_init_err_zone, &
                           norm, n_a, n_b, zone_info_var_a

  implicit none

  type(plotfile) :: pf_a, pf_b
  character (len=256) :: plotfile_a, plotfile_b
  character (len=256) :: diffvar, zone_info_var_name
  integer :: unit_a, unit_b
  logical :: zone_info

  real(kind=dp_t), pointer :: p_a(:,:,:,:), p_b(:,:,:,:), p(:,:,:,:)
  real(kind=dp_t), pointer :: mp(:,:,:,:)

  integer :: lo_a(MAX_SPACEDIM), hi_a(MAX_SPACEDIM)

  integer :: flo_b(MAX_SPACEDIM), fhi_b(MAX_SPACEDIM)

  real(kind=dp_t) :: dx_a(MAX_SPACEDIM), dx_b(MAX_SPACEDIM)

  integer :: nboxes_a, nboxes_b

  integer, allocatable :: ivar_b(:)
  integer :: save_var_a

  real(kind=dp_t), allocatable :: aerror(:), rerror(:), rerror_denom(:)

  logical, allocatable :: has_nan_a(:), has_nan_b(:)
  logical :: any_nans, all_variables_found

  integer :: narg, farg
  character (len=256) :: fname

  integer :: i, j
  integer :: ii, jj, kk

  integer ir_a, ir_b, ng

  integer :: itest

  real(kind=dp_t) :: global_error
  real(kind=dp_t) :: pa, pb, pd, aerr, rerr

  integer :: dm
  type(box) :: bx_a, bx_b

  type(multifab), allocatable :: mf_array(:)
  type(layout) :: la
  type(boxarray) :: ba
  type(list_box) :: bl
  integer, allocatable :: ref_ratio(:)

  character(len=20), allocatable :: plot_names(:)

  logical :: do_ghost, gc_warn

  type(zone_t) :: err_zone

  !$omp declare reduction(err_reduce: real(kind=dp_t): fort_error_reduce(omp_in, omp_out))

  !$omp declare reduction(err_zone_reduce: zone_t: fort_err_zone_reduce(omp_in, omp_out)) &
  !$omp   initializer(fort_init_err_zone(omp_priv))

  !---------------------------------------------------------------------------
  ! process the command line arguments

  narg = command_argument_count()

  ! defaults
  norm = 0
  plotfile_a = ""
  plotfile_b = ""

  diffvar = ""
  do_ghost = .false.
  zone_info = .false.
  zone_info_var_name = ""

  global_error = 0.d0

  allocate(plot_names(1))

  farg = 1
  do while (farg <= narg)
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('--infile1')
        farg = farg + 1
        call get_command_argument(farg, value = plotfile_a)

     case ('--infile2')
        farg = farg + 1
        call get_command_argument(farg, value = plotfile_b)

     case ('-n','--norm')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) norm

     case ('-g','--ghost')
        farg = farg + 1
        do_ghost = .true.

     case ('-z','--zone_info')
        farg = farg + 1
        call get_command_argument(farg, value = zone_info_var_name)
        zone_info = .true.

     case ('--diffvar')
        farg = farg + 1
        call get_command_argument(farg, value = diffvar)
        plot_names(1) = trim(diffvar)

     case default
        exit

     end select
     farg = farg + 1
  enddo

  if (len_trim(plotfile_a) == 0) then
     call get_command_argument(farg, value = plotfile_a)
     farg = farg + 1
  endif

  if (len_trim(plotfile_b) == 0) then
     call get_command_argument(farg, value = plotfile_b)
     farg = farg + 1
  endif

  if (len_trim(plotfile_a) == 0 .OR. len_trim(plotfile_b) == 0) then
     print *, " "
     print *, "Compare two plotfiles, zone by zone, to machine precision"
     print *, "and report the maximum absolute and relative errors for each"
     print *, "variable."
     print *, " "
     print *, "usage:"
     print *, "   fcompare [-g|--ghost] [-n|--norm num] [--diffvar var] [-z|--zone_info var] file1 file2"
     print *, " "
     print *, "optional arguments:"
     print *, "   -g|--ghost         : compare the ghost cells too (if stored)"
     print *, "   -n|--norm num      : what norm to use (default is 0 for inf norm)"
     print *, "   --diffvar var      : output a plotfile showing the differences for "
     print *, "                        variable var"
     print *, "   -z|--zone_info var : output the information for a zone corresponding"
     print *, "                        to the maximum error for the given variable"
     print *, " "
     stop
  endif

  !---------------------------------------------------------------------------
  ! build the plotfiles and do initial comparisons
  !---------------------------------------------------------------------------

  unit_a = unit_new()
  call build(pf_a, plotfile_a, unit_a)

  unit_b = unit_new()
  call build(pf_b, plotfile_b, unit_b)

  dm = pf_a%dim

  ! check if they are the same dimensionality
  if (pf_a%dim /= pf_b%dim) then
     call bl_error("ERROR: plotfiles have different numbers of spatial dimensions")
  endif

  ! check if they have the same number of levels
  if (pf_a%flevel /= pf_b%flevel) then
     call bl_error("ERROR: number of levels do not match")
  endif

  if (pf_a%flevel /= 0) then
     ! check if the finest domains are the same size
     bx_a = plotfile_get_pd_box(pf_a, pf_a%flevel)
     bx_b = plotfile_get_pd_box(pf_b, pf_b%flevel)
     
     if (.not. box_equal(bx_a, bx_b)) then
        call bl_error("ERROR: grids do not match")
     endif
  end if

  ! check if they have the same number of variables
  if (pf_a%nvars /= pf_b%nvars) then
     print *, " "
     print *, "WARNING: number of variables do not match"
  endif

  allocate(aerror(pf_a%nvars))
  allocate(rerror(pf_a%nvars))
  allocate(rerror_denom(pf_a%nvars))

  allocate(has_nan_a(pf_a%nvars))
  allocate(has_nan_b(pf_a%nvars))

  any_nans = .false.
  all_variables_found = .true.

  save_var_a = -1
  zone_info_var_a = -1

  ! in case the variables are not in the same order, figure out the
  ! mapping between pf_a and pf_b variables
  allocate(ivar_b(pf_a%nvars))

  do n_a = 1, pf_a%nvars

     ivar_b(n_a) = -1
     do n_b = 1, pf_b%nvars

        if (pf_a%names(n_a) == pf_b%names(n_b)) then
           ivar_b(n_a) = n_b
           exit
        endif

     enddo

     if (ivar_b(n_a) == -1) then
        print *, "WARNING: variable ", trim(pf_a%names(n_a)), &
                 " not found in plotfile 2"
        all_variables_found = .false.
     endif

     if (.not. diffvar == "") then
        if (pf_a%names(n_a) == trim(diffvar)) then
           save_var_a = n_a
        endif
     endif

     if (.not. zone_info_var_name == "") then
        if (pf_a%names(n_a) == trim(zone_info_var_name)) then
           zone_info_var_a = n_a
        endif
     endif

  enddo


  ! also print out, as a diagnostic, those variables in plotfile 1 that
  ! are not in plotfile 2
  do n_b = 1, pf_b%nvars
     itest = -1
     do n_a = 1, pf_a%nvars

        if (pf_a%names(n_a) == pf_b%names(n_b)) then
           itest = n_a
           exit
        endif

     enddo

     if (itest == -1) then
        print *, "WARNING: variable ", trim(pf_b%names(n_b)), &
                 " not found in plotfile 1"
        all_variables_found = .false.
     endif

  enddo


  !---------------------------------------------------------------------------
  ! create a multifab to store the difference for output, if desired
  !---------------------------------------------------------------------------

  if (save_var_a > 0) then

     allocate(mf_array(plotfile_nlevels(pf_a)))
     allocate(ref_ratio(plotfile_nlevels(pf_a)-1))

     ! define ref_ratio
     do i = 1, pf_a%flevel-1
        ref_ratio(i) = pf_a%refrat(i,1)
     enddo

     ! loop over levels and grids and define the boxes needed to build the
     ! multifab
     do i = 1, pf_a%flevel

        do j = 1, nboxes(pf_a, i)
           call push_back(bl, get_box(pf_a,i,j))
        enddo

        call build(ba,bl,sort=.false.)
        call layout_build_ba(la,ba,plotfile_get_pd_box(pf_a,1))

        ! destroy the list and boxarray so we start over next level
        call destroy(bl)
        call destroy(ba)

        ! create a new multifab with 0 ghost cells and 1 component
        call multifab_build(mf_array(i),la,1,0)

     enddo

  endif

  !---------------------------------------------------------------------------
  ! go level-by-level and patch-by-patch and compare the data
  !---------------------------------------------------------------------------

998 format(1x,a24,2x,a24,   2x,a24)
999 format(1x,70("-"))

  write (*,*) " "
  write (*,998) "variable name", "absolute error",  "relative error"
  write (*,998) "",              "(||A - B||)",     "(||A - B||/||A||)"
  write (*,999)

  gc_warn = .false.

  err_zone%max_abs_err = -1.d33

  do i = 1, pf_a%flevel

     aerror(:) = ZERO
     rerror(:) = ZERO
     rerror_denom(:) = ZERO

     has_nan_a(:) = .false.
     has_nan_b(:) = .false.

     ! make sure the dx's agree
     dx_a = 0.0_dp_t
     dx_b = 0.0_dp_t
     dx_a(1:dm) = plotfile_get_dx(pf_a, i)
     dx_b(1:dm) = plotfile_get_dx(pf_b, i)

     if ((dx_a(1) /= dx_b(1)) .OR. &
         (pf_a%dim >= 2 .AND. dx_a(2) /= dx_b(2)) .OR. &
         (pf_a%dim == 3 .AND. dx_a(3) /= dx_b(3))) then
        call bl_error("ERROR: grid dx does not match")
     endif

     ! make sure the number of boxes agree
     nboxes_a = nboxes(pf_a, i)
     nboxes_b = nboxes(pf_b, i)

     if (nboxes_a /= nboxes_b) then
        call bl_error("ERROR: number of boxes do not match")
     endif

     do j = 1, nboxes_a

        ! make sure that the grids match
        bx_a = get_box(pf_a, i, j)

        lo_a = 1
        hi_a = 1
        lo_a(1:dm) = lwb(bx_a)
        hi_a(1:dm) = upb(bx_a)

        bx_b = get_box(pf_b, i, j)

        if (.not. box_equal(bx_a, bx_b)) then
           call bl_error("ERROR: grids do not match")
        endif

        if (.not. nghost(pf_a, i, j) == nghost(pf_b, i, j)) then
           if (.not. gc_warn) then
              call bl_warn("WARNING: grids have different numbers of ghost cells")
              gc_warn = .true.
           endif
           ng = 0
        else
           if (do_ghost) then
              ng = nghost(pf_a, i, j)
           else
              ng = 0
           endif
        endif

        ! loop over the variables.  Take plotfile_a to be the one defining
        ! the list of variables, and bind them one-by-one.  Don't assume that
        ! the variables are in the same order in plotfile_b.

        do n_a = 1, pf_a%nvars

           n_b = ivar_b(n_a)
           if (n_b == -1) cycle

           call fab_bind_comp_vec(pf_a, i, j, (/ n_a /) )
           call fab_bind_comp_vec(pf_b, i, j, (/ n_b /) )

           p_a => dataptr(pf_a, i, j)
           p_b => dataptr(pf_b, i, j)

           ! are we storing the diff?
           if (n_a == save_var_a) then
              mp => dataptr(mf_array(i), j, get_box(mf_array(i), j))

              mp(:,:,:,1) = abs(p_a(:,:,:,1) - p_b(:,:,:,1))
           endif

           ! check for NaNs -- comparisons don't work when they are present
           ! note: regardless of do_ghost, this will check the ghostcells
           ! too if they are present.
           call fab_contains_nan(p_a, volume(get_pbox(pf_a, i, j)), ir_a)
           if (ir_a == 1) has_nan_a(n_a) = .true.

           call fab_contains_nan(p_b, volume(get_pbox(pf_b, i, j)), ir_b)
           if (ir_b == 1) has_nan_b(n_a) = .true.

           if (has_nan_a(n_a) .or. has_nan_b(n_a)) cycle

           !$omp parallel do collapse(2) default(none) &
           !$omp   shared(lo_a, hi_a, ng, p_a, p_b, norm, n_a, zone_info_var_a, i, j) &
           !$omp   private(pa, pb, pd) &
           !$omp   reduction(err_reduce:aerror, rerror, rerror_denom) &
           !$omp   reduction(err_zone_reduce:err_zone)
           do kk = lo_a(3)-ng, hi_a(3)+ng
              do jj = lo_a(2)-ng, hi_a(2)+ng
                 do ii = lo_a(1)-ng, hi_a(1)+ng

                    pa = abs(p_a(ii,jj,kk,1))
                    pb = abs(p_b(ii,jj,kk,1))
                    pd = abs(p_a(ii,jj,kk,1) - p_b(ii,jj,kk,1))

                    if (norm == 0) then
                       aerror(n_a) = max(aerror(n_a), pd)

                       rerror(n_a) = max(rerror(n_a), pd)
                       rerror_denom(n_a) = max(rerror_denom(n_a), pa)
                    else
                       aerror(n_a) = aerror(n_a) + pd**norm

                       rerror(n_a) = rerror(n_a) + pd**norm
                       rerror_denom(n_a) = rerror_denom(n_a) + pa**norm
                    endif

                    if (n_a == zone_info_var_a .and. pd > err_zone%max_abs_err) then
                       err_zone % max_abs_err = pd
                       err_zone % level = i
                       err_zone % box = j
                       err_zone % i = ii
                       err_zone % j = jj
                       err_zone % k = kk
                    endif

                 enddo
              enddo
           enddo

           call fab_unbind(pf_a, i, j)
           call fab_unbind(pf_b, i, j)

        enddo  ! variable loop

     enddo  ! boxes loop

     ! Normalize
     if (norm > 0) then

        do n_a = 1, pf_a%nvars
           aerror(n_a) = aerror(n_a)*product(dx_a(1:pf_a%dim))
           aerror(n_a) = aerror(n_a)**(ONE/real(norm,dp_t))

           ! since we are taking the ratio of two norms, no grid normalization
           ! is needed
           rerror(n_a) = (rerror(n_a)/rerror_denom(n_a))**(ONE/real(norm,dp_t))
        enddo

     else

        do n_a = 1, pf_a%nvars
           rerror(n_a) = rerror(n_a)/rerror_denom(n_a)
        enddo

     endif


     !------------------------------------------------------------------------
     ! print out the comparison report for this level
     !------------------------------------------------------------------------

1000 format(1x,"level = ", i2)
1001 format(1x,a24,2x,g24.10,2x,g24.10)
1002 format(1x,a24,2x,a50)
1003 format(1x,a24,2x,g24.10)

     write (*,1000) i

     do n_a = 1, pf_a%nvars
        if (ivar_b(n_a) == -1) then
           write (*,1002) pf_a%names(n_a), "< variable not present in both files > "

        else if (has_nan_a(n_a) .or. has_nan_b(n_a)) then
           write (*,1002) pf_a%names(n_a), "< NaN present > "

        else
           if (aerror(n_a) > 0.0d0) then
              aerr = min(max(aerror(n_a), 1.d-99), 1.d98)
           else
              aerr = 0.0d0
           endif

           if (rerror(n_a) > 0.0d0) then
              rerr = min(max(rerror(n_a), 1.d-99), 1.d98)
           else
              rerr = 0.0d0
           endif

           write (*,1001) pf_a%names(n_a), aerr, rerr
        endif
     enddo

     if (i == 1) then
        global_error = maxval(aerror(:))
     else
        global_error = max(global_error, maxval(aerror(:)))
     endif

     any_nans = any(has_nan_a .or. has_nan_b) .or. any_nans

  enddo  ! level loop

  if (save_var_a > 0) then
     call fabio_ml_multifab_write_d(mf_array, ref_ratio, &
                                    "diffs", plot_names, &
                                    plotfile_get_pd_box(pf_a,1), &
                                    pf_a%plo, pf_a%phi, &
                                    pf_a%tm, &
                                    plotfile_get_dx(pf_a,1))
  endif


  !------------------------------------------------------------------------
  ! print out the zone info for max abs error (if desired)
  !------------------------------------------------------------------------
  if (zone_info) then

     call fab_bind(pf_a, err_zone % level, err_zone % box)
     p => dataptr(pf_a, err_zone % level, err_zone % box)

     print *, " "
     print *, "maximum error in ", trim(zone_info_var_name)
     print *, "level = ", err_zone % level, " (i,j,k) = ", &
          err_zone % i, err_zone % j, err_zone % k
     do n_a = 1, pf_a%nvars
        write (*, 1003) trim(pf_a%names(n_a)), &
             p(err_zone % i, err_zone % j, err_zone % k, n_a)
     enddo

  endif

  call destroy(pf_a)
  call destroy(pf_b)

  deallocate(plot_names)

  deallocate(aerror)
  deallocate(rerror)
  deallocate(rerror_denom)

  deallocate(has_nan_a)
  deallocate(has_nan_b)

  deallocate(ivar_b)

  if (global_error == ZERO .and. .not. any_nans) then
     if (.not. all_variables_found) then
        print *, "WARNING: not all variables present in both files"
     endif
     print *, "PLOTFILES AGREE"
     call send_success_return_code()
  else
     call send_fail_return_code()
  endif

end program fcompare
