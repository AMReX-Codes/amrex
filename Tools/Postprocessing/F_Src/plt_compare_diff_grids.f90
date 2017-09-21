!
! Take 2 plotfiles as input and compare each level zone by zone for differences.
! The grids at each level may not be identical, but must have the same problem domain.
!

program fcompare

  use f2kcli
  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use sort_d_module

  implicit none

  type(plotfile) :: pf_a, pf_b
  character (len=256) :: plotfile_a, plotfile_b
  integer :: unit_a, unit_b

  real(kind=dp_t), pointer :: p_a(:,:,:,:), p_b(:,:,:,:)

  integer :: lo_a(MAX_SPACEDIM), hi_a(MAX_SPACEDIM)
  integer :: lo_b(MAX_SPACEDIM), hi_b(MAX_SPACEDIM)

  integer :: flo_a(MAX_SPACEDIM), fhi_a(MAX_SPACEDIM)
  integer :: flo_b(MAX_SPACEDIM), fhi_b(MAX_SPACEDIM)

  real(kind=dp_t) :: dx_a(MAX_SPACEDIM), dx_b(MAX_SPACEDIM)

  integer :: nboxes_a, nboxes_b

  integer :: n_a, n_b
  integer, allocatable :: ivar_b(:)

  real(kind=dp_t), allocatable :: aerror(:), rerror(:), rerror_denom(:)

  logical, allocatable :: has_nan_a(:), has_nan_b(:)

  integer :: norm

  integer :: narg, farg
  character (len=256) :: fname

  integer :: i, ja, jb
  integer :: ii, jj, kk

  integer ir_a, ir_b

  integer :: itest

  real(kind=dp_t) :: pa, pb, pd

  integer :: dm
  type(box) :: bx_a, bx_b, bx_i

  integer :: lo_i(MAX_SPACEDIM), hi_i(MAX_SPACEDIM)

  !---------------------------------------------------------------------------
  ! process the command line arguments

  narg = command_argument_count()

  ! defaults
  norm = 0
  plotfile_a = ""
  plotfile_b = ""

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

     case default
        exit

     end select
     farg = farg + 1
  enddo

  if (len_trim(plotfile_a) == 0 .OR. len_trim(plotfile_b) == 0) then
     print *, " "
     print *, "Compare two plotfiles, zone by zone, to machine precision"
     print *, "and report the maximum absolute and relative errors for each"
     print *, "variable."
     print *, " "
     print *, "usage:"
     print *, "   fcompare --infile1 file1 --infile2 file2"
     print *, " "
     stop
  endif

  !---------------------------------------------------------------------------
  ! build the plotfiles and do initial comparisons
  
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

  ! check if the problem domains are the same size at each level
  do i = 1,  pf_a%flevel
     bx_a = plotfile_get_pd_box(pf_a, i)
     flo_a = 1
     fhi_a = 1
     flo_a(1:dm) = lwb(bx_a)
     fhi_a(1:dm) = upb(bx_a)

     bx_b = plotfile_get_pd_box(pf_b, i)
     flo_b(1:dm) = lwb(bx_b)
     fhi_b(1:dm) = upb(bx_b)

     if ( (flo_a(1) /= flo_b(1) .OR. fhi_a(1) /= fhi_b(1)) .OR. &
         ((flo_a(2) /= flo_b(2) .OR. fhi_a(2) /= fhi_b(2)) .AND. pf_a%dim >= 2) .OR. &
         ((flo_a(3) /= flo_b(3) .OR. fhi_a(3) /= fhi_b(3)) .AND. pf_a%dim == 3) ) then
        call bl_error("ERROR: problem domains do not match")
     endif
  enddo

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
     endif

  enddo

  !---------------------------------------------------------------------------
  ! go level-by-level and patch-by-patch and compare the data
  
998 format(1x,a24,2x,a20,   2x,a20)
999 format(1x,70("-"))

  write (*,*) " "
  write (*,998) "variable name", "absolute error", "relative error"
  write (*,999)

   do i = 1,  pf_a%flevel

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

     nboxes_a = nboxes(pf_a, i)
     nboxes_b = nboxes(pf_b, i)
     
     do ja = 1, nboxes_a

        bx_a = get_box(pf_a, i, ja)
        lo_a = 1
        hi_a = 1
        lo_a(1:dm) = lwb(bx_a)
        hi_a(1:dm) = upb(bx_a)

        do jb = 1, nboxes_b
           bx_b = get_box(pf_b, i, jb)
           lo_b = 1
           hi_b = 1
           lo_b(1:dm) = lwb(bx_b)
           hi_b(1:dm) = upb(bx_b)

           if (.not. box_intersects(bx_a,bx_b)) cycle

           bx_i = box_intersection(bx_a,bx_b)
           lo_i = 1
           hi_i = 1
           lo_i(1:dm) = lwb(bx_i)
           hi_i(1:dm) = upb(bx_i)

           ! loop over the variables.  Take plotfile_a to be the one defining
           ! the list of variables, and bind them one-by-one.  Don't assume that
           ! the variables are in the same order in plotfile_b.

           do n_a = 1, pf_a%nvars

              n_b = ivar_b(n_a)
              if (n_b == -1) cycle

              call fab_bind_comp_vec(pf_a, i, ja, (/ n_a /) )
              call fab_bind_comp_vec(pf_b, i, jb, (/ n_b /) )

              p_a => dataptr(pf_a, i, ja)
              p_b => dataptr(pf_b, i, jb)
              
              ! check for NaNs -- comparisons don't work when they are present
              call fab_contains_nan(p_a, &
                   (hi_a(3)-lo_a(3)+1)* &
                   (hi_a(2)-lo_a(2)+1)* &
                   (hi_a(1)-lo_a(1)+1),  ir_a)
              
              if (ir_a == 1) has_nan_a(n_a) = .true.

              call fab_contains_nan(p_b, &
                   (hi_b(3)-lo_b(3)+1)* &
                   (hi_b(2)-lo_b(2)+1)* &
                   (hi_b(1)-lo_b(1)+1),  ir_b)
              
              if (ir_b == 1) has_nan_b(n_a) = .true.
              
              if (has_nan_a(n_a) .or. has_nan_b(n_a)) cycle
              

              do kk = lo_i(3), hi_i(3)
                 do jj = lo_i(2), hi_i(2)
                    do ii = lo_i(1), hi_i(1)

                       pa = abs(p_a(ii,jj,kk,1))
                       pb = abs(p_b(ii,jj,kk,1))
                       pd = abs(p_a(ii,jj,kk,1) - p_b(ii,jj,kk,1))

                       if (norm == 0) then

                          aerror(n_a) = max(aerror(n_a),pd)
                          rerror(n_a) = max(rerror(n_a), pd)
                          rerror_denom(n_a) = max(rerror_denom(n_a), pa)

                       else

                          aerror(n_a) = aerror(n_a) + pd**norm

                          rerror(n_a) = rerror(n_a) + pd**norm
                          rerror_denom(n_a) = rerror_denom(n_a) + pa**norm

                       endif
                       
                    enddo
                 enddo
              enddo

              call fab_unbind(pf_a, i, ja)
              call fab_unbind(pf_b, i, jb)

           enddo  ! variable loop

        end do ! boxes loop b

     enddo  ! boxes loop a

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
     
1000 format(1x,"level = ", i2)
1001 format(1x,a24,2x,g20.10,2x,g20.10)
1002 format(1x,a24,2x,a42)

     write (*,1000) i

     do n_a = 1, pf_a%nvars
        if (ivar_b(n_a) == -1) then
           write (*,1002) pf_a%names(n_a), "< variable not present in both files > "
        else if (has_nan_a(n_a) .or. has_nan_b(n_a)) then
           write (*,1002) pf_a%names(n_a), "< NaN present > "
        else
           write (*,1001) pf_a%names(n_a), aerror(n_a), rerror(n_a)
        endif
     enddo

   enddo  ! level loop

end program fcompare
