! Take 2 plotfiles as input and compare them point by point for
! differences.
!
! For files with identical grids, use fcompare instead.  This version
! can work on files with different grid structures -- it maps
! everything onto a uniform grid first.
!

program fmlcompare

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use sort_d_module
  use filler_module

  implicit none

  type(plotfile) :: pf_a, pf_b
  character (len=256) :: plotfile_a, plotfile_b
  integer :: unit_a, unit_b

  integer :: tlo(MAX_SPACEDIM), thi(MAX_SPACEDIM)

  integer :: flo_a(MAX_SPACEDIM), fhi_a(MAX_SPACEDIM)
  integer :: flo_b(MAX_SPACEDIM), fhi_b(MAX_SPACEDIM)

  real(kind=dp_t) :: dx_a(MAX_SPACEDIM), dx_b(MAX_SPACEDIM)

  integer :: n_a, n_b
  integer, allocatable :: ivar_b(:)

  real(kind=dp_t), allocatable :: aerror(:), rerror(:)
  integer :: norm

  integer :: narg, farg
  character (len=256) :: fname

  integer :: i, j, k
  integer :: ii, jj, kk

  real(kind=dp_t) :: x, y, z
  real(kind=dp_t) :: xmin, xmax, ymin, ymax, zmin, zmax

  integer :: comp(1)
  real(kind=dp_t), allocatable :: fab_a_2d(:,:,:),   fab_b_2d(:,:,:)
  real(kind=dp_t), allocatable :: fab_a_3d(:,:,:,:), fab_b_3d(:,:,:,:)


  !---------------------------------------------------------------------------
  ! process the command line arguments

  narg = command_argument_count()

  ! defaults
  norm = 0
  plotfile_a = ""
  plotfile_b = ""

  xmin = -1
  xmax = -1
  ymin = -1
  ymax = -1
  zmin = -1
  zmax = -1

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

     case ('-x','--xmin')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) xmin

     case ('-X','--xmax')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) xmax

     case ('-y','--ymin')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) ymin

     case ('-Y','--ymax')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) ymax

     case ('-z','--zmin')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) zmin

     case ('-Z','--zmax')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname, *) zmax


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

  
  ! check if they are the same dimensionality
  if (pf_a%dim /= pf_b%dim) then
     call bl_error("ERROR: plotfiles have different numbers of spatial dimensions")
  endif


  ! make sure the finest level dx's agree
  dx_a = plotfile_get_dx(pf_a, pf_a%flevel)
  dx_b = plotfile_get_dx(pf_b, pf_b%flevel)  


  if ((dx_a(1) /= dx_b(1)) .OR. &
       (pf_a%dim >= 2 .AND. dx_a(2) /= dx_b(2)) .OR. &
       (pf_a%dim == 3 .AND. dx_a(3) /= dx_b(3))) then
     call bl_error("ERROR: grid dx does not match")
  endif


  ! get the bounds of the finest level
  flo_a = lwb(plotfile_get_pd_box(pf_a, pf_a%flevel))
  fhi_a = upb(plotfile_get_pd_box(pf_a, pf_a%flevel))

  flo_b = lwb(plotfile_get_pd_box(pf_b, pf_b%flevel))
  fhi_b = upb(plotfile_get_pd_box(pf_b, pf_b%flevel))

  if ( (flo_a(1) /= flo_b(1) .OR. fhi_a(1) /= fhi_b(1)) .OR. &
      ((flo_a(2) /= flo_b(2) .OR. fhi_a(2) /= fhi_b(2)) .AND. pf_a%dim >= 2) .OR. &
      ((flo_a(3) /= flo_b(3) .OR. fhi_a(3) /= fhi_b(3)) .AND. pf_a%dim == 3) ) then
     call bl_error("ERROR: grids do not match")
  endif

  print *, 'file 1 finest grid: '
  print *, '  lo = ', flo_a(1:pf_a%dim)
  print *, '  hi = ', fhi_a(1:pf_a%dim)
  print *, ' '
  print *, 'file 2 finest grid: '
  print *, '  lo = ', flo_b(1:pf_b%dim)
  print *, '  hi = ', fhi_b(1:pf_b%dim)
  print *, ' '

  ! determine the limits in index space of the subdomain over which we
  ! will compare (tlo(:), thi(:))
  if (xmin >= ZERO) then
     do i = flo_a(1), fhi_a(1)
        x = pf_a%plo(1) + (i + HALF)*dx_a(1)
        if (x >= xmin) then
           tlo(1) = i
           exit
        endif
     enddo
  else
     tlo(1) = flo_a(1)
  endif

  if (xmax >= ZERO) then
     do i = flo_a(1), fhi_a(1)
        x = pf_a%plo(1) + (i + HALF)*dx_a(1)
        if (x >= xmax) then
           thi(1) = i-1
           exit
        endif
     enddo
  else
     thi(1) = fhi_a(1)
  endif

  if (ymin >= ZERO) then
     do j = flo_a(2), fhi_a(2)
        y = pf_a%plo(2) + (j + HALF)*dx_a(2)
        if (y >= ymin) then
           tlo(2) = j
           exit
        endif
     enddo
  else
     tlo(2) = flo_a(2)
  endif

  if (ymax >= ZERO) then
     do j = flo_a(2), fhi_a(2)
        y = pf_a%plo(2) + (j + HALF)*dx_a(2)
        if (y >= ymax) then
           thi(2) = j-1
           exit
        endif
     enddo
  else
     thi(2) = fhi_a(2)
  endif

  if (zmin >= ZERO) then
     do k = flo_a(3), fhi_a(3)
        z = pf_a%plo(3) + (k + HALF)*dx_a(3)
        if (z >= zmin) then
           tlo(3) = k
           exit
        endif
     enddo
  else
     tlo(3) = flo_a(3)
  endif

  if (zmax >= ZERO) then
     do k = flo_a(3), fhi_a(3)
        z = pf_a%plo(3) + (k + HALF)*dx_a(3)
        if (z >= zmax) then
           thi(3) = k-1
           exit
        endif
     enddo
  else
     thi(3) = fhi_a(3)
  endif


  print *, "comparison grid: "
  print *, '  lo = ', tlo(1:pf_a%dim)
  print *, '  hi = ', thi(1:pf_a%dim)
  print *, ' '

  ! allocate the space for a single, uniformly gridded variable
  select case (pf_a%dim)

  case (1)
     call bl_error("ERROR: 1-d not implemented")

  case (2)
     allocate(fab_a_2d(tlo(1):thi(1),tlo(2):thi(2),1))
     allocate(fab_b_2d(tlo(1):thi(1),tlo(2):thi(2),1))

  case (3)
     allocate(fab_a_3d(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),1))
     allocate(fab_b_3d(tlo(1):thi(1),tlo(2):thi(2),tlo(3):thi(3),1))
     
  end select
  

  ! check if they have the same number of variables
  if (pf_a%nvars /= pf_b%nvars) then
     print *, "WARNING: number of variables do not match"
  endif

  allocate(aerror(pf_a%nvars))
  allocate(rerror(pf_a%nvars))

  
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


  aerror(:) = ZERO
  rerror(:) = ZERO

  
  ! loop over the variables.  Take plotfile_a to be the one defining
  ! the list of variables, and bind them one-by-one.  Don't assume that
  ! the variables are in the same order in plotfile_b.

  do n_a = 1, pf_a%nvars

     n_b = ivar_b(n_a)
     if (n_b == -1) cycle


     select case (pf_a%dim)

     case (2)

        ! put the data onto a uniform grid that covers the specified
        ! subdomain
        comp(1) = n_a
        call blow_out_to_sub_fab(fab_a_2d, tlo, thi, pf_a, comp, pf_a%flevel)

        comp(1) = n_b
        call blow_out_to_sub_fab(fab_b_2d, tlo, thi, pf_b, comp, pf_b%flevel)

        do jj = tlo(2), thi(2)
           do ii = tlo(1), thi(1)

              if (norm == 0) then
                 aerror(n_a) = max(aerror(n_a), &
                      abs(fab_a_2d(ii,jj,1) - fab_b_2d(ii,jj,1)))
                    
                 rerror(n_a) = max(rerror(n_a), &
                      abs(fab_a_2d(ii,jj,1) - fab_b_2d(ii,jj,1)) / &
                      abs(fab_a_2d(ii,jj,1)))
              else
                 aerror(n_a) = aerror(n_a) + &
                      (abs(fab_a_2d(ii,jj,1) - fab_b_2d(ii,jj,1)))**norm

                 rerror(n_a) = rerror(n_a) + &
                      (abs(fab_a_2d(ii,jj,1) - fab_b_2d(ii,jj,1)) / &
                       abs(fab_a_2d(ii,jj,1)))**norm
              endif

           enddo
        enddo

     case (3)

        ! put the data onto a uniform grid that covers the specified
        ! subdomain
        comp(1) = n_a
        call blow_out_to_sub_fab(fab_a_3d, tlo, thi, pf_a, comp, pf_a%flevel)

        comp(1) = n_b
        call blow_out_to_sub_fab(fab_b_3d, tlo, thi, pf_b, comp, pf_b%flevel)

        do kk = tlo(3), thi(3)
           do jj = tlo(2), thi(2)
              do ii = tlo(1), thi(1)

                 if (norm == 0) then
                    aerror(n_a) = max(aerror(n_a), &
                         abs(fab_a_3d(ii,jj,kk,1) - fab_b_3d(ii,jj,kk,1)))
                    
                    rerror(n_a) = max(rerror(n_a), &
                         abs(fab_a_3d(ii,jj,kk,1) - fab_b_3d(ii,jj,kk,1)) / &
                         abs(fab_a_3d(ii,jj,kk,1)))
                 else
                    aerror(n_a) = aerror(n_a) + &
                         (abs(fab_a_3d(ii,jj,kk,1) - fab_b_3d(ii,jj,kk,1)))**norm

                    rerror(n_a) = rerror(n_a) + &
                         (abs(fab_a_3d(ii,jj,kk,1) - fab_b_3d(ii,jj,kk,1)) / &
                          abs(fab_a_3d(ii,jj,kk,1)))**norm
                 endif

              enddo
           enddo
        enddo

     end select

  enddo  ! variable loop


  ! normalize
  if (norm > 0) then

     do n_a = 1, pf_a%nvars
        aerror(n_a) = aerror(n_a)*product(dx_a(1:pf_a%dim))
        aerror(n_a) = aerror(n_a)**(ONE/real(norm,dp_t))

        rerror(n_a) = rerror(n_a)*product(dx_a(1:pf_a%dim))
        rerror(n_a) = rerror(n_a)**(ONE/real(norm,dp_t))
     enddo

  endif


  !------------------------------------------------------------------------
  ! print out the comparison report for this level
     
  do n_a = 1, pf_a%nvars
     print *, pf_a%names(n_a), aerror(n_a), rerror(n_a)
  enddo

end program fmlcompare
