! Process a group of 1-d plotfiles from the dustcollapse problem
! and output the position of the interface as a function of time.

program fdustcollapse1d

  use f2kcli
  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use sort_d_module

  implicit none

  type(plotfile) pf
  integer :: unit
  integer :: i, j, ii
  integer :: f

  integer :: rr, r1

  integer :: cnt
  integer :: max_points

  real(kind=dp_t) :: dx(MAX_SPACEDIM)

  integer :: index

  real(kind=dp_t), pointer :: p(:,:,:,:)

  real(kind=dp_t), allocatable :: rcoord(:), dens(:)
  integer, allocatable :: isort(:)

  integer :: dens_comp

  real(kind=dp_t) :: dens_thres

  logical, allocatable :: imask(:)
  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  real(kind=dp_t) :: rmin

  integer :: narg, farg
  character(len=256) :: fname

  unit = unit_new()

  ! set the defaults
  dens_thres = 1.d8

  narg = command_argument_count()

  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('-d', '--dens_thres')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) dens_thres

     case default
        exit

     end select
     farg = farg + 1
  end do

  print *, '# density threshold = ', dens_thres


  do f = farg, narg

     call get_command_argument(f, value = fname)

     call build(pf, fname, unit)

     if (pf%dim /= 1) then
        print *, 'ERROR: not a 1-d file'
        stop
     endif

     rmin = pf%plo(1)

  
     ! density index
     dens_comp = plotfile_var_index(pf, "density")

     if (dens_comp < 0) then
        call bl_error("ERROR: variable(s) not defined")
     endif

     ! get the index bounds and dx for the coarse level.  Note, lo and hi are
     ! ZERO based indicies
     lo = lwb(plotfile_get_pd_box(pf, 1))
     hi = upb(plotfile_get_pd_box(pf, 1))

     dx = plotfile_get_dx(pf, 1)

     ! get the index bounds for the finest level
     flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
     fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

     ! compute the maximum number of zones, as if we were completely refined
     max_points = fhi(1) - flo(1) + 1


     ! imask will be set to false if we've already output the data.
     ! Note, imask is defined in terms of the finest level.  As we loop
     ! over levels, we will compare to the finest level index space to
     ! determine if we've already output here
     allocate(imask(flo(1):fhi(1)))

     imask(:) = .true.

     ! allocate storage for the data 
     allocate(rcoord(max_points))
     allocate(  dens(max_points))
     allocate( isort(max_points))

     rcoord(:) = ZERO
     dens(:)   = ZERO
     isort(:)  = ZERO

     ! loop over the data, starting at the finest grid, and if we haven't
     ! already store data in that grid location (according to imask),
     ! store it. 

     cnt = 0


     ! r1 is the factor between the current level grid spacing and the
     ! FINEST level
     r1  = 1

     do i = pf%flevel, 1, -1

        ! rr is the factor between the COARSEST level grid spacing and
        ! the current level
        rr = product(pf%refrat(1:i-1,1))

        do j = 1, nboxes(pf, i)

           ! read in the data 1 patch at a time
           call fab_bind(pf, i, j)

           lo = lwb(get_box(pf, i, j))
           hi = upb(get_box(pf, i, j))

           ! get a pointer to the current patch
           p => dataptr(pf, i, j)
        
           ! loop over all of the zones in the patch.  Here, we convert
           ! the cell-centered indices at the current level into the
           ! corresponding RANGE on the finest level, and test if we've
           ! stored data in any of those locations.  If we haven't then
           ! we store this level's data and mark that range as filled.
           do ii = lbound(p,dim=1), ubound(p,dim=1)

              if ( any(imask(ii*r1:(ii+1)*r1-1) ) ) then

                 cnt = cnt + 1

                 rcoord(cnt) = rmin + (ii + HALF)*dx(1)/rr
                 dens(cnt)   = p(ii,1,1,dens_comp)

                 imask(ii*r1:(ii+1)*r1-1) = .false.
                 
              end if

           enddo

           call fab_unbind(pf, i, j)
        end do

        ! adjust r1 for the next lowest level
        if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
     end do

     ! sort the data based on the coordinates
     call sort(rcoord(1:cnt),isort(1:cnt))


     ! loop over the solution, from r = 0 outward, and find the first
     ! place where the density drops below the threshold density
     index = -1
     do i = 1, cnt
        if (dens(isort(i)) < dens_thres) then
           index = i
           exit
        endif
     enddo

     if (index < 0) then
        print *, 'ERROR: density never fell below threshold'
        stop
     endif

     ! output
     print *, pf%tm, rcoord(index)

     ! clean-up
     deallocate(rcoord)
     deallocate(dens)
     deallocate(isort)
     deallocate(imask)

     call destroy(pf)
     
  enddo

end program fdustcollapse1d
