! average a sequence of plotfiles
program faverplot

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module
  use multifab_module

  implicit none

  type(plotfile) :: pf_a, pf
  character (len=256) :: plotfile_a, fname, outfile
  character (len=256) :: diffvar
  integer :: unit_a, unit

  real(kind=dp_t), pointer :: p_a(:,:,:,:), p_b(:,:,:,:), p(:,:,:,:)
  real(kind=dp_t), pointer :: mp(:,:,:,:)

  integer :: lo_a(MAX_SPACEDIM), hi_a(MAX_SPACEDIM)

  real(kind=dp_t) :: dx_a(MAX_SPACEDIM)

  real(kind=dp_t) :: files

  integer :: nboxes_a

  integer :: save_var_a, n_a

  logical :: all_variables_found

  integer :: narg, farg, nfiles

  integer :: f
  integer :: i, j
  integer :: ii, jj, kk

  integer ir_a, ir_b, ng

  integer :: itest

  real(kind=dp_t) :: pa, pd

  integer :: dm
  type(box) :: bx_a

  type(multifab), allocatable :: mf_array(:)
  type(layout) :: la
  type(boxarray) :: ba
  type(list_box) :: bl
  integer, allocatable :: ref_ratio(:)

  character(len=20), allocatable :: plot_names(:)

  logical :: do_ghost, gc_warn


  !$omp declare reduction(err_reduce: real(kind=dp_t): fort_error_reduce(omp_in, omp_out))


  !---------------------------------------------------------------------------
  ! process the command line arguments

   print *,' Caution:  no checking is done and current version requires diffvar be set'

  narg = command_argument_count()

  if ( narg .eq. 0 ) then
     print *, ''
     print *, 'Usage:'
     print *, '     faverage --diffvar <variable> <pltfile_list>'
     print *, ''
     print *, 'Description:'
     print *, '     This program takes a whitespace-separated list of plotfiles and'
     print *, '     returns the average for variable for each plotfile.'
     stop
  endif


  ! defaults
  plotfile_a = ""

  diffvar = ""

  allocate(plot_names(1))

  farg = 1
  do while (farg <= narg)
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('--diffvar')
        farg = farg + 1
        call get_command_argument(farg, value = diffvar)
        plot_names(1) = trim(diffvar)

     case default
        exit

     end select
     farg = farg + 1
  enddo

     nfiles = 0

     unit_a = unit_new()

  call get_command_argument(farg, value = plotfile_a)
  write(6,*)"Opening ",fname
  call build(pf_a, plotfile_a, unit_a)
  nfiles = nfiles + 1
  farg = farg + 1

  dm = pf_a%dim
  bx_a = plotfile_get_pd_box(pf_a, pf_a%flevel)
  save_var_a = -1
   
  do n_a = 1, pf_a%nvars

     if (.not. diffvar == "") then
        if (pf_a%names(n_a) == trim(diffvar)) then
           save_var_a = n_a
        endif
     endif

  enddo

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

  do i = 1, pf_a%flevel

     nboxes_a = nboxes(pf_a, i)

     do j = 1, nboxes_a

        ! make sure that the grids match
        bx_a = get_box(pf_a, i, j)

        lo_a = 1
        hi_a = 1
        lo_a(1:dm) = lwb(bx_a)
        hi_a(1:dm) = upb(bx_a)

        n_a = save_var_a


        call fab_bind_comp_vec(pf_a, i, j, (/ n_a /) )

        p_a => dataptr(pf_a, i, j)

        ! are we storing the diff?
        mp => dataptr(mf_array(i), j, get_box(mf_array(i), j))

        mp(:,:,:,1) = p_a(:,:,:,1) 

        call fab_unbind(pf_a, i, j)


     enddo  ! boxes loop

  enddo  ! level loop

do f = farg, narg

  unit = unit_new()

  call get_command_argument(f, value = fname)
  write(6,*)"Opening ",fname
  call build(pf, fname, unit)
  nfiles = nfiles + 1

  do i = 1, pf_a%flevel

     nboxes_a = nboxes(pf_a, i)

     do j = 1, nboxes_a

        ! make sure that the grids match
        bx_a = get_box(pf_a, i, j)

        lo_a = 1
        hi_a = 1
        lo_a(1:dm) = lwb(bx_a)
        hi_a(1:dm) = upb(bx_a)

        n_a = save_var_a

        call fab_bind_comp_vec(pf, i, j, (/ n_a /) )

        p => dataptr(pf, i, j)

        ! are we storing the diff?
        mp => dataptr(mf_array(i), j, get_box(mf_array(i), j))

        mp(:,:,:,1) = p(:,:,:,1) + mp(:,:,:,1)

        call fab_unbind(pf, i, j)


     enddo  ! boxes loop

  enddo  ! level loop



enddo !  file loop

  print *,'Average of ',nfiles, 'files'
  files = nfiles

  do i = 1, pf_a%flevel

     nboxes_a = nboxes(pf_a, i)

     do j = 1, nboxes_a

        mp => dataptr(mf_array(i), j, get_box(mf_array(i), j))
        mp(:,:,:,1) =  mp(:,:,:,1)/files

     enddo  ! boxes loop

  enddo  ! level loop

  outfile = trim("plt_average_" // diffvar)
  print *, 'Writing average to ',outfile

  if (save_var_a > 0) then
     call fabio_ml_multifab_write_d(mf_array, ref_ratio, &
                                    outfile, plot_names, &
                                    plotfile_get_pd_box(pf_a,1), &
                                    pf_a%plo, pf_a%phi, &
                                    pf_a%tm, &
                                    plotfile_get_dx(pf_a,1))
  endif


  call destroy(pf_a)

  deallocate(plot_names)

end program faverplot
