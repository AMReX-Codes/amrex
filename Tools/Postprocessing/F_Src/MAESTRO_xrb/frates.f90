program frates

  use network
  use rpar_indices
  use bl_space, only: MAX_SPACEDIM
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use bl_types
  use plotfile_module
  use multifab_module

  implicit none

  ! argument variables
  character(len=256) :: pltfile, outputfile
  logical :: use_tfromp

  ! f2kcli variables
  integer :: narg, farg
  character(len=256) :: fname

  ! local variables
  integer :: uin, dim, nlevs, i, j, ii, jj
  type(plotfile) :: pf
  type(layout) :: la
  type(boxarray) :: ba
  type(list_box) :: bl
  type(box) :: bx,pd
  type(multifab), allocatable :: rates(:)
  integer :: dens_comp, temp_comp, spec_comp
  integer, dimension(MAX_SPACEDIM) :: lo, hi
  integer, allocatable :: rr(:,:)
  real(kind=dp_t), pointer :: p(:,:,:,:), r(:,:,:,:)
  real(kind=dp_t) :: dens, t9, ymol(nspec)
  real(kind=dp_t), parameter :: T2T9 = 1e-9_dp_t
  real(kind=dp_t), allocatable :: rpar(:)
  character(len=20) :: plot_names(nrat)

  uin = unit_new()

  ! defaults
  pltfile =''
  outputfile = 'rates'
  use_tfromp = .false.
  
  ! parse arguments
  narg = command_argument_count()

  farg = 1
  do while (farg<=narg)
     call get_command_argument(farg,value=fname)

     select case(fname)

     case ('-i', '--input')
        farg = farg + 1
        call get_command_argument(farg,value=pltfile)
     case ('-o', '--output')
        farg = farg + 1
        call get_command_argument(farg,value=outputfile)
     case ('--tfromp')
        use_tfromp = .true.
     case default
        exit
     end select
     farg = farg + 1
  end do

  ! sanity check
  if (pltfile == '') then
     call print_usage()
     stop
  end if

  print *, 'working on pltfile: ', trim(pltfile)
  print *, 'saving to pltfile: ', trim(outputfile)
  if (use_tfromp) print *, 'using tfromp instead of tfromh'

  call network_init()
  allocate(rpar(n_rpar_comps))
  
  ! build the input plotfile
  call build(pf,pltfile,uin)

  nlevs = plotfile_nlevels(pf)
  dim = plotfile_dim(pf)

  allocate(rr(nlevs,dim),rates(nlevs))
  rr = plotfile_refrat(pf)

  do i = 1, nrat
     plot_names(i) = trim(reac_names(i))
  end do

  dens_comp = plotfile_var_index(pf,"density")
  if (use_tfromp) then
     temp_comp = plotfile_var_index(pf,"tfromp")
  else
     temp_comp = plotfile_var_index(pf,"tfromh")
  end if
  spec_comp = plotfile_var_index(pf,"X(" // trim(short_spec_names(1)) // ")")

  if (dens_comp < 0 .or. spec_comp < 0 .or. temp_comp < 0) then
     print *, dens_comp, temp_comp, spec_comp
     call bl_error("Variables not found")
  endif

  do i = 1, nlevs
     do j = 1, nboxes(pf,i)
        call push_back(bl,get_box(pf,i,j))
     end do

     call build(ba,bl)
     call build(la,ba,plotfile_get_pd_box(pf,i))
     call destroy(bl)
     call destroy(ba)

     ! build the multifab with 0 ghost cells and nrat components
     call multifab_build(rates(i),la,nrat,0)
  end do

  ! loop over the plotfile data starting at the finest
  do i = nlevs, 1, -1
     ! loop over each box at this level
     do j = 1, nboxes(pf,i)
        ! read in the data 1 patch at a time
        call fab_bind(pf,i,j)

        lo(1:dim) = lwb(get_box(pf,i,j))
        hi(1:dim) = upb(get_box(pf,i,j))

        p => dataptr(pf,i,j)
        r => dataptr(rates(i),j)

        do jj = lo(2), hi(2)
           do ii = lo(1), hi(1)
              
              dens = p(ii,jj,1,dens_comp)
              t9 = p(ii,jj,1,temp_comp)*T2T9
              ymol(:) = p(ii,jj,1,spec_comp:spec_comp+nspec-1)/aion(:)

              call make_rates(t9,dens,ymol,rpar)

              r(ii,jj,1,:) = rpar(irp_rates:irp_rates+nrat-1)
              
           end do
        end do

        call fab_unbind(pf,i,j)
        
     end do
  end do

  call fabio_ml_multifab_write_d(rates, rr(:,1), &
                                 trim(outputfile), &
                                 plot_names, plotfile_get_pd_box(pf,1), &
                                 pf%plo, pf%phi, plotfile_time(pf), &
                                 plotfile_get_dx(pf,1))
  call destroy(pf)


contains
  subroutine print_usage()
    implicit none
    
    print *,""
    print *, "This program takes a plotfile, calls a reaction network, and "
    print *, "dumps a new plotfile containing the reaction rates in each zone."
    print *, ""
    print *, "usage: "
    print *, " *frates* -i|--input <pltfile in> [-o|--output <pltfile out>]", &
         " [--tfromp]"    
    print *, ""
    print *, "    -i|--input: <pltfile in>"
    print *, "        Specify which plotfile to work on. (Required)"
    print *, "    -o|--output:"
    print *, "        Name of the out new plotfile to create. (Default: 'rates')"
    print *, "    --tfromp:"
    print *, "        Toggles the use of 'temperature' to be tfromp instead", &
         " of tfromh."
    print *, "        (Default: use tfromh)"
    print *, ""

  end subroutine print_usage
end program frates
