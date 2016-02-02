module make_plotfile_module

  use bl_types
  use fabio_module
  use multifab_module
  use variables_module

  use chemistry_module, only : nspecies, spec_names
  use probin_module, only : nOutFiles, lUsingNFiles, prob_lo, prob_hi


  implicit none

  ! the total number of plot components
  integer, save :: n_plot_comps = 0
  integer, save :: icomp_rho=0, icomp_vel=0, icomp_pres=0, icomp_temp=0, icomp_Y=0

  private
  public :: n_plot_comps, get_plot_names, make_plotfile, init_plot_variables

contains

  function get_next_plot_index(num) result (next)

    ! return the next starting index for a plotfile quantity,
    ! and increment the counter of plotfile quantities by num
    integer :: num, next

    next = n_plot_comps + 1
    n_plot_comps = n_plot_comps + num

    return
  end function get_next_plot_index


  subroutine init_plot_variables()

    icomp_rho  = get_next_plot_index(1)
    icomp_vel  = get_next_plot_index(3)
    icomp_pres = get_next_plot_index(1)
    icomp_temp = get_next_plot_index(1)
    icomp_Y = get_next_plot_index(nspecies)

  end subroutine init_plot_variables


  subroutine get_plot_names(plot_names)
    character(len=20), intent(inout) :: plot_names(:)

    integer :: i

    plot_names(icomp_rho  ) = "density"
    plot_names(icomp_vel  ) = "x_vel"
    plot_names(icomp_vel+1) = "y_vel"
    plot_names(icomp_vel+2) = "z_vel"
    plot_names(icomp_pres) = "pressure"
    plot_names(icomp_temp) = "temperature"

    do i=1,nspecies
       plot_names(icomp_Y+i-1) = "Y("//trim(spec_names(i))//")"
    end do

  end subroutine get_plot_names


  subroutine make_plotfile(dirname, la, U, plot_names, time, dx, write_pf_time)

    character(len=*) , intent(in   ) :: dirname
    type(layout)     , intent(in   ) :: la
    type(multifab)   , intent(inout) :: U
    character(len=20), intent(in   ) :: plot_names(:)
    real(dp_t)       , intent(in   ) :: time, dx(U%dim)
    real(dp_t)       , intent(  out) :: write_pf_time

    ! dimensioned as an array of size 1 for fabio_ml_multifab_write_d
    type(multifab) :: plotdata(1), Q

    ! dimensioned as an array of size 0 for fabio_ml_multifab_write_d
    integer :: rr(0), prec
    real(dp_t) :: writetime1, writetime2

    prec = FABIO_DOUBLE

    call multifab_build(Q,la,nprim, 0)
    call multifab_setval(Q, 0.d0, all=.true.)

    call ctoprim(U, Q, 0)

    call multifab_build(plotdata(1),la,n_plot_comps,0)

    ! copy up to temperature
    call multifab_copy_c(plotdata(1),1,Q,1, qtemp)
    call multifab_copy_c(plotdata(1),icomp_Y, Q,qy1, nspecies)       

    call multifab_destroy(Q)

    if (parallel_IOProcessor()) then
       write(6,*) ''
       write(6,*) 'Writing state to plotfile ',trim(dirname)
    end if

    writetime1 = parallel_wtime()

    call fabio_ml_multifab_write_d(plotdata, rr, dirname, plot_names, &
         la%lap%pd, prob_lo, prob_hi, time, dx, &
         nOutFiles = nOutFiles, &
         lUsingNFiles = lUsingNFiles, prec = prec)

    writetime2 = parallel_wtime() - writetime1
    call parallel_reduce(writetime1, writetime2, MPI_MAX, proc=parallel_IOProcessorNode())
    if (parallel_IOProcessor()) then
       print*,'Time to write plotfile: ',writetime1,' seconds'
       print*,''
    end if

    write_pf_time = writetime1

    call multifab_destroy(plotdata(1))

  end subroutine make_plotfile

end module make_plotfile_module
