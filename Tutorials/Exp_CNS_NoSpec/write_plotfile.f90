module write_plotfile_module

  use layout_module
  use multifab_module
  use fabio_module

  implicit none

contains
  
  subroutine write_plotfile(data,istep,dx,time,prob_lo,prob_hi)

    type(multifab), intent(in   ) :: data
    integer       , intent(in   ) :: istep
    real(dp_t)    , intent(in   ) :: dx(data%dim),time
    real(dp_t)    , intent(in   ) :: prob_lo(data%dim), prob_hi(data%dim)

    character(len=2)               :: varname
    character(len=8)               :: plotfile_name
    character(len=20), allocatable :: variable_names(:)
    type(layout)                   :: la

    ! dimensioned as an array of size 1 for fabio_ml_multifab_write_d
    type(multifab) :: plotdata(1)

    ! dimensioned as an array of size 0 for fabio_ml_multifab_write_d
    integer :: rr(0), i, nc

    la = get_layout(data)

    nc = ncomp(data)

    allocate(variable_names(nc))

    do i = 1, nc
       write(unit=varname,fmt='(i2.2)') i
       variable_names(i) = 'Variable_' // varname
    end do

    call multifab_build(plotdata(1),la,nc,0)

    call multifab_copy_c(plotdata(1),1,data,1,nc)

    write(unit=plotfile_name,fmt='("plt",i5.5)') istep

    if ( parallel_IOProcessor() ) print*, 'Writing ' // plotfile_name // ' ...'

    call fabio_ml_multifab_write_d(plotdata, rr, plotfile_name, variable_names, &
                                   la%lap%pd, prob_lo, prob_hi, time, dx)

    call multifab_destroy(plotdata(1))

  end subroutine write_plotfile

end module write_plotfile_module
