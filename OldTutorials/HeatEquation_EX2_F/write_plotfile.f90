module write_plotfile_module

  use layout_module
  use multifab_module
  use fabio_module

  implicit none

contains
  
  subroutine write_plotfile(la,phi,istep,dx,time,prob_lo,prob_hi)

    type(layout)   , intent(in   ) :: la
    type(multifab) , intent(in   ) :: phi
    integer        , intent(in   ) :: istep
    real(dp_t)     , intent(in   ) :: dx,time
    real(dp_t)     , intent(in   ) :: prob_lo(phi%dim), prob_hi(phi%dim)

    ! local variables
    character(len=8)  :: plotfile_name
    character(len=20) :: variable_names(1)

    ! dimensioned as an array of size 1 for fabio_ml_multifab_write_d
    type(multifab) :: plotdata(1)

    ! dimensioned as an array of size 0 for fabio_ml_multifab_write_d
    integer :: rr(0)

    ! dimensioned as an array with size dm for fabio_ml_multifab_write_d
    real(dp_t) :: dx_vec(phi%dim)

    dx_vec(:) = dx

    variable_names(1) = "phi"

    ! build plotdata with 1 component and 0 ghost cells
    call multifab_build(plotdata(1),la,1,0)
    ! copy the state into plotdata
    call multifab_copy_c(plotdata(1),1,phi,1,1)

    ! define the name of the plotfile that will be written
    write(unit=plotfile_name,fmt='("plt",i5.5)') istep

    ! write the plotfile
    call fabio_ml_multifab_write_d(plotdata, rr, plotfile_name, variable_names, &
                                   la%lap%pd, prob_lo, prob_hi, time, dx_vec)

    ! make sure to destroy the multifab or you'll leak memory
    call multifab_destroy(plotdata(1))

  end subroutine write_plotfile

end module write_plotfile_module
