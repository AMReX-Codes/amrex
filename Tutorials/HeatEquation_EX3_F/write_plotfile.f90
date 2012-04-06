module write_plotfile_module

  use layout_module
  use multifab_module
  use fabio_module

  implicit none

contains
  
  subroutine write_plotfile(nlevs,ref_ratio,la,data,istep,dx,time,prob_lo,prob_hi)

    integer        , intent(in   ) :: nlevs
    integer        , intent(in   ) :: ref_ratio(:)
    type(layout)   , intent(in   ) :: la(:)
    type(multifab) , intent(in   ) :: data(:)
    integer        , intent(in   ) :: istep
    real(dp_t)     , intent(in   ) :: dx(:),time
    real(dp_t)     , intent(in   ) :: prob_lo(data(1)%dim), prob_hi(data(1)%dim)

    ! local variables
    character(len=8)  :: plotfile_name
    character(len=20) :: variable_names(1)
    integer           :: n

    ! dimensioned as an array of size 1 for fabio_ml_multifab_write_d
    type(multifab) :: plotdata(nlevs)

    ! dimensioned as an array with size dm for fabio_ml_multifab_write_d
    real(dp_t) :: dx_vec(nlevs,data(1)%dim)

    do n = 1,nlevs
       dx_vec(n,:) = dx(n) 
    end do

    variable_names(1) = "Variable 1"

    ! build plotdata with 1 component and 0 ghost cells
    ! copy the state into plotdata
    do n = 1,nlevs
       call multifab_build(plotdata(n),la(n),multifab_ncomp(data(1)),0)
       call multifab_copy_c(plotdata(n),1,data(n),1,multifab_ncomp(data(1)))
    end do

    ! define the name of the plotfile that will be written
    write(unit=plotfile_name,fmt='("plt",i5.5)') istep

    ! write the plotfile
    call fabio_ml_multifab_write_d(plotdata, ref_ratio, plotfile_name, variable_names, &
                                   la(1)%lap%pd, prob_lo, prob_hi, time, dx_vec(1,:))

    ! make sure to destroy the multifab or you'll leak memory
    do n = 1,nlevs
       call multifab_destroy(plotdata(n))
    end do

  end subroutine write_plotfile

end module write_plotfile_module
