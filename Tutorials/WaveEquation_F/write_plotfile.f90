module write_plotfile_module

  use ml_layout_module
  use multifab_module
  use fabio_module

  implicit none

contains
  
  subroutine write_plotfile(mla,data,istep_to_write,dx,time,prob_lo,prob_hi)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: data(mla%dim)
    integer        , intent(in   ) :: istep_to_write
    real(dp_t)     , intent(in   ) :: dx(mla%nlevel,mla%dim), time
    real(dp_t)     , intent(in   ) :: prob_lo(mla%dim), prob_hi(mla%dim)

    type(multifab) :: plotdata(mla%dim)

    integer           :: n, nlevs
    character(len=8)  :: sd_name
    character(len=20) :: plot_names(2)

    nlevs = mla%nlevel

    plot_names(1) = "Variable 1"
    plot_names(2) = "Variable 2"

    do n=1,nlevs
       call multifab_build(plotdata(n),mla%la(n),2,0)
       call multifab_copy_c(plotdata(n),1,data(n),1,2)
    end do

    write(unit=sd_name,fmt='("plt",i5.5)') istep_to_write

    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), sd_name, plot_names, &
                                   mla%mba%pd(1), prob_lo, prob_hi, time, dx(1,:))

    do n = 1,nlevs
       call multifab_destroy(plotdata(n))
    end do

  end subroutine write_plotfile

end module write_plotfile_module
