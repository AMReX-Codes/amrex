module return_test

  implicit none

contains

  AMREX_CUDA_FORT_DEVICE function fort_return_test (mallocd) result(returned) &
       bind(c,name='fort_return_test')
    use iso_c_binding, only : c_long
    integer(c_long) :: mallocd
    integer(c_long) :: returned

    mallocd  = 1
    returned = 1
  end function fort_return_test

end module return_test
