module tagging_params_module

  double precision, save :: phierr(0:2) = (/ 1.01d0, 1.1d0, 1.5d0 /)

  integer, save :: max_phierr_lev = 2

end module tagging_params_module
