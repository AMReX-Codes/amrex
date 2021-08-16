module tagging_params_module

  use ISO_C_BINDING

  integer(C_INT), bind(C), save :: max_phierr_lev, max_phigrad_lev
  real(C_DOUBLE), bind(C), save :: phierr(0:15), phigrad(0:15)

end module tagging_params_module
