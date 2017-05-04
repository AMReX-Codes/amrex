
module amrex_base_module

  use iso_c_binding

  use amrex_error_module
  use amrex_fort_module
  use amrex_box_module
  use amrex_boxarray_module
  use amrex_geometry_module
  use amrex_distromap_module
  use amrex_multifab_module
  use amrex_multifabutil_module
  use amrex_omp_module
  use amrex_parallel_module
  use amrex_parmparse_module
  use amrex_string_module
  use amrex_plotfile_module
  use amrex_physbc_module

  use amrex_bc_types_module
  use amrex_filcc_module

  use amrex_fab_module

  use mempool_module

contains

  subroutine amrex_fi_init () bind(c)
  end subroutine amrex_fi_init


  subroutine amrex_fi_finalize () bind(c)
  end subroutine amrex_fi_finalize

end module amrex_base_module
