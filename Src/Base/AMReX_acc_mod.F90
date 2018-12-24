
module amrex_acc_module

  implicit none

contains

  subroutine amrex_initialize_acc (id) bind(c,name='amrex_initialize_acc')
#ifdef AMREX_USE_ACC
    use openacc, only : acc_init, acc_set_device_num, acc_device_nvidia
#endif
    integer, intent(in), value :: id
#ifdef AMREX_USE_ACC
    call acc_init(acc_device_nvidia)
    call acc_set_device_num(id, acc_device_nvidia)
#endif
  end subroutine amrex_initialize_acc

  subroutine amrex_finalize_acc () bind(c,name='amrex_finalize_acc')
#ifdef AMREX_USE_ACC
    use openacc, only: acc_shutdown, acc_device_nvidia
    call acc_shutdown(acc_device_nvidia)
#endif
  end subroutine amrex_finalize_acc

end module amrex_acc_module
