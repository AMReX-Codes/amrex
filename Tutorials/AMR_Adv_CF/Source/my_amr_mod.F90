module my_amr_module

  use amrex_famrcore_module

  implicit none

  type, public, extends(amrex_famrcore) :: my_amr
   contains
     procedure  ::            my_amr_build
     procedure  ::            my_amr_destroy
     generic    :: build   => my_amr_build
     generic    :: destroy => my_amr_destroy
  end type my_amr

contains

  subroutine my_amr_build (this)
    class(my_amr) :: this
    print *, 'in my_amr_build'
    call this%famrcore_build()
  end subroutine my_amr_build

  subroutine my_amr_destroy (this)
    class(my_amr) :: this
    print *, 'in my_amr_destroy'
    call this%famrcore_destroy()
  end subroutine my_amr_destroy

end module my_amr_module
