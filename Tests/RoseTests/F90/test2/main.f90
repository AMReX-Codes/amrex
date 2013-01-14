program main

  implicit none

  type :: t
     integer, pointer :: d(:) => Null()
  end type t

  type(t) :: x
  type(t) :: d

  ! This is fine.
  allocate(x%d(3))  
  if (associated(x%d)) then
     print *, 'x%d is allocated and the size is', size(x%d)
     deallocate(x%d)
  end if

  ! ROSE does not like this.
  allocate(d%d(3))  
  if (associated(d%d)) then
     print *, 'd%d is allocated and the size is', size(d%d)
     deallocate(d%d)
  end if

end program main
