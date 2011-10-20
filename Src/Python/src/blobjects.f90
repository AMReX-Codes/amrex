module blobjects
  use multifab_module
  use layout_module
  integer, parameter :: PYBL_MAX_STORE = 256
  
type lmultifab_store
  integer :: oid = -1
  type(lmultifab), pointer :: ptr
end type lmultifab_store
type(lmultifab_store), save :: pybl_lmultifab_store(PYBL_MAX_STORE)
integer, save :: pybl_lmultifab_count = 0


type multifab_store
  integer :: oid = -1
  type(multifab), pointer :: ptr
end type multifab_store
type(multifab_store), save :: pybl_multifab_store(PYBL_MAX_STORE)
integer, save :: pybl_multifab_count = 0


type layout_store
  integer :: oid = -1
  type(layout), pointer :: ptr
end type layout_store
type(layout_store), save :: pybl_layout_store(PYBL_MAX_STORE)
integer, save :: pybl_layout_count = 0


type boxarray_store
  integer :: oid = -1
  type(boxarray), pointer :: ptr
end type boxarray_store
type(boxarray_store), save :: pybl_boxarray_store(PYBL_MAX_STORE)
integer, save :: pybl_boxarray_count = 0

contains
  
subroutine pybl_lmultifab_get(oid,object)
  integer, intent(in) :: oid
  type(lmultifab), pointer, intent(out) :: object

  object => pybl_lmultifab_store(oid)%ptr
end subroutine pybl_lmultifab_get

subroutine pybl_lmultifab_new(oid,object)
  integer, intent(out) :: oid
  type(lmultifab), pointer, intent(out) :: object

  pybl_lmultifab_count = pybl_lmultifab_count + 1
  oid = pybl_lmultifab_count

  allocate(pybl_lmultifab_store(oid)%ptr)
  object => pybl_lmultifab_store(oid)%ptr
end subroutine pybl_lmultifab_new


subroutine pybl_multifab_get(oid,object)
  integer, intent(in) :: oid
  type(multifab), pointer, intent(out) :: object

  object => pybl_multifab_store(oid)%ptr
end subroutine pybl_multifab_get

subroutine pybl_multifab_new(oid,object)
  integer, intent(out) :: oid
  type(multifab), pointer, intent(out) :: object

  pybl_multifab_count = pybl_multifab_count + 1
  oid = pybl_multifab_count

  allocate(pybl_multifab_store(oid)%ptr)
  object => pybl_multifab_store(oid)%ptr
end subroutine pybl_multifab_new


subroutine pybl_layout_get(oid,object)
  integer, intent(in) :: oid
  type(layout), pointer, intent(out) :: object

  object => pybl_layout_store(oid)%ptr
end subroutine pybl_layout_get

subroutine pybl_layout_new(oid,object)
  integer, intent(out) :: oid
  type(layout), pointer, intent(out) :: object

  pybl_layout_count = pybl_layout_count + 1
  oid = pybl_layout_count

  allocate(pybl_layout_store(oid)%ptr)
  object => pybl_layout_store(oid)%ptr
end subroutine pybl_layout_new


subroutine pybl_boxarray_get(oid,object)
  integer, intent(in) :: oid
  type(boxarray), pointer, intent(out) :: object

  object => pybl_boxarray_store(oid)%ptr
end subroutine pybl_boxarray_get

subroutine pybl_boxarray_new(oid,object)
  integer, intent(out) :: oid
  type(boxarray), pointer, intent(out) :: object

  pybl_boxarray_count = pybl_boxarray_count + 1
  oid = pybl_boxarray_count

  allocate(pybl_boxarray_store(oid)%ptr)
  object => pybl_boxarray_store(oid)%ptr
end subroutine pybl_boxarray_new

end module blobjects
