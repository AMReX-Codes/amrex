"""BoxLib/Python object store for various F90 types."""

types = [ 'lmultifab', 'multifab', 'ml_layout', 'layout', 'boxarray' ]

module = '''\
module blobjects
  use multifab_module
  use layout_module
  use ml_layout_module
contains
  {routines}
end module blobjects
'''

get = '''
subroutine pybl_{type}_get(cptr,fptr)
  use iso_c_binding
  type(c_ptr), intent(in), value :: cptr
  type({type}), pointer, intent(out) :: fptr

  call c_f_pointer(cptr, fptr)
end subroutine pybl_{type}_get'''

new = '''
subroutine pybl_{type}_new(cptr,fptr)
  use iso_c_binding
  type(c_ptr), intent(out) :: cptr
  type({type}), pointer, intent(out) :: fptr

  allocate(fptr)
  cptr = c_loc(fptr)
end subroutine pybl_{type}_new
'''

# stores = []
routines = []

for t in types:
    # stores.append(store.format(type=t))
    routines.append(get.format(type=t))
    routines.append(new.format(type=t))

with open('blobjects.f90', 'w') as f:
    f.write(module.format(
        # stores='\n'.join(stores),
        routines='\n'.join(routines)))
