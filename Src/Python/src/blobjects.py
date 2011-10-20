"""BoxLib/Python object store for various F90 types.

For each BoxLib type in *types*, an object store is created along with
routines for adding and retrieving objects.  Objects in each store are
identified by an object id(oid).

The current implementation for a store is simply a fixed sized array.
"""

types = [ 'lmultifab', 'multifab', 'layout', 'boxarray' ]

# XXX: probably better to use arrays of pointers here!

module = '''\
module blobjects
  use multifab_module
  use layout_module
  integer, parameter :: PYBL_MAX_STORE = 256
  {stores}
contains
  {routines}
end module blobjects
'''

store = '''
type {type}_store
  integer :: oid = -1
  type({type}), pointer :: ptr
end type {type}_store
type({type}_store), save :: pybl_{type}_store(PYBL_MAX_STORE)
integer, save :: pybl_{type}_count = 0
'''

get = '''
subroutine pybl_{type}_get(oid,object)
  integer, intent(in) :: oid
  type({type}), pointer, intent(out) :: object

  object => pybl_{type}_store(oid)%ptr
end subroutine pybl_{type}_get'''

new = '''
subroutine pybl_{type}_new(oid,object)
  integer, intent(out) :: oid
  type({type}), pointer, intent(out) :: object

  pybl_{type}_count = pybl_{type}_count + 1
  oid = pybl_{type}_count

  allocate(pybl_{type}_store(oid)%ptr)
  object => pybl_{type}_store(oid)%ptr
end subroutine pybl_{type}_new
'''

stores = []
routines = []

for t in types:
    stores.append(store.format(type=t))
    routines.append(get.format(type=t))
    routines.append(new.format(type=t))

with open('blobjects.f90', 'w') as f:
    f.write(module.format(
        stores='\n'.join(stores),
        routines='\n'.join(routines)))
