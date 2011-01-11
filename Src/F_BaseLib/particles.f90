module particle_module

  use bl_space
  use ml_layout_module

  implicit none

  type particle

     integer :: id   = -1
     integer :: cpu  = -1
     integer :: lev  = -1
     integer :: grd  = -1

     integer :: cell(MAX_SPACEDIM)

     double precision :: pos(MAX_SPACEDIM)

  end type particle

  interface print
     module procedure particle_print
  end interface print

  type particle_vector
     private
     integer :: size = 0
     integer :: bump = 2
     type(particle), pointer :: d(:) => NULL()
  end type particle_vector

  interface capacity
    module procedure particle_vector_capacity
  end interface capacity

  interface size
     module procedure particle_vector_size
  end interface size

  interface build
     module procedure particle_vector_build
  end interface build

  interface destroy
     module procedure particle_vector_destroy
  end interface destroy

  interface empty
     module procedure particle_vector_empty
  end interface empty

  interface at
     module procedure particle_vector_at
  end interface at
  interface set_at
     module procedure particle_vector_at_set
  end interface set_at

  interface push_back
     module procedure particle_vector_push_back
  end interface push_back

  interface pop_back
     module procedure particle_vector_pop_back
  end interface pop_back

  interface erase
     module procedure particle_vector_erase
  end interface erase

  interface back
     module procedure particle_vector_back
  end interface back

  interface reserve
    module procedure particle_vector_reserve
  end interface reserve

  interface clear
     module procedure particle_vector_clear
  end interface clear

  interface print
     module procedure particle_vector_print
  end interface print

contains

  subroutine particle_print(p)

    type(particle), intent(in) :: p

    print*, 'id   = ', p%id
    print*, 'cpu  = ', p%cpu
    print*, 'lev  = ', p%lev
    print*, 'grd  = ', p%grd
    print*, 'cell = ', p%cell
    print*, 'pos  = ', p%pos

  end subroutine particle_print

  subroutine particle_vector_build(d)
    type(particle_vector), intent(out) :: d
    allocate(d%d(d%size))
  end subroutine particle_vector_build

  subroutine particle_vector_destroy(d)
    type(particle_vector), intent(inout) :: d
    call clear(d)
  end subroutine particle_vector_destroy

  pure function particle_vector_empty(d) result(r)
    logical :: r
    type(particle_vector), intent(in) :: d
    r = (d%size == 0)
  end function particle_vector_empty

  pure function particle_vector_size(d) result(r)
    integer :: r
    type(particle_vector), intent(in) :: d
    r = d%size
  end function particle_vector_size

  pure function particle_vector_capacity(d) result(r)
    integer :: r
    type(particle_vector), intent(in) :: d
    r = size(d%d)
  end function particle_vector_capacity

  pure function particle_vector_at(d, i) result(r)
    type(particle) :: r
    integer, intent(in) :: i
    type(particle_vector), intent(in) :: d
    r = d%d(i)
  end function particle_vector_at

  subroutine particle_vector_at_set(d, i, v)
    type(particle),        intent(in   ) :: v
    integer,               intent(in   ) :: i
    type(particle_vector), intent(inout) :: d
    d%d(i) = v
  end subroutine particle_vector_at_set

  subroutine particle_vector_reserve(d, size)
    type(particle_vector), intent(inout) :: d
    integer,               intent(in   ) :: size
    type(particle), pointer :: np(:)
    if ( size <= particle_vector_capacity(d) ) return
    allocate(np(size))
    np(1:d%size) = d%d(1:d%size)
    if ( associated(d%d) ) then
       deallocate(d%d)
    end if
    d%d => np
  end subroutine particle_vector_reserve

  subroutine particle_vector_push_back(d,v)
    type(particle_vector), intent(inout) :: d
    type(particle),        intent(in   ) :: v
    if ( d%size >= particle_vector_capacity(d) ) then
       call particle_vector_reserve(d,max(d%size+1,d%size*d%bump))
    end if
    d%size      = d%size + 1
    d%d(d%size) = v
  end subroutine particle_vector_push_back

  subroutine particle_vector_pop_back(d)
    type(particle_vector), intent(inout) :: d
    d%size = d%size - 1
  end subroutine particle_vector_pop_back

  pure function particle_vector_back(d) result(r)
    type(particle_vector), intent(in) :: d
    type(particle)                    :: r
    r = d%d(d%size)
  end function particle_vector_back

  subroutine particle_vector_erase(d,i)
    type(particle_vector), intent(inout) :: d
    integer,               intent(in   ) :: i
    d%d(i:d%size-1) = d%d(i+1:d%size)
    d%size = d%size - 1
  end subroutine particle_vector_erase

  subroutine particle_vector_clear(d)
    type(particle_vector), intent(inout) :: d
    d%size = 0
    if ( associated(d%d) ) then
       deallocate(d%d)
    end if
  end subroutine particle_vector_clear

  subroutine particle_vector_print(d, str)
    type(particle_vector), intent(in) :: d
    character (len=*), intent(in), optional :: str
    integer :: i
    if ( present(str) ) print*, str
    if ( empty(d) ) then
       print*, '"Empty"'
    else
       do i = 1, d%size
          call print(d%d(i))
       end do
    end if
  end subroutine particle_vector_print

end module particle_module
