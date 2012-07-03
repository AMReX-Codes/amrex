
module cutcell_module

  implicit none

  type cutcell
     !
     ! For now we just contain our cell index.
     !
     integer :: cell(3)
     double precision :: centroid(3)

  end type cutcell

  interface build
     module procedure cutcell_build
  end interface build

  interface print
     module procedure cutcell_print
  end interface print

  type cutcell_container
     private
     integer                :: size = 0
     type(cutcell), pointer :: d(:) => NULL()
  end type cutcell_container

  interface build
     module procedure cutcell_container_build
  end interface build

  interface destroy
     module procedure cutcell_container_destroy
  end interface destroy

  interface size
     module procedure cutcell_container_size
  end interface size

  interface empty
     module procedure cutcell_container_empty
  end interface empty

  interface print
     module procedure cutcell_container_print
  end interface print
  !
  ! Gives direct access to the underlying vector of cutcells.
  !
  ! You can change/update cutcells using this interface.
  !
  interface dataptr
    module procedure cutcell_container_dataptr
  end interface dataptr

contains

  subroutine cutcell_build(d,a,c)
    type(cutcell),    intent(inout) :: d
    integer,          intent(in   ) :: a(:)
    double precision, intent(in   ) :: c(:)
    d%cell = a
    d%centroid = c
  end subroutine cutcell_build

  subroutine cutcell_print(p)
    type(cutcell), intent(in) :: p
    print*, '  cell = ', p%cell, 'centroid = ', p%centroid
  end subroutine cutcell_print

  subroutine cutcell_container_build(d,sz)
    type(cutcell_container), intent(inout) :: d
    integer,                 intent(in   ) :: sz
    if ( .not. associated(d%d) ) then
       d%size = sz
       allocate(d%d(d%size))
    end if
  end subroutine cutcell_container_build

  subroutine cutcell_container_destroy(d)
    type(cutcell_container), intent(inout) :: d
    if ( associated(d%d) ) then
       deallocate(d%d)
       d%size  = 0
       d%d     => Null()
    end if
  end subroutine cutcell_container_destroy

  pure function cutcell_container_empty(d) result(r)
    logical :: r
    type(cutcell_container), intent(in) :: d
    r = (d%size == 0)
  end function cutcell_container_empty

  pure function cutcell_container_size(d) result(r)
    integer :: r
    type(cutcell_container), intent(in) :: d
    r = d%size
  end function cutcell_container_size

  subroutine cutcell_container_print(d, str)
    type(cutcell_container), intent(in)           :: d
    character (len=*),        intent(in), optional :: str

    integer :: i

    if ( present(str) ) print*, str
    if ( empty(d) ) then
       print*, '"Empty"'
    else
       do i = 1, d%size
          call print(d%d(i))
       end do
    end if
  end subroutine cutcell_container_print

  function cutcell_container_dataptr(d) result(r)
    type(cutcell), pointer :: r(:)
    type(cutcell_container), intent(in) :: d
    r => d%d
  end function cutcell_container_dataptr

end module cutcell_module

! program main

!   use cutcell_module

!   implicit none

!   integer, parameter :: ngrids = 16

!   logical, allocatable :: own(:)

!   type(cutcell_container), allocatable :: array_of_containers(:)

!   type(cutcell), pointer :: r(:)

!   type(cutcell) :: cell

!   integer i,j

!   allocate(own(ngrids))
!   allocate(array_of_containers(ngrids))

!   own = .false.
!   !
!   ! We only own grids 3, 7 and 9.
!   !
!   own(3) = .true.
!   own(7) = .true.
!   own(9) = .true.

!   do i = 1, ngrids
!      if (.not. empty(array_of_containers(i))) print*, 'cutcell container should be empty'

!      if (own(i)) then
!         !
!         ! Allocate space for "i" cutcells in array_of_containers(i).
!         !
!         call build(array_of_containers(i), 2*i)

!         if (empty(array_of_containers(i))) print*, 'cutcell container should NOT be empty'

!         r => dataptr(array_of_containers(i))

!         do j = 1,size(array_of_containers(i))
!            !
!            ! Populate the container.
!            !
!            call build(r(j), (/j,j,j/) , (/ 1.d0,1.d0,1.d0 /))
!            ! cell%cell(1) = j
!            ! cell%cell(2) = j
!            ! cell%cell(3) = j

!            ! r(j) = cell
!         end do

!         print*, 'cutcell container @ index = ', i

!         call print(array_of_containers(i))

!      end if

!   end do
  
!   do i = 1, ngrids
!      call destroy(array_of_containers(i))
!   end do

!   deallocate(own)
!   deallocate(array_of_containers)

! end program main
