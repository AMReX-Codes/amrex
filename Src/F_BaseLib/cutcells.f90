
module cutcell_module

  use bl_error_module

  implicit none

  ! for flagarray:
  integer, parameter, public :: FLOW        = -9
  integer, parameter, public :: IRRFLOW     = -3
  integer, parameter, public :: SOLID       = -8
  ! cut cells have index >=1

  type cutcell
     !
     ! For now we just contain our cell index.
     !
     integer          :: cutNumber        ! its index in the cutcell list
     integer          :: cutDim           ! Dimension of cutcell         
     integer          :: cutIndex(3)      ! its Cartesian index (i,j,k)
     double precision :: centroid(3)      ! centroid coordinates (x,y,z)
     double precision :: volume           ! volume of cell
     double precision :: bdyArea          ! area of boundary face (the actual cut face)
     double precision :: bdyNormal(3)     ! normalized normal vector for boundary face
     double precision :: bdyCentroid(3)   ! centroid of boundary face
     integer          :: splitIndex       ! -99999 if unsplit
                                          ! if split index points to split kids
     logical          :: FaceExist(3,2)  ! 1st comp: x,y,z face
                                          ! 2nd comp: lo,hi face
                                          ! says whether the face of the cutcell exists or not
     double precision :: FaceCentroid(3,2,3)  ! 1st and 2nd comp: as for face_exist
                                               ! 3rd component: actual centroid coord
     double precision :: FaceArea(3,2)   ! 1st and 2nd comp as for face_centroid
                                          ! area of face

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

  subroutine cutcell_build(dcut,aaa,aab,a,b,c,e,f,g,h,arr_a,arr_b,arr_c)
    type(cutcell),    intent(inout) :: dcut
    integer,          intent(in   ) :: aaa    ! cutNumber
    integer,          intent(in   ) :: aab    ! cutDim
    integer,          intent(in   ) :: a(:)   ! cutIndex
    double precision, intent(in   ) :: b(:)   ! centroid
    double precision, intent(in   ) :: c      ! volume
    double precision, intent(in   ) :: e      ! bdyArea
    double precision, intent(in   ) :: f(:)   ! bdyNormal
    double precision, intent(in   ) :: g(:)   ! bdyCentroid
    integer,          intent(in   ) :: h      ! splitIndex
    logical,          intent(in   ) :: arr_a(:,:)   ! FaceExist(3,2) 
    double precision, intent(in   ) :: arr_b(:,:,:) ! FaceCentroid(3,2,3)
    double precision, intent(in   ) :: arr_c(:,:)   ! FaceArea(3,2)

    integer i,j,k
    double precision :: my_z_cent

    if (aab == 3) then   ! for 3d problem

       dcut%cutNumber = aaa
       dcut%cutDim = aab
       dcut%cutIndex  = a
       dcut%centroid  = b
       dcut%volume    = c
       dcut%bdyArea   = e
       dcut%bdyNormal = f
       dcut%bdyCentroid = g
       dcut%splitIndex = h

       do j=1,2
          do i=1,3
             dcut%FaceExist(i,j) = arr_a(i,j)
             dcut%FaceArea(i,j)  = arr_c(i,j)
          end do
       end do
       do k=1,3
          do j=1,2
             do i=1,3
                dcut%FaceCentroid(i,j,k) = arr_b(i,j,k)
             end do
          end do
       end do

    elseif (aab == 2) then ! 2d problem

       dcut%cutNumber = aaa
       dcut%cutDim = aab
       if (a(3) .ne. 0) call bl_error("cutcells.f90: cutIndex(3) not equal 0!")
       dcut%cutIndex  = a
       dcut%centroid  = b
       dcut%volume    = c     ! if we assume Delta z = 1
       dcut%bdyArea   = e     ! if we assume Delta z = 1
       if (abs(f(3)) > 1.E-10) call bl_error("cutcells.f90: bdyNormal(3) not equal 0.")
       dcut%bdyNormal = f
       dcut%bdyCentroid = g
       dcut%splitIndex = h

       ! do centroid check!!!!!
       do j=1,2
          do i=1,2
             dcut%FaceExist(i,j) = arr_a(i,j)
             dcut%FaceArea(i,j)  = arr_c(i,j)
          end do
          dcut%FaceExist(3,j) = .FALSE.
          dcut%FaceArea(3,j)  = 0.0d0
       end do
       do k=1,3
          do j=1,2
             do i=1,2
                dcut%FaceCentroid(i,j,k) = arr_b(i,j,k)
             end do
             dcut%FaceCentroid(3,j,k) = 0.0d0
          end do
       end do    

       ! centroid sanity check:
       my_z_cent = dcut%centroid(3)
       if (abs(dcut%bdyCentroid(3)-my_z_cent) > 1.E-10) then
          call bl_error("cutcells.f90: bdyCentroid and centroid don't have same z-comp")
       end if
       if (abs(arr_b(3,1,1)-arr_b(3,2,1)) > 1.E-10 .OR. &
            abs(arr_b(3,1,2)-arr_b(3,2,2)) > 1.E-10) then
          call bl_error("cutcells.f90: centroids from upper and lower z-face don't match")
       end if
       if (abs(arr_b(1,1,3)-my_z_cent) > 1.E-10 .AND. dcut%FaceExist(1,1)) then
          call bl_error("cutcells.f90: facecentroid and centroid don't have same z-comp")
       end if


    else
       print *,'Problem in cutcell_build'
       print *,'Need dim =3 or =2 but dim = ', aab
       call bl_error("Stop")
    end if
 
  end subroutine cutcell_build

  subroutine cutcell_print(p)
    type(cutcell), intent(in) :: p
    print*, ' >>>>>>>>>> Cut cell ', p%cutNumber, ' <<<<<<<<<< ' 
    print*, ' cutIndex  = ', p%cutIndex
    print*, ' centroid  = ', p%centroid
    print*, ' volume    = ', p%volume
    print*, ' bdyArea   = ', p%bdyArea
    print*, ' bdyNormal = ', p%bdyNormal
    print*, ' bdyCentroid = ', p%bdyCentroid
    print*, ' splitIndex = ', p%splitIndex
    print*, ' x face:'
    print*, ' face exist = ',p%FaceExist(1,1), p%FaceExist(1,2)
    write(*,'(A,es12.2,es12.2,es12.2)') '  face centroid = ', p%FaceCentroid(1,1,1), &
                             p%FaceCentroid(1,1,2), p%FaceCentroid(1,1,3)
    write(*,'(A,es12.2,es12.2,es12.2)') '  face centroid = ',p%FaceCentroid(1,2,1), &
                             p%FaceCentroid(1,2,2),p%FaceCentroid(1,2,3)
    write(*,'(A,es12.2,es12.2)') '  face area = ',p%FaceArea(1,1), p%FaceArea(1,2)
    print*, ' y face:'
    print*, ' face exist = ',p%FaceExist(2,1), p%FaceExist(2,2)
    write(*,'(A,es12.2,es12.2,es12.2)') '  face centroid = ', p%FaceCentroid(2,1,1), &
                             p%FaceCentroid(2,1,2), p%FaceCentroid(2,1,3)
    write(*,'(A,es12.2,es12.2,es12.2)') '  face centroid = ',p%FaceCentroid(2,2,1), &
                             p%FaceCentroid(2,2,2),p%FaceCentroid(2,2,3)
    write(*,'(A,es12.2,es12.2)') '  face area = ',p%FaceArea(2,1), p%FaceArea(2,2)
    print*, ' z face:'
    print*, ' face exist = ',p%FaceExist(3,1), p%FaceExist(3,2)

    write(*,'(A,es12.2,es12.2,es12.2)') '  face centroid = ', p%FaceCentroid(3,1,1), &
                             p%FaceCentroid(3,1,2), p%FaceCentroid(3,1,3)
    write(*,'(A,es12.2,es12.2,es12.2)') '  face centroid = ',p%FaceCentroid(3,2,1), &
                             p%FaceCentroid(3,2,2),p%FaceCentroid(3,2,3)
    write(*,'(A,es12.2,es12.2)') '  face area = ',p%FaceArea(3,1), p%FaceArea(3,2)
    print*,''
  end subroutine cutcell_print

  subroutine cutcell_container_build(ccc,sz)
    type(cutcell_container), intent(inout) :: ccc
    integer,                 intent(in   ) :: sz
    if ( .not. associated(ccc%d) ) then
       ccc%size = sz
       allocate(ccc%d(ccc%size))
    end if
  end subroutine cutcell_container_build

  subroutine cutcell_container_destroy(ccc)
    type(cutcell_container), intent(inout) :: ccc
    if ( associated(ccc%d) ) then
       deallocate(ccc%d)
       ccc%size  = 0
       ccc%d     => Null()
    end if
  end subroutine cutcell_container_destroy

  pure function cutcell_container_empty(ccc) result(r)
    logical :: r
    type(cutcell_container), intent(in) :: ccc
    r = (ccc%size == 0)
  end function cutcell_container_empty

  pure function cutcell_container_size(ccc) result(r)
    integer :: r
    type(cutcell_container), intent(in) :: ccc
    r = ccc%size
  end function cutcell_container_size

  subroutine cutcell_container_print(ccc, str)
    type(cutcell_container), intent(in)           :: ccc
    character (len=*),        intent(in), optional :: str

    integer :: i

    if ( present(str) ) print*, str
    if ( empty(ccc) ) then
       print*, '"Empty"'
    else
       do i = 1, ccc%size
          call print(ccc%d(i))
       end do
    end if
  end subroutine cutcell_container_print

  function cutcell_container_dataptr(ccc) result(r)
    type(cutcell), pointer :: r(:)
    type(cutcell_container), intent(in) :: ccc
    r => ccc%d
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
