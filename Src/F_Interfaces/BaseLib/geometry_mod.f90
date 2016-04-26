
module geometry_module

  use iso_c_binding
  use bl_space_module, only : ndims => bl_num_dims
  use box_module

  implicit none

  private

  public :: geometry_build

  type, public :: Geometry
     logical          :: pmask(3)  = .false.
     double precision :: problo(3) = 0.0d0
     double precision :: probhi(3) = 1.0d0
     double precision :: dx(3)     = 0.0d0
     type(Box)        :: domain
     type(c_ptr)      :: p         = c_null_ptr
   contains
     final :: geometry_destroy
  end type Geometry

  ! interfaces to c++ functions

  interface
     subroutine fi_new_geometry (geom,lo,hi) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr) :: geom
       integer, intent(in) :: lo(3), hi(3)
     end subroutine fi_new_geometry

     subroutine fi_delete_geometry (geom) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: geom
     end subroutine fi_delete_geometry

     subroutine fi_geometry_get_pmask (geom,pmask) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: geom
       integer(c_int) :: pmask(3)
     end subroutine fi_geometry_get_pmask

     subroutine fi_geometry_get_probdomain (geom,problo,probhi) bind(c)
       use iso_c_binding
       implicit none
       type(c_ptr), value, intent(in) :: geom
       real(c_double) :: problo(3), probhi(3)
     end subroutine fi_geometry_get_probdomain
  end interface

contains

  subroutine geometry_build (geom, domain)
    type(Geometry) :: geom
    type(Box), intent(in) :: domain
    integer :: imask(3), i
    call fi_new_geometry(geom%p, domain%lo, domain%hi)
    imask = 0
    call fi_geometry_get_pmask(geom%p, imask)
    where (imask .eq. 1) geom%pmask = .true.
    call fi_geometry_get_probdomain(geom%p, geom%problo, geom%probhi)
    geom%domain = domain
    do i = 1, ndims
       geom%dx(i) = (geom%probhi(i)-geom%problo(i)) / dble(domain%hi(i)-domain%lo(i)+1)
    end do
  end subroutine geometry_build

  subroutine geometry_destroy (geom)
    type(Geometry) :: geom
    call fi_delete_geometry(geom%p)
  end subroutine geometry_destroy

end module geometry_module
