!! An interface to Beebe's IEEE functions
module bl_ieee_module

  implicit none

  interface bl_isden
     function disden(x) result(r)
       use bl_types
       logical :: r
       real(kind=dp_t), intent(in) :: x
     end function disden
     function sisden(x) result(r)
       use bl_types
       logical :: r
       real(kind=sp_t), intent(in) :: x
     end function sisden
  end interface

  interface bl_isinf
     function disinf(x) result(r)
       use bl_types
       logical :: r
       real(kind=dp_t), intent(in) :: x
     end function disinf
     function sisinf(x) result(r)
       use bl_types
       logical :: r
       real(kind=sp_t), intent(in) :: x
     end function sisinf
  end interface

  interface bl_isnan
     function disnan(x) result(r)
       use bl_types
       logical :: r
       real(kind=dp_t), intent(in) :: x
     end function disnan
     function sisnan(x) result(r)
       use bl_types
       logical :: r
       real(kind=sp_t), intent(in) :: x
     end function sisnan
  end interface

end module bl_ieee_module
