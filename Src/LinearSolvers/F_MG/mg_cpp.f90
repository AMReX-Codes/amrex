module cpp_mg_module
  use mg_module

  implicit none

  type mg_server
     logical        :: final = .false.
     integer        :: dim  = 0
     integer        :: han  = -1
     type(layout)   :: lay
     type(mg_tower) :: mgt
     type(box)      :: pd
     integer, pointer :: bc(:,:) => Null()
     type(multifab) :: rh
     type(multifab) :: uu
  end type mg_server

  integer, parameter :: max_mgt = 10
  type(mg_server), save :: mgts(max_mgt)

contains
  
  subroutine mgt_verify(mgt, str)
    integer, intent(in) :: mgt
    character(len=*), intent(in) :: str
    if ( mgt < 1 .or. mgt > max_mgt ) then
       call bl_error( trim(str) // ": MGT out of bounds")
    end if
    if ( mgts(mgt)%dim == 0 ) then
       call bl_error( trim(str) // ": MGT invalid DIM: not allocated")
    end if
    
  end subroutine mgt_verify

  subroutine mgt_verify_n(mgt, str, n, lo, hi)
    integer, intent(in) :: mgt, n, lo(:), hi(:)
    character(len=*), intent(in) :: str
    type(box) :: bx
    call mgt_verify(mgt, str)
    if ( n < 1 .or. n > nboxes(mgts(mgt)%lay) ) then
       call bl_error( trim(str) // ": Box out of bounds")
    end if
    bx = make_box(lo, hi)
    if ( bx /= get_box(mgts(mgt)%lay, n)) then
       call bl_error( trim(str) // ": Box no filling")
    end if
  end subroutine mgt_verify_n

  subroutine mgt_not_final(mgt, str)
    integer, intent(in) :: mgt
    character(len=*), intent(in) :: str
    call mgt_verify(mgt, str)
    if ( mgts(mgt)%final ) then
       call bl_error( trim(str) // ": Changes made to finalized solver!")
    end if
  end subroutine mgt_not_final

end module cpp_mg_module

subroutine mgt_alloc(mgt, dm, nlevel)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: dm, nlevel
  integer, intent(out) :: mgt
  integer i
  if ( nlevel /= 1 ) then
     call bl_error("MGT_ALLOC: you've caught me at a bad time")
  end if
  do i = 1, max_mgt
     if ( mgts(i)%dim == 0 ) then
        mgt = i
        mgts(i)%dim = dm
        mgts(i)%han =  i
     end if
  end do
end subroutine mgt_alloc

subroutine mgt_set_level(mgt, lev, nb, dm, lo, hi, pd_lo, pd_hi, bc)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: dm, mgt, lev, nb, lo(nb,dm), hi(nb,dm), pd_lo(dm), pd_hi(dm), bc(2,dm)
  type(box) :: bxs(nb)
  type(boxarray) :: ba
  integer   :: i
  integer   :: nc

  call mgt_verify(mgt, "MGT_SET_LEVEL")

  if ( dm /= mgts(mgt)%dim ) then
     call bl_error("MGT_SET_LEVEL: Input DIM doesn't match internal DIM")
  end if
  call build(mgts(mgt)%pd, pd_lo(1:dm), pd_hi(1:dm))
  do i = 1, nb
     bxs(i) = make_box(lo(i,:), hi(i,:))
  end do
  call build(ba, bxs)
  call build(mgts(mgt)%lay, ba)
  allocate(mgts(mgt)%bc(dm,2))
  do i = 1, dm
     mgts(mgt)%bc = transpose(bc)
  end do

  nc = 1
  call build(mgts(mgt)%uu, mgts(mgt)%lay, nc, 1)
  call build(mgts(mgt)%rh, mgts(mgt)%lay, nc, 0)
  call destroy(ba)
end subroutine mgt_set_level

subroutine mgt_finalize(mgt)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  call mgt_verify(mgt, "MGT_FINALIZE")
  call mg_tower_build(mgts(mgt)%mgt, mgts(mgt)%lay, mgts(mgt)%pd, mgts(mgt)%bc)
  mgts(mgt)%final = .true.

end subroutine mgt_finalize

subroutine mgt_set_rh_2d(mgt, lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: rh(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_SET_RH", n, lo, hi)

  rp => dataptr(mgts(mgt)%rh, n)
  rp(lo(1):hi(1), lo(2):hi(2),1,1) = rh(lo(1):hi(1), lo(2):hi(2))

end subroutine mgt_set_rh_2d

subroutine mgt_get_rh_2d(mgt, lev, n, rh, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(out) :: rh(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: rp(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_GET_RH", n, lo, hi)

  rp => dataptr(mgts(mgt)%rh, n)
  rh(lo(1):hi(1), lo(2):hi(2)) =   rp(lo(1):hi(1), lo(2):hi(2),1,1)

end subroutine mgt_get_rh_2d

subroutine mgt_set_uu_2d(mgt, lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(in) :: uu(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_SET_UU", n, lo, hi)

  up => dataptr(mgts(mgt)%uu, n)
  up(lo(1):hi(1), lo(2):hi(2),1,1) = uu(lo(1):hi(1), lo(2):hi(2))

end subroutine mgt_set_uu_2d

subroutine mgt_get_uu_2d(mgt, lev, n, uu, plo, phi, lo, hi)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, lev, n, lo(2), hi(2), plo(2), phi(2)
  real(kind=dp_t), intent(out) :: uu(plo(1):phi(1), plo(2):phi(2))
  real(kind=dp_t), pointer :: up(:,:,:,:)
  
  call mgt_verify_n(mgt, "MGT_GET_UU", n, lo, hi)

  up => dataptr(mgts(mgt)%uu, n)
  uu(lo(1):hi(1), lo(2):hi(2)) = up(lo(1):hi(1), lo(2):hi(2),1,1)

end subroutine mgt_get_uu_2d

subroutine mgt_dealloc(mgt)
  use cpp_mg_module
  implicit none
  integer, intent(inout) :: mgt
  
  call mgt_verify(mgt, "MGT_DEALLOC")

  call mg_tower_destroy(mgts(mgt)%mgt)
  call layout_destroy(mgts(mgt)%lay)
  call destroy(mgts(mgt)%rh)
  call destroy(mgts(mgt)%uu)
  mgts(mgt)%han = -1
  mgts(mgt)%dim = 0
  mgt = -1
end subroutine mgt_dealloc

subroutine mgt_solve(mgt)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt

  call mgt_verify(mgt, "MGT_SOLVE")
  if ( .not. mgts(mgt)%final ) then
     call bl_error("MGT_SOLVE: MGT not finalized")
  end if

  call mg_tower_solve(mgts(mgt)%mgt, mgts(mgt)%uu, mgts(mgt)%rh)

end subroutine mgt_solve

subroutine mgt_set_nu1(mgt, nu1)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, nu1
  call mgt_not_final(mgt, "MGT_SET_NU1")
  mgts(mgt)%mgt%nu1 = nu1
end subroutine mgt_set_nu1

subroutine mgt_set_nu2(mgt, nu2)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt, nu2
  call mgt_not_final(mgt, "MGT_SET_NU1")
  mgts(mgt)%mgt%nu2 = nu2
end subroutine mgt_set_nu2

subroutine mgt_set_eps(mgt, eps)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  real(kind=dp_t), intent(in) :: eps
  call mgt_not_final(mgt, "MGT_SET_NU1")
  mgts(mgt)%mgt%eps = eps
end subroutine mgt_set_eps

subroutine mgt_set_bottom_solver(mgt, bottom_solver)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  integer, intent(in) :: bottom_solver
  call mgt_not_final(mgt, "MGT_SET_NU1") 
  mgts(mgt)%mgt%bottom_solver = bottom_solver
end subroutine mgt_set_bottom_solver

subroutine mgt_set_bottom_max_iter(mgt, bottom_max_iter)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  integer, intent(in) :: bottom_max_iter
  call mgt_not_final(mgt, "MGT_SET_NU1")
  mgts(mgt)%mgt%bottom_max_iter = bottom_max_iter
end subroutine mgt_set_bottom_max_iter

subroutine mgt_set_bottom_solver_eps(mgt, bottom_solver_eps)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  real(kind=dp_t), intent(in) :: bottom_solver_eps
  call mgt_not_final(mgt, "MGT_SET_NU1")
  mgts(mgt)%mgt%bottom_solver_eps = bottom_solver_eps
end subroutine mgt_set_bottom_solver_eps

subroutine mgt_set_max_nlevel(mgt, max_nlevel)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  integer, intent(in) :: max_nlevel
  call mgt_not_final(mgt, "MGT_SET_NU1")
  mgts(mgt)%mgt%max_nlevel = max_nlevel
end subroutine mgt_set_max_nlevel

subroutine mgt_set_min_width(mgt, min_width)
  use cpp_mg_module
  implicit none
  integer, intent(in) :: mgt
  integer, intent(in) :: min_width
  call mgt_not_final(mgt, "MGT_SET_NU1")
  mgts(mgt)%mgt%min_width = min_width
end subroutine mgt_set_min_width
