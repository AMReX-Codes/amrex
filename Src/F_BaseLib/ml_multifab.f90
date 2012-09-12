module ml_multifab_module

  use ml_layout_module

  implicit none

  type ml_multifab
     integer                 :: dim    =  0
     integer                 :: nlevel =  0
     integer                 :: nc     =  1
     integer                 :: ng     =  0
     type(ml_layout)         :: mla
     type(multifab), pointer :: mf(:)  => Null()
  end type ml_multifab

  interface built_q
     module procedure ml_multifab_built_q
  end interface

  interface build
     module procedure ml_multifab_build
  end interface

  interface destroy
     module procedure ml_multifab_destroy
  end interface

  interface copy
     module procedure ml_multifab_copy
     module procedure ml_multifab_copy_c
  end interface

  interface setval
     module procedure ml_multifab_setval
  end interface

  interface rescale
     module procedure ml_multifab_rescale_c
     module procedure ml_multifab_rescale
  end interface

  interface dot
     module procedure ml_multifab_dot_cc
  end interface

  interface norm_l2
     module procedure ml_multifab_norm_l2
     module procedure ml_multifab_norm_l2_c
  end interface

  interface norm_inf
     module procedure ml_multifab_norm_inf
     module procedure ml_multifab_norm_inf_c
  end interface

  interface ncomp
     module procedure ml_multifab_ncomp
  end interface

  interface saxpy
     module procedure ml_multifab_saxpy_3_c
     module procedure ml_multifab_saxpy_3_cc
     module procedure ml_multifab_saxpy_3
     module procedure ml_multifab_saxpy_5
  end interface

  interface div
     module procedure ml_multifab_div_s_c
  end interface

  interface div_div
     module procedure ml_multifab_div_div
     module procedure ml_multifab_div_div_s
     module procedure ml_multifab_div_div_c_s
  end interface

  interface mult
     module procedure ml_multifab_mult_s_c
  end interface

  interface mult_mult
     module procedure ml_multifab_mult_mult
     module procedure ml_multifab_mult_mult_s
     module procedure ml_multifab_mult_mult_c_s
  end interface

  interface sub_sub
     module procedure ml_multifab_sub_sub_s
     module procedure ml_multifab_sub_sub
  end interface

  interface plus_plus
     module procedure ml_multifab_plus_plus_s
     module procedure ml_multifab_plus_plus
  end interface

  interface dataptr
     module procedure ml_multifab_dataptr
     module procedure ml_multifab_dataptr_c
     module procedure ml_multifab_dataptr_bx
     module procedure ml_multifab_dataptr_bx_c
  end interface

  interface get_box
     module procedure ml_multifab_get_box
  end interface

  interface get_ibox
     module procedure ml_multifab_get_ibox
  end interface

  interface get_pbox
     module procedure ml_multifab_get_pbox
  end interface

  interface get_boxarray
     module procedure ml_multifab_get_boxarray
  end interface

  interface nlevels
     module procedure ml_multifab_nlevels
  end interface

  interface print
     module procedure ml_multifab_print
  end interface

contains

  function ml_multifab_nlevels(mmf) result(r)
    integer :: r
    type(ml_multifab), intent(in) :: mmf
    r = mmf%nlevel
  end function ml_multifab_nlevels

  function ml_multifab_get_box(mmf, lev, n) result(r)
    type(box) :: r
    type(ml_multifab), intent(in) :: mmf
    integer, intent(in) :: lev, n
    r = get_box(mmf%mf(lev), n)
  end function ml_multifab_get_box

  function ml_multifab_get_ibox(mmf, lev, n) result(r)
    type(box) :: r
    type(ml_multifab), intent(in) :: mmf
    integer, intent(in) :: lev, n
    r = get_ibox(mmf%mf(lev), n)
  end function ml_multifab_get_ibox

  function ml_multifab_get_pbox(mmf, lev, n) result(r)
    type(box) :: r
    type(ml_multifab), intent(in) :: mmf
    integer, intent(in) :: lev, n
    r = get_pbox(mmf%mf(lev), n)
  end function ml_multifab_get_pbox

  function ml_multifab_get_boxarray(mmf, lev) result(r)
    type(boxarray) :: r
    type(ml_multifab), intent(in) :: mmf
    integer, intent(in) :: lev
    r = get_boxarray(mmf%mf(lev))
  end function ml_multifab_get_boxarray

  function ml_multifab_built_q(mmf) result(r)
    logical :: r
    type(ml_multifab), intent(in) :: mmf
    r = mmf%dim /= 0
  end function ml_multifab_built_q

  subroutine ml_multifab_saxpy_3_cc(a, ia, b1, b, ib, nc)
    real(dp_t), intent(in) :: b1
    type(ml_multifab), intent(inout) :: a
    type(ml_multifab), intent(in)  :: b
    integer, intent(in) :: ia, ib
    integer, intent(in), optional :: nc
    integer :: n
    do n = 1, a%nlevel
       call multifab_saxpy_3_cc(a%mf(n), ia, b1, b%mf(n), ib, nc)
    end do
  end subroutine ml_multifab_saxpy_3_cc
  subroutine ml_multifab_saxpy_3_c(a, ia, b1, b)
    real(dp_t), intent(in) :: b1
    type(ml_multifab), intent(inout) :: a
    type(ml_multifab), intent(in)  :: b
    integer, intent(in) :: ia
    integer :: n
    do n = 1, a%nlevel
       call multifab_saxpy_3_c(a%mf(n), ia, b1, b%mf(n))
    end do
  end subroutine ml_multifab_saxpy_3_c
  subroutine ml_multifab_saxpy_3(a, b1, b)
    real(dp_t), intent(in) :: b1
    type(ml_multifab), intent(inout) :: a
    type(ml_multifab), intent(in)  :: b
    integer :: n
    do n = 1, a%nlevel
       call multifab_saxpy_3(a%mf(n), b1, b%mf(n))
    end do
  end subroutine ml_multifab_saxpy_3
  subroutine ml_multifab_saxpy_5(a, b1, b, c1, c)
    real(dp_t), intent(in) :: b1, c1
    type(ml_multifab), intent(inout) :: a
    type(ml_multifab), intent(in)  :: b, c
    integer :: n
    do n = 1, a%nlevel
       call multifab_saxpy_5(a%mf(n), b1, b%mf(n), c1, c%mf(n))
    end do
  end subroutine ml_multifab_saxpy_5

  subroutine ml_multifab_mult_mult(a, b)
    type(ml_multifab), intent(inout) :: a
    type(ml_multifab), intent(in)  :: b
    integer :: n
    do n = 1, a%nlevel
       call mult_mult(a%mf(n), b%mf(n))
    end do
  end subroutine ml_multifab_mult_mult

  subroutine ml_multifab_mult_mult_s(a, b)
    type(ml_multifab), intent(inout) :: a
    real(kind=dp_t), intent(in)  :: b
    integer :: n
    do n = 1, a%nlevel
       call mult_mult(a%mf(n), b)
    end do
  end subroutine ml_multifab_mult_mult_s

  subroutine ml_multifab_mult_mult_c_s(a, ia, b)
    type(ml_multifab), intent(inout) :: a
    integer, intent(in) :: ia
    real(kind=dp_t), intent(in)  :: b
    integer :: n
    do n = 1, a%nlevel
       call mult_mult(a%mf(n), ia, b, ncomp(a%mf(n)), nghost(a%mf(n)))
    end do
  end subroutine ml_multifab_mult_mult_c_s
  
  subroutine ml_multifab_mult_s_c(a, ia, b, ib, val)
    type(ml_multifab), intent(inout) :: a
    type(ml_multifab), intent(in) :: b
    real(kind=dp_t), intent(in)  :: val
    integer, intent(in) :: ia, ib
    integer :: n
    do n = 1, a%nlevel
       call mult(a%mf(n), ia, b%mf(n), ib, val, ncomp(a%mf(n)), nghost(a%mf(n)))
    end do
  end subroutine ml_multifab_mult_s_c

  subroutine ml_multifab_div_div(a, b)
    type(ml_multifab), intent(inout) :: a
    type(ml_multifab), intent(in)  :: b
    integer :: n
    do n = 1, a%nlevel
       call div_div(a%mf(n), b%mf(n))
    end do
  end subroutine ml_multifab_div_div

  subroutine ml_multifab_div_div_s(a, b)
    type(ml_multifab), intent(inout) :: a
    real(kind=dp_t), intent(in)  :: b
    integer :: n
    do n = 1, a%nlevel
       call div_div(a%mf(n), b)
    end do
  end subroutine ml_multifab_div_div_s

  subroutine ml_multifab_div_div_c_s(a, ia, b)
    type(ml_multifab), intent(inout) :: a
    integer, intent(in) :: ia
    real(kind=dp_t), intent(in)  :: b
    integer :: n
    do n = 1, a%nlevel
       call div_div(a%mf(n), ia, b, ncomp(a%mf(n)), nghost(a%mf(n)))
    end do
  end subroutine ml_multifab_div_div_c_s

  subroutine ml_multifab_div_s_c(a, ia, b, ib, val)
    type(ml_multifab), intent(inout) :: a
    type(ml_multifab), intent(in) :: b
    real(kind=dp_t), intent(in)  :: val
    integer, intent(in) :: ia, ib
    integer :: n
    do n = 1, a%nlevel
       call div(a%mf(n), ia, b%mf(n), ib, val, ncomp(a%mf(n)), nghost(a%mf(n)))
    end do
  end subroutine ml_multifab_div_s_c

  subroutine ml_multifab_sub_sub_s(a, b)
    type(ml_multifab), intent(inout) :: a
    real(kind=dp_t), intent(in)  :: b
    integer :: n
    do n = 1, a%nlevel
       call sub_sub(a%mf(n), b)
    end do
  end subroutine ml_multifab_sub_sub_s
  subroutine ml_multifab_sub_sub(a, b)
    type(ml_multifab), intent(inout) :: a
    type(ml_multifab), intent(in)  :: b
    integer :: n
    do n = 1, a%nlevel
       call sub_sub(a%mf(n), b%mf(n))
    end do
  end subroutine ml_multifab_sub_sub

  subroutine ml_multifab_plus_plus_s(a, b)
    type(ml_multifab), intent(inout) :: a
    real(kind=dp_t), intent(in)  :: b
    integer :: n
    do n = 1, a%nlevel
       call plus_plus(a%mf(n), b)
    end do
  end subroutine ml_multifab_plus_plus_s
  subroutine ml_multifab_plus_plus(a, b)
    type(ml_multifab), intent(inout) :: a
    type(ml_multifab), intent(in)  :: b
    integer :: n
    do n = 1, a%nlevel
       call plus_plus(a%mf(n), b%mf(n))
    end do
  end subroutine ml_multifab_plus_plus

  subroutine ml_multifab_copy(target, source, ng)
    type(ml_multifab), intent(inout) :: target
    type(ml_multifab), intent(in) :: source
    integer, intent(in), optional :: ng
    integer :: n, lng
    lng = 0;       if ( present(ng)     ) lng = ng
    if ( target%mla /= source%mla ) call bl_error("ML_MULTFAB_COPY: not same ml_layout")
    do n = 1, source%nlevel
       call copy(target%mf(n), source%mf(n), ng = lng)
    end do
  end subroutine ml_multifab_copy

  subroutine ml_multifab_copy_c(target, tc, source, sc, nc, ng)
    type(ml_multifab), intent(inout) :: target
    type(ml_multifab), intent(in) :: source
    integer, intent(in), optional :: ng
    integer, intent(in) :: tc, sc
    integer, intent(in), optional :: nc
    integer :: n, lng
    lng = 0;       if ( present(ng)     ) lng = ng
    if ( target%mla /= source%mla ) call bl_error("ML_MULTFAB_COPY: not same ml_layout")
    do n = 1, source%nlevel
       call copy(target%mf(n), tc, source%mf(n), sc, nc, ng = lng)
    end do
  end subroutine ml_multifab_copy_c

  subroutine ml_multifab_setval(target, val, all)
    type(ml_multifab), intent(inout) :: target
    real(kind=dp_t), intent(in) :: val
    logical, intent(in), optional :: all
    integer :: n
    do n = 1, target%nlevel
       call setval(target%mf(n), val, all)
    end do
  end subroutine ml_multifab_setval

  function ml_multifab_ncomp(mmf) result(r)
    integer :: r
    type(ml_multifab), intent(in) :: mmf
    r = mmf%nc
  end function ml_multifab_ncomp

  subroutine ml_multifab_build(mmf, mla, nc, ng, maxlev, nodal)
    type(ml_multifab), intent(inout) :: mmf
    type(ml_layout), intent(inout) :: mla
    integer, intent(in), optional :: nc
    integer, intent(in), optional :: ng
    integer, intent(in), optional :: maxlev
    logical, intent(in), optional :: nodal(:)
    integer :: n
    mmf%dim = mla%dim
    if ( present(maxlev) ) then
       if ( mla%nlevel < maxlev ) then
          call bl_error("ML_MULTIFAB_BUILD: mla%nlevel < maxlev: maxlev = ", maxlev)
       end if
       mmf%nlevel = maxlev
    else
       mmf%nlevel = mla%nlevel
    end if
    mmf%mla = mla
    if ( present(nc) ) mmf%nc = nc
    if ( present(ng) ) mmf%ng = ng
    allocate(mmf%mf(mmf%nlevel))
    do n = 1, mmf%nlevel
       call multifab_build(mmf%mf(n), mla%la(n), nc, ng, nodal)
    end do
  end subroutine ml_multifab_build

  subroutine ml_multifab_destroy(mmf)
    type(ml_multifab), intent(inout) :: mmf
    integer :: n
    do n = 1, mmf%nlevel
       call destroy(mmf%mf(n))
    end do
    deallocate(mmf%mf)
    mmf%dim = 0
    mmf%nlevel = 0
    mmf%nc = 1
    mmf%ng = 0
  end subroutine ml_multifab_destroy

  function ml_multifab_dot_cc(x, compx, y, compy) result(r)
    real(kind=dp_t) :: r
    type(ml_multifab), intent(in) :: x
    type(ml_multifab), intent(in) :: y
    integer, intent(in) :: compx, compy
    integer :: n
    type(bl_prof_timer), save :: bpt
    call build(bpt, "ml_multifab_dot_cc")

    if ( x%mla /= y%mla ) call bl_error("ML_DOT: incommensurate")
    n = x%nlevel
    r = dot(x%mf(n), compx, y%mf(n), compy)
    do n = x%nlevel-1, 1, -1
       r = r/product(x%mla%mba%rr(n,:)) &
            + dot(x%mf(n), compx, y%mf(n), compy, mask=x%mla%mask(n))
    end do

    call destroy(bpt)
  end function ml_multifab_dot_cc

  function ml_multifab_sum(x) result(r)
    real(kind=dp_t) :: r
    type(ml_multifab), intent(in) :: x
    integer :: n
    n = x%nlevel
    r = multifab_sum(x%mf(n))
    do n = x%nlevel-1, 1, -1
       r = r/product(x%mla%mba%rr(n,:)) &
            + multifab_sum(x%mf(n), mask = x%mla%mask(n))
    end do
  end function ml_multifab_sum

  function ml_multifab_norm_l2(x) result(r)
    real(kind=dp_t) :: r
    type(ml_multifab) :: x
    integer :: n
    n = x%nlevel
    r = norm_l2(x%mf(n))**2
    do n = x%nlevel-1, 1, -1
       r = r/product(x%mla%mba%rr(n,:)) &
            + norm_l2(x%mf(n), mask = x%mla%mask(n))**2
    end do
    r = sqrt(r)
  end function ml_multifab_norm_l2
  function ml_multifab_norm_l2_c(x, c) result(r)
    real(kind=dp_t) :: r
    type(ml_multifab) :: x
    integer, intent(in) :: c
    integer :: n
    n = x%nlevel
    r = norm_l2(x%mf(n), c)**2
    do n = x%nlevel-1, 1, -1
       r = r/product(x%mla%mba%rr(n,:)) &
            + norm_l2(x%mf(n), c, mask = x%mla%mask(n))**2
    end do
    r = sqrt(r)
  end function ml_multifab_norm_l2_c

  function ml_multifab_norm_inf(x) result(r)
    real(kind=dp_t) :: r
    type(ml_multifab) :: x
    integer :: n
    n = x%nlevel
    r = norm_inf(x%mf(n))
    do n = x%nlevel-1, 1, -1
       r = max(r, norm_inf(x%mf(n), mask = x%mla%mask(n)))
    end do
  end function ml_multifab_norm_inf
  function ml_multifab_norm_inf_c(x, c) result(r)
    real(kind=dp_t) :: r
    type(ml_multifab) :: x
    integer, intent(in) :: c
    integer :: n
    n = x%nlevel
    r = norm_inf(x%mf(n))
    do n = x%nlevel-1, 1, -1
       r = max(r, norm_inf(x%mf(n), c, mask = x%mla%mask(n)))
    end do
  end function ml_multifab_norm_inf_c

  subroutine ml_multifab_rescale_c(x, c, val, off)
    type(ml_multifab), intent(inout) :: x
    real(dp_t), intent(in) :: val
    integer, intent(in) :: c
    real(dp_t), intent(in), optional :: off
    integer :: n
    type(bl_prof_timer), save :: bpt
    call build(bpt, "ml_multifab_rescale_c")

    do n = 1, x%nlevel
       call rescale(x%mf(n), c, val, off)
    end do

    call destroy(bpt)
  end subroutine ml_multifab_rescale_c
  subroutine ml_multifab_rescale(x, val, off)
    real(dp_t), intent(in) :: val
    real(dp_t), intent(in), optional :: off
    type(ml_multifab), intent(inout) :: x
    integer :: n
    do n = 1, x%nlevel
       call rescale(x%mf(n), val, off)
    end do
  end subroutine ml_multifab_rescale

  function ml_multifab_dataptr(mmf, lev, n) result(r)
    real(dp_t), pointer :: r(:,:,:,:)
    type(ml_multifab), intent(in) :: mmf
    integer, intent(in) :: lev, n
    r => dataptr(mmf%mf(lev), n)
  end function ml_multifab_dataptr
  function ml_multifab_dataptr_c(mmf, lev, n, c, nc) result(r)
    real(dp_t), pointer :: r(:,:,:,:)
    type(ml_multifab), intent(in) :: mmf
    integer, intent(in) :: lev, n, c
    integer, intent(in), optional :: nc
    r => dataptr(mmf%mf(lev), n, c, nc)
  end function ml_multifab_dataptr_c
  function ml_multifab_dataptr_bx(mmf, lev, n, bx) result(r)
    real(dp_t), pointer :: r(:,:,:,:)
    type(ml_multifab), intent(in) :: mmf
    integer, intent(in) :: lev, n
    type(box), intent(in) :: bx
    r => dataptr(mmf%mf(lev), n, bx)
  end function ml_multifab_dataptr_bx
  function ml_multifab_dataptr_bx_c(mmf, lev, n, bx, c, nc) result(r)
    real(dp_t), pointer :: r(:,:,:,:)
    type(ml_multifab), intent(in) :: mmf
    integer, intent(in) :: lev, n, c
    integer, intent(in), optional :: nc
    type(box), intent(in) :: bx
    r => dataptr(mmf%mf(lev), n, bx, c, nc)
  end function ml_multifab_dataptr_bx_c

  subroutine ml_multifab_print(mmf, str, unit, all, data, skip)
    use bl_IO_module
    type(ml_multifab), intent(in) :: mmf
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: all, data
    integer, intent(in), optional :: skip
    integer :: n, un
    character(len=16) :: levstr
    un = unit_stdout(unit)
    if ( parallel_IOProcessor() ) then
       call unit_skip(un, skip)
       write(unit=un, fmt='("ML_MULTIFAB ", i1)', advance = 'NO') 
       if ( present(str) ) then
          write(unit=un, fmt='(": ",A)') str
       else
          write(unit=un, fmt='()')
       end if
    end if
    call unit_skip(un, skip)
    write(unit=un, fmt='(" DIM     = ",i2)') mmf%dim
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NC      = ",i2)') mmf%nc
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NG      = ",i2)') mmf%ng
    call unit_skip(un, skip)
    write(unit=un, fmt='(" NLEVEL  = ",i2)') mmf%nlevel
    do n = 1, mmf%nlevel
       write(unit=levstr,fmt='("LEVEL ", i1)') n
       call print(mmf%mf(n), &
            str = trim(levstr), &
            unit = unit,  &
            all = all, &
            data = data, &
            skip = unit_get_skip(skip) + 2 &
            )
    end do
  end subroutine ml_multifab_print

end module ml_multifab_module

