module mempool_module

  use iso_c_binding

  implicit none

  integer, parameter, private :: szd = 8_c_size_t

  interface bl_allocate
     module procedure bl_allocate_d1
     module procedure bl_allocate_d2
     module procedure bl_allocate_d3
     module procedure bl_allocate_d4
     module procedure bl_allocate_d5
     module procedure bl_allocate_d6
  end interface

  interface bl_deallocate
     module procedure bl_deallocate_d1
     module procedure bl_deallocate_d2
     module procedure bl_deallocate_d3
     module procedure bl_deallocate_d4
     module procedure bl_deallocate_d5
     module procedure bl_deallocate_d6
  end interface

  interface 
     function bl_allocate_c (nbytes) result(p) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr) :: p
       integer(kind=c_size_t), intent(in), value :: nbytes
     end function bl_allocate_c
     
     subroutine bl_deallocate_c (p) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: p
     end subroutine bl_deallocate_c
  end interface

contains

  subroutine bl_allocate_d1(a, lo1, hi1)
    double precision, pointer, intent(inout) :: a(:)
    integer, intent(in) :: lo1, hi1
    integer :: n1
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    double precision, pointer :: fp(:)
    n1 = hi1-lo1+1 
    sz = szd * int(n1,c_size_t)
    cp = bl_allocate_c(sz)
    call c_f_pointer(cp, fp, shape=(/n1/))
    call shift_bound_d1(fp, lo1, a)
!    a(lo1:) => fp  ! some compilers may not support this
  contains
    subroutine shift_bound_d1(fp, lo1, a)
      integer, intent(in) :: lo1
      double precision, target, intent(in) :: fp(lo1:)
      double precision, pointer, intent(inout) :: a(:)
      a => fp
    end subroutine shift_bound_d1
  end subroutine bl_allocate_d1

  subroutine bl_allocate_d2(a, lo1, hi1, lo2, hi2)
    double precision, pointer, intent(inout) :: a(:,:)
    integer, intent(in) :: lo1, hi1, lo2, hi2
    integer :: n1, n2
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    double precision, pointer :: fp(:,:)
    n1 = hi1-lo1+1 
    n2 = hi2-lo2+1 
    sz = szd * int(n1,c_size_t) * int(n2,c_size_t)
    cp = bl_allocate_c(sz)
    call c_f_pointer(cp, fp, shape=(/n1,n2/))
    call shift_bound_d2(fp, lo1, lo2, a)
!    a(lo1:,lo2:) => fp
  contains
    subroutine shift_bound_d2(fp, lo1, lo2, a)
      integer, intent(in) :: lo1, lo2
      double precision, target, intent(in) :: fp(lo1:,lo2:)
      double precision, pointer, intent(inout) :: a(:,:)
      a => fp
    end subroutine shift_bound_d2
  end subroutine bl_allocate_d2

  subroutine bl_allocate_d3(a, lo1, hi1, lo2, hi2, lo3, hi3)
    double precision, pointer, intent(inout) :: a(:,:,:)
    integer, intent(in) :: lo1, hi1, lo2, hi2, lo3, hi3
    integer :: n1, n2, n3
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    double precision, pointer :: fp(:,:,:)
    n1 = hi1-lo1+1 
    n2 = hi2-lo2+1 
    n3 = hi3-lo3+1 
    sz = szd * int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t)
    cp = bl_allocate_c(sz)
    call c_f_pointer(cp, fp, shape=(/n1,n2,n3/))
    call shift_bound_d3(fp, lo1, lo2, lo3, a)
!    a(lo1:,lo2:,lo3:) => fp
  contains
    subroutine shift_bound_d3(fp, lo1, lo2, lo3, a)
      integer, intent(in) :: lo1, lo2, lo3
      double precision, target, intent(in) :: fp(lo1:,lo2:,lo3:)
      double precision, pointer, intent(inout) :: a(:,:,:)
      a => fp
    end subroutine shift_bound_d3
  end subroutine bl_allocate_d3

  subroutine bl_allocate_d4(a, lo1, hi1, lo2, hi2, lo3, hi3, lo4, hi4)
    double precision, pointer, intent(inout) :: a(:,:,:,:)
    integer, intent(in) :: lo1, hi1, lo2, hi2, lo3, hi3, lo4, hi4
    integer :: n1, n2, n3, n4
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    double precision, pointer :: fp(:,:,:,:)
    n1 = hi1-lo1+1 
    n2 = hi2-lo2+1 
    n3 = hi3-lo3+1 
    n4 = hi4-lo4+1 
    sz = szd * int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t) &
         * int(n4,c_size_t)
    cp = bl_allocate_c(sz)
    call c_f_pointer(cp, fp, shape=(/n1,n2,n3,n4/))
    call shift_bound_d4(fp, lo1, lo2, lo3, lo4, a)
!    a(lo1:,lo2:,lo3:,lo4:) => fp
  contains
    subroutine shift_bound_d4(fp, lo1, lo2, lo3, lo4, a)
      integer, intent(in) :: lo1, lo2, lo3, lo4
      double precision, target, intent(in) :: fp(lo1:,lo2:,lo3:,lo4:)
      double precision, pointer, intent(inout) :: a(:,:,:,:)
      a => fp
    end subroutine shift_bound_d4
  end subroutine bl_allocate_d4

  subroutine bl_allocate_d5(a,lo1,hi1,lo2,hi2,lo3,hi3,lo4,hi4,lo5,hi5)
    double precision, pointer, intent(inout) :: a(:,:,:,:,:)
    integer, intent(in) :: lo1,hi1,lo2,hi2,lo3,hi3,lo4,hi4,lo5,hi5
    integer :: n1, n2, n3, n4, n5
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    double precision, pointer :: fp(:,:,:,:,:)
    n1 = hi1-lo1+1 
    n2 = hi2-lo2+1 
    n3 = hi3-lo3+1 
    n4 = hi4-lo4+1 
    n5 = hi5-lo5+1 
    sz = szd * int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t) &
         * int(n4,c_size_t) * int(n5,c_size_t)
    cp = bl_allocate_c(sz)
    call c_f_pointer(cp, fp, shape=(/n1,n2,n3,n4,n5/))
    call shift_bound_d5(fp, lo1, lo2, lo3, lo4, lo5, a)
!    a(lo1:,lo2:,lo3:,lo4:,lo5:) => fp
  contains
    subroutine shift_bound_d5(fp, lo1, lo2, lo3, lo4, lo5, a)
      integer, intent(in) :: lo1, lo2, lo3, lo4, lo5
      double precision, target, intent(in) :: fp(lo1:,lo2:,lo3:,lo4:,lo5:)
      double precision, pointer, intent(inout) :: a(:,:,:,:,:)
      a => fp
    end subroutine shift_bound_d5
  end subroutine bl_allocate_d5

  subroutine bl_allocate_d6(a,lo1,hi1,lo2,hi2,lo3,hi3,lo4,hi4,lo5,hi5,lo6,hi6)
    double precision, pointer, intent(inout) :: a(:,:,:,:,:,:)
    integer, intent(in) :: lo1,hi1,lo2,hi2,lo3,hi3,lo4,hi4,lo5,hi5,lo6,hi6
    integer :: n1, n2, n3, n4, n5, n6
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    double precision, pointer :: fp(:,:,:,:,:,:)
    n1 = hi1-lo1+1 
    n2 = hi2-lo2+1 
    n3 = hi3-lo3+1 
    n4 = hi4-lo4+1 
    n5 = hi5-lo5+1 
    n6 = hi6-lo6+1 
    sz = szd * int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t) &
         * int(n4,c_size_t) * int(n5,c_size_t) * int(n6,c_size_t)
    cp = bl_allocate_c(sz)
    call c_f_pointer(cp, fp, shape=(/n1,n2,n3,n4,n5,n6/))
    call shift_bound_d6(fp, lo1, lo2, lo3, lo4, lo5, lo6, a)
!    a(lo1:,lo2:,lo3:,lo4:,lo5:,lo6) => fp
  contains
    subroutine shift_bound_d6(fp, lo1, lo2, lo3, lo4, lo5, lo6, a)
      integer, intent(in) :: lo1, lo2, lo3, lo4, lo5, lo6
      double precision, target, intent(in) :: fp(lo1:,lo2:,lo3:,lo4:,lo5:,lo6:)
      double precision, pointer, intent(inout) :: a(:,:,:,:,:,:)
      a => fp
    end subroutine shift_bound_d6
  end subroutine bl_allocate_d6

  subroutine bl_deallocate_d1(a)
    double precision, pointer, intent(inout) :: a(:)
    integer :: lo(1)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1)))
    call bl_deallocate_c(cp)
  end subroutine bl_deallocate_d1

  subroutine bl_deallocate_d2(a)
    double precision, pointer, intent(inout) :: a(:,:)
    integer :: lo(2)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1),lo(2)))
    call bl_deallocate_c(cp)
  end subroutine bl_deallocate_d2

  subroutine bl_deallocate_d3(a)
    double precision, pointer, intent(inout) :: a(:,:,:)
    integer :: lo(3)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1),lo(2),lo(3)))
    call bl_deallocate_c(cp)
  end subroutine bl_deallocate_d3

  subroutine bl_deallocate_d4(a)
    double precision, pointer, intent(inout) :: a(:,:,:,:)
    integer :: lo(4)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1),lo(2),lo(3),lo(4)))
    call bl_deallocate_c(cp)
  end subroutine bl_deallocate_d4

  subroutine bl_deallocate_d5(a)
    double precision, pointer, intent(inout) :: a(:,:,:,:,:)
    integer :: lo(5)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1),lo(2),lo(3),lo(4),lo(5)))
    call bl_deallocate_c(cp)
  end subroutine bl_deallocate_d5

  subroutine bl_deallocate_d6(a)
    double precision, pointer, intent(inout) :: a(:,:,:,:,:,:)
    integer :: lo(6)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1),lo(2),lo(3),lo(4),lo(5),lo(6)))
    call bl_deallocate_c(cp)
  end subroutine bl_deallocate_d6

end module mempool_module
