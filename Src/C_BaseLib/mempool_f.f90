module mempool_module

  use iso_c_binding

  implicit none

  integer (kind=c_size_t), parameter, private :: szd = 8_c_size_t
  integer (kind=c_size_t), parameter, private :: szi = 4_c_size_t

  interface bl_allocate
     module procedure bl_allocate_d1
     module procedure bl_allocate_d2
     module procedure bl_allocate_d3
     module procedure bl_allocate_d4
     module procedure bl_allocate_d5
     module procedure bl_allocate_d6
     module procedure bl_allocate_i1
     module procedure bl_allocate_i2
     module procedure bl_allocate_i3
     module procedure bl_allocate_d1_v
     module procedure bl_allocate_d2_v
     module procedure bl_allocate_d3_v
     module procedure bl_allocate_d1_vc
     module procedure bl_allocate_d2_vc
     module procedure bl_allocate_d3_vc
  end interface

  interface bl_deallocate
     module procedure bl_deallocate_d1
     module procedure bl_deallocate_d2
     module procedure bl_deallocate_d3
     module procedure bl_deallocate_d4
     module procedure bl_deallocate_d5
     module procedure bl_deallocate_d6
     module procedure bl_deallocate_i1
     module procedure bl_deallocate_i2
     module procedure bl_deallocate_i3
  end interface

  interface 
     function mempool_alloc (nbytes) result(p) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr) :: p
       integer(kind=c_size_t), intent(in), value :: nbytes
     end function mempool_alloc
     
     subroutine mempool_free (p) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: p
     end subroutine mempool_free

     subroutine double_array_init (p, n) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: p
       integer(kind=c_size_t), intent(in), value :: n
     end subroutine double_array_init
  end interface

contains

  subroutine bl_allocate_d1(a, lo1, hi1)
    double precision, pointer, intent(inout) :: a(:)
    integer, intent(in) :: lo1, hi1
    integer :: n1
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    double precision, pointer :: fp(:)
    n1 = max(hi1-lo1+1, 1)
    sz = int(n1,c_size_t)
    cp = mempool_alloc(szd*sz)
    call double_array_init(cp, sz)
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
    n1 = max(hi1-lo1+1, 1)
    n2 = max(hi2-lo2+1, 1)
    sz = int(n1,c_size_t) * int(n2,c_size_t) 
    cp = mempool_alloc(szd*sz)
    call double_array_init(cp, sz)
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
    n1 = max(hi1-lo1+1, 1)
    n2 = max(hi2-lo2+1, 1)
    n3 = max(hi3-lo3+1, 1)
    sz = int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t)
    cp = mempool_alloc(szd*sz)
    call double_array_init(cp, sz)
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
    n1 = max(hi1-lo1+1, 1)
    n2 = max(hi2-lo2+1, 1)
    n3 = max(hi3-lo3+1, 1)
    n4 = max(hi4-lo4+1, 1)
    sz = int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t) &
         * int(n4,c_size_t)
    cp = mempool_alloc(szd*sz)
    call double_array_init(cp, sz)
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
    n1 = max(hi1-lo1+1, 1)
    n2 = max(hi2-lo2+1, 1)
    n3 = max(hi3-lo3+1, 1)
    n4 = max(hi4-lo4+1, 1)
    n5 = max(hi5-lo5+1, 1)
    sz = int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t) &
         * int(n4,c_size_t) * int(n5,c_size_t)
    cp = mempool_alloc(szd*sz)
    call double_array_init(cp, sz)
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
    n1 = max(hi1-lo1+1, 1)
    n2 = max(hi2-lo2+1, 1)
    n3 = max(hi3-lo3+1, 1)
    n4 = max(hi4-lo4+1, 1)
    n5 = max(hi5-lo5+1, 1)
    n6 = max(hi6-lo6+1, 1)
    sz = int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t) &
         * int(n4,c_size_t) * int(n5,c_size_t) * int(n6,c_size_t)
    cp = mempool_alloc(szd*sz)
    call double_array_init(cp, sz)
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

  subroutine bl_allocate_d1_v(a, lo, hi)
    double precision, pointer, intent(inout) :: a(:)
    integer, intent(in) :: lo(1), hi(1)
    integer :: n(1)
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    double precision, pointer :: fp(:)
    n = hi - lo + 1
    sz = int(n(1),c_size_t)
    cp = mempool_alloc(szd*sz)
    call double_array_init(cp, sz)
    call c_f_pointer(cp, fp, shape=n)
    call shift_bound_d1_v(fp, lo, a)
  contains
    subroutine shift_bound_d1_v(fp, lo, a)
      integer, intent(in) :: lo(1)
      double precision, target, intent(in) :: fp(lo(1):)
      double precision, pointer, intent(inout) :: a(:)
      a => fp
    end subroutine shift_bound_d1_v
  end subroutine bl_allocate_d1_v

  subroutine bl_allocate_d2_v(a, lo, hi)
    double precision, pointer, intent(inout) :: a(:,:)
    integer, intent(in) :: lo(2), hi(2)
    integer :: n(2)
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    double precision, pointer :: fp(:,:)
    n = hi - lo + 1
    sz = int(n(1),c_size_t) * int(n(2),c_size_t)
    cp = mempool_alloc(szd*sz)
    call double_array_init(cp, sz)
    call c_f_pointer(cp, fp, shape=n)
    call shift_bound_d2_v(fp, lo, a)
  contains
    subroutine shift_bound_d2_v(fp, lo, a)
      integer, intent(in) :: lo(2)
      double precision, target, intent(in) :: fp(lo(1):,lo(2):)
      double precision, pointer, intent(inout) :: a(:,:)
      a => fp
    end subroutine shift_bound_d2_v
  end subroutine bl_allocate_d2_v

  subroutine bl_allocate_d3_v(a, lo, hi)
    double precision, pointer, intent(inout) :: a(:,:,:)
    integer, intent(in) :: lo(3), hi(3)
    integer :: n(3)
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    double precision, pointer :: fp(:,:,:)
    n = hi - lo + 1
    sz = int(n(1),c_size_t) * int(n(2),c_size_t) * int(n(3),c_size_t)
    cp = mempool_alloc(szd*sz)
    call double_array_init(cp, sz)
    call c_f_pointer(cp, fp, shape=n)
    call shift_bound_d3_v(fp, lo, a)
  contains
    subroutine shift_bound_d3_v(fp, lo, a)
      integer, intent(in) :: lo(3)
      double precision, target, intent(in) :: fp(lo(1):,lo(2):,lo(3):)
      double precision, pointer, intent(inout) :: a(:,:,:)
      a => fp
    end subroutine shift_bound_d3_v
  end subroutine bl_allocate_d3_v

  subroutine bl_allocate_d1_vc(a, lo, hi, ncomp)
    double precision, pointer, intent(inout) :: a(:,:)
    integer, intent(in) :: lo(1), hi(1), ncomp
    integer :: n(2)
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    double precision, pointer :: fp(:,:)
    n(1:1) = hi - lo + 1
    n(2) = ncomp
    sz = int(n(1),c_size_t) * int(n(2),c_size_t)
    cp = mempool_alloc(szd*sz)
    call double_array_init(cp, sz)
    call c_f_pointer(cp, fp, shape=n)
    call shift_bound_d1_vc(fp, lo, a)
  contains
    subroutine shift_bound_d1_vc(fp, lo, a)
      integer, intent(in) :: lo(1)
      double precision, target, intent(in) :: fp(lo(1):,1:)
      double precision, pointer, intent(inout) :: a(:,:)
      a => fp
    end subroutine shift_bound_d1_vc
  end subroutine bl_allocate_d1_vc

  subroutine bl_allocate_d2_vc(a, lo, hi, ncomp)
    double precision, pointer, intent(inout) :: a(:,:,:)
    integer, intent(in) :: lo(2), hi(2), ncomp
    integer :: n(3)
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    double precision, pointer :: fp(:,:,:)
    n(1:2) = hi - lo + 1
    n(3) = ncomp
    sz = int(n(1),c_size_t) * int(n(2),c_size_t) * int(n(3),c_size_t)
    cp = mempool_alloc(szd*sz)
    call double_array_init(cp, sz)
    call c_f_pointer(cp, fp, shape=n)
    call shift_bound_d2_vc(fp, lo, a)
  contains
    subroutine shift_bound_d2_vc(fp, lo, a)
      integer, intent(in) :: lo(2)
      double precision, target, intent(in) :: fp(lo(1):,lo(2):,1:)
      double precision, pointer, intent(inout) :: a(:,:,:)
      a => fp
    end subroutine shift_bound_d2_vc
  end subroutine bl_allocate_d2_vc

  subroutine bl_allocate_d3_vc(a, lo, hi, ncomp)
    double precision, pointer, intent(inout) :: a(:,:,:,:)
    integer, intent(in) :: lo(3), hi(3), ncomp
    integer :: n(4)
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    double precision, pointer :: fp(:,:,:,:)
    n(1:3) = hi - lo + 1
    n(4) = ncomp
    sz = int(n(1),c_size_t) * int(n(2),c_size_t) * int(n(3),c_size_t) &
         * int(n(4),c_size_t)
    cp = mempool_alloc(szd*sz)
    call double_array_init(cp, sz)
    call c_f_pointer(cp, fp, shape=n)
    call shift_bound_d3_vc(fp, lo, a)
  contains
    subroutine shift_bound_d3_vc(fp, lo, a)
      integer, intent(in) :: lo(3)
      double precision, target, intent(in) :: fp(lo(1):,lo(2):,lo(3):,1:)
      double precision, pointer, intent(inout) :: a(:,:,:,:)
      a => fp
    end subroutine shift_bound_d3_vc
  end subroutine bl_allocate_d3_vc

  subroutine bl_deallocate_d1(a)
    double precision, pointer, intent(inout) :: a(:)
    integer :: lo(1)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1)))
    call mempool_free(cp)
    a => Null()
  end subroutine bl_deallocate_d1

  subroutine bl_deallocate_d2(a)
    double precision, pointer, intent(inout) :: a(:,:)
    integer :: lo(2)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1),lo(2)))
    call mempool_free(cp)
    a => Null()
  end subroutine bl_deallocate_d2

  subroutine bl_deallocate_d3(a)
    double precision, pointer, intent(inout) :: a(:,:,:)
    integer :: lo(3)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1),lo(2),lo(3)))
    call mempool_free(cp)
    a => Null()
  end subroutine bl_deallocate_d3

  subroutine bl_deallocate_d4(a)
    double precision, pointer, intent(inout) :: a(:,:,:,:)
    integer :: lo(4)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1),lo(2),lo(3),lo(4)))
    call mempool_free(cp)
    a => Null()
  end subroutine bl_deallocate_d4

  subroutine bl_deallocate_d5(a)
    double precision, pointer, intent(inout) :: a(:,:,:,:,:)
    integer :: lo(5)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1),lo(2),lo(3),lo(4),lo(5)))
    call mempool_free(cp)
    a => Null()
  end subroutine bl_deallocate_d5

  subroutine bl_deallocate_d6(a)
    double precision, pointer, intent(inout) :: a(:,:,:,:,:,:)
    integer :: lo(6)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1),lo(2),lo(3),lo(4),lo(5),lo(6)))
    call mempool_free(cp)
    a => Null()
  end subroutine bl_deallocate_d6

  subroutine bl_allocate_i1(a, lo1, hi1)
    integer, pointer, intent(inout) :: a(:)
    integer, intent(in) :: lo1, hi1
    integer :: n1
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    integer, pointer :: fp(:)
    n1 = max(hi1-lo1+1, 1)
    sz = szi * int(n1,c_size_t)
    cp = mempool_alloc(sz)
    call c_f_pointer(cp, fp, shape=(/n1/))
    call shift_bound_i1(fp, lo1, a)
!    a(lo1:) => fp  ! some compilers may not support this
  contains
    subroutine shift_bound_i1(fp, lo1, a)
      integer, intent(in) :: lo1
      integer, target, intent(in) :: fp(lo1:)
      integer, pointer, intent(inout) :: a(:)
      a => fp
    end subroutine shift_bound_i1
  end subroutine bl_allocate_i1

  subroutine bl_allocate_i2(a, lo1, hi1, lo2, hi2)
    integer, pointer, intent(inout) :: a(:,:)
    integer, intent(in) :: lo1, hi1, lo2, hi2
    integer :: n1, n2
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    integer, pointer :: fp(:,:)
    n1 = max(hi1-lo1+1, 1)
    n2 = max(hi2-lo2+1, 1)
    sz = szi * int(n1,c_size_t) * int(n2,c_size_t)
    cp = mempool_alloc(sz)
    call c_f_pointer(cp, fp, shape=(/n1,n2/))
    call shift_bound_i2(fp, lo1, lo2, a)
!    a(lo1:,lo2:) => fp
  contains
    subroutine shift_bound_i2(fp, lo1, lo2, a)
      integer, intent(in) :: lo1, lo2
      integer, target, intent(in) :: fp(lo1:,lo2:)
      integer, pointer, intent(inout) :: a(:,:)
      a => fp
    end subroutine shift_bound_i2
  end subroutine bl_allocate_i2

  subroutine bl_allocate_i3(a, lo1, hi1, lo2, hi2, lo3, hi3)
    integer, pointer, intent(inout) :: a(:,:,:)
    integer, intent(in) :: lo1, hi1, lo2, hi2, lo3, hi3
    integer :: n1, n2, n3
    integer (kind=c_size_t) :: sz
    type(c_ptr) :: cp
    integer, pointer :: fp(:,:,:)
    n1 = max(hi1-lo1+1, 1)
    n2 = max(hi2-lo2+1, 1)
    n3 = max(hi3-lo3+1, 1)
    sz = szi * int(n1,c_size_t) * int(n2,c_size_t) * int(n3,c_size_t)
    cp = mempool_alloc(sz)
    call c_f_pointer(cp, fp, shape=(/n1,n2,n3/))
    call shift_bound_i3(fp, lo1, lo2, lo3, a)
!    a(lo1:,lo2:,lo3:) => fp
  contains
    subroutine shift_bound_i3(fp, lo1, lo2, lo3, a)
      integer, intent(in) :: lo1, lo2, lo3
      integer, target, intent(in) :: fp(lo1:,lo2:,lo3:)
      integer, pointer, intent(inout) :: a(:,:,:)
      a => fp
    end subroutine shift_bound_i3
  end subroutine bl_allocate_i3

  subroutine bl_deallocate_i1(a)
    integer, pointer, intent(inout) :: a(:)
    integer :: lo(1)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1)))
    call mempool_free(cp)
    a => Null()
  end subroutine bl_deallocate_i1

  subroutine bl_deallocate_i2(a)
    integer, pointer, intent(inout) :: a(:,:)
    integer :: lo(2)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1),lo(2)))
    call mempool_free(cp)
    a => Null()
  end subroutine bl_deallocate_i2

  subroutine bl_deallocate_i3(a)
    integer, pointer, intent(inout) :: a(:,:,:)
    integer :: lo(3)
    type(c_ptr) :: cp
    lo = lbound(a)
    cp = c_loc(a(lo(1),lo(2),lo(3)))
    call mempool_free(cp)
    a => Null()
  end subroutine bl_deallocate_i3

end module mempool_module
