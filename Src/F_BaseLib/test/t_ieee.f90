
subroutine t_ieee
  use bl_types
  call t_rwinfn
  call t_minmax
end subroutine t_ieee

subroutine t_rwinfn
  use bl_ieee_module
  implicit none
  character(50) line
  real arg, s, x, nan, inf

  open (unit = 1, file = 'rwinfnan.tmp', &
       status = 'replace', &
       action = 'readwrite', &
       access = 'sequential', &
       form = 'formatted')

  s = 0.0
  x = s
  inf = 1.0/s
  nan = x/s

  print *, 'inf = ', inf
  print *, 'nan = ', nan

  if ( bl_isnan(nan) ) print *, 'nan is nan'
  if ( bl_isinf(inf) ) print *, 'inf is inf'

  write (1,*) nan
  write (1,*) inf

  rewind (1)
  read (1,'(a)') line
  write (6,*) 'NaN was written as: ', line
  read (1,'(a)') line
  write (6,*) 'Inf was written as: ', line

  rewind (1)
  read (1,*) x
  if (bl_isnan(x)) then
     write (6,*) 'NaN was correctly input'
  else
     write (6,*) 'NaN was INCORRECTLY input as ', x
  endif

  read (1,*) x
  if ( bl_isinf(x) )  then
     write (6,*) 'Inf was correctly input'
  else
     write (6,*) 'Inf was INCORRECTLY input as ', x
  endif

  close(1)
contains
  logical function isnan(arg)
    real arg
    isnan  = (arg .ne. arg)
  end function isnan
  logical function isninf(arg)
    real arg
    isninf = (arg .eq. -inf)
  end function isninf
  logical function ispinf(arg)
    real arg
    ispinf = (arg .eq. inf)
  end function ispinf

end subroutine t_rwinfn

subroutine t_minmax
  !     (Minimum/Maximum of NaN)
  !     Fortran program to test the handling of IEEE 754 NaN by the
  !     min() and max() intrinsic functions.
  implicit none
  real NaN

  NaN = store(0.0)
  NaN = NaN/store(NaN)

  call test(1.0, 2.0)
  call test(1.0, NaN)
  call test(NaN, 1.0)
  call test(NaN, NaN)

contains

  real function store(x)
    real x
    store = x
    return

  end function store

  subroutine test(x,y)
    real x, y

    write (6,1000) 'max', x, y, max(x,y)
    write (6,1000) 'min', x, y, min(x,y)
1000 format (a,'(', f15.2, ',', f15.2, ') ->', f15.2)
  end subroutine test

end subroutine t_minmax


subroutine t_nantest
  !====NAN.F90 illustrates what works and what doesn't when
  !    detecting a NaN
  ! Platforms: Windows 9x/Me/NT/2000, AIX 4.3.
  ! Compilers:I)Compaq Visual Fortran 6.6a with default
  !             fpe settings (/fpe:3 /nocheck) and /NOPTIMIZE.
  !             (ISNAN is an Elemental Intrinsic Function)
  !             Options /fpe:0 /traceback will cause this
  !             program to stop with error message,
  !          II) AIX XLF90 without optimization.
  !              (ISNAN is part of a BOS Runtime C library;
  !               thus ISNAN must be declared LOGICAL.)
  !
  ! Author: hdkLESS at SPAM psu dot edu
  ! Date: March, 2002
  !
  ! Output:
  ! Y = Plus Infinity
  ! i=           1  Y= Infinity
  !
  ! Y = Minus Infinity
  ! i=           2  Y= -Infinity
  !
  ! Y = Minus Zero
  ! i=           3  Y=  0.0000000E+00
  !
  ! 1) Y is a NaN
  ! 2) Y is a NaN
  ! 3) Y is a NaN
  ! i=           4  Y= NaN
  !
  ! 1) Y is a NaN
  ! 2) Y is a NaN
  ! 3) Y is a NaN
  ! i=           5  Y= NaN
  !
  ! References:
  ! http://www.psc.edu/general/software/packages/ieee/ieee.html
  ! http://homepages.borland.com/efg2lab/Mathematics/NaN.htm
  !
  use bl_ieee_module
  implicit none
  real :: y(6)
  integer :: iy(6)
  real    ::  PInf, MInf, MZero, DivNan
  integer ::  iPInf, iMInf, iMZero, iDivNan, iNaN1, iNaN2
  equivalence(pinf,ipinf),(iminf,minf),(imzero,mzero),(idivnan,divnan)
  equivalence(iy,y)
  data iPInf/B'01111111100000000000000000000000'/    ! +Infinity
  data iMInf/B'11111111100000000000000000000000'/    ! -Infinity
  data iMZero/B'10000000000000000000000000000000'/   ! -0
  data iNaN1/B'01111111100000100000000000000000'/     ! NaN
  data iNan2/B'11111111100100010001001010101010'/     ! NaN
  integer :: i

  iy(1) = iPInf
  iy(2) = iMInf
  iy(3) = iMZero
  iy(4) = iNaN1
  iy(5) = iNaN2
  y(6) = DivNan/DivNan
  Do i = 1, 6
     if ( y(i) == PInf )  write(*,*) 'Y = Plus Infinity'
     if ( y(i) == MInf )  write(*,*) 'Y = Minus Infinity'
     if ( y(i) == Mzero ) write(*,*) 'Y = Minus Zero'
     if ( i == 6 )        write(*,*) 'Y = 0/0'

     write(*,*) 'i = ', i, ', Y = ', y(i)
     ! Check if y(i) is a NaN.
     ! EQV -> true iff both A and B are True or iff both A and B are False.
     if ( (y(i) > 0.0) .EQV. (y(i) <= 0.0) ) write(*,*) '1) Y is a NaN'
     if ( y(i) /= y(i) )  write(*,*) '2) Y is a NaN'
     ! If ISNAN is not available for a specific compiler, comment the
     ! following line.
     if ( bl_isnan(y(i)) ) write(*,*) '3) Y is a NaN'
     write(*,*) ' '
  End Do

end subroutine t_nantest
