! A C-program for MT19937, with initialization improved 2002/1/26.
! Coded by Takuji Nishimura and Makoto Matsumoto.

! Code converted to Fortran 95 by Josi Rui Faustino de Sousa
! Date: 2002-02-01

! Before using, initialize the state by using init_genrand(seed)
! or init_by_array(init_key, key_length).

! This library is free software.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

! Copyright (C) 1997, 2002 Makoto Matsumoto and Takuji Nishimura.
! Any feedback is very welcome.
! http://www.math.keio.ac.jp/matumoto/emt.html
! email: matumoto@math.keio.ac.jp

module mt19937_module

  implicit none

  private
  public :: init_genrand, init_by_array
  public :: genrand_bit
  public :: genrand_int
  public :: genrand_int32, genrand_int31
  public :: genrand_real1, genrand_real2, genrand_real3, genrand_res53

  public :: mt_init_genrand, mt_init_by_array
  public :: mt_genrand_bit
  public :: mt_genrand_int
  public :: mt_genrand_int32, mt_genrand_int31
  public :: mt_genrand_real1, mt_genrand_real2, mt_genrand_real3, mt_genrand_res53
  public :: mt_validate, build_random_boxarray

  ! genrand_init
  ! Internal RNG state initialization function accepts either an genrand_intg integer
  ! or a vector returns genrand_state if called with no arguments uses a default seed.
  !   genrand_int32
  !     Subroutine takes two arguments the state to use and a scalar or an array whose
  !     elements are random integer on [0,0xffffffff] interval. Note that this numbers
  !     will positive/negative since there are no unsigned integers in F95
  !   genrand_int31
  !     Subroutine takes two arguments the state to use and a scalar or an array whose
  !     elements are random integer on [0,0x7fffffff] interval.
  !   genrand_real1
  !     Subroutine takes two arguments the state to use and a scalar or an array whose
  !     elements are random real on [0,1] interval.
  !   genrand_real2
  !     Subroutine takes two arguments the state to use and a scalar or an array whose
  !     elements are random numbers on [0,1) interval.
  !   genrand_real3
  !     Subroutine takes two arguments the state to use and a scalar or an array whose
  !     elements are random numbers on (0,1) interval.
  !   genrand_res53
  !     Subroutine takes two arguments the state to use and a scalar or an array whose
  !     elements are random numbers on [0,1) interval with 53-bit resolution.

  ! mt_random_number returns floats,doubles in range [0,1) ala the build random_number

  interface mt_random_number
     module procedure mt_random_number_l_0, mt_random_number_l_0a
     module procedure mt_random_number_l_1, mt_random_number_l_1a
     module procedure mt_random_number_l_2, mt_random_number_l_2a
     module procedure mt_random_number_l_3, mt_random_number_l_3a
     module procedure mt_random_number_l_4, mt_random_number_l_4a
     module procedure mt_random_number_l_5, mt_random_number_l_5a
     module procedure mt_random_number_l_6, mt_random_number_l_6a
     module procedure mt_random_number_l_7, mt_random_number_l_7a

     module procedure mt_random_number_i_0, mt_random_number_i_0a
     module procedure mt_random_number_i_1, mt_random_number_i_1a
     module procedure mt_random_number_i_2, mt_random_number_i_2a
     module procedure mt_random_number_i_3, mt_random_number_i_3a
     module procedure mt_random_number_i_4, mt_random_number_i_4a
     module procedure mt_random_number_i_5, mt_random_number_i_5a
     module procedure mt_random_number_i_6, mt_random_number_i_6a
     module procedure mt_random_number_i_7, mt_random_number_i_7a

     module procedure mt_random_number_d_0, mt_random_number_d_0a
     module procedure mt_random_number_d_1, mt_random_number_d_1a
     module procedure mt_random_number_d_2, mt_random_number_d_2a
     module procedure mt_random_number_d_3, mt_random_number_d_3a
     module procedure mt_random_number_d_4, mt_random_number_d_4a
     module procedure mt_random_number_d_5, mt_random_number_d_5a
     module procedure mt_random_number_d_6, mt_random_number_d_6a
     module procedure mt_random_number_d_7, mt_random_number_d_7a

     module procedure mt_random_number_r_0, mt_random_number_r_0a
     module procedure mt_random_number_r_1, mt_random_number_r_1a
     module procedure mt_random_number_r_2, mt_random_number_r_2a
     module procedure mt_random_number_r_3, mt_random_number_r_3a
     module procedure mt_random_number_r_4, mt_random_number_r_4a
     module procedure mt_random_number_r_5, mt_random_number_r_5a
     module procedure mt_random_number_r_6, mt_random_number_r_6a
     module procedure mt_random_number_r_7, mt_random_number_r_7a
  end interface mt_random_number

  interface mt_random_seed
     module procedure mt_random_seed_0
     module procedure mt_random_seed_a
  end interface

  public mt_random_number
  public mt_random_seed

  integer, parameter  :: intg = selected_int_kind( 9 )
  integer, parameter  :: flot = selected_real_kind( 6, 37 )
  integer, parameter  :: dobl = selected_real_kind( 15, 307 )

  integer, parameter :: wi = intg
  integer, parameter :: wr = dobl
  integer, parameter :: wf = flot


  ! Period parameters
  integer( kind = wi ), parameter :: n = 624_wi
  integer( kind = wi ), parameter :: m = 397_wi
  integer, parameter :: hbs = bit_size( n ) / 2
  integer, parameter :: qbs = hbs / 2
  integer, parameter :: tbs = 3 * qbs

  type, public :: mt19937
     integer( kind = wi ) :: mt(n)
     logical( kind = wi ) :: mtinit = .false._wi
     integer( kind = wi ) :: mti = n + 1_wi
  end type mt19937

  type(mt19937), save :: the_mt

contains

  !! Places a box into the world.
  !! The center of the is uniformly distributed in the world.
  !! The size of the box is uniformly distrubed between MN and MX
  !! The resulting box is intersected with the world
  function box_random_box(world, mn, mx, mt) result(r)

    use box_module
    type(mt19937), intent(inout), optional :: mt
    type(box) :: r
    type(box), intent(in) :: world
    integer, intent(in) :: mn, mx
    real(kind=dp_t) :: aa(get_dim(world),2)
    integer ::  spot(get_dim(world))
    integer :: hwide(get_dim(world))

    if ( present(mt) ) then
       call mt_random_number(mt, aa)
    else
       call mt_random_number(aa)
    end if
    hwide = int((mn + aa(:,1)*(mx-mn))/2)
    spot = int(lwb(world) + aa(:,2)*(upb(world)-lwb(world)))
    call build(r, spot-hwide, spot+hwide)
    r = intersection(r, world)

  end function box_random_box

  !! Makes an array of random boxes
  subroutine make_random_boxes_bv(bxs, world, mn, mx, mt)

    use box_module
    type(box), dimension(:), intent(out) :: bxs
    type(mt19937), intent(inout), optional :: mt
    integer, intent(in) :: mn, mx
    type(box), intent(in) :: world
    integer :: i

    do i = 1, size(bxs)
       bxs(i) = box_random_box(world, mn, mx, mt)
    end do

  end subroutine make_random_boxes_bv

  subroutine build_random_boxarray(ba, pd, nb, mn, mx, bf)
    use boxarray_module
    implicit none
    type(boxarray), intent(out) :: ba
    integer, intent(in) :: nb, mn, mx, bf
    type(box), intent(in) :: pd
    type(box) :: tpd
    type(box) :: bxs(nb)

    integer :: tmn, tmx

    tpd = coarsen(pd, bf)
    tmn = max(mn/bf, 1)
    tmx = max(mx/bf, 1)

    call make_random_boxes_bv(bxs, tpd, tmn, tmx)
    call boxarray_build_v(ba, bxs)
    call boxarray_refine(ba, bf)

    call boxarray_to_domain(ba)

  end subroutine build_random_boxarray

  elemental function uiadd( a, b ) result( c )
    integer( kind = wi ), intent( in )  :: a, b
    integer( kind = wi )  :: c
    integer( kind = wi )  :: a1, a2, b1, b2, s1, s2

    a1 = ibits( a, 0, hbs )
    a2 = ibits( a, hbs, hbs )
    b1 = ibits( b, 0, hbs )
    b2 = ibits( b, hbs, hbs )
    s1 = a1 + b1
    s2 = a2 + b2 + ibits( s1, hbs, hbs )
    c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )
  end function uiadd

  elemental function uisub( a, b ) result( c )
    integer( kind = wi ), intent( in )  :: a, b
    integer( kind = wi )  :: c
    integer( kind = wi )  :: a1, a2, b1, b2, s1, s2

    a1 = ibits( a, 0, hbs )
    a2 = ibits( a, hbs, hbs )
    b1 = ibits( b, 0, hbs )
    b2 = ibits( b, hbs, hbs )
    s1 = a1 - b1
    s2 = a2 - b2 + ibits( s1, hbs, hbs )
    c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )
  end function uisub

  elemental function uimlt( a, b ) result( c )
    integer( kind = wi ), intent( in )  :: a, b
    integer( kind = wi )  :: c
    integer( kind = wi )  :: a0, a1, a2, a3
    integer( kind = wi )  :: b0, b1, b2, b3
    integer( kind = wi )  :: p0, p1, p2, p3

    a0 = ibits( a, 0, qbs )
    a1 = ibits( a, qbs, qbs )
    a2 = ibits( a, hbs, qbs )
    a3 = ibits( a, tbs, qbs )
    b0 = ibits( b, 0, qbs )
    b1 = ibits( b, qbs, qbs )
    b2 = ibits( b, hbs, qbs )
    b3 = ibits( b, tbs, qbs )
    p0 = a0 * b0
    p1 = a1 * b0 + a0 * b1 + ibits( p0, qbs, tbs )
    p2 = a2 * b0 + a1 * b1 + a0 * b2 + ibits( p1, qbs, tbs )
    p3 = a3 * b0 + a2 * b1 + a1 * b2 + a0 * b3 + ibits( p2, qbs, tbs )
    c  = ior( ishft( p1, qbs ), ibits( p0, 0, qbs ) )
    c  = ior( ishft( p2, hbs ), ibits( c, 0, hbs ) )
    c  = ior( ishft( p3, tbs ), ibits( c, 0, tbs ) )
  end function uimlt

  ! initializes mt[N] with a seed
  subroutine init_genrand( s )
    integer( kind = wi ), intent( in )  :: s
    call mt_init_genrand(the_mt, s)
  end subroutine init_genrand

  subroutine mt_init_genrand( mt, s )
    type(mt19937), intent(out) :: mt
    integer( kind = wi ), intent( in )  :: s
    integer( kind = wi )  :: i, mult_a

    data mult_a /z'6C078965'/

    mt%mtinit = .true._wi
    mt%mt(1) = ibits( s, 0, 32 )
    do i = 2, n, 1
       mt%mt(i) = ieor( mt%mt(i-1), ishft( mt%mt(i-1), -30 ) )
       mt%mt(i) = uimlt( mt%mt(i), mult_a )
       mt%mt(i) = uiadd( mt%mt(i),  i-1_wi )
       ! See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
       ! In the previous versions, MSBs of the seed affect
       ! only MSBs of the array mt%mt[].
       ! 2002/01/09 modified by Makoto Matsumoto
       mt%mt(i) = ibits( mt%mt(i), 0, 32 )
       ! for >32 bit machines
    end do

  end subroutine mt_init_genrand

  ! initialize by an array with array-length
  ! init_key is the array for initializing keys
  ! key_length is its length
  subroutine init_by_array( init_key )
    integer( kind = wi ), intent( in )  :: init_key(:)
    call mt_init_by_array(the_mt, init_key)
  end subroutine init_by_array

  subroutine mt_init_by_array( mt, init_key )
    type(mt19937), intent(out) :: mt
    integer( kind = wi ), intent( in )  :: init_key(:)
    integer( kind = wi )  :: i, j, k, tp, key_length
    integer( kind = wi )  :: seed_d, mult_a, mult_b, msb1_d

    data seed_d /z'12BD6AA'/
    data mult_a /z'19660D'/
    data mult_b /z'5D588B65'/
    data msb1_d /z'80000000'/

    key_length = size( init_key, dim=1 )
    call mt_init_genrand( mt, seed_d )
    i = 2_wi
    j = 1_wi
    do k = max( n, key_length ), 1, -1
       tp = ieor( mt%mt(i-1), ishft( mt%mt(i-1), -30 ) )
       tp = uimlt( tp, mult_a )
       mt%mt(i) = ieor( mt%mt(i), tp )
       mt%mt(i) = uiadd( mt%mt(i), uiadd( init_key(j), j-1_wi  ) ) ! non linear
       mt%mt(i) = ibits( mt%mt(i), 0, 32 ) ! for WORDSIZE > 32 machines
       i = i + 1_wi
       j = j + 1_wi
       if ( i > n ) then
          mt%mt(1) = mt%mt(n)
          i = 2_wi
       end if
       if ( j > key_length) j = 1_wi
    end do
    do k = n-1, 1, -1
       tp = ieor( mt%mt(i-1), ishft( mt%mt(i-1), -30 ) )
       tp = uimlt( tp, mult_b )
       mt%mt(i) = ieor( mt%mt(i), tp )
       mt%mt(i) = uisub( mt%mt(i), i-1_wi ) ! non linear
       mt%mt(i) = ibits( mt%mt(i), 0, 32 ) ! for WORDSIZE > 32 machines
       i = i + 1_wi
       if ( i > n ) then
          mt%mt(1) = mt%mt(n)
          i = 2_wi
       end if
    end do
    mt%mt(1) = msb1_d      ! MSB is 1; assuring non-zero initial array
  end subroutine mt_init_by_array

  ! generates a random number on [0,0xffffffff]-interval
  function genrand_int32( ) result( y )
    integer( kind = wi )  :: y
    y = mt_genrand_int32(the_mt)
  end function genrand_int32

  function mt_genrand_int32( mt ) result( y )
    type(mt19937), intent(inout) :: mt
    integer( kind = wi )  :: y
    integer( kind = wi )  :: kk
    integer( kind = wi )  :: seed_d, matrix_a, matrix_b, temper_a, temper_b

    data seed_d   /z'5489'/
    data matrix_a /z'9908B0DF'/
    data matrix_b /z'0'/
    data temper_a /z'9D2C5680'/
    data temper_b /z'EFC60000'/

    if ( mt%mti > n ) then      ! generate N words at one time
       ! if init_genrand() has not been called, a default initial seed is used
       if ( .not. mt%mtinit ) call mt_init_genrand( mt, seed_d )
       do kk = 1, n-m, 1
          y = ibits( mt%mt(kk+1), 0, 31 )
          call mvbits( mt%mt(kk), 31, 1, y, 31 )
          if ( btest( y, 0 ) ) then
             mt%mt(kk) = ieor( ieor( mt%mt(kk+m), ishft( y, -1 ) ), matrix_a )
          else
             mt%mt(kk) = ieor( ieor( mt%mt(kk+m), ishft( y, -1 ) ), matrix_b )
          end if
       end do
       do kk = n-m+1, n-1, 1
          y = ibits( mt%mt(kk+1), 0, 31 )
          call mvbits( mt%mt(kk), 31, 1, y, 31 )
          if ( btest( y, 0 ) ) then
             mt%mt(kk) = ieor(ieor( mt%mt(kk+m-n), ishft( y, -1 ) ), matrix_a)
          else
             mt%mt(kk) = ieor(ieor( mt%mt(kk+m-n), ishft( y, -1 ) ), matrix_b)
          end if
       end do
       y = ibits( mt%mt(1), 0, 31 )
       call mvbits( mt%mt(n), 31, 1, y, 31 )
       if ( btest( y, 0 ) ) then
          mt%mt(kk) = ieor( ieor( mt%mt(m), ishft( y, -1 ) ), matrix_a )
       else
          mt%mt(kk) = ieor( ieor( mt%mt(m), ishft( y, -1 ) ), matrix_b )
       end if
       mt%mti = 1_wi
    end if
    y = mt%mt(mt%mti)
    mt%mti = mt%mti + 1_wi
    ! Tempering
    y = ieor( y, ishft( y, -11) )
    y = ieor( y, iand( ishft( y, 7 ), temper_a ) )
    y = ieor( y, iand( ishft( y, 15 ), temper_b ) )
    y = ieor( y, ishft( y, -18 ) )
  end function mt_genrand_int32

  ! generates a random number on [0,0x7fffffff]-interval
  function genrand_int31( ) result( i )
    integer( kind = wi )  :: i
    i = mt_genrand_int31(the_mt)
  end function genrand_int31
  function mt_genrand_int31( mt ) result( i )
    type(mt19937), intent(inout) :: mt
    integer( kind = wi )  :: i
    i = ishft( mt_genrand_int32( mt ), -1 )
  end function mt_genrand_int31

  function genrand_int(natural) result(i)
    integer( kind = wi ) :: i
    logical, intent(in), optional :: natural
    i = mt_genrand_int(the_mt, natural)
  end function genrand_int
  function mt_genrand_int(mt, natural) result(i)
    type(mt19937), intent(inout) :: mt
    integer( kind = wi ) :: i
    logical, intent(in), optional :: natural
    if ( present(natural) ) then
       if ( natural ) then
          i = mt_genrand_int31(mt)
       else
          i = mt_genrand_int32(mt)
       end if
    else
       i = mt_genrand_int32(mt)
    end if
  end function mt_genrand_int

  function genrand_bit() result(l)
    logical :: l
    l = mt_genrand_bit(the_mt)
  end function genrand_bit
  function mt_genrand_bit(mt) result(l)
    type(mt19937), intent(inout) :: mt
    logical :: l
    l = mod(mt_genrand_int32(mt),2_wi) == 0
  end function mt_genrand_bit

  ! generates a random number on [0,1]-real-interval
  function genrand_real1( ) result( r )
    real( kind = wr )  :: r
    r = mt_genrand_real1(the_mt)
  end function genrand_real1

  function mt_genrand_real1( mt ) result( r )
    type(mt19937), intent(inout) :: mt
    real( kind = wr )  :: r
    integer( kind = wi )  :: a, a1, a0

    a = mt_genrand_int32( mt )
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    r = real( a0, kind = wr ) / 4294967295.0_wr
    r = real( a1, kind = wr ) * ( 65536.0_wr / 4294967295.0_wr ) + r
    ! divided by 2^32-1
  end function mt_genrand_real1

  ! generates a random number on [0,1)-real-interval
  function genrand_real2( ) result( r )
    real( kind = wr )  :: r
    r = mt_genrand_real2(the_mt)
  end function genrand_real2

  function mt_genrand_real2( mt ) result( r )
    type(mt19937), intent(inout) :: mt
    real( kind = wr )  :: r
    integer( kind = wi )  :: a, a1, a0

    a = mt_genrand_int32( mt )
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    r = real( a0, kind = wr ) / 4294967296.0_wr
    r = real( a1, kind = wr ) / 65536.0_wr + r
    ! divided by 2^32
  end function mt_genrand_real2

  ! generates a random number on (0,1)-real-interval
  function genrand_real3( ) result( r )
    real( kind = wr )  :: r
    r = mt_genrand_real3(the_mt)
  end function genrand_real3

  function mt_genrand_real3( mt ) result( r )
    type(mt19937), intent(inout) :: mt
    real( kind = wr )  :: r
    integer( kind = wi )  :: a, a1, a0

    a = mt_genrand_int32( mt )
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    r = ( real( a0, kind = wr ) + 0.5_wr ) / 4294967296.0_wr
    r = real( a1, kind = wr ) / 65536.0_wr + r
    ! divided by 2^32
  end function mt_genrand_real3

  !! returns various kinds of real results with intervals
  !! selected by lo_open, hi_open
  function genrand_real(lo_open, hi_open) result(r)
    real( kind = wr )  :: r
    logical, intent(in), optional :: lo_open, hi_open
    r = mt_genrand_real(the_mt, lo_open, hi_open)
  end function genrand_real
  function mt_genrand_real(mt, lo_open, hi_open) result(r)
    use bl_error_module
    real( kind = wr )  :: r
    type(mt19937), intent(inout) :: mt
    logical, intent(in), optional :: lo_open, hi_open
    logical :: lop, hop
    lop = .false.; if ( present(lo_open) ) lop = lo_open
    hop = .true.;  if ( present(hi_open) ) hop = hi_open
    if ( hop .and. lop ) then
       r = mt_genrand_real3(mt)
    else if ( hop .and. .not. lop ) then
       r = mt_genrand_real2(mt)
    else if ( .not. hop .and. .not. lop) then
       r = mt_genrand_real1(mt)
    else
       call bl_error("MT_GENRAND_REAL: lo_open hi_closed: fails")
    end if
  end function mt_genrand_real

  subroutine mt_random_seed_0(size, put, get)
    integer, intent(out), optional :: size
    integer, intent(in), optional :: put(:)
    integer, intent(out), optional :: get(:)
    call mt_random_seed_a(the_mt, size, put, get)
  end subroutine mt_random_seed_0

  subroutine mt_random_seed_a(mt, size, put, get)
    use bl_error_module
    type(mt19937), intent(inout) :: mt
    integer, intent(out), optional :: size
    integer, intent(in), optional :: put(:)
    integer, intent(out), optional :: get(:)
    integer(kind=wi) :: seed_d
    data seed_d /z'5489'/
    if ( .not. present(size) .and. .not. present(put) &
         .and. .not. present(get) ) then
       call mt_init_genrand(mt, seed_d)
    else if ( present(size) .and. .not. present(put) &
         .and. .not. present(get) ) then
       size = n + 1
    else if ( .not. present(size) .and. present(put) &
         .and. .not. present(get) ) then
       if ( ubound(put,dim=1) < n + 1 ) &
          call bl_error("MT_RANDOM_SEED: put array too small: ", &
          ubound(put,dim=1))
       mt%mt = put(1:n)
       mt%mti = put(n+1)
       mt%mtinit = .true._wi
    else if ( .not. present(size) .and. .not. present(put) &
         .and. present(get) ) then
       if ( .not. mt%mtinit) then
          call mt_init_genrand(mt, seed_d)
       end if
       if ( ubound(get,dim=1) < n + 1 ) &
          call bl_error("MT_RANDOM_SEED: get array too small: ", &
          ubound(get,dim=1))
       get(1:n) = mt%mt
       get(n+1) = mt%mti
    else
       call bl_error("MT_RANDOM_SEED: args either SIZE, GET, PUT")
    end if
  end subroutine mt_random_seed_a

  subroutine mt_random_number_l_0a(harvest)
    logical, intent(out) :: harvest
    call mt_random_number_l_0(the_mt, harvest)
  end subroutine mt_random_number_l_0a
  subroutine mt_random_number_l_0(mt, harvest)
    type(mt19937), intent(inout) :: mt
    logical, intent(out) :: harvest
    harvest = mt_genrand_bit(mt)
  end subroutine mt_random_number_l_0

  subroutine mt_random_number_l_1a(harvest)
    logical, intent(out) :: harvest(:)
    call mt_random_number_l_1(the_mt, harvest)
  end subroutine mt_random_number_l_1a
  subroutine mt_random_number_l_1(mt, harvest)
    type(mt19937), intent(inout) :: mt
    logical, intent(out) :: harvest(:)
    integer :: i
    do i = 1, size(harvest)
       harvest(i) = mt_genrand_bit(mt)
    end do
  end subroutine mt_random_number_l_1

  subroutine mt_random_number_l_2a(harvest)
    logical, intent(out) :: harvest(:,:)
    call mt_random_number_l_2(the_mt, harvest)
  end subroutine mt_random_number_l_2a
  subroutine mt_random_number_l_2(mt, harvest)
    type(mt19937), intent(inout) :: mt
    logical, intent(out) :: harvest(:,:)
    integer :: i,j
    do j = 1, size(harvest,dim=2)
       do i = 1, size(harvest, dim=1)
          harvest(i,j) = mt_genrand_bit(mt)
       end do
    end do
  end subroutine mt_random_number_l_2

  subroutine mt_random_number_l_3a(harvest)
    logical, intent(out) :: harvest(:,:,:)
    call mt_random_number_l_3(the_mt, harvest)
  end subroutine mt_random_number_l_3a
  subroutine mt_random_number_l_3(mt, harvest)
    type(mt19937), intent(inout) :: mt
    logical, intent(out) :: harvest(:,:,:)
    integer :: i,j,k
    do k = 1, size(harvest,dim=3)
       do j = 1, size(harvest,dim=2)
          do i = 1, size(harvest, dim=1)
             harvest(i,j,k) = mt_genrand_bit(mt)
          end do
       end do
    end do
  end subroutine mt_random_number_l_3

  subroutine mt_random_number_l_4a(harvest)
    logical, intent(out) :: harvest(:,:,:,:)
    call mt_random_number_l_4(the_mt, harvest)
  end subroutine mt_random_number_l_4a
  subroutine mt_random_number_l_4(mt, harvest)
    type(mt19937), intent(inout) :: mt
    logical, intent(out) :: harvest(:,:,:,:)
    integer :: i,j,k,l
    do l = 1, size(harvest,dim=4)
       do k = 1, size(harvest,dim=3)
          do j = 1, size(harvest,dim=2)
             do i = 1, size(harvest, dim=1)
                harvest(i,j,k,l) = mt_genrand_bit(mt)
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_l_4

  subroutine mt_random_number_l_5a(harvest)
    logical, intent(out) :: harvest(:,:,:,:,:)
    call mt_random_number_l_5(the_mt, harvest)
  end subroutine mt_random_number_l_5a
  subroutine mt_random_number_l_5(mt, harvest)
    type(mt19937), intent(inout) :: mt
    logical, intent(out) :: harvest(:,:,:,:,:)
    integer :: i,j,k,l,m
    do m = 1, size(harvest,dim=5)
       do l = 1, size(harvest,dim=4)
          do k = 1, size(harvest,dim=3)
             do j = 1, size(harvest,dim=2)
                do i = 1, size(harvest, dim=1)
                   harvest(i,j,k,l,m) = mt_genrand_bit(mt)
                end do
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_l_5

  subroutine mt_random_number_l_6a(harvest)
    logical, intent(out) :: harvest(:,:,:,:,:,:)
    call mt_random_number_l_6(the_mt, harvest)
  end subroutine mt_random_number_l_6a
  subroutine mt_random_number_l_6(mt, harvest)
    type(mt19937), intent(inout) :: mt
    logical, intent(out) :: harvest(:,:,:,:,:,:)
    integer :: i,j,k,l,m,n
    do n = 1, size(harvest,dim=6)
       do m = 1, size(harvest,dim=5)
          do l = 1, size(harvest,dim=4)
             do k = 1, size(harvest,dim=3)
                do j = 1, size(harvest,dim=2)
                   do i = 1, size(harvest, dim=1)
                      harvest(i,j,k,l,m,n) = mt_genrand_bit(mt)
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_l_6

  subroutine mt_random_number_l_7a(harvest)
    logical, intent(out) :: harvest(:,:,:,:,:,:,:)
    call mt_random_number_l_7(the_mt, harvest)
  end subroutine mt_random_number_l_7a
  subroutine mt_random_number_l_7(mt, harvest)
    type(mt19937), intent(inout) :: mt
    logical, intent(out) :: harvest(:,:,:,:,:,:,:)
    integer :: i,j,k,l,m,n,o
    do o = 1, size(harvest,dim=7)
       do n = 1, size(harvest,dim=6)
          do m = 1, size(harvest,dim=5)
             do l = 1, size(harvest,dim=4)
                do k = 1, size(harvest,dim=3)
                   do j = 1, size(harvest,dim=2)
                      do i = 1, size(harvest, dim=1)
                         harvest(i,j,k,l,m,n,o) = mt_genrand_bit(mt)
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_l_7

  subroutine mt_random_number_i_0a(harvest, natural)
    integer, intent(out) :: harvest
    logical, intent(in), optional :: natural
    call mt_random_number_i_0(the_mt, harvest, natural)
  end subroutine mt_random_number_i_0a
  subroutine mt_random_number_i_0(mt, harvest, natural)
    type(mt19937), intent(inout) :: mt
    logical, intent(in), optional :: natural
    integer, intent(out) :: harvest
    harvest = mt_genrand_int(mt, natural)
  end subroutine mt_random_number_i_0

  subroutine mt_random_number_i_1a(harvest, natural)
    integer, intent(out) :: harvest(:)
    logical, intent(in), optional :: natural
    call mt_random_number_i_1(the_mt, harvest, natural)
  end subroutine mt_random_number_i_1a
  subroutine mt_random_number_i_1(mt, harvest, natural)
    type(mt19937), intent(inout) :: mt
    logical, intent(in), optional :: natural
    integer, intent(out) :: harvest(:)
    integer :: i
    do i = 1, size(harvest)
       harvest(i) = mt_genrand_int(mt, natural)
    end do
  end subroutine mt_random_number_i_1

  subroutine mt_random_number_i_2a(harvest, natural)
    integer, intent(out) :: harvest(:,:)
    logical, intent(in), optional :: natural
    call mt_random_number_i_2(the_mt, harvest, natural)
  end subroutine mt_random_number_i_2a
  subroutine mt_random_number_i_2(mt, harvest, natural)
    type(mt19937), intent(inout) :: mt
    logical, intent(in), optional :: natural
    integer, intent(out) :: harvest(:,:)
    integer :: i,j
    do j = 1, size(harvest,dim=2)
       do i = 1, size(harvest, dim=1)
          harvest(i,j) = mt_genrand_int(mt, natural)
       end do
    end do
  end subroutine mt_random_number_i_2

  subroutine mt_random_number_i_3a(harvest, natural)
    integer, intent(out) :: harvest(:,:,:)
    logical, intent(in), optional :: natural
    call mt_random_number_i_3(the_mt, harvest, natural)
  end subroutine mt_random_number_i_3a
  subroutine mt_random_number_i_3(mt, harvest, natural)
    type(mt19937), intent(inout) :: mt
    logical, intent(in), optional :: natural
    integer, intent(out) :: harvest(:,:,:)
    integer :: i,j,k
    do k = 1, size(harvest,dim=3)
       do j = 1, size(harvest,dim=2)
          do i = 1, size(harvest, dim=1)
             harvest(i,j,k) = mt_genrand_int(mt, natural)
          end do
       end do
    end do
  end subroutine mt_random_number_i_3

  subroutine mt_random_number_i_4a(harvest, natural)
    integer, intent(out) :: harvest(:,:,:,:)
    logical, intent(in), optional :: natural
    call mt_random_number_i_4(the_mt, harvest, natural)
  end subroutine mt_random_number_i_4a
  subroutine mt_random_number_i_4(mt, harvest, natural)
    type(mt19937), intent(inout) :: mt
    logical, intent(in), optional :: natural
    integer, intent(out) :: harvest(:,:,:,:)
    integer :: i,j,k,l
    do l = 1, size(harvest,dim=4)
       do k = 1, size(harvest,dim=3)
          do j = 1, size(harvest,dim=2)
             do i = 1, size(harvest, dim=1)
                harvest(i,j,k,l) = mt_genrand_int(mt, natural)
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_i_4

  subroutine mt_random_number_i_5a(harvest, natural)
    integer, intent(out) :: harvest(:,:,:,:,:)
    logical, intent(in), optional :: natural
    call mt_random_number_i_5(the_mt, harvest, natural)
  end subroutine mt_random_number_i_5a
  subroutine mt_random_number_i_5(mt, harvest, natural)
    type(mt19937), intent(inout) :: mt
    logical, intent(in), optional :: natural
    integer, intent(out) :: harvest(:,:,:,:,:)
    integer :: i,j,k,l,m
    do m = 1, size(harvest,dim=5)
       do l = 1, size(harvest,dim=4)
          do k = 1, size(harvest,dim=3)
             do j = 1, size(harvest,dim=2)
                do i = 1, size(harvest, dim=1)
                   harvest(i,j,k,l,m) = mt_genrand_int(mt, natural)
                end do
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_i_5

  subroutine mt_random_number_i_6a(harvest, natural)
    integer, intent(out) :: harvest(:,:,:,:,:,:)
    logical, intent(in), optional :: natural
    call mt_random_number_i_6(the_mt, harvest, natural)
  end subroutine mt_random_number_i_6a
  subroutine mt_random_number_i_6(mt, harvest, natural)
    type(mt19937), intent(inout) :: mt
    logical, intent(in), optional :: natural
    integer, intent(out) :: harvest(:,:,:,:,:,:)
    integer :: i,j,k,l,m,n
    do n = 1, size(harvest,dim=6)
       do m = 1, size(harvest,dim=5)
          do l = 1, size(harvest,dim=4)
             do k = 1, size(harvest,dim=3)
                do j = 1, size(harvest,dim=2)
                   do i = 1, size(harvest, dim=1)
                      harvest(i,j,k,l,m,n) = mt_genrand_int(mt, natural)
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_i_6

  subroutine mt_random_number_i_7a(harvest, natural)
    integer, intent(out) :: harvest(:,:,:,:,:,:,:)
    logical, intent(in), optional :: natural
    call mt_random_number_i_7(the_mt, harvest, natural)
  end subroutine mt_random_number_i_7a
  subroutine mt_random_number_i_7(mt, harvest, natural)
    type(mt19937), intent(inout) :: mt
    logical, intent(in), optional :: natural
    integer, intent(out) :: harvest(:,:,:,:,:,:,:)
    integer :: i,j,k,l,m,n,o
    do o = 1, size(harvest,dim=7)
       do n = 1, size(harvest,dim=6)
          do m = 1, size(harvest,dim=5)
             do l = 1, size(harvest,dim=4)
                do k = 1, size(harvest,dim=3)
                   do j = 1, size(harvest,dim=2)
                      do i = 1, size(harvest, dim=1)
                         harvest(i,j,k,l,m,n,o) = mt_genrand_int(mt, natural)
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_i_7

  subroutine mt_random_number_d_0a(harvest)
    real(kind = wr), intent(out) :: harvest
    call mt_random_number_d_0(the_mt, harvest)
  end subroutine mt_random_number_d_0a
  subroutine mt_random_number_d_0(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wr), intent(out) :: harvest
    harvest = mt_genrand_real2(mt)
  end subroutine mt_random_number_d_0

  subroutine mt_random_number_d_1a(harvest)
    real(kind = wr), intent(out) :: harvest(:)
    call mt_random_number_d_1(the_mt, harvest)
  end subroutine mt_random_number_d_1a
  subroutine mt_random_number_d_1(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wr), intent(out) :: harvest(:)
    integer :: i
    do i = 1, size(harvest)
       harvest(i) = mt_genrand_real2(mt)
    end do
  end subroutine mt_random_number_d_1

  subroutine mt_random_number_d_2a(harvest)
    real(kind = wr), intent(out) :: harvest(:,:)
    call mt_random_number_d_2(the_mt, harvest)
  end subroutine mt_random_number_d_2a

  subroutine mt_random_number_d_2(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wr), intent(out) :: harvest(:,:)
    integer :: i,j
    do j = 1, size(harvest,dim=2)
       do i = 1, size(harvest, dim=1)
          harvest(i,j) = mt_genrand_real2(mt)
       end do
    end do
  end subroutine mt_random_number_d_2

  subroutine mt_random_number_d_3a(harvest)
    real(kind = wr), intent(out) :: harvest(:,:,:)
    call mt_random_number_d_3(the_mt, harvest)
  end subroutine mt_random_number_d_3a

  subroutine mt_random_number_d_3(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wr), intent(out) :: harvest(:,:,:)
    integer :: i,j,k
    do k = 1, size(harvest,dim=3)
       do j = 1, size(harvest,dim=2)
          do i = 1, size(harvest, dim=1)
             harvest(i,j,k) = mt_genrand_real2(mt)
          end do
       end do
    end do
  end subroutine mt_random_number_d_3

  subroutine mt_random_number_d_4a(harvest)
    real(kind = wr), intent(out) :: harvest(:,:,:,:)
    call mt_random_number_d_4(the_mt, harvest)
  end subroutine mt_random_number_d_4a

  subroutine mt_random_number_d_4(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wr), intent(out) :: harvest(:,:,:,:)
    integer :: i,j,k,l
    do l = 1, size(harvest,dim=4)
       do k = 1, size(harvest,dim=3)
          do j = 1, size(harvest,dim=2)
             do i = 1, size(harvest, dim=1)
                harvest(i,j,k,l) = mt_genrand_real2(mt)
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_d_4

  subroutine mt_random_number_d_5a(harvest)
    real(kind = wr), intent(out) :: harvest(:,:,:,:,:)
    call mt_random_number_d_5(the_mt, harvest)
  end subroutine mt_random_number_d_5a

  subroutine mt_random_number_d_5(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wr), intent(out) :: harvest(:,:,:,:,:)
    integer :: i,j,k,l,m
    do m = 1, size(harvest,dim=5)
       do l = 1, size(harvest,dim=4)
          do k = 1, size(harvest,dim=3)
             do j = 1, size(harvest,dim=2)
                do i = 1, size(harvest, dim=1)
                   harvest(i,j,k,l,m) = mt_genrand_real2(mt)
                end do
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_d_5

  subroutine mt_random_number_d_6a(harvest)
    real(kind = wr), intent(out) :: harvest(:,:,:,:,:,:)
    call mt_random_number_d_6(the_mt, harvest)
  end subroutine mt_random_number_d_6a

  subroutine mt_random_number_d_6(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wr), intent(out) :: harvest(:,:,:,:,:,:)
    integer :: i,j,k,l,m,n
    do n = 1, size(harvest,dim=6)
       do m = 1, size(harvest,dim=5)
          do l = 1, size(harvest,dim=4)
             do k = 1, size(harvest,dim=3)
                do j = 1, size(harvest,dim=2)
                   do i = 1, size(harvest, dim=1)
                      harvest(i,j,k,l,m,n) = mt_genrand_real2(mt)
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_d_6

  subroutine mt_random_number_d_7a(harvest)
    real(kind = wr), intent(out) :: harvest(:,:,:,:,:,:,:)
    call mt_random_number_d_7(the_mt, harvest)
  end subroutine mt_random_number_d_7a

  subroutine mt_random_number_d_7(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wr), intent(out) :: harvest(:,:,:,:,:,:,:)
    integer :: i,j,k,l,m,n,o
    do o = 1, size(harvest,dim=7)
       do n = 1, size(harvest,dim=6)
          do m = 1, size(harvest,dim=5)
             do l = 1, size(harvest,dim=4)
                do k = 1, size(harvest,dim=3)
                   do j = 1, size(harvest,dim=2)
                      do i = 1, size(harvest, dim=1)
                         harvest(i,j,k,l,m,n,o) = mt_genrand_real2(mt)
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_d_7

  subroutine mt_random_number_r_0a(harvest)
    real(kind = wf), intent(out) :: harvest
    call mt_random_number_r_0(the_mt, harvest)
  end subroutine mt_random_number_r_0a

  subroutine mt_random_number_r_0(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wf), intent(out) :: harvest
    harvest = real(mt_genrand_real2(mt),kind=wf)
  end subroutine mt_random_number_r_0

  subroutine mt_random_number_r_1a(harvest)
    real(kind = wf), intent(out) :: harvest(:)
    call mt_random_number_r_1(the_mt, harvest)
  end subroutine mt_random_number_r_1a

  subroutine mt_random_number_r_1(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wf), intent(out) :: harvest(:)
    integer :: i
    do i = 1, size(harvest)
       harvest(i) = real(mt_genrand_real2(mt),kind=wf)
    end do
  end subroutine mt_random_number_r_1

  subroutine mt_random_number_r_2a(harvest)
    real(kind = wf), intent(out) :: harvest(:,:)
    call mt_random_number_r_2(the_mt, harvest)
  end subroutine mt_random_number_r_2a

  subroutine mt_random_number_r_2(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wf), intent(out) :: harvest(:,:)
    integer :: i,j
    do j = 1, size(harvest,dim=2)
       do i = 1, size(harvest, dim=1)
          harvest(i,j) = real(mt_genrand_real2(mt),kind=wf)
       end do
    end do
  end subroutine mt_random_number_r_2

  subroutine mt_random_number_r_3a(harvest)
    real(kind = wf), intent(out) :: harvest(:,:,:)
    call mt_random_number_r_3(the_mt, harvest)
  end subroutine mt_random_number_r_3a

  subroutine mt_random_number_r_3(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wf), intent(out) :: harvest(:,:,:)
    integer :: i,j,k
    do k = 1, size(harvest,dim=3)
       do j = 1, size(harvest,dim=2)
          do i = 1, size(harvest, dim=1)
             harvest(i,j,k) = real(mt_genrand_real2(mt),kind=wf)
          end do
       end do
    end do
  end subroutine mt_random_number_r_3

  subroutine mt_random_number_r_4a(harvest)
    real(kind = wf), intent(out) :: harvest(:,:,:,:)
    call mt_random_number_r_4(the_mt, harvest)
  end subroutine mt_random_number_r_4a

  subroutine mt_random_number_r_4(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wf), intent(out) :: harvest(:,:,:,:)
    integer :: i,j,k,l
    do l = 1, size(harvest,dim=4)
       do k = 1, size(harvest,dim=3)
          do j = 1, size(harvest,dim=2)
             do i = 1, size(harvest, dim=1)
                harvest(i,j,k,l) = real(mt_genrand_real2(mt),kind=wf)
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_r_4

  subroutine mt_random_number_r_5a(harvest)
    real(kind = wf), intent(out) :: harvest(:,:,:,:,:)
    call mt_random_number_r_5(the_mt, harvest)
  end subroutine mt_random_number_r_5a

  subroutine mt_random_number_r_5(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wf), intent(out) :: harvest(:,:,:,:,:)
    integer :: i,j,k,l,m
    do m = 1, size(harvest,dim=5)
       do l = 1, size(harvest,dim=4)
          do k = 1, size(harvest,dim=3)
             do j = 1, size(harvest,dim=2)
                do i = 1, size(harvest, dim=1)
                   harvest(i,j,k,l,m) = real(mt_genrand_real2(mt),kind=wf)
                end do
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_r_5

  subroutine mt_random_number_r_6a(harvest)
    real(kind = wf), intent(out) :: harvest(:,:,:,:,:,:)
    call mt_random_number_r_6(the_mt, harvest)
  end subroutine mt_random_number_r_6a

  subroutine mt_random_number_r_6(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wf), intent(out) :: harvest(:,:,:,:,:,:)
    integer :: i,j,k,l,m,n
    do n = 1, size(harvest,dim=6)
       do m = 1, size(harvest,dim=5)
          do l = 1, size(harvest,dim=4)
             do k = 1, size(harvest,dim=3)
                do j = 1, size(harvest,dim=2)
                   do i = 1, size(harvest, dim=1)
                      harvest(i,j,k,l,m,n) = real(mt_genrand_real2(mt),kind=wf)
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_r_6

  subroutine mt_random_number_r_7a(harvest)
    real(kind = wf), intent(out) :: harvest(:,:,:,:,:,:,:)
    call mt_random_number_r_7(the_mt, harvest)
  end subroutine mt_random_number_r_7a

  subroutine mt_random_number_r_7(mt, harvest)
    type(mt19937), intent(inout) :: mt
    real(kind = wf), intent(out) :: harvest(:,:,:,:,:,:,:)
    integer :: i,j,k,l,m,n,o
    do o = 1, size(harvest,dim=7)
       do n = 1, size(harvest,dim=6)
          do m = 1, size(harvest,dim=5)
             do l = 1, size(harvest,dim=4)
                do k = 1, size(harvest,dim=3)
                   do j = 1, size(harvest,dim=2)
                      do i = 1, size(harvest, dim=1)
                         harvest(i,j,k,l,m,n,o) = real(mt_genrand_real2(mt),kind=wf)
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine mt_random_number_r_7

  ! generates a random number on [0,1) with 53-bit resolution
  function genrand_res53( )  result( r )
    real( kind = wr )  :: r
    r = mt_genrand_res53(the_mt)
  end function genrand_res53

  function mt_genrand_res53( mt )  result( r )
    type(mt19937), intent(inout) :: mt
    real( kind = wr )  :: r
    integer( kind = wi )  :: a, a0, a1
    integer( kind = wi )  :: b, b0, b1

    a = ishft( mt_genrand_int32( mt ), -5 )
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    b = ishft( mt_genrand_int32( mt ), -6 )
    b0 = ibits( b, 0, hbs )
    b1 = ibits( b, hbs, hbs )
    r = real( a1, kind = wr ) / 2048.0_wr
    r = real( a0, kind = wr ) / 134217728.0_wr + r
    r = real( b1, kind = wr ) / 137438953472.0_wr + r
    r = real( b0, kind = wr ) / 9007199254740992.0_wr + r
  end function mt_genrand_res53
  ! These real versions are due to Isaku Wada, 2002/01/09 added

  !! Returns whether mt rn returns same number as original source.
  function mt_validate()
    logical :: mt_validate
    integer :: i
    integer(kind=wi) :: init(4), ia
    data init /Z'123', Z'234', Z'345', Z'456'/
    type(mt19937) :: mt

    call mt_init_by_array(mt,init)
    do i = 1, 1000
       ia = mt_genrand_int32(mt);
    end do
    mt_validate = ia == (-834941650)
  end function mt_validate

end module mt19937_module

