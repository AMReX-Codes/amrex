
!> Module to create quadrature matrices and accompanying routines for SDC
module SDCquadrature_mod
  use iso_c_binding
  implicit none

  integer,  parameter :: qp = c_long_double   
  integer,  parameter :: dp = c_double
  real(qp), parameter :: eps = 1.0e-23_qp

  !> Quadrature node types
  integer, parameter :: SDC_GAUSS_LOBATTO   = 1
  integer, parameter :: SDC_GAUSS_RADAU     = 2
  integer, parameter :: SDC_CLENSHAW_CURTIS = 3
  integer, parameter :: SDC_UNIFORM         = 4
  integer, parameter :: SDC_GAUSS_LEGENDRE  = 5
  integer, parameter :: SDC_PROPER_NODES    = 2**8
  integer, parameter :: SDC_COMPOSITE_NODES = 2**9
  integer, parameter :: SDC_NO_LEFT         = 2**10
  integer, parameter :: amrex_real = c_double
  private :: qsort_partition

  interface poly_eval
     module procedure poly_eval
     module procedure poly_eval_complex
  end interface

contains

  !>  Subroutine to create quadrature matrices
  subroutine SDC_quadrature(qtype_in,nnodes, nnodes0, nodes, nflags, qmats) bind(C, name="SDC_quadrature")

   integer,    intent(in)  :: qtype_in, nnodes, nnodes0
    real(amrex_real), intent(inout) :: nodes(nnodes)
    real(amrex_real), intent(inout) :: qmats(nnodes,nnodes-1,4)
    integer,    intent(out) :: nflags(nnodes)

    real(amrex_real) :: qmatFE(nnodes-1,nnodes), qmatBE(nnodes-1,nnodes),qmatLU(nnodes-1,nnodes)
    real(amrex_real) :: smat(nnodes-1,nnodes), qmat(nnodes-1,nnodes)
    real(c_long_double) :: qnodes0(nnodes0), qnodes(nnodes), dt
    real(amrex_real)          :: qmat0(nnodes0-1,nnodes0), smat0(nnodes0-1,nnodes0)
    integer             :: flags0(nnodes0)

    integer :: qtype, i, r, refine,m,n
    logical :: composite, proper, no_left

    proper    = btest(qtype_in, 8)
    composite = btest(qtype_in, 9)
    no_left   = btest(qtype_in, 10)
    print *,'proper,composite,no_left',proper,composite,no_left

    qmat = 0
    smat = 0
    flags0 = 0
    nflags = 0

    qtype = qtype_in
    if (proper)    qtype = qtype - SDC_PROPER_NODES
    if (composite) qtype = qtype - SDC_COMPOSITE_NODES
    if (no_left)   qtype = qtype - SDC_NO_LEFT


    if (composite) then

       ! nodes are given by repeating the coarsest set of nodes.  note
       ! that in this case nnodes0 corresponds to the coarsest number
       ! of nodes.

       refine = (nnodes - 1) / (nnodes0 - 1)

       call sdc_qnodes(qnodes0, flags0, qtype, nnodes0)
       call sdc_qmats(qmat0, smat0, qnodes0, qnodes0, flags0, nnodes0, nnodes0)

       dt = 1.q0 / refine
       do i = 1, refine
          r = (i-1)*(nnodes0-1)+1
          qnodes(r:r+nnodes0) = dt * ((i-1) + qnodes0)
          smat(r:r+nnodes0,r:r+nnodes0-1) = smat0 / refine
       end do

       nodes = real(qnodes, amrex_real)

    else if (proper) then

       ! nodes are given by proper quadrature rules

       call sdc_qnodes(qnodes, nflags, qtype, nnodes)
       nodes = real(qnodes, amrex_real)

       call sdc_qmats(qmat, smat, qnodes, qnodes, nflags, nnodes, nnodes)

    else

       ! nodes are given by refining the finest set of nodes.  note
       ! that in this case nnodes0 corresponds to the finest number of
       ! nodes.

       refine = (nnodes0 - 1) / (nnodes - 1)
       call sdc_qnodes(qnodes0, flags0, qtype, nnodes0)

       qnodes = qnodes0(::refine)
       nodes  = real(qnodes, amrex_real)
       nflags = flags0(::refine)

       if (no_left) nflags(1) = 0

       call sdc_qmats(qmat, smat, qnodes, qnodes, nflags, nnodes, nnodes)

    end if


    !  Make forward and backward Euler matrices
    qmatFE=0.0d0
    qmatBE=0.0d0
    do m = 1, nnodes-1
       do n = 1,m
          qmatBE(m,n+1) =  qnodes(n+1)-qnodes(n)
       end do
    end do
    ! Explicit matrix
    do m = 1, nnodes-1
       do n = 1,m
          qmatFE(m,n)   =  qnodes(n+1)-qnodes(n)
       end do
    end do


    if (all(nodes == 0.0d0)) then
       stop 'ERROR: pf_quadrature: invalid SDC nnodes.'
    end if
!!$    print *,'printing matrices in quadrature'
!!$    do m = 1,nnodes-1
!!$       print *,qmat(m,:)
!!$       print *,qmatFE(m,:)
!!$       print *,qmatBE(m,:)
!!$    end do
!!$    print *,'transpose'
    do n = 1,nnodes-1
       do m = 1,nnodes
          qmats(m,n,1)=qmat(n,m)
          qmats(m,n,2)=qmatFE(n,m)
          qmats(m,n,3)=qmatBE(n,m)
          qmats(m,n,4)=qmatLU(n,m)
       end do
    end do

    
  end subroutine SDC_quadrature


  !>  Function to decide if the restriction of the nodes is pointwise, e.g. coarse nodes are every other fine node
  logical function not_proper(flags, node)
    integer(c_int), intent(in) :: flags(:)
    integer,        intent(in) :: node

    not_proper = .not. btest(flags(node), 0)
  end function not_proper



  !> Subroutine to compute high precision quadrature nodes.
  subroutine sdc_qnodes(qnodes, flags, qtype, nnodes) bind(c)
    integer(c_int),       intent(in), value  :: nnodes          !<  Number of nodes
    integer(c_int),       intent(in), value  :: qtype           !<  Type of nodes (see pf_dtype)
    real(c_long_double),  intent(out)        :: qnodes(nnodes)  !<  The computed quadrature nodes
    integer(c_int),       intent(out)        :: flags(nnodes)   !<

    integer :: j, degree
    real(qp), allocatable :: roots(:)
    real(qp), allocatable :: coeffs(:), coeffs2(:)

    real(qp), parameter :: pi = 3.141592653589793115997963468544185161590576171875_qp

    flags = 0

    select case(qtype)

    case (SDC_GAUSS_LEGENDRE)

       degree = nnodes-2
       allocate(roots(degree))
       allocate(coeffs(degree+1))

       call poly_legendre(coeffs, degree)
       call poly_roots(roots, coeffs, degree)

       qnodes(1) = 0.0_qp
       qnodes(2:nnodes-1) = 0.5_qp * (1.0_qp + roots)
       qnodes(nnodes) = 1.0_qp

       deallocate(coeffs)
       deallocate(roots)

       do j = 2, nnodes-1
          flags(j) = ibset(flags(j), 0)
       end do

    case (SDC_GAUSS_LOBATTO)

       degree = nnodes - 1
       allocate(roots(degree-1))
       allocate(coeffs(degree+1))

       call poly_legendre(coeffs, degree)
       call poly_diff(coeffs, degree)
       call poly_roots(roots, coeffs(:degree), degree-1)

!       print*, 'roots = ', roots
!       print*, 'coefs = ', coeffs

       qnodes(1)          = 0.0_qp
       qnodes(2:nnodes-1) = 0.5_qp * (1.0_qp + roots)
       qnodes(nnodes)     = 1.0_qp

       deallocate(coeffs)
       deallocate(roots)

       do j = 1, nnodes
          flags(j) = ibset(flags(j), 0)
       end do

    case (SDC_GAUSS_RADAU)

       degree = nnodes - 1
       allocate(roots(degree))
       allocate(coeffs(degree+1))


       call poly_legendre(coeffs, degree)
       call poly_legendre(coeffs2, degree-1)
       coeffs(:degree) = coeffs(:degree) + coeffs2
       call poly_roots(roots, coeffs, degree)

       qnodes(1)      = 0.0_qp
       do j = 2, nnodes-1
          qnodes(j) = 0.5_qp * (1.0_qp - roots(nnodes+1-j))
       end do
       qnodes(nnodes) = 1.0_qp

       deallocate(coeffs2)
       deallocate(coeffs)
       deallocate(roots)

       do j = 2, nnodes
          flags(j) = ibset(flags(j), 0)
       end do

    case (SDC_CLENSHAW_CURTIS)

       do j = 0, nnodes-1
          qnodes(j+1) = 0.5_qp * (1.0_qp - cos(j * pi / (nnodes-1)))
       end do

       do j = 1, nnodes
          flags(j) = ibset(flags(j), 0)
       end do

    case (SDC_UNIFORM)

       do j = 0, nnodes-1
          qnodes(j+1) = j * (1.0_qp / (nnodes-1))
       end do

       do j = 1, nnodes
          flags(j) = ibset(flags(j), 0)
       end do

    case default
       print *,'qtype = ',qtype
       stop "ERROR: Invalid qtype in sdc_quadrature.f90."

    end select

  end subroutine sdc_qnodes

  !>  Subroutine to compute the quadrature matrices 
  subroutine sdc_qmats(qmat, smat, dst, src, flags, ndst, nsrc) bind(c)
    integer(c_int),      intent(in), value  :: ndst   !<  Number of destination points
    integer(c_int),      intent(in), value  :: nsrc   !<  Number of source points
    real(c_long_double), intent(in)  :: dst(ndst)     !<  Destination points
    real(c_long_double), intent(in)  :: src(nsrc)     !<  Source points
    real(amrex_real),      intent(out) :: qmat(ndst-1, nsrc)  !<  O to dst quadrature weights
    real(amrex_real),      intent(out) :: smat(ndst-1, nsrc)  !< dst(m) to dst(m+1) quadrature weights
    integer(c_int),      intent(in)  :: flags(nsrc)     

    integer  :: i, j, m
    real(qp) :: q, s, den, p(0:nsrc)

    qmat = 0.0_dp
    smat = 0.0_dp

    ! construct qmat and smat
    do i = 1, nsrc

       if (not_proper(flags, i)) cycle

       ! construct interpolating polynomial coefficients
       p    = 0.0_qp
       p(0) = 1.0_qp
       do m = 1, nsrc
          if (not_proper(flags, m) .or. m == i) cycle
          p = eoshift(p, -1) - src(m) * p
       end do

       den = poly_eval(p, nsrc, src(i))

       call poly_int(p, nsrc)

       ! evaluate integrals
       do j = 2, ndst
          q = poly_eval(p, nsrc, dst(j)) - poly_eval(p, nsrc,   0.0_qp)
          s = poly_eval(p, nsrc, dst(j)) - poly_eval(p, nsrc, dst(j-1))

          qmat(j-1,i) = real(q / den, dp)
          smat(j-1,i) = real(s / den, dp)
       end do
    end do

  end subroutine sdc_qmats


  !> Polynomial manipulation routines.
  !!
  !! A polynomial p
  !!
  !!   p(x) = a_n x^n + ... + a_2 x^2 + a_1 x + a_0
  !!
  !! is stored as a Fortran array p(0:n) according to
  !!
  !!   p = [ a_0, a_1, ..., a_n ].
  !!
  
  !> Function to evaluate real polynomial
  real(qp) function poly_eval(p, n, x) result(v) bind(c)
    integer, intent(in), value :: n
    real(qp),       intent(in)        :: p(0:n), x

    integer :: j

    v = p(n)
    do j = n-1, 0, -1
       v = x * v + p(j)
    end do
  end function

  !> Function to evaluate complex polynomial
  complex(qp) function poly_eval_complex(p, n, x) result(v)
    integer, intent(in), value :: n
    real(qp),       intent(in)        :: p(0:n)
    complex(qp),    intent(in)        :: x

    integer :: j

    v = p(n)
    do j = n-1, 0, -1
       v = x * v + p(j)
    end do
  end function



  !> Subroutine to differentiate polynomial (in place)
  subroutine poly_diff(p, n) bind(c)
    integer, intent(in),   value :: n
    real(qp),       intent(inout) :: p(0:n)

    integer  :: j
    real(qp) :: pp(0:n)

    pp = 0.0_qp

    do j = 1, n
       pp(j-1) = j * p(j)
    end do

    p = pp
  end subroutine poly_diff


  !> Subroutine to integrate polynomial (in place)
  subroutine poly_int(p, n) bind(c)
    integer, intent(in),   value :: n
    real(qp),       intent(inout) :: p(0:n)

    integer  :: j
    real(qp) :: pp(0:n)

    pp = 0.0_qp

    do j = 0, n-1
       pp(j+1) = p(j) / (j+1)
    end do

    p = pp
  end subroutine poly_int


  !> Subroutine to compute Legendre polynomial coefficients using Bonnet's recursion formula.
  subroutine poly_legendre(p, n) bind(c)
    integer, intent(in), value :: n
    real(qp),       intent(out)       :: p(0:n)

    real(qp), dimension(0:n) :: p0, p1, p2
    integer :: j, m

    if (n == 0) then
       p = [ 1.0_qp ]
       return
    end if

    if (n == 1) then
       p = [ 0.0_qp, 1.0_qp ]
       return
    end if

    p0 = 0.0_qp; p1 = 0.0_qp; p2 = 0.0_qp

    p0(0) = 1.0_qp
    p1(1) = 1.0_qp

    ! (n + 1) P_{n+1} = (2n + 1) x P_{n} - n P_{n-1}
    do m = 1, n-1
       do j = 1, n
          p2(j) = ( (2*m + 1) * p1(j-1) - m * p0(j) ) / (m + 1)
       end do
       p2(0) = - m * p0(0) / (m + 1)

       p0 = p1
       p1 = p2
    end do

    p = p2
  end subroutine poly_legendre

  !> Subroutine to compute polynomial roots using the Durand-Kerner algorithm.
  !! The roots are assumed to be real.
  subroutine poly_roots(roots, p0, n) bind(c)
    integer,  intent(in), value   :: n
    real(qp),        intent(out)  :: roots(n)
    real(qp),        intent(in)   :: p0(0:n)

    integer     :: i, j, k
    complex(qp) :: num, den, z0(n), z1(n)
    real(qp)    :: p(0:n)

    p = p0 / p0(n)

    ! initial guess
    do i = 1, n
       z0(i) = (0.4_qp, 0.9_qp)**i
    end do

    ! durand-kerner-weierstrass iterations
    z1 = z0
    do k = 1, 100
       do i = 1, n

          ! evaluate poly at z0(i)
          num = poly_eval(p, n, z0(i))

          ! evaluate denominator
          den = 1.0_qp
          do j = 1, n
             if (j == i) cycle
             den = den * (z0(i) - z0(j))
          end do

          ! update
          z0(i) = z0(i) - num / den
       end do

       ! converged?
       if (sum(abs(z0 - z1)) < eps) exit

       z1 = z0
    end do

    roots = real(z0)
    where (abs(roots) < eps) roots = 0.0_qp
    call qsort(roots)

  end subroutine poly_roots

  !> Subroutine to sort (inplace) using the quick sort algorithm.
  !> Adapted from http://www.fortran.com/qsort_c.f95.
  recursive subroutine qsort(a)
    real(qp), intent(inout) :: a(:)
    integer :: iq

    if (size(a) > 1) then
       call qsort_partition(a, iq)
       call qsort(a(:iq-1))
       call qsort(a(iq:))
    end if
  end subroutine qsort

  subroutine qsort_partition(a, marker)
    real(qp), intent(inout) :: a(:)
    integer,  intent(out)   :: marker

    integer  :: i, j
    real(qp) :: temp, x

    x = a(1)
    i = 0
    j = size(a) + 1

    do
       j = j-1
       do
          if (a(j) <= x) exit
          j = j-1
       end do

       i = i+1
       do
          if (a(i) >= x) exit
          i = i+1
       end do

       if (i < j) then
          temp = a(i)
          a(i) = a(j)
          a(j) = temp
       else if (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do
  end subroutine qsort_partition


end module SDCquadrature_mod






