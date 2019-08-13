module variables_module

  use iso_c_binding
  use chemistry_module, only : nspecies

  implicit none

  ! Indices
  integer, save :: irho, imx, imy, imz, iene, iry1
  integer, save :: qrho, qu, qv, qw, qpres, qtemp, qe, qy1, qx1, qh1

  ! Number of conserved and primitive variables
  integer, save :: ncons, nprim

  ! Indices for S3D style first-derivatives
  integer, parameter :: idu=1, idv=2, idw=3, idT=4, idp=5, idX1=6 

  ! Arithmetic constants 
  double precision, parameter :: Zero          = 0.d0
  double precision, parameter :: One           = 1.d0
  double precision, parameter :: OneThird      = 1.d0/3.d0
  double precision, parameter :: TwoThirds     = 2.d0/3.d0
  double precision, parameter :: FourThirds    = 4.d0/3.d0
  double precision, parameter :: OneQuarter    = 1.d0/4.d0
  double precision, parameter :: ThreeQuarters = 3.d0/4.d0

  integer, save, public :: iN2, iO2, iH2, iias, iry_ias

contains

  ! 
  ! Initialize various indices
  !
  subroutine variables_init() bind(c,name='variables_init')

    irho = 1
    imx = 2
    imy = 3
    imz = 4
    iene = 5
    iry1 = iene+1

    ncons = iry1-1 + nspecies

    qrho = 1
    qu = 2
    qv = 3
    qw = 4
    qpres = 5
    qtemp = 6
    qe    = 7
    qy1   = 8
    qx1 = qy1 + nspecies
    qh1 = qx1 + nspecies

    nprim = qh1-1 + nspecies

    iH2 = 1
    iO2 = 2
    iN2 = 9

    iias = iN2
    iry_ias = iry1 + iias - 1
    
  end subroutine variables_init

  function get_num_cons () result(r) bind(c,name='get_num_cons')
    integer :: r
    r = ncons
  end function get_num_cons

  function get_num_prim () result(r) bind(c,name='get_num_prim')
    integer :: r
    r = nprim
  end function get_num_prim

  subroutine ctoprim_3d(tlo, thi, lo, hi, u, q, ngu, ngq) bind(c,name='ctoprim_3d')
    integer, intent(in) :: tlo(3), thi(3), lo(3), hi(3), ngu, ngq
    double precision, intent(in ) :: u(lo(1)-ngu:hi(1)+ngu,lo(2)-ngu:hi(2)+ngu,lo(3)-ngu:hi(3)+ngu,ncons)
    double precision, intent(out) :: q(lo(1)-ngq:hi(1)+ngq,lo(2)-ngq:hi(2)+ngq,lo(3)-ngq:hi(3)+ngq,nprim)
    
    integer :: i, j, k, n, iwrk, ierr
    double precision :: rho, rhoinv, rwrk, X(nspecies), Y(nspecies), h(nspecies), ei, Tt, Pt

    do k = tlo(3),thi(3)
       do j = tlo(2),thi(2)
          do i = tlo(1),thi(1)
             
             rho = u(i,j,k,irho)
             rhoinv = 1.d0/rho
             q(i,j,k,qrho) = rho
             q(i,j,k,qu) = u(i,j,k,imx) * rhoinv
             q(i,j,k,qv) = u(i,j,k,imy) * rhoinv
             q(i,j,k,qw) = u(i,j,k,imz) * rhoinv

             do n=1,nspecies
                Y(n) = u(i,j,k,iry1+n-1) * rhoinv
                q(i,j,k,qy1+n-1) = Y(n)
             end do

             call ckytx(Y, iwrk, rwrk, X)

             do n=1,nspecies
                q(i,j,k,qx1+n-1) = X(n)
             end do

             ei = rhoinv*u(i,j,k,iene) - 0.5d0*(q(i,j,k,qu)**2+q(i,j,k,qv)**2+q(i,j,k,qw)**2)
             q(i,j,k,qe) = ei

             Tt = q(i,j,k,qtemp)
             call get_T_given_ey(ei, Y, iwrk, rwrk, Tt, ierr)
             q(i,j,k,qtemp) = Tt

             call CKPY(rho, Tt, Y, iwrk, rwrk, Pt)
             q(i,j,k,qpres) = Pt

             call ckhms(Tt, iwrk, rwrk, h)

             do n=1,nspecies
                q(i,j,k,qh1+n-1) = h(n)
             end do
          enddo
       enddo
    enddo

  end subroutine ctoprim_3d


  subroutine reset_rho_3d(tlo, thi, lo, hi, u, ng) bind(c,name='reset_rho_3d')
    integer, intent(in) :: tlo(3), thi(3), lo(3), hi(3), ng
    double precision, intent(inout) :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,ncons)

    integer :: i, j, k, n, iryn
    double precision :: rho
    double precision, parameter :: eps = -1.0d-16
    integer          :: idom
    double precision :: rhoy_dom, rhoy_under

    do k = tlo(3),thi(3)
       do j = tlo(2),thi(2)
          do i = tlo(1),thi(1)

             rho = 0.d0
             do n=1, nspecies
                rho = rho + U(i,j,k,iry1+n-1)
             end do
             U(i,j,k,irho) = rho

             !
             ! Enforce nonnegative species
             !
             rhoy_under = 0.d0
             do n = 1, nspecies
                iryn = iry1+n-1
                if (U(i,j,k,iryn) .lt. 0.d0) then
                   rhoy_under = rhoy_under + U(i,j,k,iryn)
                   U(i,j,k,iryn) = 0.d0
                end if
             end do

             if (rhoy_under .lt. rho*eps) then
                !
                ! Find the dominant species.
                !
                idom = 1
                rhoy_dom = U(i,j,k,iry1)
                do n = 2, nspecies
                   iryn = iry1+n-1
                   if (U(i,j,k,iryn) .gt. rhoy_dom) then
                      idom = n
                      rhoy_dom = U(i,j,k,iryn)
                   end if
                end do
                !
                ! Take enough from the dominant species to fill the negative one.
                !
                !iryn = iry1+idom-1
                !U(i,j,k,iryn) = U(i,j,k,iryn) + rhoy_under
                !if (U(i,j,k,iryn) .lt. 0.d0) then
                !   print *,'Just made dominant species',idom, &
                !        'negative', U(i,j,k,iryn)/rho, 'at ',i,j,k 
                !   call bl_error("Error:: variables :: reset_rho_3d")
                !end if
             end if

          end do
       end do
    end do

  end subroutine reset_rho_3d

end module variables_module
