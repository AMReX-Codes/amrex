module variables_module

  use chemistry_module, only : nspecies
  use multifab_module

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
  subroutine init_variables()

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
    
  end subroutine init_variables

  !
  ! Convert conserved variables U to primitive variables Q
  !
  subroutine ctoprim(U, Q, ng)
    type(multifab), intent(in   ) :: U
    type(multifab), intent(inout) :: Q
    integer, optional, intent(in) :: ng

    integer :: ngu, ngq, ngto
    integer :: n, lo(3), hi(3)
    double precision, pointer, dimension(:,:,:,:) :: up, qp

    ngu = nghost(U)
    ngq = nghost(Q)

    if (present(ng)) then
       ngto = ng
    else
       ngto = min(ngu, ngq)
    end if

    do n=1,nfabs(Q)
       up => dataptr(U,n)
       qp => dataptr(Q,n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       call ctoprim_3d(lo,hi,up,qp,ngu,ngq,ngto)
    end do

  end subroutine ctoprim

  subroutine ctoprim_3d(lo, hi, u, q, ngu, ngq, ngto)
    integer, intent(in) :: lo(3), hi(3), ngu, ngq, ngto
    double precision, intent(in ) :: u(lo(1)-ngu:hi(1)+ngu,lo(2)-ngu:hi(2)+ngu,lo(3)-ngu:hi(3)+ngu,ncons)
    double precision, intent(out) :: q(lo(1)-ngq:hi(1)+ngq,lo(2)-ngq:hi(2)+ngq,lo(3)-ngq:hi(3)+ngq,nprim)
    
    integer :: i, j, k, n, iwrk, ierr
    double precision :: rho, rhoinv, rwrk, X(nspecies), Y(nspecies), h(nspecies), ei, Tt, Pt
    integer :: llo(3), lhi(3)

    do i=1,3
       llo(i) = lo(i)-ngto
       lhi(i) = hi(i)+ngto
    end do

    !$omp parallel private(i, j, k, n, iwrk, rho, rhoinv, rwrk) &
    !$omp private(X, Y, h, ei, Tt, Pt, ierr)
    !$omp do collapse(3)
    do k = llo(3),lhi(3)
       do j = llo(2),lhi(2)
          do i = llo(1),lhi(1)
             
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
    !$omp end do
    !$omp end parallel

  end subroutine ctoprim_3d

  !
  ! Compute total density
  !
  subroutine reset_density(U)
    type(multifab), intent(inout) :: U
    
    integer :: ng
    integer :: n, lo(3), hi(3)
    double precision, pointer, dimension(:,:,:,:) :: up

    ng = nghost(U)
    
    do n=1,nfabs(U)
       up => dataptr(U, n)
       
       lo = lwb(get_box(U,n))
       hi = upb(get_box(U,n))

       call reset_rho_3d(lo,hi,ng,up)
    end do

  end subroutine reset_density

  subroutine reset_rho_3d(lo, hi, ng, u)
    integer, intent(in) :: lo(3), hi(3), ng
    double precision, intent(inout) :: u(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,ncons)

    integer :: i, j, k, n, iryn
    double precision :: rho
    double precision, parameter :: eps = -1.0d-16
    integer          :: idom
    double precision :: rhoy_dom, rhoy_under

    !$omp parallel private(i,j,k,n,iryn,rho) &
    !$omp private(idom, rhoy_dom, rhoy_under)
    !$omp do collapse (3)
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

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
                iryn = iry1+idom-1
                U(i,j,k,iryn) = U(i,j,k,iryn) + rhoy_under
                if (U(i,j,k,iryn) .lt. 0.d0) then
                   print *,'Just made dominant species',idom, &
                        'negative', U(i,j,k,iryn)/rho, 'at ',i,j,k 
                   call bl_error("Error:: variables :: reset_rho_3d")
                end if
             end if

          end do
       end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine reset_rho_3d

end module variables_module
