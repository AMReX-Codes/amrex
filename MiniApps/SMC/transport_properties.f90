module transport_properties

  use chemistry_module
  use multifab_module
  use variables_module

  implicit none

  double precision, parameter :: Pr = 0.72d0   ! Prandtl number
  double precision, parameter :: Sc = 0.72d0   ! Schmidt number
  double precision, parameter :: C_S = 1.458d-5  ! constant in Sutherland's law
  double precision, parameter :: T_S = 110.4d0   ! Sutherland temperature

  private

  public get_transport_properties

contains

  subroutine get_transport_properties(Q, mu, xi, lam, Ddiag)
    type(multifab), intent(in   ) :: Q
    type(multifab), intent(inout) :: mu, xi, lam, Ddiag
 
    integer :: ng, n, lo(3), hi(3), wlo(3), whi(3)
    double precision, pointer, dimension(:,:,:,:) :: qp, mup, xip, lamp, dp

    ng = nghost(Q)

    do n=1,nfabs(Q)
       
       qp => dataptr(Q,n)
       mup => dataptr(mu,n)
       xip => dataptr(xi,n)
       lamp => dataptr(lam,n)
       dp => dataptr(Ddiag,n)

       lo = lwb(get_box(Q,n))
       hi = upb(get_box(Q,n))

       wlo = lo - ng
       whi = hi + ng

       call get_trans_prop_3d(lo,hi,ng,qp,mup,xip,lamp,dp,wlo,whi)

    end do

  end subroutine get_transport_properties

  subroutine get_trans_prop_3d(lo,hi,ng,q,mu,xi,lam,Ddiag,wlo,whi)
    integer, intent(in) :: lo(3), hi(3), ng, wlo(3), whi(3)
    double precision,intent(in )::    q(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,nprim)
    double precision,intent(out)::   mu(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(out)::   xi(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(out)::  lam(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(out)::Ddiag(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,nspecies)

    integer :: i, j, k, n, iwrk
    double precision :: rwrk, cp, Wbar, Wbarinv, Yt(nspecies)
    double precision :: Prinv, Scinv

    !$omp parallel private(i,j,k,n,iwrk,rwrk,cp,Wbar,Wbarinv,Yt,Prinv,Scinv)
       
    Prinv = 1.d0/Pr
    Scinv = 1.d0/Sc

    !$omp do collapse(2)
    do k=wlo(3),whi(3)
       do j=wlo(2),whi(2)
          do i=wlo(1),whi(1)
          
             mu(i,j,k) = C_S * sqrt(q(i,j,k,qtemp)) * q(i,j,k,qtemp) &
                  / (q(i,j,k,qtemp) + T_S)

             xi(i,j,k) = 0.d0

             do n=1,nspecies
                Yt(n) = q(i,j,k,qy1+n-1)
             end do

             call ckmmwy(Yt, iwrk, rwrk, Wbar)
             Wbarinv = 1.d0 / Wbar

             call ckcpbs(q(i,j,k,qtemp), Yt, iwrk, rwrk, cp)

             lam(i,j,k) = mu(i,j,k) * cp * Prinv

             do n=1,nspecies
                Ddiag(i,j,k,n) = mu(i,j,k) * molecular_weight(n) * Wbarinv * Scinv
             end do

          end do
       end do
    end do
    !$omp end do
       
    !$omp end parallel

  end subroutine get_trans_prop_3d

end module transport_properties
