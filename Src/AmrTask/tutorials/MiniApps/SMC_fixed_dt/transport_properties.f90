module transport_properties

  use iso_c_binding
  use chemistry_module
  use variables_module

  implicit none

  double precision, parameter :: Pr = 0.72d0   ! Prandtl number
  double precision, parameter :: Sc = 0.72d0   ! Schmidt number
  double precision, parameter :: C_S = 1.458d-5  ! constant in Sutherland's law
  double precision, parameter :: T_S = 110.4d0   ! Sutherland temperature

  private

  public get_trans_prop_3d

contains

  subroutine get_trans_prop_3d(tlo,thi,lo,hi,q,mu,xi,lam,Ddiag,ng) bind(c,name='get_trans_prop_3d')
    integer, intent(in) :: lo(3), hi(3), ng, tlo(3), thi(3)
    double precision,intent(in   )::    q(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,nprim)
    double precision,intent(inout)::   mu(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(inout)::   xi(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(inout)::  lam(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng)
    double precision,intent(inout)::Ddiag(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng,nspecies)

    integer :: i, j, k, n, iwrk
    double precision :: rwrk, cp, Wbar, Wbarinv, Yt(nspecies)
    double precision :: Prinv, Scinv

    Prinv = 1.d0/Pr
    Scinv = 1.d0/Sc

    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
          
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

  end subroutine get_trans_prop_3d

end module transport_properties
