module init_data_module

contains

  subroutine init_data_3d(tlo,thi,lo,hi,cons,ng,dx,phlo,phhi) bind(c,name='init_data_3d')
    
    use variables_module, only : irho, imx,imy,imz,iene,iry1,ncons, iH2, iO2, iN2
    use chemistry_module, only : nspecies, patm
    
    implicit none
    
    integer,          intent(in   ) :: tlo(3),thi(3),lo(3),hi(3),ng
    double precision, intent(in   ) :: dx(3),phlo(3),phhi(3)
    double precision, intent(inout) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)
    
    double precision, parameter :: Pi = 3.141592653589793238462643383279502884197d0
    double precision, parameter :: rfire = 0.01d0
    double precision, parameter :: Lperiodic = 0.1d0
    
    integer          :: i,j,k,n
    double precision :: x, y, z, r
    
    double precision Xt(nspecies), Yt(nspecies)
    double precision rhot,u1t,u2t,u3t,Tt,et,Pt
    integer :: iwrk
    double precision :: rwrk, expfac, kx, ky, kz
  
    kx = 2.d0*Pi/(phhi(1) - phlo(1))
    ky = 2.d0*Pi/(phhi(2) - phlo(2))
    kz = 2.d0*Pi/(phhi(3) - phlo(3))
    
    do k=tlo(3),thi(3)
       z = mod(dx(3)*(k+0.5d0), Lperiodic) - 0.5d0*Lperiodic 
       do j=tlo(2),thi(2)
          y = mod(dx(2)*(j+0.5d0), Lperiodic) - 0.5d0*Lperiodic
          do i=tlo(1),thi(1)
             x = mod(dx(1)*(i+0.5d0), Lperiodic) - 0.5d0*Lperiodic
             
             r = sqrt(x**2+y**2+z**2)
             
             Pt = Patm
             Tt = 300.0d0
             
             Xt = 0.0d0
             Xt(iH2) = 0.10d0
             Xt(iO2) = 0.25d0
             
             expfac = exp(-(r / rfire)**2)
             Pt      = Pt      + 0.1d0*patm * expfac
             Tt      = Tt      + 1100.0d0   * expfac
             Xt(iH2) = Xt(iH2) + 0.025d0    * expfac
             Xt(iO2) = Xt(iO2) - 0.050d0    * expfac
             
             Xt(iN2) = 1.0d0 - Xt(iH2) - Xt(iO2)
             
             u1t =  sin(kx*x)*cos(ky*y)*cos(kz*z) * 300.d0
             u2t = -cos(kx*x)*sin(ky*y)*cos(kz*z) * 300.d0
             u3t = 0.d0
             
             CALL CKXTY (Xt, IWRK, RWRK, Yt)
             CALL CKRHOY(Pt,Tt,Yt,IWRK,RWRK,rhot)
             call CKUBMS(Tt,Yt,IWRK,RWRK,et)
             
             cons(i,j,k,irho) = rhot
             cons(i,j,k,imx)  = rhot*u1t
             cons(i,j,k,imy)  = rhot*u2t
             cons(i,j,k,imz)  = rhot*u3t
             cons(i,j,k,iene) = rhot*(et + 0.5d0*(u1t**2 + u2t**2 + u3t**2))
             
             do n=1,nspecies
                cons(i,j,k,iry1-1+n) = Yt(n)*rhot
             end do
             
          enddo
       enddo
    enddo
    
  end subroutine init_data_3d
  
end module init_data_module
