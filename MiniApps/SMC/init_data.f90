module init_data_module

  use multifab_module

  implicit none

  double precision, parameter :: Pi = 3.141592653589793238462643383279502884197d0

  private

  public :: init_data

contains

  subroutine init_data(data,dx,plo,phi)

    type(multifab),   intent(inout) :: data
    double precision, intent(in   ) :: dx(data%dim)
    double precision, intent(in   ) :: plo(data%dim), phi(data%dim)
    
    integer                   :: lo(data%dim), hi(data%dim), ng, i
    integer                   :: tlo(data%dim), thi(data%dim)
    type(mfiter)              :: mfi
    type(box)                 :: tbx
    double precision, pointer :: dp(:,:,:,:)

    ng = data%ng

    !$omp parallel private(lo,hi,i,tlo,thi,mfi,tbx,dp)

    call mfiter_build(mfi,data,tiling=.true.)

    do while (next_tile(mfi,i))
       
       tbx = get_tilebox(mfi)
       tlo = lwb(tbx)
       thi = upb(tbx)
       
       dp => dataptr(data,i)
       lo = lwb(get_box(data,i))
       hi = upb(get_box(data,i))
       
       select case(data%dim)
       case (2)
          call bl_error('We only support 3-D')
       case (3)
          call init_data_3d(tlo,thi,lo,hi,ng,dx,dp,plo,phi)
       end select
    end do
    !$omp end parallel

  end subroutine init_data

  subroutine init_data_3d(lo,hi,clo,chi,ng,dx,cons,phlo,phhi)

    use variables_module, only : irho, imx,imy,imz,iene,iry1,ncons, iH2, iO2, iN2
    use chemistry_module, only : nspecies, patm
    use probin_module,    only : rfire

    integer,          intent(in   ) :: lo(3),hi(3),clo(3),chi(3),ng
    double precision, intent(in   ) :: dx(3),phlo(3),phhi(3)
    double precision, intent(inout) :: cons(-ng+clo(1):chi(1)+ng,-ng+clo(2):chi(2)+ng,-ng+clo(3):chi(3)+ng,ncons)

    integer          :: i,j,k,n
    double precision :: x, y, z, r

    double precision Xt(nspecies), Yt(nspecies)
    double precision rhot,u1t,u2t,u3t,Tt,et,Pt
    integer :: iwrk
    double precision :: rwrk, expfac, kx, ky, kz

    kx = 2.d0*Pi/(phhi(1) - phlo(1))
    ky = 2.d0*Pi/(phhi(2) - phlo(2))
    kz = 2.d0*Pi/(phhi(3) - phlo(3))

    do k=lo(3),hi(3)
       z = phlo(3) + dx(3)*(k + 0.5d0)
       do j=lo(2),hi(2)
          y = phlo(2) + dx(2)*(j + 0.5d0)
          do i=lo(1),hi(1)
             x = phlo(1) + dx(1)*(i + 0.5d0)

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
