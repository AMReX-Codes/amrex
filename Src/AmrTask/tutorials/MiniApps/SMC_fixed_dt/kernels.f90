module kernels_module
  use chemistry_module, only : nspecies, molecular_weight, Ru
  use derivative_stencil_module, only : stencil_ng, M8, M8T, D8 
  use variables_module
  use iso_c_binding

  implicit none

  private

  public :: hypterm_3d, narrow_diffterm_3d, chemterm_3d, comp_courno_3d

contains

  subroutine hypterm_3d (lo,hi,dx,cons,clo,chi,q,qlo,qhi,rhs,rlo,rhi) bind(c,name='hypterm_3d')

    integer,         intent(in)   :: lo(3),hi(3),clo(3),chi(3),qlo(3),qhi(3),rlo(3),rhi(3)
    double precision,intent(in)   :: dx(3)
    double precision,intent(in)   ::cons(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),ncons)
    double precision,intent(in)   ::   q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)
    double precision,intent(inout):: rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3),ncons)

    integer          :: i,j,k,n
    double precision :: dxinv(3)

    double precision :: tmpx(lo(1)-4:hi(1)+4)
    double precision :: tmpy(lo(1)  :hi(1)  ,lo(2)-4:hi(2)+4)
    double precision :: tmpz(lo(1)  :hi(1)  ,lo(2)  :hi(2)  ,lo(3)-4:hi(3)+4)

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    ! ------- BEGIN x-direction -------

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=lo(1),hi(1)
             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                ( D8(1)*(cons(i+1,j,k,imx)-cons(i-1,j,k,imx)) &
                + D8(2)*(cons(i+2,j,k,imx)-cons(i-2,j,k,imx)) &
                + D8(3)*(cons(i+3,j,k,imx)-cons(i-3,j,k,imx)) &
                + D8(4)*(cons(i+4,j,k,imx)-cons(i-4,j,k,imx)) )
          end do

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = cons(i,j,k,imx)*q(i,j,k,qu)+q(i,j,k,qpres)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                ( D8(1)*(tmpx(i+1)-tmpx(i-1)) &
                + D8(2)*(tmpx(i+2)-tmpx(i-2)) &
                + D8(3)*(tmpx(i+3)-tmpx(i-3)) &
                + D8(4)*(tmpx(i+4)-tmpx(i-4)) )
          end do

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = cons(i,j,k,imy)*q(i,j,k,qu)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                ( D8(1)*(tmpx(i+1)-tmpx(i-1)) &
                + D8(2)*(tmpx(i+2)-tmpx(i-2)) &
                + D8(3)*(tmpx(i+3)-tmpx(i-3)) &
                + D8(4)*(tmpx(i+4)-tmpx(i-4)) )
          end do

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = cons(i,j,k,imz)*q(i,j,k,qu)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                ( D8(1)*(tmpx(i+1)-tmpx(i-1)) &
                + D8(2)*(tmpx(i+2)-tmpx(i-2)) &
                + D8(3)*(tmpx(i+3)-tmpx(i-3)) &
                + D8(4)*(tmpx(i+4)-tmpx(i-4)) )
          end do

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = (cons(i,j,k,iene)+q(i,j,k,qpres))*q(i,j,k,qu)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                ( D8(1)*(tmpx(i+1)-tmpx(i-1)) &
                + D8(2)*(tmpx(i+2)-tmpx(i-2)) &
                + D8(3)*(tmpx(i+3)-tmpx(i-3)) &
                + D8(4)*(tmpx(i+4)-tmpx(i-4)) )
          end do

       end do
    end do

    do n = iry1, iry1+nspecies-1
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
    
             do i=lo(1)-4,hi(1)+4
                tmpx(i) = cons(i,j,k,n)*q(i,j,k,qu)
             end do
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                   ( D8(1)*(tmpx(i+1)-tmpx(i-1)) &
                   + D8(2)*(tmpx(i+2)-tmpx(i-2)) &
                   + D8(3)*(tmpx(i+3)-tmpx(i-3)) &
                   + D8(4)*(tmpx(i+4)-tmpx(i-4)) )
             end do

          enddo
       enddo
    enddo

    ! ------- END x-direction -------

    ! ------- BEGIN y-direction -------

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(2) * &
                ( D8(1)*(cons(i,j+1,k,imy)-cons(i,j-1,k,imy)) &
                + D8(2)*(cons(i,j+2,k,imy)-cons(i,j-2,k,imy)) &
                + D8(3)*(cons(i,j+3,k,imy)-cons(i,j-3,k,imy)) &
                + D8(4)*(cons(i,j+4,k,imy)-cons(i,j-4,k,imy)) )
          enddo
       enddo

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = cons(i,j,k,imx)*q(i,j,k,qv)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(2) * &
                ( D8(1)*(tmpy(i,j+1)-tmpy(i,j-1)) &
                + D8(2)*(tmpy(i,j+2)-tmpy(i,j-2)) &
                + D8(3)*(tmpy(i,j+3)-tmpy(i,j-3)) &
                + D8(4)*(tmpy(i,j+4)-tmpy(i,j-4)) )
          enddo
       enddo

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = cons(i,j,k,imy)*q(i,j,k,qv)+q(i,j,k,qpres)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(2) * &
                ( D8(1)*(tmpy(i,j+1)-tmpy(i,j-1)) &
                + D8(2)*(tmpy(i,j+2)-tmpy(i,j-2)) &
                + D8(3)*(tmpy(i,j+3)-tmpy(i,j-3)) &
                + D8(4)*(tmpy(i,j+4)-tmpy(i,j-4)) )
          enddo
       enddo

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = cons(i,j,k,imz)*q(i,j,k,qv)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(2) * &
                ( D8(1)*(tmpy(i,j+1)-tmpy(i,j-1)) &
                + D8(2)*(tmpy(i,j+2)-tmpy(i,j-2)) &
                + D8(3)*(tmpy(i,j+3)-tmpy(i,j-3)) &
                + D8(4)*(tmpy(i,j+4)-tmpy(i,j-4)) )
          enddo
       enddo

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = (cons(i,j,k,iene)+q(i,j,k,qpres))*q(i,j,k,qv)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(2) * &
                ( D8(1)*(tmpy(i,j+1)-tmpy(i,j-1)) &
                + D8(2)*(tmpy(i,j+2)-tmpy(i,j-2)) &
                + D8(3)*(tmpy(i,j+3)-tmpy(i,j-3)) &
                + D8(4)*(tmpy(i,j+4)-tmpy(i,j-4)) )
          enddo
       enddo
    enddo

    do n = iry1, iry1+nspecies-1
       do k=lo(3),hi(3)
          do j=lo(2)-4,hi(2)+4
             do i=lo(1),hi(1)
                tmpy(i,j) = cons(i,j,k,n)*q(i,j,k,qv)
             end do
          end do
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                   ( D8(1)*(tmpy(i,j+1)-tmpy(i,j-1)) &
                   + D8(2)*(tmpy(i,j+2)-tmpy(i,j-2)) &
                   + D8(3)*(tmpy(i,j+3)-tmpy(i,j-3)) &
                   + D8(4)*(tmpy(i,j+4)-tmpy(i,j-4)) )
             end do
          enddo
       enddo
    enddo

    ! ------- END y-direction -------

    ! ------- BEGIN z-direction -------

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(3) * &
                ( D8(1)*(cons(i,j,k+1,imz)-cons(i,j,k-1,imz)) &
                + D8(2)*(cons(i,j,k+2,imz)-cons(i,j,k-2,imz)) &
                + D8(3)*(cons(i,j,k+3,imz)-cons(i,j,k-3,imz)) &
                + D8(4)*(cons(i,j,k+4,imz)-cons(i,j,k-4,imz)) )
          end do
       end do
    end do

    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = cons(i,j,k,imx) * q(i,j,k,qw)
          end do
       end do
    end do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(3) * &
                ( D8(1)*(tmpz(i,j,k+1)-tmpz(i,j,k-1)) &
                + D8(2)*(tmpz(i,j,k+2)-tmpz(i,j,k-2)) &
                + D8(3)*(tmpz(i,j,k+3)-tmpz(i,j,k-3)) &
                + D8(4)*(tmpz(i,j,k+4)-tmpz(i,j,k-4)) )
          end do
       end do
    end do

    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = cons(i,j,k,imy) * q(i,j,k,qw)
          end do
       end do
    end do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(3) * &
                ( D8(1)*(tmpz(i,j,k+1)-tmpz(i,j,k-1)) &
                + D8(2)*(tmpz(i,j,k+2)-tmpz(i,j,k-2)) &
                + D8(3)*(tmpz(i,j,k+3)-tmpz(i,j,k-3)) &
                + D8(4)*(tmpz(i,j,k+4)-tmpz(i,j,k-4)) )
          end do
       end do
    end do

    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = cons(i,j,k,imz)*q(i,j,k,qw) + q(i,j,k,qpres)
          end do
       end do
    end do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(3) * &
                ( D8(1)*(tmpz(i,j,k+1)-tmpz(i,j,k-1)) &
                + D8(2)*(tmpz(i,j,k+2)-tmpz(i,j,k-2)) &
                + D8(3)*(tmpz(i,j,k+3)-tmpz(i,j,k-3)) &
                + D8(4)*(tmpz(i,j,k+4)-tmpz(i,j,k-4)) )
          end do
       end do
    end do

    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = (cons(i,j,k,iene)+q(i,j,k,qpres))*q(i,j,k,qw)
          end do
       end do
    end do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(3) * &
                ( D8(1)*(tmpz(i,j,k+1)-tmpz(i,j,k-1)) &
                + D8(2)*(tmpz(i,j,k+2)-tmpz(i,j,k-2)) &
                + D8(3)*(tmpz(i,j,k+3)-tmpz(i,j,k-3)) &
                + D8(4)*(tmpz(i,j,k+4)-tmpz(i,j,k-4)) )
          end do
       end do
    end do

    do n = iry1, iry1+nspecies-1
       do k=lo(3)-4,hi(3)+4
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                tmpz(i,j,k) = cons(i,j,k,n)*q(i,j,k,qw)
             end do
          end do
       end do
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                   ( D8(1)*(tmpz(i,j,k+1)-tmpz(i,j,k-1)) &
                   + D8(2)*(tmpz(i,j,k+2)-tmpz(i,j,k-2)) &
                   + D8(3)*(tmpz(i,j,k+3)-tmpz(i,j,k-3)) &
                   + D8(4)*(tmpz(i,j,k+4)-tmpz(i,j,k-4)) )
             end do
          enddo
       enddo
    enddo

  end subroutine hypterm_3d


  subroutine narrow_diffterm_3d (lo,hi,dx,q,qlo,qhi,rhs_g,glo,ghi,mu,xi,lam,dxy) &
       bind(c,name='narrow_diffterm_3d')

    integer,         intent(in):: lo(3),hi(3),qlo(3),qhi(3),glo(3),ghi(3)
    double precision,intent(in):: dx(3)
    double precision,intent(in)   ::  q  (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)
    double precision,intent(in)   ::  mu (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in)   ::  xi (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in)   ::  lam(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in)   ::  dxy(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nspecies)
    double precision,intent(inout)::rhs_g(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),ncons)

    integer :: i, dlo(3), dhi(3)
    double precision :: dxinv(3), dx2inv(3)

    double precision :: rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons)

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
       dx2inv(i) = dxinv(i)**2
    end do

    dlo = lo - stencil_ng
    dhi = hi + stencil_ng

    rhs = 0.d0

    call diffterm_1(q,qlo,qhi,rhs,lo,hi,mu,xi, lo,hi,dlo,dhi,dxinv)

    call diffterm_2(q,qlo,qhi,rhs,lo,hi,mu,xi,lam,dxy, lo,hi,dlo,dhi,dxinv,dx2inv)

    rhs_g(     lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = &
         rhs_g(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) &
         + rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)

  end subroutine narrow_diffterm_3d

  subroutine diffterm_1(q,qlo,qhi,rhs,rlo,rhi,mu,xi, lo,hi,dlo,dhi,dxinv)
    integer,         intent(in):: lo(3),hi(3),dlo(3),dhi(3)
    integer,         intent(in):: qlo(3),qhi(3),rlo(3),rhi(3)
    double precision,intent(in):: dxinv(3)
    double precision,intent(in)   :: q (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)
    double precision,intent(in)   :: mu(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in)   :: xi(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(inout)::rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3),ncons)

    !
    ! local variables
    !
    integer          :: i,j,k
    double precision, dimension(lo(1):hi(1)) :: tauxx,tauyy,tauzz,divu

    double precision :: ux  ( lo(1): hi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    double precision :: vx  ( lo(1): hi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    double precision :: wx  ( lo(1): hi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    double precision :: uy  (dlo(1):dhi(1), lo(2): hi(2),dlo(3):dhi(3))
    double precision :: vy  (dlo(1):dhi(1), lo(2): hi(2),dlo(3):dhi(3))
    double precision :: wy  (dlo(1):dhi(1), lo(2): hi(2),dlo(3):dhi(3))
    double precision :: uz  (dlo(1):dhi(1),dlo(2):dhi(2), lo(3): hi(3))
    double precision :: vz  (dlo(1):dhi(1),dlo(2):dhi(2), lo(3): hi(3))
    double precision :: wz  (dlo(1):dhi(1),dlo(2):dhi(2), lo(3): hi(3))
    double precision :: vsm (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    double precision :: tmpx(dlo(1):dhi(1))
    double precision :: tmpy( lo(1): hi(1),dlo(2):dhi(2))
    double precision :: tmpz( lo(1): hi(1), lo(2): hi(2),dlo(3):dhi(3))

    do k=dlo(3),dhi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             vsm(i,j,k) = xi(i,j,k) -  TwoThirds*mu(i,j,k)
          enddo
       enddo
    enddo

    do k=dlo(3),dhi(3)
       do j=dlo(2),dhi(2)
          do i=lo(1),hi(1)
             ux(i,j,k) = dxinv(1) * &
                ( D8(1)*(q(i+1,j,k,qu)-q(i-1,j,k,qu)) &
                + D8(2)*(q(i+2,j,k,qu)-q(i-2,j,k,qu)) &
                + D8(3)*(q(i+3,j,k,qu)-q(i-3,j,k,qu)) &
                + D8(4)*(q(i+4,j,k,qu)-q(i-4,j,k,qu)) )
             vx(i,j,k) = dxinv(1) * &
                ( D8(1)*(q(i+1,j,k,qv)-q(i-1,j,k,qv)) &
                + D8(2)*(q(i+2,j,k,qv)-q(i-2,j,k,qv)) &
                + D8(3)*(q(i+3,j,k,qv)-q(i-3,j,k,qv)) &
                + D8(4)*(q(i+4,j,k,qv)-q(i-4,j,k,qv)) )
             wx(i,j,k) = dxinv(1) * &
                ( D8(1)*(q(i+1,j,k,qw)-q(i-1,j,k,qw)) &
                + D8(2)*(q(i+2,j,k,qw)-q(i-2,j,k,qw)) &
                + D8(3)*(q(i+3,j,k,qw)-q(i-3,j,k,qw)) &
                + D8(4)*(q(i+4,j,k,qw)-q(i-4,j,k,qw)) )
          enddo
       enddo
    enddo

    do k=dlo(3),dhi(3)
       do j=lo(2),hi(2)   

          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2) * &
                ( D8(1)*(q(i,j+1,k,qu)-q(i,j-1,k,qu)) &
                + D8(2)*(q(i,j+2,k,qu)-q(i,j-2,k,qu)) &
                + D8(3)*(q(i,j+3,k,qu)-q(i,j-3,k,qu)) &
                + D8(4)*(q(i,j+4,k,qu)-q(i,j-4,k,qu)) )
             vy(i,j,k) = dxinv(2) * &
                ( D8(1)*(q(i,j+1,k,qv)-q(i,j-1,k,qv)) &
                + D8(2)*(q(i,j+2,k,qv)-q(i,j-2,k,qv)) &
                + D8(3)*(q(i,j+3,k,qv)-q(i,j-3,k,qv)) &
                + D8(4)*(q(i,j+4,k,qv)-q(i,j-4,k,qv)) )
             wy(i,j,k) = dxinv(2) * &
                ( D8(1)*(q(i,j+1,k,qw)-q(i,j-1,k,qw)) &
                + D8(2)*(q(i,j+2,k,qw)-q(i,j-2,k,qw)) &
                + D8(3)*(q(i,j+3,k,qw)-q(i,j-3,k,qw)) &
                + D8(4)*(q(i,j+4,k,qw)-q(i,j-4,k,qw)) )
          enddo

       enddo
    enddo

    do k=lo(3),hi(3)
       do j=dlo(2),dhi(2)    
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3) * &
                ( D8(1)*(q(i,j,k+1,qu)-q(i,j,k-1,qu)) &
                + D8(2)*(q(i,j,k+2,qu)-q(i,j,k-2,qu)) &
                + D8(3)*(q(i,j,k+3,qu)-q(i,j,k-3,qu)) &
                + D8(4)*(q(i,j,k+4,qu)-q(i,j,k-4,qu)) )
             vz(i,j,k) = dxinv(3) * &
                ( D8(1)*(q(i,j,k+1,qv)-q(i,j,k-1,qv)) &
                + D8(2)*(q(i,j,k+2,qv)-q(i,j,k-2,qv)) &
                + D8(3)*(q(i,j,k+3,qv)-q(i,j,k-3,qv)) &
                + D8(4)*(q(i,j,k+4,qv)-q(i,j,k-4,qv)) )
             wz(i,j,k) = dxinv(3) * &
                ( D8(1)*(q(i,j,k+1,qw)-q(i,j,k-1,qw)) &
                + D8(2)*(q(i,j,k+2,qw)-q(i,j,k-2,qw)) &
                + D8(3)*(q(i,j,k+3,qw)-q(i,j,k-3,qw)) &
                + D8(4)*(q(i,j,k+4,qw)-q(i,j,k-4,qw)) )
          end do      
       end do
    end do

    !----- mx -----

    !----- mx : d()/dx -----
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = vsm(i,j,k)*(vy(i,j,k)+wz(i,j,k))
          end do

          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(1) * &
                ( D8(1)*(tmpx(i+1)-tmpx(i-1)) &
                + D8(2)*(tmpx(i+2)-tmpx(i-2)) &
                + D8(3)*(tmpx(i+3)-tmpx(i-3)) &
                + D8(4)*(tmpx(i+4)-tmpx(i-4)) )
          end do

       end do
    end do

    !----- mx : d()/dy -----
    do k=lo(3),hi(3)

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = mu(i,j,k)*vx(i,j,k)
          end do
       end do

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(2) * &
                ( D8(1)*(tmpy(i,j+1)-tmpy(i,j-1)) &
                + D8(2)*(tmpy(i,j+2)-tmpy(i,j-2)) &
                + D8(3)*(tmpy(i,j+3)-tmpy(i,j-3)) &
                + D8(4)*(tmpy(i,j+4)-tmpy(i,j-4)) )
          end do
       end do

    end do

    !----- mx : d()/dz -----
    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = mu(i,j,k)*wx(i,j,k)
          end do
       end do
    end do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(3) * &
                ( D8(1)*(tmpz(i,j,k+1)-tmpz(i,j,k-1)) &
                + D8(2)*(tmpz(i,j,k+2)-tmpz(i,j,k-2)) &
                + D8(3)*(tmpz(i,j,k+3)-tmpz(i,j,k-3)) &
                + D8(4)*(tmpz(i,j,k+4)-tmpz(i,j,k-4)) )
          end do
       end do
    end do

    !----- my -----

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = mu(i,j,k)*uy(i,j,k)
          end do

          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(1) * &
                ( D8(1)*(tmpx(i+1)-tmpx(i-1)) &
                + D8(2)*(tmpx(i+2)-tmpx(i-2)) &
                + D8(3)*(tmpx(i+3)-tmpx(i-3)) &
                + D8(4)*(tmpx(i+4)-tmpx(i-4)) )
          end do

       end do
    end do

    do k=lo(3),hi(3)

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = vsm(i,j,k)*(ux(i,j,k)+wz(i,j,k))
          end do
       end do

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(2) * &
                ( D8(1)*(tmpy(i,j+1)-tmpy(i,j-1)) &
                + D8(2)*(tmpy(i,j+2)-tmpy(i,j-2)) &
                + D8(3)*(tmpy(i,j+3)-tmpy(i,j-3)) &
                + D8(4)*(tmpy(i,j+4)-tmpy(i,j-4)) )
          end do
       end do

    end do

    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = mu(i,j,k)*wy(i,j,k)
          end do
       end do
    end do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(3) * &
                ( D8(1)*(tmpz(i,j,k+1)-tmpz(i,j,k-1)) &
                + D8(2)*(tmpz(i,j,k+2)-tmpz(i,j,k-2)) &
                + D8(3)*(tmpz(i,j,k+3)-tmpz(i,j,k-3)) &
                + D8(4)*(tmpz(i,j,k+4)-tmpz(i,j,k-4)) )
          end do
       end do
    end do

    !----- mz -----

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = mu(i,j,k)*uz(i,j,k)
          end do

          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(1) * &
                ( D8(1)*(tmpx(i+1)-tmpx(i-1)) &
                + D8(2)*(tmpx(i+2)-tmpx(i-2)) &
                + D8(3)*(tmpx(i+3)-tmpx(i-3)) &
                + D8(4)*(tmpx(i+4)-tmpx(i-4)) )
          end do

       end do
    end do

    do k=lo(3),hi(3)

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = mu(i,j,k)*vz(i,j,k)
          end do
       end do

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(2) * &
                ( D8(1)*(tmpy(i,j+1)-tmpy(i,j-1)) &
                + D8(2)*(tmpy(i,j+2)-tmpy(i,j-2)) &
                + D8(3)*(tmpy(i,j+3)-tmpy(i,j-3)) &
                + D8(4)*(tmpy(i,j+4)-tmpy(i,j-4)) )
          end do
       end do

    end do

    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = vsm(i,j,k)*(ux(i,j,k)+vy(i,j,k))
          end do
       end do
    end do
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(3) * &
                ( D8(1)*(tmpz(i,j,k+1)-tmpz(i,j,k-1)) &
                + D8(2)*(tmpz(i,j,k+2)-tmpz(i,j,k-2)) &
                + D8(3)*(tmpz(i,j,k+3)-tmpz(i,j,k-3)) &
                + D8(4)*(tmpz(i,j,k+4)-tmpz(i,j,k-4)) )
          end do
       end do
    end do

    !----- energy -----

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             divu(i) = (ux(i,j,k)+vy(i,j,k)+wz(i,j,k))*vsm(i,j,k)
             tauxx(i) = 2.d0*mu(i,j,k)*ux(i,j,k) + divu(i)
             tauyy(i) = 2.d0*mu(i,j,k)*vy(i,j,k) + divu(i)
             tauzz(i) = 2.d0*mu(i,j,k)*wz(i,j,k) + divu(i)
             
             ! change in internal energy
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + &
                  tauxx(i)*ux(i,j,k) + tauyy(i)*vy(i,j,k) + tauzz(i)*wz(i,j,k) &
                  + mu(i,j,k)*((uy(i,j,k)+vx(i,j,k))**2 &
                  &          + (wx(i,j,k)+uz(i,j,k))**2 &
                  &          + (vz(i,j,k)+wy(i,j,k))**2 )

          end do
       end do
    end do

  end subroutine diffterm_1

  subroutine diffterm_2(q,qlo,qhi,rhs,rlo,rhi,mu,xi,lam,dxy,lo,hi,dlo,dhi,dxinv,dx2inv)
    integer,         intent(in):: lo(3),hi(3),dlo(3),dhi(3)
    integer,         intent(in):: qlo(3),qhi(3),rlo(3),rhi(3)
    double precision,intent(in):: dxinv(3),dx2inv(3)
    double precision,intent(in)   ::  q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)
    double precision,intent(in)   :: mu(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in)   :: xi(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in)   ::lam(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in)   ::dxy(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nspecies)
    double precision,intent(inout)::rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3),ncons)

    !
    ! local variables
    !
    integer          :: i,j,k,n, qxn, qyn, qhn, iryn
    double precision :: mmtmp(8,lo(1):hi(1)+1)
    double precision :: ene_c

    double precision :: vsp(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))
    double precision :: dpy(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies)
    double precision :: dxe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies)
    double precision :: dpe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3))

    double precision :: Hg(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1,2:ncons)

    double precision :: M8p (8,lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)

    double precision :: sumdrY(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: sumrYv(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    double precision :: gradp (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

    do k=dlo(3),dhi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             vsp(i,j,k) = xi(i,j,k) + FourThirds*mu(i,j,k)
          enddo
       enddo
    enddo

    dpe = 0.d0

    do n=1,nspecies
       qxn = qx1+n-1
       qyn = qy1+n-1
       qhn = qh1+n-1
       do k=dlo(3),dhi(3)
          do j=dlo(2),dhi(2)
             do i=dlo(1),dhi(1)
                dpy(i,j,k,n) = dxy(i,j,k,n)/q(i,j,k,qpres)*(q(i,j,k,qxn)-q(i,j,k,qyn))
                dxe(i,j,k,n) = dxy(i,j,k,n)*q(i,j,k,qhn)
                dpe(i,j,k) = dpe(i,j,k) + dpy(i,j,k,n)*q(i,j,k,qhn)
             end do
          end do
       end do
    end do

    ! ------- BEGIN x-direction -------

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
!EXPAND             mmtmp(1:8,i) = matmul(vsp(i-4:i+3,j,k), M8)
             mmtmp(1,i) = vsp(i-4,j,k) * M8(1,1) &
                        + vsp(i-3,j,k) * M8(2,1) &
                        + vsp(i-2,j,k) * M8(3,1) &
                        + vsp(i-1,j,k) * M8(4,1) &
                        + vsp(i  ,j,k) * M8(5,1)
             mmtmp(2,i) = vsp(i-4,j,k) * M8(1,2) &
                        + vsp(i-3,j,k) * M8(2,2) &
                        + vsp(i-2,j,k) * M8(3,2) &
                        + vsp(i-1,j,k) * M8(4,2) &
                        + vsp(i  ,j,k) * M8(5,2) &
                        + vsp(i+1,j,k) * M8(6,2)
             mmtmp(3,i) = vsp(i-4,j,k) * M8(1,3) &
                        + vsp(i-3,j,k) * M8(2,3) &
                        + vsp(i-2,j,k) * M8(3,3) &
                        + vsp(i-1,j,k) * M8(4,3) &
                        + vsp(i  ,j,k) * M8(5,3) &
                        + vsp(i+1,j,k) * M8(6,3) &
                        + vsp(i+2,j,k) * M8(7,3)
             mmtmp(4,i) = vsp(i-4,j,k) * M8(1,4) &
                        + vsp(i-3,j,k) * M8(2,4) &
                        + vsp(i-2,j,k) * M8(3,4) &
                        + vsp(i-1,j,k) * M8(4,4) &
                        + vsp(i  ,j,k) * M8(5,4) &
                        + vsp(i+1,j,k) * M8(6,4) &
                        + vsp(i+2,j,k) * M8(7,4) &
                        + vsp(i+3,j,k) * M8(8,4)
             mmtmp(5,i) = vsp(i-4,j,k) * M8(1,5) &
                        + vsp(i-3,j,k) * M8(2,5) &
                        + vsp(i-2,j,k) * M8(3,5) &
                        + vsp(i-1,j,k) * M8(4,5) &
                        + vsp(i  ,j,k) * M8(5,5) &
                        + vsp(i+1,j,k) * M8(6,5) &
                        + vsp(i+2,j,k) * M8(7,5) &
                        + vsp(i+3,j,k) * M8(8,5)
             mmtmp(6,i) = vsp(i-3,j,k) * M8(2,6) &
                        + vsp(i-2,j,k) * M8(3,6) &
                        + vsp(i-1,j,k) * M8(4,6) &
                        + vsp(i  ,j,k) * M8(5,6) &
                        + vsp(i+1,j,k) * M8(6,6) &
                        + vsp(i+2,j,k) * M8(7,6) &
                        + vsp(i+3,j,k) * M8(8,6)
             mmtmp(7,i) = vsp(i-2,j,k) * M8(3,7) &
                        + vsp(i-1,j,k) * M8(4,7) &
                        + vsp(i  ,j,k) * M8(5,7) &
                        + vsp(i+1,j,k) * M8(6,7) &
                        + vsp(i+2,j,k) * M8(7,7) &
                        + vsp(i+3,j,k) * M8(8,7)
             mmtmp(8,i) = vsp(i-1,j,k) * M8(4,8) &
                        + vsp(i  ,j,k) * M8(5,8) &
                        + vsp(i+1,j,k) * M8(6,8) &
                        + vsp(i+2,j,k) * M8(7,8) &
                        + vsp(i+3,j,k) * M8(8,8)
!EXPAND             Hg(i,j,k,imx) = dot_product(mmtmp(1:8,i), q(i-4:i+3,j,k,qu))
             Hg(i,j,k,imx) =  &
                ( mmtmp(1,i)*q(i-4,j,k,qu) + mmtmp(2,i)*q(i-3,j,k,qu) &
                + mmtmp(3,i)*q(i-2,j,k,qu) + mmtmp(4,i)*q(i-1,j,k,qu) &
                + mmtmp(5,i)*q(i  ,j,k,qu) + mmtmp(6,i)*q(i+1,j,k,qu) &
                + mmtmp(7,i)*q(i+2,j,k,qu) + mmtmp(8,i)*q(i+3,j,k,qu) )
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
!EXPAND             mmtmp(1:8,i) = matmul(mu(i-4:i+3,j,k), M8)
             mmtmp(1,i) = mu(i-4,j,k) * M8(1,1) &
                        + mu(i-3,j,k) * M8(2,1) &
                        + mu(i-2,j,k) * M8(3,1) &
                        + mu(i-1,j,k) * M8(4,1) &
                        + mu(i  ,j,k) * M8(5,1)
             mmtmp(2,i) = mu(i-4,j,k) * M8(1,2) &
                        + mu(i-3,j,k) * M8(2,2) &
                        + mu(i-2,j,k) * M8(3,2) &
                        + mu(i-1,j,k) * M8(4,2) &
                        + mu(i  ,j,k) * M8(5,2) &
                        + mu(i+1,j,k) * M8(6,2)
             mmtmp(3,i) = mu(i-4,j,k) * M8(1,3) &
                        + mu(i-3,j,k) * M8(2,3) &
                        + mu(i-2,j,k) * M8(3,3) &
                        + mu(i-1,j,k) * M8(4,3) &
                        + mu(i  ,j,k) * M8(5,3) &
                        + mu(i+1,j,k) * M8(6,3) &
                        + mu(i+2,j,k) * M8(7,3)
             mmtmp(4,i) = mu(i-4,j,k) * M8(1,4) &
                        + mu(i-3,j,k) * M8(2,4) &
                        + mu(i-2,j,k) * M8(3,4) &
                        + mu(i-1,j,k) * M8(4,4) &
                        + mu(i  ,j,k) * M8(5,4) &
                        + mu(i+1,j,k) * M8(6,4) &
                        + mu(i+2,j,k) * M8(7,4) &
                        + mu(i+3,j,k) * M8(8,4)
             mmtmp(5,i) = mu(i-4,j,k) * M8(1,5) &
                        + mu(i-3,j,k) * M8(2,5) &
                        + mu(i-2,j,k) * M8(3,5) &
                        + mu(i-1,j,k) * M8(4,5) &
                        + mu(i  ,j,k) * M8(5,5) &
                        + mu(i+1,j,k) * M8(6,5) &
                        + mu(i+2,j,k) * M8(7,5) &
                        + mu(i+3,j,k) * M8(8,5)
             mmtmp(6,i) = mu(i-3,j,k) * M8(2,6) &
                        + mu(i-2,j,k) * M8(3,6) &
                        + mu(i-1,j,k) * M8(4,6) &
                        + mu(i  ,j,k) * M8(5,6) &
                        + mu(i+1,j,k) * M8(6,6) &
                        + mu(i+2,j,k) * M8(7,6) &
                        + mu(i+3,j,k) * M8(8,6)
             mmtmp(7,i) = mu(i-2,j,k) * M8(3,7) &
                        + mu(i-1,j,k) * M8(4,7) &
                        + mu(i  ,j,k) * M8(5,7) &
                        + mu(i+1,j,k) * M8(6,7) &
                        + mu(i+2,j,k) * M8(7,7) &
                        + mu(i+3,j,k) * M8(8,7)
             mmtmp(8,i) = mu(i-1,j,k) * M8(4,8) &
                        + mu(i  ,j,k) * M8(5,8) &
                        + mu(i+1,j,k) * M8(6,8) &
                        + mu(i+2,j,k) * M8(7,8) &
                        + mu(i+3,j,k) * M8(8,8)
!EXPAND             Hg(i,j,k,imy) = dot_product(mmtmp(1:8,i), q(i-4:i+3,j,k,qv))
             Hg(i,j,k,imy) =  &
                ( mmtmp(1,i)*q(i-4,j,k,qv) + mmtmp(2,i)*q(i-3,j,k,qv) &
                + mmtmp(3,i)*q(i-2,j,k,qv) + mmtmp(4,i)*q(i-1,j,k,qv) &
                + mmtmp(5,i)*q(i  ,j,k,qv) + mmtmp(6,i)*q(i+1,j,k,qv) &
                + mmtmp(7,i)*q(i+2,j,k,qv) + mmtmp(8,i)*q(i+3,j,k,qv) )
!EXPAND             Hg(i,j,k,imz) = dot_product(mmtmp(1:8,i), q(i-4:i+3,j,k,qw))
             Hg(i,j,k,imz) =  &
                ( mmtmp(1,i)*q(i-4,j,k,qw) + mmtmp(2,i)*q(i-3,j,k,qw) &
                + mmtmp(3,i)*q(i-2,j,k,qw) + mmtmp(4,i)*q(i-1,j,k,qw) &
                + mmtmp(5,i)*q(i  ,j,k,qw) + mmtmp(6,i)*q(i+1,j,k,qw) &
                + mmtmp(7,i)*q(i+2,j,k,qw) + mmtmp(8,i)*q(i+3,j,k,qw) )
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
!EXPAND             mmtmp(1:8,i) = matmul(lam(i-4:i+3,j,k), M8)
             mmtmp(1,i) = lam(i-4,j,k) * M8(1,1) &
                        + lam(i-3,j,k) * M8(2,1) &
                        + lam(i-2,j,k) * M8(3,1) &
                        + lam(i-1,j,k) * M8(4,1) &
                        + lam(i  ,j,k) * M8(5,1)
             mmtmp(2,i) = lam(i-4,j,k) * M8(1,2) &
                        + lam(i-3,j,k) * M8(2,2) &
                        + lam(i-2,j,k) * M8(3,2) &
                        + lam(i-1,j,k) * M8(4,2) &
                        + lam(i  ,j,k) * M8(5,2) &
                        + lam(i+1,j,k) * M8(6,2)
             mmtmp(3,i) = lam(i-4,j,k) * M8(1,3) &
                        + lam(i-3,j,k) * M8(2,3) &
                        + lam(i-2,j,k) * M8(3,3) &
                        + lam(i-1,j,k) * M8(4,3) &
                        + lam(i  ,j,k) * M8(5,3) &
                        + lam(i+1,j,k) * M8(6,3) &
                        + lam(i+2,j,k) * M8(7,3)
             mmtmp(4,i) = lam(i-4,j,k) * M8(1,4) &
                        + lam(i-3,j,k) * M8(2,4) &
                        + lam(i-2,j,k) * M8(3,4) &
                        + lam(i-1,j,k) * M8(4,4) &
                        + lam(i  ,j,k) * M8(5,4) &
                        + lam(i+1,j,k) * M8(6,4) &
                        + lam(i+2,j,k) * M8(7,4) &
                        + lam(i+3,j,k) * M8(8,4)
             mmtmp(5,i) = lam(i-4,j,k) * M8(1,5) &
                        + lam(i-3,j,k) * M8(2,5) &
                        + lam(i-2,j,k) * M8(3,5) &
                        + lam(i-1,j,k) * M8(4,5) &
                        + lam(i  ,j,k) * M8(5,5) &
                        + lam(i+1,j,k) * M8(6,5) &
                        + lam(i+2,j,k) * M8(7,5) &
                        + lam(i+3,j,k) * M8(8,5)
             mmtmp(6,i) = lam(i-3,j,k) * M8(2,6) &
                        + lam(i-2,j,k) * M8(3,6) &
                        + lam(i-1,j,k) * M8(4,6) &
                        + lam(i  ,j,k) * M8(5,6) &
                        + lam(i+1,j,k) * M8(6,6) &
                        + lam(i+2,j,k) * M8(7,6) &
                        + lam(i+3,j,k) * M8(8,6)
             mmtmp(7,i) = lam(i-2,j,k) * M8(3,7) &
                        + lam(i-1,j,k) * M8(4,7) &
                        + lam(i  ,j,k) * M8(5,7) &
                        + lam(i+1,j,k) * M8(6,7) &
                        + lam(i+2,j,k) * M8(7,7) &
                        + lam(i+3,j,k) * M8(8,7)
             mmtmp(8,i) = lam(i-1,j,k) * M8(4,8) &
                        + lam(i  ,j,k) * M8(5,8) &
                        + lam(i+1,j,k) * M8(6,8) &
                        + lam(i+2,j,k) * M8(7,8) &
                        + lam(i+3,j,k) * M8(8,8)
!EXPAND             Hg(i,j,k,iene) = dot_product(mmtmp(1:8,i), q(i-4:i+3,j,k,qtemp))
             Hg(i,j,k,iene) =  &
                ( mmtmp(1,i)*q(i-4,j,k,qtemp) + mmtmp(2,i)*q(i-3,j,k,qtemp) &
                + mmtmp(3,i)*q(i-2,j,k,qtemp) + mmtmp(4,i)*q(i-1,j,k,qtemp) &
                + mmtmp(5,i)*q(i  ,j,k,qtemp) + mmtmp(6,i)*q(i+1,j,k,qtemp) &
                + mmtmp(7,i)*q(i+2,j,k,qtemp) + mmtmp(8,i)*q(i+3,j,k,qtemp) )
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
!EXPAND             mmtmp(1:8,i) = matmul(M8, q(i-4:i+3,j,k,qpres))
             mmtmp(1,i) = M8T(1,1) * q(i-4,j,k,qpres) &
                        + M8T(2,1) * q(i-3,j,k,qpres) &
                        + M8T(3,1) * q(i-2,j,k,qpres) &
                        + M8T(4,1) * q(i-1,j,k,qpres) &
                        + M8T(5,1) * q(i  ,j,k,qpres)
             mmtmp(2,i) = M8T(1,2) * q(i-4,j,k,qpres) &
                        + M8T(2,2) * q(i-3,j,k,qpres) &
                        + M8T(3,2) * q(i-2,j,k,qpres) &
                        + M8T(4,2) * q(i-1,j,k,qpres) &
                        + M8T(5,2) * q(i  ,j,k,qpres) &
                        + M8T(6,2) * q(i+1,j,k,qpres)
             mmtmp(3,i) = M8T(1,3) * q(i-4,j,k,qpres) &
                        + M8T(2,3) * q(i-3,j,k,qpres) &
                        + M8T(3,3) * q(i-2,j,k,qpres) &
                        + M8T(4,3) * q(i-1,j,k,qpres) &
                        + M8T(5,3) * q(i  ,j,k,qpres) &
                        + M8T(6,3) * q(i+1,j,k,qpres) &
                        + M8T(7,3) * q(i+2,j,k,qpres)
             mmtmp(4,i) = M8T(1,4) * q(i-4,j,k,qpres) &
                        + M8T(2,4) * q(i-3,j,k,qpres) &
                        + M8T(3,4) * q(i-2,j,k,qpres) &
                        + M8T(4,4) * q(i-1,j,k,qpres) &
                        + M8T(5,4) * q(i  ,j,k,qpres) &
                        + M8T(6,4) * q(i+1,j,k,qpres) &
                        + M8T(7,4) * q(i+2,j,k,qpres) &
                        + M8T(8,4) * q(i+3,j,k,qpres)
             mmtmp(5,i) = M8T(1,5) * q(i-4,j,k,qpres) &
                        + M8T(2,5) * q(i-3,j,k,qpres) &
                        + M8T(3,5) * q(i-2,j,k,qpres) &
                        + M8T(4,5) * q(i-1,j,k,qpres) &
                        + M8T(5,5) * q(i  ,j,k,qpres) &
                        + M8T(6,5) * q(i+1,j,k,qpres) &
                        + M8T(7,5) * q(i+2,j,k,qpres) &
                        + M8T(8,5) * q(i+3,j,k,qpres)
             mmtmp(6,i) = M8T(2,6) * q(i-3,j,k,qpres) &
                        + M8T(3,6) * q(i-2,j,k,qpres) &
                        + M8T(4,6) * q(i-1,j,k,qpres) &
                        + M8T(5,6) * q(i  ,j,k,qpres) &
                        + M8T(6,6) * q(i+1,j,k,qpres) &
                        + M8T(7,6) * q(i+2,j,k,qpres) &
                        + M8T(8,6) * q(i+3,j,k,qpres)
             mmtmp(7,i) = M8T(3,7) * q(i-2,j,k,qpres) &
                        + M8T(4,7) * q(i-1,j,k,qpres) &
                        + M8T(5,7) * q(i  ,j,k,qpres) &
                        + M8T(6,7) * q(i+1,j,k,qpres) &
                        + M8T(7,7) * q(i+2,j,k,qpres) &
                        + M8T(8,7) * q(i+3,j,k,qpres)
             mmtmp(8,i) = M8T(4,8) * q(i-1,j,k,qpres) &
                        + M8T(5,8) * q(i  ,j,k,qpres) &
                        + M8T(6,8) * q(i+1,j,k,qpres) &
                        + M8T(7,8) * q(i+2,j,k,qpres) &
                        + M8T(8,8) * q(i+3,j,k,qpres)
!EXPAND             Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dpe(i-4:i+3,j,k), mmtmp(1:8,i))
             Hg(i,j,k,iene) = Hg(i,j,k,iene)+ &
                ( dpe(i-4,j,k)*mmtmp(1,i) + dpe(i-3,j,k)*mmtmp(2,i) &
                + dpe(i-2,j,k)*mmtmp(3,i) + dpe(i-1,j,k)*mmtmp(4,i) &
                + dpe(i  ,j,k)*mmtmp(5,i) + dpe(i+1,j,k)*mmtmp(6,i) &
                + dpe(i+2,j,k)*mmtmp(7,i) + dpe(i+3,j,k)*mmtmp(8,i) )
          end do
          do i=lo(1),hi(1)+1
             M8p(:,i,j,k) = mmtmp(:,i)             
          end do
       end do
    end do

    do n = 1, nspecies

       if (n .eq. iias) cycle  ! inactive speices

       qxn = qx1+n-1
       iryn = iry1+n-1

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1 
!EXPAND                Hg(i,j,k,iryn) = dot_product(dpy(i-4:i+3,j,k,n), M8p(:,i,j,k))
                Hg(i,j,k,iryn) =  &
                   ( dpy(i-4,j,k,n)*M8p(1,i,j,k) + dpy(i-3,j,k,n)*M8p(2,i,j,k) &
                   + dpy(i-2,j,k,n)*M8p(3,i,j,k) + dpy(i-1,j,k,n)*M8p(4,i,j,k) &
                   + dpy(i  ,j,k,n)*M8p(5,i,j,k) + dpy(i+1,j,k,n)*M8p(6,i,j,k) &
                   + dpy(i+2,j,k,n)*M8p(7,i,j,k) + dpy(i+3,j,k,n)*M8p(8,i,j,k) )
             end do
          end do
       end do

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)    
             do i=lo(1),hi(1)+1
!EXPAND                mmtmp(1:8,i) = matmul(M8, q(i-4:i+3,j,k,qxn))
                mmtmp(1,i) = M8T(1,1) * q(i-4,j,k,qxn) &
                           + M8T(2,1) * q(i-3,j,k,qxn) &
                           + M8T(3,1) * q(i-2,j,k,qxn) &
                           + M8T(4,1) * q(i-1,j,k,qxn) &
                           + M8T(5,1) * q(i  ,j,k,qxn)
                mmtmp(2,i) = M8T(1,2) * q(i-4,j,k,qxn) &
                           + M8T(2,2) * q(i-3,j,k,qxn) &
                           + M8T(3,2) * q(i-2,j,k,qxn) &
                           + M8T(4,2) * q(i-1,j,k,qxn) &
                           + M8T(5,2) * q(i  ,j,k,qxn) &
                           + M8T(6,2) * q(i+1,j,k,qxn)
                mmtmp(3,i) = M8T(1,3) * q(i-4,j,k,qxn) &
                           + M8T(2,3) * q(i-3,j,k,qxn) &
                           + M8T(3,3) * q(i-2,j,k,qxn) &
                           + M8T(4,3) * q(i-1,j,k,qxn) &
                           + M8T(5,3) * q(i  ,j,k,qxn) &
                           + M8T(6,3) * q(i+1,j,k,qxn) &
                           + M8T(7,3) * q(i+2,j,k,qxn)
                mmtmp(4,i) = M8T(1,4) * q(i-4,j,k,qxn) &
                           + M8T(2,4) * q(i-3,j,k,qxn) &
                           + M8T(3,4) * q(i-2,j,k,qxn) &
                           + M8T(4,4) * q(i-1,j,k,qxn) &
                           + M8T(5,4) * q(i  ,j,k,qxn) &
                           + M8T(6,4) * q(i+1,j,k,qxn) &
                           + M8T(7,4) * q(i+2,j,k,qxn) &
                           + M8T(8,4) * q(i+3,j,k,qxn)
                mmtmp(5,i) = M8T(1,5) * q(i-4,j,k,qxn) &
                           + M8T(2,5) * q(i-3,j,k,qxn) &
                           + M8T(3,5) * q(i-2,j,k,qxn) &
                           + M8T(4,5) * q(i-1,j,k,qxn) &
                           + M8T(5,5) * q(i  ,j,k,qxn) &
                           + M8T(6,5) * q(i+1,j,k,qxn) &
                           + M8T(7,5) * q(i+2,j,k,qxn) &
                           + M8T(8,5) * q(i+3,j,k,qxn)
                mmtmp(6,i) = M8T(2,6) * q(i-3,j,k,qxn) &
                           + M8T(3,6) * q(i-2,j,k,qxn) &
                           + M8T(4,6) * q(i-1,j,k,qxn) &
                           + M8T(5,6) * q(i  ,j,k,qxn) &
                           + M8T(6,6) * q(i+1,j,k,qxn) &
                           + M8T(7,6) * q(i+2,j,k,qxn) &
                           + M8T(8,6) * q(i+3,j,k,qxn)
                mmtmp(7,i) = M8T(3,7) * q(i-2,j,k,qxn) &
                           + M8T(4,7) * q(i-1,j,k,qxn) &
                           + M8T(5,7) * q(i  ,j,k,qxn) &
                           + M8T(6,7) * q(i+1,j,k,qxn) &
                           + M8T(7,7) * q(i+2,j,k,qxn) &
                           + M8T(8,7) * q(i+3,j,k,qxn)
                mmtmp(8,i) = M8T(4,8) * q(i-1,j,k,qxn) &
                           + M8T(5,8) * q(i  ,j,k,qxn) &
                           + M8T(6,8) * q(i+1,j,k,qxn) &
                           + M8T(7,8) * q(i+2,j,k,qxn) &
                           + M8T(8,8) * q(i+3,j,k,qxn)
!EXPAND                Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dxe(i-4:i+3,j,k,n), mmtmp(1:8,i))
                Hg(i,j,k,iene) = Hg(i,j,k,iene)+ &
                   ( dxe(i-4,j,k,n)*mmtmp(1,i) + dxe(i-3,j,k,n)*mmtmp(2,i) &
                   + dxe(i-2,j,k,n)*mmtmp(3,i) + dxe(i-1,j,k,n)*mmtmp(4,i) &
                   + dxe(i  ,j,k,n)*mmtmp(5,i) + dxe(i+1,j,k,n)*mmtmp(6,i) &
                   + dxe(i+2,j,k,n)*mmtmp(7,i) + dxe(i+3,j,k,n)*mmtmp(8,i) )
!EXPAND                Hg(i,j,k,iryn) = Hg(i,j,k,iryn) &
!EXPAND                     + dot_product(dxy(i-4:i+3,j,k,n), mmtmp(1:8,i))
                Hg(i,j,k,iryn) = Hg(i,j,k,iryn)+ &
                   ( dxy(i-4,j,k,n)*mmtmp(1,i) + dxy(i-3,j,k,n)*mmtmp(2,i) &
                   + dxy(i-2,j,k,n)*mmtmp(3,i) + dxy(i-1,j,k,n)*mmtmp(4,i) &
                   + dxy(i  ,j,k,n)*mmtmp(5,i) + dxy(i+1,j,k,n)*mmtmp(6,i) &
                   + dxy(i+2,j,k,n)*mmtmp(7,i) + dxy(i+3,j,k,n)*mmtmp(8,i) )
             end do
          end do
       end do
    end do

    ! add x-direction rhs

    do n=2,iene
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hg(i+1,j,k,n) - Hg(i,j,k,n)) * dx2inv(1)
             end do
          end do
       end do
    end do
       
    sumdrY = 0.d0
    do n=iry1,ncons

       if (n.eq.iry_ias) cycle

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                sumdry(i,j,k) = sumdry(i,j,k) + (Hg(i+1,j,k,n) - Hg(i,j,k,n)) * dx2inv(1)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hg(i+1,j,k,n) - Hg(i,j,k,n)) * dx2inv(1)
             end do
          end do
       end do

    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
!EXPAND             gradp(i,j,k) = dxinv(1) * first_deriv_8(q(i-4:i+4,j,k,qpres))
             gradp(i,j,k) = dxinv(1) * &
                ( D8(1)*(q(i+1,j,k,qpres)-q(i-1,j,k,qpres)) &
                + D8(2)*(q(i+2,j,k,qpres)-q(i-2,j,k,qpres)) &
                + D8(3)*(q(i+3,j,k,qpres)-q(i-3,j,k,qpres)) &
                + D8(4)*(q(i+4,j,k,qpres)-q(i-4,j,k,qpres)) )
          end do
       end do
    end do
       
    sumryv = 0.d0
    do n = 1, nspecies

       if (n.eq.iias) cycle

       qxn = qx1+n-1
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
!EXPAND                sumryv(i,j,k) = sumryv(i,j,k) + dpy(i,j,k,n)*gradp(i,j,k)  &
!EXPAND                     + dxy(i,j,k,n)*dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qxn))
                sumryv(i,j,k) = sumryv(i,j,k)+dpy(i,j,k,n)*gradp(i,j,k)+dxy(i,j,k,n) * dxinv(1) * &
                   ( D8(1)*(q(i+1,j,k,qxn)-q(i-1,j,k,qxn)) &
                   + D8(2)*(q(i+2,j,k,qxn)-q(i-2,j,k,qxn)) &
                   + D8(3)*(q(i+3,j,k,qxn)-q(i-3,j,k,qxn)) &
                   + D8(4)*(q(i+4,j,k,qxn)-q(i-4,j,k,qxn)) )
             end do
          end do
       end do

    end do

    n = iias
    qhn = qh1+n-1
    iryn = iry1+n-1
    
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ene_c = sumdry(i,j,k)*q(i,j,k,qhn)+sumryv(i,j,k) * dxinv(1) * &
                  ( D8(1)*(q(i+1,j,k,qhn)-q(i-1,j,k,qhn)) &
                  + D8(2)*(q(i+2,j,k,qhn)-q(i-2,j,k,qhn)) &
                  + D8(3)*(q(i+3,j,k,qhn)-q(i-3,j,k,qhn)) &
                  + D8(4)*(q(i+4,j,k,qhn)-q(i-4,j,k,qhn)) )
             rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
             rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdry(i,j,k)
          end do
       end do
    end do

    ! ------- END x-direction -------

    ! ------- BEGIN y-direction -------

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)             
!EXPAND             mmtmp(1:8,i) = matmul(mu(i,j-4:j+3,k), M8)
             mmtmp(1,i) = mu(i,j-4,k) * M8(1,1) &
                        + mu(i,j-3,k) * M8(2,1) &
                        + mu(i,j-2,k) * M8(3,1) &
                        + mu(i,j-1,k) * M8(4,1) &
                        + mu(i,j  ,k) * M8(5,1)
             mmtmp(2,i) = mu(i,j-4,k) * M8(1,2) &
                        + mu(i,j-3,k) * M8(2,2) &
                        + mu(i,j-2,k) * M8(3,2) &
                        + mu(i,j-1,k) * M8(4,2) &
                        + mu(i,j  ,k) * M8(5,2) &
                        + mu(i,j+1,k) * M8(6,2)
             mmtmp(3,i) = mu(i,j-4,k) * M8(1,3) &
                        + mu(i,j-3,k) * M8(2,3) &
                        + mu(i,j-2,k) * M8(3,3) &
                        + mu(i,j-1,k) * M8(4,3) &
                        + mu(i,j  ,k) * M8(5,3) &
                        + mu(i,j+1,k) * M8(6,3) &
                        + mu(i,j+2,k) * M8(7,3)
             mmtmp(4,i) = mu(i,j-4,k) * M8(1,4) &
                        + mu(i,j-3,k) * M8(2,4) &
                        + mu(i,j-2,k) * M8(3,4) &
                        + mu(i,j-1,k) * M8(4,4) &
                        + mu(i,j  ,k) * M8(5,4) &
                        + mu(i,j+1,k) * M8(6,4) &
                        + mu(i,j+2,k) * M8(7,4) &
                        + mu(i,j+3,k) * M8(8,4)
             mmtmp(5,i) = mu(i,j-4,k) * M8(1,5) &
                        + mu(i,j-3,k) * M8(2,5) &
                        + mu(i,j-2,k) * M8(3,5) &
                        + mu(i,j-1,k) * M8(4,5) &
                        + mu(i,j  ,k) * M8(5,5) &
                        + mu(i,j+1,k) * M8(6,5) &
                        + mu(i,j+2,k) * M8(7,5) &
                        + mu(i,j+3,k) * M8(8,5)
             mmtmp(6,i) = mu(i,j-3,k) * M8(2,6) &
                        + mu(i,j-2,k) * M8(3,6) &
                        + mu(i,j-1,k) * M8(4,6) &
                        + mu(i,j  ,k) * M8(5,6) &
                        + mu(i,j+1,k) * M8(6,6) &
                        + mu(i,j+2,k) * M8(7,6) &
                        + mu(i,j+3,k) * M8(8,6)
             mmtmp(7,i) = mu(i,j-2,k) * M8(3,7) &
                        + mu(i,j-1,k) * M8(4,7) &
                        + mu(i,j  ,k) * M8(5,7) &
                        + mu(i,j+1,k) * M8(6,7) &
                        + mu(i,j+2,k) * M8(7,7) &
                        + mu(i,j+3,k) * M8(8,7)
             mmtmp(8,i) = mu(i,j-1,k) * M8(4,8) &
                        + mu(i,j  ,k) * M8(5,8) &
                        + mu(i,j+1,k) * M8(6,8) &
                        + mu(i,j+2,k) * M8(7,8) &
                        + mu(i,j+3,k) * M8(8,8)
!EXPAND             Hg(i,j,k,imx) = dot_product(mmtmp(1:8,i), q(i,j-4:j+3,k,qu))
             Hg(i,j,k,imx) =  &
                ( mmtmp(1,i)*q(i,j-4,k,qu) + mmtmp(2,i)*q(i,j-3,k,qu) &
                + mmtmp(3,i)*q(i,j-2,k,qu) + mmtmp(4,i)*q(i,j-1,k,qu) &
                + mmtmp(5,i)*q(i,j  ,k,qu) + mmtmp(6,i)*q(i,j+1,k,qu) &
                + mmtmp(7,i)*q(i,j+2,k,qu) + mmtmp(8,i)*q(i,j+3,k,qu) )
!EXPAND             Hg(i,j,k,imz) = dot_product(mmtmp(1:8,i), q(i,j-4:j+3,k,qw))
             Hg(i,j,k,imz) =  &
                ( mmtmp(1,i)*q(i,j-4,k,qw) + mmtmp(2,i)*q(i,j-3,k,qw) &
                + mmtmp(3,i)*q(i,j-2,k,qw) + mmtmp(4,i)*q(i,j-1,k,qw) &
                + mmtmp(5,i)*q(i,j  ,k,qw) + mmtmp(6,i)*q(i,j+1,k,qw) &
                + mmtmp(7,i)*q(i,j+2,k,qw) + mmtmp(8,i)*q(i,j+3,k,qw) )
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
!EXPAND             mmtmp(1:8,i) = matmul(vsp(i,j-4:j+3,k), M8)
             mmtmp(1,i) = vsp(i,j-4,k) * M8(1,1) &
                        + vsp(i,j-3,k) * M8(2,1) &
                        + vsp(i,j-2,k) * M8(3,1) &
                        + vsp(i,j-1,k) * M8(4,1) &
                        + vsp(i,j  ,k) * M8(5,1)
             mmtmp(2,i) = vsp(i,j-4,k) * M8(1,2) &
                        + vsp(i,j-3,k) * M8(2,2) &
                        + vsp(i,j-2,k) * M8(3,2) &
                        + vsp(i,j-1,k) * M8(4,2) &
                        + vsp(i,j  ,k) * M8(5,2) &
                        + vsp(i,j+1,k) * M8(6,2)
             mmtmp(3,i) = vsp(i,j-4,k) * M8(1,3) &
                        + vsp(i,j-3,k) * M8(2,3) &
                        + vsp(i,j-2,k) * M8(3,3) &
                        + vsp(i,j-1,k) * M8(4,3) &
                        + vsp(i,j  ,k) * M8(5,3) &
                        + vsp(i,j+1,k) * M8(6,3) &
                        + vsp(i,j+2,k) * M8(7,3)
             mmtmp(4,i) = vsp(i,j-4,k) * M8(1,4) &
                        + vsp(i,j-3,k) * M8(2,4) &
                        + vsp(i,j-2,k) * M8(3,4) &
                        + vsp(i,j-1,k) * M8(4,4) &
                        + vsp(i,j  ,k) * M8(5,4) &
                        + vsp(i,j+1,k) * M8(6,4) &
                        + vsp(i,j+2,k) * M8(7,4) &
                        + vsp(i,j+3,k) * M8(8,4)
             mmtmp(5,i) = vsp(i,j-4,k) * M8(1,5) &
                        + vsp(i,j-3,k) * M8(2,5) &
                        + vsp(i,j-2,k) * M8(3,5) &
                        + vsp(i,j-1,k) * M8(4,5) &
                        + vsp(i,j  ,k) * M8(5,5) &
                        + vsp(i,j+1,k) * M8(6,5) &
                        + vsp(i,j+2,k) * M8(7,5) &
                        + vsp(i,j+3,k) * M8(8,5)
             mmtmp(6,i) = vsp(i,j-3,k) * M8(2,6) &
                        + vsp(i,j-2,k) * M8(3,6) &
                        + vsp(i,j-1,k) * M8(4,6) &
                        + vsp(i,j  ,k) * M8(5,6) &
                        + vsp(i,j+1,k) * M8(6,6) &
                        + vsp(i,j+2,k) * M8(7,6) &
                        + vsp(i,j+3,k) * M8(8,6)
             mmtmp(7,i) = vsp(i,j-2,k) * M8(3,7) &
                        + vsp(i,j-1,k) * M8(4,7) &
                        + vsp(i,j  ,k) * M8(5,7) &
                        + vsp(i,j+1,k) * M8(6,7) &
                        + vsp(i,j+2,k) * M8(7,7) &
                        + vsp(i,j+3,k) * M8(8,7)
             mmtmp(8,i) = vsp(i,j-1,k) * M8(4,8) &
                        + vsp(i,j  ,k) * M8(5,8) &
                        + vsp(i,j+1,k) * M8(6,8) &
                        + vsp(i,j+2,k) * M8(7,8) &
                        + vsp(i,j+3,k) * M8(8,8)
!EXPAND             Hg(i,j,k,imy) = dot_product(mmtmp(1:8,i), q(i,j-4:j+3,k,qv))
             Hg(i,j,k,imy) =  &
                ( mmtmp(1,i)*q(i,j-4,k,qv) + mmtmp(2,i)*q(i,j-3,k,qv) &
                + mmtmp(3,i)*q(i,j-2,k,qv) + mmtmp(4,i)*q(i,j-1,k,qv) &
                + mmtmp(5,i)*q(i,j  ,k,qv) + mmtmp(6,i)*q(i,j+1,k,qv) &
                + mmtmp(7,i)*q(i,j+2,k,qv) + mmtmp(8,i)*q(i,j+3,k,qv) )
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
!EXPAND             mmtmp(1:8,i) = matmul(lam(i,j-4:j+3,k), M8)
             mmtmp(1,i) = lam(i,j-4,k) * M8(1,1) &
                        + lam(i,j-3,k) * M8(2,1) &
                        + lam(i,j-2,k) * M8(3,1) &
                        + lam(i,j-1,k) * M8(4,1) &
                        + lam(i,j  ,k) * M8(5,1)
             mmtmp(2,i) = lam(i,j-4,k) * M8(1,2) &
                        + lam(i,j-3,k) * M8(2,2) &
                        + lam(i,j-2,k) * M8(3,2) &
                        + lam(i,j-1,k) * M8(4,2) &
                        + lam(i,j  ,k) * M8(5,2) &
                        + lam(i,j+1,k) * M8(6,2)
             mmtmp(3,i) = lam(i,j-4,k) * M8(1,3) &
                        + lam(i,j-3,k) * M8(2,3) &
                        + lam(i,j-2,k) * M8(3,3) &
                        + lam(i,j-1,k) * M8(4,3) &
                        + lam(i,j  ,k) * M8(5,3) &
                        + lam(i,j+1,k) * M8(6,3) &
                        + lam(i,j+2,k) * M8(7,3)
             mmtmp(4,i) = lam(i,j-4,k) * M8(1,4) &
                        + lam(i,j-3,k) * M8(2,4) &
                        + lam(i,j-2,k) * M8(3,4) &
                        + lam(i,j-1,k) * M8(4,4) &
                        + lam(i,j  ,k) * M8(5,4) &
                        + lam(i,j+1,k) * M8(6,4) &
                        + lam(i,j+2,k) * M8(7,4) &
                        + lam(i,j+3,k) * M8(8,4)
             mmtmp(5,i) = lam(i,j-4,k) * M8(1,5) &
                        + lam(i,j-3,k) * M8(2,5) &
                        + lam(i,j-2,k) * M8(3,5) &
                        + lam(i,j-1,k) * M8(4,5) &
                        + lam(i,j  ,k) * M8(5,5) &
                        + lam(i,j+1,k) * M8(6,5) &
                        + lam(i,j+2,k) * M8(7,5) &
                        + lam(i,j+3,k) * M8(8,5)
             mmtmp(6,i) = lam(i,j-3,k) * M8(2,6) &
                        + lam(i,j-2,k) * M8(3,6) &
                        + lam(i,j-1,k) * M8(4,6) &
                        + lam(i,j  ,k) * M8(5,6) &
                        + lam(i,j+1,k) * M8(6,6) &
                        + lam(i,j+2,k) * M8(7,6) &
                        + lam(i,j+3,k) * M8(8,6)
             mmtmp(7,i) = lam(i,j-2,k) * M8(3,7) &
                        + lam(i,j-1,k) * M8(4,7) &
                        + lam(i,j  ,k) * M8(5,7) &
                        + lam(i,j+1,k) * M8(6,7) &
                        + lam(i,j+2,k) * M8(7,7) &
                        + lam(i,j+3,k) * M8(8,7)
             mmtmp(8,i) = lam(i,j-1,k) * M8(4,8) &
                        + lam(i,j  ,k) * M8(5,8) &
                        + lam(i,j+1,k) * M8(6,8) &
                        + lam(i,j+2,k) * M8(7,8) &
                        + lam(i,j+3,k) * M8(8,8)
!EXPAND             Hg(i,j,k,iene) = dot_product(mmtmp(1:8,i), q(i,j-4:j+3,k,qtemp))
             Hg(i,j,k,iene) =  &
                ( mmtmp(1,i)*q(i,j-4,k,qtemp) + mmtmp(2,i)*q(i,j-3,k,qtemp) &
                + mmtmp(3,i)*q(i,j-2,k,qtemp) + mmtmp(4,i)*q(i,j-1,k,qtemp) &
                + mmtmp(5,i)*q(i,j  ,k,qtemp) + mmtmp(6,i)*q(i,j+1,k,qtemp) &
                + mmtmp(7,i)*q(i,j+2,k,qtemp) + mmtmp(8,i)*q(i,j+3,k,qtemp) )
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
!EXPAND             mmtmp(1:8,i) = matmul(M8, q(i,j-4:j+3,k,qpres))
             mmtmp(1,i) = M8T(1,1) * q(i,j-4,k,qpres) &
                        + M8T(2,1) * q(i,j-3,k,qpres) &
                        + M8T(3,1) * q(i,j-2,k,qpres) &
                        + M8T(4,1) * q(i,j-1,k,qpres) &
                        + M8T(5,1) * q(i,j  ,k,qpres)
             mmtmp(2,i) = M8T(1,2) * q(i,j-4,k,qpres) &
                        + M8T(2,2) * q(i,j-3,k,qpres) &
                        + M8T(3,2) * q(i,j-2,k,qpres) &
                        + M8T(4,2) * q(i,j-1,k,qpres) &
                        + M8T(5,2) * q(i,j  ,k,qpres) &
                        + M8T(6,2) * q(i,j+1,k,qpres)
             mmtmp(3,i) = M8T(1,3) * q(i,j-4,k,qpres) &
                        + M8T(2,3) * q(i,j-3,k,qpres) &
                        + M8T(3,3) * q(i,j-2,k,qpres) &
                        + M8T(4,3) * q(i,j-1,k,qpres) &
                        + M8T(5,3) * q(i,j  ,k,qpres) &
                        + M8T(6,3) * q(i,j+1,k,qpres) &
                        + M8T(7,3) * q(i,j+2,k,qpres)
             mmtmp(4,i) = M8T(1,4) * q(i,j-4,k,qpres) &
                        + M8T(2,4) * q(i,j-3,k,qpres) &
                        + M8T(3,4) * q(i,j-2,k,qpres) &
                        + M8T(4,4) * q(i,j-1,k,qpres) &
                        + M8T(5,4) * q(i,j  ,k,qpres) &
                        + M8T(6,4) * q(i,j+1,k,qpres) &
                        + M8T(7,4) * q(i,j+2,k,qpres) &
                        + M8T(8,4) * q(i,j+3,k,qpres)
             mmtmp(5,i) = M8T(1,5) * q(i,j-4,k,qpres) &
                        + M8T(2,5) * q(i,j-3,k,qpres) &
                        + M8T(3,5) * q(i,j-2,k,qpres) &
                        + M8T(4,5) * q(i,j-1,k,qpres) &
                        + M8T(5,5) * q(i,j  ,k,qpres) &
                        + M8T(6,5) * q(i,j+1,k,qpres) &
                        + M8T(7,5) * q(i,j+2,k,qpres) &
                        + M8T(8,5) * q(i,j+3,k,qpres)
             mmtmp(6,i) = M8T(2,6) * q(i,j-3,k,qpres) &
                        + M8T(3,6) * q(i,j-2,k,qpres) &
                        + M8T(4,6) * q(i,j-1,k,qpres) &
                        + M8T(5,6) * q(i,j  ,k,qpres) &
                        + M8T(6,6) * q(i,j+1,k,qpres) &
                        + M8T(7,6) * q(i,j+2,k,qpres) &
                        + M8T(8,6) * q(i,j+3,k,qpres)
             mmtmp(7,i) = M8T(3,7) * q(i,j-2,k,qpres) &
                        + M8T(4,7) * q(i,j-1,k,qpres) &
                        + M8T(5,7) * q(i,j  ,k,qpres) &
                        + M8T(6,7) * q(i,j+1,k,qpres) &
                        + M8T(7,7) * q(i,j+2,k,qpres) &
                        + M8T(8,7) * q(i,j+3,k,qpres)
             mmtmp(8,i) = M8T(4,8) * q(i,j-1,k,qpres) &
                        + M8T(5,8) * q(i,j  ,k,qpres) &
                        + M8T(6,8) * q(i,j+1,k,qpres) &
                        + M8T(7,8) * q(i,j+2,k,qpres) &
                        + M8T(8,8) * q(i,j+3,k,qpres)
!EXPAND             Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dpe(i,j-4:j+3,k), mmtmp(1:8,i))
             Hg(i,j,k,iene) = Hg(i,j,k,iene)+ &
                ( dpe(i,j-4,k)*mmtmp(1,i) + dpe(i,j-3,k)*mmtmp(2,i) &
                + dpe(i,j-2,k)*mmtmp(3,i) + dpe(i,j-1,k)*mmtmp(4,i) &
                + dpe(i,j  ,k)*mmtmp(5,i) + dpe(i,j+1,k)*mmtmp(6,i) &
                + dpe(i,j+2,k)*mmtmp(7,i) + dpe(i,j+3,k)*mmtmp(8,i) )
          end do
          do i=lo(1),hi(1)
             M8p(:,i,j,k) = mmtmp(:,i)
          end do
       end do
    end do

    do n = 1, nspecies

       if (n .eq. iias) cycle  ! inactive speices

       qxn = qx1+n-1
       iryn = iry1+n-1

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1) 
!EXPAND                Hg(i,j,k,iryn) = dot_product(dpy(i,j-4:j+3,k,n), M8p(:,i,j,k))
                Hg(i,j,k,iryn) =  &
                   ( dpy(i,j-4,k,n)*M8p(1,i,j,k) + dpy(i,j-3,k,n)*M8p(2,i,j,k) &
                   + dpy(i,j-2,k,n)*M8p(3,i,j,k) + dpy(i,j-1,k,n)*M8p(4,i,j,k) &
                   + dpy(i,j  ,k,n)*M8p(5,i,j,k) + dpy(i,j+1,k,n)*M8p(6,i,j,k) &
                   + dpy(i,j+2,k,n)*M8p(7,i,j,k) + dpy(i,j+3,k,n)*M8p(8,i,j,k) )
             end do
          end do
       end do
       
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
!EXPAND                mmtmp(1:8,i) = matmul(M8, q(i,j-4:j+3,k,qxn))
                mmtmp(1,i) = M8T(1,1) * q(i,j-4,k,qxn) &
                           + M8T(2,1) * q(i,j-3,k,qxn) &
                           + M8T(3,1) * q(i,j-2,k,qxn) &
                           + M8T(4,1) * q(i,j-1,k,qxn) &
                           + M8T(5,1) * q(i,j  ,k,qxn)
                mmtmp(2,i) = M8T(1,2) * q(i,j-4,k,qxn) &
                           + M8T(2,2) * q(i,j-3,k,qxn) &
                           + M8T(3,2) * q(i,j-2,k,qxn) &
                           + M8T(4,2) * q(i,j-1,k,qxn) &
                           + M8T(5,2) * q(i,j  ,k,qxn) &
                           + M8T(6,2) * q(i,j+1,k,qxn)
                mmtmp(3,i) = M8T(1,3) * q(i,j-4,k,qxn) &
                           + M8T(2,3) * q(i,j-3,k,qxn) &
                           + M8T(3,3) * q(i,j-2,k,qxn) &
                           + M8T(4,3) * q(i,j-1,k,qxn) &
                           + M8T(5,3) * q(i,j  ,k,qxn) &
                           + M8T(6,3) * q(i,j+1,k,qxn) &
                           + M8T(7,3) * q(i,j+2,k,qxn)
                mmtmp(4,i) = M8T(1,4) * q(i,j-4,k,qxn) &
                           + M8T(2,4) * q(i,j-3,k,qxn) &
                           + M8T(3,4) * q(i,j-2,k,qxn) &
                           + M8T(4,4) * q(i,j-1,k,qxn) &
                           + M8T(5,4) * q(i,j  ,k,qxn) &
                           + M8T(6,4) * q(i,j+1,k,qxn) &
                           + M8T(7,4) * q(i,j+2,k,qxn) &
                           + M8T(8,4) * q(i,j+3,k,qxn)
                mmtmp(5,i) = M8T(1,5) * q(i,j-4,k,qxn) &
                           + M8T(2,5) * q(i,j-3,k,qxn) &
                           + M8T(3,5) * q(i,j-2,k,qxn) &
                           + M8T(4,5) * q(i,j-1,k,qxn) &
                           + M8T(5,5) * q(i,j  ,k,qxn) &
                           + M8T(6,5) * q(i,j+1,k,qxn) &
                           + M8T(7,5) * q(i,j+2,k,qxn) &
                           + M8T(8,5) * q(i,j+3,k,qxn)
                mmtmp(6,i) = M8T(2,6) * q(i,j-3,k,qxn) &
                           + M8T(3,6) * q(i,j-2,k,qxn) &
                           + M8T(4,6) * q(i,j-1,k,qxn) &
                           + M8T(5,6) * q(i,j  ,k,qxn) &
                           + M8T(6,6) * q(i,j+1,k,qxn) &
                           + M8T(7,6) * q(i,j+2,k,qxn) &
                           + M8T(8,6) * q(i,j+3,k,qxn)
                mmtmp(7,i) = M8T(3,7) * q(i,j-2,k,qxn) &
                           + M8T(4,7) * q(i,j-1,k,qxn) &
                           + M8T(5,7) * q(i,j  ,k,qxn) &
                           + M8T(6,7) * q(i,j+1,k,qxn) &
                           + M8T(7,7) * q(i,j+2,k,qxn) &
                           + M8T(8,7) * q(i,j+3,k,qxn)
                mmtmp(8,i) = M8T(4,8) * q(i,j-1,k,qxn) &
                           + M8T(5,8) * q(i,j  ,k,qxn) &
                           + M8T(6,8) * q(i,j+1,k,qxn) &
                           + M8T(7,8) * q(i,j+2,k,qxn) &
                           + M8T(8,8) * q(i,j+3,k,qxn)
!EXPAND                Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dxe(i,j-4:j+3,k,n), mmtmp(1:8,i))
                Hg(i,j,k,iene) = Hg(i,j,k,iene)+ &
                   ( dxe(i,j-4,k,n)*mmtmp(1,i) + dxe(i,j-3,k,n)*mmtmp(2,i) &
                   + dxe(i,j-2,k,n)*mmtmp(3,i) + dxe(i,j-1,k,n)*mmtmp(4,i) &
                   + dxe(i,j  ,k,n)*mmtmp(5,i) + dxe(i,j+1,k,n)*mmtmp(6,i) &
                   + dxe(i,j+2,k,n)*mmtmp(7,i) + dxe(i,j+3,k,n)*mmtmp(8,i) )
!EXPAND                Hg(i,j,k,iryn) = Hg(i,j,k,iryn) &
!EXPAND                     + dot_product(dxy(i,j-4:j+3,k,n), mmtmp(1:8,i))
                Hg(i,j,k,iryn) = Hg(i,j,k,iryn)+ &
                   ( dxy(i,j-4,k,n)*mmtmp(1,i) + dxy(i,j-3,k,n)*mmtmp(2,i) &
                   + dxy(i,j-2,k,n)*mmtmp(3,i) + dxy(i,j-1,k,n)*mmtmp(4,i) &
                   + dxy(i,j  ,k,n)*mmtmp(5,i) + dxy(i,j+1,k,n)*mmtmp(6,i) &
                   + dxy(i,j+2,k,n)*mmtmp(7,i) + dxy(i,j+3,k,n)*mmtmp(8,i) )
             end do
          end do
       end do

    end do
       
    ! add y-direction rhs

    do n=2,iene
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hg(i,j+1,k,n) - Hg(i,j,k,n)) * dx2inv(2)
             end do
          end do
       end do
    end do

    sumdrY = 0.d0
    do n=iry1,ncons

       if (n.eq.iry_ias) cycle

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                sumdry(i,j,k) = sumdry(i,j,k) + (Hg(i,j+1,k,n) - Hg(i,j,k,n)) * dx2inv(2)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hg(i,j+1,k,n) - Hg(i,j,k,n)) * dx2inv(2)
             end do
          end do
       end do

    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             gradp(i,j,k) = dxinv(2) * &
                ( D8(1)*(q(i,j+1,k,qpres)-q(i,j-1,k,qpres)) &
                + D8(2)*(q(i,j+2,k,qpres)-q(i,j-2,k,qpres)) &
                + D8(3)*(q(i,j+3,k,qpres)-q(i,j-3,k,qpres)) &
                + D8(4)*(q(i,j+4,k,qpres)-q(i,j-4,k,qpres)) )
          end do
       end do
    end do
    
    sumryv = 0.d0
    do n = 1, nspecies

       if (n.eq.iias) cycle

       qxn = qx1+n-1
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
!EXPAND                sumryv(i,j,k) = sumryv(i,j,k) + dpy(i,j,k,n)*gradp(i,j,k)  &
!EXPAND                     + dxy(i,j,k,n)*dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qxn))
                sumryv(i,j,k) = sumryv(i,j,k)+dpy(i,j,k,n)*gradp(i,j,k)+dxy(i,j,k,n) * dxinv(2) * &
                   ( D8(1)*(q(i,j+1,k,qxn)-q(i,j-1,k,qxn)) &
                   + D8(2)*(q(i,j+2,k,qxn)-q(i,j-2,k,qxn)) &
                   + D8(3)*(q(i,j+3,k,qxn)-q(i,j-3,k,qxn)) &
                   + D8(4)*(q(i,j+4,k,qxn)-q(i,j-4,k,qxn)) )
             end do
          end do
       end do

    end do

    n = iias
    qhn = qh1+n-1
    iryn = iry1+n-1
    
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ene_c = sumdry(i,j,k)*q(i,j,k,qhn)+sumryv(i,j,k) * dxinv(2) * &
                  ( D8(1)*(q(i,j+1,k,qhn)-q(i,j-1,k,qhn)) &
                  + D8(2)*(q(i,j+2,k,qhn)-q(i,j-2,k,qhn)) &
                  + D8(3)*(q(i,j+3,k,qhn)-q(i,j-3,k,qhn)) &
                  + D8(4)*(q(i,j+4,k,qhn)-q(i,j-4,k,qhn)) )
             rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
             rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdry(i,j,k)
          end do
       end do
    end do

    ! ------- END y-direction -------

    ! ------- BEGIN z-direction -------

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
!EXPAND             mmtmp(1:8,i) = matmul(mu(i,j,k-4:k+3), M8)
             mmtmp(1,i) = mu(i,j,k-4) * M8(1,1) &
                        + mu(i,j,k-3) * M8(2,1) &
                        + mu(i,j,k-2) * M8(3,1) &
                        + mu(i,j,k-1) * M8(4,1) &
                        + mu(i,j,k  ) * M8(5,1)
             mmtmp(2,i) = mu(i,j,k-4) * M8(1,2) &
                        + mu(i,j,k-3) * M8(2,2) &
                        + mu(i,j,k-2) * M8(3,2) &
                        + mu(i,j,k-1) * M8(4,2) &
                        + mu(i,j,k  ) * M8(5,2) &
                        + mu(i,j,k+1) * M8(6,2)
             mmtmp(3,i) = mu(i,j,k-4) * M8(1,3) &
                        + mu(i,j,k-3) * M8(2,3) &
                        + mu(i,j,k-2) * M8(3,3) &
                        + mu(i,j,k-1) * M8(4,3) &
                        + mu(i,j,k  ) * M8(5,3) &
                        + mu(i,j,k+1) * M8(6,3) &
                        + mu(i,j,k+2) * M8(7,3)
             mmtmp(4,i) = mu(i,j,k-4) * M8(1,4) &
                        + mu(i,j,k-3) * M8(2,4) &
                        + mu(i,j,k-2) * M8(3,4) &
                        + mu(i,j,k-1) * M8(4,4) &
                        + mu(i,j,k  ) * M8(5,4) &
                        + mu(i,j,k+1) * M8(6,4) &
                        + mu(i,j,k+2) * M8(7,4) &
                        + mu(i,j,k+3) * M8(8,4)
             mmtmp(5,i) = mu(i,j,k-4) * M8(1,5) &
                        + mu(i,j,k-3) * M8(2,5) &
                        + mu(i,j,k-2) * M8(3,5) &
                        + mu(i,j,k-1) * M8(4,5) &
                        + mu(i,j,k  ) * M8(5,5) &
                        + mu(i,j,k+1) * M8(6,5) &
                        + mu(i,j,k+2) * M8(7,5) &
                        + mu(i,j,k+3) * M8(8,5)
             mmtmp(6,i) = mu(i,j,k-3) * M8(2,6) &
                        + mu(i,j,k-2) * M8(3,6) &
                        + mu(i,j,k-1) * M8(4,6) &
                        + mu(i,j,k  ) * M8(5,6) &
                        + mu(i,j,k+1) * M8(6,6) &
                        + mu(i,j,k+2) * M8(7,6) &
                        + mu(i,j,k+3) * M8(8,6)
             mmtmp(7,i) = mu(i,j,k-2) * M8(3,7) &
                        + mu(i,j,k-1) * M8(4,7) &
                        + mu(i,j,k  ) * M8(5,7) &
                        + mu(i,j,k+1) * M8(6,7) &
                        + mu(i,j,k+2) * M8(7,7) &
                        + mu(i,j,k+3) * M8(8,7)
             mmtmp(8,i) = mu(i,j,k-1) * M8(4,8) &
                        + mu(i,j,k  ) * M8(5,8) &
                        + mu(i,j,k+1) * M8(6,8) &
                        + mu(i,j,k+2) * M8(7,8) &
                        + mu(i,j,k+3) * M8(8,8)
!EXPAND             Hg(i,j,k,imx) = dot_product(mmtmp(1:8,i), q(i,j,k-4:k+3,qu))
             Hg(i,j,k,imx) =  &
                ( mmtmp(1,i)*q(i,j,k-4,qu) + mmtmp(2,i)*q(i,j,k-3,qu) &
                + mmtmp(3,i)*q(i,j,k-2,qu) + mmtmp(4,i)*q(i,j,k-1,qu) &
                + mmtmp(5,i)*q(i,j,k  ,qu) + mmtmp(6,i)*q(i,j,k+1,qu) &
                + mmtmp(7,i)*q(i,j,k+2,qu) + mmtmp(8,i)*q(i,j,k+3,qu) )
!EXPAND             Hg(i,j,k,imy) = dot_product(mmtmp(1:8,i), q(i,j,k-4:k+3,qv))
             Hg(i,j,k,imy) =  &
                ( mmtmp(1,i)*q(i,j,k-4,qv) + mmtmp(2,i)*q(i,j,k-3,qv) &
                + mmtmp(3,i)*q(i,j,k-2,qv) + mmtmp(4,i)*q(i,j,k-1,qv) &
                + mmtmp(5,i)*q(i,j,k  ,qv) + mmtmp(6,i)*q(i,j,k+1,qv) &
                + mmtmp(7,i)*q(i,j,k+2,qv) + mmtmp(8,i)*q(i,j,k+3,qv) )
          end do
       end do
    end do

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
!EXPAND             mmtmp(1:8,i) = matmul(vsp(i,j,k-4:k+3), M8)
             mmtmp(1,i) = vsp(i,j,k-4) * M8(1,1) &
                        + vsp(i,j,k-3) * M8(2,1) &
                        + vsp(i,j,k-2) * M8(3,1) &
                        + vsp(i,j,k-1) * M8(4,1) &
                        + vsp(i,j,k  ) * M8(5,1)
             mmtmp(2,i) = vsp(i,j,k-4) * M8(1,2) &
                        + vsp(i,j,k-3) * M8(2,2) &
                        + vsp(i,j,k-2) * M8(3,2) &
                        + vsp(i,j,k-1) * M8(4,2) &
                        + vsp(i,j,k  ) * M8(5,2) &
                        + vsp(i,j,k+1) * M8(6,2)
             mmtmp(3,i) = vsp(i,j,k-4) * M8(1,3) &
                        + vsp(i,j,k-3) * M8(2,3) &
                        + vsp(i,j,k-2) * M8(3,3) &
                        + vsp(i,j,k-1) * M8(4,3) &
                        + vsp(i,j,k  ) * M8(5,3) &
                        + vsp(i,j,k+1) * M8(6,3) &
                        + vsp(i,j,k+2) * M8(7,3)
             mmtmp(4,i) = vsp(i,j,k-4) * M8(1,4) &
                        + vsp(i,j,k-3) * M8(2,4) &
                        + vsp(i,j,k-2) * M8(3,4) &
                        + vsp(i,j,k-1) * M8(4,4) &
                        + vsp(i,j,k  ) * M8(5,4) &
                        + vsp(i,j,k+1) * M8(6,4) &
                        + vsp(i,j,k+2) * M8(7,4) &
                        + vsp(i,j,k+3) * M8(8,4)
             mmtmp(5,i) = vsp(i,j,k-4) * M8(1,5) &
                        + vsp(i,j,k-3) * M8(2,5) &
                        + vsp(i,j,k-2) * M8(3,5) &
                        + vsp(i,j,k-1) * M8(4,5) &
                        + vsp(i,j,k  ) * M8(5,5) &
                        + vsp(i,j,k+1) * M8(6,5) &
                        + vsp(i,j,k+2) * M8(7,5) &
                        + vsp(i,j,k+3) * M8(8,5)
             mmtmp(6,i) = vsp(i,j,k-3) * M8(2,6) &
                        + vsp(i,j,k-2) * M8(3,6) &
                        + vsp(i,j,k-1) * M8(4,6) &
                        + vsp(i,j,k  ) * M8(5,6) &
                        + vsp(i,j,k+1) * M8(6,6) &
                        + vsp(i,j,k+2) * M8(7,6) &
                        + vsp(i,j,k+3) * M8(8,6)
             mmtmp(7,i) = vsp(i,j,k-2) * M8(3,7) &
                        + vsp(i,j,k-1) * M8(4,7) &
                        + vsp(i,j,k  ) * M8(5,7) &
                        + vsp(i,j,k+1) * M8(6,7) &
                        + vsp(i,j,k+2) * M8(7,7) &
                        + vsp(i,j,k+3) * M8(8,7)
             mmtmp(8,i) = vsp(i,j,k-1) * M8(4,8) &
                        + vsp(i,j,k  ) * M8(5,8) &
                        + vsp(i,j,k+1) * M8(6,8) &
                        + vsp(i,j,k+2) * M8(7,8) &
                        + vsp(i,j,k+3) * M8(8,8)
!EXPAND             Hg(i,j,k,imz) = dot_product(mmtmp(1:8,i), q(i,j,k-4:k+3,qw))
             Hg(i,j,k,imz) =  &
                ( mmtmp(1,i)*q(i,j,k-4,qw) + mmtmp(2,i)*q(i,j,k-3,qw) &
                + mmtmp(3,i)*q(i,j,k-2,qw) + mmtmp(4,i)*q(i,j,k-1,qw) &
                + mmtmp(5,i)*q(i,j,k  ,qw) + mmtmp(6,i)*q(i,j,k+1,qw) &
                + mmtmp(7,i)*q(i,j,k+2,qw) + mmtmp(8,i)*q(i,j,k+3,qw) )
          end do
       end do
    end do

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
!EXPAND             mmtmp(1:8,i) = matmul(lam(i,j,k-4:k+3), M8)
             mmtmp(1,i) = lam(i,j,k-4) * M8(1,1) &
                        + lam(i,j,k-3) * M8(2,1) &
                        + lam(i,j,k-2) * M8(3,1) &
                        + lam(i,j,k-1) * M8(4,1) &
                        + lam(i,j,k  ) * M8(5,1)
             mmtmp(2,i) = lam(i,j,k-4) * M8(1,2) &
                        + lam(i,j,k-3) * M8(2,2) &
                        + lam(i,j,k-2) * M8(3,2) &
                        + lam(i,j,k-1) * M8(4,2) &
                        + lam(i,j,k  ) * M8(5,2) &
                        + lam(i,j,k+1) * M8(6,2)
             mmtmp(3,i) = lam(i,j,k-4) * M8(1,3) &
                        + lam(i,j,k-3) * M8(2,3) &
                        + lam(i,j,k-2) * M8(3,3) &
                        + lam(i,j,k-1) * M8(4,3) &
                        + lam(i,j,k  ) * M8(5,3) &
                        + lam(i,j,k+1) * M8(6,3) &
                        + lam(i,j,k+2) * M8(7,3)
             mmtmp(4,i) = lam(i,j,k-4) * M8(1,4) &
                        + lam(i,j,k-3) * M8(2,4) &
                        + lam(i,j,k-2) * M8(3,4) &
                        + lam(i,j,k-1) * M8(4,4) &
                        + lam(i,j,k  ) * M8(5,4) &
                        + lam(i,j,k+1) * M8(6,4) &
                        + lam(i,j,k+2) * M8(7,4) &
                        + lam(i,j,k+3) * M8(8,4)
             mmtmp(5,i) = lam(i,j,k-4) * M8(1,5) &
                        + lam(i,j,k-3) * M8(2,5) &
                        + lam(i,j,k-2) * M8(3,5) &
                        + lam(i,j,k-1) * M8(4,5) &
                        + lam(i,j,k  ) * M8(5,5) &
                        + lam(i,j,k+1) * M8(6,5) &
                        + lam(i,j,k+2) * M8(7,5) &
                        + lam(i,j,k+3) * M8(8,5)
             mmtmp(6,i) = lam(i,j,k-3) * M8(2,6) &
                        + lam(i,j,k-2) * M8(3,6) &
                        + lam(i,j,k-1) * M8(4,6) &
                        + lam(i,j,k  ) * M8(5,6) &
                        + lam(i,j,k+1) * M8(6,6) &
                        + lam(i,j,k+2) * M8(7,6) &
                        + lam(i,j,k+3) * M8(8,6)
             mmtmp(7,i) = lam(i,j,k-2) * M8(3,7) &
                        + lam(i,j,k-1) * M8(4,7) &
                        + lam(i,j,k  ) * M8(5,7) &
                        + lam(i,j,k+1) * M8(6,7) &
                        + lam(i,j,k+2) * M8(7,7) &
                        + lam(i,j,k+3) * M8(8,7)
             mmtmp(8,i) = lam(i,j,k-1) * M8(4,8) &
                        + lam(i,j,k  ) * M8(5,8) &
                        + lam(i,j,k+1) * M8(6,8) &
                        + lam(i,j,k+2) * M8(7,8) &
                        + lam(i,j,k+3) * M8(8,8)
!EXPAND             Hg(i,j,k,iene) = dot_product(mmtmp(1:8,i), q(i,j,k-4:k+3,qtemp))
             Hg(i,j,k,iene) =  &
                ( mmtmp(1,i)*q(i,j,k-4,qtemp) + mmtmp(2,i)*q(i,j,k-3,qtemp) &
                + mmtmp(3,i)*q(i,j,k-2,qtemp) + mmtmp(4,i)*q(i,j,k-1,qtemp) &
                + mmtmp(5,i)*q(i,j,k  ,qtemp) + mmtmp(6,i)*q(i,j,k+1,qtemp) &
                + mmtmp(7,i)*q(i,j,k+2,qtemp) + mmtmp(8,i)*q(i,j,k+3,qtemp) )
          end do
       end do
    end do

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
!EXPAND             mmtmp(1:8,i) = matmul(M8, q(i,j,k-4:k+3,qpres))
             mmtmp(1,i) = M8T(1,1) * q(i,j,k-4,qpres) &
                        + M8T(2,1) * q(i,j,k-3,qpres) &
                        + M8T(3,1) * q(i,j,k-2,qpres) &
                        + M8T(4,1) * q(i,j,k-1,qpres) &
                        + M8T(5,1) * q(i,j,k  ,qpres)
             mmtmp(2,i) = M8T(1,2) * q(i,j,k-4,qpres) &
                        + M8T(2,2) * q(i,j,k-3,qpres) &
                        + M8T(3,2) * q(i,j,k-2,qpres) &
                        + M8T(4,2) * q(i,j,k-1,qpres) &
                        + M8T(5,2) * q(i,j,k  ,qpres) &
                        + M8T(6,2) * q(i,j,k+1,qpres)
             mmtmp(3,i) = M8T(1,3) * q(i,j,k-4,qpres) &
                        + M8T(2,3) * q(i,j,k-3,qpres) &
                        + M8T(3,3) * q(i,j,k-2,qpres) &
                        + M8T(4,3) * q(i,j,k-1,qpres) &
                        + M8T(5,3) * q(i,j,k  ,qpres) &
                        + M8T(6,3) * q(i,j,k+1,qpres) &
                        + M8T(7,3) * q(i,j,k+2,qpres)
             mmtmp(4,i) = M8T(1,4) * q(i,j,k-4,qpres) &
                        + M8T(2,4) * q(i,j,k-3,qpres) &
                        + M8T(3,4) * q(i,j,k-2,qpres) &
                        + M8T(4,4) * q(i,j,k-1,qpres) &
                        + M8T(5,4) * q(i,j,k  ,qpres) &
                        + M8T(6,4) * q(i,j,k+1,qpres) &
                        + M8T(7,4) * q(i,j,k+2,qpres) &
                        + M8T(8,4) * q(i,j,k+3,qpres)
             mmtmp(5,i) = M8T(1,5) * q(i,j,k-4,qpres) &
                        + M8T(2,5) * q(i,j,k-3,qpres) &
                        + M8T(3,5) * q(i,j,k-2,qpres) &
                        + M8T(4,5) * q(i,j,k-1,qpres) &
                        + M8T(5,5) * q(i,j,k  ,qpres) &
                        + M8T(6,5) * q(i,j,k+1,qpres) &
                        + M8T(7,5) * q(i,j,k+2,qpres) &
                        + M8T(8,5) * q(i,j,k+3,qpres)
             mmtmp(6,i) = M8T(2,6) * q(i,j,k-3,qpres) &
                        + M8T(3,6) * q(i,j,k-2,qpres) &
                        + M8T(4,6) * q(i,j,k-1,qpres) &
                        + M8T(5,6) * q(i,j,k  ,qpres) &
                        + M8T(6,6) * q(i,j,k+1,qpres) &
                        + M8T(7,6) * q(i,j,k+2,qpres) &
                        + M8T(8,6) * q(i,j,k+3,qpres)
             mmtmp(7,i) = M8T(3,7) * q(i,j,k-2,qpres) &
                        + M8T(4,7) * q(i,j,k-1,qpres) &
                        + M8T(5,7) * q(i,j,k  ,qpres) &
                        + M8T(6,7) * q(i,j,k+1,qpres) &
                        + M8T(7,7) * q(i,j,k+2,qpres) &
                        + M8T(8,7) * q(i,j,k+3,qpres)
             mmtmp(8,i) = M8T(4,8) * q(i,j,k-1,qpres) &
                        + M8T(5,8) * q(i,j,k  ,qpres) &
                        + M8T(6,8) * q(i,j,k+1,qpres) &
                        + M8T(7,8) * q(i,j,k+2,qpres) &
                        + M8T(8,8) * q(i,j,k+3,qpres)
!EXPAND             Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dpe(i,j,k-4:k+3), mmtmp(1:8,i))
             Hg(i,j,k,iene) = Hg(i,j,k,iene)+ &
                ( dpe(i,j,k-4)*mmtmp(1,i) + dpe(i,j,k-3)*mmtmp(2,i) &
                + dpe(i,j,k-2)*mmtmp(3,i) + dpe(i,j,k-1)*mmtmp(4,i) &
                + dpe(i,j,k  )*mmtmp(5,i) + dpe(i,j,k+1)*mmtmp(6,i) &
                + dpe(i,j,k+2)*mmtmp(7,i) + dpe(i,j,k+3)*mmtmp(8,i) )
          end do
          do i=lo(1),hi(1)          
             M8p(:,i,j,k) = mmtmp(:,i) 
          end do
       end do
    end do

    do n = 1, nspecies

       if (n .eq. iias) cycle  ! inactive speices

       qxn = qx1+n-1
       iryn = iry1+n-1

       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
!EXPAND                Hg(i,j,k,iryn) = dot_product(dpy(i,j,k-4:k+3,n), M8p(:,i,j,k))
                Hg(i,j,k,iryn) =  &
                   ( dpy(i,j,k-4,n)*M8p(1,i,j,k) + dpy(i,j,k-3,n)*M8p(2,i,j,k) &
                   + dpy(i,j,k-2,n)*M8p(3,i,j,k) + dpy(i,j,k-1,n)*M8p(4,i,j,k) &
                   + dpy(i,j,k  ,n)*M8p(5,i,j,k) + dpy(i,j,k+1,n)*M8p(6,i,j,k) &
                   + dpy(i,j,k+2,n)*M8p(7,i,j,k) + dpy(i,j,k+3,n)*M8p(8,i,j,k) )
             end do
          end do
       end do

       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
!EXPAND                mmtmp(1:8,i) = matmul(M8, q(i,j,k-4:k+3,qxn))
                mmtmp(1,i) = M8T(1,1) * q(i,j,k-4,qxn) &
                           + M8T(2,1) * q(i,j,k-3,qxn) &
                           + M8T(3,1) * q(i,j,k-2,qxn) &
                           + M8T(4,1) * q(i,j,k-1,qxn) &
                           + M8T(5,1) * q(i,j,k  ,qxn)
                mmtmp(2,i) = M8T(1,2) * q(i,j,k-4,qxn) &
                           + M8T(2,2) * q(i,j,k-3,qxn) &
                           + M8T(3,2) * q(i,j,k-2,qxn) &
                           + M8T(4,2) * q(i,j,k-1,qxn) &
                           + M8T(5,2) * q(i,j,k  ,qxn) &
                           + M8T(6,2) * q(i,j,k+1,qxn)
                mmtmp(3,i) = M8T(1,3) * q(i,j,k-4,qxn) &
                           + M8T(2,3) * q(i,j,k-3,qxn) &
                           + M8T(3,3) * q(i,j,k-2,qxn) &
                           + M8T(4,3) * q(i,j,k-1,qxn) &
                           + M8T(5,3) * q(i,j,k  ,qxn) &
                           + M8T(6,3) * q(i,j,k+1,qxn) &
                           + M8T(7,3) * q(i,j,k+2,qxn)
                mmtmp(4,i) = M8T(1,4) * q(i,j,k-4,qxn) &
                           + M8T(2,4) * q(i,j,k-3,qxn) &
                           + M8T(3,4) * q(i,j,k-2,qxn) &
                           + M8T(4,4) * q(i,j,k-1,qxn) &
                           + M8T(5,4) * q(i,j,k  ,qxn) &
                           + M8T(6,4) * q(i,j,k+1,qxn) &
                           + M8T(7,4) * q(i,j,k+2,qxn) &
                           + M8T(8,4) * q(i,j,k+3,qxn)
                mmtmp(5,i) = M8T(1,5) * q(i,j,k-4,qxn) &
                           + M8T(2,5) * q(i,j,k-3,qxn) &
                           + M8T(3,5) * q(i,j,k-2,qxn) &
                           + M8T(4,5) * q(i,j,k-1,qxn) &
                           + M8T(5,5) * q(i,j,k  ,qxn) &
                           + M8T(6,5) * q(i,j,k+1,qxn) &
                           + M8T(7,5) * q(i,j,k+2,qxn) &
                           + M8T(8,5) * q(i,j,k+3,qxn)
                mmtmp(6,i) = M8T(2,6) * q(i,j,k-3,qxn) &
                           + M8T(3,6) * q(i,j,k-2,qxn) &
                           + M8T(4,6) * q(i,j,k-1,qxn) &
                           + M8T(5,6) * q(i,j,k  ,qxn) &
                           + M8T(6,6) * q(i,j,k+1,qxn) &
                           + M8T(7,6) * q(i,j,k+2,qxn) &
                           + M8T(8,6) * q(i,j,k+3,qxn)
                mmtmp(7,i) = M8T(3,7) * q(i,j,k-2,qxn) &
                           + M8T(4,7) * q(i,j,k-1,qxn) &
                           + M8T(5,7) * q(i,j,k  ,qxn) &
                           + M8T(6,7) * q(i,j,k+1,qxn) &
                           + M8T(7,7) * q(i,j,k+2,qxn) &
                           + M8T(8,7) * q(i,j,k+3,qxn)
                mmtmp(8,i) = M8T(4,8) * q(i,j,k-1,qxn) &
                           + M8T(5,8) * q(i,j,k  ,qxn) &
                           + M8T(6,8) * q(i,j,k+1,qxn) &
                           + M8T(7,8) * q(i,j,k+2,qxn) &
                           + M8T(8,8) * q(i,j,k+3,qxn)
!EXPAND                Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dxe(i,j,k-4:k+3,n), mmtmp(1:8,i))
                Hg(i,j,k,iene) = Hg(i,j,k,iene)+ &
                   ( dxe(i,j,k-4,n)*mmtmp(1,i) + dxe(i,j,k-3,n)*mmtmp(2,i) &
                   + dxe(i,j,k-2,n)*mmtmp(3,i) + dxe(i,j,k-1,n)*mmtmp(4,i) &
                   + dxe(i,j,k  ,n)*mmtmp(5,i) + dxe(i,j,k+1,n)*mmtmp(6,i) &
                   + dxe(i,j,k+2,n)*mmtmp(7,i) + dxe(i,j,k+3,n)*mmtmp(8,i) )
!EXPAND                Hg(i,j,k,iryn) = Hg(i,j,k,iryn) &
!EXPAND                     + dot_product(dxy(i,j,k-4:k+3,n), mmtmp(1:8,i))
                Hg(i,j,k,iryn) = Hg(i,j,k,iryn)+ &
                   ( dxy(i,j,k-4,n)*mmtmp(1,i) + dxy(i,j,k-3,n)*mmtmp(2,i) &
                   + dxy(i,j,k-2,n)*mmtmp(3,i) + dxy(i,j,k-1,n)*mmtmp(4,i) &
                   + dxy(i,j,k  ,n)*mmtmp(5,i) + dxy(i,j,k+1,n)*mmtmp(6,i) &
                   + dxy(i,j,k+2,n)*mmtmp(7,i) + dxy(i,j,k+3,n)*mmtmp(8,i) )
             end do
          end do
       end do
    end do

    ! add z-direction rhs

    do n=2,iene
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hg(i,j,k+1,n) - Hg(i,j,k,n)) * dx2inv(3)
             end do
          end do
       end do
    end do
    
    sumdrY = 0.d0
    do n=iry1,ncons

       if (n.eq.iry_ias) cycle

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                sumdry(i,j,k) = sumdry(i,j,k) + (Hg(i,j,k+1,n) - Hg(i,j,k,n)) * dx2inv(3)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hg(i,j,k+1,n) - Hg(i,j,k,n)) * dx2inv(3)
             end do
          end do
       end do

    end do
    
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             gradp(i,j,k) = dxinv(3) * &
                ( D8(1)*(q(i,j,k+1,qpres)-q(i,j,k-1,qpres)) &
                + D8(2)*(q(i,j,k+2,qpres)-q(i,j,k-2,qpres)) &
                + D8(3)*(q(i,j,k+3,qpres)-q(i,j,k-3,qpres)) &
                + D8(4)*(q(i,j,k+4,qpres)-q(i,j,k-4,qpres)) )
          end do
       end do
    end do
    
    sumryv = 0.d0
    do n = 1, nspecies

       if (n.eq.iias) cycle

       qxn = qx1+n-1
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
!EXPAND                sumryv(i,j,k) = sumryv(i,j,k) + dpy(i,j,k,n)*gradp(i,j,k)  &
!EXPAND                     + dxy(i,j,k,n)*dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qxn))
                sumryv(i,j,k) = sumryv(i,j,k)+dpy(i,j,k,n)*gradp(i,j,k)+dxy(i,j,k,n) * dxinv(3) * &
                   ( D8(1)*(q(i,j,k+1,qxn)-q(i,j,k-1,qxn)) &
                   + D8(2)*(q(i,j,k+2,qxn)-q(i,j,k-2,qxn)) &
                   + D8(3)*(q(i,j,k+3,qxn)-q(i,j,k-3,qxn)) &
                   + D8(4)*(q(i,j,k+4,qxn)-q(i,j,k-4,qxn)) )
             end do
          end do
       end do

    end do

    n = iias
    qhn = qh1+n-1
    iryn = iry1+n-1
    
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             ene_c = sumdry(i,j,k)*q(i,j,k,qhn)+sumryv(i,j,k) * dxinv(3) * &
                  ( D8(1)*(q(i,j,k+1,qhn)-q(i,j,k-1,qhn)) &
                  + D8(2)*(q(i,j,k+2,qhn)-q(i,j,k-2,qhn)) &
                  + D8(3)*(q(i,j,k+3,qhn)-q(i,j,k-3,qhn)) &
                  + D8(4)*(q(i,j,k+4,qhn)-q(i,j,k-4,qhn)) )
             rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
             rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdry(i,j,k)
          end do
       end do
    end do

    ! ------- END z-direction -------
    
    ! add kinetic energy
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) &
                  + rhs(i,j,k,imx)*q(i,j,k,qu) &
                  + rhs(i,j,k,imy)*q(i,j,k,qv) &
                  + rhs(i,j,k,imz)*q(i,j,k,qw)
          end do
       end do
    end do

  end subroutine diffterm_2


  subroutine chemterm_3d(lo,hi,q,qlo,qhi,up,uplo,uphi) bind(c,name='chemterm_3d')
    integer,         intent(in):: lo(3),hi(3),qlo(3),qhi(3),uplo(3),uphi(3)
    double precision,intent(in):: q ( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),nprim)
    double precision           :: up(uplo(1):uphi(1),uplo(2):uphi(2),uplo(3):uphi(3),ncons)

    integer :: iwrk, i,j,k,n,np
    double precision :: Yt(lo(1):hi(1),nspecies), wdot(lo(1):hi(1),nspecies), rwrk

    np = hi(1) - lo(1) + 1

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do n=1, nspecies
             do i=lo(1),hi(1)
                Yt(i,n) = q(i,j,k,qy1+n-1)
             end do
          end do

          call vckwyr(np, q(lo(1),j,k,qrho), q(lo(1),j,k,qtemp), Yt, iwrk, rwrk, wdot)

          do n=1, nspecies
             do i=lo(1),hi(1)
                up(i,j,k,iry1+n-1) = up(i,j,k,iry1+n-1) + wdot(i,n) * molecular_weight(n)
             end do
          end do
          
       end do
    end do

  end subroutine chemterm_3d


  subroutine comp_courno_3d(lo,hi,dx,Q,qlo,qhi,courno) bind(c,name='comp_courno_3d')
    integer, intent(in) :: lo(3), hi(3), qlo(3), qhi(3)
    double precision, intent(in) :: dx(3)
    double precision, intent(in) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)
    double precision, intent(inout) :: courno

    integer :: i,j,k, iwrk
    double precision :: dxinv(3), c, rwrk, Cv, Cp
    double precision :: Tt, X(nspecies), gamma
    double precision :: courx, coury, courz

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             Tt = q(i,j,k,qtemp)
             X  = q(i,j,k,qx1:qx1+nspecies-1)
             call ckcvbl(Tt, X, iwrk, rwrk, Cv)
             Cp = Cv + Ru
             gamma = Cp / Cv
             c = sqrt(gamma*q(i,j,k,qpres)/q(i,j,k,qrho))

             courx = (c+abs(q(i,j,k,qu))) * dxinv(1)
             coury = (c+abs(q(i,j,k,qv))) * dxinv(2)
             courz = (c+abs(q(i,j,k,qw))) * dxinv(3)
             
             courno = max( courx, coury, courz , courno )

          end do
       end do
    end do

  end subroutine comp_courno_3d

end module kernels_module
