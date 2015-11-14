
! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine advect(time,lo,hi,domlo,domhi, &
                        uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                        uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                        ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                        ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2, &
                        flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                        flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                        dx,dt)

      implicit none

      integer lo(2),hi(2)
      integer domlo(2),domhi(2)
      integer uin_l1,uin_l2,uin_h1,uin_h2
      integer uout_l1,uout_l2,uout_h1,uout_h2
      integer ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2
      integer ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2
      integer flux1_l1,flux1_l2,flux1_h1,flux1_h2
      integer flux2_l1,flux2_l2,flux2_h1,flux2_h2

      double precision uin(uin_l1:uin_h1,uin_l2:uin_h2)
      double precision uout(uout_l1:uout_h1,uout_l2:uout_h2)
      double precision ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2)
      double precision ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2)
      double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2)
      double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2)
      double precision dx(2),dt,time

      call set_ugdv(ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                    ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2)

!     Compute hyperbolic fluxes using unsplit Godunov
      call umeth2d(uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                   uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                   lo,hi,dx,dt, &
                   flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                   flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                   ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                   ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2)

      end subroutine advect

! ::: ---------------------------------------------------------------
! ::: :: UMETH2D     Compute fluxes using unsplit second
! ::: ::               order Godunov integrator.
! ::: ----------------------------------------------------------------

      subroutine umeth2d(uin, uin_l1, uin_l2, uin_h1, uin_h2,&
                         uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                         lo, hi, dx, dt, &
                         flux1, fd1_l1, fd1_l2, fd1_h1, fd1_h2, &
                         flux2, fd2_l1, fd2_l2, fd2_h1, fd2_h2, &
                         ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                         ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2)

      use meth_params_module, only : UFA, NVAR

      implicit none

      integer uin_l1, uin_l2, uin_h1, uin_h2
      integer uout_l1,uout_l2,uout_h1,uout_h2
      integer fd1_l1, fd1_l2, fd1_h1, fd1_h2
      integer fd2_l1, fd2_l2, fd2_h1, fd2_h2
      integer ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2
      integer ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2
      integer lo(2),hi(2)

      double precision dx(2),dt
      double precision  uin( uin_l1: uin_h1, uin_l2: uin_h2,NVAR)
      double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
      double precision ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2)
      double precision ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2)
      double precision flux1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,NVAR)
      double precision flux2(fd2_l1:fd2_h1,fd2_l2:fd2_h2,NVAR)

!     Left and right state arrays (edge centered, cell centered)
      double precision, allocatable:: dq(:,:,:),  qm(:,:,:),   qp(:,:,:)
      double precision, allocatable::qxm(:,:,:),qym(:,:,:)
      double precision, allocatable::qxp(:,:,:),qyp(:,:,:)

!     Temporary flux arrays
      double precision, allocatable::   fx(:,:,:),  fy(:,:,:)

!     Local scalar variables
      double precision :: hdtdx, hdt, hdtdy, qavg, eps
      integer          :: i,j,n,iadv

      allocate (  dq(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NVAR) )
      allocate (  qm(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NVAR) )
      allocate (  qp(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NVAR) )
      allocate ( qxm(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NVAR) )
      allocate ( qxp(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NVAR) )
      allocate ( qym(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NVAR) )
      allocate ( qyp(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,NVAR) )
      allocate (  fx(lo(1)  :hi(1)+1,lo(2)-1:hi(2)+1,NVAR))
      allocate (  fy(lo(1)-1:hi(1)+1,lo(2)  :hi(2)+1,NVAR))

!     Local constants
      hdtdx = 0.5d0*dt/dx(1)
      hdtdy = 0.5d0*dt/dx(2)
      hdt = 0.5d0*dt
      eps = 1.d-8

!     Trace to edges w/o transverse flux correction terms
      call trace(uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                 dq,qxm,qxp,qym,qyp,lo(1)-1,lo(2)-1,hi(1)+2,hi(2)+2, &
                 lo(1),lo(2),hi(1),hi(2),dx(1),dx(2),dt)

      do iadv = 1, NVAR
          n  = UFA + iadv - 1

          ! Upwind in x-direction to create temporary fluxes on x-edges
          do j = lo(2)  ,hi(2)
          do i = lo(1)-1,hi(1)+1
               if (abs(ugdx(i,j)) .lt. eps) then
                  qavg = 0.5d0 * (qxm(i,j,n) + qxp(i,j,n))
                  fx(i,j,n) = qavg*ugdx(i,j)
               else if (ugdx(i,j) .gt. 0.d0) then
                  fx(i,j,n) = qxm(i,j,n)*ugdx(i,j)
               else 
                  fx(i,j,n) = qxp(i,j,n)*ugdx(i,j)
               endif
          enddo
          enddo

          ! Upwind in y-direction to create temporary fluxes on y-edges
          do j = lo(2)-1,hi(2)+1
          do i = lo(1)  ,hi(1)
               if (abs(ugdy(i,j)) .lt. eps) then
                  qavg = 0.5d0 * (qym(i,j,n) + qyp(i,j,n))
                  fy(i,j,n) = qavg*ugdy(i,j)
               else if (ugdy(i,j) .gt. 0.d0) then
                  fy(i,j,n) = qym(i,j,n)*ugdy(i,j)
               else 
                  fy(i,j,n) = qyp(i,j,n)*ugdy(i,j)
               endif
          enddo
          enddo

          ! Use the fluxes on y-edges to add transverse contributions on x-edges
          do j = lo(2)  ,hi(2)
          do i = lo(1)-1,hi(1)+1
 
              qp(i  ,j,n) = qxp(i  ,j,n) - hdtdy*(fy(i,j+1,n)-fy(i,j,n))
              qm(i+1,j,n) = qxm(i+1,j,n) - hdtdy*(fy(i,j+1,n)-fy(i,j,n))
 
          enddo
          enddo

          ! Upwind in x-direction to create final fluxes on x-edges
          do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
               if (abs(ugdx(i,j)) .lt. eps) then
                  qavg = 0.5d0 * (qm(i,j,n) + qp(i,j,n))
                  flux1(i,j,n) = qavg*ugdx(i,j)
               else if (ugdx(i,j) .gt. 0.d0) then
                  flux1(i,j,n) = qm(i,j,n)*ugdx(i,j)
               else 
                  flux1(i,j,n) = qp(i,j,n)*ugdx(i,j)
               endif
          enddo
          enddo

          ! Use the fluxes on x-edges to add transverse contributions on y-edges
          do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
 
              qp(i,j  ,n) = qyp(i,j  ,n) - hdtdx*(fx(i+1,j,n)-fx(i,j,n))
              qm(i,j+1,n) = qym(i,j+1,n) - hdtdx*(fx(i+1,j,n)-fx(i,j,n))
 
          enddo
          enddo

          ! Upwind in y-direction to create final fluxes on y-edges
          do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
               if (abs(ugdy(i,j)) .lt. eps) then
                  qavg = 0.5d0 * (qm(i,j,n) + qp(i,j,n))
                  flux2(i,j,n) = qavg*ugdy(i,j)
               else if (ugdy(i,j) .gt. 0.d0) then
                  flux2(i,j,n) = qm(i,j,n)*ugdy(i,j)
               else 
                  flux2(i,j,n) = qp(i,j,n)*ugdy(i,j)
               endif
          enddo
          enddo

      enddo

      deallocate(dq,qm,qp,qxm,qxp,qym,qyp)
      deallocate(fx,fy)

      ! Do a conservative update
      do iadv = 1, NVAR
         n = UFA + iadv - 1
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               uout(i,j,n) = uin(i,j,n) + dt * &
                  ( (flux1(i,j,n) - flux1(i+1,j,n))/dx(1) &
                   +(flux2(i,j,n) - flux2(i,j+1,n))/dx(2) )
            enddo
         enddo
      enddo

      ! Scale by face area in order to correctly reflux
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)+1
        flux1(i,j,1:NVAR) = dt * flux1(i,j,1:NVAR) * dx(2)
      enddo
      enddo

      ! Scale by face area in order to correctly reflux
      do j = lo(2),hi(2)+1 
      do i = lo(1),hi(1)
        flux2(i,j,1:NVAR) = dt * flux2(i,j,1:NVAR) * dx(1)
      enddo
      enddo

      end subroutine umeth2d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine set_ugdv(ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                          ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2)

      use probdata_module, only : uadv, vadv

      integer ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2
      integer ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2

      double precision ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2)
      double precision ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2)

      do j = ugdx_l2,ugdx_h2
      do i = ugdx_l1,ugdx_h1
         ugdx(i,j) = uadv
      end do
      end do

      do j = ugdy_l2,ugdy_h2
      do i = ugdy_l1,ugdy_h1
         ugdy(i,j) = vadv
      end do
      end do

      end subroutine set_ugdv
