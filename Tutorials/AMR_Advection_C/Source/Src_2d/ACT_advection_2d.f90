
! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine advect(time,lo,hi,domlo,domhi, &
                        uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                        uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                        ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                        ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2, &
                        src,src_l1,src_l2,src_h1,src_h2, &
                        flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                        flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                        delta,dt,verbose)

      use meth_params_module, only : URHO, QVAR, NVAR, NHYP, normalize_species

      implicit none

      integer lo(2),hi(2),verbose
      integer domlo(2),domhi(2)
      integer uin_l1,uin_l2,uin_h1,uin_h2
      integer uout_l1,uout_l2,uout_h1,uout_h2
      integer ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2
      integer ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2
      integer flux1_l1,flux1_l2,flux1_h1,flux1_h2
      integer flux2_l1,flux2_l2,flux2_h1,flux2_h2
      integer src_l1,src_l2,src_h1,src_h2

      double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
      double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
      double precision ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2)
      double precision ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2)
      double precision src(src_l1:src_h1,src_l2:src_h2,NVAR)
      double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
      double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
      double precision delta(2),dt,time

!     Automatic arrays for workspace
      double precision, allocatable:: q(:,:,:)
      double precision, allocatable:: div(:,:)
      double precision, allocatable:: srcQ(:,:,:)

      integer ngq,ngf
!     integer i_c,j_c

      double precision dx,dy

      allocate(     q(uin_l1:uin_h1,uin_l2:uin_h2,QVAR))
      allocate(  srcQ(src_l1:src_h1,src_l2:src_h2,QVAR))
      allocate(   div(lo(1)  :hi(1)+1,lo(2)  :hi(2)+1))

      dx = delta(1)
      dy = delta(2)

      ngq = NHYP
      ngf = 1

!     Translate to primitive variables, compute sound speeds
      call ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                   q,uin_l1,uin_l2,uin_h1,uin_h2, &
                   src,srcQ,src_l1,src_l2,src_h1,src_h2, &
                   dx,dy,dt,ngq,ngf)

!     Compute hyperbolic fluxes using unsplit Godunov
      call umeth2d(q,uin_l1,uin_l2,uin_h1,uin_h2, &
                   srcQ, src_l1, src_l2, src_h1, src_h2,  &
                   lo(1),lo(2),hi(1),hi(2),dx,dy,dt, &
                   flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                   flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                   ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                   ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2)

!     Compute divergence of velocity field (on surroundingNodes(lo,hi))
      call divu(lo,hi,q,uin_l1,uin_l2,uin_h1,uin_h2, &
                delta,div,lo(1),lo(2),hi(1)+1,hi(2)+1)

!     Conservative update
      call consup(uin,    uin_l1,  uin_l2,  uin_h1,  uin_h2, &
                  uout,  uout_l1, uout_l2, uout_h1, uout_h2, &
                  src,    src_l1,  src_l2,  src_h1,  src_h2, &
                  flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                  flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                  div,lo,hi,delta,dt)

      ! Enforce the species >= 0
      call enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi)

      ! Normalize the species 
      if (normalize_species .eq. 1) &
         call normalize_new_species(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi)

      deallocate(q,div)

      end subroutine advect

! ::: ---------------------------------------------------------------
! ::: :: UMETH2D     Compute hyperbolic fluxes using unsplit second
! ::: ::               order Godunov integrator.
! ::: :: 
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
! ::: :: src         => (const)  source
! ::: :: nx          => (const)  number of cells in X direction
! ::: :: ny          => (const)  number of cells in Y direction
! ::: :: dx          => (const)  grid spacing in X direction
! ::: :: dy          => (const)  grid spacing in Y direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux1      <=  (modify) flux in X direction on X edges
! ::: :: flux2      <=  (modify) flux in Y direction on Y edges
! ::: ----------------------------------------------------------------

      subroutine umeth2d(q, qd_l1, qd_l2, qd_h1, qd_h2,&
                         srcQ, src_l1, src_l2, src_h1, src_h2, &
                         ilo1, ilo2, ihi1, ihi2, dx, dy, dt, &
                         flux1, fd1_l1, fd1_l2, fd1_h1, fd1_h2, &
                         flux2, fd2_l1, fd2_l2, fd2_h1, fd2_h2, &
                         ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                         ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2)

      use network, only : nspec
      use meth_params_module, only : QU, QV, QVAR, NVAR

      implicit none

      integer qd_l1, qd_l2, qd_h1, qd_h2
      integer src_l1, src_l2, src_h1, src_h2
      integer fd1_l1, fd1_l2, fd1_h1, fd1_h2
      integer fd2_l1, fd2_l2, fd2_h1, fd2_h2
      integer ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2
      integer ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2
      integer ilo1, ilo2, ihi1, ihi2

      double precision dx, dy, dt
      double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision  srcQ(src_l1:src_h1,src_l2:src_h2)
      double precision ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2)
      double precision ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2)
      double precision flux1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,NVAR)
      double precision flux2(fd2_l1:fd2_h1,fd2_l2:fd2_h2,NVAR)

!     Left and right state arrays (edge centered, cell centered)
      double precision, allocatable:: dq(:,:,:),  qm(:,:,:),   qp(:,:,:)
      double precision, allocatable::qxm(:,:,:),qym(:,:,:)
      double precision, allocatable::qxp(:,:,:),qyp(:,:,:)

!     Work arrays to hold 3 planes of riemann state and conservative fluxes
      double precision, allocatable::   fx(:,:,:),  fy(:,:,:)

!     Local scalar variables
      double precision :: hdtdx, hdt, hdtdy
      integer          :: i,j

      allocate (  dq(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
      allocate (  qm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
      allocate (  qp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
      allocate ( qxm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
      allocate ( qxp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
      allocate ( qym(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
      allocate ( qyp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
      allocate (  fx(ilo1  :ihi1+1,ilo2-1:ihi2+1,NVAR))
      allocate (  fy(ilo1-1:ihi1+1,ilo2  :ihi2+1,NVAR))

!     Local constants
      hdtdx = 0.5d0*dt/dx
      hdtdy = 0.5d0*dt/dy
      hdt = 0.5d0*dt

      do j = ugdx_l2,ugdx_h2
      do i = ugdx_l1,ugdx_h1
         ugdx(i,j) = 0.5d0 * (q(i,j,QU) + q(i-1,j,QU))
      end do
      end do

      do j = ugdy_l2,ugdy_h2
      do i = ugdy_l1,ugdy_h1
         ugdy(i,j) = 0.5d0 * (q(i,j,QV) + q(i,j-1,QV))
      end do
      end do

!     Trace to edges w/o transverse flux correction terms
      call trace(q,qd_l1,qd_l2,qd_h1,qd_h2, &
                 dq,qxm,qxp,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                 ilo1,ilo2,ihi1,ihi2,dx,dy,dt)

!     Upwind on x-edges
      call upwind(qxm, qxp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                  fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
                  ugdx, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                  1, ilo1, ihi1, ilo2-1, ihi2+1)

!     Upwind on y-edges
      call upwind(qym, qyp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                  fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
                  ugdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                  2, ilo1-1, ihi1+1, ilo2, ihi2)

!     Use edge states to create transverse derivative in y-direction
      call transy(qxm, qm, qxp, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                  fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
                  srcQ, src_l1, src_l2, src_h1, src_h2, &
                  hdt, hdtdy, ilo1-1, ihi1+1, ilo2, ihi2)

!     Upwind on x-edges to create final fluxes
      call upwind(qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                  flux1, fd1_l1, fd1_l2, fd1_h1, fd1_h2, &
                  ugdx, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                  1, ilo1, ihi1, ilo2, ihi2)
      
!     Use edge states to create transverse derivative in x-direction
      call transx(qym, qm,qyp,qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                  fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
                  srcQ,  src_l1,  src_l2,  src_h1,  src_h2, &
                  hdt, hdtdx, ilo1, ihi1, ilo2-1, ihi2+1)

!     Upwind on y-edges to create final fluxes
      call upwind(qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                  flux2, fd2_l1, fd2_l2, fd2_h1, fd2_h2, &
                  ugdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                  2, ilo1, ihi1, ilo2, ihi2)
      
      deallocate(dq,qm,qp,qxm,qxp,qym,qyp)
      deallocate(fx,fy)

      end subroutine umeth2d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ctoprim(lo,hi, &
                         uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                         q,q_l1,q_l2,q_h1,q_h2, &
                         src,srcQ,src_l1,src_l2,src_h1,src_h2, &
                         dx,dy,dt,ngp,ngf)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UX, UY, UFA, UFS, &
                                     QVAR, QRHO, QU, QV, QFA, QFS, &
                                     nadv

      implicit none

      double precision, parameter:: small = 1.d-8

      integer lo(2), hi(2)
      integer uin_l1,uin_l2,uin_h1,uin_h2
      integer q_l1,q_l2,q_h1,q_h2
      integer src_l1,src_l2,src_h1,src_h2

      double precision :: uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
      double precision :: q(q_l1:q_h1,q_l2:q_h2,QVAR)
      double precision :: src (src_l1:src_h1,src_l2:src_h2,NVAR)
      double precision :: srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
      double precision :: dx, dy, dt

      integer          :: i, j
      integer          :: ngp, ngf, loq(2), hiq(2)
      integer          :: iadv, ispec, n, nq
      double precision :: courx, coury, courmx, courmy

      do i=1,2
         loq(i) = lo(i)-ngp
         hiq(i) = hi(i)+ngp
      enddo

      do j = loq(2),hiq(2)
         do i = loq(1),hiq(1)
            q(i,j,QRHO) = uin(i,j,URHO)
            q(i,j,QU  ) = uin(i,j,UX  )
            q(i,j,QV  ) = uin(i,j,UY  )
         enddo
      enddo

!    Load advected quatities, c, into q, assuming they arrived in uin as rho.c
     do iadv = 1, nadv
         n  = UFA + iadv - 1
         nq = QFA + iadv - 1
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
               q(i,j,nq) = uin(i,j,n)/q(i,j,QRHO)
            enddo
         enddo
      enddo
      
!     Load chemical species, c, into q, assuming they arrived in uin as rho.c
      do ispec = 1, nspec
         n  = UFS + ispec - 1
         nq = QFS + ispec - 1
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
               q(i,j,nq) = uin(i,j,n)/q(i,j,QRHO)
            enddo
         enddo
      enddo
      
!     Compute sources in terms of Q
      do j = lo(2)-1, hi(2)+1
        do i = lo(1)-1, hi(1)+1

           do ispec = 1,nspec
              srcQ(i,j,QFS+ispec-1) = src(i,j,UFS+ispec-1)/q(i,j,QRHO)
           enddo

           do iadv = 1,nadv
              srcQ(i,j,QFA+iadv-1) = src(i,j,UFA+iadv-1)/q(i,j,QRHO)
           enddo

        end do
      end do

      end subroutine ctoprim

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine trace(q,qd_l1,qd_l2,qd_h1,qd_h2, &
                       dq,qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                       ilo1,ilo2,ihi1,ihi2,dx,dy,dt)

      use network           , only : nspec
      use meth_params_module, only : iorder, QVAR, QRHO, QFA, QFS, nadv

      implicit none

      integer ilo1,ilo2,ihi1,ihi2
      integer qd_l1,qd_l2,qd_h1,qd_h2
      integer qpd_l1,qpd_l2,qpd_h1,qpd_h2

      double precision dx, dy, dt
      double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision  dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)

      ! Local variables
      integer i, j
      integer n, iadv, ispec

      ! Default dq=0.
      dq(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:QVAR) = 0.d0

      ! Compute slopes in x-direction
      if (iorder .ne. 1) then

         call uslope( q, qd_l1, qd_l2, qd_h1, qd_h2, &
                     dq, qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                     ilo1,ilo2,ihi1,ihi2,QVAR,1)

      endif

      n = QRHO
      call trace_x(n,q,qd_l1,qd_l2,qd_h1,qd_h2, &
                   dq,qxm,qxp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                   ilo1,ilo2,ihi1,ihi2,dx,dt)

      do ispec = 1, nspec
         n = QFS + ispec - 1
         call trace_x(n,q,qd_l1,qd_l2,qd_h1,qd_h2, &
                      dq,qxm,qxp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                      ilo1,ilo2,ihi1,ihi2,dx,dt)
      end do

      do iadv = 1, nadv
         n = QFA + iadv - 1
         call trace_x(n,q,qd_l1,qd_l2,qd_h1,qd_h2, &
                      dq,qxm,qxp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                      ilo1,ilo2,ihi1,ihi2,dx,dt)
      end do

!     ------------------------------------------------------------------

      ! Compute slopes in y-direction
      if (iorder .ne. 1) then

         call uslope(q,qd_l1,qd_l2,qd_h1,qd_h2, &
                     dq,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                     ilo1,ilo2,ihi1,ihi2,QVAR,2)

      endif

      n = QRHO
      call trace_y(n,q,qd_l1,qd_l2,qd_h1,qd_h2, &
                   dq,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                   ilo1,ilo2,ihi1,ihi2,dy,dt)

      do ispec = 1, nspec
         n = QFS + ispec - 1
         call trace_y(n,q,qd_l1,qd_l2,qd_h1,qd_h2, &
                        dq,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                        ilo1,ilo2,ihi1,ihi2,dy,dt)
      end do

      do iadv = 1, nadv
         n = QFA + iadv - 1
         call trace_y(n,q,qd_l1,qd_l2,qd_h1,qd_h2, &
                      dq,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                      ilo1,ilo2,ihi1,ihi2,dy,dt)
      end do

      end subroutine trace

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine trace_x(n,q,qd_l1,qd_l2,qd_h1,qd_h2, &
                         dq,qxm,qxp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                         ilo1,ilo2,ihi1,ihi2,dx,dt)

      use network           , only : nspec
      use meth_params_module, only : iorder, QVAR, QRHO, QU, QFA, QFS, nadv

      implicit none

      integer n
      integer ilo1,ilo2,ihi1,ihi2
      integer qd_l1,qd_l2,qd_h1,qd_h2
      integer qpd_l1,qpd_l2,qpd_h1,qpd_h2

      double precision :: dx, dt
      double precision ::   q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision ::  dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision :: qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision :: qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)

      ! Local variables
      integer          :: i, j
      double precision :: dtdx
      double precision :: u, spzero, acmprght, acmpleft

      dtdx = dt/dx

      do j = ilo2-1, ihi2+1

         ! Right state
         do i = ilo1, ihi1+1
            u = q(i,j,QU)
            if (u .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = u*dtdx
            endif
            acmprght = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
            qxp(i,j,n) = q(i,j,n) + acmprght
         enddo

         ! Left state
         do i = ilo1-1, ihi1
            u = q(i,j,QU)
            if (u .ge. 0.d0) then
               spzero = u*dtdx
            else
               spzero = 1.d0
            endif
            acmpleft = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
            qxm(i+1,j,n) = q(i,j,n) + acmpleft
         enddo

      enddo

      end subroutine trace_x

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine trace_y(n,q,qd_l1,qd_l2,qd_h1,qd_h2, &
                         dq,qym,qyp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                         ilo1,ilo2,ihi1,ihi2,dy,dt)

      use network, only : nspec
      use meth_params_module, only : iorder, QVAR, QRHO, QV, QFA, QFS, nadv

      implicit none

      integer n
      integer ilo1,ilo2,ihi1,ihi2
      integer qd_l1,qd_l2,qd_h1,qd_h2
      integer qpd_l1,qpd_l2,qpd_h1,qpd_h2

      double precision ::  dy, dt
      double precision ::   q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision ::  dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision :: qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision :: qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)

      ! Local variables
      integer          :: i, j
      double precision :: v, spzero, acmptop, acmpbot
      double precision :: dtdy

      dtdy = dt/dy

      do i = ilo1-1, ihi1+1

         ! Top state
         do j = ilo2, ihi2+1
            v = q(i,j,QV)
            if (v .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = v*dtdy
            endif
            acmptop = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
            qyp(i,j,n) = q(i,j,n) + acmptop
         enddo

         ! Bottom state
         do j = ilo2-1, ihi2
            v = q(i,j,QV)
            if (v .ge. 0.d0) then
               spzero = v*dtdy
            else
               spzero = 1.d0
            endif
            acmpbot = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
            qym(i,j+1,n) = q(i,j,n) + acmpbot
         enddo

      enddo

      end subroutine trace_y

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine consup( uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                        uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                        src , src_l1, src_l2, src_h1, src_h2, &
                        flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                        flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                        div,lo,hi,dx,dt)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UX, UY, UFA, UFS, normalize_species, nadv

      implicit none

      integer lo(2), hi(2)
      integer uin_l1,uin_l2,uin_h1,uin_h2
      integer uout_l1,uout_l2,uout_h1,uout_h2
      integer   src_l1,  src_l2,  src_h1,  src_h2
      integer flux1_l1,flux1_l2,flux1_h1,flux1_h2
      integer flux2_l1,flux2_l2,flux2_h1,flux2_h2

      double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
      double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
      double precision   src(  src_l1:  src_h1,  src_l2:  src_h2,NVAR)
      double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
      double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
      double precision div(lo(1):hi(1)+1,lo(2):hi(2)+1)
      double precision dx(2), dt

      integer :: i, j, n
      integer :: ispec, iadv

      double precision div1
      double precision SrU, SrV
      double precision rho, Up, Vp, SrE

      ! Normalize the species fluxes
      if (normalize_species .eq. 1) &
         call normalize_species_fluxes( &
                  flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                  flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                  lo,hi)

      n = URHO
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            uout(i,j,n) = uin(i,j,n) + dt * &
                   ( (flux1(i,j,n) - flux1(i+1,j,n)) / dx(1) &
                 +   (flux2(i,j,n) - flux2(i,j+1,n)) / dx(2) ) &
                 +   dt * src(i,j,n)
         enddo
      enddo

      do iadv = 1, nadv
         n = UFA + iadv - 1
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               uout(i,j,n) = uin(i,j,n) + dt * &
                      ( (flux1(i,j,n) - flux1(i+1,j,n)) / dx(1) &
                    +   (flux2(i,j,n) - flux2(i,j+1,n)) / dx(2) ) &
                    +   dt * src(i,j,n)
            enddo
         enddo
      enddo

      do ispec = 1, nspec
         n = UFS + ispec -1
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               uout(i,j,n) = uin(i,j,n) + dt * &
                      ( (flux1(i,j,n) - flux1(i+1,j,n)) / dx(1) &
                    +   (flux2(i,j,n) - flux2(i,j+1,n)) / dx(2) ) &
                    +   dt * src(i,j,n)
            enddo
         enddo
      enddo

      ! Scale by face area in order to correctly reflux
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)+1
        flux1(i,j,1:NVAR) = dt * flux1(i,j,1:NVAR) * dx(2)
        flux1(i,j,UX:UY)  = 0.d0
      enddo
      enddo

      ! Scale by face area in order to correctly reflux
      do j = lo(2),hi(2)+1 
      do i = lo(1),hi(1)
        flux2(i,j,1:NVAR) = dt * flux2(i,j,1:NVAR) * dx(1)
        flux2(i,j,UX:UY)  = 0.d0
      enddo
      enddo

      end subroutine consup

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine upwind(ql, qr, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                        flx, flx_l1, flx_l2, flx_h1, flx_h2, &
                        ugd, ugd_l1, ugd_l2, ugd_h1, ugd_h2, &
                        idir, ilo1, ihi1, ilo2, ihi2)

      use network, only : nspec
      use prob_params_module, only : physbc_lo,physbc_hi,Symmetry
      use meth_params_module, only : QVAR, NVAR, QRHO, QFA, QFS, &
                                     URHO, UFA, UFS, nadv

      implicit none

      double precision, parameter:: small = 1.d-8

      integer qpd_l1, qpd_l2, qpd_h1, qpd_h2
      integer gd_l1, gd_l2, gd_h1, gd_h2
      integer flx_l1, flx_l2, flx_h1, flx_h2
      integer ugd_l1, ugd_l2, ugd_h1, ugd_h2
      integer idir, ilo1, ihi1, ilo2, ihi2
      integer ilo,ihi,jlo,jhi

      double precision ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision flx(flx_l1:flx_h1,flx_l2:flx_h2,NVAR)
      double precision ugd(ugd_l1:ugd_h1,ugd_l2:ugd_h2)
      double precision ustar, qavg

      integer iadv, ispec, n, nq
      integer i, j

!************************************************************
!  set min/max based on normal direction
      if(idir.eq.1) then
         ilo = ilo1
         ihi = ihi1 + 1
         jlo = ilo2
         jhi = ihi2
      else
         ilo = ilo1
         ihi = ihi1
         jlo = ilo2
         jhi = ihi2+1
      endif

      do j = jlo, jhi
         do i = ilo, ihi

            ustar = ugd(i,j)

            ! Density
            if (ustar .gt. 0.d0) then
               flx(i,j,URHO) = ql(i,j,QRHO) * ugd(i,j)
            else if (ustar .lt. 0.d0) then
               flx(i,j,URHO) = qr(i,j,QRHO) * ugd(i,j)
            else 
               qavg = 0.5d0 * (ql(i,j,QRHO) + qr(i,j,QRHO))
               flx(i,j,URHO) = qavg * ugd(i,j)
            endif

            do iadv = 1, nadv
               n = UFA + iadv - 1
               nq = QFA + iadv - 1
               if (ustar .gt. 0.d0) then
                  flx(i,j,n) = ql(i,j,nq)
               else if (ustar .lt. 0.d0) then
                  flx(i,j,n) = qr(i,j,nq)
               else 
                  qavg = 0.5d0 * (ql(i,j,nq) + qr(i,j,nq))
                  flx(i,j,n) = qavg
               endif
               flx(i,j,n) = flx(i,j,URHO)*flx(i,j,n)
            enddo

            do ispec = 1, nspec
               n  = UFS + ispec - 1
               nq = QFS + ispec - 1
               if (ustar .gt. 0.d0) then
                  flx(i,j,n) = ql(i,j,nq)
               else if (ustar .lt. 0.d0) then
                  flx(i,j,n) = qr(i,j,nq)
               else 
                  qavg = 0.5d0 * (ql(i,j,nq) + qr(i,j,nq))
                  flx(i,j,n) = qavg
               endif
               flx(i,j,n) = flx(i,j,URHO)*flx(i,j,n)
            enddo

         enddo
      enddo

      end subroutine upwind

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine uslope(q,qd_l1,qd_l2,qd_h1,qd_h2, &
                        dq,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                        ilo1,ilo2,ihi1,ihi2,nv,idir)

      implicit none

      integer ilo,ihi
      integer qd_l1,qd_l2,qd_h1,qd_h2
      integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
      integer ilo1,ilo2,ihi1,ihi2,nv,idir

      double precision     q( qd_l1: qd_h1, qd_l2: qd_h2,nv)
      double precision    dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,nv)

!     Local arrays
      double precision, allocatable::dsgn(:),dlim(:),df(:),dcen(:)

      integer i, j, n
      double precision dlft, drgt, slop, dq1
      double precision four3rd, sixth

      four3rd = 4.d0/3.d0
      sixth = 1.d0/6.d0

      ilo = MIN(ilo1,ilo2)
      ihi = MAX(ihi1,ihi2)

      allocate (dsgn(ilo-2:ihi+2))
      allocate (dlim(ilo-2:ihi+2))
      allocate (  df(ilo-2:ihi+2))
      allocate (dcen(ilo-2:ihi+2))

      do n = 1, nv 
          if (idir .eq. 1) then

             ! slopes in first coordinate direction
             do j = ilo2-1, ihi2+1

                ! first compute Fromm slopes
                do i = ilo1-2, ihi1+2
                      dlft = 2.d0*(q(i  ,j,n) - q(i-1,j,n))
                      drgt = 2.d0*(q(i+1,j,n) - q(i  ,j,n))
                      dcen(i) = .25d0 * (dlft+drgt)
                      dsgn(i) = sign(1.d0, dcen(i))
                      slop = min( abs(dlft), abs(drgt) )
!                      dlim(i) = cvmgp( slop, 0.d0, dlft*drgt )
                      if (dlft*drgt .ge. 0.d0) then
                         dlim(i) = slop
                      else
                         dlim(i) = 0.d0
                      endif
                      df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
                  enddo

!                 Now limited fourth order slopes
                  do i = ilo1-1, ihi1+1
                      dq1 = four3rd*dcen(i) - sixth*(df(i+1) + df(i-1))
                      dq(i,j,n) = dsgn(i)*min(dlim(i),abs(dq1))
                  enddo
              enddo

          else


!            Compute slopes in second coordinate direction
             do i = ilo1-1, ihi1+1

!               First compute Fromm slopes for this column
                do j = ilo2-2, ihi2+2
                      dlft = 2.d0*(q(i,j  ,n) - q(i,j-1,n))
                      drgt = 2.d0*(q(i,j+1,n) - q(i,j  ,n))
                      dcen(j) = .25d0 * (dlft+drgt)
                      dsgn(j) = sign( 1.d0, dcen(j) )
                      slop = min( abs(dlft), abs(drgt) )
                      if (dlft*drgt .ge. 0.d0) then
                         dlim(j) = slop
                      else
                         dlim(j) = 0.d0
                      endif
                      df(j) = dsgn(j)*min( dlim(j),abs(dcen(j)) )
                  enddo

!                 Now compute limited fourth order slopes
                  do j = ilo2-1, ihi2+1
                      dq1 = four3rd*dcen(j) - &
                           sixth*( df(j+1) + df(j-1) )
                      dq(i,j,n) = dsgn(j)*min(dlim(j),abs(dq1))
                  enddo
              enddo

          endif
      enddo

      deallocate(dsgn,dlim,df,dcen)

      end subroutine uslope

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine divu(lo,hi,q,q_l1,q_l2,q_h1,q_h2,delta, &
                      div,div_l1,div_l2,div_h1,div_h2)

      use meth_params_module, only : QU, QV

      implicit none

      integer          :: lo(2),hi(2)
      integer          :: q_l1,q_l2,q_h1,q_h2
      integer          :: div_l1,div_l2,div_h1,div_h2
      double precision :: q(q_l1:q_h1,q_l2:q_h2,*)
      double precision :: div(div_l1:div_h1,div_l2:div_h2)
      double precision :: delta(2)

      integer          :: i, j
      double precision :: ux,vy

      do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)+1
         ux = 0.5d0*(q(i,j,QU)-q(i-1,j,QU)+q(i,j-1,QU)-q(i-1,j-1,QU))/delta(1)
         vy = 0.5d0*(q(i,j,QV)-q(i,j-1,QV)+q(i-1,j,QV)-q(i-1,j-1,QV))/delta(2)
         div(i,j) = ux + vy
      enddo
      enddo

      end subroutine divu

