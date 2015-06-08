
! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine advect(time,lo,hi,domlo,domhi, &
                        uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                        uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                        ugdx,ugdx_l1,ugdx_l2,ugdx_l3,ugdx_h1,ugdx_h2,ugdx_h3, &
                        ugdy,ugdy_l1,ugdy_l2,ugdy_l3,ugdy_h1,ugdy_h2,ugdy_h3, &
                        ugdz,ugdz_l1,ugdz_l2,ugdz_l3,ugdz_h1,ugdz_h2,ugdz_h3, &
                        src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                        flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                        flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                        flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                        delta,dt,verbose)

      use meth_params_module, only : QVAR, NVAR, NHYP, normalize_species

      implicit none

      integer lo(3),hi(3),verbose
      integer domlo(3),domhi(3)
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
      integer ugdx_l1,ugdx_l2,ugdx_l3,ugdx_h1,ugdx_h2,ugdx_h3
      integer ugdy_l1,ugdy_l2,ugdy_l3,ugdy_h1,ugdy_h2,ugdy_h3
      integer ugdz_l1,ugdz_l2,ugdz_l3,ugdz_h1,ugdz_h2,ugdz_h3
      integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
      integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
      integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3

      double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
      double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
      double precision ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2,ugdx_l3:ugdx_h3)
      double precision ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2,ugdy_l3:ugdy_h3)
      double precision ugdz(ugdz_l1:ugdz_h1,ugdz_l2:ugdz_h2,ugdz_l3:ugdz_h3)
      double precision src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
      double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
      double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
      double precision flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
      double precision delta(3),dt,time

!     Automatic arrays for workspace
      double precision, allocatable:: q(:,:,:,:)
      double precision, allocatable:: div(:,:,:)
      double precision, allocatable:: srcQ(:,:,:,:)

      integer ngq,ngf
!     integer i_c,j_c,k_c

      double precision dx,dy,dz

      allocate(     q(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,QVAR))
      allocate(  srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR))
      allocate(   div(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1))

      dx = delta(1)
      dy = delta(2)
      dz = delta(3)

      ngq = NHYP
      ngf = 1

!     Translate to primitive variables, compute sound speeds
      call ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                   q,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                   src,srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                   dx,dy,dz,dt,ngq,ngf)

!     Compute hyperbolic fluxes using unsplit Godunov
      call umeth3d(q,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                   srcQ, src_l1, src_l2, src_l3, src_h1, src_h2,src_h3, &
                   lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),dx,dy,dz,dt, &
                   flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                   flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                   flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                   ugdx,ugdx_l1,ugdx_l2,ugdx_l3,ugdx_h1,ugdx_h2,ugdx_h3, &
                   ugdy,ugdy_l1,ugdy_l2,ugdy_l3,ugdy_h1,ugdy_h2,ugdy_h3, &
                   ugdz,ugdz_l1,ugdz_l2,ugdz_l3,ugdz_h1,ugdz_h2,ugdz_h3)

!     Compute divergence of velocity field (on surroundingNodes(lo,hi))
      call divu(lo,hi,q,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                delta,div,lo(1),lo(2),lo(3),hi(1)+1,hi(2)+1,hi(3)+1)

!     Conservative update
      call consup(uin,    uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2, uin_h3, &
                  uout,  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3, &
                  src,    src_l1,  src_l2,  src_l3,  src_h1,  src_h2, src_h3, &
                  flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2, flux1_h3, &
                  flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2, flux2_h3, &
                  flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2, flux3_h3, &
                  div,lo,hi,delta,dt)

      ! Enforce the species >= 0
      call enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3,lo,hi)

      ! Normalize the species 
      if (normalize_species .eq. 1) &
         call normalize_new_species(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3,lo,hi)

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

      subroutine umeth3d(q, qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3, &
                         srcQ, src_l1, src_l2, src_l3, src_h1, src_h2, src_h3,  &
                         ilo1, ilo2, ilo3, ihi1, ihi2, ihi3, dx, dy, dz, dt, &
                         flux1, fd1_l1, fd1_l2, fd1_l3, fd1_h1, fd1_h2, fd1_h3, &
                         flux2, fd2_l1, fd2_l2, fd2_l3, fd2_h1, fd2_h2, fd2_h3, &
                         flux3, fd3_l1, fd3_l2, fd3_l3, fd3_h1, fd3_h2, fd3_h3, &
                         ugdx,ugdx_l1,ugdx_l2,ugdx_l3,ugdx_h1,ugdx_h2,ugdx_h3, &
                         ugdy,ugdy_l1,ugdy_l2,ugdy_l3,ugdy_h1,ugdy_h2,ugdy_h3, &
                         ugdz,ugdz_l1,ugdz_l2,ugdz_l3,ugdz_h1,ugdz_h2,ugdz_h3)

      use network, only : nspec
      use meth_params_module, only : QU, QV, QW, QVAR, NVAR

      implicit none

      integer qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3
      integer src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
      integer fd1_l1, fd1_l2, fd1_l3, fd1_h1, fd1_h2, fd1_h3
      integer fd2_l1, fd2_l2, fd2_l3, fd2_h1, fd2_h2, fd2_h3
      integer fd3_l1, fd3_l2, fd3_l3, fd3_h1, fd3_h2, fd3_h3
      integer ugdx_l1,ugdx_l2,ugdx_l3, ugdx_h1,ugdx_h2, ugdx_h3
      integer ugdy_l1,ugdy_l2,ugdy_l3, ugdy_h1,ugdy_h2, ugdy_h3
      integer ugdz_l1,ugdz_l2,ugdz_l3, ugdz_h1,ugdz_h2, ugdz_h3
      integer ilo1, ilo2, ilo3, ihi1, ihi2, ihi3

      double precision dx, dy, dz, dt
      double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision  srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3)
      double precision ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2,ugdx_l3:ugdx_h3)
      double precision ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2,ugdy_l3:ugdy_h3)
      double precision ugdz(ugdz_l1:ugdz_h1,ugdz_l2:ugdz_h2,ugdz_l3:ugdz_h3)
      double precision flux1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,fd1_l3:fd1_h3,NVAR)
      double precision flux2(fd2_l1:fd2_h1,fd2_l2:fd2_h2,fd2_l3:fd2_h3,NVAR)
      double precision flux3(fd3_l1:fd3_h1,fd3_l2:fd3_h2,fd3_l3:fd3_h3,NVAR)

!     Left and right state arrays (edge centered, cell centered)
      double precision, allocatable:: dq(:,:,:,:),  qm(:,:,:,:),   qp(:,:,:,:)
      double precision, allocatable::qxm(:,:,:,:),qym(:,:,:,:), qzm(:,:,:,:)
      double precision, allocatable::qxp(:,:,:,:),qyp(:,:,:,:), qzp(:,:,:,:)

!     Work arrays to hold 3 planes of riemann state and conservative fluxes
      double precision, allocatable::   fx(:,:,:,:),  fy(:,:,:,:), fz(:,:,:,:)

!     Local scalar variables
      double precision :: dtdx
      double precision :: hdtdx, hdt, hdtdy, hdtdz
      integer          :: i,j,k

      allocate (  dq(ilo1-1:ihi1+2,ilo2-1:ihi2+2,ilo3-1:ihi3+2,QVAR) )
      allocate (  qm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,ilo3-1:ihi3+2,QVAR) )
      allocate (  qp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,ilo3-1:ihi3+2,QVAR) )
      allocate ( qxm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,ilo3-1:ihi3+2,QVAR) )
      allocate ( qxp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,ilo3-1:ihi3+2,QVAR) )
      allocate ( qym(ilo1-1:ihi1+2,ilo2-1:ihi2+2,ilo3-1:ihi3+2,QVAR) )
      allocate ( qyp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,ilo3-1:ihi3+2,QVAR) )
      allocate ( qzm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,ilo3-1:ihi3+2,QVAR) )
      allocate ( qzp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,ilo3-1:ihi3+2,QVAR) )
      allocate (  fx(ilo1  :ihi1+1,ilo2-1:ihi2+1,ilo3-1:ihi3+1,NVAR) )
      allocate (  fy(ilo1-1:ihi1+1,ilo2  :ihi2+1,ilo3-1:ihi3+1,NVAR) )
      allocate (  fz(ilo1-1:ihi1+1,ilo2-1:ihi2+1,ilo3  :ihi3+1,NVAR) )

!     Local constants
      dtdx = dt/dx
      hdtdx = 0.5d0*dtdx
      hdtdy = 0.5d0*dt/dy
      hdtdz = 0.5d0*dt/dz
      hdt = 0.5d0*dt

      do k = ugdx_l3,ugdx_h3
      do j = ugdx_l2,ugdx_h2
      do i = ugdx_l1,ugdx_h1
         ugdx(i,j,k) = 0.5d0 * (q(i,j,k,QU) + q(i-1,j,k,QU))
      end do
      end do
      end do

      do k = ugdy_l3,ugdy_h3
      do j = ugdy_l2,ugdy_h2
      do i = ugdy_l1,ugdy_h1
         ugdy(i,j,k) = 0.5d0 * (q(i,j,k,QV) + q(i,j-1,k,QV))
      end do
      end do
      end do

      do k = ugdz_l3,ugdz_h3
      do j = ugdz_l2,ugdz_h2
      do i = ugdz_l1,ugdz_h1
         ugdz(i,j,k) = 0.5d0 * (q(i,j,k,QW) + q(i,j,k-1,QW))
      end do
      end do
      end do

!     Trace to edges w/o transverse flux correction terms
      call trace(q,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                 dq,qxm,qxp,qym,qyp,qzm,qzp,ilo1-1,ilo2-1,ilo3-1,ihi1+2,ihi2+2,ihi3+2, &
                 ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dx,dy,dz,dt)

!     Upwind on x-edges
      call upwind(qxm, qxp, ilo1-1, ilo2-1, ilo3-1, ihi1+2, ihi2+2, ihi3+2, &
                  fx,       ilo1  , ilo2-1, ilo3-1, ihi1+1, ihi2+1, ihi3+1, &
                  ugdx, ugdx_l1, ugdx_l2, ugdx_l3, ugdx_h1, ugdx_h2, ugdx_h3, &
                  1, ilo1, ihi1, ilo2-1, ihi2+1, ilo3-1, ihi3+1)

!     Upwind on y-edges
      call upwind(qym, qyp, ilo1-1, ilo2-1, ilo3-1, ihi1+2, ihi2+2, ihi3+2, &
                  fy,       ilo1-1, ilo2  , ilo3-1, ihi1+1, ihi2+1, ihi3+1, &
                  ugdy, ugdy_l1, ugdy_l2, ugdy_l3, ugdy_h1, ugdy_h2, ugdy_h3, &
                  2, ilo1-1, ihi1+1, ilo2, ihi2, ilo3-1, ihi3+1)

!     Upwind on z-edges
      call upwind(qzm, qzp, ilo1-1, ilo2-1, ilo3-1, ihi1+2, ihi2+2, ihi3+2, &
                  fz,       ilo1-1, ilo2-1, ilo3  , ihi1+1, ihi2+1, ihi3+1, &
                  ugdz, ugdz_l1, ugdz_l2, ugdz_l3, ugdz_h1, ugdz_h2, ugdz_h3, &
                  3, ilo1-1, ihi1+1, ilo2-1, ihi2+1, ilo3, ihi3)

!     Use edge states to create transverse derivatives in y- and z-directions
      call transyz(qxm, qm, qxp, qp, ilo1-1, ilo2-1, ilo3-1, ihi1+2, ihi2+2, ihi3+2, &
                   fy,               ilo1-1, ilo2  , ilo3-1, ihi1+1, ihi2+1, ihi3+1, &
                   fz,               ilo1-1, ilo2-1, ilo3  , ihi1+1, ihi2+1, ihi3+1, &
                   srcQ,             src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                   hdt, hdtdy, hdtdz,ilo1-1, ihi1+1, ilo2, ihi2, ilo3, ihi3)

!     Upwind on x-edges to create final fluxes
      call upwind(qm, qp, ilo1-1, ilo2-1, ilo3-1, ihi1+2, ihi2+2, ihi3+2,&
                  flux1, fd1_l1, fd1_l2, fd1_l3, fd1_h1, fd1_h2, fd1_h3, &
                  ugdx, ugdx_l1, ugdx_l2, ugdx_l3, ugdx_h1, ugdx_h2, ugdx_h3, &
                  1, ilo1, ihi1, ilo2, ihi2, ilo3, ihi3)
      
!     Use edge states to create transverse derivative in x- and z-directions
      call transxz(qym, qm, qyp,qp, ilo1-1, ilo2-1, ilo3-1, ihi1+2, ihi2+2, ihi3+2, &
                   fx,              ilo1  , ilo2-1, ilo3-1, ihi1+1, ihi2+1, ihi3+1, &
                   fz,              ilo1-1, ilo2-1, ilo3  , ihi1+1, ihi2+1, ihi3+1, &
                   srcQ,            src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                   hdt,hdtdx,hdtdz, ilo1,   ihi1,   ilo2-1, ihi2+1, ilo3, ihi3)

!     Upwind on y-edges to create final fluxes
      call upwind(qm, qp, ilo1-1, ilo2-1, ilo3-1, ihi1+2, ihi2+2, ihi3+2, &
                  flux2,  fd2_l1, fd2_l2, fd2_l3, fd2_h1, fd2_h2, fd2_h3, &
                  ugdy, ugdy_l1, ugdy_l2, ugdy_l3, ugdy_h1, ugdy_h2, ugdy_h3, &
                  2, ilo1, ihi1, ilo2, ihi2, ilo3, ihi3)
      
!     Use edge states to create transverse derivative in x- and y-directions
      call transxy(qzm,qm,qzp,qp,ilo1-1, ilo2-1, ilo3-1, ihi1+2, ihi2+2, ihi3+2, &
                   fx,           ilo1  , ilo2-1, ilo3-1, ihi1+1, ihi2+1, ihi3+1, &
                   fy,           ilo1-1, ilo2  , ilo3-1, ihi1+1, ihi2+1, ihi3+1, &
                   srcQ,         src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                   hdt, hdtdz, ilo1, ihi1, ilo2, ihi2, ilo3-1, ihi3+1)

!     Upwind on z-edges to create final fluxes
      call upwind(qm, qp, ilo1-1, ilo2-1, ilo3-1, ihi1+2, ihi2+2, ihi3+2, &
                  flux3, fd3_l1, fd3_l2, fd3_l3, fd3_h1, fd3_h2, fd3_h3, &
                  ugdz, ugdz_l1, ugdz_l2, ugdz_l3, ugdz_h1, ugdz_h2, ugdz_h3, &
                  3, ilo1, ihi1, ilo2, ihi2, ilo3, ihi3)
      
      deallocate(dq,qm,qp,qxm,qxp,qym,qyp,qzm,qzp)
      deallocate(fx,fy,fz)

      end subroutine umeth3d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ctoprim(lo,hi, &
                         uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                         q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3, &
                         src,srcQ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                         dx,dy,dz,dt,ngp,ngf)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UX, UY, UZ, UFA, UFS, &
                                     QVAR, QRHO, QU, QV, QW, QFA, QFS, &
                                     nadv

      implicit none

      double precision, parameter:: small = 1.d-8

      integer lo(3), hi(3)
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3

      double precision :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
      double precision :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
      double precision :: src (src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,NVAR)
      double precision :: srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      double precision :: dx, dy, dz, dt

      integer          :: i, j, k
      integer          :: ngp, ngf, loq(3), hiq(3)
      integer          :: iadv, ispec, n, nq

      do i=1,3
         loq(i) = lo(i)-ngp
         hiq(i) = hi(i)+ngp
      enddo

      do k = loq(3),hiq(3)
      do j = loq(2),hiq(2)
      do i = loq(1),hiq(1)
         q(i,j,k,QRHO) = uin(i,j,k,URHO)
         q(i,j,k,QU  ) = uin(i,j,k,UX  )
         q(i,j,k,QV  ) = uin(i,j,k,UY  )
         q(i,j,k,QW  ) = uin(i,j,k,UZ  )
      enddo
      enddo
      enddo

!    Load advected quatities, c, into q, assuming they arrived in uin as rho.c
     do iadv = 1, nadv
         n  = UFA + iadv - 1
         nq = QFA + iadv - 1
         do k = loq(3),hiq(3)
         do j = loq(2),hiq(2)
         do i = loq(1),hiq(1)
            q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
         enddo
         enddo
         enddo
      enddo
      
!     Load chemical species, c, into q, assuming they arrived in uin as rho.c
      do ispec = 1, nspec
         n  = UFS + ispec - 1
         nq = QFS + ispec - 1
         do k = loq(3),hiq(3)
         do j = loq(2),hiq(2)
         do i = loq(1),hiq(1)
            q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
         enddo
         enddo
         enddo
      enddo
      
!     Compute sources in terms of Q
      do k = lo(3)-1, hi(3)+1
      do j = lo(2)-1, hi(2)+1
      do i = lo(1)-1, hi(1)+1

        do ispec = 1,nspec
           srcQ(i,j,k,QFS+ispec-1) = src(i,j,k,UFS+ispec-1)/q(i,j,k,QRHO)
        enddo

        do iadv = 1,nadv
           srcQ(i,j,k,QFA+iadv-1) = src(i,j,k,UFA+iadv-1)/q(i,j,k,QRHO)
        enddo

      end do
      end do
      end do

     end subroutine ctoprim

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine consup( uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3,&
                        uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3,&
                        src , src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                        flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                        flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                        flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                        div,lo,hi,dx,dt)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UX, UZ, UFA, UFS, normalize_species, nadv

      implicit none

      integer lo(3), hi(3)
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
      integer   src_l1,  src_l2,  src_l3, src_h1,  src_h2, src_h3
      integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
      integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
      integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3

      double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l1:uin_h1,NVAR)
      double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
      double precision   src(  src_l1:  src_h1,  src_l2:  src_h2,src_l3:src_h3,NVAR)
      double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
      double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
      double precision flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
      double precision div(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
      double precision dx(3), dt

      integer :: i, j, k, n
      integer :: ispec, iadv

      ! Normalize the species fluxes
      if (normalize_species .eq. 1) &
         call normalize_species_fluxes( &
                  flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                  flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                  flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                  lo,hi)

      n = URHO
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            uout(i,j,k,n) = uin(i,j,k,n) + dt * &
                   ( (flux1(i,j,k,n) - flux1(i+1,j,k,n)) / dx(1) &
                 +   (flux2(i,j,k,n) - flux2(i,j+1,k,n)) / dx(2) &
                 +   (flux3(i,j,k,n) - flux3(i,j,k+1,n)) / dx(3) ) &
                 +   dt * src(i,j,k,n)
         enddo
      enddo
      enddo

      do iadv = 1, nadv
         n = UFA + iadv - 1
         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               uout(i,j,k,n) = uin(i,j,k,n) + dt * &
                      ( (flux1(i,j,k,n) - flux1(i+1,j,k,n)) / dx(1) &
                    +   (flux2(i,j,k,n) - flux2(i,j+1,k,n)) / dx(2) &
                    +   (flux3(i,j,k,n) - flux3(i,j,k+1,n)) / dx(3) ) &
                    +   dt * src(i,j,k,n)
            enddo
         enddo
         enddo
      enddo

      do ispec = 1, nspec
         n = UFS + ispec -1
         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               uout(i,j,k,n) = uin(i,j,k,n) + dt * &
                      ( (flux1(i,j,k,n) - flux1(i+1,j,k,n)) / dx(1) &
                    +   (flux2(i,j,k,n) - flux2(i,j+1,k,n)) / dx(2) &
                    +   (flux3(i,j,k,n) - flux3(i,j,k+1,n)) / dx(3) ) &
                    +   dt * src(i,j,k,n)
            enddo
         enddo
         enddo
      enddo

      ! Scale by face area in order to correctly reflux
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)+1
        flux1(i,j,k,1:NVAR) = dt * flux1(i,j,k,1:NVAR) * dx(2) * dx(3)
        flux1(i,j,k,UX:UZ)  = 0.d0
      enddo
      enddo
      enddo

      ! Scale by face area in order to correctly reflux
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)+1 
      do i = lo(1),hi(1)
        flux2(i,j,k,1:NVAR) = dt * flux2(i,j,k,1:NVAR) * dx(1) * dx(3)
        flux2(i,j,k,UX:UZ)  = 0.d0
      enddo
      enddo
      enddo

      ! Scale by face area in order to correctly reflux
      do k = lo(3),hi(3)+1 
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        flux3(i,j,k,1:NVAR) = dt * flux3(i,j,k,1:NVAR) * dx(1) * dx(2)
        flux3(i,j,k,UX:UZ)  = 0.d0
      enddo
      enddo
      enddo

      end subroutine consup

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine upwind(ql, qr, qpd_l1, qpd_l2, qpd_l3, qpd_h1, qpd_h2, qpd_h3, &
                        flx, flx_l1, flx_l2, flx_l3, flx_h1, flx_h2, flx_h3, &
                        ugd, ugd_l1, ugd_l2, ugd_l3, ugd_h1, ugd_h2, ugd_h3, &
                        idir, ilo1, ihi1, ilo2, ihi2, ilo3, ihi3)

      use network, only : nspec
      use meth_params_module, only : QVAR, NVAR, QRHO, QFA, QFS, &
                                     URHO, UFA, UFS, nadv

      implicit none

      double precision, parameter:: small = 1.d-8

      integer qpd_l1, qpd_l2, qpd_l3, qpd_h1, qpd_h2, qpd_h3
      integer flx_l1, flx_l2, flx_l3, flx_h1, flx_h2, flx_h3
      integer ugd_l1, ugd_l2, ugd_l3, ugd_h1, ugd_h2, ugd_h3
      integer idir, ilo1, ihi1, ilo2, ihi2, ilo3, ihi3
      integer ilo,ihi,jlo,jhi,klo,khi

      double precision ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      double precision flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,NVAR)
      double precision ugd(ugd_l1:ugd_h1,ugd_l2:ugd_h2,ugd_l3:ugd_h3)
      double precision ustar, qavg

      integer iadv, ispec, n, nq
      integer i, j, k

!************************************************************
!  set min/max based on normal direction
      if (idir .eq. 1) then
         ilo = ilo1
         ihi = ihi1 + 1
         jlo = ilo2
         jhi = ihi2
         klo = ilo3
         khi = ihi3
      else if (idir .eq. 2) then
         ilo = ilo1
         ihi = ihi1
         jlo = ilo2
         jhi = ihi2+1
         klo = ilo3
         khi = ihi3
      else if (idir .eq. 3) then
         ilo = ilo1
         ihi = ihi1
         jlo = ilo2
         jhi = ihi2
         klo = ilo3
         khi = ihi3+1
      endif

      do k = klo, khi
      do j = jlo, jhi
         do i = ilo, ihi

            ustar = ugd(i,j,k)

            ! Density
            if (ustar .gt. 0.d0) then
               flx(i,j,k,URHO) = ql(i,j,k,QRHO) * ugd(i,j,k)
            else if (ustar .lt. 0.d0) then
               flx(i,j,k,URHO) = qr(i,j,k,QRHO) * ugd(i,j,k)
            else 
               qavg = 0.5d0 * (ql(i,j,k,QRHO) + qr(i,j,k,QRHO))
               flx(i,j,k,URHO) = qavg * ugd(i,j,k)
            endif

            do iadv = 1, nadv
               n = UFA + iadv - 1
               nq = QFA + iadv - 1
               if (ustar .gt. 0.d0) then
                  flx(i,j,k,n) = ql(i,j,k,nq)
               else if (ustar .lt. 0.d0) then
                  flx(i,j,k,n) = qr(i,j,k,nq)
               else 
                  qavg = 0.5d0 * (ql(i,j,k,nq) + qr(i,j,k,nq))
                  flx(i,j,k,n) = qavg
               endif
               flx(i,j,k,n) = flx(i,j,k,URHO)*flx(i,j,k,n)
            enddo

            do ispec = 1, nspec
               n  = UFS + ispec - 1
               nq = QFS + ispec - 1
               if (ustar .gt. 0.d0) then
                  flx(i,j,k,n) = ql(i,j,k,nq)
               else if (ustar .lt. 0.d0) then
                  flx(i,j,k,n) = qr(i,j,k,nq)
               else 
                  qavg = 0.5d0 * (ql(i,j,k,nq) + qr(i,j,k,nq))
                  flx(i,j,k,n) = qavg
               endif
               flx(i,j,k,n) = flx(i,j,k,URHO)*flx(i,j,k,n)
            enddo

         enddo
      enddo
      enddo

      end subroutine upwind

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine divu(lo,hi,q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,delta, &
                  div,div_l1,div_l2,div_l3,div_h1,div_h2,div_h3)
    
    use meth_params_module, only : QU, QV, QW
    use bl_constants_module
    
    implicit none

    integer          :: lo(3),hi(3)
    integer          :: q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
    integer          :: div_l1,div_l2,div_l3,div_h1,div_h2,div_h3
    double precision :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,*)
    double precision :: div(div_l1:div_h1,div_l2:div_h2,div_l3:div_h3)
    double precision :: delta(3)

    integer          :: i, j, k
    double precision :: ux, vy, wz

    do k=lo(3),hi(3)+1
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)+1
             
             ux = FOURTH*( &
                    + q(i  ,j  ,k  ,QU) - q(i-1,j  ,k  ,QU) &
                    + q(i  ,j  ,k-1,QU) - q(i-1,j  ,k-1,QU) &
                    + q(i  ,j-1,k  ,QU) - q(i-1,j-1,k  ,QU) &
                    + q(i  ,j-1,k-1,QU) - q(i-1,j-1,k-1,QU) )/ delta(1)

             vy = FOURTH*( &
                    + q(i  ,j  ,k  ,QV) - q(i  ,j-1,k  ,QV) &
                    + q(i  ,j  ,k-1,QV) - q(i  ,j-1,k-1,QV) &
                    + q(i-1,j  ,k  ,QV) - q(i-1,j-1,k  ,QV) &
                    + q(i-1,j  ,k-1,QV) - q(i-1,j-1,k-1,QV) )/ delta(2)

             wz = FOURTH*( &
                    + q(i  ,j  ,k  ,QW) - q(i  ,j  ,k-1,QW) &
                    + q(i  ,j-1,k  ,QW) - q(i  ,j-1,k-1,QW) &
                    + q(i-1,j  ,k  ,QW) - q(i-1,j  ,k-1,QW) &
                    + q(i-1,j-1,k  ,QW) - q(i-1,j-1,k-1,QW) )/ delta(3)

             div(i,j,k) = ux + vy + wz

          enddo
       enddo
    enddo
    
  end subroutine divu
