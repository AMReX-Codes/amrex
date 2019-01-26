! ::: ---------------------------------------------------------------
! ::: :: UMETH3D     Compute hyperbolic fluxes using unsplit second
! ::: ::               order Godunov integrator.
! ::: ::
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
! ::: :: c           => (const)  sound speed
! ::: :: csml        => (const)  local small c val
! ::: :: flatn       => (const)  flattening parameter
! ::: :: src         => (const)  source
! ::: :: nx          => (const)  number of cells in X direction
! ::: :: ny          => (const)  number of cells in Y direction
! ::: :: nz          => (const)  number of cells in Z direction
! ::: :: dx          => (const)  grid spacing in X direction
! ::: :: dy          => (const)  grid spacing in Y direction
! ::: :: dz          => (const)  grid spacing in Z direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux1      <=  (modify) flux in X direction on X edges
! ::: :: flux2      <=  (modify) flux in Y direction on Y edges
! ::: :: flux3      <=  (modify) flux in Z direction on Z edges
! L:: ----------------------------------------------------------------

      subroutine umeth3d(q, c, csml, flatn, qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3, &
                         srcQ, srcq_l1, srcq_l2, srcq_l3, srcq_h1, srcq_h2, srcq_h3, &
                         ilo1, ilo2, ilo3, ihi1, ihi2, ihi3, dx, dy, dz, dt, &
                         flux1, fd1_l1, fd1_l2, fd1_l3, fd1_h1, fd1_h2, fd1_h3, &
                         flux2, fd2_l1, fd2_l2, fd2_l3, fd2_h1, fd2_h2, fd2_h3, &
                         flux3, fd3_l1, fd3_l2, fd3_l3, fd3_h1, fd3_h2, fd3_h3, &
                         grav ,  gv_l1,  gv_l2,  gv_l3,  gv_h1,  gv_h2,  gv_h3, &
                         ugdnvx_out,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3, &
                         ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
                         ugdnvy_out,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3, &
                         ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
                         ugdnvz_out,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3, &
                         ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
                         divu_cc, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                         a_old,a_new,print_fortran_warnings)

      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
      use amrex_constants_module
      use meth_params_module, only : QVAR, NVAR, QU, ppm_type, &
                                     use_colglaz, corner_coupling, &
                                     version_2
      use slope_module
      use trace_ppm_module
      use trace_src_module
      use transverse_module
      use ppm_module

      implicit none

      integer qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3
      integer srcq_l1, srcq_l2, srcq_l3, srcq_h1, srcq_h2, srcq_h3
      integer ilo1, ilo2, ilo3, ihi1, ihi2, ihi3
      integer   d_l1,   d_l2,   d_l3,   d_h1,   d_h2,   d_h3
      integer fd1_l1, fd1_l2, fd1_l3, fd1_h1, fd1_h2, fd1_h3
      integer fd2_l1, fd2_l2, fd2_l3, fd2_h1, fd2_h2, fd2_h3
      integer fd3_l1, fd3_l2, fd3_l3, fd3_h1, fd3_h2, fd3_h3
      integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
      integer ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3
      integer ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3
      integer ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3
      integer km,kc,kt,k3d,n
      integer print_fortran_warnings
      integer i,j

      real(rt)     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      real(rt)     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      real(rt)  csml(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      real(rt) flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      real(rt)  srcQ(srcq_l1:srcq_h1,srcq_l2:srcq_h2,srcq_l3:srcq_h3,QVAR)
      real(rt) flux1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,fd1_l3:fd1_h3,NVAR)
      real(rt) flux2(fd2_l1:fd2_h1,fd2_l2:fd2_h2,fd2_l3:fd2_h3,NVAR)
      real(rt) flux3(fd3_l1:fd3_h1,fd3_l2:fd3_h2,fd3_l3:fd3_h3,NVAR)
      real(rt)  grav( gv_l1:gv_h1,  gv_l2:gv_h2,   gv_l3:gv_h3,   3)
      real(rt) ugdnvx_out(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3)
      real(rt) ugdnvy_out(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3)
      real(rt) ugdnvz_out(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3)
      real(rt) divu_cc(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3)
      real(rt) dx, dy, dz, dt
      real(rt) dtdx, dtdy, dtdz, hdt
      real(rt) cdtdx, cdtdy, cdtdz
      real(rt) hdtdx, hdtdy, hdtdz
      real(rt) a_old,a_new,a_half

      ! Left and right state arrays (edge centered, cell centered)
      real(rt), pointer:: dqx(:,:,:,:), dqy(:,:,:,:), dqz(:,:,:,:)
      real(rt), pointer::qxm(:,:,:,:),qym(:,:,:,:), qzm(:,:,:,:)
      real(rt), pointer::qxp(:,:,:,:),qyp(:,:,:,:), qzp(:,:,:,:)

      real(rt), pointer::qmxy(:,:,:,:),qpxy(:,:,:,:)
      real(rt), pointer::qmxz(:,:,:,:),qpxz(:,:,:,:)

      real(rt), pointer::qmyx(:,:,:,:),qpyx(:,:,:,:)
      real(rt), pointer::qmyz(:,:,:,:),qpyz(:,:,:,:)

      real(rt), pointer::qmzx(:,:,:,:),qpzx(:,:,:,:)
      real(rt), pointer::qmzy(:,:,:,:),qpzy(:,:,:,:)

      real(rt), pointer::qxl(:,:,:,:),qxr(:,:,:,:)
      real(rt), pointer::qyl(:,:,:,:),qyr(:,:,:,:)
      real(rt), pointer::qzl(:,:,:,:),qzr(:,:,:,:)

      ! Work arrays to hold 3 planes of riemann state and conservative fluxes
      real(rt), pointer::   fx(:,:,:,:),  fy(:,:,:,:), fz(:,:,:,:)

      real(rt), pointer::fxy(:,:,:,:),fxz(:,:,:,:)
      real(rt), pointer::fyx(:,:,:,:),fyz(:,:,:,:)
      real(rt), pointer::fzx(:,:,:,:),fzy(:,:,:,:)

      real(rt), pointer:: pgdnvx(:,:,:), ugdnvx(:,:,:)
      real(rt), pointer:: pgdnvxf(:,:,:), ugdnvxf(:,:,:)
      real(rt), pointer:: pgdnvtmpx(:,:,:), ugdnvtmpx(:,:,:)

      real(rt), pointer:: pgdnvy(:,:,:), ugdnvy(:,:,:)
      real(rt), pointer:: pgdnvyf(:,:,:), ugdnvyf(:,:,:)
      real(rt), pointer:: pgdnvtmpy(:,:,:), ugdnvtmpy(:,:,:)

      real(rt), pointer:: pgdnvz(:,:,:), ugdnvz(:,:,:)
      real(rt), pointer:: pgdnvtmpz1(:,:,:), ugdnvtmpz1(:,:,:)
      real(rt), pointer:: pgdnvtmpz2(:,:,:), ugdnvtmpz2(:,:,:)

      real(rt), pointer:: pgdnvzf(:,:,:), ugdnvzf(:,:,:)

      real(rt), pointer:: Ip(:,:,:,:,:,:), Im(:,:,:,:,:,:)
      real(rt), pointer:: Ip_g(:,:,:,:,:,:), Im_g(:,:,:,:,:,:)

      integer :: dnv_lo(3), dnv_hi(3)
      integer :: q_lo(3), q_hi(3)
      integer :: fx_lo(3), fx_hi(3)
      integer :: fy_lo(3), fy_hi(3)
      integer :: fz_lo(3), fz_hi(3)

      dnv_lo = [ilo1-1, ilo2-1, 1]
      dnv_hi = [ihi1+2, ihi2+2, 2]

      q_lo = [ilo1-1, ilo2-1, 1]
      q_hi = [ihi1+2, ihi2+2, 2]

      fx_lo = [ilo1  , ilo2-1, 1]
      fx_hi = [ihi1+1, ihi2+1, 2]

      fy_lo = [ilo1-1, ilo2  , 1]
      fy_hi = [ihi1+1, ihi2+1, 2]


      call amrex_allocate ( pgdnvx   ,dnv_lo, dnv_hi)
      call amrex_allocate ( ugdnvx   ,dnv_lo, dnv_hi)
      call amrex_allocate ( pgdnvxf  ,dnv_lo, dnv_hi)
      call amrex_allocate ( ugdnvxf  ,dnv_lo, dnv_hi)
      call amrex_allocate ( pgdnvtmpx,dnv_lo, dnv_hi)
      call amrex_allocate ( ugdnvtmpx,dnv_lo, dnv_hi)

      call amrex_allocate ( pgdnvy   ,dnv_lo, dnv_hi)
      call amrex_allocate ( ugdnvy   ,dnv_lo, dnv_hi)
      call amrex_allocate ( pgdnvyf  ,dnv_lo, dnv_hi)
      call amrex_allocate ( ugdnvyf  ,dnv_lo, dnv_hi)
      call amrex_allocate ( pgdnvtmpy,dnv_lo, dnv_hi)
      call amrex_allocate ( ugdnvtmpy,dnv_lo, dnv_hi)

      call amrex_allocate ( pgdnvz    ,dnv_lo, dnv_hi)
      call amrex_allocate ( ugdnvz    ,dnv_lo, dnv_hi)
      call amrex_allocate ( pgdnvtmpz1,dnv_lo, dnv_hi)
      call amrex_allocate ( ugdnvtmpz1,dnv_lo, dnv_hi)
      call amrex_allocate ( pgdnvtmpz2,dnv_lo, dnv_hi)
      call amrex_allocate ( ugdnvtmpz2,dnv_lo, dnv_hi)
      call amrex_allocate ( pgdnvzf   ,dnv_lo, dnv_hi)
      call amrex_allocate ( ugdnvzf   ,dnv_lo, dnv_hi)

      call amrex_allocate ( dqx, q_lo, q_hi, QVAR)
      call amrex_allocate ( dqy, q_lo, q_hi, QVAR)
      call amrex_allocate ( dqz, q_lo, q_hi, QVAR)

      ! One-sided states on x-edges
      call amrex_allocate ( qxm , q_lo, q_hi, QVAR)
      call amrex_allocate ( qxp , q_lo, q_hi, QVAR)
      call amrex_allocate ( qmxy, q_lo, q_hi, QVAR)
      call amrex_allocate ( qpxy, q_lo, q_hi, QVAR)
      call amrex_allocate ( qmxz, q_lo, q_hi, QVAR)
      call amrex_allocate ( qpxz, q_lo, q_hi, QVAR)
      call amrex_allocate ( qxl , q_lo, q_hi, QVAR)
      call amrex_allocate ( qxr , q_lo, q_hi, QVAR)

      ! One-sided states on y-edges
      call amrex_allocate ( qym , q_lo, q_hi, QVAR)
      call amrex_allocate ( qyp , q_lo, q_hi, QVAR)
      call amrex_allocate ( qmyx, q_lo, q_hi, QVAR)
      call amrex_allocate ( qpyx, q_lo, q_hi, QVAR)
      call amrex_allocate ( qmyz, q_lo, q_hi, QVAR)
      call amrex_allocate ( qpyz, q_lo, q_hi, QVAR)
      call amrex_allocate ( qyl , q_lo, q_hi, QVAR)
      call amrex_allocate ( qyr , q_lo, q_hi, QVAR)

      ! One-sided states on z-edges
      call amrex_allocate ( qzm , q_lo, q_hi, QVAR)
      call amrex_allocate ( qzp , q_lo, q_hi, QVAR)
      call amrex_allocate ( qmzx, q_lo, q_hi, QVAR)
      call amrex_allocate ( qpzx, q_lo, q_hi, QVAR)
      call amrex_allocate ( qmzy, q_lo, q_hi, QVAR)
      call amrex_allocate ( qpzy, q_lo, q_hi, QVAR)
      call amrex_allocate ( qzl , q_lo, q_hi, QVAR)
      call amrex_allocate ( qzr , q_lo, q_hi, QVAR)

      ! Output of cmpflx on x-edges
      call amrex_allocate ( fx , fx_lo, fx_hi, NVAR)
      call amrex_allocate ( fxy, fx_lo, fx_hi, NVAR)
      call amrex_allocate ( fxz, fx_lo, fx_hi, NVAR)

      ! Output of cmpflx on y-edges
      call amrex_allocate ( fy , fy_lo, fy_hi, NVAR)
      call amrex_allocate ( fyx, fy_lo, fy_hi, NVAR)
      call amrex_allocate ( fyz, fy_lo, fy_hi, NVAR)

      ! Output of cmpflx on z-edges
      fz_lo = [ilo1-1, ilo2-1, 1]
      fz_hi = [ihi1+1, ihi2+1, 2]
      call amrex_allocate ( fz , fz_lo, fz_hi, NVAR)
      fz_lo = [ilo1, ilo2-1, 1]
      fz_hi = [ihi1, ihi2+1, 2]
      call amrex_allocate ( fzx, fz_lo, fz_hi, NVAR)
      fz_lo = [ilo1-1, ilo2, 1]
      fz_hi = [ihi1+1, ihi2, 2]
      call amrex_allocate ( fzy, fz_lo, fz_hi, NVAR)

      ! x-index, y-index, z-index, dim, characteristics, variables
      call amrex_allocate ( Ip,ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,QVAR)
      call amrex_allocate ( Im,ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,QVAR)

      call amrex_allocate (Ip_g,ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,3)
      call amrex_allocate (Im_g,ilo1-1,ihi1+1,ilo2-1,ihi2+1,1,2,1,3,1,3,1,3)

      a_half = HALF * (a_old + a_new)

      ! Local constants
      dtdx = dt/dx
      dtdy = dt/dy
      dtdz = dt/dz
      hdt = HALF*dt

      ! Note that we divide by a_half here instead of in every loop.
      hdtdx = HALF*dtdx/a_half
      hdtdy = HALF*dtdy/a_half
      hdtdz = HALF*dtdz/a_half

      ! Note that we divide by a_half here instead of in every loop.
      cdtdx = THIRD*dtdx/a_half
      cdtdy = THIRD*dtdy/a_half
      cdtdz = THIRD*dtdz/a_half

      ! Initialize kc (current k-level) and km (previous k-level)
      kc = 1
      km = 2

      do k3d = ilo3-1, ihi3+1

         ! Swap pointers to levels
         kt = km
         km = kc
         kc = kt

         if (ppm_type .gt. 0) then

            do n=1,QVAR
               call ppm(q(:,:,:,n),qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        q(:,:,:,QU:),c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        Ip(:,:,:,:,:,n),Im(:,:,:,:,:,n), &
                        ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc,a_old)
            end do

            if (version_2 .eq. 2) then
               do n=1,3
                  call ppm(srcQ(:,:,:,QU+n-1),srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                           q(:,:,:,QU:),c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           Ip_g(:,:,:,:,:,n),Im_g(:,:,:,:,:,n), &
                           ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc,a_old)
               end do
            else
               Ip_g(:,:,:,:,:,:) = ZERO
               Im_g(:,:,:,:,:,:) = ZERO
            end if

            ! Compute U_x and U_y at kc (k3d)
            call tracexy_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                             Ip,Im,Ip_g,Im_g, &
                             qxm,qxp,qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                             srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                             ilo1,ilo2,ihi1,ihi2,dt,a_old,kc,k3d)

         else if (ppm_type .eq. 0) then

            ! Compute all slopes at kc (k3d)
            call uslope(q,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        dqx,dqy,dqz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        ilo1,ilo2,ihi1,ihi2,kc,k3d,QVAR)

            ! Compute U_x and U_y at kc (k3d)
            if (use_colglaz .eq. 1) then
               call tracexy_cg(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                               dqx,dqy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                               qxm,qxp,qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                               ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d,a_old)
            else
               call tracexy(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                            dqx,dqy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                            qxm,qxp,qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                            ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d,a_old)
            end if

         else 
            print *,'>>> ... we only support ppm_type >= 0, not: ',ppm_type 
            call amrex_error("Error:: Nyx_advection_3d.f90 :: umeth3d")
         end if

         ! On x-edges -- choose state fx based on qxm, qxp
         call cmpflx(qxm,qxp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                     fx,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                     ugdnvx,pgdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                     csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                     1,ilo1,ihi1+1,ilo2-1,ihi2+1,kc,kc,k3d,print_fortran_warnings)

         ! On y-edges -- choose state fy based on qym, qyp
         call cmpflx(qym,qyp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                     fy,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                     ugdnvy,pgdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                     csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                     2,ilo1-1,ihi1+1,ilo2,ihi2+1,kc,kc,k3d,print_fortran_warnings)

 
         if (corner_coupling .eq. 1) then

             ! On x-edges
             ! qxm + d/dy (fy) --> qmxy
             ! qxp + d/dy (fy) --> qpxy
             call transy1(qxm,qmxy,qxp,qpxy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                          fy,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                          ugdnvy,pgdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                          cdtdy,ilo1-1,ihi1+1,ilo2,ihi2,kc,k3d)

             ! On y-edges
             ! qym + d/dx (fx) --> qmyx
             ! qyp + d/dx (fx) --> qpyx
             call transx1(qym,qmyx,qyp,qpyx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                          fx,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                          ugdnvx,pgdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                          cdtdx,ilo1,ihi1,ilo2-1,ihi2+1,kc,k3d)

             ! On x-edges -- choose state fxy based on qmxy, qpxy
             call cmpflx(qmxy,qpxy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         fxy,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                         ugdnvtmpx,pgdnvtmpx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         1,ilo1,ihi1+1,ilo2,ihi2,kc,kc,k3d,print_fortran_warnings)

             ! On y-edges -- choose state fyx based on qmyx, qpyx
             call cmpflx(qmyx,qpyx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         fyx,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                         ugdnvtmpy,pgdnvtmpy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                         csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         2,ilo1,ihi1,ilo2,ihi2+1,kc,kc,k3d,print_fortran_warnings)

         end if

         if (k3d.ge.ilo3) then

            ! Compute U_z at kc (k3d)
            if (ppm_type .gt. 0) then
               call tracez_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                               Ip,Im,Ip_g,Im_g, &
                               qzm,qzp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                               srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                               ilo1,ilo2,ihi1,ihi2,dt,a_old,km,kc,k3d)
            else if (ppm_type .eq. 0) then
               if (use_colglaz .eq. 1) then
                  call tracez_cg(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                              dqz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                              qzm,qzp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                              ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d,a_old)
               else
                  call tracez(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                              dqz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                              qzm,qzp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                              ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d,a_old)
               end if
            else 
               print *,'>>> ... we only support ppm_type >= 0, not: ',ppm_type 
               call amrex_error("Error:: Nyx_advection_3d.f90 :: umeth3d")
            end if

            ! Compute \tilde{F}^z at kc (k3d)
            call cmpflx(qzm,qzp,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        fz,ilo1-1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                        ugdnvz,pgdnvz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        3,ilo1-1,ihi1+1,ilo2-1,ihi2+1,kc,kc,k3d,print_fortran_warnings)

            if (corner_coupling .eq. 1) then

                ! On z-edges
                ! qzm + d/dy (fy) --> qmzy
                ! qzp + d/dy (fy) --> qpzy
                call transy2(qzm,qmzy,qzp,qpzy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                             fy,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                             ugdnvy,pgdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                             cdtdy,ilo1-1,ihi1+1,ilo2,ihi2,kc,km,k3d)

                ! On z-edges
                ! qzm + d/dx (fx) --> qmzx
                ! qzp + d/dx (fx) --> qpzx
                call transx2(qzm,qmzx,qzp,qpzx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                             fx,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                             ugdnvx,pgdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                             cdtdx,ilo1,ihi1,ilo2-1,ihi2+1,kc,km,k3d)

                ! On z-edges -- choose state fzx based on qmzx, qpzx
                call cmpflx(qmzx,qpzx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                            fzx,ilo1,ilo2-1,1,ihi1,ihi2+1,2, &
                            ugdnvtmpz1,pgdnvtmpz1,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                            csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                            3,ilo1,ihi1,ilo2-1,ihi2+1,kc,kc,k3d,print_fortran_warnings)

                ! On z-edges -- choose state fzy based on qmzy, qpzy
                call cmpflx(qmzy,qpzy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                            fzy,ilo1-1,ilo2,1,ihi1+1,ihi2,2, &
                            ugdnvtmpz2,pgdnvtmpz2,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                            csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                            3,ilo1-1,ihi1+1,ilo2,ihi2,kc,kc,k3d,print_fortran_warnings)

                ! On z-edges
                ! qzm + d/dx (fxy) + d/dy (fyx) --> qzl
                ! qzp + d/dx (fxy) + d/dy (fyx) --> qzr
                call transxy(qzm,qzl,qzp,qzr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                             fxy,ilo1  ,ilo2-1,1,ihi1+1,ihi2+1,2, &
                             fyx,ilo1-1,ilo2  ,1,ihi1+1,ihi2+1,2, &
                             ugdnvtmpx,pgdnvtmpx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                             ugdnvtmpy,pgdnvtmpy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                             srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                             hdt,hdtdx,hdtdy,ilo1,ihi1,ilo2,ihi2,kc,km,k3d,a_old,a_new)

            else
                ! On z-edges
                ! qzm + d/dx (fx) + d/dy (fy) --> qzl
                ! qzp + d/dx (fx) + d/dy (fy) --> qzr
                call transxy(qzm,qzl,qzp,qzr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                             fx,ilo1  ,ilo2-1,1,ihi1+1,ihi2+1,2, &
                             fy,ilo1-1,ilo2  ,1,ihi1+1,ihi2+1,2, &
                             ugdnvx,pgdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                             ugdnvy,pgdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                             srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                             hdt,hdtdx,hdtdy,ilo1,ihi1,ilo2,ihi2,kc,km,k3d,a_old,a_new)
            endif

            if (version_2 .eq. 3) then
               call tracez_src(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                               qzl,qzr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                               srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                               ilo1,ilo2,ihi1,ihi2,dt,a_old,km,kc,k3d)
            endif

            ! Compute F^z at kc (k3d) -- note that flux3 is indexed by k3d, not kc
            ! On z-edges -- choose state fluxf3 based on qzl, qzr
            call cmpflx(qzl,qzr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        flux3,fd3_l1,fd3_l2,fd3_l3,fd3_h1,fd3_h2,fd3_h3, &
                        ugdnvzf,pgdnvzf,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                        csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        3,ilo1,ihi1,ilo2,ihi2,kc,k3d,k3d,print_fortran_warnings)

            do j=ilo2-1,ihi2+1
               do i=ilo1-1,ihi1+1
                  ugdnvz_out(i,j,k3d) = ugdnvzf(i,j,kc)
               end do
            end do

            if (k3d .ge. ilo3+1) then
               do j = ilo2,ihi2
                  do i = ilo1,ihi1
                     divu_cc(i,j,k3d-1) = divu_cc(i,j,k3d-1) +  &
                                (ugdnvzf(i,j,kc)-ugdnvzf(i,j,km))/dz
                  end do
               end do
            end if

            if (k3d.gt.ilo3) then

               if (corner_coupling .eq. 1) then

                   ! On x-edges:
                   ! qxm + d/dz (fz) --> qmxz
                   ! qxp + d/dz (fz) --> qpxz
                   ! On y-edges:
                   ! qym + d/dz (fz) --> qmyz 
                   ! qyp + d/dz (fz) --> qpyz
                   call transz(qxm,qmxz,qxp,qpxz, &
                               qym,qmyz,qyp,qpyz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                               fz,ilo1-1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                               ugdnvz,pgdnvz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                               cdtdz,ilo1-1,ihi1+1,ilo2-1,ihi2+1,km,kc,k3d)
    
                   ! On x-edges -- choose state fxz based on qmxz, qpxz
                   call cmpflx(qmxz,qpxz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                               fxz,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                               ugdnvx,pgdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                               csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                               1,ilo1,ihi1+1,ilo2-1,ihi2+1,km,km,k3d-1,print_fortran_warnings)

                   ! On y-edges -- choose state fyz based on qmyz, qpyz
                   call cmpflx(qmyz,qpyz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                               fyz,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                               ugdnvy,pgdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                               csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                               2,ilo1-1,ihi1+1,ilo2,ihi2+1,km,km,k3d-1,print_fortran_warnings)
    
                   ! On x-edges:
                   ! qxm + d/dy (fyz) + d/dz(fzy) --> qxl
                   ! qxp + d/dy (fyz) + d/dz(fzy) --> qxr
                   call transyz(qxm,qxl,qxp,qxr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                                fyz,ilo1-1,ilo2,1,ihi1+1,ihi2+1,2, &
                                fzy,ilo1-1,ilo2,1,ihi1+1,ihi2  ,2, &
                                ugdnvy,pgdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                                ugdnvtmpz2,pgdnvtmpz2,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                                srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                                hdt,hdtdy,hdtdz,ilo1-1,ihi1+1,ilo2,ihi2,km,kc,k3d-1,a_old,a_new)

                   ! On y-edges:
                   ! qym + d/dx (fxz) + d/dz(fzx) --> qyl
                   ! qyp + d/dx (fxz) + d/dz(fzx) --> qyr
                   call transxz(qym,qyl,qyp,qyr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                                fxz,ilo1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                                fzx,ilo1,ilo2-1,1,ihi1  ,ihi2+1,2, &
                                ugdnvx,pgdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                                ugdnvtmpz1,pgdnvtmpz1,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                                srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                                hdt,hdtdx,hdtdz,ilo1,ihi1,ilo2-1,ihi2+1,km,kc,k3d-1,a_old,a_new)

               else

                   ! On x-edges:
                   ! qxm + d/dy (fy) + d/dz(fz) --> qxl
                   ! qxp + d/dy (fy) + d/dz(fz) --> qxr
                   call transyz(qxm,qxl,qxp,qxr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                                fy,ilo1-1,ilo2  ,1,ihi1+1,ihi2+1,2, &
                                fz,ilo1-1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                                ugdnvy,pgdnvy,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                                ugdnvz,pgdnvz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                                srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                                hdt,hdtdy,hdtdz,ilo1-1,ihi1+1,ilo2,ihi2,km,kc,k3d-1,a_old,a_new)

                   ! On y-edges:
                   ! qym + d/dx (fx) + d/dz(fz) --> qyl
                   ! qyp + d/dx (fx) + d/dz(fz) --> qyr
                   call transxz(qym,qyl,qyp,qyr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                                fx,ilo1  ,ilo2-1,1,ihi1+1,ihi2+1,2, &
                                fz,ilo1-1,ilo2-1,1,ihi1+1,ihi2+1,2, &
                                ugdnvx,pgdnvx,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                                ugdnvz,pgdnvz,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                                srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                                hdt,hdtdx,hdtdz,ilo1,ihi1,ilo2-1,ihi2+1,km,kc,k3d-1,a_old,a_new)

               end if

               if (version_2 .eq. 3) then
                  call tracex_src(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                                  qxl,qxr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                                  srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                                  ilo1,ilo2,ihi1,ihi2,dt,a_old,kc,k3d)
                  call tracey_src(q,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                                  qyl,qyr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                                  srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                                  ilo1,ilo2,ihi1,ihi2,dt,a_old,kc,k3d)
               endif

               ! On x-edges -- choose state flux1 based on qxl, qxr
               call cmpflx(qxl,qxr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           flux1,fd1_l1,fd1_l2,fd1_l3,fd1_h1,fd1_h2,fd1_h3, &
                           ugdnvxf,pgdnvxf,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           1,ilo1,ihi1+1,ilo2,ihi2,km,k3d-1,k3d-1,print_fortran_warnings)

               do j=ilo2-1,ihi2+1
                  do i=ilo1-1,ihi1+2
                     ugdnvx_out(i,j,k3d-1) = ugdnvxf(i,j,km)
                  end do
               end do

               ! On y-edges -- choose state flux2 based on qyl, qyr
               call cmpflx(qyl,qyr,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           flux2,fd2_l1,fd2_l2,fd2_l3,fd2_h1,fd2_h2,fd2_h3, &
                           ugdnvyf,pgdnvyf,ilo1-1,ilo2-1,1,ihi1+2,ihi2+2,2, &
                           csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                           2,ilo1,ihi1,ilo2,ihi2+1,km,k3d-1,k3d-1,print_fortran_warnings)

               do j=ilo2-1,ihi2+2
                  do i=ilo1-1,ihi1+1
                     ugdnvy_out(i,j,k3d-1) = ugdnvyf(i,j,km)
                  end do
               end do

               do j = ilo2,ihi2
                  do i = ilo1,ihi1
                     divu_cc(i,j,k3d-1) = divu_cc(i,j,k3d-1) +  &
                          (ugdnvxf(i+1,j,km)-ugdnvxf(i,j,km))/dx + &
                          (ugdnvyf(i,j+1,km)-ugdnvyf(i,j,km))/dy
                  end do
               end do

            end if
         end if
      enddo

      call amrex_deallocate(pgdnvx)
      call amrex_deallocate(pgdnvxf)
      call amrex_deallocate(pgdnvtmpx)
      call amrex_deallocate(ugdnvx)
      call amrex_deallocate(ugdnvxf)
      call amrex_deallocate(ugdnvtmpx)

      call amrex_deallocate(pgdnvy)
      call amrex_deallocate(pgdnvyf)
      call amrex_deallocate(pgdnvtmpy)
      call amrex_deallocate(ugdnvy)
      call amrex_deallocate(ugdnvyf)
      call amrex_deallocate(ugdnvtmpy)

      call amrex_deallocate ( pgdnvz    )
      call amrex_deallocate ( ugdnvz    )
      call amrex_deallocate ( pgdnvtmpz1)
      call amrex_deallocate ( ugdnvtmpz1)
      call amrex_deallocate ( pgdnvtmpz2)
      call amrex_deallocate ( ugdnvtmpz2)
      call amrex_deallocate ( pgdnvzf   )
      call amrex_deallocate ( ugdnvzf   )

      call amrex_deallocate ( dqx)
      call amrex_deallocate ( dqy)
      call amrex_deallocate ( dqz)

      call amrex_deallocate ( qxm )
      call amrex_deallocate ( qxp )
      call amrex_deallocate ( qmxy)
      call amrex_deallocate ( qpxy)
      call amrex_deallocate ( qmxz)
      call amrex_deallocate ( qpxz)
      call amrex_deallocate ( qxl )
      call amrex_deallocate ( qxr )

      call amrex_deallocate ( qym )
      call amrex_deallocate ( qyp )
      call amrex_deallocate ( qmyx)
      call amrex_deallocate ( qpyx)
      call amrex_deallocate ( qmyz)
      call amrex_deallocate ( qpyz)
      call amrex_deallocate ( qyl )
      call amrex_deallocate ( qyr )

      call amrex_deallocate ( qzm )
      call amrex_deallocate ( qzp )
      call amrex_deallocate ( qmzx)
      call amrex_deallocate ( qpzx)
      call amrex_deallocate ( qmzy)
      call amrex_deallocate ( qpzy)
      call amrex_deallocate ( qzl )
      call amrex_deallocate ( qzr )

      call amrex_deallocate ( fx )
      call amrex_deallocate ( fxy)
      call amrex_deallocate ( fxz)

      call amrex_deallocate ( fy )
      call amrex_deallocate ( fyx)
      call amrex_deallocate ( fyz)

      call amrex_deallocate ( fz )
      call amrex_deallocate ( fzx)
      call amrex_deallocate ( fzy)

      call amrex_deallocate ( Ip)
      call amrex_deallocate ( Im)

      call amrex_deallocate (Ip_g)
      call amrex_deallocate (Im_g)


      end subroutine umeth3d

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine ctoprim(lo,hi,          uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                         q,c,csml,flatn,  q_l1,  q_l2,  q_l3,  q_h1,  q_h2,  q_h3, &
                         src,  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                         srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                         grav,gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                         dx,dy,dz,dt,ngq,ngf,a_old,a_new)
      !
      !     Will give primitive variables on lo-ngq:hi+ngq, and flatn on lo-ngf:hi+ngf
      !     if use_flattening=1.  Declared dimensions of q,c,csml,flatn are given
      !     by DIMS(q).  This declared region is assumed to encompass lo-ngq:hi+ngq.
      !     Also, uflaten call assumes ngq>=ngf+3 (ie, primitve data is used by the
      !     routine that computes flatn).
      !
      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use network, only : nspec, naux
      use eos_params_module
      use eos_module
      use flatten_module
      use amrex_constants_module
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                     UEDEN, UEINT, UFA, UFS, &
                                     QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QPRES, QFA, QFS, &
                                     nadv, small_dens, small_pres, &
                                     gamma_const, gamma_minus_1, use_flattening

      implicit none

      real(rt), parameter:: small = 1.d-8

      integer lo(3), hi(3)
      integer  uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3
      integer    q_l1,   q_l2,   q_l3,   q_h1,   q_h2,   q_h3
      integer   gv_l1,  gv_l2,  gv_l3,  gv_h1,  gv_h2,  gv_h3
      integer  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
      integer srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3

      real(rt) :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
      real(rt) ::     q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
      real(rt) ::     c(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) ::  csml(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) :: flatn(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) ::   src( src_l1: src_h1, src_l2: src_h2, src_l3: src_h3,NVAR)
      real(rt) ::  srcQ(srcq_l1:srcq_h1,srcq_l2:srcq_h2,srcq_l3:srcq_h3,QVAR)
      real(rt) :: grav( gv_l1: gv_h1, gv_l2: gv_h2, gv_l3: gv_h3,3)
      real(rt) :: dx, dy, dz, dt, a_old, a_new
      real(rt) :: dpdr, dpde

      integer          :: i, j, k
      integer          :: ngq, ngf, loq(3), hiq(3)
      integer          :: n, nq
      integer          :: iadv, ispec
      real(rt) :: courx, coury, courz
      real(rt) :: a_half, a_dot, rhoInv
      real(rt) :: dtdxaold, dtdyaold, dtdzaold, small_pres_over_dens

      do i=1,3
         loq(i) = lo(i)-ngq
         hiq(i) = hi(i)+ngq
      enddo
      !
      ! Make q (all but p), except put e in slot for rho.e, fix after eos call.
      ! The temperature is used as an initial guess for the eos call and will be overwritten.
      !
      do k = loq(3),hiq(3)
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)

               if (uin(i,j,k,URHO) .le. ZERO) then
                  !
                  ! A critical region since we usually can't write from threads.
                  !
                  print *,'   '
                  print *,'>>> Error: Nyx_advection_3d::ctoprim ',i,j,k
                  print *,'>>> ... negative density ',uin(i,j,k,URHO)
                  call amrex_error("Error:: Nyx_advection_3d.f90 :: ctoprim")
               end if

               rhoInv = ONE/uin(i,j,k,URHO)

               q(i,j,k,QRHO) = uin(i,j,k,URHO)
               q(i,j,k,QU)   = uin(i,j,k,UMX)*rhoInv
               q(i,j,k,QV)   = uin(i,j,k,UMY)*rhoInv
               q(i,j,k,QW)   = uin(i,j,k,UMZ)*rhoInv

               ! Convert "rho e" to "e"
               q(i,j,k,QREINT ) = uin(i,j,k,UEINT)*rhoInv

            enddo
         enddo
      enddo

      ! Load advected quatities, c, into q, assuming they arrived in uin as rho.c
      do iadv = 1, nadv
         n = UFA + iadv - 1
         nq = QFA + iadv - 1
         do k = loq(3),hiq(3)
            do j = loq(2),hiq(2)
               do i = loq(1),hiq(1)
                  q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
               enddo
            enddo
         enddo
      enddo

      ! Load chemical species and aux. variables, c, into q, assuming they arrived in uin as rho.c
      if (UFS .gt. 0) then
         do ispec = 1, nspec+naux
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
      end if ! UFS > 0

      small_pres_over_dens = small_pres / small_dens

      ! Get p, T, c, csml using q state
      do k = loq(3), hiq(3)
         do j = loq(2), hiq(2)
            do i = loq(1), hiq(1)

               ! If necessary, reset the energy using small_temp
               if (q(i,j,k,QREINT) .lt. ZERO) then

!                 HACK HACK HACK 
!                 call nyx_eos_given_RT(q(i,j,k,QREINT),q(i,j,k,QPRES),q(i,j,k,QRHO), &
!                                       small_temp,diag_eos(i,j,k,NE_COMP),a_old)

                  if (q(i,j,k,QREINT) .lt. ZERO) then
                     !
                     ! A critical region since we usually can't write from threads.
                     !
                     print *,'   '
                     print *,'>>> Error: Nyx_advection_3d::ctoprim ',i,j,k
                     print *,'>>> ... new e from eos_given_RT call is negative ',q(i,j,k,QREINT)
                     print *,'    '
                     call amrex_error("Error:: Nyx_advection_3d.f90 :: ctoprim")
                  end if
               end if

               ! Define the soundspeed from the EOS
               call nyx_eos_soundspeed(c(i,j,k), q(i,j,k,QRHO), q(i,j,k,QREINT))

               ! Set csmal based on small_pres and small_dens
               csml(i,j,k) = sqrt(gamma_const * small_pres_over_dens)

               ! Convert "e" back to "rho e"
               q(i,j,k,QREINT) = q(i,j,k,QREINT)*q(i,j,k,QRHO)

               ! Pressure = (gamma - 1) * rho * e
               q(i,j,k,QPRES) = gamma_minus_1 * q(i,j,k,QREINT)

            end do
         end do
      end do

      a_half = HALF * (a_old + a_new)
      a_dot   = (a_new - a_old) / dt

      ! Make sure these are initialized to zero.
      srcQ = ZERO

      ! NOTE - WE ASSUME HERE THAT src(i,j,k,URHO) = 0. --
      !        IF NOT THEN THE FORMULAE BELOW ARE INCOMPLETE.

      ! compute srcQ terms
      do k = loq(3), hiq(3)
         do j = loq(2), hiq(2)
            do i = loq(1), hiq(1)

               rhoInv = ONE/q(i,j,k,QRHO)

               srcQ(i,j,k,QRHO  ) = src(i,j,k,URHO)
               srcQ(i,j,k,QU    ) = src(i,j,k,UMX) * rhoInv - a_dot * q(i,j,k,QU) + &
                                    grav(i,j,k,1)
               srcQ(i,j,k,QV    ) = src(i,j,k,UMY) * rhoInv - a_dot * q(i,j,k,QV) + &
                                    grav(i,j,k,2)
               srcQ(i,j,k,QW    ) = src(i,j,k,UMZ) * rhoInv - a_dot * q(i,j,k,QW) + &
                                    grav(i,j,k,3)
               srcQ(i,j,k,QREINT) = src(i,j,k,UEDEN) - q(i,j,k,QU)*src(i,j,k,UMX) - &
                                                       q(i,j,k,QV)*src(i,j,k,UMY) - &
                                                       q(i,j,k,QW)*src(i,j,k,UMZ) - &
                                                       a_dot * THREE * gamma_minus_1 * q(i,j,k,QREINT)

               dpde = gamma_minus_1 * q(i,j,k,QRHO)
               dpdr = gamma_minus_1 * q(i,j,k,QREINT)/q(i,j,k,QRHO)
               srcQ(i,j,k,QPRES ) = dpde * srcQ(i,j,k,QREINT) * rhoInv &
                                  + dpdr * srcQ(i,j,k,QRHO)

               if (UFS .gt. 0) then
                  do ispec = 1,nspec+naux
                     srcQ(i,j,k,QFS+ispec-1) = src(i,j,k,UFS+ispec-1)*rhoInv
                  enddo
               end if ! UFS > 0

               do iadv = 1,nadv
                  srcQ(i,j,k,QFA+iadv-1) = src(i,j,k,UFA+iadv-1)*rhoInv
               enddo

            enddo
         enddo
      enddo

      dtdxaold = dt / dx / a_old
      dtdyaold = dt / dy / a_old
      dtdzaold = dt / dz / a_old

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               courx = ( c(i,j,k)+abs(q(i,j,k,QU)) ) * dtdxaold
               coury = ( c(i,j,k)+abs(q(i,j,k,QV)) ) * dtdyaold
               courz = ( c(i,j,k)+abs(q(i,j,k,QW)) ) * dtdzaold

               if (courx .gt. ONE) then
                  !
                  ! A critical region since we usually can't write from threads.
                  !
                  print *,'   '
                  print *,'>>> ... (u+c) * a * dt / dx > 1 ', courx
                  print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                  print *,'>>> ... u, c                ',q(i,j,k,QU), c(i,j,k)
                  print *,'>>> ... density             ',q(i,j,k,QRHO)
                  call amrex_error("Error:: Nyx_advection_3d.f90 :: CFL violation in x-dir in ctoprim")
               end if

               if (coury .gt. ONE) then
                  !
                  ! A critical region since we usually can't write from threads.
                  !
                  print *,'   '
                  print *,'>>> ... (v+c) * a * dt / dx > 1 ', coury
                  print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                  print *,'>>> ... v, c                ',q(i,j,k,QV), c(i,j,k)
                  print *,'>>> ... density             ',q(i,j,k,QRHO)
                  call amrex_error("Error:: Nyx_advection_3d.f90 :: CFL violation in y-dir in ctoprim")
               end if

               if (courz .gt. ONE) then
                  !
                  ! A critical region since we usually can't write from threads.
                  !
                  print *,'   '
                  print *,'>>> ... (w+c) * a * dt / dx > 1 ', courz
                  print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                  print *,'>>> ... w, c                ',q(i,j,k,QW), c(i,j,k)
                  print *,'>>> ... density             ',q(i,j,k,QRHO)
                  call amrex_error("Error:: Nyx_advection_3d.f90 :: CFL violation in z-dir in ctoprim")
               end if

            enddo
         enddo
      enddo

      ! Compute flattening coef for slope calculations
      if (use_flattening == 1) then
         do n=1,3
            loq(n)=lo(n)-ngf
            hiq(n)=hi(n)+ngf
         enddo
         call uflaten(loq,hiq, &
                      q(q_l1,q_l2,q_l3,QPRES), &
                      q(q_l1,q_l2,q_l3,QU), &
                      q(q_l1,q_l2,q_l3,QV), &
                      q(q_l1,q_l2,q_l3,QW), &
                      flatn,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3)
      else
         flatn = ONE
      endif

      end subroutine ctoprim
! :::
! ::: ------------------------------------------------------------------
! :::

    subroutine consup(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                      hydro_src ,hsrc_l1,hsrc_l2,hsrc_l3,hsrc_h1,hsrc_h2,hsrc_h3, &
                      flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                      flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                      flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                      divu_nd,divu_cc,d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                      lo,hi,dx,dy,dz,dt,a_old,a_new)

      use amrex_fort_module, only : rt => amrex_real
      use amrex_constants_module
      use meth_params_module, only : difmag, NVAR, URHO, UMX, UMZ, &
           UEDEN, UEINT, UFS, normalize_species, gamma_minus_1

      implicit none

      integer lo(3), hi(3)
      integer   uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3
      integer   hsrc_l1,  hsrc_l2,  hsrc_l3,  hsrc_h1,  hsrc_h2,  hsrc_h3
      integer d_l1,   d_l2,   d_l3,   d_h1,   d_h2,   d_h3
      integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
      integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
      integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3

      real(rt)  :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NVAR)
      real(rt)  :: hydro_src(hsrc_l1:hsrc_h1,hsrc_l2:hsrc_h2,hsrc_l3:hsrc_h3,NVAR)
      real(rt)  :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
      real(rt)  :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
      real(rt)  :: flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
      real(rt)  :: divu_nd(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
      real(rt)  :: divu_cc(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3)
      real(rt)  :: dx, dy, dz, dt, a_old, a_new

      real(rt) :: div1, a_half, a_oldsq, a_newsq
      real(rt) :: area1, area2, area3
      real(rt) :: vol, volinv, a_newsq_inv
      real(rt) :: a_half_inv, a_new_inv, dt_a_new
      integer  :: i, j, k, n

      a_half  = HALF * (a_old + a_new)
      a_oldsq = a_old * a_old
      a_newsq = a_new * a_new
      vol     = dx * dy * dz
      volinv  = ONE / vol

      area1   =      dy * dz
      area2   = dx      * dz
      area3   = dx * dy

      do n = 1, NVAR

            do k = lo(3),hi(3)
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)+1
                     div1 = FOURTH*(divu_nd(i,j,k) + divu_nd(i,j+1,k) + divu_nd(i,j,k+1) + divu_nd(i,j+1,k+1))
                     div1 = difmag*min(ZERO,div1)
                     flux1(i,j,k,n) = flux1(i,j,k,n) + dx*div1*(uin(i,j,k,n)-uin(i-1,j,k,n))
                     flux1(i,j,k,n) = flux1(i,j,k,n) * area1 * dt
                  enddo
               enddo
            enddo
            do k = lo(3),hi(3)
               do j = lo(2),hi(2)+1
                  do i = lo(1),hi(1)
                     div1 = FOURTH*(divu_nd(i,j,k) + divu_nd(i+1,j,k) + divu_nd(i,j,k+1) + divu_nd(i+1,j,k+1))
                     div1 = difmag*min(ZERO,div1)
                     flux2(i,j,k,n) = flux2(i,j,k,n) + dy*div1*(uin(i,j,k,n)-uin(i,j-1,k,n))
                     flux2(i,j,k,n) = flux2(i,j,k,n) * area2 * dt
                  enddo
               enddo
            enddo
            do k = lo(3),hi(3)+1
               do j = lo(2),hi(2)
                  do i = lo(1),hi(1)
                     div1 = FOURTH*(divu_nd(i,j,k) + divu_nd(i+1,j,k) + divu_nd(i,j+1,k) + divu_nd(i+1,j+1,k))
                     div1 = difmag*min(ZERO,div1)
                     flux3(i,j,k,n) = flux3(i,j,k,n) + dz*div1*(uin(i,j,k,n)-uin(i,j,k-1,n))
                     flux3(i,j,k,n) = flux3(i,j,k,n) * area3 * dt
                  enddo
               enddo
            enddo

      enddo

      if (UFS .gt. 0 .and. normalize_species .eq. 1) &
         call normalize_species_fluxes( &
                  flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                  flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                  flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                  lo,hi)

      a_half_inv  = ONE / a_half
      a_new_inv   = ONE / a_new
      dt_a_new    = dt / a_new
      a_newsq_inv = ONE / a_newsq

      do n = 1, NVAR
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)

                  ! Density
                  if (n .eq. URHO) then
                      hydro_src(i,j,k,n) = &
                          ( ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                          +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                          +   flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * volinv  ) * a_half_inv

                  ! Momentum
                  else if (n .ge. UMX .and. n .le. UMZ) then
                     hydro_src(i,j,k,n) =  &
                            ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                          +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                          +   flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * volinv

                  ! (rho E)
                  else if (n .eq. UEDEN) then
                     hydro_src(i,j,k,n) =  &
                            ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                          +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                          +   flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * a_half * volinv &
                          +   a_half * (a_new - a_old) * ( TWO - THREE * gamma_minus_1) * uin(i,j,k,UEINT)

                  ! (rho e)
                  else if (n .eq. UEINT) then

                     hydro_src(i,j,k,n) =  &
                          ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                           +flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                           +flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * a_half * volinv &
                           +a_half*(a_new - a_old) * ( TWO - THREE * gamma_minus_1) * uin(i,j,k,UEINT) &
                           -a_half*dt*(gamma_minus_1 * uin(i,j,k,n)) * divu_cc(i,j,k)

                  ! (rho X_i) and (rho adv_i) and (rho aux_i)
                  else
                     hydro_src(i,j,k,n) = &
                            ( flux1(i,j,k,n) - flux1(i+1,j,k,n) &
                          +   flux2(i,j,k,n) - flux2(i,j+1,k,n) &
                          +   flux3(i,j,k,n) - flux3(i,j,k+1,n) ) * volinv 
                  endif

                  hydro_src(i,j,k,n) = hydro_src(i,j,k,n) / dt

               enddo
            enddo
         enddo
      enddo

      do n = 1, NVAR
         if (n .eq. URHO) then
            flux1(:,:,:,n) = flux1(:,:,:,n) * a_half_inv
            flux2(:,:,:,n) = flux2(:,:,:,n) * a_half_inv
            flux3(:,:,:,n) = flux3(:,:,:,n) * a_half_inv
         else if (n.ge.UMX .and. n.le.UMZ) then
            flux1(:,:,:,n) = flux1(:,:,:,n) * a_new_inv
            flux2(:,:,:,n) = flux2(:,:,:,n) * a_new_inv
            flux3(:,:,:,n) = flux3(:,:,:,n) * a_new_inv
         else if (n.eq.UEINT .or. n.eq.UEDEN) then
            flux1(:,:,:,n) = flux1(:,:,:,n) * (a_half * a_newsq_inv)
            flux2(:,:,:,n) = flux2(:,:,:,n) * (a_half * a_newsq_inv)
            flux3(:,:,:,n) = flux3(:,:,:,n) * (a_half * a_newsq_inv)
         else
            flux1(:,:,:,n) = flux1(:,:,:,n) * a_half_inv
            flux2(:,:,:,n) = flux2(:,:,:,n) * a_half_inv
            flux3(:,:,:,n) = flux3(:,:,:,n) * a_half_inv
         end if
      end do

      end subroutine consup

! :::
! ::: ------------------------------------------------------------------
! :::

     subroutine update_state(lo,hi, &
                             uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                             uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                             src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
                             hydro_src ,hsrc_l1,hsrc_l2,hsrc_l3,hsrc_h1,hsrc_h2,hsrc_h3, &
                             divu_cc,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                             dt,a_old,a_new,print_fortran_warnings) &
                             bind(C, name="fort_update_state")
 
      use amrex_fort_module, only : rt => amrex_real
      use amrex_constants_module
      use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, &
                                     gamma_minus_1, normalize_species
      use enforce_module, only : enforce_nonnegative_species

      implicit none

      integer, intent(in) :: lo(3), hi(3)
      integer, intent(in) ::print_fortran_warnings
      integer, intent(in) ::  uin_l1,   uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3
      integer, intent(in) ::  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
      integer, intent(in) ::   src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3
      integer, intent(in) ::  hsrc_l1, hsrc_l2, hsrc_l3, hsrc_h1, hsrc_h2, hsrc_h3
      integer, intent(in) ::  d_l1,d_l2,d_l3,d_h1,d_h2,d_h3

      real(rt), intent(in)  ::       uin( uin_l1: uin_h1, uin_l2: uin_h2, uin_l3: uin_h3,NVAR)
      real(rt), intent(out) ::      uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
      real(rt), intent(in)  ::       src( src_l1: src_h1, src_l2: src_h2, src_l3: src_h3,NVAR)
      real(rt), intent(in)  :: hydro_src(hsrc_l1:hsrc_h1,hsrc_l2:hsrc_h2,hsrc_l3:hsrc_h3,NVAR)
      real(rt), intent(in)  ::   divu_cc(   d_l1:   d_h1,   d_l2:   d_h2,   d_l3:   d_h3)
      real(rt), intent(in)  ::  dt, a_old, a_new

      real(rt) :: a_half, a_oldsq, a_newsq
      real(rt) :: a_new_inv, a_newsq_inv, a_half_inv, dt_a_new
      integer  :: i, j, k, n

      a_half     = HALF * (a_old + a_new)
      a_oldsq = a_old * a_old
      a_newsq = a_new * a_new

      a_half_inv  = ONE / a_half
      a_new_inv   = ONE / a_new
      a_newsq_inv = ONE / a_newsq

      do n = 1, NVAR

         ! Actually do the update here
         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)

                  ! Density
                  if (n .eq. URHO) then
                     uout(i,j,k,n) = uin(i,j,k,n) + dt * hydro_src(i,j,k,n) &
                                    + dt *  src(i,j,k,n) * a_half_inv

                  ! Momentum
                  else if (n .ge. UMX .and. n .le. UMZ) then
                     uout(i,j,k,n) = a_old * uin(i,j,k,n) + dt * hydro_src(i,j,k,n) &
                                    + dt   * src(i,j,k,n)
                     uout(i,j,k,n) = uout(i,j,k,n) * a_new_inv

                  ! (rho E)
                  else if (n .eq. UEDEN) then
                     uout(i,j,k,n) =  a_oldsq * uin(i,j,k,n) + dt * hydro_src(i,j,k,n) &
                                    + a_half  * dt * src(i,j,k,n)  
                     uout(i,j,k,n) = uout(i,j,k,n) * a_newsq_inv

                  ! (rho e)
                  else if (n .eq. UEINT) then

                     uout(i,j,k,n) =  a_oldsq*uin(i,j,k,n) + dt * hydro_src(i,j,k,n) &
                                    + a_half * dt * src(i,j,k,n) 

                     uout(i,j,k,n) = uout(i,j,k,n) * a_newsq_inv 

!                    We don't do this here because we are adding all of this term explicitly in hydro_src
!                    uout(i,j,k,n) = uout(i,j,k,n) * a_newsq_inv / &
!                        ( ONE + a_half * dt * (HALF * gamma_minus_1 * divu_cc(i,j,k)) * a_newsq_inv )

                  ! (rho X_i) and (rho adv_i) and (rho aux_i)
                  else
                     uout(i,j,k,n) = uin(i,j,k,n) +  dt * hydro_src(i,j,k,n) &
                                    + dt * src(i,j,k,n) * a_half_inv

                  endif

               enddo
            enddo
         enddo
      enddo

      ! Enforce the density >= small_dens.  Make sure we do this immediately after consup.
      call enforce_minimum_density(uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                                   uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                                   lo,hi,print_fortran_warnings)
      
      ! Enforce species >= 0
      call enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3, &
                                       uout_h1,uout_h2,uout_h3,lo,hi,0)

      ! Re-normalize the species
      if (normalize_species .eq. 1) then
         call normalize_new_species(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                                    lo,hi)
      end if

      end subroutine update_state

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine cmpflx(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                        flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
                        ugdnv,pgdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                        csml,c,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,print_fortran_warnings)

      use amrex_fort_module, only : rt => amrex_real
      use amrex_mempool_module, only : amrex_allocate, amrex_deallocate
      use amrex_constants_module
      use meth_params_module, only : QVAR, NVAR

      implicit none

      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
      integer pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3
      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer idir,ilo,ihi,jlo,jhi
      integer i,j,kc,kflux,k3d
      integer print_fortran_warnings

      real(rt) qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,NVAR)
      real(rt) ugdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
      real(rt) pgdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
      real(rt) csml(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      real(rt)    c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

      real(rt), pointer :: smallc(:,:),cavg(:,:)

      integer :: c_lo(2), c_hi(2)

      c_lo = [ilo-1, jlo-1]
      c_hi = [ihi+1, jhi+1]

      call amrex_allocate ( smallc, c_lo, c_hi )
      call amrex_allocate (   cavg, c_lo, c_hi )

      if(idir.eq.1) then
         do j = jlo, jhi
            do i = ilo, ihi
               smallc(i,j) = max( csml(i,j,k3d), csml(i-1,j,k3d) )
               cavg(i,j) = HALF*( c(i,j,k3d) + c(i-1,j,k3d) )
            enddo
         enddo
      elseif(idir.eq.2) then
         do j = jlo, jhi
            do i = ilo, ihi
               smallc(i,j) = max( csml(i,j,k3d), csml(i,j-1,k3d) )
               cavg(i,j) = HALF*( c(i,j,k3d) + c(i,j-1,k3d) )
            enddo
         enddo
      else
         do j = jlo, jhi
            do i = ilo, ihi
               smallc(i,j) = max( csml(i,j,k3d), csml(i,j,k3d-1) )
               cavg(i,j) = HALF*( c(i,j,k3d) + c(i,j,k3d-1) )
            enddo
         enddo
      endif

      ! Solve Riemann problem
      call riemannus(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                     cavg,smallc,ilo-1,jlo-1,ihi+1,jhi+1, &
                     flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
                     ugdnv,pgdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                     idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,print_fortran_warnings)

      call amrex_deallocate(smallc)
      call amrex_deallocate(cavg)

      end subroutine cmpflx

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine riemannus(ql,qr,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                           cav,smallc,gd_l1,gd_l2,gd_h1,gd_h2, &
                           uflx,uflx_l1,uflx_l2,uflx_l3,uflx_h1,uflx_h2,uflx_h3, &
                           ugdnv,pgdnv,pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3, &
                           idir,ilo,ihi,jlo,jhi,kc,kflux,k3d,print_fortran_warnings)

      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use network, only : nspec, naux
      use amrex_constants_module
      use prob_params_module, only : physbc_lo, Symmetry
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QW, QPRES, QREINT, QFA, QFS, &
                                     URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFA, UFS, &
                                     nadv, small_dens, small_pres, gamma_const, gamma_minus_1
      use analriem_module

      implicit none
      real(rt), parameter:: small = 1.d-8
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer gd_l1,gd_l2,gd_h1,gd_h2
      integer uflx_l1,uflx_l2,uflx_l3,uflx_h1,uflx_h2,uflx_h3
      integer pg_l1,pg_l2,pg_l3,pg_h1,pg_h2,pg_h3
      integer idir,ilo,ihi,jlo,jhi
      integer i,j,kc,kflux,k3d
      integer print_fortran_warnings

      real(rt) ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt) qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
      real(rt)    cav(gd_l1:gd_h1,gd_l2:gd_h2)
      real(rt) smallc(gd_l1:gd_h1,gd_l2:gd_h2)
      real(rt) uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,uflx_l3:uflx_h3,NVAR)
      real(rt) ugdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)
      real(rt) pgdnv(pg_l1:pg_h1,pg_l2:pg_h2,pg_l3:pg_h3)

      integer n, nq
      integer iadv, ispec

      real(rt), dimension(ilo:ihi) :: rgdnv,v1gdnv,v2gdnv,regdnv,ustar
      real(rt), dimension(ilo:ihi) :: rl, ul, v1l, v2l, pl, rel
      real(rt), dimension(ilo:ihi) :: rr, ur, v1r, v2r, pr, rer
!     real(rt), dimension(ilo:ihi) :: wl, wr
      real(rt), dimension(ilo:ihi) :: rhoetot, scr
      real(rt), dimension(ilo:ihi) :: rstar, cstar, estar, pstar
      real(rt), dimension(ilo:ihi) :: ro, uo, po, reo, co, entho
      real(rt), dimension(ilo:ihi) :: sgnm, spin, spout, ushock, frac
      real(rt), dimension(ilo:ihi) :: wsmall, csmall,qavg

      do j = jlo, jhi

         rl = ql(ilo:ihi,j,kc,QRHO)

         do i = ilo, ihi
            if (rl(i) .lt. ZERO) then
               if (print_fortran_warnings .gt. 0) then
                  print *,'... setting NEG RL IN RIEMANN:IDIR ',idir,i,j,k3d,rl(i),' to ',small_dens
                  call flush(6)
               endif
               rl(i) = max(rl(i),small_dens)
            else if (rl(i) .lt. small_dens) then

               if (print_fortran_warnings .gt. 0) then
                  print *,'... setting     RL IN RIEMANN:IDIR ',idir,i,j,k3d,rl(i),' to ',small_dens
                  call flush(6)
               endif
               rl(i) = max(rl(i),small_dens)
            end if
         end do

         ! pick left velocities based on direction
         if (idir.eq.1) then
            ul  = ql(ilo:ihi,j,kc,QU)
            v1l = ql(ilo:ihi,j,kc,QV)
            v2l = ql(ilo:ihi,j,kc,QW)
         elseif (idir.eq.2) then
            ul  = ql(ilo:ihi,j,kc,QV)
            v1l = ql(ilo:ihi,j,kc,QU)
            v2l = ql(ilo:ihi,j,kc,QW)
         else
            ul  = ql(ilo:ihi,j,kc,QW)
            v1l = ql(ilo:ihi,j,kc,QU)
            v2l = ql(ilo:ihi,j,kc,QV)
         endif

         pl  = ql(ilo:ihi,j,kc,QPRES)
         rel = ql(ilo:ihi,j,kc,QREINT)

         do i = ilo, ihi
            if (ql(i,j,kc,QPRES) .lt. small_pres) then
               if (print_fortran_warnings .gt. 0) then
                  print *,'... setting PL/REL IN RIEMANN:IDIR ',idir,i,j,k3d, &
                          ql(i,j,kc,QPRES),' to ',small_pres
                  call flush(6)
               endif
               pl(i)  = max(pl(i),small_pres)
               rel(i) = pl(i) / gamma_minus_1
            end if
         end do

         rr = qr(ilo:ihi,j,kc,QRHO)

         do i = ilo, ihi
            if (rr(i) .lt. ZERO) then
               !
               ! A critical region since we usually can't write from threads.
               !
               if (print_fortran_warnings .gt. 0) then
                  print *,'... setting NEG RR IN RIEMANN:IDIR ',idir,i,j,k3d,rr(i),' to ',small_dens
                  call flush(6)
               endif
               rr(i) = max(rr(i),small_dens)
            else if (rr(i) .lt. small_dens) then
               !
               ! A critical region since we usually can't write from threads.
               !

               if (print_fortran_warnings .gt. 0) then
                  print *,'... setting     RR IN RIEMANN:IDIR ',idir,i,j,k3d,rr(i),' to ',small_dens
                  call flush(6)
               endif
               rr(i) = max(rr(i),small_dens)
            end if
         end do

         ! pick right velocities based on direction
         if (idir.eq.1) then
            ur  = qr(ilo:ihi,j,kc,QU)
            v1r = qr(ilo:ihi,j,kc,QV)
            v2r = qr(ilo:ihi,j,kc,QW)
         elseif (idir.eq.2) then
            ur  = qr(ilo:ihi,j,kc,QV)
            v1r = qr(ilo:ihi,j,kc,QU)
            v2r = qr(ilo:ihi,j,kc,QW)
         else
            ur  = qr(ilo:ihi,j,kc,QW)
            v1r = qr(ilo:ihi,j,kc,QU)
            v2r = qr(ilo:ihi,j,kc,QV)
         endif

         pr  = qr(ilo:ihi,j,kc,QPRES)
         rer = qr(ilo:ihi,j,kc,QREINT)

         do i = ilo, ihi
            if (pr(i) .lt. small_pres) then
               if (print_fortran_warnings .gt. 0) then
                 print *,'... setting PR/RER IN RIEMANN:IDIR ',idir,i,j,k3d,qr(i,j,kc,QPRES),' to ',small_pres
                 call flush(6)
               endif
               pr(i) = max(pr(i),small_pres)
               rer(i) = pr(i) / gamma_minus_1
            end if
         end do

         csmall = smallc(ilo:ihi,j)
         wsmall = small_dens*csmall

         ! We keep these here in case we want to use these instead of calling the analytic solver
         ! wl = max(wsmall,sqrt(abs(gamma_const*pl*rl)))
         ! wr = max(wsmall,sqrt(abs(gamma_const*pr*rr)))
         ! pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
         ! ustar = ((wl*ul + wr*ur) +       (pl - pr))/(wl + wr)

         ! Call analytic Riemann solver
         call analriem(ilo,ihi, &
                       gamma_const, &
                       pl(ilo:ihi), &
                       rl(ilo:ihi), &
                       ul(ilo:ihi), &
                       pr(ilo:ihi), &
                       rr(ilo:ihi), &
                       ur(ilo:ihi), &
                       small_pres, &
                       pstar(ilo:ihi), &
                       ustar(ilo:ihi))

         ! This loop has more conditions and won't vectorize. But it's cheap
         ! compared to the loop that calls the analytic Riemann solver.

         do i = ilo, ihi

            if (ustar(i) .gt. ZERO) then
               ro(i) = rl(i)
               uo(i) = ul(i)
               po(i) = pl(i)
               reo(i) = rel(i)
            else if (ustar(i) .lt. ZERO) then
               ro(i) = rr(i)
               uo(i) = ur(i)
               po(i) = pr(i)
               reo(i) = rer(i)
            else
               ro(i) = HALF*(rl(i)+rr(i))
               uo(i) = HALF*(ul(i)+ur(i))
               po(i) = HALF*(pl(i)+pr(i))
               reo(i) = HALF*(rel(i)+rer(i))
            endif

         end do

         do i = ilo, ihi

            if (ro(i) .lt. small_dens) then
               !
               ! A critical region since we usually can't write from threads.
               !
               if (print_fortran_warnings .gt. 0) then
                  print *,'... setting RO     IN RIEMANN:IDIR ',idir,i,j,k3d,ro(i),' to ',small_dens
                  call flush(6)
               endif
               ro(i) = max(ro(i),small_dens)
            end if

         end do

         co = sqrt(abs(gamma_const*po/ro))
         co = max(csmall,co)
         entho = ((reo + po)/ro)/(co*co)

         rstar = ro  + (pstar - po)/(co*co)
         estar = reo + (pstar - po)*entho

         do i = ilo, ihi
            if (rstar(i) .lt. small_dens) then
               !
               ! A critical region since we usually can't write from threads.
               !
               if (print_fortran_warnings .gt. 0) then
                  print *,'... setting RSTAR  IN RIEMANN:IDIR ',idir,i,j,k3d,rstar(i),' to ',small_dens
                  call flush(6)
               endif
               rstar(i) = max(rstar(i),small_dens)
            end if
         end do

         cstar = sqrt(abs(gamma_const*pstar/rstar))
         cstar = max(cstar,csmall)

         sgnm = sign(ONE,ustar)
         spout = co - sgnm*uo
         spin = cstar - sgnm*ustar
         ushock = HALF*(spin + spout)

         do i = ilo, ihi
            if (pstar(i)-po(i) .ge. ZERO) then
               spin(i) = ushock(i)
               spout(i) = ushock(i)
            endif
            if (spout(i)-spin(i) .eq. ZERO) then
               scr(i) = small*cav(i,j)
            else
               scr(i) = spout(i)-spin(i)
            endif
         end do

         frac = (ONE + (spout + spin)/scr)*HALF
         frac = max(ZERO,min(ONE,frac))

         do i = ilo, ihi
            if (ustar(i) .gt. ZERO) then
               v1gdnv(i) = v1l(i)
               v2gdnv(i) = v2l(i)
            else if (ustar(i) .lt. ZERO) then
               v1gdnv(i) = v1r(i)
               v2gdnv(i) = v2r(i)
            else
               v1gdnv(i) = HALF*(v1l(i)+v1r(i))
               v2gdnv(i) = HALF*(v2l(i)+v2r(i))
            endif
         end do

         rgdnv = frac*rstar + (ONE - frac)*ro

         do i = ilo, ihi
            if (rgdnv(i) .lt. small_dens) then
               !
               ! A critical region since we usually can't write from threads.
               !
               print *,'SMALL RGDNV IN RIEMANN ',idir,i,j,k3d,rgdnv(i)
               print *,'LEFT ',ul(i),rl(i),pl(i)
               print *,'RGHT ',ur(i),rr(i),pr(i)
               call amrex_error("Error:: Nyx_advection_3d.f90 :: riemannus")
            end if
         end do

         ugdnv(ilo:ihi,j,kc) = frac*ustar + (ONE - frac)*uo
         pgdnv(ilo:ihi,j,kc) = frac*pstar + (ONE - frac)*po

         regdnv = frac*estar + (ONE - frac)*reo

         do i = ilo, ihi
            if (spout(i) .lt. ZERO) then
               rgdnv(i) = ro(i)
               ugdnv(i,j,kc) = uo(i)
               pgdnv(i,j,kc) = po(i)
               regdnv(i) = reo(i)
            endif
            if (spin(i) .ge. ZERO) then
               rgdnv(i) = rstar(i)
               ugdnv(i,j,kc) = ustar(i)
               pgdnv(i,j,kc) = pstar(i)
               regdnv(i) = estar(i)
            endif
         end do

         do i = ilo, ihi
            if (pgdnv(i,j,kc) .lt. small_pres) then
               !
               ! A critical region since we usually can't write from threads.
               !
               print *,'SMALL P ',i,j,k3d,pgdnv(i,j,kc)
               print *,'WITH IDIR ',idir
               print *,'PSTAR PO ',pstar(i), po
               print *,'SPIN SPOUT ',spin, spout
               print *,'FRAC ',frac
               print *,'LEFT ',ul(i),rl(i),pl(i)
               print *,'RGHT ',ur(i),rr(i),pr(i)
               call amrex_error("Error:: Nyx_advection_3d.f90 :: riemannus")
            end if
         end do

         pgdnv(ilo:ihi,j,kc) = max(pgdnv(ilo:ihi,j,kc),small_pres)

         ! NOTE: Here we assume constant gamma.
         regdnv        = pgdnv(ilo:ihi,j,kc) / gamma_minus_1

         do i = ilo, ihi
            ! Enforce that fluxes through a symmetry plane are hard zero.
             if (i     .eq. 0 .and. physbc_lo(1) .eq. Symmetry .and. idir .eq. 1) &
                  ugdnv(i,j,kc) = ZERO
             if (j     .eq. 0 .and. physbc_lo(2) .eq. Symmetry .and. idir .eq. 2) &
                  ugdnv(i,j,kc) = ZERO
             if (kflux .eq. 0 .and. physbc_lo(3) .eq. Symmetry .and. idir .eq. 3) &
                  ugdnv(i,j,kc) = ZERO
         end do

         ! Compute fluxes, order as conserved state (not q)
         uflx(ilo:ihi,j,kflux,URHO) = rgdnv*ugdnv(ilo:ihi,j,kc)

         if (idir.eq.1) then
            uflx(ilo:ihi,j,kflux,UMX) = uflx(ilo:ihi,j,kflux,URHO)*ugdnv(ilo:ihi,j,kc) + pgdnv(ilo:ihi,j,kc)
            uflx(ilo:ihi,j,kflux,UMY) = uflx(ilo:ihi,j,kflux,URHO)*v1gdnv
            uflx(ilo:ihi,j,kflux,UMZ) = uflx(ilo:ihi,j,kflux,URHO)*v2gdnv
         elseif (idir.eq.2) then
            uflx(ilo:ihi,j,kflux,UMX) = uflx(ilo:ihi,j,kflux,URHO)*v1gdnv
            uflx(ilo:ihi,j,kflux,UMY) = uflx(ilo:ihi,j,kflux,URHO)*ugdnv(ilo:ihi,j,kc) + pgdnv(ilo:ihi,j,kc)
            uflx(ilo:ihi,j,kflux,UMZ) = uflx(ilo:ihi,j,kflux,URHO)*v2gdnv
         else
            uflx(ilo:ihi,j,kflux,UMX) = uflx(ilo:ihi,j,kflux,URHO)*v1gdnv
            uflx(ilo:ihi,j,kflux,UMY) = uflx(ilo:ihi,j,kflux,URHO)*v2gdnv
            uflx(ilo:ihi,j,kflux,UMZ) = uflx(ilo:ihi,j,kflux,URHO)*ugdnv(ilo:ihi,j,kc) + pgdnv(ilo:ihi,j,kc)
         endif

         rhoetot = regdnv + HALF*rgdnv*(ugdnv(ilo:ihi,j,kc)**2 + v1gdnv**2 + v2gdnv**2)

         uflx(ilo:ihi,j,kflux,UEDEN) = ugdnv(ilo:ihi,j,kc)*(rhoetot + pgdnv(ilo:ihi,j,kc))
         uflx(ilo:ihi,j,kflux,UEINT) = ugdnv(ilo:ihi,j,kc)*regdnv

         do iadv = 1, nadv
            n  = UFA + iadv - 1
            nq = QFA + iadv - 1
            do i = ilo, ihi
               if (ustar(i) .gt. ZERO) then
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
               else if (ustar(i) .lt. ZERO) then
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
               else
                  qavg(i) = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                  uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg(i)
               endif
            end do
         enddo

         if (UFS .gt. 0) then
            do ispec = 1, nspec+naux
               n  = UFS + ispec - 1
               nq = QFS + ispec - 1
               do i = ilo, ihi
                  if (ustar(i) .gt. ZERO) then
                     uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*ql(i,j,kc,nq)
                  else if (ustar(i) .lt. ZERO) then
                     uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qr(i,j,kc,nq)
                  else
                     qavg(i) = HALF * (ql(i,j,kc,nq) + qr(i,j,kc,nq))
                     uflx(i,j,kflux,n) = uflx(i,j,kflux,URHO)*qavg(i)
                  endif
               end do
            enddo
        end if ! UFS > 0
      enddo

      end subroutine riemannus

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine make_divu_nd(lo,hi,q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,dx,dy,dz, &
                              divu_nd,div_l1,div_l2,div_l3,div_h1,div_h2,div_h3)

      use amrex_fort_module, only : rt => amrex_real
      use amrex_constants_module
      use meth_params_module, only : QU, QV, QW

      implicit none

      integer          :: lo(3),hi(3)
      integer          :: q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
      integer          :: div_l1,div_l2,div_l3,div_h1,div_h2,div_h3
      real(rt) :: dx, dy, dz
      real(rt) :: divu_nd(div_l1:div_h1,div_l2:div_h2,div_l3:div_h3)
      real(rt) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,*)

      integer          :: i, j, k
      real(rt) :: ux, vy, wz, dxinv, dyinv, dzinv

      dxinv = ONE/dx
      dyinv = ONE/dy
      dzinv = ONE/dz

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)+1

               ux = ( q(i  ,j  ,k  ,QU) - q(i-1,j  ,k  ,QU) &
                    + q(i  ,j  ,k-1,QU) - q(i-1,j  ,k-1,QU) &
                    + q(i  ,j-1,k  ,QU) - q(i-1,j-1,k  ,QU) &
                    + q(i  ,j-1,k-1,QU) - q(i-1,j-1,k-1,QU) ) * dxinv

               vy = ( q(i  ,j  ,k  ,QV) - q(i  ,j-1,k  ,QV) &
                    + q(i  ,j  ,k-1,QV) - q(i  ,j-1,k-1,QV) &
                    + q(i-1,j  ,k  ,QV) - q(i-1,j-1,k  ,QV) &
                    + q(i-1,j  ,k-1,QV) - q(i-1,j-1,k-1,QV) ) * dyinv

               wz = ( q(i  ,j  ,k  ,QW) - q(i  ,j  ,k-1,QW) &
                    + q(i  ,j-1,k  ,QW) - q(i  ,j-1,k-1,QW) &
                    + q(i-1,j  ,k  ,QW) - q(i-1,j  ,k-1,QW) &
                    + q(i-1,j-1,k  ,QW) - q(i-1,j-1,k-1,QW) ) * dzinv

               divu_nd(i,j,k) = FOURTH*(ux + vy + wz)

            enddo
         enddo
      enddo

      end subroutine make_divu_nd

