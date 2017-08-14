module advection_module implicit none private public umeth3d, consup contains
! ::: ---------------------------------------------------------------
! ::: :: hyp_gam_3d     Compute hyperbolic fluxes using unsplit second
! ::: ::               order Godunov integrator. 
! ::: ::
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
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
! ::: ----------------------------------------------------------------

  subroutine hyp_mol_gam_eb_3d(q, qd_lo, qd_hi, &
                     lo, hi, dx, dt, &
                     flux1, fd1_lo, fd1_hi, &
                     flux2, fd2_lo, fd2_hi, &
                     flux3, fd3_lo, fd3_hi, &
                     domlo, domhi )

    use mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : QVAR, NVAR, QPRES, QRHO, QU, QV, QW, QC

    use ebslope_module, only : ebslopex, ebslopey, ebslopez

    use eos_type_module
    use eos_module, only : eos_t, eos_input_rt, eos
    use bl_constants_module
    implicit none

    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: fd1_lo(3), fd1_hi(3)
    integer, intent(in) :: fd2_lo(3), fd2_hi(3)
    integer, intent(in) :: fd3_lo(3), fd3_hi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    double precision, intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),QVAR)

    double precision, intent(inout) :: flux1(fd1_lo(1):fd1_hi(1),fd1_lo(2):fd1_hi(2),fd1_lo(3):fd1_hi(3),NVAR)
    double precision, intent(inout) :: flux2(fd2_lo(1):fd2_hi(1),fd2_lo(2):fd2_hi(2),fd2_lo(3):fd2_hi(3),NVAR)
    double precision, intent(inout) :: flux3(fd3_lo(1):fd3_hi(1),fd3_lo(2):fd3_hi(2),fd3_lo(3):fd3_hi(3),NVAR)
    

    double precision, intent(in) :: dx(3), dt

    double precision :: qtempl(QVAR),qtempr(QVAR), fluxtemp(NVAR)

    double precision :: dxinv, dyinv, dzinv

    integer :: n
    integer :: i, j, k,  iwave, idim
    integer :: qtlo(3), qthi(3)

    ! Left and right state arrays (edge centered, cell centered)
    double precision, pointer :: dqxal(:,:,:,:), dqyal(:,:,:,:), dqzal(:,:,:,:)

    type (eos_t) :: eos_state


    call build(eos_state)

    
    qtlo(1) = lo(1) - 1 
    qtlo(2) = lo(2) - 1 
    qtlo(3) = lo(3) - 1 
    qthi(1) = hi(1) + 1 
    qthi(2) = hi(2) + 1 
    qthi(3) = hi(3) + 1 

    call bl_allocate ( qxal1, qt_lo(1), qt_hi(1), qt_lo(2), qt_hi(2), qt_lo(3), qt_hi(3), 1, QVAR)
    call bl_allocate ( qyal1, qt_lo(1), qt_hi(1), qt_lo(2), qt_hi(2), qt_lo(3), qt_hi(3), 1, QVAR)
    call bl_allocate ( qzal1, qt_lo(1), qt_hi(1), qt_lo(2), qt_hi(2), qt_lo(3), qt_hi(3), 1, QVAR)

    ! Local constants
    dxinv = ONE/dx(1)
    dyinv = ONE/dx(2)
    dzinv = ONE/dx(3)

    ! Initialize pdivu to zero


          ! Compute all slopes at kc (k3d)
       call ebslopex(q,qd_lo,qd_hi, &
                      dqxal,qt_lo,qt_hi, &
                      lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),QVAR)

       do k = lo(3), hi(3)
       do j = lo(2), hi(2)

          do i = lo(1), hi(1)+1


!  alphas   1,2,3 correspond to u-c, u, u+c repsectively, 4 and 5 are transverse velocities

! right eigenvectors are rho, u, p, v, w

!   in qtemp, 1 is rho, 2 is u, 3 is q , 4 is v and 5 is w

              cspeed = q(i-1,j,k,QC)
              qtempl(1) = q(i-1,j,k,QRHO) + 0.5d0 * ( (dqxal(i-1,j,k,1)+dqxal(i-1,j,k,3))/cspeed + dqxal(i-1,j,k,2))

              qtempl(2) = q(i-1,j,k,QU) + 0.5d0 * ( (dqxal(i-1,j,k,3)-dqxal(i-1,j,k,1))/q(i-1,j,k,QRHO))
     
              qtempl(3)=  q(i-1,j,k,QP) + 0.5d0 *  (dqxal(i-1,j,k,1)+dqxal(i-1,j,k,3))*cspeed 

              qtempl(4) = q(i-1,j,k,QV) + 0.5d0 * dqxal(i-1,j,k,4)
              qtempl(5) = q(i-1,j,k,Qw) + 0.5d0 * dqxal(i-1,j,k,5)

              cspeed = q(i,j,k,QC)
              qtempr(1) = q(i,j,k,QRHO) - 0.5d0 * ( (dqxal(i,j,k,1)+dqxal(i,j,k,3))/cspeed + dqxal(i,j,k,2))

              qtempr(2) = q(i,j,k,QU) - 0.5d0 * ( (dqxal(i,j,k,3)-dqxal(i,j,k,1))/q(i,j,k,QRHO))
     
              qtempr(3)=  q(i,j,k,QP) - 0.5d0 *  (dqxal(i,j,k,1)+dqxal(i,j,k,3))*cspeed 

              qtempr(4) = q(i,j,k,QV) - 0.5d0 * dqxal(i,j,k,4)

              qtempr(5) = q(i,j,k,Qw) - 0.5d0 *  dqxal(i,j,k,5)

              call analriem(gamma,qtempl, qtempr,  smallp,smallr,fluxtemp,  debug)

              flux1(i,j,k,URHO) = fluxtemp(1)
              flux1(i,j,k,UMX) = fluxtemp(2)
              flux1(i,j,k,UMY) = fluxtemp(4)
              flux1(i,j,k,UMZ) = fluxtemp(5)
              flux1(i,j,k,UEDEN) = fluxtemp(3)



          enddo

       enddo
       enddo

          ! Compute all slopes at kc (k3d)
       call ebslopey(q,qd_lo,qd_hi, &
                      dqyal,qt_lo,qt_hi, &
                      lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),QVAR)

       do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1

          do i = lo(1), hi(1)


!     1,2,3 correspond to u-c, u, u+c repsectively

! right eigenvectors are rho, v, p, u, w

!   in qtemp, 1 is rho, 2 is v, 3 is q , 4 is u and 5 is w


              cspeed = q(i,j-1,k,QC)
              qtempl(1) = q(i,j-1,k,QRHO) + 0.5d0 * ( (dqyal(i,j-1,k,1)+dqyal(i,j-1,k,3))/cspeed + dqyal(i,j-1,k,2))

              qtempl(2) = q(i,j-1,k,QV) + 0.5d0 * ( (dqyal(i,j-1,k,3)-dqyal(i,j-1,k,1))/q(i,j-1,k,QRHO))
     
              qtempl(3)=  q(i,j-1,k,QP) + 0.5d0 *  (dqyal(i,j-1,k,1)+dqyal(i,j-1,k,3))*cspeed 

              qtempl(4) = q(i,j-1,k,QU) + 0.5d0 * dqyal(i,j-1,k,4)
              qtempl(5) = q(i,j-1,k,Qw) + 0.5d0 * dqyal(i,j-1,k,5)

              cspeed = q(i,j,k,QC)
              qtempr(1) = q(i,j,k,QRHO) - 0.5d0 * ( (dqyal(i,j,k,1)+dqyal(i,j,k,3))/cspeed + dqyal(i,j,k,2))

              qtempr(2) = q(i,j,k,QV) - 0.5d0 * ( (dqyal(i,j,k,3)-dqyal(i,j,k,1))/q(i,j,k,QRHO))
     
              qtempr(3)=  q(i,j,k,QP) - 0.5d0 *  (dqyal(i,j,k,1)+dqyal(i,j,k,3))*cspeed 

              qtempr(4) = q(i,j,k,QU) - 0.5d0 * dqyal(i,j,k,4)

              qtempr(5) = q(i,j,k,Qw) - 0.5d0 *  dqyal(i,j,k,5)

              call analriem(gamma,qtempl, qtempr,  smallp,smallr,fluxtemp,  debug)


              flux2(i,j,k,URHO) = fluxtemp(1)
              flux2(i,j,k,UMX) = fluxtemp(4)
              flux2(i,j,k,UMY) = fluxtemp(2)
              flux2(i,j,k,UMZ) = fluxtemp(5)
              flux2(i,j,k,UEDEN) = fluxtemp(3)


          enddo

       enddo
       enddo

          ! Compute all slopes at kc (k3d)
       call ebslopez(q,qd_lo,qd_hi, &
                      dqxal,qt_lo,qt_hi, &
                      lo(1),lo(2),lo(3),hi(1),hi(2),hi(3),QVAR)

       do k = lo(3), hi(3)+1
       do j = lo(2), hi(2)

          do i = lo(1), hi(1)


!     1,2,3 correspond to u-c, u, u+c repsectively

! right eigenvectors are rho, v, p, u, w

!   in qtemp, 1 is rho, 2 is v, 3 is q , 4 is u and 5 is w


              cspeed = q(i,j,k-1,QC)
              qtempl(1) = q(i,j,k-1,QRHO) + 0.5d0 * ( (dqzal(i,j,k-1,1)+dqzal(i,j,k-1,3))/cspeed + dqzal(i,j,k-1,2))

              qtempl(2) = q(i,j,k-1,QW) + 0.5d0 * ( (dqzal(i,j,k-1,3)-dqzal(i,j,k-1,1))/q(i,j,k-1,QRHO))
    
              qtempl(3)=  q(i,j,k-1,QP) + 0.5d0 *  (dqzal(i,j,k-1,1)+dqzal(i,j,k-1,3))*cspeed 

              qtempl(4) = q(i,j,k-1,QU) + 0.5d0 * dqzal(i,j,k-1,4)
              qtempl(5) = q(i,j,k-1,QV) + 0.5d0 * dqzal(i,j,k-1,5)

              cspeed = q(i,j,k,QC)
              qtempr(1) = q(i,j,k,QRHO) - 0.5d0 * ( (dqqzl(i,j,k,1)+dqzal(i,j,k,3))/cspeed + dqzal(i,j,k,2))

              qtempr(2) = q(i,j,k,QW) - 0.5d0 * ( (dqzal(i,j,k,3)-dqzal(i,j,k,1))/q(i,j,k,QRHO))
     
              qtempr(3)=  q(i,j,k,QP) - 0.5d0 *  (dqzal(i,j,k,1)+dqzal(i,j,k,3))*cspeed 

              qtempr(4) = q(i,j,k,QU) - 0.5d0 * dqzal(i,j,k,4)

              qtempr(5) = q(i,j,k,QV) - 0.5d0 *  dqzal(i,j,k,5)

              call analriem(gamma,qtempl, qtempr,  smallp,smallr,fluxtemp,  debug)

              flux3(i,j,k,URHO) = fluxtemp(1)
              flux3(i,j,k,UMX) = fluxtemp(4)
              flux3(i,j,k,UMY) = fluxtemp(5)
              flux3(i,j,k,UMZ) = fluxtemp(2)
              flux3(i,j,k,UDEDN) = fluxtemp(3)

          enddo

       enddo
       enddo

       call bl_dealloate(dqxal)
       call bl_dealloate(dqyal)
       call bl_dealloate(dqzal)

  end subroutine hyp_gam_eb_3d

end module advection_module
