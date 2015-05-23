! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine transxz(qm, qmo, qp, qpo, qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3, &
                         fx, fx_l1, fx_l2, fx_l3, fx_h1, fx_h2, fx_h3, &
                         fz, fz_l1, fz_l2, fz_l3, fz_h1, fz_h2, fz_h3, &
                         srcQ, src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                         hdt, hdtdx, hdtdz, &
                         ilo, ihi, jlo, jhi, klo, khi)

      use network, only : nspec
      use meth_params_module, only : QVAR, NVAR, QRHO, QFA, QFS, &
                                     URHO, UFA, UFS, nadv
      implicit none

      integer qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3
      integer gc_l1, gc_l2, gc_l3, gc_h1, gc_h2, gc_h3
      integer fx_l1, fx_l2, fx_l3, fx_h1, fx_h2, fx_h3
      integer fz_l1, fz_l2, fz_l3, fz_h1, fz_h2, fz_h3
      integer src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
      integer ilo, ihi, jlo, jhi, klo, khi

      double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fx(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
      double precision fz(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
      double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      double precision hdt, hdtdx, hdtdz

      integer i, j, k
      integer n, nq, iadv
      integer ispec
      double precision rr, rrnew, compo, compn

      ! Update density
      do k = klo, khi
      do j = jlo, jhi
      do i = ilo, ihi

          qpo(i,j  ,k,QRHO) = qp(i,j  ,k,QRHO) - hdtdx*(fx(i+1,j,k,URHO)-fx(i  ,j,k,URHO)) &
                                               - hdtdz*(fz(i,j,k+1,URHO)-fz(i  ,j,k,URHO))
          qmo(i,j+1,k,QRHO) = qm(i,j+1,k,QRHO) - hdtdx*(fx(i+1,j,k,URHO)-fx(i  ,j,k,URHO)) &
                                               - hdtdz*(fz(i,j,k+1,URHO)-fz(i  ,j,k,URHO))

      enddo
      enddo
      enddo

      do iadv = 1, nadv
          n  = UFA + iadv - 1
          nq = QFA + iadv - 1
          do k = klo, khi 
          do j = jlo, jhi 
          do i = ilo, ihi 

              rr = qp(i,j,k,QRHO)
              rrnew = rr - hdtdx*(fx(i+1,j,k,URHO)-fx(i,j,k,URHO)) &
                         - hdtdz*(fz(i,j,k+1,URHO)-fz(i,j,k,URHO))

              compo = rr*qp(i,j,k,nq)
              compn = compo - hdtdx*(fx(i+1,j,k,n)-fx(i,j,k,n)) &
                            - hdtdz*(fz(i,j,k+1,n)-fz(i,j,k,n))

              qpo(i,j,k,nq) = compn/rrnew + hdt*srcQ(i,j,k,nq)

              rr = qm(i,j+1,i,QRHO)
              rrnew = rr - hdtdx*(fx(i+1,j,k,URHO)-fx(i,j,k,URHO)) &
                         - hdtdz*(fz(i,j,k+1,URHO)-fz(i,j,k,URHO))

              compo = rr*qm(i,j+1,k,nq)
              compn = compo - hdtdx*(fx(i+1,j,k,n)-fx(i,j,k,n)) &
                            - hdtdz*(fz(i,j,k+1,n)-fz(i,j,k,n))

              qmo(i,j+1,k,nq) = compn/rrnew + hdt*srcQ(i,j,k,nq)

          enddo
          enddo
          enddo
      enddo

      do ispec = 1, nspec
          n  = UFS + ispec - 1
          nq = QFS + ispec - 1
          do k = klo, khi 
          do j = jlo, jhi 
          do i = ilo, ihi 

              rr = qp(i,j,k,QRHO)
              rrnew = rr - hdtdx*(fx(i+1,j,k,URHO)-fx(i,j,k,URHO)) &
                         - hdtdz*(fz(i,j,k+1,URHO)-fz(i,j,k,URHO))

              compo = rr*qp(i,j,k,nq)
              compn = compo - hdtdx*(fx(i+1,j,k,n)-fx(i,j,k,n)) &
                            - hdtdz*(fz(i,j,k+1,n)-fz(i,j,k,n))

              qpo(i,j,k,nq) = compn/rrnew + hdt*srcQ(i,j,k,nq)

              rr = qm(i,j+1,i,QRHO)
              rrnew = rr - hdtdx*(fx(i+1,j,k,URHO)-fx(i,j,k,URHO)) &
                         - hdtdz*(fz(i,j,k+1,URHO)-fz(i,j,k,URHO))

              compo = rr*qm(i,j+1,k,nq)
              compn = compo - hdtdx*(fx(i+1,j,k,n)-fx(i,j,k,n)) &
                            - hdtdz*(fz(i,j,k+1,n)-fz(i,j,k,n))

              qmo(i,j+1,k,nq) = compn/rrnew + hdt*srcQ(i,j,k,nq)

          enddo
          enddo
          enddo
      enddo

      end subroutine transxz

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine transyz(qm, qmo, qp, qpo, qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3, &
                         fy, fy_l1, fy_l2, fy_l3, fy_h1, fy_h2, fy_h3, &
                         fz, fz_l1, fz_l2, fz_l3, fz_h1, fz_h2, fz_h3, &
                         srcQ, src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                         hdt, hdtdy, hdtdz, &
                         ilo, ihi, jlo, jhi, klo, khi)

      use network, only : nspec
      use meth_params_module, only : QVAR, NVAR, QRHO, QFS, QFA, &
                                     URHO, UFA, UFS, nadv
      implicit none

      integer qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3
      integer gc_l1, gc_l2, gc_l3, gc_h1, gc_h2, gc_h3
      integer fy_l1, fy_l2, fy_l3, fy_h1, fy_h2, fy_h3
      integer fz_l1, fz_l2, fz_l3, fz_h1, fz_h2, fz_h3
      integer src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
      integer ilo, ihi, jlo, jhi, klo, khi

      double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fy(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
      double precision fz(fz_l1:fz_h1,fz_l2:fz_h2,fz_l3:fz_h3,NVAR)
      double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      double precision hdt, hdtdy, hdtdz

      integer i, j, k
      integer n, nq, iadv, ispec

      double precision rr,rrnew,compo,compn

      ! Update density -- this assumes no sources for density
      do k = klo, khi
      do j = jlo, jhi
      do i = ilo, ihi

          qpo(i  ,j,k,QRHO) = qp(i  ,j,k,QRHO) - hdtdy*(fy(i,j+1,k,URHO)-fy(i,j,k,URHO)) &
                                               - hdtdz*(fz(i,j,k+1,URHO)-fz(i,j,k,URHO))
          qmo(i+1,j,k,QRHO) = qm(i+1,j,k,QRHO) - hdtdy*(fy(i,j+1,k,URHO)-fy(i,j,k,URHO)) &
                                               - hdtdz*(fz(i,j,k+1,URHO)-fz(i,j,k,URHO))

      enddo
      enddo
      enddo

      do iadv = 1, nadv
          n  = UFA + iadv - 1
          nq = QFA + iadv - 1
          do k = klo, khi 
          do j = jlo, jhi 
          do i = ilo, ihi 

              rr = qp(i,j,k,QRHO)
              rrnew = rr - hdtdy*(fy(i,j+1,k,URHO)-fy(i,j,k,URHO)) - hdtdz*(fz(i,j,k+1,URHO)-fz(i,j,k,URHO)) 

              compo = rr*qp(i,j,k,nq)
              compn = compo - hdtdy*(fy(i,j+1,k,n)-fy(i,j,k,n)) - hdtdz*(fz(i,j,k+1,n)-fz(i,j,k,n)) 

              qpo(i,j,k,nq) = compn/rrnew + hdt*srcQ(i,j,k,nq)

              rr = qm(i+1,j,k,QRHO)
              rrnew = rr - hdtdy*(fy(i,j+1,k,URHO)-fy(i,j,k,URHO))  - hdtdz*(fz(i,j,k+1,URHO)-fz(i,j,k,URHO)) 

              compo = rr*qm(i+1,j,k,nq)
              compn = compo - hdtdy*(fy(i,j+1,k,n)-fy(i,j,k,n))  - hdtdz*(fz(i,j,k+1,n)-fz(i,j,k,n)) 

              qmo(i+1,j,k,nq) = compn/rrnew + hdt*srcQ(i,j,k,nq)

          enddo
          enddo
          enddo
      enddo

      do ispec = 1, nspec 
          n  = UFS + ispec - 1
          nq = QFS + ispec - 1
          do k = klo, khi 
          do j = jlo, jhi 
          do i = ilo, ihi 

              rr = qp(i,j,k,QRHO)
              rrnew = rr - hdtdy*(fy(i,j+1,k,URHO)-fy(i,j,k,URHO)) - &
                           hdtdz*(fz(i,j,k+1,URHO)-fz(i,j,k,URHO)) 

              compo = rr*qp(i,j,k,nq)
              compn = compo - hdtdy*(fy(i,j+1,k,n)-fy(i,j,k,n)) &
                            - hdtdz*(fz(i,j,k+1,n)-fz(i,j,k,n)) 

              qpo(i,j,k,nq) = compn/rrnew + hdt*srcQ(i,j,k,nq)

              rr = qm(i+1,j,k,QRHO)
              rrnew = rr - hdtdy*(fy(i,j+1,k,URHO)-fy(i,j,k,URHO)) - &
                           hdtdz*(fz(i,j,k+1,URHO)-fz(i,j,k,URHO)) 

              compo = rr*qm(i+1,j,k,nq)
              compn = compo - hdtdy*(fy(i,j+1,k,n)-fy(i,j,k,n)) &
                            - hdtdz*(fz(i,j,k+1,n)-fz(i,j,k,n)) 

              qmo(i+1,j,k,nq) = compn/rrnew + hdt*srcQ(i,j,k,nq)

          enddo
          enddo
          enddo
      enddo

      end subroutine transyz

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine transxy(qm, qmo, qp, qpo, qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3, &
                         fx, fx_l1, fx_l2, fx_l3, fx_h1, fx_h2, fx_h3, &
                         fy, fy_l1, fy_l2, fy_l3, fy_h1, fy_h2, fy_h3, &
                         srcQ, src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                         hdt, hdtdx, hdtdy, &
                         ilo, ihi, jlo, jhi, klo, khi)

      use network, only : nspec
      use meth_params_module, only : QVAR, NVAR, QRHO, QFA, QFS, &
                                     URHO, UFA, UFS, nadv
      implicit none

      integer qd_l1, qd_l2, qd_l3, qd_h1, qd_h2, qd_h3
      integer gc_l1, gc_l2, gc_l3, gc_h1, gc_h2, gc_h3
      integer fx_l1, fx_l2, fx_l3, fx_h1, fx_h2, fx_h3
      integer fy_l1, fy_l2, fy_l3, fy_h1, fy_h2, fy_h3
      integer src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
      integer ilo, ihi, jlo, jhi, klo, khi

      double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
      double precision fx(fx_l1:fx_h1,fx_l2:fx_h2,fx_l3:fx_h3,NVAR)
      double precision fy(fy_l1:fy_h1,fy_l2:fy_h2,fy_l3:fy_h3,NVAR)
      double precision srcQ(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,QVAR)
      double precision hdt, hdtdx, hdtdy

      integer i, j, k
      integer n, nq, iadv
      integer ispec
      double precision rr, rrnew, compo, compn

      ! Update density
      do k = klo, khi
      do j = jlo, jhi
      do i = ilo, ihi

          qpo(i,j  ,k,QRHO) = qp(i,j  ,k,QRHO) - hdtdx*(fx(i+1,j,k,URHO)-fx(i,j,k,URHO)) &
                                               - hdtdy*(fy(i,j+1,k,URHO)-fy(i,j,k,URHO))
          qmo(i,j+1,k,QRHO) = qm(i,j+1,k,QRHO) - hdtdx*(fx(i+1,j,k,URHO)-fx(i,j,k,URHO)) &
                                               - hdtdy*(fy(i,j+1,k,URHO)-fy(i,j,k,URHO))

      enddo
      enddo
      enddo

      do iadv = 1, nadv
          n  = UFA + iadv - 1
          nq = QFA + iadv - 1
          do k = klo, khi 
          do j = jlo, jhi 
          do i = ilo, ihi 

              rr = qp(i,j,k,QRHO)
              rrnew = rr - hdtdx*(fx(i+1,j,k,URHO)-fx(i,j,k,URHO)) &
                         - hdtdy*(fy(i,j+1,k,URHO)-fy(i,j,k,URHO))

              compo = rr*qp(i,j,k,nq)
              compn = compo - hdtdx*(fx(i+1,j,k,n)-fx(i,j,k,n)) &
                            - hdtdy*(fy(i,j+1,k,n)-fy(i,j,k,n))

              qpo(i,j,k,nq) = compn/rrnew + hdt*srcQ(i,j,k,nq)

              rr = qm(i,j+1,i,QRHO)
              rrnew = rr - hdtdx*(fx(i+1,j,k,URHO)-fx(i,j,k,URHO)) &
                         - hdtdy*(fy(i,j+1,k,URHO)-fy(i,j,k,URHO))

              compo = rr*qm(i,j+1,k,nq)
              compn = compo - hdtdx*(fx(i+1,j,k,n)-fx(i,j,k,n)) &
                            - hdtdy*(fy(i,j+1,k,n)-fy(i,j,k,n))

              qmo(i,j+1,k,nq) = compn/rrnew + hdt*srcQ(i,j,k,nq)

          enddo
          enddo
          enddo
      enddo

      do ispec = 1, nspec
          n  = UFS + ispec - 1
          nq = QFS + ispec - 1
          do k = klo, khi 
          do j = jlo, jhi 
          do i = ilo, ihi 

              rr = qp(i,j,k,QRHO)
              rrnew = rr - hdtdx*(fx(i+1,j,k,URHO)-fx(i,j,k,URHO)) &
                         - hdtdy*(fy(i,j+1,k,URHO)-fy(i,j,k,URHO))

              compo = rr*qp(i,j,k,nq)
              compn = compo - hdtdx*(fx(i+1,j,k,n)-fx(i,j,k,n)) &
                            - hdtdy*(fy(i,j+1,k,n)-fy(i,j,k,n))

              qpo(i,j,k,nq) = compn/rrnew + hdt*srcQ(i,j,k,nq)

              rr = qm(i,j+1,i,QRHO)
              rrnew = rr - hdtdx*(fx(i+1,j,k,URHO)-fx(i,j,k,URHO)) &
                         - hdtdy*(fy(i,j+1,k,URHO)-fy(i,j,k,URHO))

              compo = rr*qm(i,j+1,k,nq)
              compn = compo - hdtdx*(fx(i+1,j,k,n)-fx(i,j,k,n)) &
                            - hdtdy*(fy(i,j+1,k,n)-fy(i,j,k,n))

              qmo(i,j+1,k,nq) = compn/rrnew + hdt*srcQ(i,j,k,nq)

          enddo
          enddo
          enddo
      enddo

      end subroutine transxy
