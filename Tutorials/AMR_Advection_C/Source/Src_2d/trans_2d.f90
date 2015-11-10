      subroutine transx(qm, qmo, qp, qpo, qd_l1, qd_l2, qd_h1, qd_h2, &
                        fx, fx_l1, fx_l2, fx_h1, fx_h2, &
                        srcQ, src_l1, src_l2, src_h1, src_h2, &
                        hdt, hdtdx,  &
                        ilo, ihi, jlo, jhi)

      use network, only : nspec
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QFA, QFS, &
                                     URHO, UX, UY, UFA, UFS, nadv
      implicit none

      integer qd_l1, qd_l2, qd_h1, qd_h2
      integer gc_l1, gc_l2, gc_h1, gc_h2
      integer fx_l1, fx_l2, fx_h1, fx_h2
      integer src_l1, src_l2, src_h1, src_h2
      integer ilo, ihi, jlo, jhi

      double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision fx(fx_l1:fx_h1,fx_l2:fx_h2,NVAR)
      double precision srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
      double precision hdt, hdtdx

      integer i, j
      integer n, nq, iadv
      integer ispec
      double precision rr, rrnew, compo, compn

      ! Update density
      do j = jlo, jhi
          do i = ilo, ihi

              qpo(i,j  ,QRHO) = qp(i,j  ,QRHO) - hdtdx*(fx(i+1,j,URHO)-fx(i  ,j,URHO))
              qmo(i,j+1,QRHO) = qm(i,j+1,QRHO) - hdtdx*(fx(i+1,j,URHO)-fx(i  ,j,URHO))

          enddo
      enddo

      do iadv = 1, nadv
          n  = UFA + iadv - 1
          nq = QFA + iadv - 1
          do j = jlo, jhi 
              do i = ilo, ihi 

                  rr = qp(i,j,  QRHO)
                  rrnew = rr - hdtdx*(fx(i+1,j,URHO)-fx(i  ,j,URHO))

                  compo = rr*qp(i,j  ,nq)
                  compn = compo - hdtdx*(fx(i+1,j,n)-fx(i  ,j,n))

                  qpo(i,j  ,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

                  rr = qm(i,j+1,QRHO)
                  rrnew = rr - hdtdx*(fx(i+1,j,URHO)-fx(i  ,j,URHO))

                  compo = rr*qm(i,j+1,nq)
                  compn = compo - hdtdx*(fx(i+1,j,n)-fx(i  ,j,n))

                  qmo(i,j+1,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

              enddo
          enddo
      enddo

      do ispec = 1, nspec
          n  = UFS + ispec - 1
          nq = QFS + ispec - 1
          do j = jlo, jhi 
              do i = ilo, ihi 

                  rr = qp(i  ,j,QRHO)
                  rrnew = rr - hdtdx*(fx(i+1,j,URHO)-fx(i  ,j,URHO))

                  compo = rr*qp(i,j  ,nq)
                  compn = compo - hdtdx*(fx(i+1,j,n)-fx(i  ,j,n))

                  qpo(i,j  ,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

                  rr = qm(i,j+1,QRHO)
                  rrnew = rr - hdtdx*(fx(i+1,j,URHO)-fx(i  ,j,URHO))

                  compo  = rr*qm(i,j+1,nq)
                  compn = compo - hdtdx*(fx(i+1,j,n)-fx(i  ,j,n))

                  qmo(i,j+1,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

              enddo
          enddo
      enddo

      end subroutine transx

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine transy(qm, qmo, qp, qpo, qd_l1, qd_l2, qd_h1, qd_h2, &
                        fy,fy_l1,fy_l2,fy_h1,fy_h2, &
                        srcQ, src_l1, src_l2, src_h1, src_h2, &
                        hdt, hdtdy, ilo, ihi, jlo, jhi)

      use network, only : nspec
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QFA, QFS, &
                                     URHO, UX, UY, UFA, UFS, nadv
      implicit none

      integer qd_l1, qd_l2, qd_h1, qd_h2
      integer gc_l1, gc_l2, gc_h1, gc_h2
      integer fy_l1, fy_l2, fy_h1, fy_h2
      integer src_l1, src_l2, src_h1, src_h2
      integer ilo, ihi, jlo, jhi

      double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision fy(fy_l1:fy_h1,fy_l2:fy_h2,NVAR)
      double precision srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
      double precision hdt, hdtdy

      integer i, j
      integer n, nq, iadv, ispec

      double precision rr,rrnew,compo,compn

      ! Update density -- this assumes no sources for density
      do j = jlo, jhi
      do i = ilo, ihi

          if (i.eq.12 .and.j.eq.0) print *,'TRANS BEF ',qp(i,j,QRHO), qm(i+1,j,QRHO)
          if (i.eq.12 .and.j.eq.0) print *,'TRANS FY  ',fy(i,j,URHO),fy(i,j+1,URHO)
          qpo(i  ,j,QRHO) = qp(i  ,j,QRHO) - hdtdy*(fy(i,j+1,URHO)-fy(i,j  ,URHO))
          qmo(i+1,j,QRHO) = qm(i+1,j,QRHO) - hdtdy*(fy(i,j+1,URHO)-fy(i,j  ,URHO))
          if (i.eq.12 .and.j.eq.0) print *,'TRANS ',qpo(i,j,QRHO), qmo(i+1,j,QRHO)

      enddo
      enddo

      do iadv = 1, nadv
          n  = UFA + iadv - 1
          nq = QFA + iadv - 1
          do j = jlo, jhi 
          do i = ilo, ihi 

              rr = qp(i,j,QRHO)
              rrnew = rr - hdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 

              compo = rr*qp(i,j,nq)
              compn = compo - hdtdy*(fy(i,j+1,n)-fy(i,j,n)) 

              qpo(i,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

              rr = qm(i+1,j,QRHO)
              rrnew = rr - hdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 

              compo = rr*qm(i+1,j,nq)
              compn = compo - hdtdy*(fy(i,j+1,n)-fy(i,j,n)) 

              qmo(i+1,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

          enddo
          enddo
      enddo

      do ispec = 1, nspec 
          n  = UFS + ispec - 1
          nq = QFS + ispec - 1
          do j = jlo, jhi 
          do i = ilo, ihi 

              rr = qp(i,j,QRHO)
              rrnew = rr - hdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 

              compo = rr*qp(i,j,nq)
              compn = compo - hdtdy*(fy(i,j+1,n)-fy(i,j,n)) 

              qpo(i,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

              rr = qm(i+1,j,QRHO)
              rrnew = rr - hdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 

              compo = rr*qm(i+1,j,nq)
              compn = compo - hdtdy*(fy(i,j+1,n)-fy(i,j,n)) 

              qmo(i+1,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

          enddo
          enddo
      enddo

      end subroutine transy

