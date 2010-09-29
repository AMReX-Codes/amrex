
c-----------------------------------------------------------------------
c Unrolled indexing in these 3 routines uses the fact that each array
c has a border of width 1
c-----------------------------------------------------------------------
c Works for NODE-based data.
      subroutine hgip(
     & v0, v1, mask,
     &     regl0,regh0,regl1,regh1,regl2,regh2,
     & sum)
      integer regl0,regh0,regl1,regh1,regl2,regh2
      double precision v0(*)
      double precision v1(*)
      double precision mask(*)
      double precision sum
      integer i, jdiff, kdiff
      jdiff = regh0 - regl0 + 1
      kdiff = (regh1 - regl1 + 1) * jdiff
!$omp parallel do reduction(+ : sum)
      do i = kdiff + jdiff + 2, kdiff * (regh2 - regl2) - jdiff - 1
         sum = sum + mask(i) * v0(i) * v1(i)
      end do
!$omp end parallel do
      end
c-----------------------------------------------------------------------
      subroutine hgcg1(
     & r, p, z, x, w, c, mask,
     &     regl0,regh0,regl1,regh1,regl2,regh2,
     & alpha, rho)
      integer regl0,regh0,regl1,regh1,regl2,regh2
      double precision r(*)
      double precision p(*)
      double precision z(*)
      double precision x(*)
      double precision w(*)
      double precision c(*)
      double precision mask(*)
      double precision alpha, rho
      integer i, jdiff, kdiff
      jdiff = regh0 - regl0 + 1
      kdiff = (regh1 - regl1 + 1) * jdiff
!$omp parallel do reduction(+ : rho)
      do i = kdiff + jdiff + 2, kdiff * (regh2 - regl2) - jdiff - 1
         r(i) = r(i) - alpha * w(i)
         x(i) = x(i) + alpha * p(i)
         z(i) = r(i) * c(i)
         rho = rho + mask(i) * z(i) * r(i)
      end do
!$omp end parallel do
      end
c-----------------------------------------------------------------------
      subroutine hgcg2(
     & p, z,
     &     regl0,regh0,regl1,regh1,regl2,regh2,
     & alpha)
      integer regl0,regh0,regl1,regh1,regl2,regh2
      double precision p(*)
      double precision z(*)
      double precision alpha
      integer i, jdiff, kdiff
      jdiff = regh0 - regl0 + 1
      kdiff = (regh1 - regl1 + 1) * jdiff
!$omp parallel do private(i)
      do i = kdiff + jdiff + 2, kdiff * (regh2 - regl2) - jdiff - 1
         p(i) = alpha * p(i) + z(i)
      end do
!$omp end parallel do
      end
