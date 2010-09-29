
c-----------------------------------------------------------------------
c Unrolled indexing in these 3 routines uses the fact that each array
c has a border of width 1
c-----------------------------------------------------------------------

c Works for NODE-based data.
      subroutine hgip(
     & v0, v1, mask,
     &     regl0, regh0, regl1, regh1,
     & sum)
      integer regl0, regh0, regl1, regh1
      double precision v0(*)
      double precision v1(*)
      double precision mask(*)
      double precision sum
      integer i, idiff
c      do 10 i = 1, (regh0 - regl0 + 1) * (regh1 - regl1 + 1)
      idiff = regh0 - regl0 + 1
      do i = idiff + 2, idiff * (regh1 - regl1) - 1
         sum = sum + mask(i) * v0(i) * v1(i)
      end do
      end
c-----------------------------------------------------------------------
      subroutine hgcg1(
     & r, p, z, x, w, c, mask,
     &     regl0, regh0, regl1, regh1,
     & alpha, rho)
      integer regl0, regh0, regl1, regh1
      double precision r(*)
      double precision p(*)
      double precision z(*)
      double precision x(*)
      double precision w(*)
      double precision c(*)
      double precision mask(*)
      double precision alpha, rho
      integer i, idiff
c      do 10 i = 1, (regh0 - regl0 + 1) * (regh1 - regl1 + 1)
      idiff = regh0 - regl0 + 1
      do i = idiff + 2, idiff * (regh1 - regl1) - 1
         r(i) = r(i) - alpha * w(i)
         x(i) = x(i) + alpha * p(i)
         z(i) = r(i) * c(i)
         rho = rho + mask(i) * z(i) * r(i)
      end do
      end
c-----------------------------------------------------------------------
      subroutine hgcg2(p, z,
     &     regl0, regh0, regl1, regh1,
     & alpha)
      integer regl0, regh0, regl1, regh1
      double precision p(*)
      double precision z(*)
      double precision alpha
      integer i, idiff
c      do 10 i = 1, (regh0 - regl0 + 1) * (regh1 - regl1 + 1)
      idiff = regh0 - regl0 + 1
      do i = idiff + 2, idiff * (regh1 - regl1) - 1
         p(i) = alpha * p(i) + z(i)
      end do
      end

