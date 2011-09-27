#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <REAL.H>
#include "ArrayLim.H"
#include "LO_BCTYPES.H"

#include "HypreABec_F.H"

c-----------------------------------------------------------------------

      subroutine FORT_HPACOEF(mat, 
     @     a, a_l1, a_l2, a_l3, a_h1, a_h2, a_h3,
     @     lo, hi, sa)

      implicit none

      integer lo(3), hi(3)
      integer a_l1, a_l2, a_l3, a_h1, a_h2, a_h3
      double precision a(a_l1:a_h1, a_l2:a_h2, a_l3:a_h3)
      double precision mat(0:6,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      double precision sa

      integer i, j, k

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               mat(6,i,j,k) = sa * a(i,j,k)
            enddo
         enddo
      enddo

      end subroutine

c-----------------------------------------------------------------------

      subroutine FORT_HPBCOEF(mat, 
     @     b, b_l1, b_l2, b_l3, b_h1, b_h2, b_h3,
     @     lo, hi, sb, dx, idim)

      implicit none

      integer lo(3), hi(3), idim
      integer b_l1, b_l2, b_l3, b_h1, b_h2, b_h3
      double precision b(b_l1:b_h1, b_l2:b_h2, b_l3:b_h3)
      double precision mat(0:6,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      double precision sb, dx(3)

      integer i, j, k
      double precision fac

      fac = sb / dx(idim+1)**2 

      if (idim .eq. 0) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  mat(0,i,j,k) = - fac * b(i,j,k)
                  mat(3,i,j,k) = - fac * b(i+1,j,k)
                  mat(6,i,j,k) = mat(6,i,j,k) + fac * (b(i,j,k) + b(i+1,j,k))
               enddo
            enddo
         enddo
      else if (idim .eq. 1) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  mat(1,i,j,k) = - fac * b(i,j,k)
                  mat(4,i,j,k) = - fac * b(i,j+1,k)
                  mat(6,i,j,k) = mat(6,i,j,k) + fac * (b(i,j,k) + b(i,j+1,k))
               enddo
            enddo
         enddo
      else
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  mat(2,i,j,k) = - fac * b(i,j,k)
                  mat(5,i,j,k) = - fac * b(i,j,k+1)
                  mat(6,i,j,k) = mat(6,i,j,k) + fac * (b(i,j,k) + b(i,j,k+1))
               enddo
            enddo
         enddo
      endif

      end subroutine

c-----------------------------------------------------------------------

      subroutine FORT_HPBMAT3(mat, 
     @     b, b_l1, b_l2, b_l3, b_h1, b_h2, b_h3,
     @     lo, hi, sb, dx, cdir, bct, bcl,
     @     mask, m_l1, m_l2, m_l3, m_h1, m_h2, m_h3)

      implicit none

      integer lo(3), hi(3), cdir, bct
      integer b_l1, b_l2, b_l3, b_h1, b_h2, b_h3
      integer m_l1, m_l2, m_l3, m_h1, m_h2, m_h3
      double precision b(b_l1:b_h1, b_l2:b_h2, b_l3:b_h3)
      integer mask(m_l1:m_h1, m_l2:m_h2, m_l3:m_h3)
      double precision mat(0:6,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      double precision sb, dx(3), bcl

      integer i, j, k
      double precision fac, h, bfm, bfv

      if (cdir .eq. 0 .or. cdir .eq. 3) then
         h = dx(1)
      elseif (cdir .eq. 1 .or. cdir .eq. 4) then
         h = dx(2)
      else
         h = dx(3)
      endif
      fac = sb / (h**2)
      if (cdir .eq. 0) then
         i = lo(1)
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               if (mask(i-1,j,k) .gt. 0) then
                  if (bct .eq. LO_DIRICHLET) then
                     bfv = fac * h / (0.5d0 * h + bcl)
                     bfm = bfv * b(i,j,k)
                  else if (bct .eq. LO_NEUMANN) then
                     bfm = 0.d0
                  else
                     print *, "hpmat3: unsupported boundary type"
                     stop
                  endif
                  mat(0,i,j,k) = 0.d0
                  mat(6,i,j,k) = mat(6,i,j,k) + bfm - fac * b(i,j,k)
               endif
            enddo
         enddo
      else if (cdir .eq. 3) then
         i = hi(1)
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               if (mask(i+1,j,k) .gt. 0) then
                  if (bct .eq. LO_DIRICHLET) then
                     bfv = fac * h / (0.5d0 * h + bcl)
                     bfm = bfv * b(i+1,j,k)
                  else if (bct .eq. LO_NEUMANN) then
                     bfm = 0.d0
                  else
                     print *, "hpmat3: unsupported boundary type"
                     stop
                  endif
                  mat(3,i,j,k) = 0.d0
                  mat(6,i,j,k) = mat(6,i,j,k) + bfm - fac * b(i+1,j,k)
               endif
            enddo
         enddo
      else if (cdir .eq. 1) then
         j = lo(2)
         do k = lo(3), hi(3)
            do i = lo(1), hi(1)
               if (mask(i,j-1,k) .gt. 0) then
                  if (bct .eq. LO_DIRICHLET) then
                     bfv = fac * h / (0.5d0 * h + bcl)
                     bfm = bfv * b(i,j,k)
                  else if (bct .eq. LO_NEUMANN) then
                     bfm = 0.d0
                  else
                     print *, "hpmat3: unsupported boundary type"
                     stop
                  endif
                  mat(1,i,j,k) = 0.d0
                  mat(6,i,j,k) = mat(6,i,j,k) + bfm - fac * b(i,j,k)
               endif
            enddo
         enddo
      else if (cdir .eq. 4) then
         j = hi(2)
         do k = lo(3), hi(3)
            do i = lo(1), hi(1)
               if (mask(i,j+1,k) .gt. 0) then
                  if (bct .eq. LO_DIRICHLET) then
                     bfv = fac * h / (0.5d0 * h + bcl)
                     bfm = bfv * b(i,j+1,k)
                  else if (bct .eq. LO_NEUMANN) then
                     bfm = 0.d0
                  else
                     print *, "hpmat3: unsupported boundary type"
                     stop
                  endif
                  mat(4,i,j,k) = 0.d0
                  mat(6,i,j,k) = mat(6,i,j,k) + bfm - fac * b(i,j+1,k)
               endif
            enddo
         enddo
      else if (cdir .eq. 2) then
         k = lo(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (mask(i,j,k-1) .gt. 0) then
                  if (bct .eq. LO_DIRICHLET) then
                     bfv = fac * h / (0.5d0 * h + bcl)
                     bfm = bfv * b(i,j,k)
                  else if (bct .eq. LO_NEUMANN) then
                     bfm = 0.d0
                  else
                     print *, "hpmat3: unsupported boundary type"
                     stop
                  endif
                  mat(2,i,j,k) = 0.d0
                  mat(6,i,j,k) = mat(6,i,j,k) + bfm - fac * b(i,j,k)
               endif
            enddo
         enddo
      else if (cdir .eq. 5) then
         k = hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (mask(i,j,k+1) .gt. 0) then
                  if (bct .eq. LO_DIRICHLET) then
                     bfv = fac * h / (0.5d0 * h + bcl)
                     bfm = bfv * b(i,j,k+1)
                  else if (bct .eq. LO_NEUMANN) then
                     bfm = 0.d0
                  else
                     print *, "hpmat3: unsupported boundary type"
                     stop
                  endif
                  mat(5,i,j,k) = 0.d0
                  mat(6,i,j,k) = mat(6,i,j,k) + bfm - fac * b(i,j,k+1)
               endif
            enddo
         enddo
      else
         print *, "hpmat3: impossible face orientation"
      endif

      end subroutine

c-----------------------------------------------------------------------

      subroutine FORT_HPBVEC3(vec, 
     @     b, b_l1, b_l2, b_l3, b_h1, b_h2, b_h3,
     @     lo, hi, sb, dx, cdir, bct, bcl,
     @     mask, m_l1, m_l2, m_l3, m_h1, m_h2, m_h3, 
     @     bcv, bv_l1, bv_l2, bv_l3, bv_h1, bv_h2, bv_h3)

      implicit none

      integer lo(3), hi(3), cdir, bct
      integer b_l1, b_l2, b_l3, b_h1, b_h2, b_h3
      integer m_l1, m_l2, m_l3, m_h1, m_h2, m_h3
      integer bv_l1, bv_l2, bv_l3, bv_h1, bv_h2, bv_h3
      double precision b(b_l1:b_h1, b_l2:b_h2, b_l3:b_h3)
      integer mask(m_l1:m_h1, m_l2:m_h2, m_l3:m_h3)
      double precision bcv(bv_l1:bv_h1, bv_l2:bv_h2, bv_l3:bv_h3)
      double precision vec(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      double precision sb, dx(3), bcl

      integer i, j, k
      double precision h, bfv

      if (cdir .eq. 0 .or. cdir .eq. 3) then
         h = dx(1)
      elseif (cdir .eq. 1 .or. cdir .eq. 4) then
         h = dx(2)
      else
         h = dx(3)
      endif
      if (cdir .eq. 0) then
         i = lo(1)
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               if (mask(i-1,j,k) .gt. 0) then
                  if (bct .eq. LO_DIRICHLET) then
                     bfv = (sb / h) / (0.5d0 * h + bcl)
                     bfv = bfv * b(i,j,k)
                  else if (bct .eq. LO_NEUMANN) then
                     bfv = sb / h
                  else
                     print *, "hpbvec3: unsupported boundary type"
                     stop
                  endif
                  vec(i,j,k) = vec(i,j,k) + bfv * bcv(i-1,j,k)
               endif
            enddo
         enddo
      else if (cdir .eq. 3) then
         i = hi(1)
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               if (mask(i+1,j,k) .gt. 0) then
                  if (bct .eq. LO_DIRICHLET) then
                     bfv = (sb / h) / (0.5d0 * h + bcl)
                     bfv = bfv * b(i+1,j,k)
                  else if (bct .eq. LO_NEUMANN) then
                     bfv = sb / h
                  else
                     print *, "hpbvec3: unsupported boundary type"
                     stop
                  endif
                  vec(i,j,k) = vec(i,j,k) + bfv * bcv(i+1,j,k)
               endif
            enddo
         enddo
      else if (cdir .eq. 1) then
         j = lo(2)
         do k = lo(3), hi(3)
            do i = lo(1), hi(1)
               if (mask(i,j-1,k) .gt. 0) then
                  if (bct .eq. LO_DIRICHLET) then
                     bfv = (sb / h) / (0.5d0 * h + bcl)
                     bfv = bfv * b(i,j,k)
                  else if (bct .eq. LO_NEUMANN) then
                     bfv = sb / h
                  else
                     print *, "hpbvec3: unsupported boundary type"
                     stop
                  endif
                  vec(i,j,k) = vec(i,j,k) + bfv * bcv(i,j-1,k)
               endif
            enddo
         enddo
      else if (cdir .eq. 4) then
         j = hi(2)
         do k = lo(3), hi(3)
            do i = lo(1), hi(1)
               if (mask(i,j+1,k) .gt. 0) then
                  if (bct .eq. LO_DIRICHLET) then
                     bfv = (sb / h) / (0.5d0 * h + bcl)
                     bfv = bfv * b(i,j+1,k)
                  else if (bct .eq. LO_NEUMANN) then
                     bfv = sb / h
                  else
                     print *, "hpbvec3: unsupported boundary type"
                     stop
                  endif
                  vec(i,j,k) = vec(i,j,k) + bfv * bcv(i,j+1,k)
               endif
            enddo
         enddo
      else if (cdir .eq. 2) then
         k = lo(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (mask(i,j,k-1) .gt. 0) then
                  if (bct .eq. LO_DIRICHLET) then
                     bfv = (sb / h) / (0.5d0 * h + bcl)
                     bfv = bfv * b(i,j,k)
                  else if (bct .eq. LO_NEUMANN) then
                     bfv = sb / h
                  else
                     print *, "hpbvec3: unsupported boundary type"
                     stop
                  endif
                  vec(i,j,k) = vec(i,j,k) + bfv * bcv(i,j,k-1)
               endif
            enddo
         enddo
      else if (cdir .eq. 5) then
         k = hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (mask(i,j,k+1) .gt. 0) then
                  if (bct .eq. LO_DIRICHLET) then
                     bfv = (sb / h) / (0.5d0 * h + bcl)
                     bfv = bfv * b(i,j,k+1)
                  else if (bct .eq. LO_NEUMANN) then
                     bfv = sb / h
                  else
                     print *, "hpbvec3: unsupported boundary type"
                     stop
                  endif
                  vec(i,j,k) = vec(i,j,k) + bfv * bcv(i,j,k+1)
               endif
            enddo
         enddo
      else
         print *, "hpbvec3: impossible face orientation"
      endif

      end subroutine

c-----------------------------------------------------------------------
