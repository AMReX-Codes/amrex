
      subroutine apply_metric(lo, hi, &
                              rhs, rl1, rl2, rh1, rh2,  &
                              ecx, ecxl1, ecxl2, ecxh1, ecxh2, &
                              ecy, ecyl1, ecyl2, ecyh1, ecyh2, dx, coord_type)

      implicit none
      integer lo(2),hi(2)
      integer rl1, rl2, rh1, rh2
      integer ecxl1, ecxl2, ecxh1, ecxh2
      integer ecyl1, ecyl2, ecyh1, ecyh2
      integer coord_type
      double precision rhs(rl1:rh1,rl2:rh2)
      double precision ecx(ecxl1:ecxh1,ecxl2:ecxh2)
      double precision ecy(ecyl1:ecyh1,ecyl2:ecyh2)
      double precision dx(2)

      double precision r
      integer i,j

      ! r-z
      if (coord_type .eq. 1) then

         ! At centers
         do i=lo(1),hi(1)
            r = (dble(i)+0.5d0) * dx(1)
            do j=lo(2),hi(2)
               rhs(i,j) = rhs(i,j) * r
            enddo
         enddo

         ! On x-edges
         do i=lo(1),hi(1)+1
            r = dble(i)*dx(1)
            do j=lo(2),hi(2)
               ecx(i,j) = ecx(i,j) * r
            enddo
         enddo

         ! On y-edges
         do i=lo(1),hi(1)
            r = (dble(i)+0.5d0) * dx(1)
            do j=lo(2),hi(2)+1
               ecy(i,j) = ecy(i,j) * r
            enddo
         enddo

      else 
         print *,'Bogus coord_type in apply_metric ' ,coord_type
         call bl_error("Error:: MGutils_2d.f90 :: apply_metric")
      end if

      end subroutine apply_metric

!-----------------------------------------------------------------------

      subroutine unweight_cc(lo, hi, &
                             cc, cl1, cl2, ch1, ch2,  &
                             dx, coord_type)

      implicit none
      integer lo(2),hi(2)
      integer cl1, cl2, ch1, ch2
      integer coord_type
      double precision cc(cl1:ch1,cl2:ch2)
      double precision dx(2)

      double precision r
      integer i,j

      ! r-z
      if (coord_type .eq. 1) then

         ! At centers
         do i=lo(1),hi(1)
            r = (dble(i)+0.5d0) * dx(1)
            do j=lo(2),hi(2)
               cc(i,j) = cc(i,j) / r
            enddo
         enddo

      else 
         print *,'Bogus coord_type in unweight_cc ' ,coord_type
         call bl_error("Error:: MGutils_2d.f90 :: unweight_cc")
      end if

      end subroutine unweight_cc

!-----------------------------------------------------------------------

      subroutine unweight_edges(lo, hi, &
                                ecx, ecxl1, ecxl2, ecxh1, ecxh2, &
                                ecy, ecyl1, ecyl2, ecyh1, ecyh2, dx, coord_type)

      implicit none
      integer lo(2),hi(2)
      integer ecxl1, ecxl2, ecxh1, ecxh2
      integer ecyl1, ecyl2, ecyh1, ecyh2
      integer coord_type
      double precision ecx(ecxl1:ecxh1,ecxl2:ecxh2)
      double precision ecy(ecyl1:ecyh1,ecyl2:ecyh2)
      double precision dx(2)

      double precision :: r
      integer          :: i,j

      ! r-z
      if (coord_type .eq. 1) then

         ! On x-edges
         do i = lo(1), hi(1)+1
            if (i .ne. 0) then
               r = dble(i)*dx(1)
               do j = lo(2),hi(2)
                  ecx(i,j) = ecx(i,j) / r
               enddo
            end if
         enddo

         ! On y-edges
         do i = lo(1), hi(1)
            r = (dble(i)+0.5d0) * dx(1)
            do j = lo(2),hi(2)+1
               ecy(i,j) = ecy(i,j) / r
            enddo
         enddo

      else 
         print *,'Bogus coord_type in unweight_edges ' ,coord_type
         call bl_error("Error:: MGutils_2d.f90 :: unweight_edges")
      end if

      end subroutine unweight_edges
