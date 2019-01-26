! ::
! :: ----------------------------------------------------------
! ::

      subroutine normalize_species_fluxes(flux1,flux1_l1,flux1_l2,flux1_l3, &
                                          flux1_h1,flux1_h2,flux1_h3, &
                                          flux2,flux2_l1,flux2_l2,flux2_l3, &
                                          flux2_h1,flux2_h2,flux2_h3, &
                                          flux3,flux3_l1,flux3_l2,flux3_l3, &
                                          flux3_h1,flux3_h2,flux3_h3, &
                                          lo,hi)

      use amrex_fort_module, only : rt => amrex_real
      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(3),hi(3)
      integer          :: flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
      integer          :: flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
      integer          :: flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
      real(rt) :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
      real(rt) :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
      real(rt) :: flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)

      ! Local variables
      integer          :: i,j,k,n
      real(rt) :: sum,fac

      if (UFS .gt. 0) then

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)+1
                   sum = 0.d0
                   do n = UFS, UFS+nspec-1
                      sum = sum + flux1(i,j,k,n)
                   end do
                   if (sum .ne. 0.d0) then
                      fac = flux1(i,j,k,URHO) / sum
                   else
                      fac = 1.d0
                   end if
                   do n = UFS, UFS+nspec-1
                      flux1(i,j,k,n) = flux1(i,j,k,n) * fac
                   end do
                end do
             end do
          end do

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)+1
                do i = lo(1),hi(1)
                   sum = 0.d0
                   do n = UFS, UFS+nspec-1
                      sum = sum + flux2(i,j,k,n)
                   end do
                   if (sum .ne. 0.d0) then
                      fac = flux2(i,j,k,URHO) / sum
                   else
                      fac = 1.d0
                   end if
                   do n = UFS, UFS+nspec-1
                      flux2(i,j,k,n) = flux2(i,j,k,n) * fac
                   end do
                end do
             end do
          end do

          do k = lo(3),hi(3)+1
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   sum = 0.d0
                   do n = UFS, UFS+nspec-1
                       sum = sum + flux3(i,j,k,n)
                   end do
                   if (sum .ne. 0.d0) then
                      fac = flux3(i,j,k,URHO) / sum
                   else
                      fac = 1.d0
                   end if
                   do n = UFS, UFS+nspec-1
                      flux3(i,j,k,n) = flux3(i,j,k,n) * fac
                   end do
                end do
             end do
          end do

      end if

      end subroutine normalize_species_fluxes

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine normalize_new_species(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi)

      use amrex_fort_module, only : rt => amrex_real
      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
      real(rt) :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)

      ! Local variables
      integer          :: i,j,k,n
      real(rt) :: fac,sum

      if (UFS .gt. 0) then

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            sum = 0.d0
            do n = UFS, UFS+nspec-1
               sum = sum + u(i,j,k,n)
            end do
            if (sum .ne. 0.d0) then
               fac = u(i,j,k,URHO) / sum
            else
               fac = 1.d0
            end if
            do n = UFS, UFS+nspec-1
               u(i,j,k,n) = u(i,j,k,n) * fac
            end do
         end do
      end do
      end do

      end if ! UFS > 0

      end subroutine normalize_new_species

