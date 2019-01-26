      subroutine fort_init_e_from_rhoe(state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ns, &
                                       lo,hi,a_old) bind(C, name="fort_init_e_from_rhoe")

      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT, UEDEN
      use  eos_params_module

      implicit none

      integer  ,intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ns
      real(rt) ,intent(inout) :: state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ns)
      integer  ,intent(in   ) :: lo(3), hi(3)
      real(rt) ,intent(in   ) :: a_old

      ! Local variables
      integer          :: i,j,k
      !
      ! Compute energy from the EOS
      !
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                  0.5d0 * (state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2) / state(i,j,k,URHO)

            enddo
         enddo
      enddo

      end subroutine fort_init_e_from_rhoe
