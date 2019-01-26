      subroutine fort_init_e_from_t(state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ns, &
                                    diag,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3,nd, &
                                   lo,hi,a_old) bind(C, name="fort_init_e_from_t")


      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT, UEDEN, TEMP_COMP, NE_COMP
      use  eos_params_module

      implicit none

      integer , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,ns
      integer , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3,nd
      real(rt), intent(inout) :: state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,ns)
      real(rt), intent(in   ) ::  diag(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,nd)
      integer , intent(in   ) :: lo(3), hi(3)
      real(rt), intent(in   ) :: a_old

      ! Local variables
      real(rt) :: e, pres, T
      integer  :: i,j,k
      !
      ! Compute energy from the EOS
      !
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               T  = diag(i,j,k,TEMP_COMP)

               ! Call EOS to get the internal energy
               call nyx_eos_given_RT(e, pres, state(i,j,k,URHO), diag(i,j,k,TEMP_COMP), &
                                     diag(i,j,k,NE_COMP), a_old)

               state(i,j,k,UEINT) = state(i,j,k,URHO) * e

               state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                  0.5d0 * (state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2) / state(i,j,k,URHO)

            enddo
         enddo
      enddo

      end subroutine fort_init_e_from_t
