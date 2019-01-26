      subroutine fort_check_initial_species(lo,hi,&
                             state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3) &
      bind(C, name="fort_check_initial_species")

      use network           , only : nspec
      use meth_params_module, only : NVAR, URHO, UFS
      use  eos_params_module

      use amrex_fort_module, only : rt => amrex_real
      use amrex_error_module, only : amrex_abort
      implicit none

      integer  :: lo(3), hi(3)
      integer  :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
      real(rt) :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

      ! Local variables
      integer  :: i,j,k,n
      real(rt) :: sum
      character(len=256) :: errmsg_pt1, errmsg_pt2

      if (UFS .gt. 0) then

          !$OMP PARALLEL DO PRIVATE(i,j,k,n,sum)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
    
                   sum = 0.d0
                   do n = 1, nspec
                      sum = sum + state(i,j,k,UFS+n-1)
                   end do

                   if (abs(state(i,j,k,URHO)-sum).gt. 1.d-8 * state(i,j,k,URHO)) then
                      write(errmsg_pt1, *) 'Sum of (rho X)_i vs rho at (i,j,k): ', &
                        i,j,k,sum,state(i,j,k,URHO)
                      write(errmsg_pt2, *) trim(errmsg_pt1) // new_line('a') // &
                        'Failed check of initial species summing to 1'
                      call amrex_abort(errmsg_pt2)
                   end if
    
                enddo
             enddo
          enddo
      !$OMP END PARALLEL DO

      else
          sum = h_species + he_species
          if (abs(sum-1.d0).gt. 1.d-8) &
              call bl_error("Error:: Failed check of initial species summing to 1")
      end if

      end subroutine fort_check_initial_species
