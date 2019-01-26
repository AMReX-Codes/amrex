   subroutine fort_enforce_consistent_e(lo,hi,state, &
                                        state_l1,state_l2,state_l3,state_h1,state_h2,state_h3) & 
     bind(C,name="fort_enforce_consistent_e")

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT

     implicit none

     integer          :: lo(3), hi(3)
     integer          :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
     real(rt) :: state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)

     ! Local variables
     integer          :: i,j,k
     real(rt) :: u, v, w, rhoInv

     ! 
     ! Make sure to enforce (rho E) = (rho e) + 1/2 rho (u^2 +_ v^2 + w^2)
     !
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              rhoInv = 1.0d0 / state(i,j,k,URHO)

              u = state(i,j,k,UMX) * rhoInv
              v = state(i,j,k,UMY) * rhoInv
              w = state(i,j,k,UMZ) * rhoInv

              state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                     0.5d0 * state(i,j,k,URHO) * (u*u + v*v + w*w)

           end do
        end do
     end do

   end subroutine fort_enforce_consistent_e
