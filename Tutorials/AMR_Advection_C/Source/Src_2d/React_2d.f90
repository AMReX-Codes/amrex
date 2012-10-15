
       subroutine react_state(lo,hi, &
                              s_in ,so_l1,so_l2,so_h1,so_h2,&
                              s_out,sn_l1,sn_l2,sn_h1,sn_h2,&
                              time,dt_react)

      use network           , only : nspec
      use meth_params_module, only : NVAR, URHO, UX, UY, UFS
      use burner_module

      implicit none

      integer lo(2),hi(2)
      integer so_l1,so_h1,so_l2,so_h2
      integer sn_l1,sn_h1,sn_l2,sn_h2
      double precision s_in (so_l1:so_h1,so_l2:so_h2,NVAR)
      double precision s_out(sn_l1:sn_h1,sn_l2:sn_h2,NVAR)
      double precision time,dt_react

      integer          :: i,j,n
      double precision :: rho, rhoInv, u, v, T_in
      double precision :: x_in(nspec), x_out(nspec)

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)

           ! Make sure all variables are copied before reactions.
           s_out(i,j,1:NVAR)          = s_in(i,j,1:NVAR)

           rho           = s_in(i,j,URHO)
           rhoInv        = 1.d0 / rho

!          HACK -- how do we want to set the temperature?  Are the reaction rates
!                  temperature-sensitive?
           T_in          = 0.d0

           x_in(1:nspec) = s_in(i,j,UFS:UFS+nspec-1) * rhoInv

           call burner(rho, T_in, x_in, dt_react, time, x_out)

           ! Make sure that species emerge in the proper range: [0,1]
           do n = 1, nspec
             x_out(n) = max(min(x_out(n),1.d0),0.d0)
           end do

           s_out(i,j,UFS:UFS+nspec-1) = rho * x_out(1:nspec)

      end do
      end do

      end subroutine react_state

