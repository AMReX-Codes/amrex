! :::
! ::: ------------------------------------------------------------------
! :::

    subroutine fort_add_grav_source(lo,hi,&
                               uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                               uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                               grav, gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                               dt,a_old,a_new) &
                               bind(C, name="fort_add_grav_source")

      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
           UEDEN, grav_source_type

      implicit none

      integer lo(3), hi(3)
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
      integer  gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3

      real(rt)  uin( uin_l1: uin_h1, uin_l2: uin_h2, uin_l3: uin_h3,NVAR)
      real(rt) uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
      real(rt) grav(  gv_l1:  gv_h1,  gv_l2:  gv_h2,  gv_l3:  gv_h3,3)
      real(rt) dt
      real(rt) a_old, a_new

      real(rt) :: a_half, a_oldsq, a_newsq, a_newsq_inv
      real(rt) :: rho
      real(rt) :: SrU, SrV, SrW, SrE
      real(rt) :: rhoInv, dt_a_new
      real(rt) :: old_rhoeint, new_rhoeint, old_ke, new_ke
      integer          :: i, j, k

      a_half  = 0.5d0 * (a_old + a_new)
      a_oldsq = a_old * a_old
      a_newsq = a_new * a_new
      a_newsq_inv = 1.d0 / a_newsq

      dt_a_new    = dt / a_new

      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 3: Puts all gravitational work into KE, not (rho e)

      ! Add gravitational source terms
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ! **** Start Diagnostics ****
               old_ke = 0.5d0 * (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) / &
                                 uout(i,j,k,URHO) 
               old_rhoeint = uout(i,j,k,UEDEN) - old_ke
               ! ****   End Diagnostics ****

               rho    = uin(i,j,k,URHO)
               rhoInv = 1.0d0 / rho

               SrU = rho * grav(i,j,k,1)
               SrV = rho * grav(i,j,k,2)
               SrW = rho * grav(i,j,k,3)

               ! We use a_new here because we think of d/dt(a rho u) = ... + (rho g)
               uout(i,j,k,UMX)   = uout(i,j,k,UMX) + SrU * dt_a_new
               uout(i,j,k,UMY)   = uout(i,j,k,UMY) + SrV * dt_a_new
               uout(i,j,k,UMZ)   = uout(i,j,k,UMZ) + SrW * dt_a_new

               if (grav_source_type .eq. 1) then

                   ! This does work (in 1-d)
                   ! Src = rho u dot g, evaluated with all quantities at t^n
                   SrE = uin(i,j,k,UMX) * grav(i,j,k,1) + &
                         uin(i,j,k,UMY) * grav(i,j,k,2) + &
                         uin(i,j,k,UMZ) * grav(i,j,k,3)
                   uout(i,j,k,UEDEN) = (a_newsq*uout(i,j,k,UEDEN) + SrE * (dt*a_half)) * a_newsq_inv

               else if (grav_source_type .eq. 3) then

                   new_ke = 0.5d0 * (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) / &
                                     uout(i,j,k,URHO) 
                   uout(i,j,k,UEDEN) = old_rhoeint + new_ke
               else 
                  call amrex_error("Error:: Nyx_advection_3d.f90 :: bogus grav_source_type")
               end if

            enddo
         enddo
      enddo

      end subroutine fort_add_grav_source
