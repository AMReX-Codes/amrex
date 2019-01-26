
! :::
! ::: ------------------------------------------------------------------
! :::

      !===========================================================================
      ! This is called from the C++ so the threading happens here...
      !===========================================================================
      subroutine fort_correct_gsrc(lo,hi, &
                              gold,gold_l1,gold_l2,gold_l3,gold_h1,gold_h2,gold_h3, &
                              gnew,gnew_l1,gnew_l2,gnew_l3,gnew_h1,gnew_h2,gnew_h3, &
                              uold,uold_l1,uold_l2,uold_l3,uold_h1,uold_h2,uold_h3, &
                              unew,unew_l1,unew_l2,unew_l3,unew_h1,unew_h2,unew_h3, &
                              a_old,a_new,dt) &
                              bind(C, name="fort_correct_gsrc")

      use amrex_error_module
      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, grav_source_type

      implicit none

      integer lo(3),hi(3)
      integer gold_l1,gold_l2,gold_l3,gold_h1,gold_h2,gold_h3
      integer gnew_l1,gnew_l2,gnew_l3,gnew_h1,gnew_h2,gnew_h3
      integer uold_l1,uold_l2,uold_l3,uold_h1,uold_h2,uold_h3
      integer unew_l1,unew_l2,unew_l3,unew_h1,unew_h2,unew_h3
      real(rt)   gold(gold_l1:gold_h1,gold_l2:gold_h2,gold_l3:gold_h3,3)
      real(rt)   gnew(gnew_l1:gnew_h1,gnew_l2:gnew_h2,gnew_l3:gnew_h3,3)
      real(rt)  uold(uold_l1:uold_h1,uold_l2:uold_h2,uold_l3:uold_h3,NVAR)
      real(rt)  unew(unew_l1:unew_h1,unew_l2:unew_h2,unew_l3:unew_h3,NVAR)
      real(rt)  a_old,a_new,dt

      integer i,j,k
      real(rt) SrU_old, SrV_old, SrW_old
      real(rt) SrU_new, SrV_new, SrW_new
      real(rt) SrUcorr, SrVcorr, SrWcorr, SrEcorr
      real(rt) rhoo, Upo, Vpo, Wpo
      real(rt) rhon, Upn, Vpn, Wpn

      real(rt) a_half, a_newsq, rhooinv, rhoninv, a_new_inv
      real(rt) old_ke, old_rhoeint
      real(rt) new_ke, new_rhoeint

      a_half    = 0.5d0 * (a_old + a_new)
      a_newsq   = a_new*a_new
      a_new_inv = 1.0d0 / a_new

      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 3: Puts all gravitational work into KE, not (rho e)

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ! **** Start Diagnostics ****
               old_ke = 0.5d0 * (unew(i,j,k,UMX)**2 + unew(i,j,k,UMY)**2 + unew(i,j,k,UMZ)**2) / &
                                 unew(i,j,k,URHO) 
               old_rhoeint = unew(i,j,k,UEDEN) - old_ke
               ! ****   End Diagnostics ****

               rhoo    = uold(i,j,k,URHO)
               rhooinv = 1.0d0 / uold(i,j,k,URHO)
               Upo     = uold(i,j,k,UMX) * rhooinv
               Vpo     = uold(i,j,k,UMY) * rhooinv
               Wpo     = uold(i,j,k,UMZ) * rhooinv

               ! Define old source terms
               SrU_old = rhoo * gold(i,j,k,1)
               SrV_old = rhoo * gold(i,j,k,2)
               SrW_old = rhoo * gold(i,j,k,3)

               rhon    = unew(i,j,k,URHO)
               rhoninv = 1.0d0 / unew(i,j,k,URHO)
               Upn     = unew(i,j,k,UMX) * rhoninv
               Vpn     = unew(i,j,k,UMY) * rhoninv
               Wpn     = unew(i,j,k,UMZ) * rhoninv

               ! Define new source terms
               SrU_new = rhon * gnew(i,j,k,1)
               SrV_new = rhon * gnew(i,j,k,2)
               SrW_new = rhon * gnew(i,j,k,3)

               ! Define corrections to source terms
               SrUcorr = 0.5d0*(SrU_new - SrU_old)
               SrVcorr = 0.5d0*(SrV_new - SrV_old)
               SrWcorr = 0.5d0*(SrW_new - SrW_old)

               ! This does work (in 1-d)
               if (grav_source_type .eq. 1) then
                   SrEcorr =  0.5d0 * ( (SrU_new * Upn - SrU_old * Upo) + &
                                        (SrV_new * Vpn - SrV_old * Vpo) + &
                                        (SrW_new * Wpn - SrW_old * Wpo) )
               end if

               ! Correct state with correction terms
               unew(i,j,k,UMX)   = unew(i,j,k,UMX)   + SrUcorr*dt * a_new_inv
               unew(i,j,k,UMY)   = unew(i,j,k,UMY)   + SrVcorr*dt * a_new_inv
               unew(i,j,k,UMZ)   = unew(i,j,k,UMZ)   + SrWcorr*dt * a_new_inv

               if (grav_source_type .eq. 1) then
                   unew(i,j,k,UEDEN) = unew(i,j,k,UEDEN) + SrEcorr*dt * (a_half / a_newsq)
               else if (grav_source_type .eq. 3) then
                   new_ke = 0.5d0 * (unew(i,j,k,UMX)**2 + unew(i,j,k,UMY)**2 + unew(i,j,k,UMZ)**2) / &
                                     unew(i,j,k,URHO) 
                   unew(i,j,k,UEDEN) = old_rhoeint + new_ke
               else 
                  call amrex_error("Error:: Nyx_advection_3d.f90 :: bogus grav_source_type")
               end if

            enddo
         enddo
      enddo

      end subroutine fort_correct_gsrc
! :::
! ::: ------------------------------------------------------------------
! :::
      subroutine fort_syncgsrc(lo,hi, &
                              gphi,gphi_l1,gphi_l2,gphi_l3,gphi_h1,gphi_h2,gphi_h3, &
                              gdphi,gdphi_l1,gdphi_l2,gdphi_l3,gdphi_h1,gdphi_h2,gdphi_h3, &
                              state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                              dstate,dstate_l1,dstate_l2,dstate_l3, &
                              dstate_h1,dstate_h2,dstate_h3, &
                              sync_src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,a_new,dt) &
                              bind(C, name="fort_syncgsrc")

      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ

      implicit none
 
      integer lo(3),hi(3)
      integer gphi_l1,gphi_l2,gphi_l3,gphi_h1,gphi_h2,gphi_h3
      integer gdphi_l1,gdphi_l2,gdphi_l3,gdphi_h1,gdphi_h2,gdphi_h3
      integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
      integer dstate_l1,dstate_l2,dstate_l3,dstate_h1,dstate_h2,dstate_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      real(rt)   gphi(gphi_l1:gphi_h1,gphi_l2:gphi_h2,gphi_l3:gphi_h3,3)
      real(rt)  gdphi(gdphi_l1:gdphi_h1,gdphi_l2:gdphi_h2,gdphi_l3:gdphi_h3,3)
      real(rt)  state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
      real(rt) dstate(dstate_l1:dstate_h1,dstate_l2:dstate_h2,dstate_l3:dstate_h3,3+1)
      real(rt) sync_src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,3+1)
      real(rt) a_new,dt
 
      !    Note that dstate is drho and drhoU, state is the entire state, and src
      !    is S_rhoU and S_rhoE
 
      integer          :: i,j,k
      real(rt) :: rho_pre, rhoU_pre, rhoV_pre, rhoW_pre
      real(rt) :: gx, gy, gz, dgx, dgy, dgz, SrU, SrV, SrW, SrE, a_new_inv
 
      a_new_inv = 1.0d0 / a_new

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               rho_pre  = state(i,j,k,URHO) - dstate(i,j,k,1)
               rhoU_pre = state(i,j,k,UMX)  - dstate(i,j,k,2)
               rhoV_pre = state(i,j,k,UMY)  - dstate(i,j,k,3)
               rhoW_pre = state(i,j,k,UMZ)  - dstate(i,j,k,4)

               gx  = gphi(i,j,k,1)
               gy  = gphi(i,j,k,2)
               gz  = gphi(i,j,k,3)

               dgx = gdphi(i,j,k,1)
               dgy = gdphi(i,j,k,2)
               dgz = gdphi(i,j,k,3)

               SrU = dstate(i,j,k,1)*gx + rho_pre*dgx
               SrV = dstate(i,j,k,1)*gy + rho_pre*dgy
               SrW = dstate(i,j,k,1)*gz + rho_pre*dgz

               SrE = ( SrU * (rhoU_pre + (0.5d0*dt)*SrU) + &
                       SrV * (rhoV_pre + (0.5d0*dt)*SrV) + &
                       SrW * (rhoW_pre + (0.5d0*dt)*SrW) ) / rho_pre

               sync_src(i,j,k,1) = SrU * a_new_inv
               sync_src(i,j,k,2) = SrV * a_new_inv
               sync_src(i,j,k,3) = SrW * a_new_inv
               sync_src(i,j,k,4) = SrE * a_new_inv

            enddo
         enddo
      enddo

      end subroutine fort_syncgsrc
