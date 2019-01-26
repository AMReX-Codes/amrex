
! ::: -----------------------------------------------------------

      subroutine hypfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                         adv_h3,domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="hypfill")

      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only : NVAR

      implicit none
      include 'AMReX_bc_types.fi'

      integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      real(rt) delta(3), xlo(3), time
      real(rt) adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3,NVAR)

      real(rt) state(NVAR)
      real(rt) staten(NVAR)

      integer i, j, k, n, lo(3), hi(3)
      real(rt) x, y, z
      logical rho_only

      lo(1) = adv_l1
      lo(2) = adv_l2
      lo(3) = adv_l3
      hi(1) = adv_h1
      hi(2) = adv_h2
      hi(3) = adv_h3

      do n = 1,NVAR
         call filcc(adv(adv_l1,adv_l2,adv_l3,n), &
              adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
              domlo,domhi,delta,xlo,bc(1,1,n))
      enddo

!     The strategy here is to set Dirichlet condition for inflow and
!     outflow boundaries, and let the Riemann solver sort out the
!     proper upwinding.  However, this decision makes this routine
!     look somewhat non-orthodox, in that we need to set external
!     values in either case....how do we know it's Outflow?  We have
!     to assume that the setup routines converted Outflow to FOEXTRAP.

!     Set flag for bc function
      rho_only = .FALSE.

!     XLO
      if ( (bc(1,1,1).eq.EXT_DIR.or.bc(1,1,1).eq.FOEXTRAP).and.adv_l1.lt.domlo(1)) then
         do i = adv_l1, domlo(1)-1
            do j = adv_l2, adv_h2
               do k = adv_l3, adv_h3
                  do n=1,NVAR
                     state(n) = adv(domlo(1),j,k,n)
                  enddo
                  call bcnormal(state,staten,1,+1,rho_only)
                  do n=1,NVAR
                     adv(i,j,k,n) = staten(n)
                  enddo
               end do
            end do
         end do
      end if            

!     XHI
      if ( (bc(1,2,1).eq.EXT_DIR.or.bc(1,2,1).eq.FOEXTRAP).and.adv_h1.gt.domhi(1)) then
         do i = domhi(1)+1, adv_h1
            do j = adv_l2, adv_h2
               do k = adv_l3, adv_h3
                  do n=1,NVAR
                     state(n) = adv(domhi(1),j,k,n)
                  enddo
                  call bcnormal(state,staten,1,-1,rho_only)
                  do n=1,NVAR
                     adv(i,j,k,n) = staten(n)
                  enddo
               end do
            end do
	 end do
      end if            

!     YLO
      if ( (bc(2,1,1).eq.EXT_DIR.or.bc(2,1,1).eq.FOEXTRAP).and.adv_l2.lt.domlo(2)) then
         do i = adv_l1, adv_h1
            do j = adv_l2, domlo(2)-1
               do k = adv_l3, adv_h3
                  do n=1,NVAR
                     state(n) = adv(i,domlo(2),k,n)
                  enddo
                  call bcnormal(state,staten,2,+1,rho_only)
                  do n=1,NVAR
                     adv(i,j,k,n) = staten(n)
                  enddo
               end do
            end do
         end do
      end if            

!     YHI
      if ( (bc(2,2,1).eq.EXT_DIR.or.bc(2,2,1).eq.FOEXTRAP).and.adv_h2.gt.domhi(2)) then
         do i = adv_l1, adv_h1
            do j = domhi(2)+1, adv_h2
               do k = adv_l3, adv_h3
                  do n=1,NVAR
                     state(n) = adv(i,domhi(2),k,n)
                  enddo
                  call bcnormal(state,staten,2,-1,rho_only)
                  do n=1,NVAR
                     adv(i,j,k,n) = staten(n)
                  enddo
               end do
            end do
	 end do
      end if            

!     ZLO
      if ( (bc(3,1,1).eq.EXT_DIR.or.bc(3,1,1).eq.FOEXTRAP).and.adv_l3.lt.domlo(3)) then
         do i = adv_l1, adv_h1
            do j = adv_l2, adv_h2
               do k = adv_l3, domlo(3)-1
                  do n=1,NVAR
                     state(n) = adv(i,j,domlo(3),n)
                  enddo
                  call bcnormal(state,staten,3,+1,rho_only)
                  do n=1,NVAR
                     adv(i,j,k,n) = staten(n)
                  enddo
               end do
            end do
         end do
      end if            

!     ZHI
      if ( (bc(3,2,1).eq.EXT_DIR.or.bc(3,2,1).eq.FOEXTRAP).and.adv_h3.gt.domhi(3)) then
         do i = adv_l1, adv_h1
            do j = adv_l2, adv_h2
               do k = domhi(3)+1, adv_h3
                  do n=1,NVAR
                     state(n) = adv(i,j,domhi(3),n)
                  enddo
                  call bcnormal(state,staten,3,-1,rho_only)
                  do n=1,NVAR
                     adv(i,j,k,n) = staten(n)
                  enddo
               end do
            end do
	 end do
      end if            

      end subroutine hypfill

! ::: -----------------------------------------------------------

      subroutine denfill(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2, &
                         adv_h3,domlo,domhi,delta,xlo,time,bc) &
                         bind(C, name="denfill")

      use amrex_fort_module, only : rt => amrex_real
      implicit none
      include 'AMReX_bc_types.fi'
      integer adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      real(rt) delta(3), xlo(3), time
      real(rt) adv(adv_l1:adv_h1,adv_l2:adv_h2,adv_l3:adv_h3)
      logical rho_only
      integer i,j,k

!     Note: this function should not be needed, technically, but is provided
!     to filpatch because there are many times in the algorithm when just
!     the density is needed.  We try to rig up the filling so that the same
!     function is called here and in hypfill where all the states are filled.

      call filcc(adv,adv_l1,adv_l2,adv_l3,adv_h1,adv_h2,adv_h3, &
           domlo,domhi,delta,xlo,bc)

      rho_only = .TRUE.

!     XLO
      if ( (bc(1,1,1).eq.EXT_DIR.or.bc(1,1,1).eq.FOEXTRAP).and.adv_l1.lt.domlo(1)) then
         do i = adv_l1, domlo(1)-1
            do j = adv_l2, adv_h2
               do k = adv_l3, adv_h3
                  call bcnormal(adv(domlo(1),j,k),adv(i,j,k),1,+1,rho_only)
               end do
            end do
         end do
      end if            

!     XHI
      if ( (bc(1,2,1).eq.EXT_DIR.or.bc(1,2,1).eq.FOEXTRAP).and.adv_h1.gt.domhi(1)) then
         do i = domhi(1)+1, adv_h1
            do j = adv_l2, adv_h2
               do k = adv_l3, adv_h3
                  call bcnormal(adv(domhi(1),j,k),adv(i,j,k),1,-1,rho_only)
               end do
            end do
	 end do
      end if            

!     YLO
      if ( (bc(2,1,1).eq.EXT_DIR.or.bc(2,1,1).eq.FOEXTRAP).and.adv_l2.lt.domlo(2)) then
         do i = adv_l1, adv_h1
            do j = adv_l2, domlo(2)-1
               do k = adv_l3, adv_h3
                  call bcnormal(adv(i,domlo(2),k),adv(i,j,k),2,+1,rho_only)
               end do
            end do
         end do
      end if            

!     YHI
      if ( (bc(2,2,1).eq.EXT_DIR.or.bc(2,2,1).eq.FOEXTRAP).and.adv_h2.gt.domhi(2)) then
         do i = adv_l1, adv_h1
            do j = domhi(2)+1, adv_h2
               do k = adv_l3, adv_h3
                  call bcnormal(adv(i,domhi(2),k),adv(i,j,k),2,-1,rho_only)
               end do
            end do
	 end do
      end if            

!     ZLO
      if ( (bc(3,1,1).eq.EXT_DIR.or.bc(3,1,1).eq.FOEXTRAP).and.adv_l3.lt.domlo(3)) then
         do i = adv_l1, adv_h1
            do j = adv_l2, adv_h2
               do k = adv_l3, domlo(3)-1
                  call bcnormal(adv(i,j,domlo(3)),adv(i,j,k),3,+1,rho_only)
               end do
            end do
         end do
      end if            

!     ZHI
      if ( (bc(3,2,1).eq.EXT_DIR.or.bc(3,2,1).eq.FOEXTRAP).and.adv_h3.gt.domhi(3)) then
         do i = adv_l1, adv_h1
            do j = adv_l2, adv_h2
               do k = domhi(3)+1, adv_h3
                  call bcnormal(adv(i,j,domhi(3)),adv(i,j,k),3,-1,rho_only)
               end do
            end do
	 end do
      end if            

      end subroutine denfill

! ::: -----------------------------------------------------------

      subroutine bcnormal(u_int,u_ext,dir,sgn,rho_only)

      use amrex_fort_module, only : rt => amrex_real
      use probdata_module
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, gamma_minus_1
      implicit none

      real(rt) u_int(*),u_ext(*)
      logical rho_only
      integer dir,sgn
      real(rt) rho, rhou(3), eden, T, Y
      integer n,t1,t2,i

!     for the Sedov problem, we will always set the state to the ambient conditions

      if (rho_only .EQV. .TRUE. ) then

         u_ext(1) = dens_ambient

      else

!     First set everything from internal data (this is probably a bad thing to do...)
!     That is, we should have explicit boundary data for advected fields and species

         do i=1,NVAR
            u_ext(i) = u_int(i)
         enddo

         u_ext(URHO)   = dens_ambient
         u_ext(UMX)    = 0.d0
         u_ext(UMY)    = 0.d0
         u_ext(UMZ)    = 0.d0
         u_ext(UEDEN)  = p_ambient/gamma_minus_1
         u_ext(UEINT)  = u_ext(UEDEN)

      endif

      end subroutine bcnormal
      
! ::: -----------------------------------------------------------

      subroutine generic_fill(var,var_l1,var_l2,var_l3,var_h1,var_h2,var_h3, &
                              domlo,domhi,delta,xlo,time,bc) &
                              bind(C, name="generic_fill")

      use amrex_fort_module, only : rt => amrex_real
      implicit none
      include 'AMReX_bc_types.fi'
      integer var_l1,var_l2,var_l3,var_h1,var_h2,var_h3
      integer bc(3,2,*)
      integer domlo(3), domhi(3)
      real(rt) delta(3), xlo(3), time
      real(rt) var(var_l1:var_h1,var_l2:var_h2,var_l3:var_h3)

      call filcc(var,var_l1,var_l2,var_l3,var_h1,var_h2,var_h3,domlo,domhi,delta,xlo,bc)

      end subroutine generic_fill

