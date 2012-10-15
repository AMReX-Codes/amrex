
! :: ----------------------------------------------------------
! :: Volume-weight average the fine grid data onto the coarse
! :: grid.  Overlap is given in coarse grid coordinates.
! ::
! :: INPUTS / OUTPUTS:
! ::  crse      <=  coarse grid data
! ::  clo,chi    => index limits of crse array interior
! ::  ngc        => number of ghost cells in coarse array
! ::  nvar	 => number of components in arrays
! ::  fine       => fine grid data
! ::  flo,fhi    => index limits of fine array interior
! ::  ngf        => number of ghost cells in fine array
! ::  rfine      => (ignore) used in 2-D RZ calc
! ::  lo,hi      => index limits of overlap (crse grid)
! ::  lrat       => refinement ratio
! ::
! :: NOTE:
! ::  Assumes all data cell centered
! :: ----------------------------------------------------------
! ::
      subroutine avgdown(crse,c_l1,c_l2,c_h1,c_h2,nvar, &
                         cv,cv_l1,cv_l2,cv_h1,cv_h2, &
                         fine,f_l1,f_l2,f_h1,f_h2, &
                         fv,fv_l1,fv_l2,fv_h1,fv_h2,lo,hi,lrat)
      implicit none
      integer c_l1,c_l2,c_h1,c_h2
      integer cv_l1,cv_l2,cv_h1,cv_h2
      integer f_l1,f_l2,f_h1,f_h2
      integer fv_l1,fv_l2,fv_h1,fv_h2
      integer lo(2), hi(2)
      integer nvar, lrat(2)
      double precision crse(c_l1:c_h1,c_l2:c_h2,nvar)
      double precision cv(cv_l1:cv_h1,cv_l2:cv_h2)
      double precision fine(f_l1:f_h1,f_l2:f_h2,nvar)
      double precision fv(fv_l1:fv_h1,fv_l2:fv_h2)

      integer i, j, n, ic, jc, ioff, joff
      integer lenx, leny, mxlen
      integer lratx, lraty

      lratx = lrat(1)
      lraty = lrat(2)
      lenx = hi(1)-lo(1)+1
      leny = hi(2)-lo(2)+1
      mxlen = max(lenx,leny)

      if (lenx .eq. mxlen) then
         do n = 1, nvar
 
!           Set coarse grid to zero on overlap
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,n) = 0.d0
               enddo
            enddo

!           Sum fine data
            do joff = 0, lraty-1
               do jc = lo(2), hi(2)
                  j = jc*lraty + joff
                  do ioff = 0, lratx-1
                     do ic = lo(1), hi(1)
                        i = ic*lratx + ioff
                        crse(ic,jc,n) = crse(ic,jc,n) + fv(i,j) * fine(i,j,n)
                     enddo
                  enddo
               enddo
            enddo

!           Divide out by volume weight
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,n) = crse(ic,jc,n) / cv(ic,jc)
               enddo
            enddo
            
         enddo

      else

         do n = 1, nvar

!           Set coarse grid to zero on overlap
            do ic = lo(1), hi(1)
               do jc = lo(2), hi(2)
                  crse(ic,jc,n) = 0.d0
               enddo
            enddo
 
!           Sum fine data
            do ioff = 0, lratx-1
               do ic = lo(1), hi(1)
                  i = ic*lratx + ioff
                  do joff = 0, lraty-1
                     do jc = lo(2), hi(2)
                        j = jc*lraty + joff
                        crse(ic,jc,n) = crse(ic,jc,n) + fv(i,j) * fine(i,j,n)
                     enddo
                  enddo
               enddo
            enddo
             
!           Divide out by volume weight
            do ic = lo(1), hi(1)
               do jc = lo(2), hi(2)
                  crse(ic,jc,n) = crse(ic,jc,n) / cv(ic,jc)
               enddo
            enddo
            
         enddo

      end if

      end subroutine avgdown

! ::
! :: ----------------------------------------------------------
! ::

      subroutine normalize_species_fluxes(  &
                        flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                        flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                        lo,hi)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(2),hi(2)
      integer          :: flux1_l1,flux1_l2,flux1_h1,flux1_h2
      integer          :: flux2_l1,flux2_l2,flux2_h1,flux2_h2
      double precision :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
      double precision :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)

      ! Local variables
      integer          :: i,j,n
      double precision :: sum,fac

      do j = lo(2),hi(2)
            do i = lo(1),hi(1)+1
               sum = 0.d0
               do n = UFS, UFS+nspec-1
                  sum = sum + flux1(i,j,n)
      	       end do
               if (sum .ne. 0.d0) then
                  fac = flux1(i,j,URHO) / sum
               else
                  fac = 1.d0
               end if
               do n = UFS, UFS+nspec-1
                  flux1(i,j,n) = flux1(i,j,n) * fac
      	       end do
            end do
      end do
      do j = lo(2),hi(2)+1
            do i = lo(1),hi(1)
               sum = 0.d0
               do n = UFS, UFS+nspec-1
                  sum = sum + flux2(i,j,n)
      	       end do
               if (sum .ne. 0.d0) then
                  fac = flux2(i,j,URHO) / sum
               else
                  fac = 1.d0
               end if
               do n = UFS, UFS+nspec-1
                  flux2(i,j,n) = flux2(i,j,n) * fac
      	       end do
            end do
      end do

      end subroutine normalize_species_fluxes

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(2), hi(2)
      integer          :: uout_l1,uout_l2,uout_h1,uout_h2
      double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)

      ! Local variables
      integer          :: i,j,n,nn
      integer          :: int_dom_spec
      logical          :: any_negative
      double precision :: dom_spec,x,eps

      eps = -1.0d-16

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

         any_negative = .false.

         ! First deal with tiny undershoots by just setting them to zero
         do n = UFS, UFS+nspec-1
           if (uout(i,j,n) .lt. 0.d0) then
              x = uout(i,j,n)/uout(i,j,URHO)
              if (x .gt. eps) then
                 uout(i,j,n) = 0.d0
              else
                 any_negative = .true.
              end if
           end if
         end do

         ! We know there are one or more undershoots needing correction 
         if (any_negative) then

            ! Find the dominant species
            int_dom_spec  = UFS
            dom_spec      = uout(i,j,int_dom_spec)

            do n = UFS,UFS+nspec-1
              if (uout(i,j,n) .gt. dom_spec) then
                dom_spec = uout(i,j,n)
                int_dom_spec = n
              end if
            end do

           ! Now take care of undershoots greater in magnitude than 1e-16.
           do n = UFS, UFS+nspec-1

              if (uout(i,j,n) .lt. 0.d0) then

                 x = uout(i,j,n)/uout(i,j,URHO)

                 ! Here we only print the bigger negative values
                 if (x .lt. -1.d-2) then
                    print *,'At cell (i,j) = ',i,j
                    print *,'... Fixing negative species ',n           ,' with X = ',x
                    print *,'...   from dominant species ',int_dom_spec,' with X = ',&
                             uout(i,j,int_dom_spec) / uout(i,j,URHO)
                 end if

                 ! Take enough from the dominant species to fill the negative one.
                 uout(i,j,int_dom_spec) = uout(i,j,int_dom_spec) + uout(i,j,n)
   
                 ! Test that we didn't make the dominant species negative
                 if (uout(i,j,int_dom_spec) .lt. 0.d0) then 
                    print *,'Just made dominant species negative ',int_dom_spec,' at ',i,j
                    print *,'... We were fixing species ',n,' which had value ',x
                    print *,'... Dominant species became ',uout(i,j,int_dom_spec) / uout(i,j,URHO)
                    call bl_error("Error:: ACTReact_2d.f90 :: enforce_nonnegative_species")
                 end if

                 ! Now the negative species to zero
                 uout(i,j,n) = 0.d0

              end if

           enddo
         end if

      enddo
      enddo

      ! Test again to make sure all species are between 0 and 1. 
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)

         any_negative = .false.

         do n = UFS, UFS+nspec-1
           if (uout(i,j,n) .lt. eps) then
              print *,'Species ',n,' still negative at (i,j) = ',i,j,uout(i,j,n)
              any_negative = .true.
           end if
           if ( uout(i,j,n) .gt. 1.d0 .and. abs(uout(i,j,n)-1.d0) .gt. 1.e-15) then
              print *,'Species ',n,' still overshoots at (i,j) = ',i,j,(uout(i,j,n)-1.d0)
              any_negative = .true.
           end if
         end do
      enddo
      enddo

      end subroutine enforce_nonnegative_species

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine normalize_new_species(u,u_l1,u_l2,u_h1,u_h2,lo,hi)

      use network, only : nspec
      use meth_params_module, only : NVAR, URHO, UFS

      implicit none

      integer          :: lo(2), hi(2)
      integer          :: u_l1,u_l2,u_h1,u_h2
      double precision :: u(u_l1:u_h1,u_l2:u_h2,NVAR)

      ! Local variables
      integer          :: i,j,n
      double precision :: fac,sum

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            sum = 0.d0
            do n = UFS, UFS+nspec-1
               sum = sum + u(i,j,n)
            end do
            if (sum .ne. 0.d0) then
               fac = u(i,j,URHO) / sum
            else
               fac = 1.d0
            end if
            do n = UFS, UFS+nspec-1
               u(i,j,n) = u(i,j,n) * fac
            end do
         end do
      end do

      end subroutine normalize_new_species

