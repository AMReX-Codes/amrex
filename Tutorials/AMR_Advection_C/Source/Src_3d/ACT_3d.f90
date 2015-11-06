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
    
    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer          :: lo(3),hi(3)
    integer          :: flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
    integer          :: flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
    integer          :: flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
    double precision :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
    double precision :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
    double precision :: flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
    
    ! Local variables
    integer          :: i,j,k,n
    double precision :: sum,fac
    
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux1(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux1(i,j,k,URHO) / sum
             else
                fac = ONE
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
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux2(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux2(i,j,k,URHO) / sum
             else
                fac = ONE
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
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + flux3(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = flux3(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                flux3(i,j,k,n) = flux3(i,j,k,n) * fac
             end do
          end do
       end do
    end do

  end subroutine normalize_species_fluxes

! ::
! :: ----------------------------------------------------------
! ::

subroutine enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3, &
                                       uout_h1,uout_h2,uout_h3,lo,hi)

  use network, only : nspec
  use meth_params_module, only : NVAR, URHO, UFS
  use bl_constants_module
  
  implicit none
  
  integer          :: lo(3), hi(3)
  integer          :: uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
  double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
  
  ! Local variables
  integer          :: i,j,k,n
  integer          :: int_dom_spec
  logical          :: any_negative
  double precision :: dom_spec,x
  
  double precision, parameter :: eps = -1.0d-16
  
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           
           any_negative = .false.
           !
           ! First deal with tiny undershoots by just setting them to zero.
           !
           do n = UFS, UFS+nspec-1
              if (uout(i,j,k,n) .lt. ZERO) then
                 x = uout(i,j,k,n)/uout(i,j,k,URHO)
                 if (x .gt. eps) then
                    uout(i,j,k,n) = ZERO
                 else
                    any_negative = .true.
                 end if
              end if
           end do
           !
           ! We know there are one or more undershoots needing correction.
           !
           if (any_negative) then
              !
              ! Find the dominant species.
              !
              int_dom_spec = UFS
              dom_spec     = uout(i,j,k,int_dom_spec)
              
              do n = UFS,UFS+nspec-1
                 if (uout(i,j,k,n) .gt. dom_spec) then
                    dom_spec     = uout(i,j,k,n)
                    int_dom_spec = n
                 end if
              end do
              !
              ! Now take care of undershoots greater in magnitude than 1e-16.
              !
              do n = UFS, UFS+nspec-1
                 
                 if (uout(i,j,k,n) .lt. ZERO) then
                    
                    x = uout(i,j,k,n)/uout(i,j,k,URHO)
                    !
                    ! Here we only print the bigger negative values.
                    !
                    if (x .lt. -1.d-2) then
                       print *,'At cell (i,j,k) = ',i,j,k
                       print *,'...   Fixing negative species ',n           ,' with X = ',x
                       print *,'...     from dominant species ',int_dom_spec,' with X = ',&
                                uout(i,j,k,int_dom_spec) / uout(i,j,k,URHO)
                    end if

                    !
                    ! Take enough from the dominant species to fill the negative one.
                    !
                    uout(i,j,k,int_dom_spec) = uout(i,j,k,int_dom_spec) + uout(i,j,k,n)
                    !
                    ! Test that we didn't make the dominant species negative.
                    !
                    if (uout(i,j,k,int_dom_spec) .lt. ZERO) then 
                       print *,' Just made nth dominant species negative ',int_dom_spec-UFS+1,' at ',i,j,k 
                       print *,'We were fixing species ',n-UFS+1,' which had value ',x
                       print *,'Dominant species became ',uout(i,j,k,int_dom_spec) / uout(i,j,k,URHO)
                       call bl_error("Error:: Castro_3d.f90 :: ca_enforce_nonnegative_species")
                    end if
                    !
                    ! Now set the negative species to zero.
                    !
                    uout(i,j,k,n) = ZERO
                    
                 end if
                 
              enddo
           end if
        enddo
     enddo
  enddo
  
end subroutine enforce_nonnegative_species

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine normalize_new_species(u,u_l1,u_l2,u_l3,u_h1,u_h2,u_h3,lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
    double precision :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
    
    ! Local variables
    integer          :: i,j,k,n
    double precision :: fac,sum
    
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             sum = ZERO
             do n = UFS, UFS+nspec-1
                sum = sum + u(i,j,k,n)
             end do
             if (sum .ne. ZERO) then
                fac = u(i,j,k,URHO) / sum
             else
                fac = ONE
             end if
             do n = UFS, UFS+nspec-1
                u(i,j,k,n) = u(i,j,k,n) * fac
             end do
          end do
       end do
    end do
    
  end subroutine normalize_new_species
