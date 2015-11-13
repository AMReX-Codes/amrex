
! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine get_method_params(nGrowHyp)

        ! Passing data from f90 back to C++

        use meth_params_module

        implicit none 

        integer, intent(out) :: ngrowHyp

        nGrowHyp = NHYP

      end subroutine get_method_params

! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine set_method_params(dm,Density,Xvel, &
                                   FirstAdv,numadv)

        ! Passing data from C++ into f90

        use meth_params_module

        implicit none 
 
        integer, intent(in) :: dm
        integer, intent(in) :: Density, Xvel, FirstAdv
        integer, intent(in) :: numadv

        ! NTHERM: number of thermodynamic variables
        ! NVAR  : number of total variables in initial system
        ! dm refers to velocity components, '1' refers to rho
        NTHERM = dm + 1
        NVAR = NTHERM + numadv

        nadv = numadv

        ! We use these to index into the state "U"
        URHO  = Density   + 1
        UX    = Xvel      + 1
        UY   = UX + 1
        if (dm .eq. 3) UZ = UY + 1

        if (numadv .ge. 1) then
          UFA   = FirstAdv  + 1
        else 
          UFA = 1
        end if

      end subroutine set_method_params

! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine set_problem_params(dm,physbc_lo_in,physbc_hi_in,Outflow_in,Symmetry_in,coord_type_in)

        ! Passing data from C++ into f90

        use prob_params_module

        implicit none 
 
        integer, intent(in) :: dm
        integer, intent(in) :: physbc_lo_in(dm),physbc_hi_in(dm)
        integer, intent(in) :: Outflow_in
        integer, intent(in) :: Symmetry_in
        integer, intent(in) :: coord_type_in

        allocate(physbc_lo(dm))
        allocate(physbc_hi(dm))

        physbc_lo(:) = physbc_lo_in(:)
        physbc_hi(:) = physbc_hi_in(:)

        Outflow  = Outflow_in
        Symmetry = Symmetry_in

        coord_type = coord_type_in

      end subroutine set_problem_params

