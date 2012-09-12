! :::
! ::: ----------------------------------------------------------------
! :::
 
      subroutine f_network_init()
 
        use network
 
        call network_init()
 
      end subroutine f_network_init

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine get_num_spec(nspec_out)

        use network, only : nspec

        implicit none 

        integer, intent(out) :: nspec_out

        nspec_out = nspec

      end subroutine get_num_spec

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine get_spec_names(spec_names,ispec,len)

        use network, only : nspec, short_spec_names

        implicit none 

        integer, intent(in   ) :: ispec
        integer, intent(inout) :: len
        integer, intent(inout) :: spec_names(len)

        ! Local variables
        integer   :: i

        len = len_trim(short_spec_names(ispec+1))

        do i = 1,len
           spec_names(i) = ichar(short_spec_names(ispec+1)(i:i))
        end do

      end subroutine get_spec_names


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
                                   FirstAdv,FirstSpec,numadv, &
                                   normalize_species_in)

        ! Passing data from C++ into f90

        use meth_params_module
        use network, only : nspec

        implicit none 
 
        integer, intent(in) :: dm
        integer, intent(in) :: Density, Xvel, FirstAdv, FirstSpec
        integer, intent(in) :: numadv
        integer, intent(in) :: normalize_species_in

        integer             :: QLAST

        iorder = 2 
        difmag = 0.1d0

        ! NTHERM: number of thermodynamic variables
        ! NVAR  : number of total variables in initial system
        ! dm refers to velocity components, '1' refers to rho
        NTHERM = dm + 1
        NVAR = NTHERM + nspec + numadv

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

        UFS   = FirstSpec + 1

        ! QTHERM: number of primitive variables
        ! QVAR  : number of total variables in primitive form

        QTHERM = NTHERM
        QVAR = QTHERM + nspec + numadv

        ! We use these to index into the state "Q"
        QRHO  = 1

        QU    = 2
        QLAST = 2

        if (dm .ge. 2) then
           QV    = 3
           QLAST = 3
        end if

        if (dm .eq. 3) then
           QW    = 4
           QLAST = 4
        end if

        if (numadv .ge. 1) then
          QFA = QTHERM + 1
          QFS = QFA + numadv
        else 
          QFA = 1
          QFS = QTHERM + 1
        end if

        normalize_species     = normalize_species_in

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

