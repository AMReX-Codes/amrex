! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_network_init() &
        bind(C, name="fort_network_init")

        use network

        call network_init()

      end subroutine fort_network_init

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_get_num_spec(nspec_out) &
        bind(C, name="fort_get_num_spec")

        use network, only : nspec

        implicit none

        integer, intent(out) :: nspec_out

        nspec_out = nspec

      end subroutine fort_get_num_spec

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_get_num_aux(naux_out) &
        bind(C, name="fort_get_num_aux")

        use network, only : naux

        implicit none

        integer, intent(out) :: naux_out

        naux_out = naux

      end subroutine fort_get_num_aux

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_get_spec_names(spec_names,ispec,len) &
        bind(C, name="fort_get_spec_names")

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

      end subroutine fort_get_spec_names

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_get_aux_names(aux_names,iaux,len) &
        bind(C, name="fort_get_aux_names")

        use network, only : naux, short_aux_names

        implicit none

        integer, intent(in   ) :: iaux
        integer, intent(inout) :: len
        integer, intent(inout) :: aux_names(len)

        ! Local variables
        integer   :: i

        len = len_trim(short_aux_names(iaux+1))

        do i = 1,len
           aux_names(i) = ichar(short_aux_names(iaux+1)(i:i))
        end do

      end subroutine fort_get_aux_names

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_get_method_params(nGrowHyp) &
        bind(C, name = "fort_get_method_params")

        ! Passing data from f90 back to C++

        use meth_params_module

        implicit none

        integer, intent(out) :: ngrowHyp

        nGrowHyp = NHYP

      end subroutine fort_get_method_params

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_set_xhydrogen(xhydrogen_in) &
        bind(C, name="fort_set_xhydrogen")

        use amrex_error_module
        use amrex_fort_module, only : rt => amrex_real
        use atomic_rates_module

        real(rt), intent(in) :: xhydrogen_in
        if (xhydrogen_in .lt. 0.d0 .or. xhydrogen_in .gt. 1.d0) &
            call amrex_error("Bad value of xhydrogen_in")

        XHYDROGEN = xhydrogen_in
        YHELIUM   = (1.0d0-XHYDROGEN)/(4.0d0*XHYDROGEN)

      end subroutine fort_set_xhydrogen

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_set_small_values(average_dens, average_temp, a, &
                                       small_dens_inout, small_temp_inout, small_pres_inout) &
        bind(C, name="fort_set_small_values")

        use amrex_fort_module, only : rt => amrex_real
        use meth_params_module
        use eos_module

        implicit none

        real(rt), intent(in   ) :: average_dens, average_temp
        real(rt), intent(in   ) :: a
        real(rt), intent(inout) :: small_dens_inout, small_temp_inout, &
                                           small_pres_inout

        ! Local variables
        real(rt) :: frac = 1.d-6
        real(rt) :: typical_Ne
        
        if (small_dens_inout .le. 0.d0) then
           small_dens = frac * average_dens
        else
           small_dens = small_dens_inout
        end if
        if (small_temp_inout .le. 0.d0) then
           small_temp = frac * average_temp
        else
           small_temp = small_temp_inout
        end if

        ! HACK HACK HACK -- FIX ME!!!!!
        typical_Ne = 1.d0

        call eos_init_small_pres(small_dens, small_temp, typical_Ne, small_pres, a)

        small_dens_inout = small_dens
        small_temp_inout = small_temp
        small_pres_inout = small_pres

      end subroutine fort_set_small_values

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_set_method_params( &
                 dm, numadv, ndiag_in, do_hydro, ppm_type_in, ppm_ref_in, &
                 ppm_flatten_before_integrals_in, &
                 use_colglaz_in, use_flattening_in, &
                 corner_coupling_in, version_2_in, &
                 use_const_species_in, gamma_in, normalize_species_in, &
                 heat_cool_in, inhomo_reion_in) &
                 bind(C, name = "fort_set_method_params")

        ! Passing data from C++ into f90

        use amrex_error_module
        use amrex_fort_module, only : rt => amrex_real
        use meth_params_module
        use  eos_params_module
        use atomic_rates_module
        use comoving_module, only : comoving_type
        use network, only : nspec, naux
        use eos_module

        implicit none

        integer,  intent(in) :: dm
        integer,  intent(in) :: numadv
        integer,  intent(in) :: ndiag_in
        integer,  intent(in) :: do_hydro
        integer,  intent(in) :: ppm_type_in
        integer,  intent(in) :: ppm_ref_in
        integer,  intent(in) :: ppm_flatten_before_integrals_in
        integer,  intent(in) :: use_colglaz_in
        integer,  intent(in) :: use_flattening_in
        integer,  intent(in) :: version_2_in
        integer,  intent(in) :: corner_coupling_in
        real(rt), intent(in) :: gamma_in
        integer,  intent(in) :: use_const_species_in
        integer,  intent(in) :: normalize_species_in
        integer,  intent(in) :: heat_cool_in
        integer,  intent(in) :: inhomo_reion_in

        integer             :: QNEXT
        integer             :: UNEXT

        integer             :: iadv, ispec

        use_const_species = use_const_species_in

        iorder = 2
        difmag = 0.1d0

        grav_source_type = 1
        ! We may want to default to 3 when using SDC because then the gravity updates
        ! do not change the internal energy (rho e) -- but this needs further investigation
        ! grav_source_type = 3

        comoving_type = 1

        NDIAG = ndiag_in

        if (do_hydro .eq. 0) then

           NVAR = 1
           URHO  = 1
           UMX   = -1
           UMY   = -1
           UMZ   = -1
           UEDEN = -1
           UEINT = -1
           UFA   = -1
           UFS   = -1

           TEMP_COMP = -1
             NE_COMP = -1
            ZHI_COMP = -1

        else

           TEMP_COMP = 1
             NE_COMP = 2
            if (inhomo_reion_in .gt. 0) then
               ZHI_COMP = 3
            else
               ZHI_COMP = -1
            endif

           !---------------------------------------------------------------------
           ! conserved state components
           !---------------------------------------------------------------------
    
           ! NTHERM: number of thermodynamic variables
           ! NVAR  : number of total variables in initial system
           ! dm refers to momentum components, '3' refers to (rho, rhoE, rhoe)
           NTHERM = dm + 3

           if (use_const_species .eq. 1) then
              if (nspec .ne. 2 .or. naux .ne. 0) then
                  call amrex_error("Bad nspec or naux in set_method_params")
              end if
              NVAR = NTHERM + numadv
           else
              NVAR = NTHERM + nspec + naux + numadv
           end if

           nadv = numadv

           ! We use these to index into the state "U"
           URHO  = 1
           UMX   = 2
           UMY   = 3
           UMZ   = 4
           UEDEN = 5
           UEINT = 6
           UNEXT = 7

           UFA   = -1
           UFS   = -1
           if (numadv .ge. 1) then
               UFA = UNEXT 
               if (use_const_species .eq. 0) then
                   UFS = UFA + numadv
               end if
           else
             if (use_const_species .eq. 0) then
                 UFS = UNEXT
             end if
           end if

           !---------------------------------------------------------------------
           ! primitive state components
           !---------------------------------------------------------------------

           ! IMPORTANT: if use_const_species = 0, then we assume that 
           !   the auxiliary quantities immediately follow the species
           !   so we can loop over species and auxiliary quantities.
   
           ! QTHERM: number of primitive variables, which includes pressure (+1) 
           !         but not big E (-1) 
           ! QVAR  : number of total variables in primitive form

           QTHERM = NTHERM
           if (use_const_species .eq. 1) then
              QVAR = QTHERM + numadv
           else
              QVAR = QTHERM + nspec + naux + numadv
           end if

           ! We use these to index into the state "Q"
           QRHO   = 1   ! rho
           QU     = 2   ! u
           QV     = 3   ! v
           QW     = 4   ! w
           QPRES  = 5   ! p
           QREINT = 6   ! (rho e)

           QNEXT  = QREINT+1
   
           QFS = -1
           if (numadv .ge. 1) then
             QFA = QNEXT
             if (use_const_species .eq. 0) &
                 QFS = QFA + numadv
           else
             QFA = -1
             if (use_const_species .eq. 0) &
                 QFS = QNEXT
           end if

           ! constant ratio of specific heats
           if (gamma_in .gt. 0.d0) then
              gamma_const = gamma_in
           else
              gamma_const = 5.d0/3.d0
           end if
           gamma_minus_1 = gamma_const - 1.d0

           ppm_type                     = ppm_type_in
           ppm_reference                = ppm_ref_in
           ppm_flatten_before_integrals = ppm_flatten_before_integrals_in
           use_colglaz                  = use_colglaz_in
           use_flattening               = use_flattening_in
           version_2                    = version_2_in
           corner_coupling              = corner_coupling_in
           normalize_species            = normalize_species_in

           heat_cool_type               = heat_cool_in
           inhomo_reion                 = inhomo_reion_in

           ! Easy indexing for the passively advected quantities.  
           ! This lets us loop over all four groups (advected, species, aux)
           ! in a single loop.
           allocate(qpass_map(QVAR))
           allocate(upass_map(NVAR))
           npassive = 0
           do iadv = 1, nadv
              upass_map(npassive + iadv) = UFA + iadv - 1
              qpass_map(npassive + iadv) = QFA + iadv - 1
           enddo
           npassive = npassive + nadv
           if(QFS > -1) then
              do ispec = 1, nspec+naux
                 upass_map(npassive + ispec) = UFS + ispec - 1
                 qpass_map(npassive + ispec) = QFS + ispec - 1
              enddo
              npassive = npassive + nspec + naux
           endif

        end if

      end subroutine fort_set_method_params

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_set_eos_params(h_species_in, he_species_in) &
        bind(C, name="fort_set_eos_params")

        ! Passing data from C++ into f90

        use amrex_fort_module, only : rt => amrex_real
        use  eos_params_module

        implicit none

        real(rt), intent(in) :: h_species_in, he_species_in

         h_species =  h_species_in
        he_species = he_species_in

      end subroutine fort_set_eos_params

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_set_problem_params( &
        dm,physbc_lo_in,physbc_hi_in,Outflow_in,Symmetry_in,coord_type_in) &
        bind(C, name="fort_set_problem_params")

        ! Passing data from C++ into f90

        use prob_params_module

        implicit none

        integer, intent(in) :: dm
        integer, intent(in) :: physbc_lo_in(dm),physbc_hi_in(dm)
        integer, intent(in) :: Outflow_in
        integer, intent(in) :: Symmetry_in
        integer, intent(in) :: coord_type_in

        physbc_lo(1:dm) = physbc_lo_in(1:dm)
        physbc_hi(1:dm) = physbc_hi_in(1:dm)

        Outflow  = Outflow_in
        Symmetry = Symmetry_in

        coord_type = coord_type_in

      end subroutine fort_set_problem_params

! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_init_zhi(lo, hi, &
           nd, diag_eos,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
           ratio, zhi, z_l1, z_l2, z_l3, z_h1, z_h2, z_h3) &
           bind(C, name="fort_init_zhi")

        use amrex_fort_module,  only : rt => amrex_real
        use meth_params_module, only : ZHI_COMP

        implicit none

        integer ratio, nd
        integer lo(3), hi(3)
        integer d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
        integer z_l1,z_l2,z_l3,z_h1,z_h2,z_h3
        real(rt) diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,nd)
        real(rt)      zhi(z_l1:z_h1,z_l2:z_h2,z_l3:z_h3)

        integer i, j, k, ic, jc, kc

        if (ZHI_COMP .gt. -1) then

           do k = lo(3), hi(3)
              kc = k / ratio
              do j = lo(2), hi(2)
                 jc = j / ratio
                 do i = lo(1), hi(1)
                    ic = i / ratio
                    diag_eos(i,j,k, ZHI_COMP) = zhi(ic, jc, kc)
                 end do
              end do
           end do

        end if

   end subroutine fort_init_zhi


! :::
! ::: ----------------------------------------------------------------
! :::

      ! Get density component number from state data MultiFab. Useful
      ! for accessing MultiFabs in C++.

      function get_comp_urho() bind(C, name="get_comp_urho")
        use meth_params_module
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int) :: get_comp_urho
        get_comp_urho = URHO
      end function get_comp_urho

! :::
! ::: ----------------------------------------------------------------
! :::

      ! Get temperature component number from EOS data MultiFab. Useful
      ! for accessing MultiFabs in C++.

      function get_comp_temp() bind(C, name="get_comp_temp")
        use meth_params_module
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int) :: get_comp_temp
        get_comp_temp = TEMP_COMP
      end function get_comp_temp

! :::
! ::: ----------------------------------------------------------------
! :::

      ! Get internal energy component number from state data MultiFab.
      ! Useful for accessing MultiFabs in C++.

      function get_comp_e_int() bind(C, name="get_comp_e_int")
        use meth_params_module
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int) :: get_comp_e_int
        get_comp_e_int = UEINT
      end function get_comp_e_int
