module fcvode_extras_src

  implicit none

  contains

subroutine fcvode_wrapper_with_source(dt, rho_in, T_in, ne_in, e_in, neq, cvmem, &
                              sunvec_y, yvec, rho_out, T_out, ne_out, e_out, rho_src, e_src)

    use amrex_fort_module, only : rt => amrex_real
    use vode_aux_module, only: rho_vode, T_vode, ne_vode, &
                               i_vode, j_vode, k_vode, NR_vode, rho_src_vode, e_src_vode,&
                               rho_init_vode

    use eos_params_module
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_OmB

    use cvode_interface
    use fnvector_serial
    use eos_module, only: vode_rtol, vode_atol_scaled
    use, intrinsic :: iso_c_binding

    implicit none

    real(rt), intent(in   ) :: dt
    real(rt), intent(in   ) :: rho_in, T_in, ne_in, e_in
    type(c_ptr), value :: cvmem
    type(c_ptr), value :: sunvec_y
    real(rt), intent(  out) ::  rho_out, T_out,ne_out,e_out
    real(rt),  intent(in   ) ::         rho_src, e_src

    real(rt) mean_rhob

    real(c_double), pointer, dimension(:) :: atol
    real(c_double) :: rtol
    real(c_double) :: time, tout
    integer(c_long), intent(in) :: neq
    real(c_double), pointer, intent(in) :: yvec(:)
    
    integer(c_int) :: ierr
    real(c_double) :: t_soln
    type(c_ptr) :: sunvec_atol
    logical, save :: firstCall = .true.

    allocate(atol(neq))

    T_vode   = T_in
    ne_vode  = ne_in
    rho_vode = rho_in
    rho_init_vode = rho_in
    NR_vode  = 0
    rho_src_vode = rho_src
    e_src_vode = e_src

    ! Initialize the integration time
    time = 0.d0
    
    ! We will integrate "e" and "rho" in time. 
    yvec(1) = e_in
    yvec(2) = rho_in

    ! Set the tolerances.  
    atol(1) = 1.d-4 * e_in
    atol(2) = 1.d-4 * rho_in
    rtol = 1.d-4

    sunvec_atol = N_VMake_Serial(neq, atol)
    ierr = FCVodeReInit(cvmem, time, sunvec_y)
    ierr = FCVodeSVtolerances(CVmem, rtol, sunvec_atol)
    
    ierr = FCVode(CVmem, dt, sunvec_y, time, CV_NORMAL)
    
    e_out  = yvec(1)
    rho_out = yvec(2)
!    rho_out = rho_init_vode + dt * rho_src_vode
    T_out  = T_vode
    ne_out = ne_vode

    if (ierr .ne. 0) then
       print *, 'istate = ', ierr, 'at (i,j,k) ',i_vode,j_vode,k_vode
       call bl_error("ERROR in fcvode_wrapper_with_src: integration failed")
    endif

    call N_VDestroy_Serial(sunvec_atol)

end subroutine fcvode_wrapper_with_source

subroutine fcvode_wrapper_with_source_single(dt, rho_in, T_in, ne_in, e_in, neq, cvmem, &
                              sunvec_y, yvec, rho_out, T_out, ne_out, e_out, rho_src, e_src)

    use amrex_fort_module, only : rt => amrex_real
    use vode_aux_module, only: rho_vode, T_vode, ne_vode, &
                               i_vode, j_vode, k_vode, NR_vode, rho_src_vode, e_src_vode,&
                               rho_init_vode

    use eos_params_module
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_OmB

    use cvode_interface
    use fnvector_serial
    use eos_module, only: vode_rtol, vode_atol_scaled
    use, intrinsic :: iso_c_binding

    implicit none

    real(rt), intent(in   ) :: dt
    real(rt), intent(in   ) :: rho_in, T_in, ne_in, e_in
    type(c_ptr), value :: cvmem
    type(c_ptr), value :: sunvec_y
    real(rt), intent(  out) ::  rho_out, T_out,ne_out,e_out
    real(rt),  intent(in   ) ::         rho_src, e_src

    real(rt) mean_rhob

    real(c_double), pointer, dimension(:) :: atol
    real(c_double) :: rtol
    real(c_double) :: time, tout
    integer(c_long), intent(in) :: neq
    real(c_double), pointer, intent(in) :: yvec(:)
    
    integer(c_int) :: ierr, maxord
    integer(c_long) :: mxsteps
    real(c_double) :: t_soln
    type(c_ptr) :: sunvec_atol
    logical, save :: firstCall = .true.

    allocate(atol(neq))

    T_vode   = T_in
    ne_vode  = ne_in
    rho_vode = rho_in
    rho_init_vode = rho_in
    NR_vode  = 0
    rho_src_vode = rho_src
    e_src_vode = e_src

    ! Initialize the integration time
    time = 0.d0
    
    ! We will integrate "e" and "rho" in time. 
    yvec(1) = e_in
    yvec(2) = rho_in

    ! Set the tolerances.  
    atol(1) = 1.d-4 * e_in
    atol(2) = 1.d-4 * rho_in
    rtol = 1.d-4

    maxord = 1
    mxsteps = 1

    sunvec_atol = N_VMake_Serial(neq, atol)
    ierr = FCVodeReInit(cvmem, time, sunvec_y)
    ierr = FCVodeSVtolerances(CVmem, rtol, sunvec_atol)

    ierr = FCVodeSetMaxOrd(CVmem, maxord) 
    ierr = FCVodeSetMaxNumSteps(CVmem, mxsteps) 
    ierr = FCVodeSetInitStep(CVmem, dt) 
    
    ierr = FCVode(CVmem, dt, sunvec_y, time, CV_NORMAL)
    
    e_out  = yvec(1)
    rho_out = yvec(2)
!    rho_out = rho_init_vode + dt * rho_src_vode
    T_out  = T_vode
    ne_out = ne_vode

    if (ierr .ne. 0) then
       print *, 'istate = ', ierr, 'at (i,j,k) ',i_vode,j_vode,k_vode
       call bl_error("ERROR in fcvode_wrapper_with_src: integration failed")
    endif

    call N_VDestroy_Serial(sunvec_atol)

end subroutine fcvode_wrapper_with_source_single

    integer(c_int) function RhsFn_src(tn, sunvec_y, sunvec_f, user_data) &
           result(ierr) bind(C,name='RhsFn_src')

      use, intrinsic :: iso_c_binding
      use fnvector_serial
      use cvode_interface
      implicit none

      real(c_double), value :: tn
      type(c_ptr), value    :: sunvec_y, sunvec_f, user_data

      ! pointers to data in SUNDAILS vectors
      real(c_double), dimension(:), pointer :: yvec, fvec

      integer(c_long) :: neq
      real(c_double) :: energy(2)

      neq = int(2, c_long)

      ! get data arrays from SUNDIALS vectors
      call N_VGetData_Serial(sunvec_y, neq, yvec)
      call N_VGetData_Serial(sunvec_f, neq, fvec)

      call f_rhs_split(2, tn, yvec, energy, 0.0, 0)

      fvec = energy

      ierr = 0
    end function RhsFn_src
end module fcvode_extras_src

