module amrex_rungekutta_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  use amrex_multifab_module, only : amrex_multifab

  implicit none
  private

  public :: amrex_rungekutta_rk2, amrex_rungekutta_rk3, amrex_rungekutta_rk4
  public :: amrex_rungekutta_rhs_proc, amrex_rungekutta_fill_proc
  public :: amrex_rungekutta_store3_proc, amrex_rungekutta_store4_proc
  public :: amrex_rungekutta_post_stage_proc

  ! The ctx argument is an opaque user pointer passed unchanged to every
  ! callback. AMReX does not dereference it. A Fortran caller can pass
  ! c_loc(my_context) and recover the application-specific data in each
  ! callback with c_f_pointer.
  interface
     ! Compute the RK stage right-hand side. The rhs and state arguments are
     ! C pointers to MultiFabs. The stage number starts at 1. The time
     ! argument is the stage time, and dtsub is the sub-step size used by the
     ! C++ RK helper for this stage. AMR applications that reflux should also
     ! update their flux registers from this callback, passing flux-register
     ! scales consistent with dtsub for the current RK stage. The rhs itself
     ! should be the unscaled time derivative.
     subroutine amrex_rungekutta_rhs_proc (ctx, stage, rhs, state, time, dtsub) bind(c)
       import :: c_int, c_ptr, amrex_real
       implicit none
       type(c_ptr), value :: ctx, rhs, state
       integer(c_int), value :: stage
       real(amrex_real), value :: time, dtsub
     end subroutine amrex_rungekutta_rhs_proc

     ! Prepare a stage state before rhs is evaluated. This is typically where
     ! ghost cells, physical boundary conditions, and coarse/fine data from a
     ! FillPatcher are filled.
     subroutine amrex_rungekutta_fill_proc (ctx, stage, state, time) bind(c)
       import :: c_int, c_ptr, amrex_real
       implicit none
       type(c_ptr), value :: ctx, state
       integer(c_int), value :: stage
       real(amrex_real), value :: time
     end subroutine amrex_rungekutta_fill_proc

     ! Optional RK3 storage hook called after all RK3 stages are complete.
     ! The rk array contains C pointers to the temporary stage right-hand-side
     ! MultiFabs. Deep-copy any data that must live after the callback returns.
     subroutine amrex_rungekutta_store3_proc (ctx, old_state, rk) bind(c)
       import :: c_ptr
       implicit none
       type(c_ptr), value :: ctx, old_state
       type(c_ptr), intent(in) :: rk(*)
     end subroutine amrex_rungekutta_store3_proc

     ! Optional RK4 storage hook called after all RK4 stages are complete.
     ! The rk array contains C pointers to the temporary stage right-hand-side
     ! MultiFabs. Deep-copy any data that must live after the callback returns.
     subroutine amrex_rungekutta_store4_proc (ctx, old_state, rk) bind(c)
       import :: c_ptr
       implicit none
       type(c_ptr), value :: ctx, old_state
       type(c_ptr), intent(in) :: rk(*)
     end subroutine amrex_rungekutta_store4_proc

     ! Optional hook called after each RK stage has updated state.
     subroutine amrex_rungekutta_post_stage_proc (ctx, stage, state) bind(c)
       import :: c_int, c_ptr
       implicit none
       type(c_ptr), value :: ctx, state
       integer(c_int), value :: stage
     end subroutine amrex_rungekutta_post_stage_proc
  end interface

  interface
     subroutine amrex_fi_rungekutta_rk2 (old_state, new_state, time, dt, ctx, &
          rhs, fill, post_stage) bind(c)
       import :: c_funptr, c_ptr, amrex_real
       implicit none
       type(c_ptr), value :: old_state, new_state, ctx
       type(c_funptr), value :: rhs, fill, post_stage
       real(amrex_real), value :: time, dt
     end subroutine amrex_fi_rungekutta_rk2

     subroutine amrex_fi_rungekutta_rk3 (old_state, new_state, time, dt, ctx, &
          rhs, fill, store, post_stage) bind(c)
       import :: c_funptr, c_ptr, amrex_real
       implicit none
       type(c_ptr), value :: old_state, new_state, ctx
       type(c_funptr), value :: rhs, fill, store, post_stage
       real(amrex_real), value :: time, dt
     end subroutine amrex_fi_rungekutta_rk3

     subroutine amrex_fi_rungekutta_rk4 (old_state, new_state, time, dt, ctx, &
          rhs, fill, store, post_stage) bind(c)
       import :: c_funptr, c_ptr, amrex_real
       implicit none
       type(c_ptr), value :: old_state, new_state, ctx
       type(c_funptr), value :: rhs, fill, store, post_stage
       real(amrex_real), value :: time, dt
     end subroutine amrex_fi_rungekutta_rk4
  end interface

contains

  ! Advance old_state from time to time+dt with RK2 and put the result in
  ! new_state. The rhs and fill callbacks are required. The post_stage callback
  ! can be omitted if no work is needed after each stage update.
  subroutine amrex_rungekutta_rk2 (old_state, new_state, time, dt, ctx, &
       rhs, fill, post_stage)
    type(amrex_multifab), intent(inout) :: old_state, new_state
    real(amrex_real), intent(in) :: time, dt
    type(c_ptr), intent(in) :: ctx
    procedure(amrex_rungekutta_rhs_proc) :: rhs
    procedure(amrex_rungekutta_fill_proc) :: fill
    procedure(amrex_rungekutta_post_stage_proc), optional :: post_stage

    type(c_funptr) :: post_stage_ptr

    post_stage_ptr = c_null_funptr
    if (present(post_stage)) post_stage_ptr = c_funloc(post_stage)

    call amrex_fi_rungekutta_rk2(old_state%p, new_state%p, time, dt, ctx, &
         c_funloc(rhs), c_funloc(fill), post_stage_ptr)
  end subroutine amrex_rungekutta_rk2

  ! Advance old_state from time to time+dt with SSPRK3 and put the result in
  ! new_state. The optional store callback is useful for retaining the RK stage
  ! right-hand sides needed by AMR coarse/fine FillPatcher operations.
  subroutine amrex_rungekutta_rk3 (old_state, new_state, time, dt, ctx, &
       rhs, fill, store, post_stage)
    type(amrex_multifab), intent(inout) :: old_state, new_state
    real(amrex_real), intent(in) :: time, dt
    type(c_ptr), intent(in) :: ctx
    procedure(amrex_rungekutta_rhs_proc) :: rhs
    procedure(amrex_rungekutta_fill_proc) :: fill
    procedure(amrex_rungekutta_store3_proc), optional :: store
    procedure(amrex_rungekutta_post_stage_proc), optional :: post_stage

    type(c_funptr) :: store_ptr, post_stage_ptr

    store_ptr = c_null_funptr
    if (present(store)) store_ptr = c_funloc(store)

    post_stage_ptr = c_null_funptr
    if (present(post_stage)) post_stage_ptr = c_funloc(post_stage)

    call amrex_fi_rungekutta_rk3(old_state%p, new_state%p, time, dt, ctx, &
         c_funloc(rhs), c_funloc(fill), store_ptr, post_stage_ptr)
  end subroutine amrex_rungekutta_rk3

  ! Advance old_state from time to time+dt with classical RK4 and put the
  ! result in new_state. The optional store callback is useful for retaining
  ! the RK stage right-hand sides needed by AMR coarse/fine FillPatcher
  ! operations.
  subroutine amrex_rungekutta_rk4 (old_state, new_state, time, dt, ctx, &
       rhs, fill, store, post_stage)
    type(amrex_multifab), intent(inout) :: old_state, new_state
    real(amrex_real), intent(in) :: time, dt
    type(c_ptr), intent(in) :: ctx
    procedure(amrex_rungekutta_rhs_proc) :: rhs
    procedure(amrex_rungekutta_fill_proc) :: fill
    procedure(amrex_rungekutta_store4_proc), optional :: store
    procedure(amrex_rungekutta_post_stage_proc), optional :: post_stage

    type(c_funptr) :: store_ptr, post_stage_ptr

    store_ptr = c_null_funptr
    if (present(store)) store_ptr = c_funloc(store)

    post_stage_ptr = c_null_funptr
    if (present(post_stage)) post_stage_ptr = c_funloc(post_stage)

    call amrex_fi_rungekutta_rk4(old_state%p, new_state%p, time, dt, ctx, &
         c_funloc(rhs), c_funloc(fill), store_ptr, post_stage_ptr)
  end subroutine amrex_rungekutta_rk4

end module amrex_rungekutta_module
