module fmg_cycle_mod

contains

    recursive subroutine fmg_cycle(ct,lvl)

        use constants
        use traverse_mod
        use cycle_tower_mod
        use v_cycle_mod

        implicit none

        integer, intent(in) :: lvl
        type(cycle_tower), intent(inout) :: ct

        integer :: i 

        ! Always start with initial guess of zero for the solution
        call setval(ct%uu(lvl),0.d0,all=.true.)

        call print_cycle_tower(ct,lvl,'Down FMG-cycle: ')

        if (lvl == 1) then

            if (parallel_IOProcessor()) then
                call print_cycle_tower(ct,lvl,'Before relaxation: ')
            end if

            if (ct%verbose >= 2) then
                write(*,'(a,i3,a)') '    FMG-Cycle BOTTOM -- Applying ',ct%fmg_bot,' relaxations.'
            end if

            ! Just smooth a bunch of times at the bottom; can choose a different solver
            do i = 1,ct%fmg_bot
                call gsrb(ct%uu(lvl),ct%ff(lvl),ct%bb(lvl,:),ct%dx(lvl))
            end do

            if (parallel_IOProcessor()) then
                call print_cycle_tower(ct,lvl,' After relaxation: ')
            end if

        else

            ! Coarsen and continue downwards
            call restriction(ct%ff(lvl-1),ct%res(lvl))
            call fmg_cycle(ct,lvl-1)

            ! Interpolate solution
            call prolongation(ct%uu(lvl),ct%uu(lvl-1),ct%interp_type)

            ! Apply entire V-cycle on the way up
            if (ct%verbose == 1) then
                ct%verbose = ct%verbose - 1
                call v_cycle(ct,lvl,lvl)
                ct%verbose = ct%verbose + 1
            else
                call v_cycle(ct,lvl,lvl)
            end if

        end if

        ! If at the top, correct solution by adding the error that we've just solved for
        if (lvl == ct%nlvl) then
            call plus_plus(ct%sol,ct%uu(lvl))
            call multifab_fill_boundary(ct%sol)
        end if

    end subroutine fmg_cycle

end module fmg_cycle_mod


