module v_cycle_mod

contains

    recursive subroutine v_cycle(ct,lvl,nlvl)

        use constants
        use multifab_module
        use traverse_mod
        use cycle_tower_mod

        implicit none

        ! Note: nlvl is different than ct%nlvl when in the FMG-cycle, because nlvl
        ! is the number of multigrid levels in whichever "mini" V-cycle is being called by
        ! the FMG-cycle, but ct%nlvl will always be the number of multigrid levels
        ! of the FMG-cycle.
        integer          , intent(in)    :: lvl ! Current multigrid level of this V-cycle
        integer          , intent(in)    :: nlvl ! Total multigrid levels of this V-cycle
        type(cycle_tower), intent(inout) :: ct

        integer :: i

        ! If not at the top, set solution guess to zero; use existing solution guess otherwise
        if (lvl /= nlvl) then
            call setval(ct%uu(lvl),0.d0,all=.true.)
        end if

        if (lvl == 1) then

            if (parallel_IOProcessor()) then
                call print_cycle_tower(ct,lvl,'Before relaxation: ')
                if (ct%verbose >= 2) then
                    write(*,'(a,i3,a)') '    V Cycle BOTTOM -- Applying ',ct%v_bot,' relaxations.'
                end if
            end if

            ! Just smooth a bunch of times at the bottom; can choose different solver here
            do i = 1,ct%v_bot
                call gsrb(ct%uu(lvl),ct%ff(lvl),ct%bb(lvl,:),ct%dx(lvl))
            end do

            if (parallel_IOProcessor()) then
                call print_cycle_tower(ct,lvl,' After relaxation: ')
            end if

        else

            if (parallel_IOProcessor()) then
                call print_cycle_tower(ct,lvl,'Before relaxation: ')
            end if

            ! Smooth on the way down
            do i = 1,ct%s1
                call gsrb(ct%uu(lvl),ct%ff(lvl),ct%bb(lvl,:),ct%dx(lvl))
            end do

            if (parallel_IOProcessor()) then
                call print_cycle_tower(ct,lvl,' After relaxation: ')
            end if

            ! Coarsen and continue downwards
            call residual(ct%res(lvl),ct%uu(lvl),ct%ff(lvl),ct%bb(lvl,:),ct%dx(lvl))
            call restriction(ct%ff(lvl-1),ct%res(lvl))
            call v_cycle(ct,lvl-1,nlvl)

            if (parallel_IOProcessor()) then
                call print_cycle_tower(ct,lvl,'Before relaxation: ')
            end if

            ! Smooth on the way up
            do i = 1,ct%s2
                call gsrb(ct%uu(lvl),ct%ff(lvl),ct%bb(lvl,:),ct%dx(lvl))
            end do

            if (parallel_IOProcessor()) then
                call print_cycle_tower(ct,lvl,' After relaxation: ')
            end if

        end if

        ! If not at the top of this V-cycle, interpolate up and apply correction (i.e. prolongation)
        if (lvl /= nlvl) then
            call prolongation(ct%uu(lvl+1),ct%uu(lvl),ct%interp_type)
        end if

        ! If at the top of this V-cycle, and if this V-cycle is NOT part of an FMG-cycle,
        ! then save solution and fill ghost cells.
        if (lvl == nlvl .and. ct%cycle_type /= FMG_CYCLE_TYPE) then
            call copy(ct%sol,ct%uu(lvl))
            call multifab_fill_boundary(ct%sol)
        end if

    end subroutine v_cycle

end module v_cycle_mod


