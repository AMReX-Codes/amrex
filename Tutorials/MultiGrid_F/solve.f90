module solve_mod

contains

    subroutine solve(ct)

        use constants
        use cycle_tower_mod
        use traverse_mod
        use v_cycle_mod
        use fmg_cycle_mod

        implicit none

        type(cycle_tower), intent(inout) :: ct

        integer :: iter, nlvl
        double precision :: res0_norm, res_norm

        nlvl = ct%nlvl

        if (ct%cycle_type == V_CYCLE_TYPE) then

            ! Calculate initial residual
            call residual(ct%res(nlvl),ct%sol,ct%rhs,ct%bb(nlvl,:),ct%dx(nlvl))
            res_norm = norm_inf(ct%res(nlvl))

            if (parallel_IOProcessor()) then
                if (ct%verbose >= 1) then
                    print *, ''
                    write(*,'(a,i4,a,es28.18)') &
                    '    Starting V-Cycle Solve at n =', ct%n(nlvl),'.  norm(residual) =', res_norm
                end if
            end if

        else if (ct%cycle_type == FMG_CYCLE_TYPE) then

            ! Calculate initial residual
            call residual(ct%ff(nlvl),ct%sol,ct%rhs,ct%bb(nlvl,:),ct%dx(nlvl))
            res_norm = norm_inf(ct%ff(nlvl))

            if (parallel_IOProcessor()) then
                if (ct%verbose >= 1) then
                    write(*,*) ''
                    write(*,'(a,i4,a,es28.18)')&
                    '    Starting FMG-Cycle Solve at n=',ct%n(nlvl),'.   norm(residual) = ', res_norm
                end if
            end if

        end if

        res0_norm = res_norm

        do iter = 1,ct%max_iter

            if (ct%cycle_type == V_CYCLE_TYPE) then

                ! Approximate solution by executing one V-cycle
                call v_cycle(ct,nlvl,nlvl)

                ! Calculate residual after V-cycle
                call residual(ct%res(nlvl),ct%sol,ct%rhs,ct%bb(nlvl,:),ct%dx(nlvl))
                res_norm = norm_inf(ct%res(nlvl))

                if (parallel_IOProcessor()) then
                    if (ct%verbose >= 1) then
                        write(*,'(a,i3,a,i4,a,es28.18)')&
                        '   ', iter,' V-Cycle(s) completed at n =',ct%n(nlvl),'.  norm(residual) =', res_norm
                    end if
                end if

            else if (ct%cycle_type == FMG_CYCLE_TYPE) then

                ! Approximate error to solution by executing one FMG-cycle
                call fmg_cycle(ct,nlvl)

                ! Calculate residual of corrected solution and set it to RHS
                call residual(ct%ff(nlvl),ct%sol,ct%rhs,ct%bb(nlvl,:),ct%dx(nlvl))
                res_norm = norm_inf(ct%ff(nlvl))

                if (parallel_IOProcessor()) then
                    if (ct%verbose > 0) then
                        write(*,'(a,i2,a,i4,a,es28.18)')&
                        '   ',iter,' FMG-Cycle(s) completed at n=',ct%n(nlvl), '.  norm(residual) = ', res_norm
                    end if
                end if

            end if

            ! Check for convergence
            if (res_norm/res0_norm < ct%eps) then
                exit
            end if

        end do

        iter = iter - 1

        if (parallel_IOProcessor()) then
            print *, '' 
            print *, '******************************'
            if (ct%cycle_type == V_CYCLE_TYPE) then
                write(*,'(a,$)') ' V-Cycle'
            else if (ct%cycle_type == FMG_CYCLE_TYPE) then
                write(*,'(a,$)') ' FMG-Cycle'
            end if

            if (ct%dim == 2) then
                write(*,'(a,$)') ' (2D)'
            else if (ct%dim == 3) then
                write(*,'(a,$)') ' (3D)'
            end if

            if (ct%interp_type == LIN_INTERP) then
                write(*,'(a,$)') ', piecewise linear interpolation'
            else if (ct%interp_type == CONST_INTERP) then
                write(*,'(a,$)') ', piecewise constant interpolation'
            end if
            write(*,'(a,i4)') ', n = ', ct%n(ct%nlvl)

            write(*,'(a,i3)') ' Number of Cycles = ', iter
            write(*,'(a,es28.18)') ' Final residual = ', res_norm
            write(*,'(a,es28.18)') ' Final res/res0 =', res_norm/res0_norm
            write(*,'(a)') '******************************'
            write(*,'(a)') ''
        end if

    end subroutine solve


end module solve_mod

