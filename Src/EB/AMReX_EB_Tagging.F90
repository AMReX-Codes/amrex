module amrex_eb_tagging_module

  implicit none

  public

contains

  ! -----------------------------------------------------------
  !> This routine will tag high error cells based on the state
  !!
  !! \param tag        <=  integer tag array
  !! \param tag_lo,hi   => index extent of tag array
  !! \param state       => state array
  !! \param state_lo,hi => index extent of state array
  !! \param set         => integer value to tag cell for refinement
  !! \param clear       => integer value to untag cell
  !! \param lo,hi       => work region we are allowed to change
  !! \param dx          => cell size
  !! \param problo      => phys loc of lower left corner of prob domain
  !! \param time        => problem evolution time
  !!        level       => refinement level of this array
  ! -----------------------------------------------------------

  subroutine amrex_eb_levelset_error(tag,   tag_lo,   tag_hi,   &
                                     state, state_lo, state_hi, &
                                     set,   clear,              &
                                     lo,    hi,                 &
                                     dx,    problo,             &
                                     time,  phierr)             &
    bind(C, name="amrex_eb_levelset_error")

    use amrex_fort_module, only : amrex_real
    use amrex_eb_levelset_module, only: amrex_eb_interp_levelset

    implicit none

    integer, dimension(3), intent(in   ) :: lo, hi, state_lo, state_hi, tag_lo, tag_hi
    real(amrex_real),      intent(in   ) :: problo(3), dx(3), time, phierr
    integer,               intent(in   ) :: set, clear

    real(amrex_real),      intent(in   ) :: state(state_lo(1):state_hi(1), &
                                                  state_lo(2):state_hi(2), &
                                                  state_lo(3):state_hi(3))
    integer,               intent(  out) :: tag(tag_lo(1):tag_hi(1), &
                                                tag_lo(2):tag_hi(2), &
                                                tag_lo(3):tag_hi(3))

    integer          :: i, j, k
    real(amrex_real) :: pos(3), plo(3), ls_val

    plo = (/ 0., 0., 0. /)

    ! Tag on regions of high phi
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)

             pos = (/ i, j, k /)*dx + 0.5*dx
             call amrex_eb_interp_levelset ( &
                  pos, plo, 1,               &
                  state, state_lo, state_hi, &
                  dx, ls_val                 )

             if (abs(ls_val) .le. phierr) then
                tag(i, j, k) = set
             endif

          enddo
       enddo
    enddo

  end subroutine amrex_eb_levelset_error

end module amrex_eb_tagging_module
