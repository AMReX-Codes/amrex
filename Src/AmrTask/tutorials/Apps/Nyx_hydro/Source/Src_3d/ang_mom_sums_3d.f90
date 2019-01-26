! ::
! :: ----------------------------------------------------------
! ::
        subroutine volweight_sum_angularmomentum(lo,hi,problo,dx,fab, &
                                              fab_l1,fab_l2,fab_l3, &
                                              fab_h1,fab_h2,fab_h3, &
                                              dir,amom)
        use amrex_fort_module, only : rt => amrex_real
        use probdata_module, only: center
        implicit none
        integer, intent(in) :: lo(3),hi(3),fab_l1,fab_h1,fab_l2, &
                               fab_h2,fab_l3,fab_h3,dir
        real(rt), intent(in) :: problo(3),dx(3)
        real(rt), intent(in) :: fab(fab_l1:fab_h1, &
                                            fab_l2:fab_h2, &
                                            fab_l3:fab_h3)
        integer :: i,j,k
        real(rt) :: pos(3)
        real(rt), intent(out) :: amom(3)
        do i= 1,3
          amom(i) = 0.0d0
        enddo
        if(dir.eq.1) then
          do k=lo(3),hi(3)
            pos(3) = problo(3) + (DBLE(k)+0.5d0)*dx(3) - center(3)
            do j=lo(2),hi(2)
              pos(2) = problo(2) + (DBLE(j)+0.5d0)*dx(2) - center(2)
              do i=lo(1),hi(1)
                amom(2) = amom(2) - fab(i,j,k)*pos(3)
                amom(3) = amom(3) + fab(i,j,k)*pos(2)
              enddo
            enddo
          enddo
        elseif(dir.eq.2) then
          do k=lo(3),hi(3)
            pos(3) = problo(3) + (DBLE(k)+0.5d0)*dx(3) - center(3)
            do j=lo(2),hi(2)
              do i=lo(1),hi(1)
                pos(1) = problo(1) + (DBLE(i)+0.5d0)*dx(1) - center(1)
                amom(1) = amom(1) + fab(i,j,k)*pos(3)
                amom(3) = amom(3) - fab(i,j,k)*pos(1)
              enddo
            enddo
          enddo
        elseif(dir.eq.3) then
          do k=lo(3),hi(3)
            do j=lo(2),hi(2)
              pos(2) = problo(2) + (DBLE(j)+0.5d0)*dx(2) - center(2)
              do i=lo(1),hi(1)
                pos(1) = problo(1) + (DBLE(i)+0.5d0)*dx(1) - center(1)
                amom(2) = amom(2) + fab(i,j,k)*pos(1)
                amom(1) = amom(1) - fab(i,j,k)*pos(2)
              enddo
            enddo
          enddo
        endif
        do i=1,3
          amom(i) = amom(i)*(dx(1)*dx(2)*dx(3))
        enddo
!        write(*,*) '### amom ',amom
        return
        end subroutine volweight_sum_angularmomentum
