
module mytest
  use amrex_fort_module, only : amrex_real
  use amrex_constants_module, only : half, zero, one, three, seven, fifteen, m_pi
  use amrex_ebcellflag_module, only : is_regular_cell, is_covered_cell, is_single_valued_cell
  implicit none

contains

#if (AMREX_SPACEDIM == 2)

  subroutine mytest_set_phi_reg(lo, hi, xlo, xhi, ylo, yhi, phie, elo, ehi, &
       rhs, rlo, rhi, bx, bxlo, bxhi, by, bylo, byhi, dx, prob_type) &
       bind(c,name='mytest_set_phi_reg')
    integer, dimension(2), intent(in) :: lo, hi, elo, ehi, rlo, rhi, xlo, xhi, ylo, yhi, &
         bxlo, bxhi, bylo, byhi
    integer, intent(in) :: prob_type
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(inout) ::  phie( elo(1): ehi(1), elo(2): ehi(2))
    real(amrex_real), intent(inout) ::  rhs ( rlo(1): rhi(1), rlo(2): rhi(2))
    real(amrex_real), intent(inout) ::  bx  (bxlo(1):bxhi(1),bxlo(2):bxhi(2))
    real(amrex_real), intent(inout) ::  by  (bylo(1):byhi(1),bylo(2):byhi(2))

    integer :: i,j
    real(amrex_real) :: x, y, r2, theta

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          x = (i+half)*dx(1) - half
          y = (j+half)*dx(2) - half
          theta = atan2(x,y) + half*m_pi
          r2 = x**2 + y**2
          phie(i,j) = r2**2 * cos(three*theta)
          if (prob_type .eq. 1) then
             rhs(i,j) = -seven * r2 * cos(three*theta)
          else
             rhs(i,j) = -(seven * r2 - fifteen * r2**2) * cos(three*theta)
          end if
       end do
    end do

    if (prob_type .eq. 2) then
       do    j = xlo(2), xhi(2)
          do i = xlo(1), xhi(1)
             x = (i     )*dx(1) - half
             y = (j+half)*dx(2) - half
             theta = atan2(x,y) + half*m_pi
             r2 = x**2 + y**2
             bx(i,j) = one - r2
          end do
       end do

       do    j = ylo(2), yhi(2)
          do i = ylo(1), yhi(1)
             x = (i+half)*dx(1) - half
             y = (j     )*dx(2) - half
             theta = atan2(x,y) + half*m_pi
             r2 = x**2 + y**2
             by(i,j) = one - r2
          end do
       end do
    end if
  end subroutine mytest_set_phi_reg

  subroutine mytest_set_phi_eb(lo, hi, xlo, xhi, ylo, yhi, phie, elo, ehi, &
       phib, blo, bhi, rhs, rlo, rhi, &
       bx, bxlo, bxhi, by, bylo, byhi, bb, bblo, bbhi, &
       flag, flo, fhi, cent, tlo, thi, bcent, clo, chi, dx, prob_type) &
       bind(c,name='mytest_set_phi_eb')
    integer, dimension(2), intent(in) :: lo, hi, xlo, xhi, ylo, yhi, elo, ehi, blo, bhi, &
         rlo, rhi, bxlo, bxhi, bylo, byhi, bblo, bbhi, flo, fhi, tlo, thi, clo, chi
    integer, intent(in) :: prob_type
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(inout) ::  phie(elo(1):ehi(1),elo(2):ehi(2))
    real(amrex_real), intent(inout) ::  phib(blo(1):bhi(1),blo(2):bhi(2))
    real(amrex_real), intent(inout) ::  rhs (rlo(1):rhi(1),rlo(2):rhi(2))
    real(amrex_real), intent(inout) ::  bx  (bxlo(1):bxhi(1),bxlo(2):bxhi(2))
    real(amrex_real), intent(inout) ::  by  (bylo(1):byhi(1),bylo(2):byhi(2))
    real(amrex_real), intent(inout) ::  bb  (bblo(1):bbhi(1),bblo(2):bbhi(2))
    real(amrex_real), intent(in   ) ::  cent(tlo(1):thi(1),tlo(2):thi(2),2)
    real(amrex_real), intent(in   ) :: bcent(clo(1):chi(1),clo(2):chi(2),2)
    integer         , intent(in   ) ::  flag(flo(1):fhi(1),flo(2):fhi(2))

    integer :: i,j
    real(amrex_real) :: x, y, r2, theta

    do    j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (is_covered_cell(flag(i,j))) then
             phie(i,j) = zero
             phib(i,j) = zero
          else
             x = (i+half)*dx(1) - half
             y = (j+half)*dx(2) - half
             theta = atan2(x,y) + half*m_pi
             r2 = x**2 + y**2
             phie(i,j) = r2**2 * cos(three*theta)

             x = (i+half+cent(i,j,1))*dx(1) - half
             y = (j+half+cent(i,j,2))*dx(2) - half
             theta = atan2(x,y) + half*m_pi
             r2 = x**2 + y**2
             if (prob_type .eq. 1) then
                rhs(i,j) = -seven * r2 * cos(three*theta)
             else
                rhs(i,j) = -(seven * r2 - fifteen * r2**2) * cos(three*theta)
             end if

             if (is_single_valued_cell(flag(i,j))) then
                x = (i+half+bcent(i,j,1))*dx(1) - half
                y = (j+half+bcent(i,j,2))*dx(2) - half
                theta = atan2(x,y) + half*m_pi
                r2 = x**2 + y**2
                phib(i,j) = r2**2 * cos(three*theta)
                if (prob_type .eq. 2) then
                   bb(i,j) = one - r2
                end if
             else
                phib(i,j) = zero
             end if
          end if
       end do
    end do

    if (prob_type .eq. 2) then
       do    j = xlo(2), xhi(2)
          do i = xlo(1), xhi(1)
             x = (i     )*dx(1) - half
             y = (j+half)*dx(2) - half
             theta = atan2(x,y) + half*m_pi
             r2 = x**2 + y**2
             bx(i,j) = one - r2
          end do
       end do

       do    j = ylo(2), yhi(2)
          do i = ylo(1), yhi(1)
             x = (i+half)*dx(1) - half
             y = (j     )*dx(2) - half
             theta = atan2(x,y) + half*m_pi
             r2 = x**2 + y**2
             by(i,j) = one - r2
          end do
       end do
    end if

  end subroutine mytest_set_phi_eb

  subroutine mytest_set_phi_boundary (glo, ghi, dlo, dhi, phi, phlo, phhi, dx) &
       bind(c,name='mytest_set_phi_boundary')
    integer, dimension(2), intent(in) :: dlo, dhi, glo, ghi, phlo, phhi
    real(amrex_real), intent(in) :: dx(2)
    real(amrex_real), intent(inout) :: phi(phlo(1):phhi(1),phlo(2):phhi(2))

    integer :: i,j
    real(amrex_real) :: x, y, r2, theta

    if (glo(1) < dlo(1)) then
       i = dlo(1)-1
       x = (i+1)*dx(1) - half
       do j = glo(2), ghi(2)
          y = (j+half)*dx(2) - half
          theta = atan2(x,y) + half*m_pi
          r2 = x**2 + y**2
          phi(i,j) = r2**2 * cos(three*theta)
       end do
    end if

    if (ghi(1) > dhi(1)) then
       i = dhi(1)+1
       x = (i)*dx(1) - half
       do j = glo(2), ghi(2)
          y = (j+half)*dx(2) - half
          theta = atan2(x,y) + half*m_pi
          r2 = x**2 + y**2
          phi(i,j) = r2**2 * cos(three*theta)
       end do
    end if

    if (glo(2) < dlo(2)) then
       j = dlo(2)-1
       y = (j+1)*dx(2) - half
       do i = glo(1), ghi(1)
          x = (i+half)*dx(1) - half
          theta = atan2(x,y) + half*m_pi
          r2 = x**2 + y**2
          phi(i,j) = r2**2 * cos(three*theta)
       end do
    end if

    if (ghi(2) > dhi(2)) then
       j = dhi(2)+1
       y = (j)*dx(2) - half
       do i = glo(1), ghi(1)
          x = (i+half)*dx(1) - half
          theta = atan2(x,y) + half*m_pi
          r2 = x**2 + y**2
          phi(i,j) = r2**2 * cos(three*theta)
       end do
    end if    
  end subroutine mytest_set_phi_boundary

#else

  subroutine mytest_set_phi_reg(lo, hi, xlo, xhi, ylo, yhi, zlo, zhi, phie, elo, ehi, &
       rhs, rlo, rhi, bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, dx, prob_type) &
       bind(c,name='mytest_set_phi_reg')
    integer, dimension(3), intent(in) :: lo, hi, elo, ehi, rlo, rhi, xlo, xhi, ylo, yhi, &
         zlo, zhi, bxlo, bxhi, bylo, byhi, bzlo, bzhi
    integer, intent(in) :: prob_type
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(inout) ::  phie( elo(1): ehi(1), elo(2): ehi(2), elo(3): ehi(3))
    real(amrex_real), intent(inout) ::  rhs ( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(inout) ::  bx  (bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
    real(amrex_real), intent(inout) ::  by  (bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
    real(amrex_real), intent(inout) ::  bz  (bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))

    integer :: i,j,k
    real(amrex_real) :: x, y, r2, theta

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             x = (i+half)*dx(1) - half
             y = (j+half)*dx(2) - half
             theta = atan2(x,y) + half*m_pi
             r2 = x**2 + y**2
             phie(i,j,k) = r2**2 * cos(three*theta)
             if (prob_type .eq. 1) then
                rhs(i,j,k) = -seven * r2 * cos(three*theta)
             else
                rhs(i,j,k) = -(seven * r2 - fifteen * r2**2) * cos(three*theta)
             end if
          end do
       end do
    end do

    if (prob_type .eq. 2) then
       do       k = xlo(3), xhi(3)
          do    j = xlo(2), xhi(2)
             do i = xlo(1), xhi(1)
                x = (i     )*dx(1) - half
                y = (j+half)*dx(2) - half
                theta = atan2(x,y) + half*m_pi
                r2 = x**2 + y**2
                bx(i,j,k) = one - r2
             end do
          end do
       end do

       do       k = ylo(3), yhi(3)
          do    j = ylo(2), yhi(2)
             do i = ylo(1), yhi(1)
                x = (i+half)*dx(1) - half
                y = (j     )*dx(2) - half
                theta = atan2(x,y) + half*m_pi
                r2 = x**2 + y**2
                by(i,j,k) = one - r2
             end do
          end do
       end do

       do       k = zlo(3), zhi(3)
          do    j = zlo(2), zhi(2)
             do i = zlo(1), zhi(1)
                x = (i+half)*dx(1) - half
                y = (j+half)*dx(2) - half
                theta = atan2(x,y) + half*m_pi
                r2 = x**2 + y**2
                bz(i,j,k) = one - r2
             end do
          end do
       end do
    end if
  end subroutine mytest_set_phi_reg

  subroutine mytest_set_phi_eb(lo, hi, xlo, xhi, ylo, yhi, zlo, zhi, phie, elo, ehi, &
       phib, blo, bhi, rhs, rlo, rhi, &
       bx, bxlo, bxhi, by, bylo, byhi, bz, bzlo, bzhi, bb, bblo, bbhi, &
       flag, flo, fhi, cent, tlo, thi, bcent, clo, chi, dx, prob_type) &
       bind(c,name='mytest_set_phi_eb')
    integer, dimension(3), intent(in) :: lo, hi, xlo, xhi, ylo, yhi, zlo, zhi, elo, ehi, blo, bhi, &
         rlo, rhi, bxlo, bxhi, bylo, byhi, bzlo, bzhi, bblo, bbhi, flo, fhi, tlo, thi, clo, chi
    integer, intent(in) :: prob_type
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(inout) ::  phie( elo(1): ehi(1), elo(2): ehi(2), elo(3): ehi(3))
    real(amrex_real), intent(inout) ::  phib( blo(1): bhi(1), blo(2): bhi(2), blo(3): bhi(3))
    real(amrex_real), intent(inout) ::  rhs ( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3))
    real(amrex_real), intent(inout) ::  bx  (bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))
    real(amrex_real), intent(inout) ::  by  (bylo(1):byhi(1),bylo(2):byhi(2),bylo(3):byhi(3))
    real(amrex_real), intent(inout) ::  bz  (bzlo(1):bzhi(1),bzlo(2):bzhi(2),bzlo(3):bzhi(3))
    real(amrex_real), intent(inout) ::  bb  (bblo(1):bbhi(1),bblo(2):bbhi(2),bblo(3):bbhi(3))
    real(amrex_real), intent(in   ) ::  cent( tlo(1): thi(1), tlo(2): thi(2), tlo(3): thi(3),3)
    real(amrex_real), intent(in   ) :: bcent( clo(1): chi(1), clo(2): chi(2), clo(3): chi(3),3)
    integer         , intent(in   ) ::  flag( flo(1): fhi(1), flo(2): fhi(2), flo(3): fhi(3))

    integer :: i,j,k
    real(amrex_real) :: x, y, r2, theta

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (is_covered_cell(flag(i,j,k))) then
                phie(i,j,k) = zero
                phib(i,j,k) = zero
             else
                x = (i+half)*dx(1) - half
                y = (j+half)*dx(2) - half
                theta = atan2(x,y) + half*m_pi
                r2 = x**2 + y**2
                phie(i,j,k) = r2**2 * cos(three*theta)
                
                x = (i+half+cent(i,j,k,1))*dx(1) - half
                y = (j+half+cent(i,j,k,2))*dx(2) - half
                theta = atan2(x,y) + half*m_pi
                r2 = x**2 + y**2
                if (prob_type .eq. 1) then
                   rhs(i,j,k) = -seven * r2 * cos(three*theta)
                else
                   rhs(i,j,k) = -(seven * r2 - fifteen * r2**2) * cos(three*theta)
                end if
                
                if (is_single_valued_cell(flag(i,j,k))) then
                   x = (i+half+bcent(i,j,k,1))*dx(1) - half
                   y = (j+half+bcent(i,j,k,2))*dx(2) - half
                   theta = atan2(x,y) + half*m_pi
                   r2 = x**2 + y**2
                   phib(i,j,k) = r2**2 * cos(three*theta)
                   if (prob_type .eq. 3) then
                      bb(i,j,k) = one - r2
                   end if
                else
                   phib(i,j,k) = zero
                end if
             end if
          end do
       end do
    end do

    if (prob_type .eq. 2) then
       do       k = xlo(3), xhi(3)
          do    j = xlo(2), xhi(2)
             do i = xlo(1), xhi(1)
                x = (i     )*dx(1) - half
                y = (j+half)*dx(2) - half
                theta = atan2(x,y) + half*m_pi
                r2 = x**2 + y**2
                bx(i,j,k) = one - r2
             end do
          end do
       end do

       do       k = ylo(3), yhi(3)
          do    j = ylo(2), yhi(2)
             do i = ylo(1), yhi(1)
                x = (i+half)*dx(1) - half
                y = (j     )*dx(2) - half
                theta = atan2(x,y) + half*m_pi
                r2 = x**2 + y**2
                by(i,j,k) = one - r2
             end do
          end do
       end do

       do       k = zlo(3), zhi(3)
          do    j = zlo(2), zhi(2)
             do i = zlo(1), zhi(1)
                x = (i+half)*dx(1) - half
                y = (j+half)*dx(2) - half
                theta = atan2(x,y) + half*m_pi
                r2 = x**2 + y**2
                bz(i,j,k) = one - r2
             end do
          end do
       end do
    end if
  end subroutine mytest_set_phi_eb

  subroutine mytest_set_phi_boundary (glo, ghi, dlo, dhi, phi, phlo, phhi, dx) &
       bind(c,name='mytest_set_phi_boundary')
    integer, dimension(3), intent(in) :: dlo, dhi, glo, ghi, phlo, phhi
    real(amrex_real), intent(in) :: dx(3)
    real(amrex_real), intent(inout) :: phi(phlo(1):phhi(1),phlo(2):phhi(2),phlo(3):phhi(3))

    integer :: i,j,k
    real(amrex_real) :: x, y, r2, theta

    if (glo(1) < dlo(1)) then
       i = dlo(1)-1
       x = (i+1)*dx(1) - half
       do    k = glo(3), ghi(3)
          do j = glo(2), ghi(2)
             y = (j+half)*dx(2) - half
             theta = atan2(x,y) + half*m_pi
             r2 = x**2 + y**2
             phi(i,j,k) = r2**2 * cos(three*theta)
          end do
       end do
    end if

    if (ghi(1) > dhi(1)) then
       i = dhi(1)+1
       x = (i)*dx(1) - half
       do    k = glo(3), ghi(3)
          do j = glo(2), ghi(2)
             y = (j+half)*dx(2) - half
             theta = atan2(x,y) + half*m_pi
             r2 = x**2 + y**2
             phi(i,j,k) = r2**2 * cos(three*theta)
          end do
       end do
    end if

    if (glo(2) < dlo(2)) then
       j = dlo(2)-1
       y = (j+1)*dx(2) - half
       do    k = glo(3), ghi(3)
          do i = glo(1), ghi(1)
             x = (i+half)*dx(1) - half
             theta = atan2(x,y) + half*m_pi
             r2 = x**2 + y**2
             phi(i,j,k) = r2**2 * cos(three*theta)
          end do
       end do
    end if

    if (ghi(2) > dhi(2)) then
       j = dhi(2)+1
       y = (j)*dx(2) - half
       do    k = glo(3), ghi(3)
          do i = glo(1), ghi(1)
             x = (i+half)*dx(1) - half
             theta = atan2(x,y) + half*m_pi
             r2 = x**2 + y**2
             phi(i,j,k) = r2**2 * cos(three*theta)
          end do
       end do
    end if

    if (glo(3) < dlo(3)) then
       k = dlo(3)-1
       do    j = glo(2), ghi(2)
          do i = glo(1), ghi(1)
             x = (i+half)*dx(1) - half
             y = (j+half)*dx(2) - half
             theta = atan2(x,y) + half*m_pi
             r2 = x**2 + y**2
             phi(i,j,k) = r2**2 * cos(three*theta)
          end do
       end do
    end if

    if (ghi(3) > dhi(3)) then
       k = dhi(3)+1
       do    j = glo(2), ghi(2)
          do i = glo(1), ghi(1)
             x = (i+half)*dx(1) - half
             y = (j+half)*dx(2) - half
             theta = atan2(x,y) + half*m_pi
             r2 = x**2 + y**2
             phi(i,j,k) = r2**2 * cos(three*theta)
          end do
       end do
    end if

  end subroutine mytest_set_phi_boundary

#endif

end module mytest
