
module analriem_module

  contains

    subroutine analriem(ilo,ihi,gamma,pl,rl,ul,pr,rr,ur,smallp,pstar,ustar)

      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer, intent(in) :: ilo, ihi

      real(rt), intent(in ), dimension(ilo:ihi) :: pl,rl,ul,pr,rr,ur
      real(rt), intent(in ) :: gamma, smallp
      real(rt), intent(out), dimension(ilo:ihi) :: pstar,ustar

      ! Local variables
      real(rt), dimension(ilo:ihi) :: pstnm1,ustarp,ustarm
      real(rt), dimension(ilo:ihi) :: wl,wr,wlsq,wrsq
      real(rt), dimension(ilo:ihi) :: cleft,cright
      real(rt), dimension(ilo:ihi) :: dpditer,zp,zm,denom
      real(rt), dimension(ilo:ihi) :: ustnm1,ustnp1
      integer          :: iter

      real(rt), parameter :: weakwv = 1.d-3
      real(rt), parameter :: small  = 1.d-6

      integer :: i

      wl = sqrt(gamma*pl*rl)
      wr = sqrt(gamma*pr*rr)

      cleft  = wl/rl
      cright = wr/rr

      pstar=(wl*pr+wr*pl-wr*wl*(ur-ul))/(wl+wr)

      pstar=max(pstar,smallp)
      pstnm1 = pstar

      wlsq = (.5d0*(gamma-1.d0)*(pstar+pl)+pstar) * rl
      wrsq = (.5d0*(gamma-1.d0)*(pstar+pr)+pstar) * rr

      wl = sqrt(wlsq)
      wr = sqrt(wrsq)

      ustarp = ul - (pstar-pl)/wl
      ustarm = ur + (pstar-pr)/wr

      pstar = (wl*pr+wr*pl-wr*wl*(ur-ul))/(wl+wr)

      pstar=max(pstar,smallp)

      do iter = 1,3

        wlsq = (.5d0*(gamma-1.d0)*(pstar+pl)+pstar) * rl
        wrsq = (.5d0*(gamma-1.d0)*(pstar+pr)+pstar) * rr

        wl = 1.d0/sqrt(wlsq)
        wr = 1.d0/sqrt(wrsq)

        ustnm1 = ustarm
        ustnp1 = ustarp

        ustarm = ur - (pr-pstar)*wr
        ustarp = ul + (pl-pstar)*wl

        dpditer = abs(pstnm1-pstar)
        zp      = abs(ustarp-ustnp1)

        do i = ilo, ihi
          if (zp(i)-weakwv*cleft(i) .lt. 0.0d0 ) then
            zp(i) = dpditer(i)*wl(i)
          endif
        end do

        zm=abs(ustarm-ustnm1)

        do i = ilo, ihi
          if (zm(i)-weakwv*cright(i) .lt. 0.0d0 ) then
            zm(i) = dpditer(i)*wr(i)
          endif
        end do

        denom  = dpditer/max(zp+zm,small*(cleft+cright))
        pstnm1 = pstar

        pstar = pstar - denom*(ustarm-ustarp)

        pstar = max(pstar,smallp)

        ustar = 0.5d0*(ustarm+ustarp)

     end do

    end subroutine analriem

  end module analriem_module
