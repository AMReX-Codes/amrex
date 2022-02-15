module eb_diffusion_module
  use amrex_fort_module, only : rt=>amrex_real
  use cns_module, only : urho, umx, umy, umz, ueden, ueint, utemp, nvar, &
       qrho,qu,qv,qw,qp,qc,qeint,qtemp,qvar
  implicit none
  private

  public :: eb_diff_mol_3d

contains

  subroutine eb_diff_mol_3d (lo,hi, dx, &
       q,qd_lo,qd_hi, &
       lam, mu, xi, clo, chi, &
       flux1, fd1_lo, fd1_hi, &
       flux2, fd2_lo, fd2_hi, &
       flux3, fd3_lo, fd3_hi, &
       flag, fglo, fghi)
    use amrex_ebcellflag_module, only : get_neighbor_cells, get_neighbor_cells_int_single
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: lo(3), hi(3), fglo(3), fghi(3)
    integer, intent(in) :: clo(3), chi(3)
    integer, intent(in) :: fd1_lo(3), fd1_hi(3)
    integer, intent(in) :: fd2_lo(3), fd2_hi(3)
    integer, intent(in) :: fd3_lo(3), fd3_hi(3)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(in   ) ::     q( qd_lo(1): qd_hi(1), qd_lo(2): qd_hi(2), qd_lo(3): qd_hi(3),QVAR)
    real(rt), intent(in   ) :: lam  (   clo(1):   chi(1),   clo(2):   chi(2),   clo(3):   chi(3))
    real(rt), intent(in   ) :: mu   (   clo(1):   chi(1),   clo(2):   chi(2),   clo(3):   chi(3))
    real(rt), intent(in   ) :: xi   (   clo(1):   chi(1),   clo(2):   chi(2),   clo(3):   chi(3))
    real(rt), intent(  out) :: flux1(fd1_lo(1):fd1_hi(1),fd1_lo(2):fd1_hi(2),fd1_lo(3):fd1_hi(3),5)
    real(rt), intent(  out) :: flux2(fd2_lo(1):fd2_hi(1),fd2_lo(2):fd2_hi(2),fd2_lo(3):fd2_hi(3),5)
    real(rt), intent(  out) :: flux3(fd3_lo(1):fd3_hi(1),fd3_lo(2):fd3_hi(2),fd3_lo(3):fd3_hi(3),5)
    integer, intent(in) :: flag(fglo(1):fghi(1),fglo(2):fghi(2),fglo(3):fghi(3))

    integer :: i,j,k
    real(rt) :: dxinv(3)
    real(rt) :: tauxx, tauyy, tauzz, tauxy, tauxz, tauyz, muf, xif
    real(rt) :: dudx, dudy, dudz, dvdx, dvdy, dwdx, dwdz, divu, dTdx
    real(rt) :: dvdz, dwdy, dTdy, dTdz
    integer  :: nbrlo(-1:1,-1:1,-1:1), nbrhi(-1:1,-1:1,-1:1)
    integer  :: ihip, ihim, ilop, ilom, jhip, jhim, jlop, jlom, khip, khim, klop, klom
    real(rt) :: wlo, whi
    real(rt), parameter :: weights(0:2) = [0.d0, 1.d0, 0.5d0]
    real(rt), parameter :: twoThirds = 2.d0/3.d0
    integer, parameter :: nextra = 2

    dxinv = 1.d0/dx

    ! x-direction
    do       k = lo(3)-nextra-1, hi(3)+nextra+1
       do    j = lo(2)-nextra-1, hi(2)+nextra+1
          do i = lo(1)-nextra  , hi(1)+nextra+1
             dTdx = (q(i,j,k,qtemp)-q(i-1,j,k,qtemp))*dxinv(1)
             dudx = (q(i,j,k,qu)-q(i-1,j,k,qu))*dxinv(1)
             dvdx = (q(i,j,k,qv)-q(i-1,j,k,qv))*dxinv(1)
             dwdx = (q(i,j,k,qw)-q(i-1,j,k,qw))*dxinv(1)

             jhip = j + get_neighbor_cells_int_single(flag(i  ,j,k),0, 1,0)
             jhim = j - get_neighbor_cells_int_single(flag(i  ,j,k),0,-1,0)
             jlop = j + get_neighbor_cells_int_single(flag(i-1,j,k),0, 1,0)
             jlom = j - get_neighbor_cells_int_single(flag(i-1,j,k),0,-1,0)
             whi = weights(jhip-jhim)
             wlo = weights(jlop-jlom)
             dudy = (0.5d0*dxinv(2)) * &
                  ((q(i  ,jhip,k,qu)-q(i  ,jhim,k,qu))*whi &
                  +(q(i-1,jlop,k,qu)-q(i-1,jlom,k,qu))*wlo)
             dvdy = (0.5d0*dxinv(2)) * &
                  ((q(i  ,jhip,k,qv)-q(i  ,jhim,k,qv))*whi &
                  +(q(i-1,jlop,k,qv)-q(i-1,jlom,k,qv))*wlo)

             khip = k + get_neighbor_cells_int_single(flag(i  ,j,k),0,0, 1)
             khim = k - get_neighbor_cells_int_single(flag(i  ,j,k),0,0,-1)
             klop = k + get_neighbor_cells_int_single(flag(i-1,j,k),0,0, 1)
             klom = k - get_neighbor_cells_int_single(flag(i-1,j,k),0,0,-1)
             whi = weights(khip-khim)
             wlo = weights(klop-klom)
             dudz = (0.5d0*dxinv(3)) * &
                  ((q(i  ,j,khip,qu)-q(i  ,j,khim,qu))*whi &
                  +(q(i-1,j,klop,qu)-q(i-1,j,klom,qu))*wlo)
             dwdz = (0.5d0*dxinv(3)) * &
                  ((q(i  ,j,khip,qw)-q(i  ,j,khim,qw))*whi &
                  +(q(i-1,j,klop,qw)-q(i-1,j,klom,qw))*wlo)

             divu = dudx + dvdy + dwdz
             muf = 0.5d0*(mu(i,j,k)+mu(i-1,j,k))
             xif = 0.5d0*(xi(i,j,k)+xi(i-1,j,k))
             tauxx = muf*(2.d0*dudx-twoThirds*divu) + xif*divu
             tauxy = muf*(dudy+dvdx)
             tauxz = muf*(dudz+dwdx)
             flux1(i,j,k,urho)  = 0.d0
             flux1(i,j,k,umx)   = -tauxx
             flux1(i,j,k,umy)   = -tauxy
             flux1(i,j,k,umz)   = -tauxz
             flux1(i,j,k,ueden) = -0.5d0*((q(i,j,k,qu)+q(i-1,j,k,qu))*tauxx &
                  &                      +(q(i,j,k,qv)+q(i-1,j,k,qv))*tauxy &
                  &                      +(q(i,j,k,qw)+q(i-1,j,k,qw))*tauxz &
                  &                    +(lam(i,j,k) +lam(i-1,j,k))*dTdx)
          end do
       end do
    end do

    ! y-direction
    do       k = lo(3)-nextra-1, hi(3)+nextra+1
       do    j = lo(2)-nextra  , hi(2)+nextra+1
          do i = lo(1)-nextra-1, hi(1)+nextra+1
             dTdy = (q(i,j,k,qtemp)-q(i,j-1,k,qtemp))*dxinv(2)
             dudy = (q(i,j,k,qu)-q(i,j-1,k,qu))*dxinv(2)
             dvdy = (q(i,j,k,qv)-q(i,j-1,k,qv))*dxinv(2)
             dwdy = (q(i,j,k,qw)-q(i,j-1,k,qw))*dxinv(2)

             ihip = i + get_neighbor_cells_int_single(flag(i,j  ,k), 1,0,0)
             ihim = i - get_neighbor_cells_int_single(flag(i,j  ,k),-1,0,0)
             ilop = i + get_neighbor_cells_int_single(flag(i,j-1,k), 1,0,0)
             ilom = i - get_neighbor_cells_int_single(flag(i,j-1,k),-1,0,0)
             whi = weights(ihip-ihim)
             wlo = weights(ilop-ilom)
             dudx = (0.5d0*dxinv(1)) * &
                  ((q(ihip,j  ,k,qu)-q(ihim,j  ,k,qu))*whi &
                  +(q(ilop,j-1,k,qu)-q(ilom,j-1,k,qu))*wlo)
             dvdx = (0.5d0*dxinv(1)) * &
                  ((q(ihip,j  ,k,qv)-q(ihim,j  ,k,qv))*whi &
                  +(q(ilop,j-1,k,qv)-q(ilom,j-1,k,qv))*wlo)

             khip = k + get_neighbor_cells_int_single(flag(i,j  ,k),0,0, 1)
             khim = k - get_neighbor_cells_int_single(flag(i,j  ,k),0,0,-1)
             klop = k + get_neighbor_cells_int_single(flag(i,j-1,k),0,0, 1)
             klom = k - get_neighbor_cells_int_single(flag(i,j-1,k),0,0,-1)
             whi = weights(khip-khim)
             wlo = weights(klop-klom)
             dvdz = (0.5d0*dxinv(3)) * &
                  ((q(i,j  ,khip,qv)-q(i,j  ,khim,qv))*whi &
                  +(q(i,j-1,klop,qv)-q(i,j-1,klom,qv))*wlo)
             dwdz = (0.5d0*dxinv(3)) * &
                  ((q(i,j  ,khip,qw)-q(i,j  ,khim,qw))*whi &
                  +(q(i,j-1,klop,qw)-q(i,j-1,klom,qw))*wlo)

             divu = dudx + dvdy + dwdz
             muf = 0.5d0*(mu(i,j,k)+mu(i,j-1,k))
             xif = 0.5d0*(xi(i,j,k)+xi(i,j-1,k))
             tauyy = muf*(2.d0*dvdy-twoThirds*divu) + xif*divu
             tauxy = muf*(dudy+dvdx)
             tauyz = muf*(dwdy+dvdz)
             flux2(i,j,k,urho)  = 0.d0
             flux2(i,j,k,umx)   = -tauxy
             flux2(i,j,k,umy)   = -tauyy
             flux2(i,j,k,umz)   = -tauyz
             flux2(i,j,k,ueden) = -0.5d0*((q(i,j,k,qu)+q(i,j-1,k,qu))*tauxy &
                  &                      +(q(i,j,k,qv)+q(i,j-1,k,qv))*tauyy &
                  &                      +(q(i,j,k,qw)+q(i,j-1,k,qw))*tauyz &
                  &                    +(lam(i,j,k) +lam(i,j-1,k))*dTdy)
          end do
       end do
    end do

    ! z-direction
    do       k = lo(3)-nextra  , hi(3)+nextra+1
       do    j = lo(2)-nextra-1, hi(2)+nextra+1
          do i = lo(1)-nextra-1, hi(1)+nextra+1
             dTdz = (q(i,j,k,qtemp)-q(i,j,k-1,qtemp))*dxinv(3)
             dudz = (q(i,j,k,qu)-q(i,j,k-1,qu))*dxinv(3)
             dvdz = (q(i,j,k,qv)-q(i,j,k-1,qv))*dxinv(3)
             dwdz = (q(i,j,k,qw)-q(i,j,k-1,qw))*dxinv(3)

             ihip = i + get_neighbor_cells_int_single(flag(i,j,k  ), 1,0,0)
             ihim = i - get_neighbor_cells_int_single(flag(i,j,k  ),-1,0,0)
             ilop = i + get_neighbor_cells_int_single(flag(i,j,k-1), 1,0,0)
             ilom = i - get_neighbor_cells_int_single(flag(i,j,k-1),-1,0,0)
             whi = weights(ihip-ihim)
             wlo = weights(ilop-ilom)
             dudx = (0.5d0*dxinv(1)) * &
                  ((q(ihip,j,k  ,qu)-q(ihim,j,k  ,qu))*whi &
                  +(q(ilop,j,k-1,qu)-q(ilom,j,k-1,qu))*wlo)
             dwdx = (0.5d0*dxinv(1)) * &
                  ((q(ihip,j,k  ,qw)-q(ihim,j,k  ,qw))*whi &
                  +(q(ilop,j,k-1,qw)-q(ilom,j,k-1,qw))*wlo)

             jhip = j + get_neighbor_cells_int_single(flag(i,j,k  ),0 ,1,0)
             jhim = j - get_neighbor_cells_int_single(flag(i,j,k  ),0,-1,0)
             jlop = j + get_neighbor_cells_int_single(flag(i,j,k-1),0 ,1,0)
             jlom = j - get_neighbor_cells_int_single(flag(i,j,k-1),0,-1,0)
             whi = weights(jhip-jhim)
             wlo = weights(jlop-jlom)
             dvdy = (0.5d0*dxinv(2)) * &
                  ((q(i,jhip,k  ,qv)-q(i,jhim,k  ,qv))*whi &
                  +(q(i,jlop,k-1,qv)-q(i,jlom,k-1,qv))*wlo)
             dwdy = (0.5d0*dxinv(2)) * &
                  ((q(i,jhip,k  ,qw)-q(i,jhim,k  ,qw))*whi &
                  +(q(i,jlop,k-1,qw)-q(i,jlom,k-1,qw))*wlo)

             divu = dudx + dvdy + dwdz
             muf = 0.5d0*(mu(i,j,k)+mu(i,j,k-1))
             xif = 0.5d0*(xi(i,j,k)+xi(i,j,k-1))
             tauxz = muf*(dudz+dwdx)
             tauyz = muf*(dvdz+dwdy)
             tauzz = muf*(2.d0*dwdz-twoThirds*divu) + xif*divu
             flux3(i,j,k,urho)  = 0.d0
             flux3(i,j,k,umx)   = -tauxz
             flux3(i,j,k,umy)   = -tauyz
             flux3(i,j,k,umz)   = -tauzz
             flux3(i,j,k,ueden) = -0.5d0*((q(i,j,k,qu)+q(i,j,k-1,qu))*tauxz &
                  &                      +(q(i,j,k,qv)+q(i,j,k-1,qv))*tauyz &
                  &                      +(q(i,j,k,qw)+q(i,j,k-1,qw))*tauzz &
                  &                      +(lam(i,j,k) +lam(i,j,k-1))*dTdz)
          end do
       end do
    end do
  end subroutine eb_diff_mol_3d

end module eb_diffusion_module
