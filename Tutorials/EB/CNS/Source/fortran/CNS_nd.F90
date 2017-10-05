module cns_nd_module

  use amrex_fort_module, only : rt=>amrex_real
  use cns_module, only : urho, umx, umy, umz, ueden, ueint, utemp, nvar, &
       qrho,qu,qv,qw,qp,qc,qeint,qtemp,qvar
  implicit none
  private

  public :: cns_compute_temperature, cns_estdt, ctoprim, cns_eb_fixup_geom, cns_nullfill

contains

  subroutine cns_compute_temperature (lo,hi,u,ulo,uhi) bind(c,name='cns_compute_temperature')
    use cns_physics_module, only : cv
    integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3)
    real(rt), intent(inout) :: u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),nvar)

    integer :: i,j,k
    real(rt) :: rhoinv

    do       k = lo(3),hi(3)
       do    j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhoinv = 1.d0/u(i,j,k,urho)
             u(i,j,k,ueint) = u(i,j,k,ueden) - 0.5d0*rhoinv* &
                  (u(i,j,k,umx)**2 + u(i,j,k,umy)**2 + u(i,j,k,umz)**2)
             u(i,j,k,utemp) = rhoinv*u(i,j,k,ueint)*(1.d0/cv)
          end do
       end do
    end do
  end subroutine cns_compute_temperature

  subroutine cns_estdt (lo,hi,u,ulo,uhi,dx,dt) bind(c,name='cns_estdt')
    use cns_physics_module, only : gamma
    use cns_module, only : smallp, smallr
    integer, intent(in) :: lo(3), hi(3), ulo(3), uhi(3)
    real(rt), intent(in) :: u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),nvar)
    real(rt), intent(in) :: dx(3)
    real(rt), intent(inout) :: dt

    integer :: i,j,k
    real(rt) :: rhoinv, vx, vy, vz, p, cs

    do       k = lo(3),hi(3)
       do    j = lo(2),hi(2)
          do i = lo(1),hi(1)
             rhoinv = 1.d0/max(smallr,u(i,j,k,urho))
             vx = u(i,j,k,umx)*rhoinv
             vy = u(i,j,k,umy)*rhoinv
             vz = u(i,j,k,umz)*rhoinv
             p = max((gamma-1.d0)*u(i,j,k,ueint),smallp)
             cs = sqrt(gamma*p*rhoinv)
             dt = min(dt,dx(1)/(abs(vx)+cs),dx(2)/(abs(vy)+cs),dx(3)/(abs(vz)+cs))
          end do
       end do
    end do
  end subroutine cns_estdt

  subroutine ctoprim (lo,hi,u,ulo,uhi,q,qlo,qhi)
    use cns_physics_module, only : gamma, cv
    use cns_module, only : smallp, smallr
    integer, intent(in) :: lo(3),hi(3),ulo(3),uhi(3),qlo(3),qhi(3)
    real(rt), intent(in   ) :: u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),nvar)
    real(rt), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),qvar)
    
    integer :: i,j,k
    real(rt) :: rhoinv, kineng
    
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             q(i,j,k,qrho) = max(smallr,u(i,j,k,urho))
             rhoinv = 1.d0/q(i,j,k,qrho)
             q(i,j,k,qu) = u(i,j,k,umx)*rhoinv
             q(i,j,k,qv) = u(i,j,k,umy)*rhoinv
             q(i,j,k,qw) = u(i,j,k,umz)*rhoinv
             kineng = 0.5d0*q(i,j,k,qrho)*(q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2)
             q(i,j,k,qeint) = (u(i,j,k,ueden)-kineng) * rhoinv
             if (q(i,j,k,qeint) .le. 0.d0) then
                q(i,j,k,qeint) = u(i,j,k,ueint) * rhoinv
             end if
             q(i,j,k,qp) = max(smallp,(gamma-1.d0)*q(i,j,k,qeint)*q(i,j,k,qrho))
             q(i,j,k,qc) = sqrt(gamma*q(i,j,k,qp)*rhoinv)
             q(i,j,k,qtemp) = q(i,j,k,qeint) * (1.d0/cv)
          end do
       end do
    end do
  end subroutine ctoprim

  subroutine cns_eb_fixup_geom(lo, hi, flag, flo, fhi, vfrac, vlo, vhi, domlo, domhi) &
       bind(c,name='cns_eb_fixup_geom')
    use amrex_ebcellflag_module, only : is_regular_cell, is_single_valued_cell, get_neighbor_cells, &
         set_regular_cell
    integer, dimension(3), intent(in) :: lo, hi, flo, fhi, vlo, vhi, domlo, domhi
    integer, intent(inout) :: flag(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3))
    real(rt), intent(in)  :: vfrac(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))

    integer :: i,j,k, nbr(-1:1,-1:1,-1:1)
    real(rt), parameter :: almostone = 1.d0 - 1.d-13

    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             if (is_single_valued_cell(flag(i,j,k)) .and. vfrac(i,j,k) .gt. almostone) then
                call get_neighbor_cells(flag(i,j,k),nbr)
                if ( (nbr(-1,0,0).eq.1 .or. i.eq.domlo(1)) .and. &
                     (nbr( 1,0,0).eq.1 .or. i.eq.domhi(1)) .and. &
                     (nbr(0,-1,0).eq.1 .or. j.eq.domlo(2)) .and. &
                     (nbr(0, 1,0).eq.1 .or. j.eq.domhi(2)) .and. &
                     (nbr(0,0,-1).eq.1 .or. k.eq.domlo(3)) .and. &
                     (nbr(0,0, 1).eq.1 .or. k.eq.domhi(3)) ) then
                   call set_regular_cell(flag(i,j,k))
                end if
             end if
          end do
       end do
    end do
  end subroutine cns_eb_fixup_geom


  subroutine cns_nullfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="cns_nullfill")
    use amrex_fort_module, only: dim=>amrex_spacedim
    use amrex_error_module, only : amrex_error
    implicit none
    include 'AMReX_bc_types.fi'
    integer          :: adv_lo(3),adv_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    call amrex_error("How did this happen?")
  end subroutine cns_nullfill

end module cns_nd_module
