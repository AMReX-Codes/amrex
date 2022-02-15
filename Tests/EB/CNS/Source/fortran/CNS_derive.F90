
module cns_derive

  use amrex_fort_module, only : rt=>amrex_real
  implicit none
  private

  public :: cns_derpres, cns_dervel

contains

    ! This routine will derive pressure from internal energy density
    subroutine cns_derpres(p,p_lo,p_hi,ncomp_p, &
                           u,u_lo,u_hi,ncomp_u,lo,hi,domlo, &
                           domhi,dx,xlo,time,dt,bc,level,grid_no) &
                           bind(C, name="cns_derpres")

    use cns_module, only: URHO, UEINT
    use cns_physics_module, only : gamma

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: p_lo(3), p_hi(3), ncomp_p
    integer, intent(in) :: u_lo(3), u_hi(3), ncomp_u
    integer, intent(in) :: domlo(3), domhi(3)
    real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))
    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3))
    real(rt), intent(in) :: dx(3), xlo(3), time, dt
    integer, intent(in) :: bc(3,2,ncomp_u), level, grid_no

    integer  :: i, j, k
    real(rt) :: gm1

    gm1 = gamma-1.d0

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             p(i,j,k) = gm1*u(i,j,k)
          enddo
       enddo
    enddo

  end subroutine cns_derpres


  ! This routine will derive the velocity from the momentum.
  subroutine cns_dervel(vel,v_lo,v_hi,nv, &
                       dat,d_lo,d_hi,nc,lo,hi,domlo, &
                       domhi,delta,xlo,time,dt,bc,level,grid_no) &
                       bind(C, name="cns_dervel")
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: v_lo(3), v_hi(3), nv
    integer, intent(in) :: d_lo(3), d_hi(3), nc
    integer, intent(in) :: domlo(3), domhi(3)
    integer, intent(in) :: bc(3,2,nc)
    real(rt), intent(in) :: delta(3), xlo(3), time, dt
    real(rt), intent(inout) :: vel(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3))
    real(rt), intent(in   ) :: dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer, intent(in) :: level, grid_no

    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             vel(i,j,k) = dat(i,j,k,2) / dat(i,j,k,1)
          end do
       end do
    end do

  end subroutine cns_dervel

end module cns_derive
