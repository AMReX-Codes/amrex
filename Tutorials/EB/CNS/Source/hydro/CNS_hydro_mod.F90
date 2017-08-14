
module cns_hydro_module

  use amrex_fort_module, only : rt=>amrex_real
  use cns_module, only : urho, umx, umy, umz, ueden, ueint, utemp, nvar, &
       qrho,qu,qv,qw,qp,qc,qeint,qtemp,qvar
  use mempool_module, only : amrex_allocate, amrex_deallocate
  implicit none
  private

  public:: cns_compute_hydro_flux

contains

  subroutine cns_compute_hydro_flux (lo,hi,u,ulo,uhi,fx,fxlo,fxhi,fy,fylo,fyhi,fz,fzlo,fzhi,dx) &
       bind(c,name='cns_compute_hydro_flux')
    use advection_module, only : hyp_mol_gam_3d
    integer, intent(in) :: lo(3),hi(3),ulo(3),uhi(3),fxlo(3),fxhi(3),fylo(3),fyhi(3),fzlo(3),fzhi(3)
    real(rt), intent(in   ) :: u ( ulo(1): uhi(1), ulo(2): uhi(2), ulo(3): uhi(3),nvar)
    real(rt), intent(inout) :: fx(fxlo(1):fxhi(1),fxlo(2):fxhi(2),fxlo(3):fxhi(3),nvar)
    real(rt), intent(inout) :: fy(fylo(1):fyhi(1),fylo(2):fyhi(2),fylo(3):fyhi(3),nvar)
    real(rt), intent(inout) :: fz(fzlo(1):fzhi(1),fzlo(2):fzhi(2),fzlo(3):fzhi(3),nvar)
    real(rt), intent(in) :: dx(3)

    integer :: qlo(3), qhi(3)
    real(rt), pointer, contiguous :: q(:,:,:,:)

    qlo = lo - 4
    qhi = hi + 4
    call amrex_allocate(q, qlo(1),qhi(1), qlo(2),qhi(2), qlo(3),qhi(3), 1,qvar)
    
    call ctoprim(qlo, qhi, u, ulo, uhi, q, qlo, qhi)
    
    call hyp_mol_gam_3d(q, qlo, qhi, lo, hi, dx, fx, fxlo, fxhi, fy, fylo, fyhi, fz, fzlo, fzhi)

    call amrex_deallocate(q)
  end subroutine cns_compute_hydro_flux


  subroutine ctoprim (lo,hi,u,ulo,uhi,q,qlo,qhi)
    use cns_physics_module, only : gamma, cv
    integer, intent(in) :: lo(3),hi(3),ulo(3),uhi(3),qlo(3),qhi(3)
    real(rt), intent(in   ) :: u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),nvar)
    real(rt), intent(inout) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),qvar)
    
    integer :: i,j,k
    real(rt) :: rhoinv, kineng
    
    do       k = lo(3), hi(3)
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhoinv = 1.d0/u(i,j,k,urho)
             q(i,j,k,qrho) = u(i,j,k,urho)
             q(i,j,k,qu) = u(i,j,k,umx)*rhoinv
             q(i,j,k,qv) = u(i,j,k,umy)*rhoinv
             q(i,j,k,qw) = u(i,j,k,umz)*rhoinv
             kineng = 0.5d0*q(i,j,k,qrho)*(q(i,j,k,QU)**2 + q(i,j,k,QV)**2 + q(i,j,k,QW)**2)
             q(i,j,k,qeint) = (u(i,j,k,ueden)-kineng) * rhoinv
             if (q(i,j,k,qeint) .le. 0.d0) then
                q(i,j,k,qeint) = u(i,j,k,ueint) * rhoinv
             end if
             q(i,j,k,qp) = (gamma-1.d0)*q(i,j,k,qeint)*q(i,j,k,qrho)
             q(i,j,k,qc) = sqrt(gamma*q(i,j,k,qp)*rhoinv)
             q(i,j,k,qtemp) = q(i,j,k,qeint) * (1.d0/cv)
          end do
       end do
    end do
  end subroutine ctoprim

end module cns_hydro_module
