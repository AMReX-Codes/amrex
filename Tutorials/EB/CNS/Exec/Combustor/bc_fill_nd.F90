module bc_fill_module

  implicit none

  public

contains

  ! All subroutines in this file must be threadsafe because they are called
  ! inside OpenMP parallel regions.
  
  subroutine cns_hypfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="cns_hypfill")

    use cns_module, only: NVAR
    use amrex_fort_module, only: dim=>amrex_spacedim
    use amrex_filcc_module
    use amrex_bc_types_module, only : amrex_bc_ext_dir
    use probdata_module, only : inflow_state

    implicit none

    integer          :: adv_lo(3),adv_hi(3)
    integer          :: bc(dim,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)

    integer          :: i,j,k,n

    do n = 1,NVAR
       call amrex_filcc(adv(:,:,:,n),adv_lo(1),adv_lo(2),adv_lo(3),adv_hi(1),adv_hi(2),adv_hi(3),domlo,domhi,delta,xlo,bc(:,:,n))
    enddo

    if ( bc(3,1,1).eq.amrex_bc_ext_dir .and. adv_lo(3).lt.domlo(3)) then
       do       k = adv_lo(3),  domlo(3)-1
          do    j = adv_lo(2), adv_hi(2)
             do i = adv_lo(1), adv_hi(1)
                adv(i,j,k,:) = inflow_state(:)
             end do
          end do
       end do
    end if

  end subroutine cns_hypfill



  subroutine cns_denfill(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,time,bc) &
       bind(C, name="cns_denfill")

    use amrex_fort_module, only: dim=>amrex_spacedim
    use amrex_filcc_module
    use amrex_bc_types_module, only : amrex_bc_ext_dir
    use probdata_module, only : inflow_state

    implicit none

    include 'AMReX_bc_types.fi'

    integer          :: adv_lo(3),adv_hi(3)
    integer          :: bc(dim,2)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

    integer :: i,j,k

    call amrex_filcc(adv,adv_lo(1),adv_lo(2),adv_lo(3),adv_hi(1),adv_hi(2),adv_hi(3),domlo,domhi,delta,xlo,bc)

    if ( bc(3,1).eq.amrex_bc_ext_dir .and. adv_lo(3).lt.domlo(3)) then
       do       k = adv_lo(3),  domlo(3)-1
          do    j = adv_lo(2), adv_hi(2)
             do i = adv_lo(1), adv_hi(1)
                adv(i,j,k) = inflow_state(1)
             end do
          end do
       end do
    end if

  end subroutine cns_denfill

end module bc_fill_module
