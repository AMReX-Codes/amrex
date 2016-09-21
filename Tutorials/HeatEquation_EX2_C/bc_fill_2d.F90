module bc_fill_module

  implicit none

  public

contains

  subroutine phifill(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
                     domlo,domhi,delta,xlo,time,bc) &
                     bind(C, name="phifill")

    implicit none

    include 'bc_types.fi'

    integer          :: phi_l1,phi_l2,phi_h1,phi_h2
    integer          :: bc(2,2,*)
    integer          :: domlo(2), domhi(2)
    double precision :: delta(2), xlo(2), time
    double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    call filcc(phi,phi_l1,phi_l2,phi_h1,phi_h2,domlo,domhi,delta,xlo,bc)

  end subroutine phifill
  
end module bc_fill_module
