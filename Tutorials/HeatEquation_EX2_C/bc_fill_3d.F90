module bc_fill_module

  implicit none

  public

contains

  subroutine phifill(phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2, &
                     phi_h3,domlo,domhi,delta,xlo,time,bc) &
                     bind(C, name="phifill")

    implicit none

    include 'bc_types.fi'

    integer          :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
    integer          :: bc(3,2,*)
    integer          :: domlo(3), domhi(3)
    double precision :: delta(3), xlo(3), time
    double precision :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)

    call filcc(phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3, &
               domlo,domhi,delta,xlo,bc)

  end subroutine phifill
  
end module bc_fill_module
