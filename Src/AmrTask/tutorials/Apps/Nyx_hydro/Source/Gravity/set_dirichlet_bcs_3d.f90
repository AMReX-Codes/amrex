
     subroutine fort_set_homog_bcs(lo,hi,domlo,domhi, &
                                   phi,phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3,dx);
 
     use amrex_fort_module, only : rt => amrex_real
     implicit none

     integer         ,intent(in   ) :: lo(3),hi(3),domlo(3),domhi(3)
     integer         ,intent(in   ) :: phi_l1,phi_l2,phi_l3,phi_h1,phi_h2,phi_h3
     real(rt),intent(  out) :: phi(phi_l1:phi_h1,phi_l2:phi_h2,phi_l3:phi_h3)
     real(rt),intent(in   ) :: dx(3)

     phi = 0.d0
 
     end subroutine fort_set_homog_bcs

