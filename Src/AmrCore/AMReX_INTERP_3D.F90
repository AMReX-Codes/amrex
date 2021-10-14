
module amrex_interp_module

  use amrex_fort_module
  use amrex_constants_module
  use amrex_bc_types_module
  use amrex_error_module

  implicit none

contains

    subroutine AMREX_CQINTERP (fine, fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3, &
                              fb_l1,fb_l2,fb_l3,fb_h1,fb_h2,fb_h3, &
                              nvar, lratiox, lratioy, lratioz, crse, &
                              clo, chi, cb_l1,cb_l2,cb_l3,cb_h1,cb_h2,cb_h3, &
                              fslo, fshi, cslope, clen, fslope, fdat, &
                              flen, voff, bc, limslope, &
                              fvcx, fvcy, fvcz, cvcx, cvcy, cvcz, &
                              actual_comp, actual_state) bind(c,name='amrex_cqinterp')

      implicit none

      integer fine_l1,fine_l2,fine_l3,fine_h1,fine_h2,fine_h3
      integer fb_l1,fb_l2,fb_l3,fb_h1,fb_h2,fb_h3
      integer cb_l1,cb_l2,cb_l3,cb_h1,cb_h2,cb_h3
      integer fslo(3), fshi(3)
      integer nvar, lratiox, lratioy, lratioz
      integer bc(3,2,nvar)
      integer clen, flen, clo, chi, limslope
      integer actual_comp,actual_state
      real(amrex_real) fine(fine_l1:fine_h1,fine_l2:fine_h2,fine_l3:fine_h3,nvar)
      real(amrex_real) crse(clo:chi, nvar)
      real(amrex_real) cslope(clo:chi, 3)
      real(amrex_real) fslope(flen, 3)
      real(amrex_real) fdat(flen)
      real(amrex_real) voff(flen)
      real(amrex_real) fvcx(fb_l1:fb_h1+1)
      real(amrex_real) fvcy(fb_l2:fb_h2+1)
      real(amrex_real) fvcz(fb_l3:fb_h3+1)
      real(amrex_real) cvcx(cb_l1:cb_h1+1)
      real(amrex_real) cvcy(cb_l2:cb_h2+1)
      real(amrex_real) cvcz(cb_l3:cb_h3+1)

      call bl_abort('QUADRATIC INTERP NOT IMPLEMENTED IN 3-D')

    end subroutine AMREX_CQINTERP

end module amrex_interp_module
