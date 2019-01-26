
      subroutine fort_get_grav_const(Gconst_out) bind(C,name="fort_get_grav_const")

         use amrex_fort_module, only : rt => amrex_real
         use fundamental_constants_module, only: Gconst

         real(rt) :: Gconst_out

         Gconst_out = Gconst

      end subroutine fort_get_grav_const

