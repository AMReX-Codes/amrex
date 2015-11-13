      subroutine fill_diff_coeff(lo,hi, &
                                 coefx,cx_l1,cx_l2,cx_l3,cx_h1,cx_h2,cx_h3,&
                                 coefy,cy_l1,cy_l2,cy_l3,cy_h1,cy_h2,cy_h3,&
                                 coefz,cz_l1,cz_l2,cz_l3,cz_h1,cz_h2,cz_h3,dx)

      use probdata_module, only: diff_coeff

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: cx_l1,cx_l2,cx_l3,cx_h1,cx_h2,cx_h3
      integer         , intent(in   ) :: cy_l1,cy_l2,cy_l3,cy_h1,cy_h2,cy_h3
      integer         , intent(in   ) :: cz_l1,cz_l2,cz_l3,cz_h1,cz_h2,cz_h3
      double precision, intent(inout) :: coefx(cx_l1:cx_h1,cx_l2:cx_h2,cx_l3:cx_h3)
      double precision, intent(inout) :: coefy(cy_l1:cy_h1,cy_l2:cy_h2,cy_l3:cy_h3)
      double precision, intent(inout) :: coefz(cz_l1:cz_h1,cz_l2:cz_h2,cz_l3:cz_h3)
      double precision, intent(in   ) :: dx(3)

      coefx = diff_coeff
      coefy = diff_coeff
      coefz = diff_coeff

      end subroutine fill_diff_coeff
