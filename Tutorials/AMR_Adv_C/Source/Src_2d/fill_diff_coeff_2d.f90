      subroutine fill_diff_coeff(lo,hi, &
                                 coefx,coefx_l1,coefx_l2,coefx_h1,coefx_h2,&
                                 coefy,coefy_l1,coefy_l2,coefy_h1,coefy_h2,dx)

      use probdata_module, only: diff_coeff

      implicit none
      integer         , intent(in   ) :: lo(2),hi(2)
      integer         , intent(in   ) :: coefx_l1,coefx_l2,coefx_h1,coefx_h2
      integer         , intent(in   ) :: coefy_l1,coefy_l2,coefy_h1,coefy_h2
      double precision, intent(inout) :: coefx(coefx_l1:coefx_h1,coefx_l2:coefx_h2)
      double precision, intent(inout) :: coefy(coefy_l1:coefy_h1,coefy_l2:coefy_h2)
      double precision, intent(in   ) :: dx(2)

      ! Local variables
      integer          :: i,j

      do j = lo(2),hi(2)
        do i = lo(1),hi(1)+1
          coefx(i,j) = diff_coeff
        end do
      end do

      do j = lo(2),hi(2)+1
        do i = lo(1),hi(1)
          coefy(i,j) = diff_coeff
        end do
      end do

      end subroutine fill_diff_coeff
