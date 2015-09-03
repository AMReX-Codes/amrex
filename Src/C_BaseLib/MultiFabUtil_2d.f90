subroutine bl_avg_fc_to_cc (lo, hi, &
     cc, ccl1, ccl2, cch1, cch2, &
     fx, fxl1, fxl2, fxh1, fxh2, &
     fy, fyl1, fyl2, fyh1, fyh2, &
     dx, problo, coord_type)

  implicit none
  integer          :: lo(2),hi(2), coord_type
  integer          :: ccl1, ccl2, cch1, cch2
  integer          :: fxl1, fxl2, fxh1, fxh2
  integer          :: fyl1, fyl2, fyh1, fyh2
  double precision :: cc(ccl1:cch1, ccl2:cch2, 2)
  double precision :: fx(fxl1:fxh1, fxl2:fxh2)
  double precision :: fy(fyl1:fyh1, fyl2:fyh2)
  double precision :: dx(2), problo(2)

  ! Local variables
  integer          :: i,j
  
  do j=lo(2),hi(2)
     do i=lo(1),hi(1)
        cc(i,j,1) = 0.5d0 * ( fx(i,j) + fx(i+1,j) )
        cc(i,j,2) = 0.5d0 * ( fy(i,j) + fy(i,j+1) )
     enddo
  enddo

end subroutine bl_avg_fc_to_cc
