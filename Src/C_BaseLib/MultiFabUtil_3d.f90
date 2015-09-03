subroutine bl_avg_fc_to_cc (lo, hi, &
     cc, ccl1, ccl2, ccl3, cch1, cch2, cch3, &
     fx, fxl1, fxl2, fxl3, fxh1, fxh2, fxh3, &
     fy, fyl1, fyl2, fyl3, fyh1, fyh2, fyh3, &
     fz, fzl1, fzl2, fzl3, fzh1, fzh2, fzh3, &
     dx, problo, coord_type)

  implicit none
  integer          :: lo(3),hi(3), coord_type
  integer          :: ccl1, ccl2, ccl3, cch1, cch2, cch3
  integer          :: fxl1, fxl2, fxl3, fxh1, fxh2, fxh3
  integer          :: fyl1, fyl2, fyl3, fyh1, fyh2, fyh3
  integer          :: fzl1, fzl2, fzl3, fzh1, fzh2, fzh3
  double precision :: cc(ccl1:cch1, ccl2:cch2, ccl3:cch3, 3)
  double precision :: fx(fxl1:fxh1, fxl2:fxh2, fxl3:fxh3)
  double precision :: fy(fyl1:fyh1, fyl2:fyh2, fyl3:fyh3)
  double precision :: fz(fzl1:fzh1, fzl2:fzh2, fzl3:fzh3)
  double precision :: dx(3), problo(3)

  ! Local variables
  integer          :: i,j,k
  
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        do i=lo(1),hi(1)
           cc(i,j,k,1) = 0.5d0 * ( fx(i,j,k) + fx(i+1,j,k) )
           cc(i,j,k,2) = 0.5d0 * ( fy(i,j,k) + fy(i,j+1,k) )
           cc(i,j,k,3) = 0.5d0 * ( fz(i,j,k) + fz(i,j,k+1) )
        enddo
     enddo
  enddo

end subroutine bl_avg_fc_to_cc
