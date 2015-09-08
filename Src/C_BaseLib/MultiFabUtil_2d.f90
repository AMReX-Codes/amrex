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

subroutine bl_avgdown_faces (lo, hi, &
     f, f_l1, f_l2, f_h1, f_h2, &
     c, c_l1, c_l2, c_h1, c_h2, &
     ratio, idir)

  implicit none
  integer          :: lo(2),hi(2)
  integer          :: f_l1, f_l2, f_h1, f_h2
  integer          :: c_l1, c_l2, c_h1, c_h2
  integer          :: ratio(3), idir
  double precision :: f(f_l1:f_h1, f_l2:f_h2)
  double precision :: c(c_l1:c_h1, c_l2:c_h2)

  ! Local variables
  integer i,j,n,facx,facy

  facx = ratio(1)
  facy = ratio(2)

  if (idir .eq. 0) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               c(i,j) = 0.d0
               do n = 0,facy-1
                  c(i,j) = c(i,j) + f(facx*i,facy*j+n)
               end do
               c(i,j) = c(i,j) / facy
            end do
         end do
  else 
         do i = lo(1), hi(1)
            do j = lo(2), hi(2)
               c(i,j) = 0.d0
               do n = 0,facx-1
                  c(i,j) = c(i,j) + f(facx*i+n,facy*j)
               end do
               c(i,j) = c(i,j) / facx
            end do
         end do
  end if

end subroutine bl_avgdown_faces
