subroutine bl_avg_fc_to_cc (lo, hi, &
     cc, ccl1, cch1, &
     fx, fxl1, fxh1, &
     dx, problo, coord_type)

  implicit none
  integer          :: lo(1),hi(1), coord_type
  integer          :: ccl1, cch1
  integer          :: fxl1, fxh1
  integer          :: fyl1, fyh1
  double precision :: cc(ccl1:cch1)
  double precision :: fx(fxl1:fxh1)
  double precision :: dx(1), problo(1)

  ! Local variables
  integer          :: i
  double precision :: rlo,rhi,rcen

  if (coord_type .eq. 0) then ! Cartesian

     do i=lo(1),hi(1)
        cc(i) = 0.5d0 * ( fx(i) + fx(i+1) )
     enddo
  
  else if (coord_type .eq. 1) then ! Cylindrical

     do i=lo(1),hi(1)
        rlo = problo(1) + (dble(i)  )*dx(1)
        rhi = problo(1) + (dble(i+1))*dx(1)
        rcen = 0.5d0 * (rlo + rhi)
        cc(i) = 0.5d0 * ( rlo*fx(i) + rhi*fx(i+1) ) / rcen
     enddo
  
  else  ! Spherical

     do i=lo(1),hi(1)
        rlo = problo(1) + (dble(i)  )*dx(1)
        rhi = problo(1) + (dble(i+1))*dx(1)
        rcen = 0.5d0 * (rlo + rhi)
        cc(i) = 0.5d0 * ( rlo**2 * fx(i) + rhi**2 * fx(i+1) ) / rcen**2
     enddo

  end if

end subroutine bl_avg_fc_to_cc

subroutine bl_avg_cc_to_fc (lo, hi, &
     fx, fxl1, fxh1, &
     cc, ccl1, cch1, &
     dx, problo, coord_type)

  implicit none
  integer          :: lo(1),hi(1), coord_type
  integer          :: ccl1, cch1
  integer          :: fxl1, fxh1
  integer          :: fyl1, fyh1
  double precision :: cc(ccl1:cch1)
  double precision :: fx(fxl1:fxh1)
  double precision :: dx(1), problo(1)

  ! Local variables
  integer          :: i
  double precision :: rlo,rhi,rcen

  if (coord_type .eq. 0) then ! Cartesian

     do i = lo(1), hi(1)
        fx(i) = 0.5d0 * (cc(i-1) + cc(i))
     end do

  else if (coord_type .eq. 1) then ! Cylindrical

     do i = lo(1), hi(1)
        rlo = problo(1) + (dble(i-1)+0.5d0)*dx(1)
        rhi = problo(1) + (dble(i  )+0.5d0)*dx(1)
        rcen = 0.5d0 * (rlo + rhi)
        fx(i) = 0.5d0 * (rlo*cc(i-1) + rhi*cc(i)) / rcen
     end do

  else  ! Spherical

     do i = lo(1), hi(1)
        rlo = problo(1) + (dble(i-1)+0.5d0)*dx(1)
        rhi = problo(1) + (dble(i  )+0.5d0)*dx(1)
        rcen = 0.5d0 * (rlo + rhi)
        fx(i) = 0.5d0 * (rlo**2 * cc(i-1) + rhi**2 * cc(i)) / rcen**2
     end do

  end if

end subroutine bl_avg_cc_to_fc

subroutine bl_avgdown_faces (lo, hi, &
     f, f_l1, f_h1, &
     c, c_l1, c_h1, &
     ratio,n)

  implicit none
  integer          :: lo(1),hi(1)
  integer          :: f_l1, f_h1
  integer          :: c_l1, c_h1
  integer          :: ratio(1), n
  double precision :: f(f_l1:f_h1)
  double precision :: c(c_l1:c_h1)

  ! Local variables
  integer i,facx

  facx = ratio(1)

   ! lo(1)..hi(1) are edge base indices
   do i = lo(1), hi(1)
      c(i) = f(facx*i)
   end do

end subroutine bl_avgdown_faces

