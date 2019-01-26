! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the Laplacian.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: var       => array of data
! ::: nd        => number of components in var array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine tag_laplac_error(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                                  set,clear, &
                                  var,varl1,varl2,varl3,varh1,varh2,varh3, &
                                  lo,hi,nd,domlo,domhi, &
                                  delta,xlo,problo,time,level)
      use probdata_module
      use amrex_fort_module, only : rt => amrex_real
      implicit none

      integer          :: set, clear, nd, level
      integer          :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer          :: varl1,varl2,varl3,varh1,varh2,varh3
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) :: var(varl1:varh1,varl2:varh2,varl3:varh3)
      real(rt) :: delta(3), xlo(3), problo(3), time
      integer          :: i,j,k

      real(rt) ::  delu(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3)
      real(rt) :: delua(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3)
      real(rt) :: delu2(9), delu3(9), delu4(9)
      real(rt) :: num, denom, error

      ! This value is  taken from FLASH
      real(rt), parameter :: ctore=0.8
      real(rt), parameter:: epsil=0.02

      ! adapted from ref_marking.f90 in FLASH2.5

      ! d/dx
      do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(i,j,k,1) =     var(i+1,j,k)  -     var(i-1,j,k) 
         delua(i,j,k,1) = abs(var(i+1,j,k)) + abs(var(i-1,j,k))
      end do
      end do
      end do

      ! d/dy
      do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(i,j,k,2) =     var(i,j+1,k)  -     var(i,j-1,k) 
         delua(i,j,k,2) = abs(var(i,j+1,k)) + abs(var(i,j-1,k))
      end do
      end do
      end do

      ! d/dz
      do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
      do i=lo(1)-1,hi(1)+1
          delu(i,j,k,3) =     var(i,j,k+1)  -     var(i,j,k-1)
         delua(i,j,k,3) = abs(var(i,j,k+1)) + abs(var(i,j,k-1))
      end do
      end do
      end do

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         

         ! d/dxdx
          delu2(1) =     delu(i+1,j,k,1)  -     delu(i-1,j,k,1)
          delu3(1) = abs(delu(i+1,j,k,1)) + abs(delu(i-1,j,k,1))
          delu4(1) =    delua(i+1,j,k,1)  +    delua(i-1,j,k,1)

          ! d/dydx
          delu2(2) =     delu(i,j+1,k,1)  -     delu(i,j-1,k,1)
          delu3(2) = abs(delu(i,j+1,k,1)) + abs(delu(i,j-1,k,1))
          delu4(2) =    delua(i,j+1,k,1)  +    delua(i,j-1,k,1)

          ! d/dxdy
          delu2(3) =     delu(i+1,j,k,2)  -     delu(i-1,j,k,2)
          delu3(3) = abs(delu(i+1,j,k,2)) + abs(delu(i-1,j,k,2))
          delu4(3) =    delua(i+1,j,k,2)  +    delua(i-1,j,k,2)
                                       
          ! d/dydy                     
          delu2(4) =     delu(i,j+1,k,2)  -     delu(i,j-1,k,2)
          delu3(4) = abs(delu(i,j+1,k,2)) + abs(delu(i,j-1,k,2))
          delu4(4) =    delua(i,j+1,k,2)  +    delua(i,j-1,k,2)
                                                              
          ! d/dzdx                                            
          delu2(5) =     delu(i,j,k+1,1)  -     delu(i,j,k-1,1)
          delu3(5) = abs(delu(i,j,k+1,1)) + abs(delu(i,j,k-1,1))
          delu4(5) =    delua(i,j,k+1,1)  +    delua(i,j,k-1,1)
                                                              
          ! d/dzdy                                            
          delu2(6) =     delu(i,j,k+1,2)  -     delu(i,j,k-1,2)
          delu3(6) = abs(delu(i,j,k+1,2)) + abs(delu(i,j,k-1,2))
          delu4(6) =    delua(i,j,k+1,2)  +    delua(i,j,k-1,2)
                                                              
          ! d/dxdz                                            
          delu2(7) =     delu(i+1,j,k,3)  -     delu(i-1,j,k,3)
          delu3(7) = abs(delu(i+1,j,k,3)) + abs(delu(i-1,j,k,3))
          delu4(7) =    delua(i+1,j,k,3)  +    delua(i-1,j,k,3)
                                                              
          ! d/dydz                                            
          delu2(8) =     delu(i,j+1,k,3)  -     delu(i,j-1,k,3)
          delu3(8) = abs(delu(i,j+1,k,3)) + abs(delu(i,j-1,k,3))
          delu4(8) =    delua(i,j+1,k,3)  +    delua(i,j-1,k,3)
                                                              
          ! d/dzdz                                            
          delu2(9) =     delu(i,j,k+1,3)  -     delu(i,j,k-1,3)
          delu3(9) = abs(delu(i,j,k+1,3)) + abs(delu(i,j,k-1,3))
          delu4(9) =    delua(i,j,k+1,3)  +    delua(i,j,k-1,3)

         ! compute the error
         num   =  delu2(1)**2 + delu2(2)**2 + delu2(3)**2 + delu2(4)**2 &
                 +delu2(5)**2 + delu2(6)**2 + delu2(7)**2 + delu2(8)**2 &
                 +delu2(9)**2

         denom = (delu3(1) + (epsil*delu4(1)+1.d-99))**2 + &
                 (delu3(2) + (epsil*delu4(2)+1.d-99))**2 + &
                 (delu3(3) + (epsil*delu4(3)+1.d-99))**2 + &
                 (delu3(4) + (epsil*delu4(4)+1.d-99))**2 + &
                 (delu3(5) + (epsil*delu4(5)+1.d-99))**2 + &
                 (delu3(6) + (epsil*delu4(6)+1.d-99))**2 + &
                 (delu3(7) + (epsil*delu4(7)+1.d-99))**2 + &
                 (delu3(8) + (epsil*delu4(8)+1.d-99))**2 + &
                 (delu3(9) + (epsil*delu4(9)+1.d-99))**2

         error = sqrt(num/denom)

         if (error .gt. ctore) tag(i,j,k)=set

      end do
      end do
      end do

      end subroutine tag_laplac_error

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the density
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: den       => density array
! ::: nd        => number of components in den array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine tag_denerror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                              set,clear, &
                              den,denl1,denl2,denl3,denh1,denh2,denh3, &
                              lo,hi,nd,domlo,domhi, &
                              delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer          :: set, clear, nd, level
      integer          :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer          :: denl1,denl2,denl3,denh1,denh2,denh3
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) :: den(denl1:denh1,denl2:denh2,denl3:denh3,nd)
      real(rt) :: delta(3), xlo(3), problo(3), time

      integer          :: i,j,k
      real(rt) :: ax,ay,az

!     Tag on regions of high density
      if (level .lt. max_denerr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (den(i,j,k,1) .ge. denerr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high density gradient
      if (level .lt. max_dengrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(den(i+1,j,k,1) - den(i,j,k,1))
                  ay = ABS(den(i,j+1,k,1) - den(i,j,k,1))
                  az = ABS(den(i,j,k+1,1) - den(i,j,k,1))
                  ax = MAX(ax,ABS(den(i,j,k,1) - den(i-1,j,k,1)))
                  ay = MAX(ay,ABS(den(i,j,k,1) - den(i,j-1,k,1)))
                  az = MAX(az,ABS(den(i,j,k,1) - den(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. dengrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end subroutine tag_denerror

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the temperature
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: temp      => temperature array
! ::: np        => number of components in temp array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------

      subroutine tag_temperror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                               set,clear, &
                               temp,templ1,templ2,templ3,temph1, &
                               temph2,temph3, &
                               lo,hi,np,domlo,domhi, &
                               delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer          :: set, clear, np, level
      integer          :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer          :: templ1,templ2,templ3,temph1,temph2,temph3
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) :: temp(templ1:temph1,templ2:temph2, &
                            templ3:temph3,np)
      real(rt) :: delta(3), xlo(3), problo(3), time

      integer          :: i,j,k
      real(rt) :: ax,ay,az

!     Tag on regions of high temperature
      if (level .lt. max_temperr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (temp(i,j,k,1) .ge. temperr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high temperature gradient
      if (level .lt. max_tempgrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(temp(i+1,j,k,1) - temp(i,j,k,1))
                  ay = ABS(temp(i,j+1,k,1) - temp(i,j,k,1))
                  az = ABS(temp(i,j,k+1,1) - temp(i,j,k,1))
                  ax = MAX(ax,ABS(temp(i,j,k,1) - temp(i-1,j,k,1)))
                  ay = MAX(ay,ABS(temp(i,j,k,1) - temp(i,j-1,k,1)))
                  az = MAX(az,ABS(temp(i,j,k,1) - temp(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. tempgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end subroutine tag_temperror

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the pressure
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: press     => pressure array
! ::: np        => number of components in press array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine tag_presserror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                                set,clear, &
                                press,pressl1,pressl2,pressl3,pressh1, &
                                pressh2,pressh3, &
                                lo,hi,np,domlo,domhi, &
                                delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer          :: set, clear, np, level
      integer          :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer          :: pressl1,pressl2,pressl3,pressh1,pressh2,pressh3
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) :: press(pressl1:pressh1,pressl2:pressh2, &
                                pressl3:pressh3,np)
      real(rt) :: delta(3), xlo(3), problo(3), time

      integer          :: i,j,k
      real(rt) :: ax,ay,az

!     Tag on regions of high pressure
      if (level .lt. max_presserr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (press(i,j,k,1) .ge. presserr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high pressure gradient
      if (level .lt. max_pressgrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(press(i+1,j,k,1) - press(i,j,k,1))
                  ay = ABS(press(i,j+1,k,1) - press(i,j,k,1))
                  az = ABS(press(i,j,k+1,1) - press(i,j,k,1))
                  ax = MAX(ax,ABS(press(i,j,k,1) - press(i-1,j,k,1)))
                  ay = MAX(ay,ABS(press(i,j,k,1) - press(i,j-1,k,1)))
                  az = MAX(az,ABS(press(i,j,k,1) - press(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. pressgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end subroutine tag_presserror

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the velocity
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: vel       => velocity array
! ::: nv        => number of components in vel array (should be 3)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine tag_velerror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                              set,clear, &
                              vel,vell1,vell2,vell3,velh1,velh2,velh3, &
                              lo,hi,nv,domlo,domhi, &
                              delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer          :: set, clear, nv, level
      integer          :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer          :: vell1,vell2,vell3,velh1,velh2,velh3
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) :: vel(vell1:velh1,vell2:velh2,vell3:velh3,nv)
      real(rt) :: delta(3), xlo(3), problo(3), time

      integer          :: i,j,k
      real(rt) :: ax,ay,az

!     Tag on regions of high velocity gradient
      if (level .lt. max_velgrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(vel(i+1,j,k,1) - vel(i,j,k,1))
                  ay = ABS(vel(i,j+1,k,1) - vel(i,j,k,1))
                  az = ABS(vel(i,j,k+1,1) - vel(i,j,k,1))
                  ax = MAX(ax,ABS(vel(i,j,k,1) - vel(i-1,j,k,1)))
                  ay = MAX(ay,ABS(vel(i,j,k,1) - vel(i,j-1,k,1)))
                  az = MAX(az,ABS(vel(i,j,k,1) - vel(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. velgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end subroutine tag_velerror

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on overdensity
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: den       => density array
! ::: nc        => number of components in density array
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine tag_overdensity(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                                 set,clear, &
                                 den,denl1,denl2,denl3,denh1,denh2,denh3, &
                                 lo,hi,nc,domlo,domhi,delta,level,avg_den)

      use probdata_module
      implicit none

      integer          :: set, clear, nc, level
      integer          :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer          :: denl1,denl2,denl3,denh1,denh2,denh3
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) :: den(denl1:denh1,denl2:denh2,denl3:denh3,nc)
      real(rt) :: delta(3), avg_den

      integer          :: i,j,k
      real(rt) :: over_den

      over_den = 1.1d0 * avg_den

!     Tag on regions of overdensity
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if ( den(i,j,k,1) .gt. over_den ) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
      enddo

      end subroutine tag_overdensity

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the number of particles.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: var       => array of data
! ::: nd        => number of components in var array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine tag_part_cnt_err(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                                  set,clear, &
                                  var,varl1,varl2,varl3,varh1,varh2,varh3, &
                                  lo,hi,nd,domlo,domhi, &
                                  delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer          :: set, clear, nd, level
      integer          :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer          :: varl1,varl2,varl3,varh1,varh2,varh3
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) :: var(varl1:varh1,varl2:varh2,varl3:varh3,nd)
      real(rt) :: delta(3), xlo(3), problo(3), time

      integer          :: ilo,jlo,klo,ihi,jhi,khi
      integer          :: i,j,k

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         if (var(i,j,k,1) .gt. max_num_part) then
            tag(i,j,k) = set
         end if
      end do
      end do
      end do

      end subroutine tag_part_cnt_err

! ::: -----------------------------------------------------------
! ::: This routine will tag a specific region.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: var       => array of data
! ::: nd        => number of components in var array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine tag_region(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                            set,lo,hi,domlo,domhi, &
                            delta,xlo,problo,level)

      use probdata_module
      implicit none

      integer          :: set,level
      integer          :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      real(rt) :: delta(3), xlo(3), problo(3), time

      integer          :: ilo,jlo,klo,ihi,jhi,khi
      integer          :: i,j,k

      if (level .eq. 0) then
         ilo = (domhi(1)+1)*1/4
         ihi = (domhi(1)+1)*3/4 - 1
         jlo = (domhi(2)+1)*1/4
         jhi = (domhi(2)+1)*3/4 - 1
         klo = (domhi(3)+1)*1/4
         khi = (domhi(3)+1)*3/4 - 1
         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            if (i.ge.ilo .and. i.le.ihi .and. &
                j.ge.jlo .and. j.le.jhi .and. &
                k.ge.klo .and. k.le.khi) then
               tag(i,j,k) = set
            end if
         end do
         end do
         end do
      else
         ilo = (domhi(1)+1)*3/8
         ihi = (domhi(1)+1)*5/8 - 1
         jlo = (domhi(2)+1)*3/8
         jhi = (domhi(2)+1)*5/8 - 1
         klo = (domhi(3)+1)*3/8
         khi = (domhi(3)+1)*5/8 - 1
         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            if (i.ge.ilo .and. i.le.ihi .and. &
                j.ge.jlo .and. j.le.jhi .and. &
                k.ge.klo .and. k.le.khi) then
               tag(i,j,k) = set
            end if
         end do
         end do
         end do
      endif

      end subroutine tag_region
