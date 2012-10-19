   ! ::: -----------------------------------------------------------
 
   subroutine statefill(state,state_l1,state_l2,state_h1,state_h2,&
                        domlo,domhi,dx,xlo,time,bc)
 
     implicit none
     !include 'bc_types.fi'
 
     integer :: state_l1,state_l2,state_h1,state_h2
     integer :: bc(2,2,*)
     integer :: domlo(2), domhi(2)
     double precision dx(2), xlo(2), time
     double precision state(state_l1:state_h1,state_l2:state_h2,*)
     integer :: i

     !do i=1,NVAR
     !   call filcc(state(state_l1,state_l2,i),state_l1,state_l2,state_h1,state_h2,domlo,domhi,dx,xlo,bc)
     !enddo
     
   end subroutine statefill


! :: ----------------------------------------------------------
! :: Volume-weight average the fine grid data onto the coarse
! :: grid.  Overlap is given in coarse grid coordinates.
! ::
! :: INPUTS / OUTPUTS:
! ::  crse      <=  coarse grid data
! ::  clo,chi    => index limits of crse array interior
! ::  ngc        => number of ghost cells in coarse array
! ::  nvar	 => number of components in arrays
! ::  fine       => fine grid data
! ::  flo,fhi    => index limits of fine array interior
! ::  ngf        => number of ghost cells in fine array
! ::  rfine      => (ignore) used in 2-D RZ calc
! ::  lo,hi      => index limits of overlap (crse grid)
! ::  lrat       => refinement ratio
! ::
! :: NOTE:
! ::  Assumes all data cell centered
! :: ----------------------------------------------------------
! ::
      subroutine avgdown(crse,c_l1,c_l2,c_h1,c_h2,nvar, &
                         cv,cv_l1,cv_l2,cv_h1,cv_h2, &
                         fine,f_l1,f_l2,f_h1,f_h2, &
                         fv,fv_l1,fv_l2,fv_h1,fv_h2,lo,hi,lrat)
      implicit none
      integer c_l1,c_l2,c_h1,c_h2
      integer cv_l1,cv_l2,cv_h1,cv_h2
      integer f_l1,f_l2,f_h1,f_h2
      integer fv_l1,fv_l2,fv_h1,fv_h2
      integer lo(2), hi(2)
      integer nvar, lrat(2)
      double precision crse(c_l1:c_h1,c_l2:c_h2,nvar)
      double precision cv(cv_l1:cv_h1,cv_l2:cv_h2)
      double precision fine(f_l1:f_h1,f_l2:f_h2,nvar)
      double precision fv(fv_l1:fv_h1,fv_l2:fv_h2)

      integer i, j, n, ic, jc, ioff, joff
      integer lenx, leny, mxlen
      integer lratx, lraty

      lratx = lrat(1)
      lraty = lrat(2)
      lenx = hi(1)-lo(1)+1
      leny = hi(2)-lo(2)+1
      mxlen = max(lenx,leny)

      if (lenx .eq. mxlen) then
         do n = 1, nvar
 
!           Set coarse grid to zero on overlap
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,n) = 0.d0
               enddo
            enddo

!           Sum fine data
            do joff = 0, lraty-1
               do jc = lo(2), hi(2)
                  j = jc*lraty + joff
                  do ioff = 0, lratx-1
                     do ic = lo(1), hi(1)
                        i = ic*lratx + ioff
                        crse(ic,jc,n) = crse(ic,jc,n) + fv(i,j) * fine(i,j,n)
                     enddo
                  enddo
               enddo
            enddo

!           Divide out by volume weight
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,n) = crse(ic,jc,n) / cv(ic,jc)
               enddo
            enddo
            
         enddo

      else

         do n = 1, nvar

!           Set coarse grid to zero on overlap
            do ic = lo(1), hi(1)
               do jc = lo(2), hi(2)
                  crse(ic,jc,n) = 0.d0
               enddo
            enddo
 
!           Sum fine data
            do ioff = 0, lratx-1
               do ic = lo(1), hi(1)
                  i = ic*lratx + ioff
                  do joff = 0, lraty-1
                     do jc = lo(2), hi(2)
                        j = jc*lraty + joff
                        crse(ic,jc,n) = crse(ic,jc,n) + fv(i,j) * fine(i,j,n)
                     enddo
                  enddo
               enddo
            enddo
             
!           Divide out by volume weight
            do ic = lo(1), hi(1)
               do jc = lo(2), hi(2)
                  crse(ic,jc,n) = crse(ic,jc,n) / cv(ic,jc)
               enddo
            enddo
            
         enddo

      end if

      end subroutine avgdown

      subroutine derstate(state,state_l1,state_l2,state_h1,state_h2,nv,&
                          dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo,&
                          domhi,delta,xlo,time,dt,bc,level,grid_no)
!
!     This routine will derive the X from the (rho X)
!
      implicit none 

      integer          lo(2), hi(2)
      integer          state_l1,state_l2,state_h1,state_h2,nv
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision delta(2), xlo(2), time, dt
      double precision state(state_l1:state_h1,state_l2:state_h2,nv)
      double precision   dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no
 
      print *,'fix me'
      call bl_pd_abort()

      end subroutine derstate
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
      subroutine laplac_error(tag,tagl1,tagl2,tagh1,tagh2, &
                              set,clear, &
                              var,varl1,varl2,varh1,varh2, &
                              lo,hi,nd,domlo,domhi, &
                              delta,xlo,problo,time,level)
      implicit none

      integer          :: set, clear, nd, level
      integer          :: tagl1,tagl2,tagh1,tagh2
      integer          :: varl1,varl2,varh1,varh2
      integer          :: lo(2), hi(2), domlo(2), domhi(2)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2)
      double precision :: var(varl1:varh1,varl2:varh2,nd)
      double precision :: delta(2), xlo(2), problo(2), time
      integer          :: i,j

      ! This value is  taken from FLASH
      double precision, parameter :: ctore=0.8

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         if (var(i,j,1) .gt. ctore) tag(i,j)=set
      end do
      end do

      end subroutine laplac_error

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the density
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: spec      => species array
! ::: nd        => number of components in den array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine state_error(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             spec,specl1,specl2,spech1,spech2, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none
      integer set, clear, nd, level
      integer tagl1,tagl2,tagh1,tagh2
      integer specl1,specl2,spech1,spech2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision spec(specl1:spech1,specl2:spech2,nd)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

!     Tag on regions of high concentration
!     In this example we only test on whether the 2nd species has conc > specerr
      if (level .lt. max_specerr_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (spec(i,j,2) .ge. specerr) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif

!     Tag on regions of high gradient
      if (level .lt. max_specgrad_lev) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(spec(i+1,j,1) - spec(i,j,1))
               ay = ABS(spec(i,j+1,1) - spec(i,j,1))
               ax = MAX(ax,ABS(spec(i,j,1) - spec(i-1,j,1)))
               ay = MAX(ay,ABS(spec(i,j,1) - spec(i,j-1,1)))
               if ( MAX(ax,ay) .ge. specgrad) then
                  tag(i,j) = set
               endif
            enddo
         enddo
      endif
      
      end subroutine state_error
