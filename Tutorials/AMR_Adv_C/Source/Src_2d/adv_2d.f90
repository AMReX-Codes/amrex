
subroutine advect(time, lo, hi, &
     &            uin , ui_lo, ui_hi, &
     &            uout, uo_lo, uo_hi, &
     &            vx  , vx_lo, vx_hi, &
     &            vy  , vy_lo, vy_hi, &
     &            flxx, fx_lo, fx_hi, &
     &            flxy, fy_lo, fy_hi, &
     &            dx,dt)
  
  use mempool_module, only : bl_allocate, bl_deallocate
  use trace_module, only : tracex, tracey
  use flux_module, only : cmpflx
  use transverse_module, only : transy, transx

  implicit none

  integer, intent(in) :: lo(2), hi(2)
  double precision, intent(in) :: dx(2), dt, time
  integer, intent(in) :: ui_lo(2), ui_hi(2)
  integer, intent(in) :: uo_lo(2), uo_hi(2)
  integer, intent(in) :: vx_lo(2), vx_hi(2)
  integer, intent(in) :: vy_lo(2), vy_hi(2)
  integer, intent(in) :: fx_lo(2), fx_hi(2)
  integer, intent(in) :: fy_lo(2), fy_hi(2)
  double precision, intent(in   ) :: uin (ui_lo(1):ui_hi(1),ui_lo(2):ui_hi(2))
  double precision, intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2))
  double precision, intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2))
  double precision, intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2))
  double precision, intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2))
  double precision, intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2))

  integer :: i, j
  integer :: glo(2), ghi(2)
  double precision :: dtdx(2), hdtdx(2)

  ! Some compiler may not support 'contiguous'.  Remove it in that case.
  double precision, dimension(:,:), pointer, contiguous :: &
       qxm, qxp, qym, qyp, fxtmp, fytmp, dq

  dtdx = dt/dx
  hdtdx = dtdx * 0.5d0

  glo = lo - 1
  ghi = hi + 1

  ! We use BoxLib's bl_allocate to allocate memeory 
  ! instead of intrinsic allocate because it is faster inside OMP.
  ! One must remember to call bl_deallocate.

  ! *x* is on x-face and y-center
  ! *y* is on y-face and x-center
  call bl_allocate(qxm  , glo(1), ghi(1), glo(2), ghi(2))
  call bl_allocate(qxp  , glo(1), ghi(1), glo(2), ghi(2))
  call bl_allocate(qym  , glo(1), ghi(1), glo(2), ghi(2))
  call bl_allocate(qyp  , glo(1), ghi(1), glo(2), ghi(2))
  call bl_allocate(fxtmp, glo(1), ghi(1), glo(2), ghi(2))
  call bl_allocate(fytmp, glo(1), ghi(1), glo(2), ghi(2))
  ! dq is cell centered
  call bl_allocate(dq   , glo(1), ghi(1), glo(2), ghi(2))

  ! Trace to x-faces
  ! scratch space: dq       on (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1) 
  ! Input        : uin      on (lo(1)-3:hi(1)+3,lo(2)-2:hi(2)+2) 
  !                vx       on (lo(1)  :hi(1)+1,lo(2)-1:hi(2)+1)
  ! Output       : qxm, qxp on (lo(1)  :hi(1)+1,lo(2)-1:hi(2)+1)
  call tracex((/lo(1),lo(2)-1/), ghi, hdtdx(1), &
              dq, glo, ghi, &
              uin, ui_lo, ui_hi, &
              vx, vx_lo, vx_hi, &
              qxm, qxp, glo, ghi)

  ! Trace to y-faces
  ! scratch space: dq       on (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)
  ! Input        : uin      on (lo(1)-2:hi(1)+2,lo(2)-3:hi(2)+3)
  !                vy       on (lo(1)-1:hi(1)+1,lo(2):hi(2)+1)
  ! Output       : qym, qyp on (lo(1)-1:hi(1)+1,lo(2):hi(2)+1)
  call tracey((/lo(1)-1,lo(2)/), ghi, hdtdx(2), &
              dq, glo, ghi, &
              uin, ui_lo, ui_hi, &
              vy, vy_lo, vy_hi, &
              qym, qyp, glo, ghi)

  ! Compute temporary x-flux
  ! Input : qxm, qxp, vx on (lo(1):hi(1)+1,lo(2)-1:hi(2)+1)
  ! Output: fxtmp        on (lo(1):hi(1)+1,lo(2)-1:hi(2)+1)
  call cmpflx((/lo(1),lo(2)-1/), ghi, &
              qxm, qxp, glo, ghi, &
              vx, vx_lo, vx_hi, &
              fxtmp, glo, ghi)

  ! Compute temporary y-flux
  ! Input : qym, qyp, vy on (lo(1)-1:hi(1)+1,lo(2):hi(2)+1)
  ! Output: fytmp        on (lo(1)-1:hi(1)+1,lo(2):hi(2)+1)
  call cmpflx((/lo(1)-1,lo(2)/), ghi, &
              qym, qyp, glo, ghi, &
              vy, vy_lo, vy_hi, &
              fytmp, glo, ghi)

  ! Evolve states on x-faces with transverse contribution from y-fluxes
  ! Input   : fytmp    on (lo(1)-1:hi(1)+1,lo(2):hi(2)+1)
  ! In & Out: qxm, qxp on (lo(1)  :hi(1)+1,lo(2):hi(2))
  call transy(lo, (/hi(1)+1,hi(2)/), hdtdx(2), &
              qxm, qxp, glo, ghi, &
              fytmp, glo, ghi)

  ! Compute final fluxes on x-faces
  ! Input : qxm, qxp, vx on (lo(1):hi(1)+1,lo(2):hi(2))
  ! Output: flxx         on (lo(1):hi(1)+1,lo(2):hi(2))
  call cmpflx(lo, (/hi(1)+1,hi(2)/), &
              qxm, qxp, glo, ghi, &
              vx, vx_lo, vx_hi, &
              flxx, fx_lo, fx_hi)

  ! Evolve states on y-faces with transverse contribution from x-fluxes
  ! Input   : fxtmp    on (lo(1):hi(1)+1,lo(2)-1:hi(2)+1)
  ! In & Out: qym, qyp on (lo(1):hi(1)  ,lo(2)  :hi(2)+1)
  call transx(lo, (/hi(1),hi(2)+1/), hdtdx(1), &
              qym, qyp, glo, ghi, &
              fxtmp, glo, ghi)

  ! Compute final fluxes on y-faces
  ! Input : qym, qyp, vy on (lo(1):hi(1),lo(2):hi(2)+1)
  ! Output: flxy         on (lo(1):hi(1),lo(2):hi(2)+1)
  call cmpflx(lo, (/hi(1),hi(2)+1/), &
              qym, qyp, glo, ghi, &
              vy, vy_lo, vy_hi, &
              flxy, fy_lo, fy_hi)

  ! Do a conservative update
  do    j = lo(2),hi(2)
     do i = lo(1),hi(1)
        uout(i,j) = uin(i,j) + &
             ( (flxx(i,j) - flxx(i+1,j)) * dtdx(1) &
             + (flxy(i,j) - flxy(i,j+1)) * dtdx(2) )
     enddo
  enddo

  ! Scale by face area in order to correctly reflx
  do    j = lo(2),hi(2)
     do i = lo(1),hi(1)+1
        flxx(i,j) = flxx(i,j) * ( dt * dx(2))
     enddo
  enddo
  
  ! Scale by face area in order to correctly reflx
  do    j = lo(2),hi(2)+1 
     do i = lo(1),hi(1)
        flxy(i,j) = flxy(i,j) * (dt * dx(1))
     enddo
  enddo

  call bl_deallocate(qxm)
  call bl_deallocate(qxp)
  call bl_deallocate(qym)
  call bl_deallocate(qyp)
  call bl_deallocate(fxtmp)
  call bl_deallocate(fytmp)
  call bl_deallocate(dq)

end subroutine advect
