#define USE_TMP 1
#include <REAL.H>

#include "DivVis_F.H"
#include <LO_BCTYPES.H>
#include "ArrayLim.H"

c :::: usage:
c :::: this .mF file is meant to be run through Mathematica.  This converts
c :::: compact symbolic expressions into fortran which is stored in a .F
c :::: file.

c----------------------------------------------------------------
c     this is the fortran support file for the the operator 
c     L(U) = alpha*a(x)*U - beta*Div( tau )
c
c     where U is the two component vector (u,v) and
c     tau is a three by three tensor
c     tau = | t_xx     t_xy 	t_xz|
c           | t_xy     t_yy 	t_yz|
c	    | t_xz     t_yz	t_zz|
c
c     t_xx = 2*mu * u_x
c     t_yy = 2*mu * v_y
c     t_zz = 2*mu * w_z
c     t_xy = mu*(u_y + v_x)
c     t_xz = mu*(u_z + w_x)
c     t_yz = mu*(v_z + w_y)

c ::: define standard replacements used by Mathematica
c ::: see file visc3d.ma
c ::: Null

c ::: interface notes:
c ::: 1) trander* ALWAYS have values in them, even if the cells are
c :::    all covered.  Might as well.  These are edge located derivative
c :::    values.  In index space, they are located the same as the mask
c :::    cells.  Which is somewhat anomalous.  These are edge values, after
c :::    all.  While the masks are cell-centered.  But it seems easier
c :::    at the moment.
c ::: 1a) trander is 3x3 for each location; first is component number, then
c :::     derivative direction.  Not all of these are used, in fact, only the
c :::     component which is normal to the face is used, and only the 
c :::     tangential directions are used.
c ::: 2) the normal derivatives are evaluated in the normal fashion.
c ::: 3) tangential derivatives which reach outside the rectangle DO have
c :::    to check the masks.  

      subroutine FORT_APPLYBC (
     $     flagden, flagbc, maxorder,
     $     u, DIMS(u),
     $     cdir, bct, bcl,
     $     bcval, DIMS(bcval),
     $     maskn, DIMS(maskn),
     $     maske, DIMS(maske),
     $     maskw, DIMS(maskw),
     $     masks, DIMS(masks),
     $     maskt, DIMS(maskt),
     $     maskb, DIMS(maskb),
     $     den, DIMS(den),
     $     exttd,DIMS(exttd),
     $     trander,DIMS(trander),
     $     lo, hi, nc,
     $     h
     $     )
c
c     If the boundary is of Neumann type, set the ghost cell value to
c     that of the outermost point in the valid data (2nd order accurate)
c     and then fill the "den" array with the value "1"
c     
c     
c     If flagbc==1:
c     
c     If the boundary is of Dirichlet type, construct a polynomial
c     interpolation through the boundary location and internal points
c     (at locations x(-1:len-2) that generates the ghost cell value (at
c     location xInt).  Then fill the ghost cell with the interpolated value.
c     If flagden==1, load the "den" array with the interpolation
c     coefficient corresponding to outermost point in the valid region
c     ( the coef(0) corresponding to the location x(0) )
c
c     Note: 
c     The bc type = LO_REFLECT_ODD is a special type of dirichlet condition,
c     in that we want a "zeroth" order interpolant to fill the ghost cell.
c     If this were treated in the normal way, then ALL boundaries would be
c     low order.

c ::: other notes since previous developers didn't bother to document
c ::: cdir is mnemonic for coordinate direction. i.e. which side is
c ::: cdir==0->west
c ::: cdir==3->east
c ::: cdir==1->south
c ::: cdir==4->north
c ::: cdir==2->bottom
c ::: cdir==5->top
c      
      integer maxorder
      integer nc, cdir, flagden, flagbc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(u)
      REAL_T u(DIMV(u),nc)
      integer DIMDEC(den)
      REAL_T den(DIMV(den),nc)
      integer DIMDEC(exttd)
      REAL_T exttd(DIMV(exttd),nc,3)
      integer DIMDEC(bcval)
      REAL_T bcval(DIMV(bcval),nc)
      integer DIMDEC(maskn)
      integer maskn(DIMV(maskn))
      integer DIMDEC(maske)
      integer maske(DIMV(maske))
      integer DIMDEC(maskw)
      integer maskw(DIMV(maskw))
      integer DIMDEC(masks)
      integer masks(DIMV(masks))
      integer DIMDEC(maskt)
      integer maskt(DIMV(maskt))
      integer DIMDEC(maskb)
      integer maskb(DIMV(maskb))
      integer DIMDEC(trander)
      REAL_T trander(DIMV(trander),3,3)
      integer bct(nc)
      REAL_T bcl
      REAL_T h(BL_SPACEDIM)
c
      integer i
      integer j
      integer k
      integer n
      logical is_dirichlet
      logical is_neumann
      logical is_odd
      REAL_T xb
      REAL_T innder,outder,innloc
      REAL_T edgloc,outloc
      REAL_T lambda
      REAL_T hx,hy,hz
c
      integer lenx
      integer leny
      integer lenz
      integer m
c...........trial change .................
      integer tmplen
      parameter(tmplen=256)
      REAL_T tmp(0:tmplen)
c...........trial change .................
c
      integer Lmaxorder
      integer maxmaxorder
      parameter(maxmaxorder=4)
      REAL_T x(-1:maxmaxorder-2)
      REAL_T coef(-1:maxmaxorder-2)
      REAL_T xInt
      logical False, True
      parameter( False=.false.)
      parameter( True=.true. )
c
      is_dirichlet(i) = ( i .eq. LO_DIRICHLET )
      is_neumann(i) = (i .eq. LO_NEUMANN)
c
      if ( maxorder .eq. -1 ) then
         Lmaxorder = maxmaxorder
      else
         Lmaxorder = MIN(maxorder,maxmaxorder)
      endif
      lenx = MIN(hi(1)-lo(1), Lmaxorder-2)
      leny = MIN(hi(2)-lo(2), Lmaxorder-2)
      lenz = MIN(hi(3)-lo(3), Lmaxorder-2)
      hx = h(1)
      hy = h(2)
      hz = h(3)
c...........trial change .................
      if( (hi(1)-lo(1)+1).gt.tmplen .or. 
     &     (hi(2)-lo(2)+1).gt.tmplen .or. 
     &     (hi(3)-lo(3)+1).gt.tmplen) then
         write(6,*)'FORT_APPLYBC: tmplen must be larger'
         stop
      endif
c...........trial change .................

c
c     TODO:
c     In order for this to work with growing multigrid, must
c     sort xa[] because it is possible for the xb value to lay
c     within this range.

c     
c     The west face of the grid
c
      if(cdir .eq. 0) then
         do n = 1,nc
            if (is_neumann(bct(n))) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     u(lo(1)-1,j,k,n) = cvmgt(
     $                    u(lo(1),j,k,n),
     $                    u(lo(1)-1,j,k,n),
     $                    maskw(lo(1)-1,j,k) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1) then
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        den(lo(1),j,k,n) = 1.0
                     enddo
                  enddo
               endif
c ::: ::: now do the tangential derivative part
c ::: ::: Null
               i = lo(1)
c ::: ::: ::: Y
               do k=lo(3),hi(3)
c ::: ::: ::: ::: interior part of west face
                  do j=lo(2)+1,hi(2)-1
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                     trander(i-1,j,k,n,2) = innder
                  enddo
                  j = lo(2)
                  if(masks(i,-1+j,k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,1+j,k,n)-U(i,2+j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif
                  trander(i-1,j,k,n,2) = innder
                  j = hi(2)
                  if(maskn(i,1+j,k).gt.0) then
                     innder = (U(i,-2+j,k,n)-4*U(i,-1+j,k,n)+3*U(i,j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif
                  trander(i-1,j,k,n,2) = innder
               enddo
c ::: ::: ::: Z
               do j=lo(2),hi(2)
                  do k=lo(3)+1,hi(3)-1
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                     trander(i-1,j,k,n,3) = innder
                  enddo
                  k = lo(3)
                  if(maskb(i,j,-1+k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,j,1+k,n)-U(i,j,2+k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  trander(i-1,j,k,n,3) = innder
                  k = hi(3)
                  if(maskt(i,j,1+k).gt.0) then
                     innder = (U(i,j,-2+k,n)-4*U(i,j,-1+k,n)+3*U(i,j,k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  trander(i-1,j,k,n,3) = innder
               enddo
            else if (is_dirichlet(bct(n))) then
               do m=0,lenx
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(1)
               xInt = - 0.5
               call polyInterpCoeff(xInt, x, lenx+2, coef)
#if USE_TMP
               do k = lo(3), hi(3)
                  if ( flagbc .eq. 1 ) then
                     do j = lo(2), hi(2)
                        tmp(j-lo(2)) = bcval(lo(1)-1,j,k,n)*coef(-1)
                     enddo
                  else
                     do j = lo(2), hi(2)
                        tmp(j-lo(2)) = 0.0
                     enddo
                  endif
                  do m = 0, lenx
                     do j = lo(2), hi(2)
                        tmp(j-lo(2)) = tmp(j-lo(2))+u(lo(1)+m,j,k,n)*coef(m)
                     enddo
                  enddo
                  do j = lo(2), hi(2)
                     u(lo(1)-1, j, k, n) = cvmgt(
     &                    tmp(j-lo(2)),
     $                    u(lo(1)-1, j,k, n),
     $                    maskw(lo(1)-1,j,k) .gt. 0)
                  enddo
               enddo
#else
               if ( flagbc .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        u(lo(1)-1, j, k, n) = cvmgt(
     $                       bcval(lo(1)-1,j,k,n)*coef(-1),
     $                       u(lo(1)-1, j,k, n),
     $                       maskw(lo(1)-1,j,k) .gt. 0)
                     enddo
                  enddo
               else
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        u(lo(1)-1, j, k, n) = cvmgt(
     $                       0.0,
     $                       u(lo(1)-1, j, k, n),
     $                       maskw(lo(1)-1,j, k) .gt. 0)
                     enddo
                  enddo
               endif
               do k = lo(3), hi(3)
                  do m = 0, lenx
                     do j = lo(2), hi(2)
                        u(lo(1)-1,j,k,n) = cvmgt(
     $                       u(lo(1)-1,j,k,n)
     $                       + u(lo(1)+m, j, k, n)*coef(m),
     $                       u(lo(1)-1,j,k,n),
     $                       maskw(lo(1)-1,j,k) .gt. 0)
                     enddo
                  enddo
               enddo
#endif
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        den(lo(1),j,k,n) = cvmgt(coef(0), 0.0,
     $                       maskw(lo(1)-1,j,k) .gt. 0)
                     enddo
                  enddo
               endif
c ::: ::: now do the tangential deriv part
c ::: ::: Null
c ::: ::: interior part of west side
               edgloc = 0.
               i = lo(1)
c ::: ::: ::: Y
               do k=lo(3),hi(3)
                  do j=lo(2)+1,hi(2)-1
                     if( maskw(-1+i,-1+j,k).eq.0.and.maskw(-1+i,1+j,k).eq.0)then
                         outloc = -0.5
                         outder = (-U(-1+i,-1+j,k,n)+U(-1+i,1+j,k,n))/(2.d0*hy)
                     elseif( maskw(-1+i,j,k).eq.0.and.maskw(-1+i,1+j,k).eq.0)then
                         outloc = -0.5
                         outder = (-3*U(-1+i,j,k,n)+4*U(-1+i,1+j,k,n)-U(-1+i,2+j,k,n))/(2.d0
     &  *hy)
                     elseif( maskw(-1+i,-1+j,k).eq.0.and.maskw(-1+i,j,k).eq.0)then
                         outloc = -0.5
                         outder = (U(-1+i,-2+j,k,n)-4*U(-1+i,-1+j,k,n)+3*U(-1+i,j,k,n))/(2.d
     &  0*hy)
                     else
                        outloc = x(-1)
                        if( flagbc .eq. 1) then
                           outder = exttd(i-1,j,k,n,2)
                        else
                           outder = 0.
                        endif
                     endif
                     innloc = 0.5
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                     lambda = (edgloc-outloc)/(innloc-outloc)
                     trander(i-1,j,k,n,2) = lambda*innder+(1-lambda)*outder
                  enddo
c ::: ::: ::: now endpoints
                  j = lo(2)
                  if( maskw(-1+i,-1+j,k).eq.0.and.maskw(-1+i,1+j,k).eq.0)then
                     outloc = -0.5
                     outder = (-U(-1+i,-1+j,k,n)+U(-1+i,1+j,k,n))/(2.d0*hy)
                  elseif( maskw(-1+i,j,k).eq.0.and.maskw(-1+i,1+j,k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(-1+i,j,k,n)+4*U(-1+i,1+j,k,n)-U(-1+i,2+j,k,n))/(2.d0
     &  *hy)
                  elseif( maskw(-1+i,-1+j,k).eq.0.and.maskw(-1+i,j,k).eq.0)then
                     outloc = -0.5
                     outder = (U(-1+i,-2+j,k,n)-4*U(-1+i,-1+j,k,n)+3*U(-1+i,j,k,n))/(2.d
     &  0*hy)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i-1,j,k,n,2)
                     else
                        outder = 0.
                     endif
                  endif
                  if(masks(i,-1+j,k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,1+j,k,n)-U(i,2+j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i-1,j,k,n,2) = lambda*innder+(1-lambda)*outder

                  j = hi(2)
                  if( maskw(-1+i,-1+j,k).eq.0.and.maskw(-1+i,1+j,k).eq.0)then
                     outloc = -0.5
                     outder = (-U(-1+i,-1+j,k,n)+U(-1+i,1+j,k,n))/(2.d0*hy)
                  elseif( maskw(-1+i,j,k).eq.0.and.maskw(-1+i,1+j,k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(-1+i,j,k,n)+4*U(-1+i,1+j,k,n)-U(-1+i,2+j,k,n))/(2.d0
     &  *hy)
                  elseif( maskw(-1+i,-1+j,k).eq.0.and.maskw(-1+i,j,k).eq.0)then
                     outloc = -0.5
                     outder = (U(-1+i,-2+j,k,n)-4*U(-1+i,-1+j,k,n)+3*U(-1+i,j,k,n))/(2.d
     &  0*hy)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i-1,j,k,n,2)
                     else
                        outder = 0.
                     endif
                  endif
                  if(maskn(i,1+j,k).gt.0) then
                     innder = (U(i,-2+j,k,n)-4*U(i,-1+j,k,n)+3*U(i,j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i-1,j,k,n,2) = lambda*innder+(1-lambda)*outder
               enddo
c ::: ::: ::: Z
               do j=lo(2),hi(2)
                  do k=lo(3)+1,hi(3)-1
                     if( maskw(-1+i,j,-1+k).eq.0.and.maskw(-1+i,j,1+k).eq.0)then
                         outloc = -0.5
                         outder = (-U(-1+i,j,-1+k,n)+U(-1+i,j,1+k,n))/(2.d0*hz)
                     elseif( maskw(-1+i,j,k).eq.0.and.maskw(-1+i,j,1+k).eq.0)then
                         outloc = -0.5
                         outder = (-3*U(-1+i,j,k,n)+4*U(-1+i,j,1+k,n)-U(-1+i,j,2+k,n))/(2.d0
     &  *hz)
                     elseif( maskw(-1+i,j,-1+k).eq.0.and.maskw(-1+i,j,k).eq.0)then
                         outloc = -0.5
                         outder = (U(-1+i,j,-2+k,n)-4*U(-1+i,j,-1+k,n)+3*U(-1+i,j,k,n))/(2.d
     &  0*hz)
                     else
                        outloc = x(-1)
                        if( flagbc .eq. 1) then
                           outder = exttd(i-1,j,k,n,3)
                        else
                           outder = 0.
                        endif
                     endif
                     innloc = 0.5
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                     lambda = (edgloc-outloc)/(innloc-outloc)
                     trander(i-1,j,k,n,3) = lambda*innder+(1-lambda)*outder
                  enddo
c ::: ::: ::: now endpoints
                  k = lo(3)
                  if( maskw(-1+i,j,-1+k).eq.0.and.maskw(-1+i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-U(-1+i,j,-1+k,n)+U(-1+i,j,1+k,n))/(2.d0*hz)
                  elseif( maskw(-1+i,j,k).eq.0.and.maskw(-1+i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(-1+i,j,k,n)+4*U(-1+i,j,1+k,n)-U(-1+i,j,2+k,n))/(2.d0
     &  *hz)
                  elseif( maskw(-1+i,j,-1+k).eq.0.and.maskw(-1+i,j,k).eq.0)then
                     outloc = -0.5
                     outder = (U(-1+i,j,-2+k,n)-4*U(-1+i,j,-1+k,n)+3*U(-1+i,j,k,n))/(2.d
     &  0*hz)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i-1,j,k,n,3)
                     else
                        outder = 0.
                     endif
                  endif
                  if(maskb(i,j,-1+k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,j,1+k,n)-U(i,j,2+k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i-1,j,k,n,3) = lambda*innder+(1-lambda)*outder

                  k = hi(3)
                  if( maskw(-1+i,j,-1+k).eq.0.and.maskw(-1+i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-U(-1+i,j,-1+k,n)+U(-1+i,j,1+k,n))/(2.d0*hz)
                  elseif( maskw(-1+i,j,k).eq.0.and.maskw(-1+i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(-1+i,j,k,n)+4*U(-1+i,j,1+k,n)-U(-1+i,j,2+k,n))/(2.d0
     &  *hz)
                  elseif( maskw(-1+i,j,-1+k).eq.0.and.maskw(-1+i,j,k).eq.0)then
                     outloc = -0.5
                     outder = (U(-1+i,j,-2+k,n)-4*U(-1+i,j,-1+k,n)+3*U(-1+i,j,k,n))/(2.d
     &  0*hz)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i-1,j,k,n,3)
                     else
                        outder = 0.
                     endif
                  endif
                  if(maskt(i,j,1+k).gt.0) then
                     innder = (U(i,j,-2+k,n)-4*U(i,j,-1+k,n)+3*U(i,j,k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i-1,j,k,n,3) = lambda*innder+(1-lambda)*outder
               enddo
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     u(lo(1)-1, j, k, n) = cvmgt(
     $                   -u(lo(1),j,k,n),
     $                    u(lo(1)-1,j,k,n),
     $                    maskw(lo(1)-1,j,k) .gt. 0)
                     trander(lo(1)-1,j,k,n,2) = 0.
                     trander(lo(1)-1,j,k,n,3) = 0.
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        den(lo(1),j,k,n) = cvmgt(-1.0, 0.0,
     $                       maskw(lo(1)-1,j,k) .gt. 0)
                     enddo
                  enddo
               endif
            else
               print *,'UNKNOWN BC ON WEST FACE IN APPLYBC'
               stop
            endif
	 enddo
      endif
c
c     The east face of the grid
c
      if(cdir .eq. 3) then
         do n = 1,nc
            if(is_neumann(bct(n))) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     u(hi(1)+1,j,k,n) = cvmgt(
     $                    u(hi(1), j, k, n),
     $                    u(hi(1)+1, j, k, n),
     $                    maske(hi(1)+1,j,k) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        den(hi(1),j,k,n) = 1.0
                     enddo
                  enddo
               endif
c ::: ::: now do the tangential derivative part
               i = hi(1)
c ::: ::: ::: Y
               do k=lo(3),hi(3)
c ::: ::: ::: ::: interior part of west face
                  do j=lo(2)+1,hi(2)-1
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                     trander(i+1,j,k,n,2) = innder
                  enddo
                  j = lo(2)
                  if(masks(i,-1+j,k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,1+j,k,n)-U(i,2+j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif
                  trander(i+1,j,k,n,2) = innder
                  j = hi(2)
                  if(maskn(i,1+j,k).gt.0) then
                     innder = (U(i,-2+j,k,n)-4*U(i,-1+j,k,n)+3*U(i,j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif
                  trander(i+1,j,k,n,2) = innder
               enddo
c ::: ::: ::: Z
               do j=lo(2),hi(2)
                  do k=lo(3)+1,hi(3)-1
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                     trander(i+1,j,k,n,3) = innder
                  enddo
                  k = lo(3)
                  if(maskb(i,j,-1+k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,j,1+k,n)-U(i,j,2+k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  trander(i+1,j,k,n,3) = innder
                  k = hi(3)
                  if(maskt(i,j,1+k).gt.0) then
                     innder = (U(i,j,-2+k,n)-4*U(i,j,-1+k,n)+3*U(i,j,k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  trander(i+1,j,k,n,3) = innder
               enddo
            else if (is_dirichlet(bct(n))) then
               do m=0,lenx
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(1)
               xInt = - 0.5
               call polyInterpCoeff(xInt, x, lenx+2, coef)
#if USE_TMP
               do k = lo(3), hi(3)
                  if ( flagbc .eq. 1 ) then
                     do j = lo(2), hi(2)
                        tmp(j-lo(2)) = bcval(hi(1)+1,j,k,n)*coef(-1)
                     enddo
                  else
                     do j = lo(2), hi(2)
                        tmp(j-lo(2)) = 0.0
                     enddo
                  endif
                  do m = 0, lenx
                     do j = lo(2), hi(2)
                        tmp(j-lo(2)) = tmp(j-lo(2))+u(hi(1)-m,j,k,n)*coef(m)
                     enddo
                  enddo
                  do j = lo(2), hi(2)
                     u(hi(1)+1,j,k,n) = cvmgt(
     $                    tmp(j-lo(2)),
     $                    u(hi(1)+1,j,k,n),
     $                    maske(hi(1)+1,j,k) .gt. 0)
                  enddo
               enddo
#else
               if ( flagbc .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        u(hi(1)+1,j,k,n) = cvmgt(
     $                       bcval(hi(1)+1,j,k,n)*coef(-1),
     $                       u(hi(1)+1,j,k,n),
     $                       maske(hi(1)+1,j,k) .gt. 0)
                     enddo
                  enddo
               else
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        u(hi(1)+1,j,k,n) = cvmgt(
     $                       0.0,
     $                       u(hi(1)+1,j,k,n),
     $                       maske(hi(1)+1,j,k) .gt. 0)
                     enddo
                  enddo
               endif
               do k = lo(3), hi(3)
                  do m = 0, lenx
                     do j = lo(2), hi(2)
                        u(hi(1)+1,j,k,n) = cvmgt(
     $                       u(hi(1)+1,j,k,n)
     $                       + u(hi(1)-m,j,k,n)*coef(m),
     $                       u(hi(1)+1,j,k,n),
     $                       maske(hi(1)+1,j,k) .gt. 0)
                     enddo
                  enddo
               enddo
#endif
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        den(hi(1),j,k,n)   = cvmgt(coef(0), 0.0,
     $                       maske(hi(1)+1,j,k) .gt. 0)
                     enddo
                  enddo
               endif
c ::: ::: now do the tangential deriv part
c ::: ::: interior part of west side
               edgloc = 0.
               i = hi(1)
c ::: ::: ::: Y
               do k=lo(3),hi(3)
                  do j=lo(2)+1,hi(2)-1
                     if( maske(1+i,-1+j,k).eq.0.and.maske(1+i,1+j,k).eq.0)then
                         outloc = -0.5
                         outder = (-U(1+i,-1+j,k,n)+U(1+i,1+j,k,n))/(2.d0*hy)
                     elseif( maske(1+i,j,k).eq.0.and.maske(1+i,1+j,k).eq.0)then
                         outloc = -0.5
                         outder = (-3*U(1+i,j,k,n)+4*U(1+i,1+j,k,n)-U(1+i,2+j,k,n))/(2.d0*hy
     &  )
                     elseif( maske(1+i,-1+j,k).eq.0.and.maske(1+i,j,k).eq.0)then
                         outloc = -0.5
                         outder = (U(1+i,-2+j,k,n)-4*U(1+i,-1+j,k,n)+3*U(1+i,j,k,n))/(2.d0*h
     &  y)
                     else
                        outloc = x(-1)
                        if( flagbc .eq. 1) then
                           outder = exttd(i+1,j,k,n,2)
                        else
                           outder = 0.
                        endif
                     endif
                     innloc = 0.5
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                     lambda = (edgloc-outloc)/(innloc-outloc)
                     trander(i+1,j,k,n,2) = lambda*innder+(1-lambda)*outder
                  enddo
c ::: ::: ::: now endpoints
                  j = lo(2)
                  if( maske(1+i,-1+j,k).eq.0.and.maske(1+i,1+j,k).eq.0)then
                     outloc = -0.5
                     outder = (-U(1+i,-1+j,k,n)+U(1+i,1+j,k,n))/(2.d0*hy)
                  elseif( maske(1+i,j,k).eq.0.and.maske(1+i,1+j,k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(1+i,j,k,n)+4*U(1+i,1+j,k,n)-U(1+i,2+j,k,n))/(2.d0*hy
     &  )
                  elseif( maske(1+i,-1+j,k).eq.0.and.maske(1+i,j,k).eq.0)then
                     outloc = -0.5
                     outder = (U(1+i,-2+j,k,n)-4*U(1+i,-1+j,k,n)+3*U(1+i,j,k,n))/(2.d0*h
     &  y)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i+1,j,k,n,2)
                     else
                        outder = 0.
                     endif
                  endif
                  if(masks(i,-1+j,k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,1+j,k,n)-U(i,2+j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i+1,j,k,n,2) = lambda*innder+(1-lambda)*outder

                  j = hi(2)
                  if( maske(1+i,-1+j,k).eq.0.and.maske(1+i,1+j,k).eq.0)then
                     outloc = -0.5
                     outder = (-U(1+i,-1+j,k,n)+U(1+i,1+j,k,n))/(2.d0*hy)
                  elseif( maske(1+i,j,k).eq.0.and.maske(1+i,1+j,k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(1+i,j,k,n)+4*U(1+i,1+j,k,n)-U(1+i,2+j,k,n))/(2.d0*hy
     &  )
                  elseif( maske(1+i,-1+j,k).eq.0.and.maske(1+i,j,k).eq.0)then
                     outloc = -0.5
                     outder = (U(1+i,-2+j,k,n)-4*U(1+i,-1+j,k,n)+3*U(1+i,j,k,n))/(2.d0*h
     &  y)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i+1,j,k,n,2)
                     else
                        outder = 0.
                     endif
                  endif
                  if(maskn(i,1+j,k).gt.0) then
                     innder = (U(i,-2+j,k,n)-4*U(i,-1+j,k,n)+3*U(i,j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i+1,j,k,n,2) = lambda*innder+(1-lambda)*outder
               enddo
c ::: ::: ::: Z
               do j=lo(2),hi(2)
                  do k=lo(3)+1,hi(3)-1
                     if( maske(1+i,j,-1+k).eq.0.and.maske(1+i,j,1+k).eq.0)then
                         outloc = -0.5
                         outder = (-U(1+i,j,-1+k,n)+U(1+i,j,1+k,n))/(2.d0*hz)
                     elseif( maske(1+i,j,k).eq.0.and.maske(1+i,j,1+k).eq.0)then
                         outloc = -0.5
                         outder = (-3*U(1+i,j,k,n)+4*U(1+i,j,1+k,n)-U(1+i,j,2+k,n))/(2.d0*hz
     &  )
                     elseif( maske(1+i,j,-1+k).eq.0.and.maske(1+i,j,k).eq.0)then
                         outloc = -0.5
                         outder = (U(1+i,j,-2+k,n)-4*U(1+i,j,-1+k,n)+3*U(1+i,j,k,n))/(2.d0*h
     &  z)
                     else
                        outloc = x(-1)
                        if( flagbc .eq. 1) then
                           outder = exttd(i+1,j,k,n,3)
                        else
                           outder = 0.
                        endif
                     endif
                     innloc = 0.5
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                     lambda = (edgloc-outloc)/(innloc-outloc)
                     trander(i+1,j,k,n,3) = lambda*innder+(1-lambda)*outder
                  enddo
c ::: ::: ::: now endpoints
                  k = lo(3)
                  if( maske(1+i,j,-1+k).eq.0.and.maske(1+i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-U(1+i,j,-1+k,n)+U(1+i,j,1+k,n))/(2.d0*hz)
                  elseif( maske(1+i,j,k).eq.0.and.maske(1+i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(1+i,j,k,n)+4*U(1+i,j,1+k,n)-U(1+i,j,2+k,n))/(2.d0*hz
     &  )
                  elseif( maske(1+i,j,-1+k).eq.0.and.maske(1+i,j,k).eq.0)then
                     outloc = -0.5
                     outder = (U(1+i,j,-2+k,n)-4*U(1+i,j,-1+k,n)+3*U(1+i,j,k,n))/(2.d0*h
     &  z)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i+1,j,k,n,3)
                     else
                        outder = 0.
                     endif
                  endif
                  if(maskb(i,j,-1+k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,j,1+k,n)-U(i,j,2+k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i+1,j,k,n,3) = lambda*innder+(1-lambda)*outder

                  k = hi(3)
                  if( maske(1+i,j,-1+k).eq.0.and.maske(1+i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-U(1+i,j,-1+k,n)+U(1+i,j,1+k,n))/(2.d0*hz)
                  elseif( maske(1+i,j,k).eq.0.and.maske(1+i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(1+i,j,k,n)+4*U(1+i,j,1+k,n)-U(1+i,j,2+k,n))/(2.d0*hz
     &  )
                  elseif( maske(1+i,j,-1+k).eq.0.and.maske(1+i,j,k).eq.0)then
                     outloc = -0.5
                     outder = (U(1+i,j,-2+k,n)-4*U(1+i,j,-1+k,n)+3*U(1+i,j,k,n))/(2.d0*h
     &  z)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i+1,j,k,n,3)
                     else
                        outder = 0.
                     endif
                  endif
                  if(maskt(i,j,1+k).gt.0) then
                     innder = (U(i,j,-2+k,n)-4*U(i,j,-1+k,n)+3*U(i,j,k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i+1,j,k,n,3) = lambda*innder+(1-lambda)*outder
               enddo
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               do k = lo(3), hi(3)
                  do j = lo(2), hi(2)
                     u(hi(1)+1, j, k, n) = cvmgt(
     $                   -u(hi(1),j,k,n),
     $                    u(hi(1)+1,j,k,n),
     $                    maske(hi(1)+1,j,k) .gt. 0)
                     trander(hi(1)+1,j,k,n,2) = 0.
                     trander(hi(1)+1,j,k,n,3) = 0.
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do j = lo(2), hi(2)
                        den(hi(1),j,k,n) = cvmgt(-1.0, 0.0,
     $                       maske(hi(1)+1,j,k) .gt. 0)
                     enddo
                  enddo
               endif
            else
               print *,'UNKNOWN BC ON RIGHT FACE IN APPLYBC'
               stop
            endif
         enddo
      endif
c
c     The south of the Grid
c
      if(cdir .eq. 1) then
	 do n = 1,nc
            if(is_neumann(bct(n))) then
               do k = lo(3), hi(3)
                  do i = lo(1),hi(1)
                     u(i,lo(2)-1,k,n) = cvmgt(
     $                    u(i,lo(2),k,n),
     $                    u(i,lo(2)-1,k,n),
     $                    masks(i,lo(2)-1,k) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1),hi(1)
                        den(i,lo(2),k,n)   = 1.0
                     enddo
                  enddo
               endif
c ::: ::: now do the tangential derivative part
               j = lo(2)
c ::: ::: ::: X
               do k=lo(3),hi(3)
c ::: ::: ::: ::: interior part of south face
                  do i=lo(1)+1,hi(1)-1
                     innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                     trander(i,j-1,k,n,1) = innder
                  enddo
                  i = lo(1)
                  if(maskw(-1+i,j,k).gt.0) then
                       innder = (-3*U(i,j,k,n)+4*U(1+i,j,k,n)-U(2+i,j,k,n))/(2.d0*hx)
                  else
                       innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif
                  trander(i,j-1,k,n,1) = innder
                  i = hi(1)
                  if(maske(1+i,j,k).gt.0) then
                     innder = (U(-2+i,j,k,n)-4*U(-1+i,j,k,n)+3*U(i,j,k,n))/(2.d0*hx)
                  else
                     innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif
                  trander(i,j-1,k,n,1) = innder
               enddo
c ::: ::: ::: Z
               do i=lo(1),hi(1)
c ::: ::: ::: ::: interior part of south face
                  do k=lo(3)+1,hi(3)-1
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                     trander(i,j-1,k,n,3) = innder
                  enddo
                  k = lo(3)
                  if(maskb(i,j,-1+k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,j,1+k,n)-U(i,j,2+k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  trander(i,j-1,k,n,3) = innder
                  k = hi(3)
                  if(maskt(i,j,1+k).gt.0) then
                     innder = (U(i,j,-2+k,n)-4*U(i,j,-1+k,n)+3*U(i,j,k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  trander(i,j-1,k,n,3) = innder
               enddo
            else if (is_dirichlet(bct(n))) then
               do m=0,leny
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(2)
               xInt = - 0.5
               call polyInterpCoeff(xInt, x, leny+2, coef)
#if USE_TMP
               do k = lo(3), hi(3)
                  if ( flagbc .eq. 1 ) then
                     do i = lo(1), hi(1)
                        tmp(i-lo(1)) = bcval(i,lo(2)-1,k,n)*coef(-1)
                     enddo
                  else
                     do i = lo(1), hi(1)
                        tmp(i-lo(1)) = 0.0
                     enddo
                  endif
                  do m = 0, leny
                     do i = lo(1), hi(1)
                        tmp(i-lo(1)) = tmp(i-lo(1))+u(i,lo(2)+m,k,n)*coef(m)
                     enddo
                  enddo
                  do i = lo(1), hi(1)
                     u(i,lo(2)-1,k,n) = cvmgt(
     $                    tmp(i-lo(1)),
     $                    u(i,lo(2)-1,k,n),
     $                    masks(i,lo(2)-1,k) .gt. 0)
                  enddo
               enddo
#else
               if ( flagbc .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        u(i,lo(2)-1,k,n) = cvmgt(
     $                       bcval(i,lo(2)-1,k,n)*coef(-1),
     $                       u(i,lo(2)-1,k,n),
     $                       masks(i,lo(2)-1,k) .gt. 0)
                     enddo
                  enddo
               else
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        u(i,lo(2)-1,k,n) = cvmgt(
     $                       0.0,
     $                       u(i,lo(2)-1,k,n),
     $                       masks(i,lo(2)-1,k) .gt. 0)
                     enddo
                  enddo
               endif
               do k = lo(3), hi(3)
                  do m = 0, leny
                     do i = lo(1), hi(1)
                        u(i, lo(2)-1, k, n) = cvmgt(
     $                       u(i, lo(2)-1,k,n)
     $                       + u(i, lo(2)+m, k,n)*coef(m),
     $                       u(i, lo(2)-1, k, n),
     $                       masks(i, lo(2)-1, k) .gt. 0)
                     enddo
                  enddo
               enddo
#endif
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        den(i, lo(2),k,n)   = cvmgt(coef(0), 0.0,
     $                       masks(i, lo(2)-1,k) .gt. 0)
                     enddo
                  enddo
               endif
c ::: ::: now do the tangential deriv part
c ::: ::: interior part first, followed by end points
               edgloc = 0.
               j=lo(2)
c ::: ::: ::: X
               do k=lo(3),hi(3)
                  do i=lo(1)+1, hi(1)-1
                     if( masks(-1+i,-1+j,k).eq.0.and.masks(1+i,-1+j,k).eq.0)then
                       outloc = -0.5
                       outder = (-U(-1+i,-1+j,k,n)+U(1+i,-1+j,k,n))/(2.d0*hx)
                    elseif( masks(i,-1+j,k).eq.0.and.masks(1+i,-1+j,k).eq.0)then
                       outloc = -0.5
                       outder = (-3*U(i,-1+j,k,n)+4*U(1+i,-1+j,k,n)-U(2+i,-1+j,k,n))/(2.d0
     &  *hx)
                    elseif( masks(-1+i,-1+j,k).eq.0.and.masks(i,-1+j,k).eq.0)then
                       outloc = -0.5
                       outder = (U(-2+i,-1+j,k,n)-4*U(-1+i,-1+j,k,n)+3*U(i,-1+j,k,n))/(2.d
     &  0*hx)
                    else
                       outloc = x(-1)
                       if( flagbc .eq. 1 ) then
                          outder = exttd(i,j-1,k,n,1)
                       else
                          outder = 0.
                       endif
                    endif
                    innloc = 0.5
                    innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                    lambda = (edgloc-outloc)/(innloc-outloc)
                    trander(i,j-1,k,n,1) = lambda*innder+(1-lambda)*outder
                  enddo

                  i = lo(1)
                  if( masks(-1+i,-1+j,k).eq.0.and.masks(1+i,-1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (-U(-1+i,-1+j,k,n)+U(1+i,-1+j,k,n))/(2.d0*hx)
                  elseif( masks(i,-1+j,k).eq.0.and.masks(1+i,-1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (-3*U(i,-1+j,k,n)+4*U(1+i,-1+j,k,n)-U(2+i,-1+j,k,n))/(2.d0
     &  *hx)
                  elseif( masks(-1+i,-1+j,k).eq.0.and.masks(i,-1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (U(-2+i,-1+j,k,n)-4*U(-1+i,-1+j,k,n)+3*U(i,-1+j,k,n))/(2.d
     &  0*hx)
                  else
                    outloc = x(-1)
                    if( flagbc .eq. 1 ) then
                       outder = exttd(i,j-1,k,n,1)
                    else
                       outder = 0.
                    endif
                  endif
                  if(maskw(-1+i,j,k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(1+i,j,k,n)-U(2+i,j,k,n))/(2.d0*hx)
                  else
                     innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j-1,k,n,1) = lambda*innder+(1-lambda)*outder

                  i = hi(1)
                  if( masks(-1+i,-1+j,k).eq.0.and.masks(1+i,-1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (-U(-1+i,-1+j,k,n)+U(1+i,-1+j,k,n))/(2.d0*hx)
                  elseif( masks(i,-1+j,k).eq.0.and.masks(1+i,-1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (-3*U(i,-1+j,k,n)+4*U(1+i,-1+j,k,n)-U(2+i,-1+j,k,n))/(2.d0
     &  *hx)
                  elseif( masks(-1+i,-1+j,k).eq.0.and.masks(i,-1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (U(-2+i,-1+j,k,n)-4*U(-1+i,-1+j,k,n)+3*U(i,-1+j,k,n))/(2.d
     &  0*hx)
                  else
                    outloc = x(-1)
                    if( flagbc .eq. 1 ) then
                       outder = exttd(i,j-1,k,n,1)
                    else
                       outder = 0.
                    endif
                  endif
                  if(maske(1+i,j,k).gt.0) then
                     innder = (U(-2+i,j,k,n)-4*U(-1+i,j,k,n)+3*U(i,j,k,n))/(2.d0*hx)
                  else
                     innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j-1,k,n,1) = lambda*innder+(1-lambda)*outder                 
               enddo
c ::: ::: ::: Z
               do i=lo(1), hi(1)
                  do k=lo(3)+1,hi(3)-1
                     if( masks(i,-1+j,-1+k).eq.0.and.masks(i,-1+j,1+k).eq.0)then
                       outloc = -0.5
                       outder = (-U(i,-1+j,-1+k,n)+U(i,-1+j,1+k,n))/(2.d0*hz)
                    elseif( masks(i,-1+j,k).eq.0.and.masks(i,-1+j,1+k).eq.0)then
                       outloc = -0.5
                       outder = (-3*U(i,-1+j,k,n)+4*U(i,-1+j,1+k,n)-U(i,-1+j,2+k,n))/(2.d0
     &  *hz)
                    elseif( masks(i,-1+j,-1+k).eq.0.and.masks(i,-1+j,k).eq.0)then
                       outloc = -0.5
                       outder = (U(i,-1+j,-2+k,n)-4*U(i,-1+j,-1+k,n)+3*U(i,-1+j,k,n))/(2.d
     &  0*hz)
                    else
                       outloc = x(-1)
                       if( flagbc .eq. 1 ) then
                          outder = exttd(i,j-1,k,n,3)
                       else
                          outder = 0.
                       endif
                    endif
                    innloc = 0.5
                    innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                    lambda = (edgloc-outloc)/(innloc-outloc)
                    trander(i,j-1,k,n,3) = lambda*innder+(1-lambda)*outder
                  enddo

                  k = lo(3)
                  if( masks(i,-1+j,-1+k).eq.0.and.masks(i,-1+j,1+k).eq.0)then
                    outloc = -0.5
                    outder = (-U(i,-1+j,-1+k,n)+U(i,-1+j,1+k,n))/(2.d0*hz)
                  elseif( masks(i,-1+j,k).eq.0.and.masks(i,-1+j,1+k).eq.0)then
                    outloc = -0.5
                    outder = (-3*U(i,-1+j,k,n)+4*U(i,-1+j,1+k,n)-U(i,-1+j,2+k,n))/(2.d0
     &  *hz)
                  elseif( masks(i,-1+j,-1+k).eq.0.and.masks(i,-1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (U(i,-1+j,-2+k,n)-4*U(i,-1+j,-1+k,n)+3*U(i,-1+j,k,n))/(2.d
     &  0*hz)
                  else
                    outloc = x(-1)
                    if( flagbc .eq. 1 ) then
                       outder = exttd(i,j-1,k,n,3)
                    else
                       outder = 0.
                    endif
                  endif
                  if(maskb(i,j,-1+k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,j,1+k,n)-U(i,j,2+k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j-1,k,n,3) = lambda*innder+(1-lambda)*outder

                  k = hi(3)
                  if( masks(i,-1+j,-1+k).eq.0.and.masks(i,-1+j,1+k).eq.0)then
                    outloc = -0.5
                    outder = (-U(i,-1+j,-1+k,n)+U(i,-1+j,1+k,n))/(2.d0*hz)
                  elseif( masks(i,-1+j,k).eq.0.and.masks(i,-1+j,1+k).eq.0)then
                    outloc = -0.5
                    outder = (-3*U(i,-1+j,k,n)+4*U(i,-1+j,1+k,n)-U(i,-1+j,2+k,n))/(2.d0
     &  *hz)
                  elseif( masks(i,-1+j,-1+k).eq.0.and.masks(i,-1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (U(i,-1+j,-2+k,n)-4*U(i,-1+j,-1+k,n)+3*U(i,-1+j,k,n))/(2.d
     &  0*hz)
                  else
                    outloc = x(-1)
                    if( flagbc .eq. 1 ) then
                       outder = exttd(i,j-1,k,n,3)
                    else
                       outder = 0.
                    endif
                  endif
                  if(maskt(i,j,1+k).gt.0) then
                     innder = (U(i,j,-2+k,n)-4*U(i,j,-1+k,n)+3*U(i,j,k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j-1,k,n,3) = lambda*innder+(1-lambda)*outder                 
               enddo
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     u(i, lo(2)-1, k, n) = cvmgt(
     $                   -u(i,lo(2),k,n),
     $                    u(i,lo(2)-1,k,n),
     $                    masks(i,lo(2)-1,k) .gt. 0)
                     trander(i,lo(2)-1,k,n,1) = 0.
                     trander(i,lo(2)-1,k,n,3) = 0.
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        den(i,lo(2),k,n) = cvmgt(-1.0, 0.0,
     $                       masks(i,lo(2)-1,k) .gt. 0)
                     enddo
                  enddo
               endif
            else
               print *,'UNKNOWN BC ON BOTTOM FACE IN APPLYBC'
               stop
            endif
         enddo
      endif
c
c     The north of the grid
c
      if (cdir .eq. 4) then
         do n = 1,nc
            if(is_neumann(bct(n))) then
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     u(i,hi(2)+1,k,n) = cvmgt(
     $                    u(i,hi(2),k,n),
     $                    u(i,hi(2)+1,k,n),
     $                    maskn(i,hi(2)+1,k) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        den(i,hi(2),k,n)   = 1.0
                     enddo
                  enddo
               endif
c ::: ::: now do the tangential derivative part
               j = hi(2)
c ::: ::: ::: X
               do k=lo(3),hi(3)
c ::: ::: ::: ::: interior part of south face
                  do i=lo(1)+1,hi(1)-1
                     innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                     trander(i,j+1,k,n,1) = innder
                  enddo
                  i = lo(1)
                  if(maskw(-1+i,j,k).gt.0) then
                       innder = (-3*U(i,j,k,n)+4*U(1+i,j,k,n)-U(2+i,j,k,n))/(2.d0*hx)
                  else
                       innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif
                  trander(i,j+1,k,n,1) = innder
                  i = hi(1)
                  if(maske(1+i,j,k).gt.0) then
                     innder = (U(-2+i,j,k,n)-4*U(-1+i,j,k,n)+3*U(i,j,k,n))/(2.d0*hx)
                  else
                     innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif
                  trander(i,j+1,k,n,1) = innder
               enddo
c ::: ::: ::: Z
               do i=lo(1),hi(1)
c ::: ::: ::: ::: interior part of south face
                  do k=lo(3)+1,hi(3)-1
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                     trander(i,j+1,k,n,3) = innder
                  enddo
                  k = lo(3)
                  if(maskb(i,j,-1+k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,j,1+k,n)-U(i,j,2+k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  trander(i,j+1,k,n,3) = innder
                  k = hi(3)
                  if(maskt(i,j,1+k).gt.0) then
                     innder = (U(i,j,-2+k,n)-4*U(i,j,-1+k,n)+3*U(i,j,k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  trander(i,j+1,k,n,3) = innder
               enddo
            else if (is_dirichlet(bct(n))) then
               do m=0,leny
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(2)
               xInt = - 0.5
               call polyInterpCoeff(xInt, x, leny+2, coef)
#if USE_TMP
               do k = lo(3), hi(3)
                  if ( flagbc .eq. 1 ) then
                     do i = lo(1), hi(1)
                        tmp(i-lo(1)) = bcval(i,hi(2)+1,k,n)*coef(-1)
                     enddo
                  else
                     do i = lo(1), hi(1)
                        tmp(i-lo(1)) = 0.0
                     enddo
                  endif
                  do m = 0, leny
                     do i = lo(1), hi(1)
                        tmp(i-lo(1)) = tmp(i-lo(1))+u(i,hi(2)-m,k,n)*coef(m)
                     enddo
                  enddo
                  do i = lo(1), hi(1)
                     u(i,hi(2)+1,k,n) = cvmgt(
     $                    tmp(i-lo(1)),
     $                    u(i,hi(2)+1,k,n),
     $                    maskn(i,hi(2)+1,k) .gt. 0)
                  enddo
               enddo
#else
               if ( flagbc .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        u(i,hi(2)+1,k,n) = cvmgt(
     $                       bcval(i,hi(2)+1,k,n)*coef(-1),
     $                       u(i,hi(2)+1,k,n),
     $                       maskn(i,hi(2)+1,k) .gt. 0)
                     enddo
                  enddo
               else
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        u(i,hi(2)+1,k,n) = cvmgt(
     $                       0.0,
     $                       u(i,hi(2)+1,k,n),
     $                       maskn(i,hi(2)+1,k) .gt. 0)
                     enddo
                  enddo
               endif
               do k = lo(3), hi(3)
                  do m = 0, leny
                     do i = lo(1), hi(1)
                        u(i, hi(2)+1,k,n) = cvmgt(
     $                       u(i,hi(2)+1,k,n)
     $                       + u(i, hi(2)-m,k,n)*coef(m),
     $                       u(i,hi(2)+1,k,n),
     $                       maskn(i,hi(2)+1,k) .gt. 0)
                     enddo
                  enddo
               enddo
#endif
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        den(i,hi(2),k,n)   = cvmgt(coef(0), 0.0,
     $                       maskn(i,hi(2)+1,k) .gt. 0)
                     enddo
                  enddo
               endif
c ::: ::: now do the tangential deriv part
c ::: ::: interior part first, followed by end points
               edgloc = 0.
               j=hi(2)
c ::: ::: ::: X
               do k=lo(3),hi(3)
                  do i=lo(1)+1, hi(1)-1
                     if( maskn(-1+i,1+j,k).eq.0.and.maskn(1+i,1+j,k).eq.0)then
                       outloc = -0.5
                       outder = (-U(-1+i,1+j,k,n)+U(1+i,1+j,k,n))/(2.d0*hx)
                    elseif( maskn(i,1+j,k).eq.0.and.maskn(1+i,1+j,k).eq.0)then
                       outloc = -0.5
                       outder = (-3*U(i,1+j,k,n)+4*U(1+i,1+j,k,n)-U(2+i,1+j,k,n))/(2.d0*hx
     &  )
                    elseif( maskn(-1+i,1+j,k).eq.0.and.maskn(i,1+j,k).eq.0)then
                       outloc = -0.5
                       outder = (U(-2+i,1+j,k,n)-4*U(-1+i,1+j,k,n)+3*U(i,1+j,k,n))/(2.d0*h
     &  x)
                    else
                       outloc = x(-1)
                       if( flagbc .eq. 1 ) then
                          outder = exttd(i,j+1,k,n,1)
                       else
                          outder = 0.
                       endif
                    endif
                    innloc = 0.5
                    innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                    lambda = (edgloc-outloc)/(innloc-outloc)
                    trander(i,j+1,k,n,1) = lambda*innder+(1-lambda)*outder
                  enddo

                  i = lo(1)
                  if( maskn(-1+i,1+j,k).eq.0.and.maskn(1+i,1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (-U(-1+i,1+j,k,n)+U(1+i,1+j,k,n))/(2.d0*hx)
                  elseif( maskn(i,1+j,k).eq.0.and.maskn(1+i,1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (-3*U(i,1+j,k,n)+4*U(1+i,1+j,k,n)-U(2+i,1+j,k,n))/(2.d0*hx
     &  )
                  elseif( maskn(-1+i,1+j,k).eq.0.and.maskn(i,1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (U(-2+i,1+j,k,n)-4*U(-1+i,1+j,k,n)+3*U(i,1+j,k,n))/(2.d0*h
     &  x)
                  else
                    outloc = x(-1)
                    if( flagbc .eq. 1 ) then
                       outder = exttd(i,j+1,k,n,1)
                    else
                       outder = 0.
                    endif
                  endif
                  if(maskw(-1+i,j,k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(1+i,j,k,n)-U(2+i,j,k,n))/(2.d0*hx)
                  else
                     innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j+1,k,n,1) = lambda*innder+(1-lambda)*outder

                  i = hi(1)
                  if( maskn(-1+i,1+j,k).eq.0.and.maskn(1+i,1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (-U(-1+i,1+j,k,n)+U(1+i,1+j,k,n))/(2.d0*hx)
                  elseif( maskn(i,1+j,k).eq.0.and.maskn(1+i,1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (-3*U(i,1+j,k,n)+4*U(1+i,1+j,k,n)-U(2+i,1+j,k,n))/(2.d0*hx
     &  )
                  elseif( maskn(-1+i,1+j,k).eq.0.and.maskn(i,1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (U(-2+i,1+j,k,n)-4*U(-1+i,1+j,k,n)+3*U(i,1+j,k,n))/(2.d0*h
     &  x)
                  else
                    outloc = x(-1)
                    if( flagbc .eq. 1 ) then
                       outder = exttd(i,j+1,k,n,1)
                    else
                       outder = 0.
                    endif
                  endif
                  if(maske(1+i,j,k).gt.0) then
                     innder = (U(-2+i,j,k,n)-4*U(-1+i,j,k,n)+3*U(i,j,k,n))/(2.d0*hx)
                  else
                     innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j+1,k,n,1) = lambda*innder+(1-lambda)*outder                 
               enddo
c ::: ::: ::: Z
               do i=lo(1), hi(1)
                  do k=lo(3)+1,hi(3)-1
                     if( maskn(i,1+j,-1+k).eq.0.and.maskn(i,1+j,1+k).eq.0)then
                       outloc = -0.5
                       outder = (-U(i,1+j,-1+k,n)+U(i,1+j,1+k,n))/(2.d0*hz)
                    elseif( maskn(i,1+j,k).eq.0.and.maskn(i,1+j,1+k).eq.0)then
                       outloc = -0.5
                       outder = (-3*U(i,1+j,k,n)+4*U(i,1+j,1+k,n)-U(i,1+j,2+k,n))/(2.d0*hz
     &  )
                    elseif( maskn(i,1+j,-1+k).eq.0.and.maskn(i,1+j,k).eq.0)then
                       outloc = -0.5
                       outder = (U(i,1+j,-2+k,n)-4*U(i,1+j,-1+k,n)+3*U(i,1+j,k,n))/(2.d0*h
     &  z)
                    else
                       outloc = x(-1)
                       if( flagbc .eq. 1 ) then
                          outder = exttd(i,j+1,k,n,3)
                       else
                          outder = 0.
                       endif
                    endif
                    innloc = 0.5
                    innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                    lambda = (edgloc-outloc)/(innloc-outloc)
                    trander(i,j+1,k,n,3) = lambda*innder+(1-lambda)*outder
                  enddo

                  k = lo(3)
                  if( maskn(i,1+j,-1+k).eq.0.and.maskn(i,1+j,1+k).eq.0)then
                    outloc = -0.5
                    outder = (-U(i,1+j,-1+k,n)+U(i,1+j,1+k,n))/(2.d0*hz)
                  elseif( maskn(i,1+j,k).eq.0.and.maskn(i,1+j,1+k).eq.0)then
                    outloc = -0.5
                    outder = (-3*U(i,1+j,k,n)+4*U(i,1+j,1+k,n)-U(i,1+j,2+k,n))/(2.d0*hz
     &  )
                  elseif( maskn(i,1+j,-1+k).eq.0.and.maskn(i,1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (U(i,1+j,-2+k,n)-4*U(i,1+j,-1+k,n)+3*U(i,1+j,k,n))/(2.d0*h
     &  z)
                  else
                    outloc = x(-1)
                    if( flagbc .eq. 1 ) then
                       outder = exttd(i,j+1,k,n,3)
                    else
                       outder = 0.
                    endif
                  endif
                  if(maskb(i,j,-1+k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,j,1+k,n)-U(i,j,2+k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j+1,k,n,3) = lambda*innder+(1-lambda)*outder

                  k = hi(3)
                  if( maskn(i,1+j,-1+k).eq.0.and.maskn(i,1+j,1+k).eq.0)then
                    outloc = -0.5
                    outder = (-U(i,1+j,-1+k,n)+U(i,1+j,1+k,n))/(2.d0*hz)
                  elseif( maskn(i,1+j,k).eq.0.and.maskn(i,1+j,1+k).eq.0)then
                    outloc = -0.5
                    outder = (-3*U(i,1+j,k,n)+4*U(i,1+j,1+k,n)-U(i,1+j,2+k,n))/(2.d0*hz
     &  )
                  elseif( maskn(i,1+j,-1+k).eq.0.and.maskn(i,1+j,k).eq.0)then
                    outloc = -0.5
                    outder = (U(i,1+j,-2+k,n)-4*U(i,1+j,-1+k,n)+3*U(i,1+j,k,n))/(2.d0*h
     &  z)
                  else
                    outloc = x(-1)
                    if( flagbc .eq. 1 ) then
                       outder = exttd(i,j+1,k,n,3)
                    else
                       outder = 0.
                    endif
                  endif
                  if(maskt(i,j,1+k).gt.0) then
                     innder = (U(i,j,-2+k,n)-4*U(i,j,-1+k,n)+3*U(i,j,k,n))/(2.d0*hz)
                  else
                     innder = (-U(i,j,-1+k,n)+U(i,j,1+k,n))/(2.d0*hz)
                  endif
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j+1,k,n,3) = lambda*innder+(1-lambda)*outder                 
               enddo
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               do k = lo(3), hi(3)
                  do i = lo(1), hi(1)
                     u(i, hi(2)+1, k, n) = cvmgt(
     $                   -u(i,hi(2),k,n),
     $                    u(i,hi(2)+1,k,n),
     $                    maskn(i,hi(2)+1,k) .gt. 0)
                     trander(i,hi(2)+1,k,n,1) = 0.
                     trander(i,hi(2)+1,k,n,3) = 0.
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do k = lo(3), hi(3)
                     do i = lo(1), hi(1)
                        den(i,hi(2),k,n) = cvmgt(-1.0, 0.0,
     $                       maskn(i,hi(2)+1,k) .gt. 0)
                     enddo
                  enddo
               endif
            else
               print *,'UNKNOWN BC ON TOP FACE IN APPLYBC'
               stop
            endif
         enddo
      endif
c
c     The bottom of the Grid
c
      if(cdir .eq. 2) then
         do n = 1,nc
           if(is_neumann(bct(n))) then
               do j = lo(2), hi(2)
                  do i = lo(1),hi(1)
                     u(i,j,lo(3)-1,n) = cvmgt(
     $                    u(i,j,lo(3),n),
     $                    u(i,j,lo(3)-1,n),
     $                    maskb(i,j,lo(3)-1) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1),hi(1)
                        den(i,j,lo(3),n)   = 1.0
                     enddo
                  enddo
               endif
c ::: ::: now do the tangential derivative part
               k = lo(3)
c ::: ::: X
               do j=lo(2),hi(2)
c ::: ::: ::: interior part of bottom face
                  do i=lo(1)+1,hi(1)-1
                     innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                     trander(i,j,k-1,n,1) = innder
                  enddo
                  i = lo(1)
                  if( maskw(-1+i,j,k).gt.0) then
                       innder = (-3*U(i,j,k,n)+4*U(1+i,j,k,n)-U(2+i,j,k,n))/(2.d0*hx)
                  else
                       innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif
                  trander(i,j,k-1,n,1) = innder
                  i = hi(1)
                  if( maske(1+i,j,k).gt.0) then
                       innder = (U(-2+i,j,k,n)-4*U(-1+i,j,k,n)+3*U(i,j,k,n))/(2.d0*hx)
                  else
                       innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif
                  trander(i,j,k-1,n,1) = innder
               enddo
c ::: ::: ::: Y
               do i=lo(1),hi(1)
c ::: ::: ::: interior part of bottom face
                  do j=lo(2)+1,hi(2)-1
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                     trander(i,j,k-1,n,2) = innder
                  enddo
                  j = lo(2)
                  if( masks(i,-1+j,k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,1+j,k,n)-U(i,2+j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif
                  trander(i,j,k-1,n,2) = innder
                  j = hi(2)
                  if( maskn(i,1+j,k).gt.0) then
                     innder = (U(i,-2+j,k,n)-4*U(i,-1+j,k,n)+3*U(i,j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif
                  trander(i,j,k-1,n,2) = innder
               enddo
            else if (is_dirichlet(bct(n))) then
               do m=0,lenz
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(3)
               xInt = - 0.5
               call polyInterpCoeff(xInt, x, lenz+2, coef)
#if USE_TMP
               do j = lo(2), hi(2)
                  if ( flagbc .eq. 1 ) then
                     do i = lo(1), hi(1)
                        tmp(i-lo(1)) = bcval(i,j,lo(3)-1,n)*coef(-1)
                     enddo
                  else
                     do i = lo(1), hi(1)
                        tmp(i-lo(1)) = 0.0
                     enddo
                  endif
                  do m = 0, lenz
                     do i = lo(1), hi(1)
                        tmp(i-lo(1)) = tmp(i-lo(1))+u(i,j,lo(3)+m,n)*coef(m)
                     enddo
                  enddo
                  do i=lo(1),hi(1)
                     u(i,j,lo(3)-1,n) = cvmgt(
     &                    tmp(i-lo(1)),
     &                    u(i,j,lo(3)-1,n),
     &                    maskb(i,j,lo(3)-1).gt. 0 )
                  enddo
               enddo
#else
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        u(i,j,lo(3)-1,n) = cvmgt(
     $                       bcval(i,j,lo(3)-1,n)*coef(-1),
     $                       u(i,j,lo(3)-1,n),
     $                       maskb(i,j,lo(3)-1) .gt. 0)
                     enddo
                  enddo
               else
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        u(i,j,lo(3)-1,n) = cvmgt(
     $                       0.0,
     $                       u(i,j,lo(3)-1,n),
     $                       maskb(i,j,lo(3)-1) .gt. 0)
                     enddo
                  enddo
               endif
               do j = lo(2), hi(2)
                  do m = 0, lenz
                     do i = lo(1), hi(1)
                        u(i, j, lo(3)-1, n) = cvmgt(
     $                       u(i, j, lo(3)-1,n)
     $                       + u(i, j, lo(3)+m, n)*coef(m),
     $                       u(i, j, lo(3)-1,n),
     $                       maskb(i, j, lo(3)-1) .gt. 0)
                     enddo
                  enddo
               enddo
#endif
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        den(i, j, lo(3),n)   = cvmgt(coef(0), 0.0,
     $                       maskb(i, j, lo(3)-1) .gt. 0)
                     enddo
                  enddo
               endif
c ::: ::: now do the tangential deriv part
c ::: ::: interior part first, followed by end points
               edgloc = 0.
               k = lo(3)
c ::: ::: ::: X
               do j=lo(2),hi(2)
                  do i=lo(1)+1, hi(1)-1
                     if( maskb(-1+i,j,-1+k).eq.0.and.maskb(1+i,j,-1+k).eq.0)then
                        outloc = -0.5
                        outder = (-U(-1+i,j,-1+k,n)+U(1+i,j,-1+k,n))/(2.d0*hx)
                     else if( maskb(i,j,-1+k).eq.0.and.maskb(1+i,j,-1+k).eq.0)then
                        outloc = -0.5
                        outder = (-3*U(i,j,-1+k,n)+4*U(1+i,j,-1+k,n)-U(2+i,j,-1+k,n))/(2.d0
     &  *hx)
                     else if( maskb(-1+i,j,-1+k).eq.0.and.maskb(i,j,-1+k).eq.0)then
                        outloc = -0.5
                        outder = (U(-2+i,j,-1+k,n)-4*U(-1+i,j,-1+k,n)+3*U(i,j,-1+k,n))/(2.d
     &  0*hx)
                     else
                        outloc = x(-1)
                        if( flagbc .eq. 1 ) then
                           outder = exttd(i,j,k-1,n,1)
                        else
                           outder = 0.
                        endif
                     endif
                     innloc = 0.5
                     innder = (-u(-1+i,j,k,n)+u(1+i,j,k,n))/(2.d0*hx)
                     lambda = (edgloc-outloc)/(innloc-outloc)
                     trander(i,j,k-1,n,1) = lambda*innder+(1-lambda)*outder
                  enddo

                  i = lo(1)
                  if( maskb(-1+i,j,-1+k).eq.0.and.maskb(1+i,j,-1+k).eq.0)then
                     outloc = -0.5
                     outder = (-U(-1+i,j,-1+k,n)+U(1+i,j,-1+k,n))/(2.d0*hx)
                  else if( maskb(i,j,-1+k).eq.0.and.maskb(1+i,j,-1+k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(i,j,-1+k,n)+4*U(1+i,j,-1+k,n)-U(2+i,j,-1+k,n))/(2.d0
     &  *hx)
                  else if( maskb(-1+i,j,-1+k).eq.0.and.maskb(i,j,-1+k).eq.0)then
                     outloc = -0.5
                     outder = (U(-2+i,j,-1+k,n)-4*U(-1+i,j,-1+k,n)+3*U(i,j,-1+k,n))/(2.d
     &  0*hx)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1 ) then
                        outder = exttd(i,j,k-1,n,1)
                     else
                        outder = 0.
                     endif
                  endif
                  if(maskw(-1+i,j,k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(1+i,j,k,n)-U(2+i,j,k,n))/(2.d0*hx)
                  else
                     innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif                     
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j,k-1,n,1) = lambda*innder+(1-lambda)*outder

                  i = hi(1)
                  if( maskb(-1+i,j,-1+k).eq.0.and.maskb(1+i,j,-1+k).eq.0)then
                     outloc = -0.5
                     outder = (-U(-1+i,j,-1+k,n)+U(1+i,j,-1+k,n))/(2.d0*hx)
                  else if( maskb(i,j,-1+k).eq.0.and.maskb(1+i,j,-1+k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(i,j,-1+k,n)+4*U(1+i,j,-1+k,n)-U(2+i,j,-1+k,n))/(2.d0
     &  *hx)
                  else if( maskb(-1+i,j,-1+k).eq.0.and.maskb(i,j,-1+k).eq.0)then
                     outloc = -0.5
                     outder = (U(-2+i,j,-1+k,n)-4*U(-1+i,j,-1+k,n)+3*U(i,j,-1+k,n))/(2.d
     &  0*hx)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1 ) then
                        outder = exttd(i,j,k-1,n,1)
                     else
                        outder = 0.
                     endif
                  endif
                  if(maske(1+i,j,k).gt.0) then
                     innder = (U(-2+i,j,k,n)-4*U(-1+i,j,k,n)+3*U(i,j,k,n))/(2.d0*hx)
                  else
                     innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif                     
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j,k-1,n,1) = lambda*innder+(1-lambda)*outder
               enddo
c ::: ::: ::: Y               
               do i=lo(1), hi(1)
                  do j=lo(2)+1,hi(2)-1
                     if( maskb(i,-1+j,-1+k).eq.0.and.maskb(i,1+j,-1+k).eq.0)then
                        outloc = -0.5
                        outder = (-U(i,-1+j,-1+k,n)+U(i,1+j,-1+k,n))/(2.d0*hy)
                     else if( maskb(i,j,-1+k).eq.0.and.maskb(i,1+j,-1+k).eq.0)then
                        outloc = -0.5
                        outder = (-3*U(i,j,-1+k,n)+4*U(i,1+j,-1+k,n)-U(i,2+j,-1+k,n))/(2.d0
     &  *hy)
                     else if( maskb(i,-1+j,-1+k).eq.0.and.maskb(i,j,-1+k).eq.0)then
                        outloc = -0.5
                        outder = (U(i,-2+j,-1+k,n)-4*U(i,-1+j,-1+k,n)+3*U(i,j,-1+k,n))/(2.d
     &  0*hy)
                     else
                        outloc = x(-1)
                        if( flagbc .eq. 1 ) then
                           outder = exttd(i,j,k-1,n,2)
                        else
                           outder = 0.
                        endif
                     endif
                     innloc = 0.5
                     innder = (-u(i,-1+j,k,n)+u(i,1+j,k,n))/(2.d0*hy)
                     lambda = (edgloc-outloc)/(innloc-outloc)
                     trander(i,j,k-1,n,2) = lambda*innder+(1-lambda)*outder
                  enddo

                  j = lo(2)
                  if( maskb(i,-1+j,-1+k).eq.0.and.maskb(i,1+j,-1+k).eq.0)then
                     outloc = -0.5
                     outder = (-U(i,-1+j,-1+k,n)+U(i,1+j,-1+k,n))/(2.d0*hy)
                  else if( maskb(i,j,-1+k).eq.0.and.maskb(i,1+j,-1+k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(i,j,-1+k,n)+4*U(i,1+j,-1+k,n)-U(i,2+j,-1+k,n))/(2.d0
     &  *hy)
                  else if( maskb(i,-1+j,-1+k).eq.0.and.maskb(i,j,-1+k).eq.0)then
                     outloc = -0.5
                     outder = (U(i,-2+j,-1+k,n)-4*U(i,-1+j,-1+k,n)+3*U(i,j,-1+k,n))/(2.d
     &  0*hy)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1 ) then
                        outder = exttd(i,j,k-1,n,2)
                     else
                        outder = 0.
                     endif
                  endif
                  if(masks(i,-1+j,k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,1+j,k,n)-U(i,2+j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif                     
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j,k-1,n,2) = lambda*innder+(1-lambda)*outder

                  j = hi(2)
                  if( maskb(i,-1+j,-1+k).eq.0.and.maskb(i,1+j,-1+k).eq.0)then
                     outloc = -0.5
                     outder = (-U(i,-1+j,-1+k,n)+U(i,1+j,-1+k,n))/(2.d0*hy)
                  else if( maskb(i,j,-1+k).eq.0.and.maskb(i,1+j,-1+k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(i,j,-1+k,n)+4*U(i,1+j,-1+k,n)-U(i,2+j,-1+k,n))/(2.d0
     &  *hy)
                  else if( maskb(i,-1+j,-1+k).eq.0.and.maskb(i,j,-1+k).eq.0)then
                     outloc = -0.5
                     outder = (U(i,-2+j,-1+k,n)-4*U(i,-1+j,-1+k,n)+3*U(i,j,-1+k,n))/(2.d
     &  0*hy)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1 ) then
                        outder = exttd(i,j,k-1,n,2)
                     else
                        outder = 0.
                     endif
                  endif
                  if(maskn(i,1+j,k).gt.0) then
                     innder = (U(i,-2+j,k,n)-4*U(i,-1+j,k,n)+3*U(i,j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif                     
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j,k-1,n,2) = lambda*innder+(1-lambda)*outder
               enddo
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     u(i, j, lo(3)-1, n) = cvmgt(
     $                   -u(i,j,lo(3),n),
     $                    u(i,j,lo(3)-1,n),
     $                    maskb(i,j,lo(3)-1) .gt. 0)
                     trander(i,j,lo(3)-1,n,1) = 0.
                     trander(i,j,lo(3)-1,n,2) = 0.
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        den(i,j,lo(3),n) = cvmgt(-1.0, 0.0,
     $                       maskb(i,j,lo(3)-1) .gt. 0)
                     enddo
                  enddo
               endif
            else
               print *,'UNKNOWN BC ON FRONT FACE IN APPLYBC'
               stop
            endif
         enddo
      endif
c
c     The top of the grid
c
      if (cdir .eq. 5) then
         do n = 1,nc
            if(is_neumann(bct(n))) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     u(i,j, hi(3)+1,n) = cvmgt(
     $                    u(i,j, hi(3),n),
     $                    u(i,j, hi(3)+1,n),
     $                    maskt(i,j, hi(3)+1) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        den(i,j, hi(3),n)   = 1.0
                     enddo
                  enddo
               endif
c ::: ::: now do the tangential derivative part
               k = hi(3)
c ::: ::: X
               do j=lo(2),hi(2)
c ::: ::: ::: interior part of bottom face
                  do i=lo(1)+1,hi(1)-1
                     innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                     trander(i,j,k+1,n,1) = innder
                  enddo
                  i = lo(1)
                  if( maskw(-1+i,j,k).gt.0) then
                       innder = (-3*U(i,j,k,n)+4*U(1+i,j,k,n)-U(2+i,j,k,n))/(2.d0*hx)
                  else
                       innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif
                  trander(i,j,k+1,n,1) = innder
                  i = hi(1)
                  if( maske(1+i,j,k).gt.0) then
                       innder = (U(-2+i,j,k,n)-4*U(-1+i,j,k,n)+3*U(i,j,k,n))/(2.d0*hx)
                  else
                       innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif
                  trander(i,j,k+1,n,1) = innder
               enddo
c ::: ::: ::: Y
               do i=lo(1),hi(1)
c ::: ::: ::: interior part of bottom face
                  do j=lo(2)+1,hi(2)-1
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                     trander(i,j,k+1,n,2) = innder
                  enddo
                  j = lo(2)
                  if( masks(i,-1+j,k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,1+j,k,n)-U(i,2+j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif
                  trander(i,j,k+1,n,2) = innder
                  j = hi(2)
                  if( maskn(i,1+j,k).gt.0) then
                     innder = (U(i,-2+j,k,n)-4*U(i,-1+j,k,n)+3*U(i,j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif
                  trander(i,j,k+1,n,2) = innder
               enddo
            else if (is_dirichlet(bct(n))) then
               do m=0,lenz
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(3)
               xInt = - 0.5
               call polyInterpCoeff(xInt, x, lenz+2, coef)
#if USE_TMP
               do j = lo(2), hi(2)
                  if ( flagbc .eq. 1 ) then
                     do i = lo(1), hi(1)
                        tmp(i-lo(1)) = bcval(i,j, hi(3)+1,n)*coef(-1)
                     enddo
                  else
                     do i = lo(1), hi(1)
                        tmp(i-lo(1)) = 0.0
                     enddo
                  endif
                  do m = 0, lenz
                     do i = lo(1), hi(1)
                        tmp(i-lo(1)) = tmp(i-lo(1))
     $                       + u(i, j, hi(3)-m,n)*coef(m)
                     enddo
                  enddo
                  do i = lo(1), hi(1)
                     u(i,j, hi(3)+1,n) = cvmgt(
     $                    tmp(i-lo(1)),
     $                    u(i,j, hi(3)+1,n),
     $                    maskt(i,j, hi(3)+1) .gt. 0)
                  enddo
               enddo
#else
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        u(i,j, hi(3)+1,n) = cvmgt(
     $                       bcval(i,j, hi(3)+1,n)*coef(-1),
     $                       u(i,j, hi(3)+1,n),
     $                       maskt(i,j, hi(3)+1) .gt. 0)
                     enddo
                  enddo
               else
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        u(i,j, hi(3)+1,n) = cvmgt(
     $                       0.0,
     $                       u(i,j, hi(3)+1,n),
     $                       maskt(i,j, hi(3)+1) .gt. 0)
                     enddo
                  enddo
               endif
               do j = lo(2), hi(2)
                  do m = 0, lenz
                     do i = lo(1), hi(1)
                        u(i, j, hi(3)+1,n) = cvmgt(
     $                       u(i,j, hi(3)+1,n)
     $                       + u(i, j, hi(3)-m,n)*coef(m),
     $                       u(i,j, hi(3)+1,n),
     $                       maskt(i,j, hi(3)+1) .gt. 0)
                     enddo
                  enddo
               enddo
#endif
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        den(i,j, hi(3),n)   = cvmgt(coef(0), 0.0,
     $                       maskt(i,j, hi(3)+1) .gt. 0)
                     enddo
                  enddo
               endif
c ::: ::: now do the tangential deriv part
c ::: ::: interior part first, followed by end points
               edgloc = 0.
               k = hi(3)
c ::: ::: ::: X
               do j=lo(2),hi(2)
                  do i=lo(1)+1, hi(1)-1
                     if( maskt(-1+i,j,1+k).eq.0.and.maskt(1+i,j,1+k).eq.0)then
                        outloc = -0.5
                        outder = (-U(-1+i,j,1+k,n)+U(1+i,j,1+k,n))/(2.d0*hx)
                     else if( maskt(i,j,1+k).eq.0.and.maskt(1+i,j,1+k).eq.0)then
                        outloc = -0.5
                        outder = (-3*U(i,j,1+k,n)+4*U(1+i,j,1+k,n)-U(2+i,j,1+k,n))/(2.d0*hx
     &  )
                     else if( maskt(-1+i,j,1+k).eq.0.and.maskt(i,j,1+k).eq.0)then
                        outloc = -0.5
                        outder = (U(-2+i,j,1+k,n)-4*U(-1+i,j,1+k,n)+3*U(i,j,1+k,n))/(2.d0*h
     &  x)
                     else
                        outloc = x(-1)
                        if( flagbc .eq. 1 ) then
                           outder = exttd(i,j,k+1,n,1)
                        else
                           outder = 0.
                        endif
                     endif
                     innloc = 0.5
                     innder = (-u(-1+i,j,k,n)+u(1+i,j,k,n))/(2.d0*hx)
                     lambda = (edgloc-outloc)/(innloc-outloc)
                     trander(i,j,k+1,n,1) = lambda*innder+(1-lambda)*outder
                  enddo

                  i = lo(1)
                  if( maskt(-1+i,j,1+k).eq.0.and.maskt(1+i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-U(-1+i,j,1+k,n)+U(1+i,j,1+k,n))/(2.d0*hx)
                  else if( maskt(i,j,1+k).eq.0.and.maskt(1+i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(i,j,1+k,n)+4*U(1+i,j,1+k,n)-U(2+i,j,1+k,n))/(2.d0*hx
     &  )
                  else if( maskt(-1+i,j,1+k).eq.0.and.maskt(i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (U(-2+i,j,1+k,n)-4*U(-1+i,j,1+k,n)+3*U(i,j,1+k,n))/(2.d0*h
     &  x)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1 ) then
                        outder = exttd(i,j,k+1,n,1)
                     else
                        outder = 0.
                     endif
                  endif
                  if(maskw(-1+i,j,k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(1+i,j,k,n)-U(2+i,j,k,n))/(2.d0*hx)
                  else
                     innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif                     
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j,k+1,n,1) = lambda*innder+(1-lambda)*outder

                  i = hi(1)
                  if( maskt(-1+i,j,1+k).eq.0.and.maskt(1+i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-U(-1+i,j,1+k,n)+U(1+i,j,1+k,n))/(2.d0*hx)
                  else if( maskt(i,j,1+k).eq.0.and.maskt(1+i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(i,j,1+k,n)+4*U(1+i,j,1+k,n)-U(2+i,j,1+k,n))/(2.d0*hx
     &  )
                  else if( maskt(-1+i,j,1+k).eq.0.and.maskt(i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (U(-2+i,j,1+k,n)-4*U(-1+i,j,1+k,n)+3*U(i,j,1+k,n))/(2.d0*h
     &  x)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1 ) then
                        outder = exttd(i,j,k+1,n,1)
                     else
                        outder = 0.
                     endif
                  endif
                  if(maske(1+i,j,k).gt.0) then
                     innder = (U(-2+i,j,k,n)-4*U(-1+i,j,k,n)+3*U(i,j,k,n))/(2.d0*hx)
                  else
                     innder = (-U(-1+i,j,k,n)+U(1+i,j,k,n))/(2.d0*hx)
                  endif                     
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j,k+1,n,1) = lambda*innder+(1-lambda)*outder
               enddo
c ::: ::: ::: Y               
               do i=lo(1), hi(1)
                  do j=lo(2)+1,hi(2)-1
                     if( maskt(i,-1+j,1+k).eq.0.and.maskt(i,1+j,1+k).eq.0)then
                        outloc = -0.5
                        outder = (-U(i,-1+j,1+k,n)+U(i,1+j,1+k,n))/(2.d0*hy)
                     else if( maskt(i,j,1+k).eq.0.and.maskt(i,1+j,1+k).eq.0)then
                        outloc = -0.5
                        outder = (-3*U(i,j,1+k,n)+4*U(i,1+j,1+k,n)-U(i,2+j,1+k,n))/(2.d0*hy
     &  )
                     else if( maskt(i,-1+j,1+k).eq.0.and.maskt(i,j,1+k).eq.0)then
                        outloc = -0.5
                        outder = (U(i,-2+j,1+k,n)-4*U(i,-1+j,1+k,n)+3*U(i,j,1+k,n))/(2.d0*h
     &  y)
                     else
                        outloc = x(-1)
                        if( flagbc .eq. 1 ) then
                           outder = exttd(i,j,k+1,n,2)
                        else
                           outder = 0.
                        endif
                     endif
                     innloc = 0.5
                     innder = (-u(i,-1+j,k,n)+u(i,1+j,k,n))/(2.d0*hy)
                     lambda = (edgloc-outloc)/(innloc-outloc)
                     trander(i,j,k+1,n,2) = lambda*innder+(1-lambda)*outder
                  enddo

                  j = lo(2)
                  if( maskt(i,-1+j,1+k).eq.0.and.maskt(i,1+j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-U(i,-1+j,1+k,n)+U(i,1+j,1+k,n))/(2.d0*hy)
                  else if( maskt(i,j,1+k).eq.0.and.maskt(i,1+j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(i,j,1+k,n)+4*U(i,1+j,1+k,n)-U(i,2+j,1+k,n))/(2.d0*hy
     &  )
                  else if( maskt(i,-1+j,1+k).eq.0.and.maskt(i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (U(i,-2+j,1+k,n)-4*U(i,-1+j,1+k,n)+3*U(i,j,1+k,n))/(2.d0*h
     &  y)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1 ) then
                        outder = exttd(i,j,k+1,n,2)
                     else
                        outder = 0.
                     endif
                  endif
                  if(masks(i,-1+j,k).gt.0) then
                     innder = (-3*U(i,j,k,n)+4*U(i,1+j,k,n)-U(i,2+j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif                     
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j,k+1,n,2) = lambda*innder+(1-lambda)*outder

                  j = hi(2)
                  if( maskt(i,-1+j,1+k).eq.0.and.maskt(i,1+j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-U(i,-1+j,1+k,n)+U(i,1+j,1+k,n))/(2.d0*hy)
                  else if( maskt(i,j,1+k).eq.0.and.maskt(i,1+j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (-3*U(i,j,1+k,n)+4*U(i,1+j,1+k,n)-U(i,2+j,1+k,n))/(2.d0*hy
     &  )
                  else if( maskt(i,-1+j,1+k).eq.0.and.maskt(i,j,1+k).eq.0)then
                     outloc = -0.5
                     outder = (U(i,-2+j,1+k,n)-4*U(i,-1+j,1+k,n)+3*U(i,j,1+k,n))/(2.d0*h
     &  y)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1 ) then
                        outder = exttd(i,j,k+1,n,2)
                     else
                        outder = 0.
                     endif
                  endif
                  if(maskn(i,1+j,k).gt.0) then
                     innder = (U(i,-2+j,k,n)-4*U(i,-1+j,k,n)+3*U(i,j,k,n))/(2.d0*hy)
                  else
                     innder = (-U(i,-1+j,k,n)+U(i,1+j,k,n))/(2.d0*hy)
                  endif                     
                  innloc = 0.5
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  trander(i,j,k+1,n,2) = lambda*innder+(1-lambda)*outder
               enddo
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               do j = lo(2), hi(2)
                  do i = lo(1), hi(1)
                     u(i, j, hi(3)+1, n) = cvmgt(
     $                   -u(i,j,hi(3),n),
     $                    u(i,j,hi(3)+1,n),
     $                    maskt(i,j,hi(3)+1) .gt. 0)
                     trander(i,j,hi(3)+1,n,1) = 0.
                     trander(i,j,hi(3)+1,n,2) = 0.
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     do i = lo(1), hi(1)
                        den(i,j,hi(3),n) = cvmgt(-1.0, 0.0,
     $                       maskt(i,j,hi(3)+1) .gt. 0)
                     enddo
                  enddo
               endif
            else
               print *,'UNKNOWN BC ON BACK FACE IN APPLYBC'
               stop
            endif
	 enddo
      endif
      end

