#include <AMReX_REAL.H>

#include "AMReX_DivVis_F.H"
#include <AMReX_LO_BCTYPES.H>
#include "AMReX_ArrayLim.H"

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
c ::: 2) the normal derivatives are evaluated in the normal fashion.
c ::: 3) tangential derivatives which reach outside the rectangle DO have
c :::    to check the masks.  

      subroutine FORT_GSRB(
     $                         u, DIMS(u),
     $                         rhs, DIMS(rhs),
     $                         alpha, beta,
     $                         a, DIMS(a),
     $                         muX, DIMS(muX),
     $                         muY, DIMS(muY),
     $                         muZ, DIMS(muZ),
     $                         maskn,DIMS(maskn),
     $                         fn, DIMS(fn),
     $                         maske,DIMS(maske),
     $                         fe, DIMS(fe),
     $                         maskw,DIMS(maskw),
     $                         fw, DIMS(fw),
     $                         masks,DIMS(masks),
     $                         fs, DIMS(fs),
     $                         maskt,DIMS(maskt),
     $                         ft, DIMS(ft),
     $                         maskb,DIMS(maskb),
     $                         fb, DIMS(fb),
     $                         trandern,DIMS(trandern),
     $                         trandere,DIMS(trandere),
     $                         tranderw,DIMS(tranderw),
     $                         tranders,DIMS(tranders),
     $                         trandert,DIMS(trandert),
     $                         tranderb,DIMS(tranderb),
     $                         lo,hi,h,ncomp,phaseflag
     $                         )
      implicit none
      REAL_T alpha, beta
      integer DIMDEC(u)
      integer DIMDEC(rhs)
      integer DIMDEC(a)
      integer DIMDEC(muX)
      integer DIMDEC(muY)
      integer DIMDEC(muZ)
      integer DIMDEC(maskn)
      integer DIMDEC(fn)
      integer DIMDEC(maske)
      integer DIMDEC(fe)
      integer DIMDEC(maskw)
      integer DIMDEC(fw)
      integer DIMDEC(masks)
      integer DIMDEC(fs)
      integer DIMDEC(maskt)
      integer DIMDEC(ft)
      integer DIMDEC(maskb)
      integer DIMDEC(fb)
      integer DIMDEC(trandern)
      integer DIMDEC(trandere)
      integer DIMDEC(tranderw)
      integer DIMDEC(tranders)
      integer DIMDEC(trandert)
      integer DIMDEC(tranderb)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      integer ncomp
      integer phaseflag
      REAL_T h(BL_SPACEDIM)
      REAL_T u(DIMV(u),3)
      REAL_T rhs(DIMV(rhs),3)
      REAL_T a(DIMV(a))
      REAL_T muX(DIMV(muX))
      REAL_T muY(DIMV(muY))
      REAL_T muZ(DIMV(muZ))
      integer  maskn(DIMV(maskn))
      REAL_T fn(DIMV(fn),3)
      integer  maske(DIMV(maske))
      REAL_T fe(DIMV(fe),3)
      integer  maskw(DIMV(maskw))
      REAL_T fw(DIMV(fw),3)
      integer  masks(DIMV(masks))
      REAL_T fs(DIMV(fs),3)
      integer  maskt(DIMV(maskt))
      REAL_T ft(DIMV(ft),3)
      integer  maskb(DIMV(maskb))
      REAL_T fb(DIMV(fb),3)

c ::: for transverse derives, first 3 is for variable, second 3 is for
c ::: direction of derivative.  Obviously not all are used, but this is easy.
      REAL_T trandern(DIMV(trandern),3,3)
      REAL_T trandere(DIMV(trandere),3,3)
      REAL_T tranderw(DIMV(tranderw),3,3)
      REAL_T tranders(DIMV(tranders),3,3)
      REAL_T trandert(DIMV(trandert),3,3)
      REAL_T tranderb(DIMV(tranderb),3,3)

      integer i,j,k
      REAL_T hx,hy,hz,invhx,invhy,invhz,ihxsq,ihysq,ihzsq,ihxyz
      REAL_T dudxe,dudxw
      REAL_T dvdyn,dvdys
      REAL_T dudye,dudyw
      REAL_T dudyn,dudys
      REAL_T dvdxe,dvdxw
      REAL_T dvdxn,dvdxs
      REAL_T dudzb,dudzt
      REAL_T dwdzt,dwdzb
      REAL_T dwdxb,dwdxt
      REAL_T dvdzb,dvdzt
      REAL_T dwdyb,dwdyt
      REAL_T dudzw,dudze
      REAL_T dwdxw,dwdxe
      REAL_T dvdzs,dvdzn
      REAL_T dwdys,dwdyn
      integer modx,mody,modz
      integer istart,jstart,kstart
      REAL_T diagu,diagv,diagw
      REAL_T operu,operv,operw

      logical False, True
      parameter( False=.false.)
      parameter( True=.true. )

      hx = h(1)
      hy = h(2)
      hz = h(3)

      invhx = 1.0d0/hx
      invhy = 1.0d0/hy
      invhz = 1.0d0/hz

      ihxsq = 1.0d0/hx**2
      ihysq = 1.0d0/hy**2
      ihzsq = 1.0d0/hz**2

      ihxyz = 1.0d0/(hx*hy*hz)

      if( ncomp .ne. 3 ) then
	write(6,*)'FORT_GSRB: bad ncomp',ncomp
	stop
      endif
c ::: convert phase flag into modx and mody and modz

      select case (phaseflag)
      case (0)
         modx = 0
         mody = 0
         modz = 0
      case (1)
         modx = 1
         mody = 0
         modz = 0
      case (2)
         modx = 0
         mody = 1
         modz = 0
      case (3)
         modx = 1
         mody = 1
         modz = 0
      case (4)
         modx = 0
         mody = 0
         modz = 1
      case (5)
         modx = 1
         mody = 0
         modz = 1
      case (6)
         modx = 0
         mody = 1
         modz = 1
      case (7)
         modx = 1
         mody = 1
         modz = 1
      case default
         write(6,*)'FORT_GSRB: bad phaseflag', phaseflag
         stop
      end select


c ::: cases:
c     ::: 0) interior.  All usual derivative expressions
c     ::: 1-6) face but not on edge.  
c     ::: 7-18) edge but not on corner.
c     ::: 19-26) corners.

c ::: case 0
c ::: Null
c ::: Null
      istart = lo(1)+1
      if( mod(istart,2) .ne. modx ) istart = istart+1
      jstart = lo(2)+1
      if( mod(jstart,2) .ne. mody ) jstart = jstart+1
      kstart = lo(3)+1
      if( mod(kstart,2) .ne. modz ) kstart = kstart+1
      do k=kstart,hi(3)-1,2
         do j=jstart,hi(2)-1,2
            do i=istart,hi(1)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy
      
c     ::: get expressions for tangential derivatives
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: in interior, diagonal elements need no corrections

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)
            enddo
         enddo
      enddo

c ::: case 1, top face
c ::: Null
c ::: Null
c ::: Null
      k=hi(3)
      if( mod(k,2) .eq. modz ) then
         do j=jstart,hi(2)-1,2
            do i=istart,hi(1)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
        dwdxt=trandert(i,j,1+k,3,1)
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
        dwdyt=trandert(i,j,1+k,3,2)
      if(maskt(-1+i,j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dudzw=((U(-1+i,j,-2+k,1)-4*U(-1+i,j,-1+k,1)+3*U(-1
     &  +i,j,k,1))*(.5d0*invhz)+(U(i,j,-2+k,1)-4*U(i,j,-1+k,1
     &  )+3*U(i,j,k,1))*(.5d0*invhz))*0.5d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
      endif
      if(maskt(i,j,1+k).gt.0.or.maskt(1+i,j,1+k).gt.0) then 
        dudze=((U(i,j,-2+k,1)-4*U(i,j,-1+k,1)+3*U(i,j,k,1)
     &  )*(.5d0*invhz)+(U(1+i,j,-2+k,1)-4*U(1+i,j,-1+k,1)+3*U
     &  (1+i,j,k,1))*(.5d0*invhz))*0.5d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
      endif
      if(maskt(i,-1+j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dvdzs=((U(i,-1+j,-2+k,2)-4*U(i,-1+j,-1+k,2)+3*U(i,
     &  -1+j,k,2))*(.5d0*invhz)+(U(i,j,-2+k,2)-4*U(i,j,-1+k,2
     &  )+3*U(i,j,k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
      endif
      if(maskt(i,j,1+k).gt.0.or.maskt(i,1+j,1+k).gt.0) then 
        dvdzn=((U(i,j,-2+k,2)-4*U(i,j,-1+k,2)+3*U(i,j,k,2)
     &  )*(.5d0*invhz)+(U(i,1+j,-2+k,2)-4*U(i,1+j,-1+k,2)+3*U
     &  (i,1+j,k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)
      endif

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskt(i,j,k+1) .gt. 0 ) then
         diagu = diagu -ft(i,j,k,1)*(beta*muZ(i,j,1+k)*ihzsq)
         diagv = diagv -ft(i,j,k,2)*(beta*muZ(i,j,1+k)*ihzsq)
         diagw = diagw -ft(i,j,k,3)*2.d0*beta*muZ(i,j,1+k)*ihzsq
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)
            enddo
         enddo
      endif

c ::: case 2, bottom face
c ::: Null
c ::: Null
c ::: Null
      k=lo(3)
      if( mod(k,2).eq.modz) then
         do j=jstart,hi(2)-1,2
            do i=istart,hi(1)-1,2
c     :::  get expressions for normal derivatives
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy
               
c     ::: get expressions for tangential derivatives
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
        dwdxb=tranderb(i,j,-1+k,3,1)
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
        dwdyb=tranderb(i,j,-1+k,3,2)
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
      if(maskb(-1+i,j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dudzw=((-3*U(-1+i,j,k,1)+4*U(-1+i,j,1+k,1)-U(-1+i,
     &  j,2+k,1))*(.5d0*invhz)+(-3*U(i,j,k,1)+4*U(i,j,1+k,1)-
     &  U(i,j,2+k,1))*(.5d0*invhz))*0.5d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
      endif
      if(maskb(i,j,-1+k).gt.0.or.maskb(1+i,j,-1+k).gt.0) then 
        dudze=((-3*U(i,j,k,1)+4*U(i,j,1+k,1)-U(i,j,2+k,1))
     &  *(.5d0*invhz)+(-3*U(1+i,j,k,1)+4*U(1+i,j,1+k,1)-U(1+i
     &  ,j,2+k,1))*(.5d0*invhz))*0.5d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
      endif
      if(maskb(i,-1+j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dvdzs=((-3*U(i,-1+j,k,2)+4*U(i,-1+j,1+k,2)-U(i,-1+
     &  j,2+k,2))*(.5d0*invhz)+(-3*U(i,j,k,2)+4*U(i,j,1+k,2)-
     &  U(i,j,2+k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
      endif
      if(maskb(i,j,-1+k).gt.0.or.maskb(i,1+j,-1+k).gt.0) then 
        dvdzn=((-3*U(i,j,k,2)+4*U(i,j,1+k,2)-U(i,j,2+k,2))
     &  *(.5d0*invhz)+(-3*U(i,1+j,k,2)+4*U(i,1+j,1+k,2)-U(i,1
     &  +j,2+k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)
      endif

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskb(i,j,k-1) .gt. 0 ) then
         diagu = diagu -fb(i,j,k,1)*(beta*muZ(i,j,k)*ihzsq)
         diagv = diagv -fb(i,j,k,2)*(beta*muZ(i,j,k)*ihzsq)
         diagw = diagw -fb(i,j,k,3)*2.d0*beta*muZ(i,j,k)*ihzsq
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)
            enddo
         enddo
      endif

c ::: case 3, west face
c ::: Null
c ::: Null
c ::: Null
      i = lo(1)
      if( mod(i,2) .eq. modx) then
      do k=kstart,hi(3)-1,2
         do j=jstart,hi(2)-1,2
c     :::  get expressions for normal 
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
        dudyw=tranderw(-1+i,j,k,1,2)
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,1+j,k).gt.0) then 
        dvdxn=((-3*U(i,j,k,2)+4*U(1+i,j,k,2)-U(2+i,j,k,2))
     &  *(.5d0*invhx)+(-3*U(i,1+j,k,2)+4*U(1+i,1+j,k,2)-U(2+i
     &  ,1+j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
      endif
      if(maskw(-1+i,-1+j,k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dvdxs=((-3*U(i,-1+j,k,2)+4*U(1+i,-1+j,k,2)-U(2+i,-
     &  1+j,k,2))*(.5d0*invhx)+(-3*U(i,j,k,2)+4*U(1+i,j,k,2)-
     &  U(2+i,j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
      endif
      if(maskw(-1+i,j,-1+k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dwdxb=((-3*U(i,j,-1+k,3)+4*U(1+i,j,-1+k,3)-U(2+i,j
     &  ,-1+k,3))*(.5d0*invhx)+(-3*U(i,j,k,3)+4*U(1+i,j,k,3)-
     &  U(2+i,j,k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
      endif
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,j,1+k).gt.0) then 
        dwdxt=((-3*U(i,j,k,3)+4*U(1+i,j,k,3)-U(2+i,j,k,3))
     &  *(.5d0*invhx)+(-3*U(i,j,1+k,3)+4*U(1+i,j,1+k,3)-U(2+i
     &  ,j,1+k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
      endif
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
        dudzw=tranderw(-1+i,j,k,1,3)
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskw(i-1,j,k) .gt. 0 ) then
         diagu = diagu -fw(i,j,k,1)*2.d0*beta*muX(i,j,k)*ihxsq
         diagv = diagv -fw(i,j,k,2)*(beta*muX(i,j,k)*ihxsq)
         diagw = diagw -fw(i,j,k,3)*(beta*muX(i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)

            enddo
         enddo
      endif

c ::: case 4, east face
c ::: Null
c ::: Null
c ::: Null
      i = hi(1)
      if( mod(i,2).eq.modx) then
      do k=kstart,hi(3)-1,2
         do j=jstart,hi(2)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
        dudye=trandere(1+i,j,k,1,2)
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
      if(maske(1+i,j,k).gt.0.or.maske(1+i,1+j,k).gt.0) then 
        dvdxn=((U(-2+i,j,k,2)-4*U(-1+i,j,k,2)+3*U(i,j,k,2)
     &  )*(.5d0*invhx)+(U(-2+i,1+j,k,2)-4*U(-1+i,1+j,k,2)+3*U
     &  (i,1+j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
      endif
      if(maske(1+i,-1+j,k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dvdxs=((U(-2+i,-1+j,k,2)-4*U(-1+i,-1+j,k,2)+3*U(i,
     &  -1+j,k,2))*(.5d0*invhx)+(U(-2+i,j,k,2)-4*U(-1+i,j,k,2
     &  )+3*U(i,j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
      endif
      if(maske(1+i,j,-1+k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dwdxb=((U(-2+i,j,-1+k,3)-4*U(-1+i,j,-1+k,3)+3*U(i,
     &  j,-1+k,3))*(.5d0*invhx)+(U(-2+i,j,k,3)-4*U(-1+i,j,k,3
     &  )+3*U(i,j,k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
      endif
      if(maske(1+i,j,k).gt.0.or.maske(1+i,j,1+k).gt.0) then 
        dwdxt=((U(-2+i,j,k,3)-4*U(-1+i,j,k,3)+3*U(i,j,k,3)
     &  )*(.5d0*invhx)+(U(-2+i,j,1+k,3)-4*U(-1+i,j,1+k,3)+3*U
     &  (i,j,1+k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
      endif
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
        dudze=trandere(1+i,j,k,1,3)
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maske(i+1,j,k) .gt. 0 ) then
         diagu = diagu -fe(i,j,k,1)*2.d0*beta*muX(1+i,j,k)*ihxsq
         diagv = diagv -fe(i,j,k,2)*(beta*muX(1+i,j,k)*ihxsq)
         diagw = diagw -fe(i,j,k,3)*(beta*muX(1+i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)
         enddo
      enddo
      endif
	

c ::: case 5, north face
c ::: Null
c ::: Null
c ::: Null
      j = hi(2)
      if( mod(j,2).eq.mody) then
      do k=kstart,hi(3)-1,2
            do i=istart,hi(1)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy
               
c     ::: get expressions for tangential derivatives
      if(maskn(i,1+j,k).gt.0.or.maskn(1+i,1+j,k).gt.0) then 
        dudye=((U(i,-2+j,k,1)-4*U(i,-1+j,k,1)+3*U(i,j,k,1)
     &  )*(.5d0*invhy)+(U(1+i,-2+j,k,1)-4*U(1+i,-1+j,k,1)+3*U
     &  (1+i,j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
      endif
      if(maskn(-1+i,1+j,k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dudyw=((U(-1+i,-2+j,k,1)-4*U(-1+i,-1+j,k,1)+3*U(-1
     &  +i,j,k,1))*(.5d0*invhy)+(U(i,-2+j,k,1)-4*U(i,-1+j,k,1
     &  )+3*U(i,j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
      endif
        dvdxn=trandern(i,1+j,k,2,1)
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
      if(maskn(i,1+j,-1+k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dwdyb=((U(i,-2+j,-1+k,3)-4*U(i,-1+j,-1+k,3)+3*U(i,
     &  j,-1+k,3))*(.5d0*invhy)+(U(i,-2+j,k,3)-4*U(i,-1+j,k,3
     &  )+3*U(i,j,k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
      endif
      if(maskn(i,1+j,k).gt.0.or.maskn(i,1+j,1+k).gt.0) then 
        dwdyt=((U(i,-2+j,k,3)-4*U(i,-1+j,k,3)+3*U(i,j,k,3)
     &  )*(.5d0*invhy)+(U(i,-2+j,1+k,3)-4*U(i,-1+j,1+k,3)+3*U
     &  (i,j,1+k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
      endif
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
        dvdzn=trandern(i,1+j,k,2,3)

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskn(i,j+1,k) .gt. 0 ) then
         diagu = diagu -fn(i,j,k,1)*(beta*muY(i,1+j,k)*ihysq)
         diagv = diagv -fn(i,j,k,2)*2.d0*beta*muY(i,1+j,k)*ihysq
         diagw = diagw -fn(i,j,k,3)*(beta*muY(i,1+j,k)*ihysq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)

             enddo
          enddo
       endif



c ::: case 6, south face
c ::: Null
c ::: Null
c ::: Null
      j = lo(2)
      if( mod(j,2).eq.mody) then
      do k=kstart,hi(3)-1,2
            do i=istart,hi(1)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
      if(masks(i,-1+j,k).gt.0.or.masks(1+i,-1+j,k).gt.0) then 
        dudye=((-3*U(i,j,k,1)+4*U(i,1+j,k,1)-U(i,2+j,k,1))
     &  *(.5d0*invhy)+(-3*U(1+i,j,k,1)+4*U(1+i,1+j,k,1)-U(1+i
     &  ,2+j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
      endif
      if(masks(-1+i,-1+j,k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dudyw=((-3*U(-1+i,j,k,1)+4*U(-1+i,1+j,k,1)-U(-1+i,
     &  2+j,k,1))*(.5d0*invhy)+(-3*U(i,j,k,1)+4*U(i,1+j,k,1)-
     &  U(i,2+j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
      endif
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
      if(masks(i,-1+j,-1+k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dwdyb=((-3*U(i,j,-1+k,3)+4*U(i,1+j,-1+k,3)-U(i,2+j
     &  ,-1+k,3))*(.5d0*invhy)+(-3*U(i,j,k,3)+4*U(i,1+j,k,3)-
     &  U(i,2+j,k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
      endif
      if(masks(i,-1+j,k).gt.0.or.masks(i,-1+j,1+k).gt.0) then 
        dwdyt=((-3*U(i,j,k,3)+4*U(i,1+j,k,3)-U(i,2+j,k,3))
     &  *(.5d0*invhy)+(-3*U(i,j,1+k,3)+4*U(i,1+j,1+k,3)-U(i,2
     &  +j,1+k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
      endif
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
        dvdzs=tranders(i,-1+j,k,2,3)
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( masks(i,j-1,k) .gt. 0 ) then
         diagu = diagu -fs(i,j,k,1)*(beta*muY(i,j,k)*ihysq)
         diagv = diagv -fs(i,j,k,2)*2.d0*beta*muY(i,j,k)*ihysq
         diagw = diagw -fs(i,j,k,3)*(beta*muY(i,j,k)*ihysq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)


            enddo
         enddo
      endif



c ::: case 7, top-north edge 
c ::: Null
c ::: Null
c ::: Null
      k=hi(3)
         j=hi(2)
         if( mod(k,2).eq.modz .and. mod(j,2).eq.mody ) then
            do i=istart,hi(1)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
      if(maskn(i,1+j,k).gt.0.or.maskn(1+i,1+j,k).gt.0) then 
        dudye=((U(i,-2+j,k,1)-4*U(i,-1+j,k,1)+3*U(i,j,k,1)
     &  )*(.5d0*invhy)+(U(1+i,-2+j,k,1)-4*U(1+i,-1+j,k,1)+3*U
     &  (1+i,j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
      endif
      if(maskn(-1+i,1+j,k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dudyw=((U(-1+i,-2+j,k,1)-4*U(-1+i,-1+j,k,1)+3*U(-1
     &  +i,j,k,1))*(.5d0*invhy)+(U(i,-2+j,k,1)-4*U(i,-1+j,k,1
     &  )+3*U(i,j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
      endif
        dvdxn=trandern(i,1+j,k,2,1)
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
        dwdxt=trandert(i,j,1+k,3,1)
      if(maskn(i,1+j,-1+k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dwdyb=((U(i,-2+j,-1+k,3)-4*U(i,-1+j,-1+k,3)+3*U(i,
     &  j,-1+k,3))*(.5d0*invhy)+(U(i,-2+j,k,3)-4*U(i,-1+j,k,3
     &  )+3*U(i,j,k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
      endif
        dwdyt=trandert(i,j,1+k,3,2)
      if(maskt(-1+i,j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dudzw=((U(-1+i,j,-2+k,1)-4*U(-1+i,j,-1+k,1)+3*U(-1
     &  +i,j,k,1))*(.5d0*invhz)+(U(i,j,-2+k,1)-4*U(i,j,-1+k,1
     &  )+3*U(i,j,k,1))*(.5d0*invhz))*0.5d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
      endif
      if(maskt(i,j,1+k).gt.0.or.maskt(1+i,j,1+k).gt.0) then 
        dudze=((U(i,j,-2+k,1)-4*U(i,j,-1+k,1)+3*U(i,j,k,1)
     &  )*(.5d0*invhz)+(U(1+i,j,-2+k,1)-4*U(1+i,j,-1+k,1)+3*U
     &  (1+i,j,k,1))*(.5d0*invhz))*0.5d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
      endif
      if(maskt(i,-1+j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dvdzs=((U(i,-1+j,-2+k,2)-4*U(i,-1+j,-1+k,2)+3*U(i,
     &  -1+j,k,2))*(.5d0*invhz)+(U(i,j,-2+k,2)-4*U(i,j,-1+k,2
     &  )+3*U(i,j,k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskt(i,j,k+1) .gt. 0 ) then
         diagu = diagu -ft(i,j,k,1)*(beta*muZ(i,j,1+k)*ihzsq)
         diagv = diagv -ft(i,j,k,2)*(beta*muZ(i,j,1+k)*ihzsq)
         diagw = diagw -ft(i,j,k,3)*2.d0*beta*muZ(i,j,1+k)*ihzsq
      endif
      if( maskn(i,j+1,k) .gt. 0 ) then
         diagu = diagu -fn(i,j,k,1)*(beta*muY(i,1+j,k)*ihysq)
         diagv = diagv -fn(i,j,k,1)*2.d0*beta*muY(i,1+j,k)*ihysq
         diagw = diagw -fn(i,j,k,1)*(beta*muY(i,1+j,k)*ihysq)
      endif
c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)

            enddo
      endif

c ::: case 8, top-south edge 
c ::: Null
      k=hi(3)
         j=lo(2)
         if( mod(k,2).eq.modz .and. mod(j,2).eq.mody) then
            do i=istart,hi(1)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
      if(masks(i,-1+j,k).gt.0.or.masks(1+i,-1+j,k).gt.0) then 
        dudye=((-3*U(i,j,k,1)+4*U(i,1+j,k,1)-U(i,2+j,k,1))
     &  *(.5d0*invhy)+(-3*U(1+i,j,k,1)+4*U(1+i,1+j,k,1)-U(1+i
     &  ,2+j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
      endif
      if(masks(-1+i,-1+j,k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dudyw=((-3*U(-1+i,j,k,1)+4*U(-1+i,1+j,k,1)-U(-1+i,
     &  2+j,k,1))*(.5d0*invhy)+(-3*U(i,j,k,1)+4*U(i,1+j,k,1)-
     &  U(i,2+j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
      endif
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
        dwdxt=trandert(i,j,1+k,3,1)
      if(masks(i,-1+j,-1+k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dwdyb=((-3*U(i,j,-1+k,3)+4*U(i,1+j,-1+k,3)-U(i,2+j
     &  ,-1+k,3))*(.5d0*invhy)+(-3*U(i,j,k,3)+4*U(i,1+j,k,3)-
     &  U(i,2+j,k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
      endif
        dwdyt=trandert(i,j,1+k,3,2)
      if(maskt(-1+i,j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dudzw=((U(-1+i,j,-2+k,1)-4*U(-1+i,j,-1+k,1)+3*U(-1
     &  +i,j,k,1))*(.5d0*invhz)+(U(i,j,-2+k,1)-4*U(i,j,-1+k,1
     &  )+3*U(i,j,k,1))*(.5d0*invhz))*0.5d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
      endif
      if(maskt(i,j,1+k).gt.0.or.maskt(1+i,j,1+k).gt.0) then 
        dudze=((U(i,j,-2+k,1)-4*U(i,j,-1+k,1)+3*U(i,j,k,1)
     &  )*(.5d0*invhz)+(U(1+i,j,-2+k,1)-4*U(1+i,j,-1+k,1)+3*U
     &  (1+i,j,k,1))*(.5d0*invhz))*0.5d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
      endif
        dvdzs=tranders(i,-1+j,k,2,3)
      if(maskt(i,j,1+k).gt.0.or.maskt(i,1+j,1+k).gt.0) then 
        dvdzn=((U(i,j,-2+k,2)-4*U(i,j,-1+k,2)+3*U(i,j,k,2)
     &  )*(.5d0*invhz)+(U(i,1+j,-2+k,2)-4*U(i,1+j,-1+k,2)+3*U
     &  (i,1+j,k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)
      endif

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskt(i,j,k+1) .gt. 0 ) then
         diagu = diagu -ft(i,j,k,1)*(beta*muZ(i,j,1+k)*ihzsq)
         diagv = diagv -ft(i,j,k,2)*(beta*muZ(i,j,1+k)*ihzsq)
         diagw = diagw -ft(i,j,k,3)*2.d0*beta*muZ(i,j,1+k)*ihzsq
      endif
      if( masks(i,j-1,k) .gt. 0 ) then
         diagu = diagu -fs(i,j,k,1)*(beta*muY(i,j,k)*ihysq)
         diagv = diagv -fs(i,j,k,2)*2.d0*beta*muY(i,j,k)*ihysq
         diagw = diagw -fs(i,j,k,3)*(beta*muY(i,j,k)*ihysq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)

            enddo
      endif


c ::: case 9, top-west edge 
c ::: Null
      k=hi(3)
      i=lo(1)
      if( mod(k,2).eq.modz .and. mod(i,2).eq.modx ) then
         do j=jstart,hi(2)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
        dudyw=tranderw(-1+i,j,k,1,2)
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,1+j,k).gt.0) then 
        dvdxn=((-3*U(i,j,k,2)+4*U(1+i,j,k,2)-U(2+i,j,k,2))
     &  *(.5d0*invhx)+(-3*U(i,1+j,k,2)+4*U(1+i,1+j,k,2)-U(2+i
     &  ,1+j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
      endif
      if(maskw(-1+i,-1+j,k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dvdxs=((-3*U(i,-1+j,k,2)+4*U(1+i,-1+j,k,2)-U(2+i,-
     &  1+j,k,2))*(.5d0*invhx)+(-3*U(i,j,k,2)+4*U(1+i,j,k,2)-
     &  U(2+i,j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
      endif
      if(maskw(-1+i,j,-1+k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dwdxb=((-3*U(i,j,-1+k,3)+4*U(1+i,j,-1+k,3)-U(2+i,j
     &  ,-1+k,3))*(.5d0*invhx)+(-3*U(i,j,k,3)+4*U(1+i,j,k,3)-
     &  U(2+i,j,k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
      endif
        dwdxt=trandert(i,j,1+k,3,1)
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
        dwdyt=trandert(i,j,1+k,3,2)
        dudzw=tranderw(-1+i,j,k,1,3)
      if(maskt(i,j,1+k).gt.0.or.maskt(1+i,j,1+k).gt.0) then 
        dudze=((U(i,j,-2+k,1)-4*U(i,j,-1+k,1)+3*U(i,j,k,1)
     &  )*(.5d0*invhz)+(U(1+i,j,-2+k,1)-4*U(1+i,j,-1+k,1)+3*U
     &  (1+i,j,k,1))*(.5d0*invhz))*0.5d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
      endif
      if(maskt(i,-1+j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dvdzs=((U(i,-1+j,-2+k,2)-4*U(i,-1+j,-1+k,2)+3*U(i,
     &  -1+j,k,2))*(.5d0*invhz)+(U(i,j,-2+k,2)-4*U(i,j,-1+k,2
     &  )+3*U(i,j,k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
      endif
      if(maskt(i,j,1+k).gt.0.or.maskt(i,1+j,1+k).gt.0) then 
        dvdzn=((U(i,j,-2+k,2)-4*U(i,j,-1+k,2)+3*U(i,j,k,2)
     &  )*(.5d0*invhz)+(U(i,1+j,-2+k,2)-4*U(i,1+j,-1+k,2)+3*U
     &  (i,1+j,k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)
      endif

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskt(i,j,k+1) .gt. 0 ) then
         diagu = diagu -ft(i,j,k,1)*(beta*muZ(i,j,1+k)*ihzsq)
         diagv = diagv -ft(i,j,k,2)*(beta*muZ(i,j,1+k)*ihzsq)
         diagw = diagw -ft(i,j,k,3)*2.d0*beta*muZ(i,j,1+k)*ihzsq
      endif
      if( maskw(i-1,j,k) .gt. 0 ) then
         diagu = diagu -fw(i,j,k,1)*2.d0*beta*muX(i,j,k)*ihxsq
         diagv = diagv -fw(i,j,k,2)*(beta*muX(i,j,k)*ihxsq)
         diagw = diagw -fw(i,j,k,3)*(beta*muX(i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)
         enddo
      endif



c ::: case 10, top-east edge 
c ::: Null
      k=hi(3)
      i=hi(1)
      if( mod(k,2).eq.modz .and. mod(i,2).eq.modx) then
         do j=jstart,hi(2)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy
               
c     ::: get expressions for tangential derivatives
        dudye=trandere(1+i,j,k,1,2)
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
      if(maske(1+i,j,k).gt.0.or.maske(1+i,1+j,k).gt.0) then 
        dvdxn=((U(-2+i,j,k,2)-4*U(-1+i,j,k,2)+3*U(i,j,k,2)
     &  )*(.5d0*invhx)+(U(-2+i,1+j,k,2)-4*U(-1+i,1+j,k,2)+3*U
     &  (i,1+j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
      endif
      if(maske(1+i,-1+j,k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dvdxs=((U(-2+i,-1+j,k,2)-4*U(-1+i,-1+j,k,2)+3*U(i,
     &  -1+j,k,2))*(.5d0*invhx)+(U(-2+i,j,k,2)-4*U(-1+i,j,k,2
     &  )+3*U(i,j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
      endif
      if(maske(1+i,j,-1+k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dwdxb=((U(-2+i,j,-1+k,3)-4*U(-1+i,j,-1+k,3)+3*U(i,
     &  j,-1+k,3))*(.5d0*invhx)+(U(-2+i,j,k,3)-4*U(-1+i,j,k,3
     &  )+3*U(i,j,k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
      endif
        dwdxt=trandert(i,j,1+k,3,1)
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
        dwdyt=trandert(i,j,1+k,3,2)
      if(maskt(-1+i,j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dudzw=((U(-1+i,j,-2+k,1)-4*U(-1+i,j,-1+k,1)+3*U(-1
     &  +i,j,k,1))*(.5d0*invhz)+(U(i,j,-2+k,1)-4*U(i,j,-1+k,1
     &  )+3*U(i,j,k,1))*(.5d0*invhz))*0.5d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
      endif
        dudze=trandere(1+i,j,k,1,3)
      if(maskt(i,-1+j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dvdzs=((U(i,-1+j,-2+k,2)-4*U(i,-1+j,-1+k,2)+3*U(i,
     &  -1+j,k,2))*(.5d0*invhz)+(U(i,j,-2+k,2)-4*U(i,j,-1+k,2
     &  )+3*U(i,j,k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
      endif
      if(maskt(i,j,1+k).gt.0.or.maskt(i,1+j,1+k).gt.0) then 
        dvdzn=((U(i,j,-2+k,2)-4*U(i,j,-1+k,2)+3*U(i,j,k,2)
     &  )*(.5d0*invhz)+(U(i,1+j,-2+k,2)-4*U(i,1+j,-1+k,2)+3*U
     &  (i,1+j,k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)
      endif
        
c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskt(i,j,k+1) .gt. 0 ) then
         diagu = diagu -ft(i,j,k,1)*(beta*muZ(i,j,1+k)*ihzsq)
         diagv = diagv -ft(i,j,k,2)*(beta*muZ(i,j,1+k)*ihzsq)
         diagw = diagw -ft(i,j,k,3)*2.d0*beta*muZ(i,j,1+k)*ihzsq
      endif
      if( maske(i+1,j,k) .gt. 0 ) then
         diagu = diagu -fe(i,j,k,1)*2.d0*beta*muX(1+i,j,k)*ihxsq
         diagv = diagv -fe(i,j,k,2)*(beta*muX(1+i,j,k)*ihxsq)
         diagw = diagw -fe(i,j,k,3)*(beta*muX(1+i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)
         enddo
      endif



c ::: case 11, bottom-north edge 
c ::: Null
      k=lo(3)
      j=hi(2)
      if( mod(k,2).eq.modz .and. mod(j,2).eq.mody) then
            do i=istart,hi(1)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
      if(maskn(i,1+j,k).gt.0.or.maskn(1+i,1+j,k).gt.0) then 
        dudye=((U(i,-2+j,k,1)-4*U(i,-1+j,k,1)+3*U(i,j,k,1)
     &  )*(.5d0*invhy)+(U(1+i,-2+j,k,1)-4*U(1+i,-1+j,k,1)+3*U
     &  (1+i,j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
      endif
      if(maskn(-1+i,1+j,k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dudyw=((U(-1+i,-2+j,k,1)-4*U(-1+i,-1+j,k,1)+3*U(-1
     &  +i,j,k,1))*(.5d0*invhy)+(U(i,-2+j,k,1)-4*U(i,-1+j,k,1
     &  )+3*U(i,j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
      endif
        dvdxn=trandern(i,1+j,k,2,1)
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
        dwdxb=tranderb(i,j,-1+k,3,1)
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
        dwdyb=tranderb(i,j,-1+k,3,2)
      if(maskn(i,1+j,k).gt.0.or.maskn(i,1+j,1+k).gt.0) then 
        dwdyt=((U(i,-2+j,k,3)-4*U(i,-1+j,k,3)+3*U(i,j,k,3)
     &  )*(.5d0*invhy)+(U(i,-2+j,1+k,3)-4*U(i,-1+j,1+k,3)+3*U
     &  (i,j,1+k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
      endif
      if(maskb(-1+i,j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dudzw=((-3*U(-1+i,j,k,1)+4*U(-1+i,j,1+k,1)-U(-1+i,
     &  j,2+k,1))*(.5d0*invhz)+(-3*U(i,j,k,1)+4*U(i,j,1+k,1)-
     &  U(i,j,2+k,1))*(.5d0*invhz))*0.5d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
      endif
      if(maskb(i,j,-1+k).gt.0.or.maskb(1+i,j,-1+k).gt.0) then 
        dudze=((-3*U(i,j,k,1)+4*U(i,j,1+k,1)-U(i,j,2+k,1))
     &  *(.5d0*invhz)+(-3*U(1+i,j,k,1)+4*U(1+i,j,1+k,1)-U(1+i
     &  ,j,2+k,1))*(.5d0*invhz))*0.5d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
      endif
      if(maskb(i,-1+j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dvdzs=((-3*U(i,-1+j,k,2)+4*U(i,-1+j,1+k,2)-U(i,-1+
     &  j,2+k,2))*(.5d0*invhz)+(-3*U(i,j,k,2)+4*U(i,j,1+k,2)-
     &  U(i,j,2+k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
      endif
        dvdzn=trandern(i,1+j,k,2,3)
        
c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskb(i,j,k-1) .gt. 0 ) then
         diagu = diagu -fb(i,j,k,1)*(beta*muZ(i,j,k)*ihzsq)
         diagv = diagv -fb(i,j,k,2)*(beta*muZ(i,j,k)*ihzsq)
         diagw = diagw -fb(i,j,k,3)*2.d0*beta*muZ(i,j,k)*ihzsq
      endif
      if( maskn(i,j+1,k) .gt. 0 ) then
         diagu = diagu -fn(i,j,k,1)*(beta*muY(i,1+j,k)*ihysq)
         diagv = diagv -fn(i,j,k,2)*2.d0*beta*muY(i,1+j,k)*ihysq
         diagw = diagw -fn(i,j,k,3)*(beta*muY(i,1+j,k)*ihysq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)
            enddo
        endif



c ::: case 12, bottom-south edge 
c ::: Null
      k=lo(3)
      j=lo(2)
      if( mod(k,2).eq.modz .and. mod(j,2).eq.mody ) then
            do i=istart,hi(1)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
      if(masks(i,-1+j,k).gt.0.or.masks(1+i,-1+j,k).gt.0) then 
        dudye=((-3*U(i,j,k,1)+4*U(i,1+j,k,1)-U(i,2+j,k,1))
     &  *(.5d0*invhy)+(-3*U(1+i,j,k,1)+4*U(1+i,1+j,k,1)-U(1+i
     &  ,2+j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
      endif
      if(masks(-1+i,-1+j,k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dudyw=((-3*U(-1+i,j,k,1)+4*U(-1+i,1+j,k,1)-U(-1+i,
     &  2+j,k,1))*(.5d0*invhy)+(-3*U(i,j,k,1)+4*U(i,1+j,k,1)-
     &  U(i,2+j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
      endif
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxb=tranderb(i,j,-1+k,3,1)
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
        dwdyb=tranderb(i,j,-1+k,3,2)
      if(masks(i,-1+j,k).gt.0.or.masks(i,-1+j,1+k).gt.0) then 
        dwdyt=((-3*U(i,j,k,3)+4*U(i,1+j,k,3)-U(i,2+j,k,3))
     &  *(.5d0*invhy)+(-3*U(i,j,1+k,3)+4*U(i,1+j,1+k,3)-U(i,2
     &  +j,1+k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
      endif
      if(maskb(-1+i,j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dudzw=((-3*U(-1+i,j,k,1)+4*U(-1+i,j,1+k,1)-U(-1+i,
     &  j,2+k,1))*(.5d0*invhz)+(-3*U(i,j,k,1)+4*U(i,j,1+k,1)-
     &  U(i,j,2+k,1))*(.5d0*invhz))*0.5d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
      endif
      if(maskb(i,j,-1+k).gt.0.or.maskb(1+i,j,-1+k).gt.0) then 
        dudze=((-3*U(i,j,k,1)+4*U(i,j,1+k,1)-U(i,j,2+k,1))
     &  *(.5d0*invhz)+(-3*U(1+i,j,k,1)+4*U(1+i,j,1+k,1)-U(1+i
     &  ,j,2+k,1))*(.5d0*invhz))*0.5d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
      endif
        dvdzs=tranders(i,-1+j,k,2,3)
      if(maskb(i,j,-1+k).gt.0.or.maskb(i,1+j,-1+k).gt.0) then 
        dvdzn=((-3*U(i,j,k,2)+4*U(i,j,1+k,2)-U(i,j,2+k,2))
     &  *(.5d0*invhz)+(-3*U(i,1+j,k,2)+4*U(i,1+j,1+k,2)-U(i,1
     &  +j,2+k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)
      endif

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskb(i,j,k-1) .gt. 0 ) then
         diagu = diagu -fb(i,j,k,1)*(beta*muZ(i,j,k)*ihzsq)
         diagv = diagv -fb(i,j,k,2)*(beta*muZ(i,j,k)*ihzsq)
         diagw = diagw -fb(i,j,k,3)*2.d0*beta*muZ(i,j,k)*ihzsq
      endif
      if( masks(i,j-1,k) .gt. 0 ) then
         diagu = diagu -fs(i,j,k,1)*(beta*muY(i,j,k)*ihysq)
         diagv = diagv -fs(i,j,k,2)*2.d0*beta*muY(i,j,k)*ihysq
         diagw = diagw -fs(i,j,k,3)*(beta*muY(i,j,k)*ihysq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)
            enddo
       endif



c ::: case 13, bottom-west edge 
c ::: Null
      k=lo(3)
      i=lo(1)
      if( mod(k,2).eq.modz .and. mod(i,2).eq.modx) then
      do j=jstart,hi(2)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
        dudyw=tranderw(-1+i,j,k,1,2)
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,1+j,k).gt.0) then 
        dvdxn=((-3*U(i,j,k,2)+4*U(1+i,j,k,2)-U(2+i,j,k,2))
     &  *(.5d0*invhx)+(-3*U(i,1+j,k,2)+4*U(1+i,1+j,k,2)-U(2+i
     &  ,1+j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
      endif
      if(maskw(-1+i,-1+j,k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dvdxs=((-3*U(i,-1+j,k,2)+4*U(1+i,-1+j,k,2)-U(2+i,-
     &  1+j,k,2))*(.5d0*invhx)+(-3*U(i,j,k,2)+4*U(1+i,j,k,2)-
     &  U(2+i,j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
      endif
        dwdxb=tranderb(i,j,-1+k,3,1)
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,j,1+k).gt.0) then 
        dwdxt=((-3*U(i,j,k,3)+4*U(1+i,j,k,3)-U(2+i,j,k,3))
     &  *(.5d0*invhx)+(-3*U(i,j,1+k,3)+4*U(1+i,j,1+k,3)-U(2+i
     &  ,j,1+k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
      endif
        dwdyb=tranderb(i,j,-1+k,3,2)
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
        dudzw=tranderw(-1+i,j,k,1,3)
      if(maskb(i,j,-1+k).gt.0.or.maskb(1+i,j,-1+k).gt.0) then 
        dudze=((-3*U(i,j,k,1)+4*U(i,j,1+k,1)-U(i,j,2+k,1))
     &  *(.5d0*invhz)+(-3*U(1+i,j,k,1)+4*U(1+i,j,1+k,1)-U(1+i
     &  ,j,2+k,1))*(.5d0*invhz))*0.5d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
      endif
      if(maskb(i,-1+j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dvdzs=((-3*U(i,-1+j,k,2)+4*U(i,-1+j,1+k,2)-U(i,-1+
     &  j,2+k,2))*(.5d0*invhz)+(-3*U(i,j,k,2)+4*U(i,j,1+k,2)-
     &  U(i,j,2+k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
      endif
      if(maskb(i,j,-1+k).gt.0.or.maskb(i,1+j,-1+k).gt.0) then 
        dvdzn=((-3*U(i,j,k,2)+4*U(i,j,1+k,2)-U(i,j,2+k,2))
     &  *(.5d0*invhz)+(-3*U(i,1+j,k,2)+4*U(i,1+j,1+k,2)-U(i,1
     &  +j,2+k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)
      endif

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskb(i,j,k-1) .gt. 0 ) then
         diagu = diagu -fb(i,j,k,1)*(beta*muZ(i,j,k)*ihzsq)
         diagv = diagv -fb(i,j,k,2)*(beta*muZ(i,j,k)*ihzsq)
         diagw = diagw -fb(i,j,k,3)*2.d0*beta*muZ(i,j,k)*ihzsq
      endif
      if( maskw(i-1,j,k) .gt. 0 ) then
         diagu = diagu -fw(i,j,k,1)*2.d0*beta*muX(i,j,k)*ihxsq
         diagv = diagv -fw(i,j,k,2)*(beta*muX(i,j,k)*ihxsq)
         diagw = diagw -fw(i,j,k,3)*(beta*muX(i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)
         enddo
         endif



c ::: case 14, bottom-east edge 
c ::: Null
      k=lo(3)
      i=hi(1)
      if( mod(k,2).eq.modz .and. mod(i,2).eq.modx ) then
         do j=jstart,hi(2)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
      if(maske(1+i,j,k).gt.0.or.maske(1+i,1+j,k).gt.0) then 
        dvdxn=((U(-2+i,j,k,2)-4*U(-1+i,j,k,2)+3*U(i,j,k,2)
     &  )*(.5d0*invhx)+(U(-2+i,1+j,k,2)-4*U(-1+i,1+j,k,2)+3*U
     &  (i,1+j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
      endif
      if(maske(1+i,-1+j,k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dvdxs=((U(-2+i,-1+j,k,2)-4*U(-1+i,-1+j,k,2)+3*U(i,
     &  -1+j,k,2))*(.5d0*invhx)+(U(-2+i,j,k,2)-4*U(-1+i,j,k,2
     &  )+3*U(i,j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
      endif
        dwdxb=tranderb(i,j,-1+k,3,1)
      if(maske(1+i,j,k).gt.0.or.maske(1+i,j,1+k).gt.0) then 
        dwdxt=((U(-2+i,j,k,3)-4*U(-1+i,j,k,3)+3*U(i,j,k,3)
     &  )*(.5d0*invhx)+(U(-2+i,j,1+k,3)-4*U(-1+i,j,1+k,3)+3*U
     &  (i,j,1+k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
      endif
        dudye=trandere(1+i,j,k,1,2)
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
        dwdyb=tranderb(i,j,-1+k,3,2)
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
      if(maskb(-1+i,j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dudzw=((-3*U(-1+i,j,k,1)+4*U(-1+i,j,1+k,1)-U(-1+i,
     &  j,2+k,1))*(.5d0*invhz)+(-3*U(i,j,k,1)+4*U(i,j,1+k,1)-
     &  U(i,j,2+k,1))*(.5d0*invhz))*0.5d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
      endif
        dudze=trandere(1+i,j,k,1,3)
      if(maskb(i,-1+j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dvdzs=((-3*U(i,-1+j,k,2)+4*U(i,-1+j,1+k,2)-U(i,-1+
     &  j,2+k,2))*(.5d0*invhz)+(-3*U(i,j,k,2)+4*U(i,j,1+k,2)-
     &  U(i,j,2+k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
      endif
      if(maskb(i,j,-1+k).gt.0.or.maskb(i,1+j,-1+k).gt.0) then 
        dvdzn=((-3*U(i,j,k,2)+4*U(i,j,1+k,2)-U(i,j,2+k,2))
     &  *(.5d0*invhz)+(-3*U(i,1+j,k,2)+4*U(i,1+j,1+k,2)-U(i,1
     &  +j,2+k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)
      endif

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskb(i,j,k-1) .gt. 0 ) then
         diagu = diagu -fb(i,j,k,1)*(beta*muZ(i,j,k)*ihzsq)
         diagv = diagv -fb(i,j,k,2)*(beta*muZ(i,j,k)*ihzsq)
         diagw = diagw -fb(i,j,k,3)*2.d0*beta*muZ(i,j,k)*ihzsq
      endif
      if( maske(i+1,j,k) .gt. 0 ) then
         diagu = diagu -fe(i,j,k,1)*2.d0*beta*muX(1+i,j,k)*ihxsq
         diagv = diagv -fe(i,j,k,2)*(beta*muX(1+i,j,k)*ihxsq)
         diagw = diagw -fe(i,j,k,3)*(beta*muX(1+i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)
         enddo
         endif



c ::: case 15, east-north edge
c ::: Null
      j=hi(2)
      i=hi(1)
      if( mod(j,2).eq.mody .and. mod(i,2).eq.modx ) then
      do k=kstart,hi(3)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if(maske(1+i,-1+j,k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dvdxs=((U(-2+i,-1+j,k,2)-4*U(-1+i,-1+j,k,2)+3*U(i,
     &  -1+j,k,2))*(.5d0*invhx)+(U(-2+i,j,k,2)-4*U(-1+i,j,k,2
     &  )+3*U(i,j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
      endif
      if(maske(1+i,j,-1+k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dwdxb=((U(-2+i,j,-1+k,3)-4*U(-1+i,j,-1+k,3)+3*U(i,
     &  j,-1+k,3))*(.5d0*invhx)+(U(-2+i,j,k,3)-4*U(-1+i,j,k,3
     &  )+3*U(i,j,k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
      endif
      if(maske(1+i,j,k).gt.0.or.maske(1+i,j,1+k).gt.0) then 
        dwdxt=((U(-2+i,j,k,3)-4*U(-1+i,j,k,3)+3*U(i,j,k,3)
     &  )*(.5d0*invhx)+(U(-2+i,j,1+k,3)-4*U(-1+i,j,1+k,3)+3*U
     &  (i,j,1+k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
      endif
        dudye=trandere(1+i,j,k,1,2)
      if(maskn(-1+i,1+j,k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dudyw=((U(-1+i,-2+j,k,1)-4*U(-1+i,-1+j,k,1)+3*U(-1
     &  +i,j,k,1))*(.5d0*invhy)+(U(i,-2+j,k,1)-4*U(i,-1+j,k,1
     &  )+3*U(i,j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
      endif
      if(maskn(i,1+j,-1+k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dwdyb=((U(i,-2+j,-1+k,3)-4*U(i,-1+j,-1+k,3)+3*U(i,
     &  j,-1+k,3))*(.5d0*invhy)+(U(i,-2+j,k,3)-4*U(i,-1+j,k,3
     &  )+3*U(i,j,k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
      endif
      if(maskn(i,1+j,k).gt.0.or.maskn(i,1+j,1+k).gt.0) then 
        dwdyt=((U(i,-2+j,k,3)-4*U(i,-1+j,k,3)+3*U(i,j,k,3)
     &  )*(.5d0*invhy)+(U(i,-2+j,1+k,3)-4*U(i,-1+j,1+k,3)+3*U
     &  (i,j,1+k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
      endif
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
        dudze=trandere(1+i,j,k,1,3)
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
        dvdzn=trandern(i,1+j,k,2,3)

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskn(i,j+1,k) .gt. 0 ) then
         diagu = diagu -fn(i,j,k,1)*(beta*muY(i,1+j,k)*ihysq)
         diagv = diagv -fn(i,j,k,2)*2.d0*beta*muY(i,1+j,k)*ihysq
         diagw = diagw -fn(i,j,k,3)*(beta*muY(i,1+j,k)*ihysq)
      endif
      if( maske(i+1,j,k) .gt. 0 ) then
         diagu = diagu -fe(i,j,k,1)*2.d0*beta*muX(1+i,j,k)*ihxsq
         diagv = diagv -fe(i,j,k,2)*(beta*muX(1+i,j,k)*ihxsq)
         diagw = diagw -fe(i,j,k,3)*(beta*muX(1+i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)

      enddo
      endif


c ::: case 16, east-south edge
c ::: Null
      j=lo(2)
      i=hi(1)
      if( mod(j,2).eq.mody .and. mod(i,2).eq.modx ) then
      do k=kstart,hi(3)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
      if(maske(1+i,j,k).gt.0.or.maske(1+i,1+j,k).gt.0) then 
        dvdxn=((U(-2+i,j,k,2)-4*U(-1+i,j,k,2)+3*U(i,j,k,2)
     &  )*(.5d0*invhx)+(U(-2+i,1+j,k,2)-4*U(-1+i,1+j,k,2)+3*U
     &  (i,1+j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
      if(maske(1+i,j,-1+k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dwdxb=((U(-2+i,j,-1+k,3)-4*U(-1+i,j,-1+k,3)+3*U(i,
     &  j,-1+k,3))*(.5d0*invhx)+(U(-2+i,j,k,3)-4*U(-1+i,j,k,3
     &  )+3*U(i,j,k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
      endif
      if(maske(1+i,j,k).gt.0.or.maske(1+i,j,1+k).gt.0) then 
        dwdxt=((U(-2+i,j,k,3)-4*U(-1+i,j,k,3)+3*U(i,j,k,3)
     &  )*(.5d0*invhx)+(U(-2+i,j,1+k,3)-4*U(-1+i,j,1+k,3)+3*U
     &  (i,j,1+k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
      endif
        dudye=trandere(1+i,j,k,1,2)
      if(masks(-1+i,-1+j,k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dudyw=((-3*U(-1+i,j,k,1)+4*U(-1+i,1+j,k,1)-U(-1+i,
     &  2+j,k,1))*(.5d0*invhy)+(-3*U(i,j,k,1)+4*U(i,1+j,k,1)-
     &  U(i,2+j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
      endif
      if(masks(i,-1+j,-1+k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dwdyb=((-3*U(i,j,-1+k,3)+4*U(i,1+j,-1+k,3)-U(i,2+j
     &  ,-1+k,3))*(.5d0*invhy)+(-3*U(i,j,k,3)+4*U(i,1+j,k,3)-
     &  U(i,2+j,k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
      endif
      if(masks(i,-1+j,k).gt.0.or.masks(i,-1+j,1+k).gt.0) then 
        dwdyt=((-3*U(i,j,k,3)+4*U(i,1+j,k,3)-U(i,2+j,k,3))
     &  *(.5d0*invhy)+(-3*U(i,j,1+k,3)+4*U(i,1+j,1+k,3)-U(i,2
     &  +j,1+k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
      endif
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
        dudze=trandere(1+i,j,k,1,3)
        dvdzs=tranders(i,-1+j,k,2,3)
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maske(i+1,j,k) .gt. 0 ) then
         diagu = diagu -fe(i,j,k,1)*2.d0*beta*muX(1+i,j,k)*ihxsq
         diagv = diagv -fe(i,j,k,2)*(beta*muX(1+i,j,k)*ihxsq)
         diagw = diagw -fe(i,j,k,3)*(beta*muX(1+i,j,k)*ihxsq)
      endif
      if( masks(i,j-1,k) .gt. 0 ) then
         diagu = diagu -fs(i,j,k,1)*(beta*muY(i,j,k)*ihysq)
         diagv = diagv -fs(i,j,k,2)*2.d0*beta*muY(i,j,k)*ihysq
         diagw = diagw -fs(i,j,k,3)*(beta*muY(i,j,k)*ihysq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)
      enddo
      endif



c ::: case 17, west-north edge
c ::: Null
      j=hi(2)
      i=lo(1)
      if( mod(j,2).eq.mody .and. mod(i,2).eq.modx ) then
      do k=kstart,hi(3)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if(maskw(-1+i,-1+j,k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dvdxs=((-3*U(i,-1+j,k,2)+4*U(1+i,-1+j,k,2)-U(2+i,-
     &  1+j,k,2))*(.5d0*invhx)+(-3*U(i,j,k,2)+4*U(1+i,j,k,2)-
     &  U(2+i,j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
      endif
      if(maskw(-1+i,j,-1+k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dwdxb=((-3*U(i,j,-1+k,3)+4*U(1+i,j,-1+k,3)-U(2+i,j
     &  ,-1+k,3))*(.5d0*invhx)+(-3*U(i,j,k,3)+4*U(1+i,j,k,3)-
     &  U(2+i,j,k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
      endif
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,j,1+k).gt.0) then 
        dwdxt=((-3*U(i,j,k,3)+4*U(1+i,j,k,3)-U(2+i,j,k,3))
     &  *(.5d0*invhx)+(-3*U(i,j,1+k,3)+4*U(1+i,j,1+k,3)-U(2+i
     &  ,j,1+k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
      endif
      if(maskn(i,1+j,k).gt.0.or.maskn(1+i,1+j,k).gt.0) then 
        dudye=((U(i,-2+j,k,1)-4*U(i,-1+j,k,1)+3*U(i,j,k,1)
     &  )*(.5d0*invhy)+(U(1+i,-2+j,k,1)-4*U(1+i,-1+j,k,1)+3*U
     &  (1+i,j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
      endif
        dudyw=tranderw(-1+i,j,k,1,2)
      if(maskn(i,1+j,-1+k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dwdyb=((U(i,-2+j,-1+k,3)-4*U(i,-1+j,-1+k,3)+3*U(i,
     &  j,-1+k,3))*(.5d0*invhy)+(U(i,-2+j,k,3)-4*U(i,-1+j,k,3
     &  )+3*U(i,j,k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
      endif
      if(maskn(i,1+j,k).gt.0.or.maskn(i,1+j,1+k).gt.0) then 
        dwdyt=((U(i,-2+j,k,3)-4*U(i,-1+j,k,3)+3*U(i,j,k,3)
     &  )*(.5d0*invhy)+(U(i,-2+j,1+k,3)-4*U(i,-1+j,1+k,3)+3*U
     &  (i,j,1+k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
      endif
        dudzw=tranderw(-1+i,j,k,1,3)
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
        dvdzn=trandern(i,1+j,k,2,3)
        
c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskn(i,j+1,k) .gt. 0 ) then
         diagu = diagu -fn(i,j,k,1)*(beta*muY(i,1+j,k)*ihysq)
         diagv = diagv -fn(i,j,k,2)*2.d0*beta*muY(i,1+j,k)*ihysq
         diagw = diagw -fn(i,j,k,3)*(beta*muY(i,1+j,k)*ihysq)
      endif
      if( maskw(i-1,j,k) .gt. 0 ) then
         diagu = diagu -fw(i,j,k,1)*2.d0*beta*muX(i,j,k)*ihxsq
         diagv = diagv -fw(i,j,k,2)*(beta*muX(i,j,k)*ihxsq)
         diagw = diagw -fw(i,j,k,3)*(beta*muX(i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)
      enddo
      endif



c ::: case 18, west-south edge
c ::: Null
      j=lo(2)
      i=lo(1)
      if( mod(j,2).eq.mody .and. mod(i,2).eq.modx ) then
      do k=kstart,hi(3)-1,2
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,1+j,k).gt.0) then 
        dvdxn=((-3*U(i,j,k,2)+4*U(1+i,j,k,2)-U(2+i,j,k,2))
     &  *(.5d0*invhx)+(-3*U(i,1+j,k,2)+4*U(1+i,1+j,k,2)-U(2+i
     &  ,1+j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
      if(maskw(-1+i,j,-1+k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dwdxb=((-3*U(i,j,-1+k,3)+4*U(1+i,j,-1+k,3)-U(2+i,j
     &  ,-1+k,3))*(.5d0*invhx)+(-3*U(i,j,k,3)+4*U(1+i,j,k,3)-
     &  U(2+i,j,k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
      endif
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,j,1+k).gt.0) then 
        dwdxt=((-3*U(i,j,k,3)+4*U(1+i,j,k,3)-U(2+i,j,k,3))
     &  *(.5d0*invhx)+(-3*U(i,j,1+k,3)+4*U(1+i,j,1+k,3)-U(2+i
     &  ,j,1+k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
      endif
      if(masks(i,-1+j,k).gt.0.or.masks(1+i,-1+j,k).gt.0) then 
        dudye=((-3*U(i,j,k,1)+4*U(i,1+j,k,1)-U(i,2+j,k,1))
     &  *(.5d0*invhy)+(-3*U(1+i,j,k,1)+4*U(1+i,1+j,k,1)-U(1+i
     &  ,2+j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
      endif
        dudyw=tranderw(-1+i,j,k,1,2)
      if(masks(i,-1+j,-1+k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dwdyb=((-3*U(i,j,-1+k,3)+4*U(i,1+j,-1+k,3)-U(i,2+j
     &  ,-1+k,3))*(.5d0*invhy)+(-3*U(i,j,k,3)+4*U(i,1+j,k,3)-
     &  U(i,2+j,k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
      endif
      if(masks(i,-1+j,k).gt.0.or.masks(i,-1+j,1+k).gt.0) then 
        dwdyt=((-3*U(i,j,k,3)+4*U(i,1+j,k,3)-U(i,2+j,k,3))
     &  *(.5d0*invhy)+(-3*U(i,j,1+k,3)+4*U(i,1+j,1+k,3)-U(i,2
     &  +j,1+k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
      endif
        dudzw=tranderw(-1+i,j,k,1,3)
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
        dvdzs=tranders(i,-1+j,k,2,3)
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskw(i-1,j,k) .gt. 0 ) then
         diagu = diagu -fw(i,j,k,1)*2.d0*beta*muX(i,j,k)*ihxsq
         diagv = diagv -fw(i,j,k,2)*(beta*muX(i,j,k)*ihxsq)
         diagw = diagw -fw(i,j,k,3)*(beta*muX(i,j,k)*ihxsq)
      endif
      if( masks(i,j-1,k) .gt. 0 ) then
         diagu = diagu -fw(i,j,k,1)*(beta*muY(i,j,k)*ihysq)
         diagv = diagv -fw(i,j,k,2)*2.d0*beta*muY(i,j,k)*ihysq
         diagw = diagw -fw(i,j,k,3)*(beta*muY(i,j,k)*ihysq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)

      enddo
      endif


c ::: case 19, top-north-east corner
c ::: Null
      k=hi(3)
      j=hi(2)
      i=hi(1)
      if( mod(k,2).eq.modz .and. mod(j,2).eq.mody .and. mod(i,2).eq.modx ) then
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if(maske(1+i,-1+j,k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dvdxs=((U(-2+i,-1+j,k,2)-4*U(-1+i,-1+j,k,2)+3*U(i,
     &  -1+j,k,2))*(.5d0*invhx)+(U(-2+i,j,k,2)-4*U(-1+i,j,k,2
     &  )+3*U(i,j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
      endif
      if(maske(1+i,j,-1+k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dwdxb=((U(-2+i,j,-1+k,3)-4*U(-1+i,j,-1+k,3)+3*U(i,
     &  j,-1+k,3))*(.5d0*invhx)+(U(-2+i,j,k,3)-4*U(-1+i,j,k,3
     &  )+3*U(i,j,k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
      endif
        dwdxt=trandert(i,j,1+k,3,1)
        dudye=trandere(1+i,j,k,1,2)
      if(maskn(-1+i,1+j,k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dudyw=((U(-1+i,-2+j,k,1)-4*U(-1+i,-1+j,k,1)+3*U(-1
     &  +i,j,k,1))*(.5d0*invhy)+(U(i,-2+j,k,1)-4*U(i,-1+j,k,1
     &  )+3*U(i,j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
      endif
      if(maskn(i,1+j,-1+k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dwdyb=((U(i,-2+j,-1+k,3)-4*U(i,-1+j,-1+k,3)+3*U(i,
     &  j,-1+k,3))*(.5d0*invhy)+(U(i,-2+j,k,3)-4*U(i,-1+j,k,3
     &  )+3*U(i,j,k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
      endif
        dwdyt=trandert(i,j,1+k,3,2)
      if(maskt(-1+i,j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dudzw=((U(-1+i,j,-2+k,1)-4*U(-1+i,j,-1+k,1)+3*U(-1
     &  +i,j,k,1))*(.5d0*invhz)+(U(i,j,-2+k,1)-4*U(i,j,-1+k,1
     &  )+3*U(i,j,k,1))*(.5d0*invhz))*0.5d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
      endif
        dudze=trandere(1+i,j,k,1,3)
      if(maskt(i,-1+j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dvdzs=((U(i,-1+j,-2+k,2)-4*U(i,-1+j,-1+k,2)+3*U(i,
     &  -1+j,k,2))*(.5d0*invhz)+(U(i,j,-2+k,2)-4*U(i,j,-1+k,2
     &  )+3*U(i,j,k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
      endif
        dvdzn=trandern(i,1+j,k,2,3)


c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskt(i,j,k+1) .gt. 0 ) then
         diagu = diagu -ft(i,j,k,1)*(beta*muZ(i,j,1+k)*ihzsq)
         diagv = diagv -ft(i,j,k,2)*(beta*muZ(i,j,1+k)*ihzsq)
         diagw = diagw -ft(i,j,k,3)*2.d0*beta*muZ(i,j,1+k)*ihzsq
      endif
      if( maskn(i,j+1,k) .gt. 0 ) then
         diagu = diagu -fn(i,j,k,1)*(beta*muY(i,1+j,k)*ihysq)
         diagv = diagv -fn(i,j,k,2)*2.d0*beta*muY(i,1+j,k)*ihysq
         diagw = diagw -fn(i,j,k,3)*(beta*muY(i,1+j,k)*ihysq)
      endif
      if( maske(i+1,j,k) .gt. 0 ) then
         diagu = diagu -fe(i,j,k,1)*2.d0*beta*muX(1+i,j,k)*ihxsq
         diagv = diagv -fe(i,j,k,2)*(beta*muX(1+i,j,k)*ihxsq)
         diagw = diagw -fe(i,j,k,3)*(beta*muX(1+i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)

       endif



c ::: case 20, top-north-west corner
c ::: Null
           i=lo(1)
           j=hi(2)
           k=hi(3)
           if( mod(i,2).eq.modx .and. mod(j,2).eq.mody .and. mod(k,2).eq.modz ) then
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if(maskw(-1+i,-1+j,k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dvdxs=((-3*U(i,-1+j,k,2)+4*U(1+i,-1+j,k,2)-U(2+i,-
     &  1+j,k,2))*(.5d0*invhx)+(-3*U(i,j,k,2)+4*U(1+i,j,k,2)-
     &  U(2+i,j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
      endif
      if(maskw(-1+i,j,-1+k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dwdxb=((-3*U(i,j,-1+k,3)+4*U(1+i,j,-1+k,3)-U(2+i,j
     &  ,-1+k,3))*(.5d0*invhx)+(-3*U(i,j,k,3)+4*U(1+i,j,k,3)-
     &  U(2+i,j,k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
      endif
        dwdxt=trandert(i,j,1+k,3,1)
      if(maskn(i,1+j,k).gt.0.or.maskn(1+i,1+j,k).gt.0) then 
        dudye=((U(i,-2+j,k,1)-4*U(i,-1+j,k,1)+3*U(i,j,k,1)
     &  )*(.5d0*invhy)+(U(1+i,-2+j,k,1)-4*U(1+i,-1+j,k,1)+3*U
     &  (1+i,j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
      endif
        dudyw=tranderw(-1+i,j,k,1,2)
      if(maskn(i,1+j,-1+k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dwdyb=((U(i,-2+j,-1+k,3)-4*U(i,-1+j,-1+k,3)+3*U(i,
     &  j,-1+k,3))*(.5d0*invhy)+(U(i,-2+j,k,3)-4*U(i,-1+j,k,3
     &  )+3*U(i,j,k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
      endif
        dwdyt=trandert(i,j,1+k,3,2)
        dudzw=tranderw(-1+i,j,k,1,3)
      if(maskt(i,j,1+k).gt.0.or.maskt(1+i,j,1+k).gt.0) then 
        dudze=((U(i,j,-2+k,1)-4*U(i,j,-1+k,1)+3*U(i,j,k,1)
     &  )*(.5d0*invhz)+(U(1+i,j,-2+k,1)-4*U(1+i,j,-1+k,1)+3*U
     &  (1+i,j,k,1))*(.5d0*invhz))*0.5d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
      endif
      if(maskt(i,-1+j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dvdzs=((U(i,-1+j,-2+k,2)-4*U(i,-1+j,-1+k,2)+3*U(i,
     &  -1+j,k,2))*(.5d0*invhz)+(U(i,j,-2+k,2)-4*U(i,j,-1+k,2
     &  )+3*U(i,j,k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskt(i,j,k+1) .gt. 0 ) then
         diagu = diagu -ft(i,j,k,1)*(beta*muZ(i,j,1+k)*ihzsq)
         diagv = diagv -ft(i,j,k,2)*(beta*muZ(i,j,1+k)*ihzsq)
         diagw = diagw -ft(i,j,k,3)*2.d0*beta*muZ(i,j,1+k)*ihzsq
      endif
      if( maskn(i,j+1,k) .gt. 0 ) then
         diagu = diagu -fn(i,j,k,1)*(beta*muY(i,1+j,k)*ihysq)
         diagv = diagv -fn(i,j,k,2)*2.d0*beta*muY(i,1+j,k)*ihysq
         diagw = diagw -fn(i,j,k,3)*(beta*muY(i,1+j,k)*ihysq)
      endif
      if( maskw(i-1,j,k) .gt. 0 ) then
         diagu = diagu -fw(i,j,k,1)*2.d0*beta*muX(i,j,k)*ihxsq
         diagv = diagv -fw(i,j,k,2)*(beta*muX(i,j,k)*ihxsq)
         diagw = diagw -fw(i,j,k,3)*(beta*muX(i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)

       endif



c ::: case 21, top-south-east corner
c ::: Null
           i=hi(1)
           j=lo(2)
           k=hi(3)
           if( mod(i,2).eq.modx .and. mod(j,2).eq.mody .and. mod(k,2).eq.modz ) then
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
      if(maske(1+i,j,k).gt.0.or.maske(1+i,1+j,k).gt.0) then 
        dvdxn=((U(-2+i,j,k,2)-4*U(-1+i,j,k,2)+3*U(i,j,k,2)
     &  )*(.5d0*invhx)+(U(-2+i,1+j,k,2)-4*U(-1+i,1+j,k,2)+3*U
     &  (i,1+j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
      if(maske(1+i,j,-1+k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dwdxb=((U(-2+i,j,-1+k,3)-4*U(-1+i,j,-1+k,3)+3*U(i,
     &  j,-1+k,3))*(.5d0*invhx)+(U(-2+i,j,k,3)-4*U(-1+i,j,k,3
     &  )+3*U(i,j,k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
      endif
        dwdxt=trandert(i,j,1+k,3,1)
        dudye=trandere(1+i,j,k,1,2)
      if(masks(-1+i,-1+j,k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dudyw=((-3*U(-1+i,j,k,1)+4*U(-1+i,1+j,k,1)-U(-1+i,
     &  2+j,k,1))*(.5d0*invhy)+(-3*U(i,j,k,1)+4*U(i,1+j,k,1)-
     &  U(i,2+j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
      endif
      if(masks(i,-1+j,-1+k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dwdyb=((-3*U(i,j,-1+k,3)+4*U(i,1+j,-1+k,3)-U(i,2+j
     &  ,-1+k,3))*(.5d0*invhy)+(-3*U(i,j,k,3)+4*U(i,1+j,k,3)-
     &  U(i,2+j,k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
      endif
        dwdyt=trandert(i,j,1+k,3,2)
      if(maskt(-1+i,j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dudzw=((U(-1+i,j,-2+k,1)-4*U(-1+i,j,-1+k,1)+3*U(-1
     &  +i,j,k,1))*(.5d0*invhz)+(U(i,j,-2+k,1)-4*U(i,j,-1+k,1
     &  )+3*U(i,j,k,1))*(.5d0*invhz))*0.5d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
      endif
        dudze=trandere(1+i,j,k,1,3)
        dvdzs=tranders(i,-1+j,k,2,3)
      if(maskt(i,j,1+k).gt.0.or.maskt(i,1+j,1+k).gt.0) then 
        dvdzn=((U(i,j,-2+k,2)-4*U(i,j,-1+k,2)+3*U(i,j,k,2)
     &  )*(.5d0*invhz)+(U(i,1+j,-2+k,2)-4*U(i,1+j,-1+k,2)+3*U
     &  (i,1+j,k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)
      endif


c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskt(i,j,k+1) .gt. 0 ) then
         diagu = diagu -ft(i,j,k,1)*(beta*muZ(i,j,1+k)*ihzsq)
         diagv = diagv -ft(i,j,k,2)*(beta*muZ(i,j,1+k)*ihzsq)
         diagw = diagw -ft(i,j,k,3)*2.d0*beta*muZ(i,j,1+k)*ihzsq
      endif
      if( masks(i,j-1,k) .gt. 0 ) then
         diagu = diagu -fs(i,j,k,1)*(beta*muY(i,j,k)*ihysq)
         diagv = diagv -fs(i,j,k,2)*2.d0*beta*muY(i,j,k)*ihysq
         diagw = diagw -fs(i,j,k,3)*(beta*muY(i,j,k)*ihysq)
      endif
      if( maske(i+1,j,k) .gt. 0 ) then
         diagu = diagu -fe(i,j,k,1)*2.d0*beta*muX(1+i,j,k)*ihxsq
         diagv = diagv -fe(i,j,k,2)*(beta*muX(1+i,j,k)*ihxsq)
         diagw = diagw -fe(i,j,k,3)*(beta*muX(1+i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)

      endif


c ::: case 22, top-south-west corner
c ::: Null
           i=lo(1)
           j=lo(2)
           k=hi(3)
           if( mod(i,2).eq.modx .and. mod(j,2).eq.mody .and. mod(k,2).eq.modz ) then
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,1+j,k).gt.0) then 
        dvdxn=((-3*U(i,j,k,2)+4*U(1+i,j,k,2)-U(2+i,j,k,2))
     &  *(.5d0*invhx)+(-3*U(i,1+j,k,2)+4*U(1+i,1+j,k,2)-U(2+i
     &  ,1+j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
      if(maskw(-1+i,j,-1+k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dwdxb=((-3*U(i,j,-1+k,3)+4*U(1+i,j,-1+k,3)-U(2+i,j
     &  ,-1+k,3))*(.5d0*invhx)+(-3*U(i,j,k,3)+4*U(1+i,j,k,3)-
     &  U(2+i,j,k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))*(.25d0*invhx)
      endif
        dwdxt=trandert(i,j,1+k,3,1)
      if(masks(i,-1+j,k).gt.0.or.masks(1+i,-1+j,k).gt.0) then 
        dudye=((-3*U(i,j,k,1)+4*U(i,1+j,k,1)-U(i,2+j,k,1))
     &  *(.5d0*invhy)+(-3*U(1+i,j,k,1)+4*U(1+i,1+j,k,1)-U(1+i
     &  ,2+j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
      endif
        dudyw=tranderw(-1+i,j,k,1,2)
      if(masks(i,-1+j,-1+k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dwdyb=((-3*U(i,j,-1+k,3)+4*U(i,1+j,-1+k,3)-U(i,2+j
     &  ,-1+k,3))*(.5d0*invhy)+(-3*U(i,j,k,3)+4*U(i,1+j,k,3)-
     &  U(i,2+j,k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))*(.25d0*invhy)
      endif
        dwdyt=trandert(i,j,1+k,3,2)
        dudzw=tranderw(-1+i,j,k,1,3)
      if(maskt(i,j,1+k).gt.0.or.maskt(1+i,j,1+k).gt.0) then 
        dudze=((U(i,j,-2+k,1)-4*U(i,j,-1+k,1)+3*U(i,j,k,1)
     &  )*(.5d0*invhz)+(U(1+i,j,-2+k,1)-4*U(1+i,j,-1+k,1)+3*U
     &  (1+i,j,k,1))*(.5d0*invhz))*0.5d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
      endif
        dvdzs=tranders(i,-1+j,k,2,3)
      if(maskt(i,j,1+k).gt.0.or.maskt(i,1+j,1+k).gt.0) then 
        dvdzn=((U(i,j,-2+k,2)-4*U(i,j,-1+k,2)+3*U(i,j,k,2)
     &  )*(.5d0*invhz)+(U(i,1+j,-2+k,2)-4*U(i,1+j,-1+k,2)+3*U
     &  (i,1+j,k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)
      endif

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskt(i,j,k+1) .gt. 0 ) then
         diagu = diagu -ft(i,j,k,1)*(beta*muZ(i,j,1+k)*ihzsq)
         diagv = diagv -ft(i,j,k,2)*(beta*muZ(i,j,1+k)*ihzsq)
         diagw = diagw -ft(i,j,k,3)*2.d0*beta*muZ(i,j,1+k)*ihzsq
      endif
      if( masks(i,j-1,k) .gt. 0 ) then
         diagu = diagu -fs(i,j,k,1)*(beta*muY(i,j,k)*ihysq)
         diagv = diagv -fs(i,j,k,2)*2.d0*beta*muY(i,j,k)*ihysq
         diagw = diagw -fs(i,j,k,3)*(beta*muY(i,j,k)*ihysq)
      endif
      if( maskw(i-1,j,k) .gt. 0 ) then
         diagu = diagu -fw(i,j,k,1)*2.d0*beta*muX(i,j,k)*ihxsq
         diagv = diagv -fw(i,j,k,2)*(beta*muX(i,j,k)*ihxsq)
         diagw = diagw -fw(i,j,k,3)*(beta*muX(i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)

        endif



c ::: case 23, bottom-north-east corner
c ::: Null
           i=hi(1)
           j=hi(2)
           k=lo(3)
           if( mod(i,2).eq.modx .and. mod(j,2).eq.mody .and. mod(k,2).eq.modz ) then
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if(maske(1+i,-1+j,k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dvdxs=((U(-2+i,-1+j,k,2)-4*U(-1+i,-1+j,k,2)+3*U(i,
     &  -1+j,k,2))*(.5d0*invhx)+(U(-2+i,j,k,2)-4*U(-1+i,j,k,2
     &  )+3*U(i,j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
      endif
        dwdxb=tranderb(i,j,-1+k,3,1)
      if(maske(1+i,j,k).gt.0.or.maske(1+i,j,1+k).gt.0) then 
        dwdxt=((U(-2+i,j,k,3)-4*U(-1+i,j,k,3)+3*U(i,j,k,3)
     &  )*(.5d0*invhx)+(U(-2+i,j,1+k,3)-4*U(-1+i,j,1+k,3)+3*U
     &  (i,j,1+k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
      endif
        dudye=trandere(1+i,j,k,1,2)
      if(maskn(-1+i,1+j,k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dudyw=((U(-1+i,-2+j,k,1)-4*U(-1+i,-1+j,k,1)+3*U(-1
     &  +i,j,k,1))*(.5d0*invhy)+(U(i,-2+j,k,1)-4*U(i,-1+j,k,1
     &  )+3*U(i,j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
      endif
        dwdyb=tranderb(i,j,-1+k,3,2)
      if(maskn(i,1+j,k).gt.0.or.maskn(i,1+j,1+k).gt.0) then 
        dwdyt=((U(i,-2+j,k,3)-4*U(i,-1+j,k,3)+3*U(i,j,k,3)
     &  )*(.5d0*invhy)+(U(i,-2+j,1+k,3)-4*U(i,-1+j,1+k,3)+3*U
     &  (i,j,1+k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
      endif
      if(maskb(-1+i,j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dudzw=((-3*U(-1+i,j,k,1)+4*U(-1+i,j,1+k,1)-U(-1+i,
     &  j,2+k,1))*(.5d0*invhz)+(-3*U(i,j,k,1)+4*U(i,j,1+k,1)-
     &  U(i,j,2+k,1))*(.5d0*invhz))*0.5d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
      endif
        dudze=trandere(1+i,j,k,1,3)
      if(maskb(i,-1+j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dvdzs=((-3*U(i,-1+j,k,2)+4*U(i,-1+j,1+k,2)-U(i,-1+
     &  j,2+k,2))*(.5d0*invhz)+(-3*U(i,j,k,2)+4*U(i,j,1+k,2)-
     &  U(i,j,2+k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
      endif
        dvdzn=trandern(i,1+j,k,2,3)


c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskb(i,j,k-1) .gt. 0 ) then
         diagu = diagu -fb(i,j,k,1)*(beta*muZ(i,j,k)*ihzsq)
         diagv = diagv -fb(i,j,k,2)*(beta*muZ(i,j,k)*ihzsq)
         diagw = diagw -fb(i,j,k,3)*2.d0*beta*muZ(i,j,k)*ihzsq
      endif
      if( maskn(i,j+1,k) .gt. 0 ) then
         diagu = diagu -fn(i,j,k,1)*(beta*muY(i,1+j,k)*ihysq)
         diagv = diagv -fn(i,j,k,2)*2.d0*beta*muY(i,1+j,k)*ihysq
         diagw = diagw -fn(i,j,k,3)*(beta*muY(i,1+j,k)*ihysq)
      endif
      if( maske(i+1,j,k) .gt. 0 ) then
         diagu = diagu -fe(i,j,k,1)*2.d0*beta*muX(1+i,j,k)*ihxsq
         diagv = diagv -fe(i,j,k,2)*(beta*muX(1+i,j,k)*ihxsq)
         diagw = diagw -fe(i,j,k,3)*(beta*muX(1+i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)


       endif


c ::: case 24, bottom-north-west corner
c ::: Null
           i=lo(1)
           j=hi(2)
           k=lo(3)
           if( mod(i,2).eq.modx .and. mod(j,2).eq.mody .and. mod(k,2).eq.modz ) then
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if(maskw(-1+i,-1+j,k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dvdxs=((-3*U(i,-1+j,k,2)+4*U(1+i,-1+j,k,2)-U(2+i,-
     &  1+j,k,2))*(.5d0*invhx)+(-3*U(i,j,k,2)+4*U(1+i,j,k,2)-
     &  U(2+i,j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))*(.25d0*invhx)
      endif
        dwdxb=tranderb(i,j,-1+k,3,1)
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,j,1+k).gt.0) then 
        dwdxt=((-3*U(i,j,k,3)+4*U(1+i,j,k,3)-U(2+i,j,k,3))
     &  *(.5d0*invhx)+(-3*U(i,j,1+k,3)+4*U(1+i,j,1+k,3)-U(2+i
     &  ,j,1+k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
      endif
      if(maskn(i,1+j,k).gt.0.or.maskn(1+i,1+j,k).gt.0) then 
        dudye=((U(i,-2+j,k,1)-4*U(i,-1+j,k,1)+3*U(i,j,k,1)
     &  )*(.5d0*invhy)+(U(1+i,-2+j,k,1)-4*U(1+i,-1+j,k,1)+3*U
     &  (1+i,j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
      endif
        dudyw=tranderw(-1+i,j,k,1,2)
        dwdyb=tranderb(i,j,-1+k,3,2)
      if(maskn(i,1+j,k).gt.0.or.maskn(i,1+j,1+k).gt.0) then 
        dwdyt=((U(i,-2+j,k,3)-4*U(i,-1+j,k,3)+3*U(i,j,k,3)
     &  )*(.5d0*invhy)+(U(i,-2+j,1+k,3)-4*U(i,-1+j,1+k,3)+3*U
     &  (i,j,1+k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
      endif
        dudzw=tranderw(-1+i,j,k,1,3)
      if(maskb(i,j,-1+k).gt.0.or.maskb(1+i,j,-1+k).gt.0) then 
        dudze=((-3*U(i,j,k,1)+4*U(i,j,1+k,1)-U(i,j,2+k,1))
     &  *(.5d0*invhz)+(-3*U(1+i,j,k,1)+4*U(1+i,j,1+k,1)-U(1+i
     &  ,j,2+k,1))*(.5d0*invhz))*0.5d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
      endif
      if(maskb(i,-1+j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dvdzs=((-3*U(i,-1+j,k,2)+4*U(i,-1+j,1+k,2)-U(i,-1+
     &  j,2+k,2))*(.5d0*invhz)+(-3*U(i,j,k,2)+4*U(i,j,1+k,2)-
     &  U(i,j,2+k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))*(.25d0*invhz)
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskb(i,j,k-1) .gt. 0 ) then
         diagu = diagu -fb(i,j,k,1)*(beta*muZ(i,j,k)*ihzsq)
         diagv = diagv -fb(i,j,k,2)*(beta*muZ(i,j,k)*ihzsq)
         diagw = diagw -fb(i,j,k,3)*2.d0*beta*muZ(i,j,k)*ihzsq
      endif
      if( maskn(i,j+1,k) .gt. 0 ) then
         diagu = diagu -fn(i,j,k,1)*(beta*muY(i,1+j,k)*ihysq)
         diagv = diagv -fn(i,j,k,2)*2.d0*beta*muY(i,1+j,k)*ihysq
         diagw = diagw -fn(i,j,k,3)*(beta*muY(i,1+j,k)*ihysq)
      endif
      if( maskw(i-1,j,k) .gt. 0 ) then
         diagu = diagu -fw(i,j,k,1)*2.d0*beta*muX(i,j,k)*ihxsq
         diagv = diagv -fw(i,j,k,2)*(beta*muX(i,j,k)*ihxsq)
         diagw = diagw -fw(i,j,k,3)*(beta*muX(i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)


       endif



c ::: case 25, bottom-south-east corner
c ::: Null
           i=hi(1)
           j=lo(2)
           k=lo(3)
           if( mod(i,2).eq.modx .and. mod(j,2).eq.mody .and. mod(k,2).eq.modz ) then
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
      if(maske(1+i,j,k).gt.0.or.maske(1+i,1+j,k).gt.0) then 
        dvdxn=((U(-2+i,j,k,2)-4*U(-1+i,j,k,2)+3*U(i,j,k,2)
     &  )*(.5d0*invhx)+(U(-2+i,1+j,k,2)-4*U(-1+i,1+j,k,2)+3*U
     &  (i,1+j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxb=tranderb(i,j,-1+k,3,1)
      if(maske(1+i,j,k).gt.0.or.maske(1+i,j,1+k).gt.0) then 
        dwdxt=((U(-2+i,j,k,3)-4*U(-1+i,j,k,3)+3*U(i,j,k,3)
     &  )*(.5d0*invhx)+(U(-2+i,j,1+k,3)-4*U(-1+i,j,1+k,3)+3*U
     &  (i,j,1+k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
      endif
        dudye=trandere(1+i,j,k,1,2)
      if(masks(-1+i,-1+j,k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dudyw=((-3*U(-1+i,j,k,1)+4*U(-1+i,1+j,k,1)-U(-1+i,
     &  2+j,k,1))*(.5d0*invhy)+(-3*U(i,j,k,1)+4*U(i,1+j,k,1)-
     &  U(i,2+j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))*(.25d0*invhy)
      endif
        dwdyb=tranderb(i,j,-1+k,3,2)
      if(masks(i,-1+j,k).gt.0.or.masks(i,-1+j,1+k).gt.0) then 
        dwdyt=((-3*U(i,j,k,3)+4*U(i,1+j,k,3)-U(i,2+j,k,3))
     &  *(.5d0*invhy)+(-3*U(i,j,1+k,3)+4*U(i,1+j,1+k,3)-U(i,2
     &  +j,1+k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
      endif
      if(maskb(-1+i,j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dudzw=((-3*U(-1+i,j,k,1)+4*U(-1+i,j,1+k,1)-U(-1+i,
     &  j,2+k,1))*(.5d0*invhz)+(-3*U(i,j,k,1)+4*U(i,j,1+k,1)-
     &  U(i,j,2+k,1))*(.5d0*invhz))*0.5d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))*(.25d0*invhz)
      endif
        dudze=trandere(1+i,j,k,1,3)
        dvdzs=tranders(i,-1+j,k,2,3)
      if(maskb(i,j,-1+k).gt.0.or.maskb(i,1+j,-1+k).gt.0) then 
        dvdzn=((-3*U(i,j,k,2)+4*U(i,j,1+k,2)-U(i,j,2+k,2))
     &  *(.5d0*invhz)+(-3*U(i,1+j,k,2)+4*U(i,1+j,1+k,2)-U(i,1
     &  +j,2+k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)
      endif

c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskb(i,j,k-1) .gt. 0 ) then
         diagu = diagu -fb(i,j,k,1)*(beta*muZ(i,j,k)*ihzsq)
         diagv = diagv -fb(i,j,k,2)*(beta*muZ(i,j,k)*ihzsq)
         diagw = diagw -fb(i,j,k,3)*2.d0*beta*muZ(i,j,k)*ihzsq
      endif
      if( masks(i,j-1,k) .gt. 0 ) then
         diagu = diagu -fs(i,j,k,1)*(beta*muY(i,j,k)*ihysq)
         diagv = diagv -fs(i,j,k,2)*2.d0*beta*muY(i,j,k)*ihysq
         diagw = diagw -fs(i,j,k,3)*(beta*muY(i,j,k)*ihysq)
      endif
      if( maske(i+1,j,k) .gt. 0 ) then
         diagu = diagu -fe(i,j,k,1)*2.d0*beta*muX(1+i,j,k)*ihxsq
         diagv = diagv -fe(i,j,k,2)*(beta*muX(1+i,j,k)*ihxsq)
         diagw = diagw -fe(i,j,k,3)*(beta*muX(1+i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)

       endif



c ::: case 26, bottom-south-west corner
c ::: Null
           i=lo(1)
           j=lo(2)
           k=lo(3)
           if( mod(i,2).eq.modx .and. mod(j,2).eq.mody .and. mod(k,2).eq.modz ) then
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))*invhx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))*invhx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))*invhy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))*invhy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))*invhz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))*invhz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))*invhy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))*invhy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))*invhx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))*invhx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))*invhz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))*invhz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))*invhz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))*invhz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))*invhx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))*invhx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))*invhy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))*invhy

c     ::: get expressions for tangential derivatives
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,1+j,k).gt.0) then 
        dvdxn=((-3*U(i,j,k,2)+4*U(1+i,j,k,2)-U(2+i,j,k,2))
     &  *(.5d0*invhx)+(-3*U(i,1+j,k,2)+4*U(1+i,1+j,k,2)-U(2+i
     &  ,1+j,k,2))*(.5d0*invhx))*0.5d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))*(.25d0*invhx)
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxb=tranderb(i,j,-1+k,3,1)
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,j,1+k).gt.0) then 
        dwdxt=((-3*U(i,j,k,3)+4*U(1+i,j,k,3)-U(2+i,j,k,3))
     &  *(.5d0*invhx)+(-3*U(i,j,1+k,3)+4*U(1+i,j,1+k,3)-U(2+i
     &  ,j,1+k,3))*(.5d0*invhx))*0.5d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))*(.25d0*invhx)
      endif
      if(masks(i,-1+j,k).gt.0.or.masks(1+i,-1+j,k).gt.0) then 
        dudye=((-3*U(i,j,k,1)+4*U(i,1+j,k,1)-U(i,2+j,k,1))
     &  *(.5d0*invhy)+(-3*U(1+i,j,k,1)+4*U(1+i,1+j,k,1)-U(1+i
     &  ,2+j,k,1))*(.5d0*invhy))*0.5d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))*(.25d0*invhy)
      endif
        dudyw=tranderw(-1+i,j,k,1,2)
        dwdyb=tranderb(i,j,-1+k,3,2)
      if(masks(i,-1+j,k).gt.0.or.masks(i,-1+j,1+k).gt.0) then 
        dwdyt=((-3*U(i,j,k,3)+4*U(i,1+j,k,3)-U(i,2+j,k,3))
     &  *(.5d0*invhy)+(-3*U(i,j,1+k,3)+4*U(i,1+j,1+k,3)-U(i,2
     &  +j,1+k,3))*(.5d0*invhy))*0.5d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))*(.25d0*invhy)
      endif
        dudzw=tranderw(-1+i,j,k,1,3)
      if(maskb(i,j,-1+k).gt.0.or.maskb(1+i,j,-1+k).gt.0) then 
        dudze=((-3*U(i,j,k,1)+4*U(i,j,1+k,1)-U(i,j,2+k,1))
     &  *(.5d0*invhz)+(-3*U(1+i,j,k,1)+4*U(1+i,j,1+k,1)-U(1+i
     &  ,j,2+k,1))*(.5d0*invhz))*0.5d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))*(.25d0*invhz)
      endif
        dvdzs=tranders(i,-1+j,k,2,3)
      if(maskb(i,j,-1+k).gt.0.or.maskb(i,1+j,-1+k).gt.0) then 
        dvdzn=((-3*U(i,j,k,2)+4*U(i,j,1+k,2)-U(i,j,2+k,2))
     &  *(.5d0*invhz)+(-3*U(i,1+j,k,2)+4*U(i,1+j,1+k,2)-U(i,1
     &  +j,2+k,2))*(.5d0*invhz))*0.5d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))*(.25d0*invhz)
      endif


c ::: get diagonal element
      diagu = alpha*a(i,j,k)+2.d0*beta*muX(i,j,k)*ihxsq+2.d0*beta*muX(1
     &  +i,j,k)*ihxsq+beta*muY(i,j,k)*ihysq+beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagv = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +2.d0*beta*muY(i,j,k)*ihysq+2.d0*beta*muY(i,1+
     &  j,k)*ihysq+beta*muZ(i,j,k)*ihzsq+beta*muZ(i,j,1+k
     &  )*ihzsq
      diagw = alpha*a(i,j,k)+beta*muX(i,j,k)*ihxsq+beta*muX(1+i,j,k)*ihxsq
     &  +beta*muY(i,j,k)*ihysq+beta*muY(i,1+j,k)*ihysq
     &  +2.d0*beta*muZ(i,j,k)*ihzsq+2.d0*beta*muZ(i,j,1+k
     &  )*ihzsq
c ::: correct diagonal elements for ghost cell contrib
      if( maskb(i,j,k-1) .gt. 0 ) then
         diagu = diagu -fb(i,j,k,1)*(beta*muZ(i,j,k)*ihzsq)
         diagv = diagv -fb(i,j,k,2)*(beta*muZ(i,j,k)*ihzsq)
         diagw = diagw -fb(i,j,k,3)*2.d0*beta*muZ(i,j,k)*ihzsq
      endif
      if( masks(i,j-1,k) .gt. 0 ) then
         diagu = diagu -fs(i,j,k,1)*(beta*muY(i,j,k)*ihysq)
         diagv = diagv -fs(i,j,k,2)*2.d0*beta*muY(i,j,k)*ihysq
         diagw = diagw -fs(i,j,k,3)*(beta*muY(i,j,k)*ihysq)
      endif
      if( maskw(i-1,j,k) .gt. 0 ) then
         diagu = diagu -fw(i,j,k,1)*2.d0*beta*muX(i,j,k)*ihxsq
         diagv = diagv -fw(i,j,k,2)*(beta*muX(i,j,k)*ihxsq)
         diagw = diagw -fw(i,j,k,3)*(beta*muX(i,j,k)*ihxsq)
      endif

c     :::  evaluate expression
      operu = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,1)
      operv = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,2)
      operw = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))*ihxyz)+alpha*a(i,
     &  j,k)*u(i,j,k,3)
c ::: do relaxation
               u(i,j,k,1) = (rhs(i,j,k,1)-operu)/diagu+u(i,j,k,1)
               u(i,j,k,2) = (rhs(i,j,k,2)-operv)/diagv+u(i,j,k,2)
               u(i,j,k,3) = (rhs(i,j,k,3)-operw)/diagw+u(i,j,k,3)

      endif





       return
       end
