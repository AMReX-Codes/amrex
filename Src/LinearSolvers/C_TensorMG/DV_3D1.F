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
c ::: 2) the normal derivatives are evaluated in the normal fashion.
c ::: 3) tangential derivatives which reach outside the rectangle DO have
c :::    to check the masks.  

      subroutine FORT_DVAPPLY(
     $     u, DIMS(u),
     $     alpha, beta,
     $     a, DIMS(a),
     $     muX, DIMS(muX),
     $     muY, DIMS(muY),
     $     muZ, DIMS(muZ),
     $     out, DIMS(out),
     $     maskn,DIMS(maskn),
     $     maske,DIMS(maske),
     $     maskw,DIMS(maskw),
     $     masks,DIMS(masks),
     $     maskt,DIMS(maskt),
     $     maskb,DIMS(maskb),
     $     trandern,DIMS(trandern),
     $     trandere,DIMS(trandere),
     $     tranderw,DIMS(tranderw),
     $     tranders,DIMS(tranders),
     $     trandert,DIMS(trandert),
     $     tranderb,DIMS(tranderb),
     $     lo,hi,h
     $     )
      REAL_T alpha, beta
      integer DIMDEC(u)
      integer DIMDEC(a)
      integer DIMDEC(muX)
      integer DIMDEC(muY)
      integer DIMDEC(muZ)
      integer DIMDEC(out)
      integer DIMDEC(maskn)
      integer DIMDEC(maske)
      integer DIMDEC(maskw)
      integer DIMDEC(masks)
      integer DIMDEC(maskt)
      integer DIMDEC(maskb)
      integer DIMDEC(trandern)
      integer DIMDEC(trandere)
      integer DIMDEC(tranderw)
      integer DIMDEC(tranders)
      integer DIMDEC(trandert)
      integer DIMDEC(tranderb)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      REAL_T h(BL_SPACEDIM)
      REAL_T u(DIMV(u),3)
      REAL_T a(DIMV(a))
      REAL_T muX(DIMV(muX))
      REAL_T muY(DIMV(muY))
      REAL_T muZ(DIMV(muZ))
      REAL_T out(DIMV(out),3)
      integer  maskn(DIMV(maskn))
      integer  maske(DIMV(maske))
      integer  maskw(DIMV(maskw))
      integer  masks(DIMV(masks))
      integer  maskt(DIMV(maskt))
      integer  maskb(DIMV(maskb))

c ::: for transverse derives, first 3 is for variable, second 3 is for
c ::: direction of derivative.  Obviously not all are used, but this is easy.
      REAL_T trandern(DIMV(trandern),3,3)
      REAL_T trandere(DIMV(trandere),3,3)
      REAL_T tranderw(DIMV(tranderw),3,3)
      REAL_T tranders(DIMV(tranders),3,3)
      REAL_T trandert(DIMV(trandert),3,3)
      REAL_T tranderb(DIMV(tranderb),3,3)

      integer i,j,k
      REAL_T tauxxw,tauxxe,tauyyn,tauyys
      REAL_T tauxye,tauxyw,tauxyn,tauxys
      REAL_T hx,hy,hz
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

      logical False, True
      parameter( False=.false.)
      parameter( True=.true. )

      hx = h(1)
      hy = h(2)
      hz = h(3)

c ::: cases:
c     ::: 0) interior.  All usual derivative expressions
c     ::: 1-6) face but not on edge.  
c     ::: 7-18) edge but not on corner.
c     ::: 19-26) corners.

c ::: case 0
c ::: Null
c ::: Null
      do k=lo(3)+1,hi(3)-1
         do j=lo(2)+1,hi(2)-1
            do i=lo(1)+1,hi(1)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
            enddo
         enddo
      enddo

c ::: case 1, top face
c ::: Null
c ::: Null
c ::: Null
      k=hi(3)
         do j=lo(2)+1,hi(2)-1
            do i=lo(1)+1,hi(1)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
        dwdxt=trandert(i,j,1+k,3,1)
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
        dwdyt=trandert(i,j,1+k,3,2)
      if(maskt(-1+i,j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dudzw=((U(-1+i,j,-2+k,1)-4*U(-1+i,j,-1+k,1)+3*U(-1
     &  +i,j,k,1))/(2.d0*hz)+(U(i,j,-2+k,1)-4*U(i,j,-1+k,1
     &  )+3*U(i,j,k,1))/(2.d0*hz))/2.d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
      endif
      if(maskt(i,j,1+k).gt.0.or.maskt(1+i,j,1+k).gt.0) then 
        dudze=((U(i,j,-2+k,1)-4*U(i,j,-1+k,1)+3*U(i,j,k,1)
     &  )/(2.d0*hz)+(U(1+i,j,-2+k,1)-4*U(1+i,j,-1+k,1)+3*U
     &  (1+i,j,k,1))/(2.d0*hz))/2.d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
      endif
      if(maskt(i,-1+j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dvdzs=((U(i,-1+j,-2+k,2)-4*U(i,-1+j,-1+k,2)+3*U(i,
     &  -1+j,k,2))/(2.d0*hz)+(U(i,j,-2+k,2)-4*U(i,j,-1+k,2
     &  )+3*U(i,j,k,2))/(2.d0*hz))/2.d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
      endif
      if(maskt(i,j,1+k).gt.0.or.maskt(i,1+j,1+k).gt.0) then 
        dvdzn=((U(i,j,-2+k,2)-4*U(i,j,-1+k,2)+3*U(i,j,k,2)
     &  )/(2.d0*hz)+(U(i,1+j,-2+k,2)-4*U(i,1+j,-1+k,2)+3*U
     &  (i,1+j,k,2))/(2.d0*hz))/2.d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)
      endif

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
            enddo
         enddo

c ::: case 2, bottom face
c ::: Null
c ::: Null
c ::: Null
      k=lo(3)
         do j=lo(2)+1,hi(2)-1
            do i=lo(1)+1,hi(1)-1
c     :::  get expressions for normal derivatives
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
        dwdxb=tranderb(i,j,-1+k,3,1)
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
        dwdyb=tranderb(i,j,-1+k,3,2)
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
      if(maskb(-1+i,j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dudzw=((-3*U(-1+i,j,k,1)+4*U(-1+i,j,1+k,1)-U(-1+i,
     &  j,2+k,1))/(2.d0*hz)+(-3*U(i,j,k,1)+4*U(i,j,1+k,1)-
     &  U(i,j,2+k,1))/(2.d0*hz))/2.d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
      endif
      if(maskb(i,j,-1+k).gt.0.or.maskb(1+i,j,-1+k).gt.0) then 
        dudze=((-3*U(i,j,k,1)+4*U(i,j,1+k,1)-U(i,j,2+k,1))
     &  /(2.d0*hz)+(-3*U(1+i,j,k,1)+4*U(1+i,j,1+k,1)-U(1+i
     &  ,j,2+k,1))/(2.d0*hz))/2.d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
      endif
      if(maskb(i,-1+j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dvdzs=((-3*U(i,-1+j,k,2)+4*U(i,-1+j,1+k,2)-U(i,-1+
     &  j,2+k,2))/(2.d0*hz)+(-3*U(i,j,k,2)+4*U(i,j,1+k,2)-
     &  U(i,j,2+k,2))/(2.d0*hz))/2.d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
      endif
      if(maskb(i,j,-1+k).gt.0.or.maskb(i,1+j,-1+k).gt.0) then 
        dvdzn=((-3*U(i,j,k,2)+4*U(i,j,1+k,2)-U(i,j,2+k,2))
     &  /(2.d0*hz)+(-3*U(i,1+j,k,2)+4*U(i,1+j,1+k,2)-U(i,1
     &  +j,2+k,2))/(2.d0*hz))/2.d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)
      endif

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
            enddo
         enddo

c ::: case 3, west face
c ::: Null
c ::: Null
c ::: Null
      i = lo(1)
      do k=lo(3)+1,hi(3)-1
         do j=lo(2)+1,hi(2)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
        dudyw=tranderw(-1+i,j,k,1,2)
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,1+j,k).gt.0) then 
        dvdxn=((-3*U(i,j,k,2)+4*U(1+i,j,k,2)-U(2+i,j,k,2))
     &  /(2.d0*hx)+(-3*U(i,1+j,k,2)+4*U(1+i,1+j,k,2)-U(2+i
     &  ,1+j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
      endif
      if(maskw(-1+i,-1+j,k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dvdxs=((-3*U(i,-1+j,k,2)+4*U(1+i,-1+j,k,2)-U(2+i,-
     &  1+j,k,2))/(2.d0*hx)+(-3*U(i,j,k,2)+4*U(1+i,j,k,2)-
     &  U(2+i,j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
      endif
      if(maskw(-1+i,j,-1+k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dwdxb=((-3*U(i,j,-1+k,3)+4*U(1+i,j,-1+k,3)-U(2+i,j
     &  ,-1+k,3))/(2.d0*hx)+(-3*U(i,j,k,3)+4*U(1+i,j,k,3)-
     &  U(2+i,j,k,3))/(2.d0*hx))/2.d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
      endif
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,j,1+k).gt.0) then 
        dwdxt=((-3*U(i,j,k,3)+4*U(1+i,j,k,3)-U(2+i,j,k,3))
     &  /(2.d0*hx)+(-3*U(i,j,1+k,3)+4*U(1+i,j,1+k,3)-U(2+i
     &  ,j,1+k,3))/(2.d0*hx))/2.d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
      endif
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
        dudzw=tranderw(-1+i,j,k,1,3)
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
         enddo
      enddo

c ::: case 4, east face
c ::: Null
c ::: Null
c ::: Null
      i = hi(1)
      do k=lo(3)+1,hi(3)-1
         do j=lo(2)+1,hi(2)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dudye=trandere(1+i,j,k,1,2)
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
      if(maske(1+i,j,k).gt.0.or.maske(1+i,1+j,k).gt.0) then 
        dvdxn=((U(-2+i,j,k,2)-4*U(-1+i,j,k,2)+3*U(i,j,k,2)
     &  )/(2.d0*hx)+(U(-2+i,1+j,k,2)-4*U(-1+i,1+j,k,2)+3*U
     &  (i,1+j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
      endif
      if(maske(1+i,-1+j,k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dvdxs=((U(-2+i,-1+j,k,2)-4*U(-1+i,-1+j,k,2)+3*U(i,
     &  -1+j,k,2))/(2.d0*hx)+(U(-2+i,j,k,2)-4*U(-1+i,j,k,2
     &  )+3*U(i,j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
      endif
      if(maske(1+i,j,-1+k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dwdxb=((U(-2+i,j,-1+k,3)-4*U(-1+i,j,-1+k,3)+3*U(i,
     &  j,-1+k,3))/(2.d0*hx)+(U(-2+i,j,k,3)-4*U(-1+i,j,k,3
     &  )+3*U(i,j,k,3))/(2.d0*hx))/2.d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
      endif
      if(maske(1+i,j,k).gt.0.or.maske(1+i,j,1+k).gt.0) then 
        dwdxt=((U(-2+i,j,k,3)-4*U(-1+i,j,k,3)+3*U(i,j,k,3)
     &  )/(2.d0*hx)+(U(-2+i,j,1+k,3)-4*U(-1+i,j,1+k,3)+3*U
     &  (i,j,1+k,3))/(2.d0*hx))/2.d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
      endif
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
        dudze=trandere(1+i,j,k,1,3)
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
         enddo
      enddo
	
c ::: case 5, north face
c ::: Null
c ::: Null
c ::: Null
      j = hi(2)
      do k=lo(3)+1,hi(3)-1
            do i=lo(1)+1,hi(1)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if(maskn(i,1+j,k).gt.0.or.maskn(1+i,1+j,k).gt.0) then 
        dudye=((U(i,-2+j,k,1)-4*U(i,-1+j,k,1)+3*U(i,j,k,1)
     &  )/(2.d0*hy)+(U(1+i,-2+j,k,1)-4*U(1+i,-1+j,k,1)+3*U
     &  (1+i,j,k,1))/(2.d0*hy))/2.d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
      endif
      if(maskn(-1+i,1+j,k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dudyw=((U(-1+i,-2+j,k,1)-4*U(-1+i,-1+j,k,1)+3*U(-1
     &  +i,j,k,1))/(2.d0*hy)+(U(i,-2+j,k,1)-4*U(i,-1+j,k,1
     &  )+3*U(i,j,k,1))/(2.d0*hy))/2.d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
      endif
        dvdxn=trandern(i,1+j,k,2,1)
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
      if(maskn(i,1+j,-1+k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dwdyb=((U(i,-2+j,-1+k,3)-4*U(i,-1+j,-1+k,3)+3*U(i,
     &  j,-1+k,3))/(2.d0*hy)+(U(i,-2+j,k,3)-4*U(i,-1+j,k,3
     &  )+3*U(i,j,k,3))/(2.d0*hy))/2.d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
      endif
      if(maskn(i,1+j,k).gt.0.or.maskn(i,1+j,1+k).gt.0) then 
        dwdyt=((U(i,-2+j,k,3)-4*U(i,-1+j,k,3)+3*U(i,j,k,3)
     &  )/(2.d0*hy)+(U(i,-2+j,1+k,3)-4*U(i,-1+j,1+k,3)+3*U
     &  (i,j,1+k,3))/(2.d0*hy))/2.d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
      endif
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
            enddo
      enddo


c ::: case 6, south face
c ::: Null
c ::: Null
c ::: Null
      j = lo(2)
      do k=lo(3)+1,hi(3)-1
            do i=lo(1)+1,hi(1)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if(masks(i,-1+j,k).gt.0.or.masks(1+i,-1+j,k).gt.0) then 
        dudye=((-3*U(i,j,k,1)+4*U(i,1+j,k,1)-U(i,2+j,k,1))
     &  /(2.d0*hy)+(-3*U(1+i,j,k,1)+4*U(1+i,1+j,k,1)-U(1+i
     &  ,2+j,k,1))/(2.d0*hy))/2.d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
      endif
      if(masks(-1+i,-1+j,k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dudyw=((-3*U(-1+i,j,k,1)+4*U(-1+i,1+j,k,1)-U(-1+i,
     &  2+j,k,1))/(2.d0*hy)+(-3*U(i,j,k,1)+4*U(i,1+j,k,1)-
     &  U(i,2+j,k,1))/(2.d0*hy))/2.d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
      endif
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
      if(masks(i,-1+j,-1+k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dwdyb=((-3*U(i,j,-1+k,3)+4*U(i,1+j,-1+k,3)-U(i,2+j
     &  ,-1+k,3))/(2.d0*hy)+(-3*U(i,j,k,3)+4*U(i,1+j,k,3)-
     &  U(i,2+j,k,3))/(2.d0*hy))/2.d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
      endif
      if(masks(i,-1+j,k).gt.0.or.masks(i,-1+j,1+k).gt.0) then 
        dwdyt=((-3*U(i,j,k,3)+4*U(i,1+j,k,3)-U(i,2+j,k,3))
     &  /(2.d0*hy)+(-3*U(i,j,1+k,3)+4*U(i,1+j,1+k,3)-U(i,2
     &  +j,1+k,3))/(2.d0*hy))/2.d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
      endif
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
        dvdzs=tranders(i,-1+j,k,2,3)
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
            enddo
      enddo

c ::: case 7, top-north edge 
c ::: Null
c ::: Null
c ::: Null
      k=hi(3)
         j=hi(2)
            do i=lo(1)+1,hi(1)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if(maskn(i,1+j,k).gt.0.or.maskn(1+i,1+j,k).gt.0) then 
        dudye=((U(i,-2+j,k,1)-4*U(i,-1+j,k,1)+3*U(i,j,k,1)
     &  )/(2.d0*hy)+(U(1+i,-2+j,k,1)-4*U(1+i,-1+j,k,1)+3*U
     &  (1+i,j,k,1))/(2.d0*hy))/2.d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
      endif
      if(maskn(-1+i,1+j,k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dudyw=((U(-1+i,-2+j,k,1)-4*U(-1+i,-1+j,k,1)+3*U(-1
     &  +i,j,k,1))/(2.d0*hy)+(U(i,-2+j,k,1)-4*U(i,-1+j,k,1
     &  )+3*U(i,j,k,1))/(2.d0*hy))/2.d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
      endif
        dvdxn=trandern(i,1+j,k,2,1)
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
        dwdxt=trandert(i,j,1+k,3,1)
      if(maskn(i,1+j,-1+k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dwdyb=((U(i,-2+j,-1+k,3)-4*U(i,-1+j,-1+k,3)+3*U(i,
     &  j,-1+k,3))/(2.d0*hy)+(U(i,-2+j,k,3)-4*U(i,-1+j,k,3
     &  )+3*U(i,j,k,3))/(2.d0*hy))/2.d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
      endif
        dwdyt=trandert(i,j,1+k,3,2)
      if(maskt(-1+i,j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dudzw=((U(-1+i,j,-2+k,1)-4*U(-1+i,j,-1+k,1)+3*U(-1
     &  +i,j,k,1))/(2.d0*hz)+(U(i,j,-2+k,1)-4*U(i,j,-1+k,1
     &  )+3*U(i,j,k,1))/(2.d0*hz))/2.d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
      endif
      if(maskt(i,j,1+k).gt.0.or.maskt(1+i,j,1+k).gt.0) then 
        dudze=((U(i,j,-2+k,1)-4*U(i,j,-1+k,1)+3*U(i,j,k,1)
     &  )/(2.d0*hz)+(U(1+i,j,-2+k,1)-4*U(1+i,j,-1+k,1)+3*U
     &  (1+i,j,k,1))/(2.d0*hz))/2.d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
      endif
      if(maskt(i,-1+j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dvdzs=((U(i,-1+j,-2+k,2)-4*U(i,-1+j,-1+k,2)+3*U(i,
     &  -1+j,k,2))/(2.d0*hz)+(U(i,j,-2+k,2)-4*U(i,j,-1+k,2
     &  )+3*U(i,j,k,2))/(2.d0*hz))/2.d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
            enddo


c ::: case 8, top-south edge 
c ::: Null
      k=hi(3)
         j=lo(2)
            do i=lo(1)+1,hi(1)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if(masks(i,-1+j,k).gt.0.or.masks(1+i,-1+j,k).gt.0) then 
        dudye=((-3*U(i,j,k,1)+4*U(i,1+j,k,1)-U(i,2+j,k,1))
     &  /(2.d0*hy)+(-3*U(1+i,j,k,1)+4*U(1+i,1+j,k,1)-U(1+i
     &  ,2+j,k,1))/(2.d0*hy))/2.d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
      endif
      if(masks(-1+i,-1+j,k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dudyw=((-3*U(-1+i,j,k,1)+4*U(-1+i,1+j,k,1)-U(-1+i,
     &  2+j,k,1))/(2.d0*hy)+(-3*U(i,j,k,1)+4*U(i,1+j,k,1)-
     &  U(i,2+j,k,1))/(2.d0*hy))/2.d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
      endif
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
        dwdxt=trandert(i,j,1+k,3,1)
      if(masks(i,-1+j,-1+k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dwdyb=((-3*U(i,j,-1+k,3)+4*U(i,1+j,-1+k,3)-U(i,2+j
     &  ,-1+k,3))/(2.d0*hy)+(-3*U(i,j,k,3)+4*U(i,1+j,k,3)-
     &  U(i,2+j,k,3))/(2.d0*hy))/2.d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
      endif
        dwdyt=trandert(i,j,1+k,3,2)
      if(maskt(-1+i,j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dudzw=((U(-1+i,j,-2+k,1)-4*U(-1+i,j,-1+k,1)+3*U(-1
     &  +i,j,k,1))/(2.d0*hz)+(U(i,j,-2+k,1)-4*U(i,j,-1+k,1
     &  )+3*U(i,j,k,1))/(2.d0*hz))/2.d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
      endif
      if(maskt(i,j,1+k).gt.0.or.maskt(1+i,j,1+k).gt.0) then 
        dudze=((U(i,j,-2+k,1)-4*U(i,j,-1+k,1)+3*U(i,j,k,1)
     &  )/(2.d0*hz)+(U(1+i,j,-2+k,1)-4*U(1+i,j,-1+k,1)+3*U
     &  (1+i,j,k,1))/(2.d0*hz))/2.d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
      endif
        dvdzs=tranders(i,-1+j,k,2,3)
      if(maskt(i,j,1+k).gt.0.or.maskt(i,1+j,1+k).gt.0) then 
        dvdzn=((U(i,j,-2+k,2)-4*U(i,j,-1+k,2)+3*U(i,j,k,2)
     &  )/(2.d0*hz)+(U(i,1+j,-2+k,2)-4*U(i,1+j,-1+k,2)+3*U
     &  (i,1+j,k,2))/(2.d0*hz))/2.d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)
      endif

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
            enddo


c ::: case 9, top-west edge 
c ::: Null
      k=hi(3)
      i=lo(1)
         do j=lo(2)+1,hi(2)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
        dudyw=tranderw(-1+i,j,k,1,2)
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,1+j,k).gt.0) then 
        dvdxn=((-3*U(i,j,k,2)+4*U(1+i,j,k,2)-U(2+i,j,k,2))
     &  /(2.d0*hx)+(-3*U(i,1+j,k,2)+4*U(1+i,1+j,k,2)-U(2+i
     &  ,1+j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
      endif
      if(maskw(-1+i,-1+j,k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dvdxs=((-3*U(i,-1+j,k,2)+4*U(1+i,-1+j,k,2)-U(2+i,-
     &  1+j,k,2))/(2.d0*hx)+(-3*U(i,j,k,2)+4*U(1+i,j,k,2)-
     &  U(2+i,j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
      endif
      if(maskw(-1+i,j,-1+k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dwdxb=((-3*U(i,j,-1+k,3)+4*U(1+i,j,-1+k,3)-U(2+i,j
     &  ,-1+k,3))/(2.d0*hx)+(-3*U(i,j,k,3)+4*U(1+i,j,k,3)-
     &  U(2+i,j,k,3))/(2.d0*hx))/2.d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
      endif
        dwdxt=trandert(i,j,1+k,3,1)
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
        dwdyt=trandert(i,j,1+k,3,2)
        dudzw=tranderw(-1+i,j,k,1,3)
      if(maskt(i,j,1+k).gt.0.or.maskt(1+i,j,1+k).gt.0) then 
        dudze=((U(i,j,-2+k,1)-4*U(i,j,-1+k,1)+3*U(i,j,k,1)
     &  )/(2.d0*hz)+(U(1+i,j,-2+k,1)-4*U(1+i,j,-1+k,1)+3*U
     &  (1+i,j,k,1))/(2.d0*hz))/2.d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
      endif
      if(maskt(i,-1+j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dvdzs=((U(i,-1+j,-2+k,2)-4*U(i,-1+j,-1+k,2)+3*U(i,
     &  -1+j,k,2))/(2.d0*hz)+(U(i,j,-2+k,2)-4*U(i,j,-1+k,2
     &  )+3*U(i,j,k,2))/(2.d0*hz))/2.d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
      endif
      if(maskt(i,j,1+k).gt.0.or.maskt(i,1+j,1+k).gt.0) then 
        dvdzn=((U(i,j,-2+k,2)-4*U(i,j,-1+k,2)+3*U(i,j,k,2)
     &  )/(2.d0*hz)+(U(i,1+j,-2+k,2)-4*U(i,1+j,-1+k,2)+3*U
     &  (i,1+j,k,2))/(2.d0*hz))/2.d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)
      endif

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
         enddo

c ::: case 10, top-east edge 
c ::: Null
      k=hi(3)
      i=hi(1)
         do j=lo(2)+1,hi(2)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dudye=trandere(1+i,j,k,1,2)
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
      if(maske(1+i,j,k).gt.0.or.maske(1+i,1+j,k).gt.0) then 
        dvdxn=((U(-2+i,j,k,2)-4*U(-1+i,j,k,2)+3*U(i,j,k,2)
     &  )/(2.d0*hx)+(U(-2+i,1+j,k,2)-4*U(-1+i,1+j,k,2)+3*U
     &  (i,1+j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
      endif
      if(maske(1+i,-1+j,k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dvdxs=((U(-2+i,-1+j,k,2)-4*U(-1+i,-1+j,k,2)+3*U(i,
     &  -1+j,k,2))/(2.d0*hx)+(U(-2+i,j,k,2)-4*U(-1+i,j,k,2
     &  )+3*U(i,j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
      endif
      if(maske(1+i,j,-1+k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dwdxb=((U(-2+i,j,-1+k,3)-4*U(-1+i,j,-1+k,3)+3*U(i,
     &  j,-1+k,3))/(2.d0*hx)+(U(-2+i,j,k,3)-4*U(-1+i,j,k,3
     &  )+3*U(i,j,k,3))/(2.d0*hx))/2.d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
      endif
        dwdxt=trandert(i,j,1+k,3,1)
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
        dwdyt=trandert(i,j,1+k,3,2)
      if(maskt(-1+i,j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dudzw=((U(-1+i,j,-2+k,1)-4*U(-1+i,j,-1+k,1)+3*U(-1
     &  +i,j,k,1))/(2.d0*hz)+(U(i,j,-2+k,1)-4*U(i,j,-1+k,1
     &  )+3*U(i,j,k,1))/(2.d0*hz))/2.d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
      endif
        dudze=trandere(1+i,j,k,1,3)
      if(maskt(i,-1+j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dvdzs=((U(i,-1+j,-2+k,2)-4*U(i,-1+j,-1+k,2)+3*U(i,
     &  -1+j,k,2))/(2.d0*hz)+(U(i,j,-2+k,2)-4*U(i,j,-1+k,2
     &  )+3*U(i,j,k,2))/(2.d0*hz))/2.d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
      endif
      if(maskt(i,j,1+k).gt.0.or.maskt(i,1+j,1+k).gt.0) then 
        dvdzn=((U(i,j,-2+k,2)-4*U(i,j,-1+k,2)+3*U(i,j,k,2)
     &  )/(2.d0*hz)+(U(i,1+j,-2+k,2)-4*U(i,1+j,-1+k,2)+3*U
     &  (i,1+j,k,2))/(2.d0*hz))/2.d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)
      endif

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
         enddo


c ::: case 11, bottom-north edge 
c ::: Null
      k=lo(3)
      j=hi(2)
            do i=lo(1)+1,hi(1)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if(maskn(i,1+j,k).gt.0.or.maskn(1+i,1+j,k).gt.0) then 
        dudye=((U(i,-2+j,k,1)-4*U(i,-1+j,k,1)+3*U(i,j,k,1)
     &  )/(2.d0*hy)+(U(1+i,-2+j,k,1)-4*U(1+i,-1+j,k,1)+3*U
     &  (1+i,j,k,1))/(2.d0*hy))/2.d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
      endif
      if(maskn(-1+i,1+j,k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dudyw=((U(-1+i,-2+j,k,1)-4*U(-1+i,-1+j,k,1)+3*U(-1
     &  +i,j,k,1))/(2.d0*hy)+(U(i,-2+j,k,1)-4*U(i,-1+j,k,1
     &  )+3*U(i,j,k,1))/(2.d0*hy))/2.d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
      endif
        dvdxn=trandern(i,1+j,k,2,1)
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
        dwdxb=tranderb(i,j,-1+k,3,1)
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
        dwdyb=tranderb(i,j,-1+k,3,2)
      if(maskn(i,1+j,k).gt.0.or.maskn(i,1+j,1+k).gt.0) then 
        dwdyt=((U(i,-2+j,k,3)-4*U(i,-1+j,k,3)+3*U(i,j,k,3)
     &  )/(2.d0*hy)+(U(i,-2+j,1+k,3)-4*U(i,-1+j,1+k,3)+3*U
     &  (i,j,1+k,3))/(2.d0*hy))/2.d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
      endif
      if(maskb(-1+i,j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dudzw=((-3*U(-1+i,j,k,1)+4*U(-1+i,j,1+k,1)-U(-1+i,
     &  j,2+k,1))/(2.d0*hz)+(-3*U(i,j,k,1)+4*U(i,j,1+k,1)-
     &  U(i,j,2+k,1))/(2.d0*hz))/2.d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
      endif
      if(maskb(i,j,-1+k).gt.0.or.maskb(1+i,j,-1+k).gt.0) then 
        dudze=((-3*U(i,j,k,1)+4*U(i,j,1+k,1)-U(i,j,2+k,1))
     &  /(2.d0*hz)+(-3*U(1+i,j,k,1)+4*U(1+i,j,1+k,1)-U(1+i
     &  ,j,2+k,1))/(2.d0*hz))/2.d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
      endif
      if(maskb(i,-1+j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dvdzs=((-3*U(i,-1+j,k,2)+4*U(i,-1+j,1+k,2)-U(i,-1+
     &  j,2+k,2))/(2.d0*hz)+(-3*U(i,j,k,2)+4*U(i,j,1+k,2)-
     &  U(i,j,2+k,2))/(2.d0*hz))/2.d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
            enddo


c ::: case 12, bottom-south edge 
c ::: Null
      k=lo(3)
      j=lo(2)
            do i=lo(1)+1,hi(1)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if(masks(i,-1+j,k).gt.0.or.masks(1+i,-1+j,k).gt.0) then 
        dudye=((-3*U(i,j,k,1)+4*U(i,1+j,k,1)-U(i,2+j,k,1))
     &  /(2.d0*hy)+(-3*U(1+i,j,k,1)+4*U(1+i,1+j,k,1)-U(1+i
     &  ,2+j,k,1))/(2.d0*hy))/2.d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
      endif
      if(masks(-1+i,-1+j,k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dudyw=((-3*U(-1+i,j,k,1)+4*U(-1+i,1+j,k,1)-U(-1+i,
     &  2+j,k,1))/(2.d0*hy)+(-3*U(i,j,k,1)+4*U(i,1+j,k,1)-
     &  U(i,2+j,k,1))/(2.d0*hy))/2.d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
      endif
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxb=tranderb(i,j,-1+k,3,1)
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
        dwdyb=tranderb(i,j,-1+k,3,2)
      if(masks(i,-1+j,k).gt.0.or.masks(i,-1+j,1+k).gt.0) then 
        dwdyt=((-3*U(i,j,k,3)+4*U(i,1+j,k,3)-U(i,2+j,k,3))
     &  /(2.d0*hy)+(-3*U(i,j,1+k,3)+4*U(i,1+j,1+k,3)-U(i,2
     &  +j,1+k,3))/(2.d0*hy))/2.d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
      endif
      if(maskb(-1+i,j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dudzw=((-3*U(-1+i,j,k,1)+4*U(-1+i,j,1+k,1)-U(-1+i,
     &  j,2+k,1))/(2.d0*hz)+(-3*U(i,j,k,1)+4*U(i,j,1+k,1)-
     &  U(i,j,2+k,1))/(2.d0*hz))/2.d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
      endif
      if(maskb(i,j,-1+k).gt.0.or.maskb(1+i,j,-1+k).gt.0) then 
        dudze=((-3*U(i,j,k,1)+4*U(i,j,1+k,1)-U(i,j,2+k,1))
     &  /(2.d0*hz)+(-3*U(1+i,j,k,1)+4*U(1+i,j,1+k,1)-U(1+i
     &  ,j,2+k,1))/(2.d0*hz))/2.d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
      endif
        dvdzs=tranders(i,-1+j,k,2,3)
      if(maskb(i,j,-1+k).gt.0.or.maskb(i,1+j,-1+k).gt.0) then 
        dvdzn=((-3*U(i,j,k,2)+4*U(i,j,1+k,2)-U(i,j,2+k,2))
     &  /(2.d0*hz)+(-3*U(i,1+j,k,2)+4*U(i,1+j,1+k,2)-U(i,1
     &  +j,2+k,2))/(2.d0*hz))/2.d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)
      endif

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
            enddo



c ::: case 13, bottom-west edge 
c ::: Null
      k=lo(3)
      i=lo(1)
      do j=lo(2)+1,hi(2)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
        dudyw=tranderw(-1+i,j,k,1,2)
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,1+j,k).gt.0) then 
        dvdxn=((-3*U(i,j,k,2)+4*U(1+i,j,k,2)-U(2+i,j,k,2))
     &  /(2.d0*hx)+(-3*U(i,1+j,k,2)+4*U(1+i,1+j,k,2)-U(2+i
     &  ,1+j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
      endif
      if(maskw(-1+i,-1+j,k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dvdxs=((-3*U(i,-1+j,k,2)+4*U(1+i,-1+j,k,2)-U(2+i,-
     &  1+j,k,2))/(2.d0*hx)+(-3*U(i,j,k,2)+4*U(1+i,j,k,2)-
     &  U(2+i,j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
      endif
        dwdxb=tranderb(i,j,-1+k,3,1)
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,j,1+k).gt.0) then 
        dwdxt=((-3*U(i,j,k,3)+4*U(1+i,j,k,3)-U(2+i,j,k,3))
     &  /(2.d0*hx)+(-3*U(i,j,1+k,3)+4*U(1+i,j,1+k,3)-U(2+i
     &  ,j,1+k,3))/(2.d0*hx))/2.d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
      endif
        dwdyb=tranderb(i,j,-1+k,3,2)
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
        dudzw=tranderw(-1+i,j,k,1,3)
      if(maskb(i,j,-1+k).gt.0.or.maskb(1+i,j,-1+k).gt.0) then 
        dudze=((-3*U(i,j,k,1)+4*U(i,j,1+k,1)-U(i,j,2+k,1))
     &  /(2.d0*hz)+(-3*U(1+i,j,k,1)+4*U(1+i,j,1+k,1)-U(1+i
     &  ,j,2+k,1))/(2.d0*hz))/2.d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
      endif
      if(maskb(i,-1+j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dvdzs=((-3*U(i,-1+j,k,2)+4*U(i,-1+j,1+k,2)-U(i,-1+
     &  j,2+k,2))/(2.d0*hz)+(-3*U(i,j,k,2)+4*U(i,j,1+k,2)-
     &  U(i,j,2+k,2))/(2.d0*hz))/2.d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
      endif
      if(maskb(i,j,-1+k).gt.0.or.maskb(i,1+j,-1+k).gt.0) then 
        dvdzn=((-3*U(i,j,k,2)+4*U(i,j,1+k,2)-U(i,j,2+k,2))
     &  /(2.d0*hz)+(-3*U(i,1+j,k,2)+4*U(i,1+j,1+k,2)-U(i,1
     &  +j,2+k,2))/(2.d0*hz))/2.d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)
      endif

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
         enddo



c ::: case 14, bottom-east edge 
c ::: Null
      k=lo(3)
      i=hi(1)
         do j=lo(2)+1,hi(2)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if(maske(1+i,j,k).gt.0.or.maske(1+i,1+j,k).gt.0) then 
        dvdxn=((U(-2+i,j,k,2)-4*U(-1+i,j,k,2)+3*U(i,j,k,2)
     &  )/(2.d0*hx)+(U(-2+i,1+j,k,2)-4*U(-1+i,1+j,k,2)+3*U
     &  (i,1+j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
      endif
      if(maske(1+i,-1+j,k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dvdxs=((U(-2+i,-1+j,k,2)-4*U(-1+i,-1+j,k,2)+3*U(i,
     &  -1+j,k,2))/(2.d0*hx)+(U(-2+i,j,k,2)-4*U(-1+i,j,k,2
     &  )+3*U(i,j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
      endif
        dwdxb=tranderb(i,j,-1+k,3,1)
      if(maske(1+i,j,k).gt.0.or.maske(1+i,j,1+k).gt.0) then 
        dwdxt=((U(-2+i,j,k,3)-4*U(-1+i,j,k,3)+3*U(i,j,k,3)
     &  )/(2.d0*hx)+(U(-2+i,j,1+k,3)-4*U(-1+i,j,1+k,3)+3*U
     &  (i,j,1+k,3))/(2.d0*hx))/2.d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
      endif
        dudye=trandere(1+i,j,k,1,2)
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
        dwdyb=tranderb(i,j,-1+k,3,2)
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
      if(maskb(-1+i,j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dudzw=((-3*U(-1+i,j,k,1)+4*U(-1+i,j,1+k,1)-U(-1+i,
     &  j,2+k,1))/(2.d0*hz)+(-3*U(i,j,k,1)+4*U(i,j,1+k,1)-
     &  U(i,j,2+k,1))/(2.d0*hz))/2.d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
      endif
        dudze=trandere(1+i,j,k,1,3)
      if(maskb(i,-1+j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dvdzs=((-3*U(i,-1+j,k,2)+4*U(i,-1+j,1+k,2)-U(i,-1+
     &  j,2+k,2))/(2.d0*hz)+(-3*U(i,j,k,2)+4*U(i,j,1+k,2)-
     &  U(i,j,2+k,2))/(2.d0*hz))/2.d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
      endif
      if(maskb(i,j,-1+k).gt.0.or.maskb(i,1+j,-1+k).gt.0) then 
        dvdzn=((-3*U(i,j,k,2)+4*U(i,j,1+k,2)-U(i,j,2+k,2))
     &  /(2.d0*hz)+(-3*U(i,1+j,k,2)+4*U(i,1+j,1+k,2)-U(i,1
     &  +j,2+k,2))/(2.d0*hz))/2.d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)
      endif

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
         enddo


c ::: case 15, east-north edge
c ::: Null
      j=hi(2)
      i=hi(1)
      do k=lo(3)+1,hi(3)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if(maske(1+i,-1+j,k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dvdxs=((U(-2+i,-1+j,k,2)-4*U(-1+i,-1+j,k,2)+3*U(i,
     &  -1+j,k,2))/(2.d0*hx)+(U(-2+i,j,k,2)-4*U(-1+i,j,k,2
     &  )+3*U(i,j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
      endif
      if(maske(1+i,j,-1+k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dwdxb=((U(-2+i,j,-1+k,3)-4*U(-1+i,j,-1+k,3)+3*U(i,
     &  j,-1+k,3))/(2.d0*hx)+(U(-2+i,j,k,3)-4*U(-1+i,j,k,3
     &  )+3*U(i,j,k,3))/(2.d0*hx))/2.d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
      endif
      if(maske(1+i,j,k).gt.0.or.maske(1+i,j,1+k).gt.0) then 
        dwdxt=((U(-2+i,j,k,3)-4*U(-1+i,j,k,3)+3*U(i,j,k,3)
     &  )/(2.d0*hx)+(U(-2+i,j,1+k,3)-4*U(-1+i,j,1+k,3)+3*U
     &  (i,j,1+k,3))/(2.d0*hx))/2.d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
      endif
        dudye=trandere(1+i,j,k,1,2)
      if(maskn(-1+i,1+j,k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dudyw=((U(-1+i,-2+j,k,1)-4*U(-1+i,-1+j,k,1)+3*U(-1
     &  +i,j,k,1))/(2.d0*hy)+(U(i,-2+j,k,1)-4*U(i,-1+j,k,1
     &  )+3*U(i,j,k,1))/(2.d0*hy))/2.d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
      endif
      if(maskn(i,1+j,-1+k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dwdyb=((U(i,-2+j,-1+k,3)-4*U(i,-1+j,-1+k,3)+3*U(i,
     &  j,-1+k,3))/(2.d0*hy)+(U(i,-2+j,k,3)-4*U(i,-1+j,k,3
     &  )+3*U(i,j,k,3))/(2.d0*hy))/2.d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
      endif
      if(maskn(i,1+j,k).gt.0.or.maskn(i,1+j,1+k).gt.0) then 
        dwdyt=((U(i,-2+j,k,3)-4*U(i,-1+j,k,3)+3*U(i,j,k,3)
     &  )/(2.d0*hy)+(U(i,-2+j,1+k,3)-4*U(i,-1+j,1+k,3)+3*U
     &  (i,j,1+k,3))/(2.d0*hy))/2.d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
      endif
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
        dudze=trandere(1+i,j,k,1,3)
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
      enddo


c ::: case 16, east-south edge
c ::: Null
      j=lo(2)
      i=hi(1)
      do k=lo(3)+1,hi(3)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if(maske(1+i,j,k).gt.0.or.maske(1+i,1+j,k).gt.0) then 
        dvdxn=((U(-2+i,j,k,2)-4*U(-1+i,j,k,2)+3*U(i,j,k,2)
     &  )/(2.d0*hx)+(U(-2+i,1+j,k,2)-4*U(-1+i,1+j,k,2)+3*U
     &  (i,1+j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
      if(maske(1+i,j,-1+k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dwdxb=((U(-2+i,j,-1+k,3)-4*U(-1+i,j,-1+k,3)+3*U(i,
     &  j,-1+k,3))/(2.d0*hx)+(U(-2+i,j,k,3)-4*U(-1+i,j,k,3
     &  )+3*U(i,j,k,3))/(2.d0*hx))/2.d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
      endif
      if(maske(1+i,j,k).gt.0.or.maske(1+i,j,1+k).gt.0) then 
        dwdxt=((U(-2+i,j,k,3)-4*U(-1+i,j,k,3)+3*U(i,j,k,3)
     &  )/(2.d0*hx)+(U(-2+i,j,1+k,3)-4*U(-1+i,j,1+k,3)+3*U
     &  (i,j,1+k,3))/(2.d0*hx))/2.d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
      endif
        dudye=trandere(1+i,j,k,1,2)
      if(masks(-1+i,-1+j,k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dudyw=((-3*U(-1+i,j,k,1)+4*U(-1+i,1+j,k,1)-U(-1+i,
     &  2+j,k,1))/(2.d0*hy)+(-3*U(i,j,k,1)+4*U(i,1+j,k,1)-
     &  U(i,2+j,k,1))/(2.d0*hy))/2.d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
      endif
      if(masks(i,-1+j,-1+k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dwdyb=((-3*U(i,j,-1+k,3)+4*U(i,1+j,-1+k,3)-U(i,2+j
     &  ,-1+k,3))/(2.d0*hy)+(-3*U(i,j,k,3)+4*U(i,1+j,k,3)-
     &  U(i,2+j,k,3))/(2.d0*hy))/2.d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
      endif
      if(masks(i,-1+j,k).gt.0.or.masks(i,-1+j,1+k).gt.0) then 
        dwdyt=((-3*U(i,j,k,3)+4*U(i,1+j,k,3)-U(i,2+j,k,3))
     &  /(2.d0*hy)+(-3*U(i,j,1+k,3)+4*U(i,1+j,1+k,3)-U(i,2
     &  +j,1+k,3))/(2.d0*hy))/2.d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
      endif
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
        dudze=trandere(1+i,j,k,1,3)
        dvdzs=tranders(i,-1+j,k,2,3)
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
      enddo


c ::: case 17, west-north edge
c ::: Null
      j=hi(2)
      i=lo(1)
      do k=lo(3)+1,hi(3)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if(maskw(-1+i,-1+j,k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dvdxs=((-3*U(i,-1+j,k,2)+4*U(1+i,-1+j,k,2)-U(2+i,-
     &  1+j,k,2))/(2.d0*hx)+(-3*U(i,j,k,2)+4*U(1+i,j,k,2)-
     &  U(2+i,j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
      endif
      if(maskw(-1+i,j,-1+k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dwdxb=((-3*U(i,j,-1+k,3)+4*U(1+i,j,-1+k,3)-U(2+i,j
     &  ,-1+k,3))/(2.d0*hx)+(-3*U(i,j,k,3)+4*U(1+i,j,k,3)-
     &  U(2+i,j,k,3))/(2.d0*hx))/2.d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
      endif
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,j,1+k).gt.0) then 
        dwdxt=((-3*U(i,j,k,3)+4*U(1+i,j,k,3)-U(2+i,j,k,3))
     &  /(2.d0*hx)+(-3*U(i,j,1+k,3)+4*U(1+i,j,1+k,3)-U(2+i
     &  ,j,1+k,3))/(2.d0*hx))/2.d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
      endif
      if(maskn(i,1+j,k).gt.0.or.maskn(1+i,1+j,k).gt.0) then 
        dudye=((U(i,-2+j,k,1)-4*U(i,-1+j,k,1)+3*U(i,j,k,1)
     &  )/(2.d0*hy)+(U(1+i,-2+j,k,1)-4*U(1+i,-1+j,k,1)+3*U
     &  (1+i,j,k,1))/(2.d0*hy))/2.d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
      endif
        dudyw=tranderw(-1+i,j,k,1,2)
      if(maskn(i,1+j,-1+k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dwdyb=((U(i,-2+j,-1+k,3)-4*U(i,-1+j,-1+k,3)+3*U(i,
     &  j,-1+k,3))/(2.d0*hy)+(U(i,-2+j,k,3)-4*U(i,-1+j,k,3
     &  )+3*U(i,j,k,3))/(2.d0*hy))/2.d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
      endif
      if(maskn(i,1+j,k).gt.0.or.maskn(i,1+j,1+k).gt.0) then 
        dwdyt=((U(i,-2+j,k,3)-4*U(i,-1+j,k,3)+3*U(i,j,k,3)
     &  )/(2.d0*hy)+(U(i,-2+j,1+k,3)-4*U(i,-1+j,1+k,3)+3*U
     &  (i,j,1+k,3))/(2.d0*hy))/2.d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
      endif
        dudzw=tranderw(-1+i,j,k,1,3)
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
      enddo



c ::: case 18, west-south edge
c ::: Null
      j=lo(2)
      i=lo(1)
      do k=lo(3)+1,hi(3)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,1+j,k).gt.0) then 
        dvdxn=((-3*U(i,j,k,2)+4*U(1+i,j,k,2)-U(2+i,j,k,2))
     &  /(2.d0*hx)+(-3*U(i,1+j,k,2)+4*U(1+i,1+j,k,2)-U(2+i
     &  ,1+j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
      if(maskw(-1+i,j,-1+k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dwdxb=((-3*U(i,j,-1+k,3)+4*U(1+i,j,-1+k,3)-U(2+i,j
     &  ,-1+k,3))/(2.d0*hx)+(-3*U(i,j,k,3)+4*U(1+i,j,k,3)-
     &  U(2+i,j,k,3))/(2.d0*hx))/2.d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
      endif
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,j,1+k).gt.0) then 
        dwdxt=((-3*U(i,j,k,3)+4*U(1+i,j,k,3)-U(2+i,j,k,3))
     &  /(2.d0*hx)+(-3*U(i,j,1+k,3)+4*U(1+i,j,1+k,3)-U(2+i
     &  ,j,1+k,3))/(2.d0*hx))/2.d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
      endif
      if(masks(i,-1+j,k).gt.0.or.masks(1+i,-1+j,k).gt.0) then 
        dudye=((-3*U(i,j,k,1)+4*U(i,1+j,k,1)-U(i,2+j,k,1))
     &  /(2.d0*hy)+(-3*U(1+i,j,k,1)+4*U(1+i,1+j,k,1)-U(1+i
     &  ,2+j,k,1))/(2.d0*hy))/2.d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
      endif
        dudyw=tranderw(-1+i,j,k,1,2)
      if(masks(i,-1+j,-1+k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dwdyb=((-3*U(i,j,-1+k,3)+4*U(i,1+j,-1+k,3)-U(i,2+j
     &  ,-1+k,3))/(2.d0*hy)+(-3*U(i,j,k,3)+4*U(i,1+j,k,3)-
     &  U(i,2+j,k,3))/(2.d0*hy))/2.d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
      endif
      if(masks(i,-1+j,k).gt.0.or.masks(i,-1+j,1+k).gt.0) then 
        dwdyt=((-3*U(i,j,k,3)+4*U(i,1+j,k,3)-U(i,2+j,k,3))
     &  /(2.d0*hy)+(-3*U(i,j,1+k,3)+4*U(i,1+j,1+k,3)-U(i,2
     &  +j,1+k,3))/(2.d0*hy))/2.d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
      endif
        dudzw=tranderw(-1+i,j,k,1,3)
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
        dvdzs=tranders(i,-1+j,k,2,3)
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)
      enddo

c ::: case 19, top-north-east corner
c ::: Null
      k=hi(3)
      j=hi(2)
      i=hi(1)
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if(maske(1+i,-1+j,k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dvdxs=((U(-2+i,-1+j,k,2)-4*U(-1+i,-1+j,k,2)+3*U(i,
     &  -1+j,k,2))/(2.d0*hx)+(U(-2+i,j,k,2)-4*U(-1+i,j,k,2
     &  )+3*U(i,j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
      endif
      if(maske(1+i,j,-1+k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dwdxb=((U(-2+i,j,-1+k,3)-4*U(-1+i,j,-1+k,3)+3*U(i,
     &  j,-1+k,3))/(2.d0*hx)+(U(-2+i,j,k,3)-4*U(-1+i,j,k,3
     &  )+3*U(i,j,k,3))/(2.d0*hx))/2.d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
      endif
        dwdxt=trandert(i,j,1+k,3,1)
        dudye=trandere(1+i,j,k,1,2)
      if(maskn(-1+i,1+j,k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dudyw=((U(-1+i,-2+j,k,1)-4*U(-1+i,-1+j,k,1)+3*U(-1
     &  +i,j,k,1))/(2.d0*hy)+(U(i,-2+j,k,1)-4*U(i,-1+j,k,1
     &  )+3*U(i,j,k,1))/(2.d0*hy))/2.d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
      endif
      if(maskn(i,1+j,-1+k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dwdyb=((U(i,-2+j,-1+k,3)-4*U(i,-1+j,-1+k,3)+3*U(i,
     &  j,-1+k,3))/(2.d0*hy)+(U(i,-2+j,k,3)-4*U(i,-1+j,k,3
     &  )+3*U(i,j,k,3))/(2.d0*hy))/2.d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
      endif
        dwdyt=trandert(i,j,1+k,3,2)
      if(maskt(-1+i,j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dudzw=((U(-1+i,j,-2+k,1)-4*U(-1+i,j,-1+k,1)+3*U(-1
     &  +i,j,k,1))/(2.d0*hz)+(U(i,j,-2+k,1)-4*U(i,j,-1+k,1
     &  )+3*U(i,j,k,1))/(2.d0*hz))/2.d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
      endif
        dudze=trandere(1+i,j,k,1,3)
      if(maskt(i,-1+j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dvdzs=((U(i,-1+j,-2+k,2)-4*U(i,-1+j,-1+k,2)+3*U(i,
     &  -1+j,k,2))/(2.d0*hz)+(U(i,j,-2+k,2)-4*U(i,j,-1+k,2
     &  )+3*U(i,j,k,2))/(2.d0*hz))/2.d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)


c ::: case 20, top-north-west corner
c ::: Null
           i=lo(1)
           j=hi(2)
           k=hi(3)
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if(maskw(-1+i,-1+j,k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dvdxs=((-3*U(i,-1+j,k,2)+4*U(1+i,-1+j,k,2)-U(2+i,-
     &  1+j,k,2))/(2.d0*hx)+(-3*U(i,j,k,2)+4*U(1+i,j,k,2)-
     &  U(2+i,j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
      endif
      if(maskw(-1+i,j,-1+k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dwdxb=((-3*U(i,j,-1+k,3)+4*U(1+i,j,-1+k,3)-U(2+i,j
     &  ,-1+k,3))/(2.d0*hx)+(-3*U(i,j,k,3)+4*U(1+i,j,k,3)-
     &  U(2+i,j,k,3))/(2.d0*hx))/2.d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
      endif
        dwdxt=trandert(i,j,1+k,3,1)
      if(maskn(i,1+j,k).gt.0.or.maskn(1+i,1+j,k).gt.0) then 
        dudye=((U(i,-2+j,k,1)-4*U(i,-1+j,k,1)+3*U(i,j,k,1)
     &  )/(2.d0*hy)+(U(1+i,-2+j,k,1)-4*U(1+i,-1+j,k,1)+3*U
     &  (1+i,j,k,1))/(2.d0*hy))/2.d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
      endif
        dudyw=tranderw(-1+i,j,k,1,2)
      if(maskn(i,1+j,-1+k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dwdyb=((U(i,-2+j,-1+k,3)-4*U(i,-1+j,-1+k,3)+3*U(i,
     &  j,-1+k,3))/(2.d0*hy)+(U(i,-2+j,k,3)-4*U(i,-1+j,k,3
     &  )+3*U(i,j,k,3))/(2.d0*hy))/2.d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
      endif
        dwdyt=trandert(i,j,1+k,3,2)
        dudzw=tranderw(-1+i,j,k,1,3)
      if(maskt(i,j,1+k).gt.0.or.maskt(1+i,j,1+k).gt.0) then 
        dudze=((U(i,j,-2+k,1)-4*U(i,j,-1+k,1)+3*U(i,j,k,1)
     &  )/(2.d0*hz)+(U(1+i,j,-2+k,1)-4*U(1+i,j,-1+k,1)+3*U
     &  (1+i,j,k,1))/(2.d0*hz))/2.d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
      endif
      if(maskt(i,-1+j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dvdzs=((U(i,-1+j,-2+k,2)-4*U(i,-1+j,-1+k,2)+3*U(i,
     &  -1+j,k,2))/(2.d0*hz)+(U(i,j,-2+k,2)-4*U(i,j,-1+k,2
     &  )+3*U(i,j,k,2))/(2.d0*hz))/2.d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)



c ::: case 21, top-south-east corner
c ::: Null
           i=hi(1)
           j=lo(2)
           k=hi(3)
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if(maske(1+i,j,k).gt.0.or.maske(1+i,1+j,k).gt.0) then 
        dvdxn=((U(-2+i,j,k,2)-4*U(-1+i,j,k,2)+3*U(i,j,k,2)
     &  )/(2.d0*hx)+(U(-2+i,1+j,k,2)-4*U(-1+i,1+j,k,2)+3*U
     &  (i,1+j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
      if(maske(1+i,j,-1+k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dwdxb=((U(-2+i,j,-1+k,3)-4*U(-1+i,j,-1+k,3)+3*U(i,
     &  j,-1+k,3))/(2.d0*hx)+(U(-2+i,j,k,3)-4*U(-1+i,j,k,3
     &  )+3*U(i,j,k,3))/(2.d0*hx))/2.d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
      endif
        dwdxt=trandert(i,j,1+k,3,1)
        dudye=trandere(1+i,j,k,1,2)
      if(masks(-1+i,-1+j,k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dudyw=((-3*U(-1+i,j,k,1)+4*U(-1+i,1+j,k,1)-U(-1+i,
     &  2+j,k,1))/(2.d0*hy)+(-3*U(i,j,k,1)+4*U(i,1+j,k,1)-
     &  U(i,2+j,k,1))/(2.d0*hy))/2.d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
      endif
      if(masks(i,-1+j,-1+k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dwdyb=((-3*U(i,j,-1+k,3)+4*U(i,1+j,-1+k,3)-U(i,2+j
     &  ,-1+k,3))/(2.d0*hy)+(-3*U(i,j,k,3)+4*U(i,1+j,k,3)-
     &  U(i,2+j,k,3))/(2.d0*hy))/2.d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
      endif
        dwdyt=trandert(i,j,1+k,3,2)
      if(maskt(-1+i,j,1+k).gt.0.or.maskt(i,j,1+k).gt.0) then 
        dudzw=((U(-1+i,j,-2+k,1)-4*U(-1+i,j,-1+k,1)+3*U(-1
     &  +i,j,k,1))/(2.d0*hz)+(U(i,j,-2+k,1)-4*U(i,j,-1+k,1
     &  )+3*U(i,j,k,1))/(2.d0*hz))/2.d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
      endif
        dudze=trandere(1+i,j,k,1,3)
        dvdzs=tranders(i,-1+j,k,2,3)
      if(maskt(i,j,1+k).gt.0.or.maskt(i,1+j,1+k).gt.0) then 
        dvdzn=((U(i,j,-2+k,2)-4*U(i,j,-1+k,2)+3*U(i,j,k,2)
     &  )/(2.d0*hz)+(U(i,1+j,-2+k,2)-4*U(i,1+j,-1+k,2)+3*U
     &  (i,1+j,k,2))/(2.d0*hz))/2.d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)
      endif

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)



c ::: case 22, top-south-west corner
c ::: Null
           i=lo(1)
           j=lo(2)
           k=hi(3)
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,1+j,k).gt.0) then 
        dvdxn=((-3*U(i,j,k,2)+4*U(1+i,j,k,2)-U(2+i,j,k,2))
     &  /(2.d0*hx)+(-3*U(i,1+j,k,2)+4*U(1+i,1+j,k,2)-U(2+i
     &  ,1+j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
      if(maskw(-1+i,j,-1+k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dwdxb=((-3*U(i,j,-1+k,3)+4*U(1+i,j,-1+k,3)-U(2+i,j
     &  ,-1+k,3))/(2.d0*hx)+(-3*U(i,j,k,3)+4*U(1+i,j,k,3)-
     &  U(2+i,j,k,3))/(2.d0*hx))/2.d0
      else
        dwdxb=(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+i,j,-1+
     &  k,3)+U(1+i,j,k,3))/(4.d0*hx)
      endif
        dwdxt=trandert(i,j,1+k,3,1)
      if(masks(i,-1+j,k).gt.0.or.masks(1+i,-1+j,k).gt.0) then 
        dudye=((-3*U(i,j,k,1)+4*U(i,1+j,k,1)-U(i,2+j,k,1))
     &  /(2.d0*hy)+(-3*U(1+i,j,k,1)+4*U(1+i,1+j,k,1)-U(1+i
     &  ,2+j,k,1))/(2.d0*hy))/2.d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
      endif
        dudyw=tranderw(-1+i,j,k,1,2)
      if(masks(i,-1+j,-1+k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dwdyb=((-3*U(i,j,-1+k,3)+4*U(i,1+j,-1+k,3)-U(i,2+j
     &  ,-1+k,3))/(2.d0*hy)+(-3*U(i,j,k,3)+4*U(i,1+j,k,3)-
     &  U(i,2+j,k,3))/(2.d0*hy))/2.d0
      else
        dwdyb=(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,1+j,-1+
     &  k,3)+U(i,1+j,k,3))/(4.d0*hy)
      endif
        dwdyt=trandert(i,j,1+k,3,2)
        dudzw=tranderw(-1+i,j,k,1,3)
      if(maskt(i,j,1+k).gt.0.or.maskt(1+i,j,1+k).gt.0) then 
        dudze=((U(i,j,-2+k,1)-4*U(i,j,-1+k,1)+3*U(i,j,k,1)
     &  )/(2.d0*hz)+(U(1+i,j,-2+k,1)-4*U(1+i,j,-1+k,1)+3*U
     &  (1+i,j,k,1))/(2.d0*hz))/2.d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
      endif
        dvdzs=tranders(i,-1+j,k,2,3)
      if(maskt(i,j,1+k).gt.0.or.maskt(i,1+j,1+k).gt.0) then 
        dvdzn=((U(i,j,-2+k,2)-4*U(i,j,-1+k,2)+3*U(i,j,k,2)
     &  )/(2.d0*hz)+(U(i,1+j,-2+k,2)-4*U(i,1+j,-1+k,2)+3*U
     &  (i,1+j,k,2))/(2.d0*hz))/2.d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)
      endif

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)



c ::: case 23, bottom-north-east corner
c ::: Null
           i=hi(1)
           j=hi(2)
           k=lo(3)
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if(maske(1+i,-1+j,k).gt.0.or.maske(1+i,j,k).gt.0) then 
        dvdxs=((U(-2+i,-1+j,k,2)-4*U(-1+i,-1+j,k,2)+3*U(i,
     &  -1+j,k,2))/(2.d0*hx)+(U(-2+i,j,k,2)-4*U(-1+i,j,k,2
     &  )+3*U(i,j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
      endif
        dwdxb=tranderb(i,j,-1+k,3,1)
      if(maske(1+i,j,k).gt.0.or.maske(1+i,j,1+k).gt.0) then 
        dwdxt=((U(-2+i,j,k,3)-4*U(-1+i,j,k,3)+3*U(i,j,k,3)
     &  )/(2.d0*hx)+(U(-2+i,j,1+k,3)-4*U(-1+i,j,1+k,3)+3*U
     &  (i,j,1+k,3))/(2.d0*hx))/2.d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
      endif
        dudye=trandere(1+i,j,k,1,2)
      if(maskn(-1+i,1+j,k).gt.0.or.maskn(i,1+j,k).gt.0) then 
        dudyw=((U(-1+i,-2+j,k,1)-4*U(-1+i,-1+j,k,1)+3*U(-1
     &  +i,j,k,1))/(2.d0*hy)+(U(i,-2+j,k,1)-4*U(i,-1+j,k,1
     &  )+3*U(i,j,k,1))/(2.d0*hy))/2.d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
      endif
        dwdyb=tranderb(i,j,-1+k,3,2)
      if(maskn(i,1+j,k).gt.0.or.maskn(i,1+j,1+k).gt.0) then 
        dwdyt=((U(i,-2+j,k,3)-4*U(i,-1+j,k,3)+3*U(i,j,k,3)
     &  )/(2.d0*hy)+(U(i,-2+j,1+k,3)-4*U(i,-1+j,1+k,3)+3*U
     &  (i,j,1+k,3))/(2.d0*hy))/2.d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
      endif
      if(maskb(-1+i,j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dudzw=((-3*U(-1+i,j,k,1)+4*U(-1+i,j,1+k,1)-U(-1+i,
     &  j,2+k,1))/(2.d0*hz)+(-3*U(i,j,k,1)+4*U(i,j,1+k,1)-
     &  U(i,j,2+k,1))/(2.d0*hz))/2.d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
      endif
        dudze=trandere(1+i,j,k,1,3)
      if(maskb(i,-1+j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dvdzs=((-3*U(i,-1+j,k,2)+4*U(i,-1+j,1+k,2)-U(i,-1+
     &  j,2+k,2))/(2.d0*hz)+(-3*U(i,j,k,2)+4*U(i,j,1+k,2)-
     &  U(i,j,2+k,2))/(2.d0*hz))/2.d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)



c ::: case 24, bottom-north-west corner
c ::: Null
           i=lo(1)
           j=hi(2)
           k=lo(3)
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if(maskw(-1+i,-1+j,k).gt.0.or.maskw(-1+i,j,k).gt.0) then 
        dvdxs=((-3*U(i,-1+j,k,2)+4*U(1+i,-1+j,k,2)-U(2+i,-
     &  1+j,k,2))/(2.d0*hx)+(-3*U(i,j,k,2)+4*U(1+i,j,k,2)-
     &  U(2+i,j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxs=(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+i,-1+j,
     &  k,2)+U(1+i,j,k,2))/(4.d0*hx)
      endif
        dwdxb=tranderb(i,j,-1+k,3,1)
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,j,1+k).gt.0) then 
        dwdxt=((-3*U(i,j,k,3)+4*U(1+i,j,k,3)-U(2+i,j,k,3))
     &  /(2.d0*hx)+(-3*U(i,j,1+k,3)+4*U(1+i,j,1+k,3)-U(2+i
     &  ,j,1+k,3))/(2.d0*hx))/2.d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
      endif
      if(maskn(i,1+j,k).gt.0.or.maskn(1+i,1+j,k).gt.0) then 
        dudye=((U(i,-2+j,k,1)-4*U(i,-1+j,k,1)+3*U(i,j,k,1)
     &  )/(2.d0*hy)+(U(1+i,-2+j,k,1)-4*U(1+i,-1+j,k,1)+3*U
     &  (1+i,j,k,1))/(2.d0*hy))/2.d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
      endif
        dudyw=tranderw(-1+i,j,k,1,2)
        dwdyb=tranderb(i,j,-1+k,3,2)
      if(maskn(i,1+j,k).gt.0.or.maskn(i,1+j,1+k).gt.0) then 
        dwdyt=((U(i,-2+j,k,3)-4*U(i,-1+j,k,3)+3*U(i,j,k,3)
     &  )/(2.d0*hy)+(U(i,-2+j,1+k,3)-4*U(i,-1+j,1+k,3)+3*U
     &  (i,j,1+k,3))/(2.d0*hy))/2.d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
      endif
        dudzw=tranderw(-1+i,j,k,1,3)
      if(maskb(i,j,-1+k).gt.0.or.maskb(1+i,j,-1+k).gt.0) then 
        dudze=((-3*U(i,j,k,1)+4*U(i,j,1+k,1)-U(i,j,2+k,1))
     &  /(2.d0*hz)+(-3*U(1+i,j,k,1)+4*U(1+i,j,1+k,1)-U(1+i
     &  ,j,2+k,1))/(2.d0*hz))/2.d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
      endif
      if(maskb(i,-1+j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dvdzs=((-3*U(i,-1+j,k,2)+4*U(i,-1+j,1+k,2)-U(i,-1+
     &  j,2+k,2))/(2.d0*hz)+(-3*U(i,j,k,2)+4*U(i,j,1+k,2)-
     &  U(i,j,2+k,2))/(2.d0*hz))/2.d0
      else
        dvdzs=(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(i,j,-1+
     &  k,2)+U(i,j,1+k,2))/(4.d0*hz)
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)



c ::: case 25, bottom-south-east corner
c ::: Null
           i=hi(1)
           j=lo(2)
           k=lo(3)
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if(maske(1+i,j,k).gt.0.or.maske(1+i,1+j,k).gt.0) then 
        dvdxn=((U(-2+i,j,k,2)-4*U(-1+i,j,k,2)+3*U(i,j,k,2)
     &  )/(2.d0*hx)+(U(-2+i,1+j,k,2)-4*U(-1+i,1+j,k,2)+3*U
     &  (i,1+j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxb=tranderb(i,j,-1+k,3,1)
      if(maske(1+i,j,k).gt.0.or.maske(1+i,j,1+k).gt.0) then 
        dwdxt=((U(-2+i,j,k,3)-4*U(-1+i,j,k,3)+3*U(i,j,k,3)
     &  )/(2.d0*hx)+(U(-2+i,j,1+k,3)-4*U(-1+i,j,1+k,3)+3*U
     &  (i,j,1+k,3))/(2.d0*hx))/2.d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
      endif
        dudye=trandere(1+i,j,k,1,2)
      if(masks(-1+i,-1+j,k).gt.0.or.masks(i,-1+j,k).gt.0) then 
        dudyw=((-3*U(-1+i,j,k,1)+4*U(-1+i,1+j,k,1)-U(-1+i,
     &  2+j,k,1))/(2.d0*hy)+(-3*U(i,j,k,1)+4*U(i,1+j,k,1)-
     &  U(i,2+j,k,1))/(2.d0*hy))/2.d0
      else
        dudyw=(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(i,-1+j,
     &  k,1)+U(i,1+j,k,1))/(4.d0*hy)
      endif
        dwdyb=tranderb(i,j,-1+k,3,2)
      if(masks(i,-1+j,k).gt.0.or.masks(i,-1+j,1+k).gt.0) then 
        dwdyt=((-3*U(i,j,k,3)+4*U(i,1+j,k,3)-U(i,2+j,k,3))
     &  /(2.d0*hy)+(-3*U(i,j,1+k,3)+4*U(i,1+j,1+k,3)-U(i,2
     &  +j,1+k,3))/(2.d0*hy))/2.d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
      endif
      if(maskb(-1+i,j,-1+k).gt.0.or.maskb(i,j,-1+k).gt.0) then 
        dudzw=((-3*U(-1+i,j,k,1)+4*U(-1+i,j,1+k,1)-U(-1+i,
     &  j,2+k,1))/(2.d0*hz)+(-3*U(i,j,k,1)+4*U(i,j,1+k,1)-
     &  U(i,j,2+k,1))/(2.d0*hz))/2.d0
      else
        dudzw=(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(i,j,-1+
     &  k,1)+U(i,j,1+k,1))/(4.d0*hz)
      endif
        dudze=trandere(1+i,j,k,1,3)
        dvdzs=tranders(i,-1+j,k,2,3)
      if(maskb(i,j,-1+k).gt.0.or.maskb(i,1+j,-1+k).gt.0) then 
        dvdzn=((-3*U(i,j,k,2)+4*U(i,j,1+k,2)-U(i,j,2+k,2))
     &  /(2.d0*hz)+(-3*U(i,1+j,k,2)+4*U(i,1+j,1+k,2)-U(i,1
     &  +j,2+k,2))/(2.d0*hz))/2.d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)
      endif

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)



c ::: case 26, bottom-south-west corner
c ::: Null
           i=lo(1)
           j=lo(2)
           k=lo(3)
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,1+j,k).gt.0) then 
        dvdxn=((-3*U(i,j,k,2)+4*U(1+i,j,k,2)-U(2+i,j,k,2))
     &  /(2.d0*hx)+(-3*U(i,1+j,k,2)+4*U(1+i,1+j,k,2)-U(2+i
     &  ,1+j,k,2))/(2.d0*hx))/2.d0
      else
        dvdxn=(-U(-1+i,j,k,2)-U(-1+i,1+j,k,2)+U(1+i,j,k,2)
     &  +U(1+i,1+j,k,2))/(4.d0*hx)
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxb=tranderb(i,j,-1+k,3,1)
      if(maskw(-1+i,j,k).gt.0.or.maskw(-1+i,j,1+k).gt.0) then 
        dwdxt=((-3*U(i,j,k,3)+4*U(1+i,j,k,3)-U(2+i,j,k,3))
     &  /(2.d0*hx)+(-3*U(i,j,1+k,3)+4*U(1+i,j,1+k,3)-U(2+i
     &  ,j,1+k,3))/(2.d0*hx))/2.d0
      else
        dwdxt=(-U(-1+i,j,k,3)-U(-1+i,j,1+k,3)+U(1+i,j,k,3)
     &  +U(1+i,j,1+k,3))/(4.d0*hx)
      endif
      if(masks(i,-1+j,k).gt.0.or.masks(1+i,-1+j,k).gt.0) then 
        dudye=((-3*U(i,j,k,1)+4*U(i,1+j,k,1)-U(i,2+j,k,1))
     &  /(2.d0*hy)+(-3*U(1+i,j,k,1)+4*U(1+i,1+j,k,1)-U(1+i
     &  ,2+j,k,1))/(2.d0*hy))/2.d0
      else
        dudye=(-U(i,-1+j,k,1)+U(i,1+j,k,1)-U(1+i,-1+j,k,1)
     &  +U(1+i,1+j,k,1))/(4.d0*hy)
      endif
        dudyw=tranderw(-1+i,j,k,1,2)
        dwdyb=tranderb(i,j,-1+k,3,2)
      if(masks(i,-1+j,k).gt.0.or.masks(i,-1+j,1+k).gt.0) then 
        dwdyt=((-3*U(i,j,k,3)+4*U(i,1+j,k,3)-U(i,2+j,k,3))
     &  /(2.d0*hy)+(-3*U(i,j,1+k,3)+4*U(i,1+j,1+k,3)-U(i,2
     &  +j,1+k,3))/(2.d0*hy))/2.d0
      else
        dwdyt=(-U(i,-1+j,k,3)-U(i,-1+j,1+k,3)+U(i,1+j,k,3)
     &  +U(i,1+j,1+k,3))/(4.d0*hy)
      endif
        dudzw=tranderw(-1+i,j,k,1,3)
      if(maskb(i,j,-1+k).gt.0.or.maskb(1+i,j,-1+k).gt.0) then 
        dudze=((-3*U(i,j,k,1)+4*U(i,j,1+k,1)-U(i,j,2+k,1))
     &  /(2.d0*hz)+(-3*U(1+i,j,k,1)+4*U(1+i,j,1+k,1)-U(1+i
     &  ,j,2+k,1))/(2.d0*hz))/2.d0
      else
        dudze=(-U(i,j,-1+k,1)+U(i,j,1+k,1)-U(1+i,j,-1+k,1)
     &  +U(1+i,j,1+k,1))/(4.d0*hz)
      endif
        dvdzs=tranders(i,-1+j,k,2,3)
      if(maskb(i,j,-1+k).gt.0.or.maskb(i,1+j,-1+k).gt.0) then 
        dvdzn=((-3*U(i,j,k,2)+4*U(i,j,1+k,2)-U(i,j,2+k,2))
     &  /(2.d0*hz)+(-3*U(i,1+j,k,2)+4*U(i,1+j,1+k,2)-U(i,1
     &  +j,2+k,2))/(2.d0*hz))/2.d0
      else
        dvdzn=(-U(i,j,-1+k,2)+U(i,j,1+k,2)-U(i,1+j,-1+k,2)
     &  +U(i,1+j,1+k,2))/(4.d0*hz)
      endif

c     :::  evaluate expression
               out(i,j,k,1) = -((beta*(hy*hz*(-2*dudxw*muX(i,j,k)+2*dudxe*muX(1+i,j,k))+
     &  hx*hz*(-((dudys+dvdxs)*muY(i,j,k))+(dudyn+dvdxn)*m
     &  uY(i,1+j,k))+hx*hy*(-((dudzb+dwdxb)*muZ(i,j,k))+(d
     &  udzt+dwdxt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,1)
               out(i,j,k,2) = -((beta*(hy*hz*(-((dudyw+dvdxw)*muX(i,j,k))+(dudye+dvdxe)*
     &  muX(1+i,j,k))+hx*hz*(-2*dvdys*muY(i,j,k)+2*dvdyn*m
     &  uY(i,1+j,k))+hx*hy*(-((dvdzb+dwdyb)*muZ(i,j,k))+(d
     &  vdzt+dwdyt)*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,2)
               out(i,j,k,3) = -((beta*(hy*hz*(-((dudzw+dwdxw)*muX(i,j,k))+(dudze+dwdxe)*
     &  muX(1+i,j,k))+hx*hz*(-((dvdzs+dwdys)*muY(i,j,k))+(
     &  dvdzn+dwdyn)*muY(i,1+j,k))+hx*hy*(-2*dwdzb*muZ(i,j
     &  ,k)+2*dwdzt*muZ(i,j,1+k))))/(hx*hy*hz))+alpha*a(i,
     &  j,k)*u(i,j,k,3)






      return
      end




