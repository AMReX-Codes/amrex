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




c -------------------------------------------------------------------
c ... a point of interpretation here: lo and hi are the covered area in
c ... cell centered coordinates.  The dimension of the x and y fluxes have
c ... to be at least one larger in the appropriate direction.
c ... we use the "usual" interpretation of fluxes: cell i has a left flux
c ... numbered i and a right flux numbered i+1

      subroutine FORT_DVFLUX(
     $     u, DIMS(u),
     $     alpha, beta,
     $     a, DIMS(a),
     $     muX, DIMS(muX),
     $     muY, DIMS(muY),
     $     muZ, DIMS(muZ),
     $     xflux,DIMS(xflux),
     $     yflux,DIMS(yflux),
     $     zflux,DIMS(zflux),
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
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(zflux)
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
      REAL_T xflux(DIMV(xflux),3)
      REAL_T yflux(DIMV(yflux),3)
      REAL_T zflux(DIMV(zflux),3)
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

#if 0
      write(6,*)'dvflux'
#else

c ::: cases:
c     ::: 0) interior.  All usual derivative expressions
c     ::: 1-6) face but not on edge.  
c     ::: 7-18) edge but not on corner.
c     ::: 19-26) corners.

c ::: case 0
c ::: Null
c ::: Null
      do k=lo(3)+1,hi(3)
         do j=lo(2)+1,hi(2)
            do i=lo(1)+1,hi(1)
c     :::  get expressions for normal derivatives
              dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

              dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

              dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

              dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

              dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

              dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

              dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

              dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

              dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxs=2.5d-1*(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+
     &  i,-1+j,k,2)+U(1+i,j,k,2))/hx
        dwdxb=2.5d-1*(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+
     &  i,j,-1+k,3)+U(1+i,j,k,3))/hx
        dudyw=2.5d-1*(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(
     &  i,-1+j,k,1)+U(i,1+j,k,1))/hy
        dwdyb=2.5d-1*(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,
     &  1+j,-1+k,3)+U(i,1+j,k,3))/hy
        dudzw=2.5d-1*(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(
     &  i,j,-1+k,1)+U(i,j,1+k,1))/hz
        dvdzs=2.5d-1*(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(
     &  i,j,-1+k,2)+U(i,j,1+k,2))/hz

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

            enddo
         enddo
      enddo


c ::: case 1, top face
c ::: Null
c ::: Null
c ::: Null
c ::: wst
      k=hi(3)
         do j=lo(2)+1,hi(2)
            do i=lo(1)+1,hi(1)
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz

               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dudyw=2.5d-1*(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(
     &  i,-1+j,k,1)+U(i,1+j,k,1))/hy
        dvdxs=2.5d-1*(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+
     &  i,-1+j,k,2)+U(1+i,j,k,2))/hx
        dwdxt=trandert(i,j,1+k,3,1)
        dwdyt=trandert(i,j,1+k,3,2)
      if((maskt(-1+i,j,1+k).gt.0).or.(maskt(i,j,1+k).gt.0)) then 
        dudzw=5.d-1*(5.d-1*(U(-1+i,j,-2+k,1)-4.d0*U(-1+i,j
     &  ,-1+k,1)+3.d0*U(-1+i,j,k,1))/hz+5.d-1*(U(i,j,-2+k,
     &  1)-4.d0*U(i,j,-1+k,1)+3.d0*U(i,j,k,1))/hz)
      else
        dudzw=2.5d-1*(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(
     &  i,j,-1+k,1)+U(i,j,1+k,1))/hz
      endif
      if((maskt(i,-1+j,1+k).gt.0).or.(maskt(i,j,1+k).gt.0)) then 
        dvdzs=5.d-1*(5.d-1*(U(i,-1+j,-2+k,2)-4.d0*U(i,-1+j
     &  ,-1+k,2)+3.d0*U(i,-1+j,k,2))/hz+5.d-1*(U(i,j,-2+k,
     &  2)-4.d0*U(i,j,-1+k,2)+3.d0*U(i,j,k,2))/hz)
      else
        dvdzs=2.5d-1*(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(
     &  i,j,-1+k,2)+U(i,j,1+k,2))/hz
      endif
c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k+1,1)=(dudzt+dwdxt)*muZ(i,j,1+k)
              zflux(i,j,k+1,2)=(dvdzt+dwdyt)*muZ(i,j,1+k)
              zflux(i,j,k+1,3)=2.d0*dwdzt*muZ(i,j,1+k)

            enddo
         enddo

c ::: case 2, bottom face
c ::: Null
c ::: Null
c ::: Null
c ::: wsb
      k=lo(3)
         do j=lo(2)+1,hi(2)
            do i=lo(1)+1,hi(1)
c     :::  get expressions for normal derivatives
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dudyw=2.5d-1*(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(
     &  i,-1+j,k,1)+U(i,1+j,k,1))/hy
        dvdxs=2.5d-1*(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+
     &  i,-1+j,k,2)+U(1+i,j,k,2))/hx
        dwdxb=tranderb(i,j,-1+k,3,1)
        dwdyb=tranderb(i,j,-1+k,3,2)
      if((maskb(-1+i,j,-1+k).gt.0).or.(maskb(i,j,-1+k).gt.0)) then 
        dudzw=5.d-1*(5.d-1*(-3.d0*U(-1+i,j,k,1)+4.d0*U(-1+
     &  i,j,1+k,1)-U(-1+i,j,2+k,1))/hz+5.d-1*(-3.d0*U(i,j,
     &  k,1)+4.d0*U(i,j,1+k,1)-U(i,j,2+k,1))/hz)
      else
        dudzw=2.5d-1*(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(
     &  i,j,-1+k,1)+U(i,j,1+k,1))/hz
      endif
      if((maskb(i,-1+j,-1+k).gt.0).or.(maskb(i,j,-1+k).gt.0)) then 
        dvdzs=5.d-1*(5.d-1*(-3.d0*U(i,-1+j,k,2)+4.d0*U(i,-
     &  1+j,1+k,2)-U(i,-1+j,2+k,2))/hz+5.d-1*(-3.d0*U(i,j,
     &  k,2)+4.d0*U(i,j,1+k,2)-U(i,j,2+k,2))/hz)
      else
        dvdzs=2.5d-1*(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(
     &  i,j,-1+k,2)+U(i,j,1+k,2))/hz
      endif

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

            enddo
         enddo


c ::: case 3, west face
c ::: Null
c ::: Null
c ::: Null
c ::: wsb
      i = lo(1)
      do k=lo(3)+1,hi(3)
         do j=lo(2)+1,hi(2)
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dudyw=tranderw(-1+i,j,k,1,2)
      if((maskw(-1+i,-1+j,k).gt.0).or.(maskw(-1+i,j,k).gt.0)) then 
        dvdxs=5.d-1*(5.d-1*(-3.d0*U(i,-1+j,k,2)+4.d0*U(1+i
     &  ,-1+j,k,2)-U(2+i,-1+j,k,2))/hx+5.d-1*(-3.d0*U(i,j,
     &  k,2)+4.d0*U(1+i,j,k,2)-U(2+i,j,k,2))/hx)
      else
        dvdxs=2.5d-1*(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+
     &  i,-1+j,k,2)+U(1+i,j,k,2))/hx
      endif
      if((maskw(-1+i,j,-1+k).gt.0).or.(maskw(-1+i,j,k).gt.0)) then 
        dwdxb=5.d-1*(5.d-1*(-3.d0*U(i,j,-1+k,3)+4.d0*U(1+i
     &  ,j,-1+k,3)-U(2+i,j,-1+k,3))/hx+5.d-1*(-3.d0*U(i,j,
     &  k,3)+4.d0*U(1+i,j,k,3)-U(2+i,j,k,3))/hx)
      else
        dwdxb=2.5d-1*(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+
     &  i,j,-1+k,3)+U(1+i,j,k,3))/hx
      endif
        dwdyb=2.5d-1*(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,
     &  1+j,-1+k,3)+U(i,1+j,k,3))/hy
        dudzw=tranderw(-1+i,j,k,1,3)
        dvdzs=2.5d-1*(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(
     &  i,j,-1+k,2)+U(i,j,1+k,2))/hz

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

         enddo
      enddo


c ::: case 4, east face
c ::: Null
c ::: Null
c ::: Null
c ::: esb
      i = hi(1)
      do k=lo(3)+1,hi(3)
         do j=lo(2)+1,hi(2)
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dudye=trandere(1+i,j,k,1,2)
      if((maske(1+i,-1+j,k).gt.0).or.(maske(1+i,j,k).gt.0)) then 
        dvdxs=5.d-1*(5.d-1*(U(-2+i,-1+j,k,2)-4.d0*U(-1+i,-
     &  1+j,k,2)+3.d0*U(i,-1+j,k,2))/hx+5.d-1*(U(-2+i,j,k,
     &  2)-4.d0*U(-1+i,j,k,2)+3.d0*U(i,j,k,2))/hx)
      else
        dvdxs=2.5d-1*(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+
     &  i,-1+j,k,2)+U(1+i,j,k,2))/hx
      endif
      if((maske(1+i,j,-1+k).gt.0).or.(maske(1+i,j,k).gt.0)) then 
        dwdxb=5.d-1*(5.d-1*(U(-2+i,j,-1+k,3)-4.d0*U(-1+i,j
     &  ,-1+k,3)+3.d0*U(i,j,-1+k,3))/hx+5.d-1*(U(-2+i,j,k,
     &  3)-4.d0*U(-1+i,j,k,3)+3.d0*U(i,j,k,3))/hx)
      else
        dwdxb=2.5d-1*(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+
     &  i,j,-1+k,3)+U(1+i,j,k,3))/hx
      endif
        dwdyb=2.5d-1*(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,
     &  1+j,-1+k,3)+U(i,1+j,k,3))/hy
        dudze=trandere(1+i,j,k,1,3)
        dvdzs=2.5d-1*(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(
     &  i,j,-1+k,2)+U(i,j,1+k,2))/hz

c     :::  evaluate expression
              xflux(i+1,j,k,1)=2.d0*dudxe*muX(1+i,j,k)
              xflux(i+1,j,k,2)=(dudye+dvdxe)*muX(1+i,j,k)
              xflux(i+1,j,k,3)=(dudze+dwdxe)*muX(1+i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

         enddo
      enddo


c ::: case 5, north face
c ::: Null
c ::: Null
c ::: Null
c ::: wnb
      j = hi(2)
      do k=lo(3)+1,hi(3)
            do i=lo(1)+1,hi(1)
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if((maskn(-1+i,1+j,k).gt.0).or.(maskn(i,1+j,k).gt.0)) then 
        dudyw=5.d-1*(5.d-1*(U(-1+i,-2+j,k,1)-4.d0*U(-1+i,-
     &  1+j,k,1)+3.d0*U(-1+i,j,k,1))/hy+5.d-1*(U(i,-2+j,k,
     &  1)-4.d0*U(i,-1+j,k,1)+3.d0*U(i,j,k,1))/hy)
      else
        dudyw=2.5d-1*(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(
     &  i,-1+j,k,1)+U(i,1+j,k,1))/hy
      endif
        dvdxn=trandern(i,1+j,k,2,1)
        dwdxb=2.5d-1*(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+
     &  i,j,-1+k,3)+U(1+i,j,k,3))/hx
      if((maskn(i,1+j,-1+k).gt.0).or.(maskn(i,1+j,k).gt.0)) then 
        dwdyb=5.d-1*(5.d-1*(U(i,-2+j,-1+k,3)-4.d0*U(i,-1+j
     &  ,-1+k,3)+3.d0*U(i,j,-1+k,3))/hy+5.d-1*(U(i,-2+j,k,
     &  3)-4.d0*U(i,-1+j,k,3)+3.d0*U(i,j,k,3))/hy)
      else
        dwdyb=2.5d-1*(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,
     &  1+j,-1+k,3)+U(i,1+j,k,3))/hy
      endif
        dudzw=2.5d-1*(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(
     &  i,j,-1+k,1)+U(i,j,1+k,1))/hz
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j+1,k,1)=(dudyn+dvdxn)*muY(i,1+j,k)
              yflux(i,j+1,k,2)=2.d0*dvdyn*muY(i,1+j,k)
              yflux(i,j+1,k,3)=(dvdzn+dwdyn)*muY(i,1+j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)
            enddo
      enddo

	
c ::: case 6, south face
c ::: Null
c ::: Null
c ::: Null
c ::: wsb
      j = lo(2)
      do k=lo(3)+1,hi(3)
            do i=lo(1)+1,hi(1)
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if((masks(-1+i,-1+j,k).gt.0).or.(masks(i,-1+j,k).gt.0)) then 
        dudyw=5.d-1*(5.d-1*(-3.d0*U(-1+i,j,k,1)+4.d0*U(-1+
     &  i,1+j,k,1)-U(-1+i,2+j,k,1))/hy+5.d-1*(-3.d0*U(i,j,
     &  k,1)+4.d0*U(i,1+j,k,1)-U(i,2+j,k,1))/hy)
      else
        dudyw=2.5d-1*(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(
     &  i,-1+j,k,1)+U(i,1+j,k,1))/hy
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxb=2.5d-1*(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+
     &  i,j,-1+k,3)+U(1+i,j,k,3))/hx
      if((masks(i,-1+j,-1+k).gt.0).or.(masks(i,-1+j,k).gt.0)) then 
        dwdyb=5.d-1*(5.d-1*(-3.d0*U(i,j,-1+k,3)+4.d0*U(i,1
     &  +j,-1+k,3)-U(i,2+j,-1+k,3))/hy+5.d-1*(-3.d0*U(i,j,
     &  k,3)+4.d0*U(i,1+j,k,3)-U(i,2+j,k,3))/hy)
      else
        dwdyb=2.5d-1*(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,
     &  1+j,-1+k,3)+U(i,1+j,k,3))/hy
      endif
        dudzw=2.5d-1*(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(
     &  i,j,-1+k,1)+U(i,j,1+k,1))/hz
        dvdzs=tranders(i,-1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)
            enddo
      enddo

c ::: case 7, top-north edge 
c ::: Null
c ::: Null
c ::: Null
c ::: tnw
      k=hi(3)
         j=hi(2)
            do i=lo(1)+1,hi(1)-1
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz

               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if((maskn(-1+i,1+j,k).gt.0).or.(maskn(i,1+j,k).gt.0)) then 
        dudyw=5.d-1*(5.d-1*(U(-1+i,-2+j,k,1)-4.d0*U(-1+i,-
     &  1+j,k,1)+3.d0*U(-1+i,j,k,1))/hy+5.d-1*(U(i,-2+j,k,
     &  1)-4.d0*U(i,-1+j,k,1)+3.d0*U(i,j,k,1))/hy)
      else
        dudyw=2.5d-1*(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(
     &  i,-1+j,k,1)+U(i,1+j,k,1))/hy
      endif
        dvdxn=trandern(i,1+j,k,2,1)
        dwdxt=trandert(i,j,1+k,3,1)
        dwdyt=trandert(i,j,1+k,3,2)
      if((maskt(-1+i,j,1+k).gt.0).or.(maskt(i,j,1+k).gt.0)) then 
        dudzw=5.d-1*(5.d-1*(U(-1+i,j,-2+k,1)-4.d0*U(-1+i,j
     &  ,-1+k,1)+3.d0*U(-1+i,j,k,1))/hz+5.d-1*(U(i,j,-2+k,
     &  1)-4.d0*U(i,j,-1+k,1)+3.d0*U(i,j,k,1))/hz)
      else
        dudzw=2.5d-1*(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(
     &  i,j,-1+k,1)+U(i,j,1+k,1))/hz
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j+1,k,1)=(dudyn+dvdxn)*muY(i,1+j,k)
              yflux(i,j+1,k,2)=2.d0*dvdyn*muY(i,1+j,k)
              yflux(i,j+1,k,3)=(dvdzn+dwdyn)*muY(i,1+j,k)

              zflux(i,j,k+1,1)=(dudzt+dwdxt)*muZ(i,j,1+k)
              zflux(i,j,k+1,2)=(dvdzt+dwdyt)*muZ(i,j,1+k)
              zflux(i,j,k+1,3)=2.d0*dwdzt*muZ(i,j,1+k)

            enddo

c ::: case 8, top-south edge 
c ::: Null
c ::: tsw
      k=hi(3)
         j=lo(2)
            do i=lo(1)+1,hi(1)-1
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz

               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if((masks(-1+i,-1+j,k).gt.0).or.(masks(i,-1+j,k).gt.0)) then 
        dudyw=5.d-1*(5.d-1*(-3.d0*U(-1+i,j,k,1)+4.d0*U(-1+
     &  i,1+j,k,1)-U(-1+i,2+j,k,1))/hy+5.d-1*(-3.d0*U(i,j,
     &  k,1)+4.d0*U(i,1+j,k,1)-U(i,2+j,k,1))/hy)
      else
        dudyw=2.5d-1*(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(
     &  i,-1+j,k,1)+U(i,1+j,k,1))/hy
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxt=trandert(i,j,1+k,3,1)
        dwdyt=trandert(i,j,1+k,3,2)
      if((maskt(-1+i,j,1+k).gt.0).or.(maskt(i,j,1+k).gt.0)) then 
        dudzw=5.d-1*(5.d-1*(U(-1+i,j,-2+k,1)-4.d0*U(-1+i,j
     &  ,-1+k,1)+3.d0*U(-1+i,j,k,1))/hz+5.d-1*(U(i,j,-2+k,
     &  1)-4.d0*U(i,j,-1+k,1)+3.d0*U(i,j,k,1))/hz)
      else
        dudzw=2.5d-1*(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(
     &  i,j,-1+k,1)+U(i,j,1+k,1))/hz
      endif
        dvdzs=tranders(i,-1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k+1,1)=(dudzt+dwdxt)*muZ(i,j,1+k)
              zflux(i,j,k+1,2)=(dvdzt+dwdyt)*muZ(i,j,1+k)
              zflux(i,j,k+1,3)=2.d0*dwdzt*muZ(i,j,1+k)
            enddo


c ::: case 9, top-west edge 
c ::: Null
c ::: tsw
      k=hi(3)
      i=lo(1)
         do j=lo(2)+1,hi(2)-1
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz

               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dudyw=tranderw(-1+i,j,k,1,2)
      if((maskw(-1+i,-1+j,k).gt.0).or.(maskw(-1+i,j,k).gt.0)) then 
        dvdxs=5.d-1*(5.d-1*(-3.d0*U(i,-1+j,k,2)+4.d0*U(1+i
     &  ,-1+j,k,2)-U(2+i,-1+j,k,2))/hx+5.d-1*(-3.d0*U(i,j,
     &  k,2)+4.d0*U(1+i,j,k,2)-U(2+i,j,k,2))/hx)
      else
        dvdxs=2.5d-1*(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+
     &  i,-1+j,k,2)+U(1+i,j,k,2))/hx
      endif
        dwdxt=trandert(i,j,1+k,3,1)
        dwdyt=trandert(i,j,1+k,3,2)
        dudzw=tranderw(-1+i,j,k,1,3)
      if((maskt(i,-1+j,1+k).gt.0).or.(maskt(i,j,1+k).gt.0)) then 
        dvdzs=5.d-1*(5.d-1*(U(i,-1+j,-2+k,2)-4.d0*U(i,-1+j
     &  ,-1+k,2)+3.d0*U(i,-1+j,k,2))/hz+5.d-1*(U(i,j,-2+k,
     &  2)-4.d0*U(i,j,-1+k,2)+3.d0*U(i,j,k,2))/hz)
      else
        dvdzs=2.5d-1*(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(
     &  i,j,-1+k,2)+U(i,j,1+k,2))/hz
      endif

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k+1,1)=(dudzt+dwdxt)*muZ(i,j,1+k)
              zflux(i,j,k+1,2)=(dvdzt+dwdyt)*muZ(i,j,1+k)
              zflux(i,j,k+1,3)=2.d0*dwdzt*muZ(i,j,1+k)

         enddo

c ::: case 10, top-east edge 
c ::: Null
c ::: tes
      k=hi(3)
      i=hi(1)
         do j=lo(2)+1,hi(2)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz

               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dudye=trandere(1+i,j,k,1,2)
      if((maske(1+i,-1+j,k).gt.0).or.(maske(1+i,j,k).gt.0)) then 
        dvdxs=5.d-1*(5.d-1*(U(-2+i,-1+j,k,2)-4.d0*U(-1+i,-
     &  1+j,k,2)+3.d0*U(i,-1+j,k,2))/hx+5.d-1*(U(-2+i,j,k,
     &  2)-4.d0*U(-1+i,j,k,2)+3.d0*U(i,j,k,2))/hx)
      else
        dvdxs=2.5d-1*(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+
     &  i,-1+j,k,2)+U(1+i,j,k,2))/hx
      endif
        dwdxt=trandert(i,j,1+k,3,1)
        dwdyt=trandert(i,j,1+k,3,2)
        dudze=trandere(1+i,j,k,1,3)
      if((maskt(i,-1+j,1+k).gt.0).or.(maskt(i,j,1+k).gt.0)) then 
        dvdzs=5.d-1*(5.d-1*(U(i,-1+j,-2+k,2)-4.d0*U(i,-1+j
     &  ,-1+k,2)+3.d0*U(i,-1+j,k,2))/hz+5.d-1*(U(i,j,-2+k,
     &  2)-4.d0*U(i,j,-1+k,2)+3.d0*U(i,j,k,2))/hz)
      else
        dvdzs=2.5d-1*(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(
     &  i,j,-1+k,2)+U(i,j,1+k,2))/hz
      endif

c     :::  evaluate expression
              xflux(i+1,j,k,1)=2.d0*dudxe*muX(1+i,j,k)
              xflux(i+1,j,k,2)=(dudye+dvdxe)*muX(1+i,j,k)
              xflux(i+1,j,k,3)=(dudze+dwdxe)*muX(1+i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k+1,1)=(dudzt+dwdxt)*muZ(i,j,1+k)
              zflux(i,j,k+1,2)=(dvdzt+dwdyt)*muZ(i,j,1+k)
              zflux(i,j,k+1,3)=2.d0*dwdzt*muZ(i,j,1+k)

         enddo
c ::: case 11, bottom-north edge 
c ::: Null
c ::: bnw
      k=lo(3)
      j=hi(2)
            do i=lo(1)+1,hi(1)-1
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if((maskn(-1+i,1+j,k).gt.0).or.(maskn(i,1+j,k).gt.0)) then 
        dudyw=5.d-1*(5.d-1*(U(-1+i,-2+j,k,1)-4.d0*U(-1+i,-
     &  1+j,k,1)+3.d0*U(-1+i,j,k,1))/hy+5.d-1*(U(i,-2+j,k,
     &  1)-4.d0*U(i,-1+j,k,1)+3.d0*U(i,j,k,1))/hy)
      else
        dudyw=2.5d-1*(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(
     &  i,-1+j,k,1)+U(i,1+j,k,1))/hy
      endif
        dvdxn=trandern(i,1+j,k,2,1)
        dwdxb=tranderb(i,j,-1+k,3,1)
        dwdyb=tranderb(i,j,-1+k,3,2)
      if((maskb(-1+i,j,-1+k).gt.0).or.(maskb(i,j,-1+k).gt.0)) then 
        dudzw=5.d-1*(5.d-1*(-3.d0*U(-1+i,j,k,1)+4.d0*U(-1+
     &  i,j,1+k,1)-U(-1+i,j,2+k,1))/hz+5.d-1*(-3.d0*U(i,j,
     &  k,1)+4.d0*U(i,j,1+k,1)-U(i,j,2+k,1))/hz)
      else
        dudzw=2.5d-1*(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(
     &  i,j,-1+k,1)+U(i,j,1+k,1))/hz
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j+1,k,1)=(dudyn+dvdxn)*muY(i,1+j,k)
              yflux(i,j+1,k,2)=2.d0*dvdyn*muY(i,1+j,k)
              yflux(i,j+1,k,3)=(dvdzn+dwdyn)*muY(i,1+j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

            enddo

c ::: case 12, bottom-south edge 
c ::: Null
c ::: bsw
      k=lo(3)
      j=lo(2)
            do i=lo(1)+1,hi(1)-1
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if((masks(-1+i,-1+j,k).gt.0).or.(masks(i,-1+j,k).gt.0)) then 
        dudyw=5.d-1*(5.d-1*(-3.d0*U(-1+i,j,k,1)+4.d0*U(-1+
     &  i,1+j,k,1)-U(-1+i,2+j,k,1))/hy+5.d-1*(-3.d0*U(i,j,
     &  k,1)+4.d0*U(i,1+j,k,1)-U(i,2+j,k,1))/hy)
      else
        dudyw=2.5d-1*(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(
     &  i,-1+j,k,1)+U(i,1+j,k,1))/hy
      endif
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxb=tranderb(i,j,-1+k,3,1)
        dwdyb=tranderb(i,j,-1+k,3,2)
      if((maskb(-1+i,j,-1+k).gt.0).or.(maskb(i,j,-1+k).gt.0)) then 
        dudzw=5.d-1*(5.d-1*(-3.d0*U(-1+i,j,k,1)+4.d0*U(-1+
     &  i,j,1+k,1)-U(-1+i,j,2+k,1))/hz+5.d-1*(-3.d0*U(i,j,
     &  k,1)+4.d0*U(i,j,1+k,1)-U(i,j,2+k,1))/hz)
      else
        dudzw=2.5d-1*(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(
     &  i,j,-1+k,1)+U(i,j,1+k,1))/hz
      endif
        dvdzs=tranders(i,-1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

            enddo


c ::: case 13, bottom-west edge 
c ::: Null
c ::: bws
      k=lo(3)
      i=lo(1)
      do j=lo(2)+1,hi(2)-1
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dudyw=tranderw(-1+i,j,k,1,2)
      if((maskw(-1+i,-1+j,k).gt.0).or.(maskw(-1+i,j,k).gt.0)) then 
        dvdxs=5.d-1*(5.d-1*(-3.d0*U(i,-1+j,k,2)+4.d0*U(1+i
     &  ,-1+j,k,2)-U(2+i,-1+j,k,2))/hx+5.d-1*(-3.d0*U(i,j,
     &  k,2)+4.d0*U(1+i,j,k,2)-U(2+i,j,k,2))/hx)
      else
        dvdxs=2.5d-1*(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+
     &  i,-1+j,k,2)+U(1+i,j,k,2))/hx
      endif
        dwdxb=tranderb(i,j,-1+k,3,1)
        dwdyb=tranderb(i,j,-1+k,3,2)
        dudzw=tranderw(-1+i,j,k,1,3)
      if((maskb(i,-1+j,-1+k).gt.0).or.(maskb(i,j,-1+k).gt.0)) then 
        dvdzs=5.d-1*(5.d-1*(-3.d0*U(i,-1+j,k,2)+4.d0*U(i,-
     &  1+j,1+k,2)-U(i,-1+j,2+k,2))/hz+5.d-1*(-3.d0*U(i,j,
     &  k,2)+4.d0*U(i,j,1+k,2)-U(i,j,2+k,2))/hz)
      else
        dvdzs=2.5d-1*(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(
     &  i,j,-1+k,2)+U(i,j,1+k,2))/hz
      endif

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

         enddo



c ::: case 14, bottom-east edge 
c ::: Null
c ::: bse
      k=lo(3)
      i=hi(1)
         do j=lo(2)+1,hi(2)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
      if((maske(1+i,-1+j,k).gt.0).or.(maske(1+i,j,k).gt.0)) then 
        dvdxs=5.d-1*(5.d-1*(U(-2+i,-1+j,k,2)-4.d0*U(-1+i,-
     &  1+j,k,2)+3.d0*U(i,-1+j,k,2))/hx+5.d-1*(U(-2+i,j,k,
     &  2)-4.d0*U(-1+i,j,k,2)+3.d0*U(i,j,k,2))/hx)
      else
        dvdxs=2.5d-1*(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+
     &  i,-1+j,k,2)+U(1+i,j,k,2))/hx
      endif
        dwdxb=tranderb(i,j,-1+k,3,1)
        dudye=trandere(1+i,j,k,1,2)
        dwdyb=tranderb(i,j,-1+k,3,2)
        dudze=trandere(1+i,j,k,1,3)
      if((maskb(i,-1+j,-1+k).gt.0).or.(maskb(i,j,-1+k).gt.0)) then 
        dvdzs=5.d-1*(5.d-1*(-3.d0*U(i,-1+j,k,2)+4.d0*U(i,-
     &  1+j,1+k,2)-U(i,-1+j,2+k,2))/hz+5.d-1*(-3.d0*U(i,j,
     &  k,2)+4.d0*U(i,j,1+k,2)-U(i,j,2+k,2))/hz)
      else
        dvdzs=2.5d-1*(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(
     &  i,j,-1+k,2)+U(i,j,1+k,2))/hz
      endif

c     :::  evaluate expression
              xflux(i+1,j,k,1)=2.d0*dudxe*muX(1+i,j,k)
              xflux(i+1,j,k,2)=(dudye+dvdxe)*muX(1+i,j,k)
              xflux(i+1,j,k,3)=(dudze+dwdxe)*muX(1+i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

         enddo



c ::: case 15, east-north edge
c ::: Null
c ::: neb
      j=hi(2)
      i=hi(1)
      do k=lo(3)+1,hi(3)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if((maske(1+i,j,-1+k).gt.0).or.(maske(1+i,j,k).gt.0)) then 
        dwdxb=5.d-1*(5.d-1*(U(-2+i,j,-1+k,3)-4.d0*U(-1+i,j
     &  ,-1+k,3)+3.d0*U(i,j,-1+k,3))/hx+5.d-1*(U(-2+i,j,k,
     &  3)-4.d0*U(-1+i,j,k,3)+3.d0*U(i,j,k,3))/hx)
      else
        dwdxb=2.5d-1*(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+
     &  i,j,-1+k,3)+U(1+i,j,k,3))/hx
      endif
        dudye=trandere(1+i,j,k,1,2)
      if((maskn(i,1+j,-1+k).gt.0).or.(maskn(i,1+j,k).gt.0)) then 
        dwdyb=5.d-1*(5.d-1*(U(i,-2+j,-1+k,3)-4.d0*U(i,-1+j
     &  ,-1+k,3)+3.d0*U(i,j,-1+k,3))/hy+5.d-1*(U(i,-2+j,k,
     &  3)-4.d0*U(i,-1+j,k,3)+3.d0*U(i,j,k,3))/hy)
      else
        dwdyb=2.5d-1*(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,
     &  1+j,-1+k,3)+U(i,1+j,k,3))/hy
      endif
        dudze=trandere(1+i,j,k,1,3)
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
              xflux(i+1,j,k,1)=2.d0*dudxe*muX(1+i,j,k)
              xflux(i+1,j,k,2)=(dudye+dvdxe)*muX(1+i,j,k)
              xflux(i+1,j,k,3)=(dudze+dwdxe)*muX(1+i,j,k)
              
              yflux(i,j+1,k,1)=(dudyn+dvdxn)*muY(i,1+j,k)
              yflux(i,j+1,k,2)=2.d0*dvdyn*muY(i,1+j,k)
              yflux(i,j+1,k,3)=(dvdzn+dwdyn)*muY(i,1+j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

      enddo


c ::: case 16, east-south edge
c ::: Null
c ::: esb
      j=lo(2)
      i=hi(1)
      do k=lo(3)+1,hi(3)-1
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxs=tranders(i,-1+j,k,2,1)
      if((maske(1+i,j,-1+k).gt.0).or.(maske(1+i,j,k).gt.0)) then 
        dwdxb=5.d-1*(5.d-1*(U(-2+i,j,-1+k,3)-4.d0*U(-1+i,j
     &  ,-1+k,3)+3.d0*U(i,j,-1+k,3))/hx+5.d-1*(U(-2+i,j,k,
     &  3)-4.d0*U(-1+i,j,k,3)+3.d0*U(i,j,k,3))/hx)
      else
        dwdxb=2.5d-1*(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+
     &  i,j,-1+k,3)+U(1+i,j,k,3))/hx
      endif
        dudye=trandere(1+i,j,k,1,2)
      if((masks(i,-1+j,-1+k).gt.0).or.(masks(i,-1+j,k).gt.0)) then 
        dwdyb=5.d-1*(5.d-1*(-3.d0*U(i,j,-1+k,3)+4.d0*U(i,1
     &  +j,-1+k,3)-U(i,2+j,-1+k,3))/hy+5.d-1*(-3.d0*U(i,j,
     &  k,3)+4.d0*U(i,1+j,k,3)-U(i,2+j,k,3))/hy)
      else
        dwdyb=2.5d-1*(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,
     &  1+j,-1+k,3)+U(i,1+j,k,3))/hy
      endif
        dudze=trandere(1+i,j,k,1,3)
        dvdzs=tranders(i,-1+j,k,2,3)

c     :::  evaluate expression
              xflux(i+1,j,k,1)=2.d0*dudxe*muX(1+i,j,k)
              xflux(i+1,j,k,2)=(dudye+dvdxe)*muX(1+i,j,k)
              xflux(i+1,j,k,3)=(dudze+dwdxe)*muX(1+i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

      enddo

c ::: case 17, west-north edge
c ::: Null
c ::: nwb
      j=hi(2)
      i=lo(1)
      do k=lo(3)+1,hi(3)-1
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if((maskw(-1+i,j,-1+k).gt.0).or.(maskw(-1+i,j,k).gt.0)) then 
        dwdxb=5.d-1*(5.d-1*(-3.d0*U(i,j,-1+k,3)+4.d0*U(1+i
     &  ,j,-1+k,3)-U(2+i,j,-1+k,3))/hx+5.d-1*(-3.d0*U(i,j,
     &  k,3)+4.d0*U(1+i,j,k,3)-U(2+i,j,k,3))/hx)
      else
        dwdxb=2.5d-1*(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+
     &  i,j,-1+k,3)+U(1+i,j,k,3))/hx
      endif
        dudyw=tranderw(-1+i,j,k,1,2)
      if((maskn(i,1+j,-1+k).gt.0).or.(maskn(i,1+j,k).gt.0)) then 
        dwdyb=5.d-1*(5.d-1*(U(i,-2+j,-1+k,3)-4.d0*U(i,-1+j
     &  ,-1+k,3)+3.d0*U(i,j,-1+k,3))/hy+5.d-1*(U(i,-2+j,k,
     &  3)-4.d0*U(i,-1+j,k,3)+3.d0*U(i,j,k,3))/hy)
      else
        dwdyb=2.5d-1*(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,
     &  1+j,-1+k,3)+U(i,1+j,k,3))/hy
      endif
        dudzw=tranderw(-1+i,j,k,1,3)
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j+1,k,1)=(dudyn+dvdxn)*muY(i,1+j,k)
              yflux(i,j+1,k,2)=2.d0*dvdyn*muY(i,1+j,k)
              yflux(i,j+1,k,3)=(dvdzn+dwdyn)*muY(i,1+j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)
      enddo


c ::: case 18, west-south edge
c ::: Null
c ::: swb
      j=lo(2)
      i=lo(1)
      do k=lo(3)+1,hi(3)-1
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxs=tranders(i,-1+j,k,2,1)
      if((maskw(-1+i,j,-1+k).gt.0).or.(maskw(-1+i,j,k).gt.0)) then 
        dwdxb=5.d-1*(5.d-1*(-3.d0*U(i,j,-1+k,3)+4.d0*U(1+i
     &  ,j,-1+k,3)-U(2+i,j,-1+k,3))/hx+5.d-1*(-3.d0*U(i,j,
     &  k,3)+4.d0*U(1+i,j,k,3)-U(2+i,j,k,3))/hx)
      else
        dwdxb=2.5d-1*(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+
     &  i,j,-1+k,3)+U(1+i,j,k,3))/hx
      endif
        dudyw=tranderw(-1+i,j,k,1,2)
      if((masks(i,-1+j,-1+k).gt.0).or.(masks(i,-1+j,k).gt.0)) then 
        dwdyb=5.d-1*(5.d-1*(-3.d0*U(i,j,-1+k,3)+4.d0*U(i,1
     &  +j,-1+k,3)-U(i,2+j,-1+k,3))/hy+5.d-1*(-3.d0*U(i,j,
     &  k,3)+4.d0*U(i,1+j,k,3)-U(i,2+j,k,3))/hy)
      else
        dwdyb=2.5d-1*(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,
     &  1+j,-1+k,3)+U(i,1+j,k,3))/hy
      endif
        dudzw=tranderw(-1+i,j,k,1,3)
        dvdzs=tranders(i,-1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)
      enddo


c ::: case 19, top-north-east corner
c ::: Null
c ::: tbnsew
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
      if((maske(1+i,-1+j,k).gt.0).or.(maske(1+i,j,k).gt.0)) then 
        dvdxs=5.d-1*(5.d-1*(U(-2+i,-1+j,k,2)-4.d0*U(-1+i,-
     &  1+j,k,2)+3.d0*U(i,-1+j,k,2))/hx+5.d-1*(U(-2+i,j,k,
     &  2)-4.d0*U(-1+i,j,k,2)+3.d0*U(i,j,k,2))/hx)
      else
        dvdxs=2.5d-1*(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+
     &  i,-1+j,k,2)+U(1+i,j,k,2))/hx
      endif
      if((maske(1+i,j,-1+k).gt.0).or.(maske(1+i,j,k).gt.0)) then 
        dwdxb=5.d-1*(5.d-1*(U(-2+i,j,-1+k,3)-4.d0*U(-1+i,j
     &  ,-1+k,3)+3.d0*U(i,j,-1+k,3))/hx+5.d-1*(U(-2+i,j,k,
     &  3)-4.d0*U(-1+i,j,k,3)+3.d0*U(i,j,k,3))/hx)
      else
        dwdxb=2.5d-1*(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+
     &  i,j,-1+k,3)+U(1+i,j,k,3))/hx
      endif
        dwdxt=trandert(i,j,1+k,3,1)
        dudye=trandere(1+i,j,k,1,2)
      if((maskn(-1+i,1+j,k).gt.0).or.(maskn(i,1+j,k).gt.0)) then 
        dudyw=5.d-1*(5.d-1*(U(-1+i,-2+j,k,1)-4.d0*U(-1+i,-
     &  1+j,k,1)+3.d0*U(-1+i,j,k,1))/hy+5.d-1*(U(i,-2+j,k,
     &  1)-4.d0*U(i,-1+j,k,1)+3.d0*U(i,j,k,1))/hy)
      else
        dudyw=2.5d-1*(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(
     &  i,-1+j,k,1)+U(i,1+j,k,1))/hy
      endif
      if((maskn(i,1+j,-1+k).gt.0).or.(maskn(i,1+j,k).gt.0)) then 
        dwdyb=5.d-1*(5.d-1*(U(i,-2+j,-1+k,3)-4.d0*U(i,-1+j
     &  ,-1+k,3)+3.d0*U(i,j,-1+k,3))/hy+5.d-1*(U(i,-2+j,k,
     &  3)-4.d0*U(i,-1+j,k,3)+3.d0*U(i,j,k,3))/hy)
      else
        dwdyb=2.5d-1*(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,
     &  1+j,-1+k,3)+U(i,1+j,k,3))/hy
      endif
        dwdyt=trandert(i,j,1+k,3,2)
      if((maskt(-1+i,j,1+k).gt.0).or.(maskt(i,j,1+k).gt.0)) then 
        dudzw=5.d-1*(5.d-1*(U(-1+i,j,-2+k,1)-4.d0*U(-1+i,j
     &  ,-1+k,1)+3.d0*U(-1+i,j,k,1))/hz+5.d-1*(U(i,j,-2+k,
     &  1)-4.d0*U(i,j,-1+k,1)+3.d0*U(i,j,k,1))/hz)
      else
        dudzw=2.5d-1*(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(
     &  i,j,-1+k,1)+U(i,j,1+k,1))/hz
      endif
        dudze=trandere(1+i,j,k,1,3)
      if((maskt(i,-1+j,1+k).gt.0).or.(maskt(i,j,1+k).gt.0)) then 
        dvdzs=5.d-1*(5.d-1*(U(i,-1+j,-2+k,2)-4.d0*U(i,-1+j
     &  ,-1+k,2)+3.d0*U(i,-1+j,k,2))/hz+5.d-1*(U(i,j,-2+k,
     &  2)-4.d0*U(i,j,-1+k,2)+3.d0*U(i,j,k,2))/hz)
      else
        dvdzs=2.5d-1*(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(
     &  i,j,-1+k,2)+U(i,j,1+k,2))/hz
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              xflux(i+1,j,k,1)=2.d0*dudxe*muX(1+i,j,k)
              xflux(i+1,j,k,2)=(dudye+dvdxe)*muX(1+i,j,k)
              xflux(i+1,j,k,3)=(dudze+dwdxe)*muX(1+i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              yflux(i,j+1,k,1)=(dudyn+dvdxn)*muY(i,1+j,k)
              yflux(i,j+1,k,2)=2.d0*dvdyn*muY(i,1+j,k)
              yflux(i,j+1,k,3)=(dvdzn+dwdyn)*muY(i,1+j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

              zflux(i,j,k+1,1)=(dudzt+dwdxt)*muZ(i,j,1+k)
              zflux(i,j,k+1,2)=(dvdzt+dwdyt)*muZ(i,j,1+k)
              zflux(i,j,k+1,3)=2.d0*dwdzt*muZ(i,j,1+k)



c ::: case 20, top-north-west corner
c ::: Null
c ::: tbnsw
           i=lo(1)
           j=hi(2)
           k=hi(3)
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if((maskw(-1+i,-1+j,k).gt.0).or.(maskw(-1+i,j,k).gt.0)) then 
        dvdxs=5.d-1*(5.d-1*(-3.d0*U(i,-1+j,k,2)+4.d0*U(1+i
     &  ,-1+j,k,2)-U(2+i,-1+j,k,2))/hx+5.d-1*(-3.d0*U(i,j,
     &  k,2)+4.d0*U(1+i,j,k,2)-U(2+i,j,k,2))/hx)
      else
        dvdxs=2.5d-1*(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+
     &  i,-1+j,k,2)+U(1+i,j,k,2))/hx
      endif
      if((maskw(-1+i,j,-1+k).gt.0).or.(maskw(-1+i,j,k).gt.0)) then 
        dwdxb=5.d-1*(5.d-1*(-3.d0*U(i,j,-1+k,3)+4.d0*U(1+i
     &  ,j,-1+k,3)-U(2+i,j,-1+k,3))/hx+5.d-1*(-3.d0*U(i,j,
     &  k,3)+4.d0*U(1+i,j,k,3)-U(2+i,j,k,3))/hx)
      else
        dwdxb=2.5d-1*(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+
     &  i,j,-1+k,3)+U(1+i,j,k,3))/hx
      endif
        dwdxt=trandert(i,j,1+k,3,1)
        dudyw=tranderw(-1+i,j,k,1,2)
      if((maskn(i,1+j,-1+k).gt.0).or.(maskn(i,1+j,k).gt.0)) then 
        dwdyb=5.d-1*(5.d-1*(U(i,-2+j,-1+k,3)-4.d0*U(i,-1+j
     &  ,-1+k,3)+3.d0*U(i,j,-1+k,3))/hy+5.d-1*(U(i,-2+j,k,
     &  3)-4.d0*U(i,-1+j,k,3)+3.d0*U(i,j,k,3))/hy)
      else
        dwdyb=2.5d-1*(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,
     &  1+j,-1+k,3)+U(i,1+j,k,3))/hy
      endif
        dwdyt=trandert(i,j,1+k,3,2)
        dudzw=tranderw(-1+i,j,k,1,3)
      if((maskt(i,-1+j,1+k).gt.0).or.(maskt(i,j,1+k).gt.0)) then 
        dvdzs=5.d-1*(5.d-1*(U(i,-1+j,-2+k,2)-4.d0*U(i,-1+j
     &  ,-1+k,2)+3.d0*U(i,-1+j,k,2))/hz+5.d-1*(U(i,j,-2+k,
     &  2)-4.d0*U(i,j,-1+k,2)+3.d0*U(i,j,k,2))/hz)
      else
        dvdzs=2.5d-1*(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(
     &  i,j,-1+k,2)+U(i,j,1+k,2))/hz
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              yflux(i,j+1,k,1)=(dudyn+dvdxn)*muY(i,1+j,k)
              yflux(i,j+1,k,2)=2.d0*dvdyn*muY(i,1+j,k)
              yflux(i,j+1,k,3)=(dvdzn+dwdyn)*muY(i,1+j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

              zflux(i,j,k+1,1)=(dudzt+dwdxt)*muZ(i,j,1+k)
              zflux(i,j,k+1,2)=(dvdzt+dwdyt)*muZ(i,j,1+k)
              zflux(i,j,k+1,3)=2.d0*dwdzt*muZ(i,j,1+k)


c ::: case 21, top-south-east corner
c ::: Null
c ::: tbsew
           i=hi(1)
           j=lo(2)
           k=hi(3)
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

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

c     ::: get expressions for tangential derivatives
        dvdxs=tranders(i,-1+j,k,2,1)
      if((maske(1+i,j,-1+k).gt.0).or.(maske(1+i,j,k).gt.0)) then 
        dwdxb=5.d-1*(5.d-1*(U(-2+i,j,-1+k,3)-4.d0*U(-1+i,j
     &  ,-1+k,3)+3.d0*U(i,j,-1+k,3))/hx+5.d-1*(U(-2+i,j,k,
     &  3)-4.d0*U(-1+i,j,k,3)+3.d0*U(i,j,k,3))/hx)
      else
        dwdxb=2.5d-1*(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+
     &  i,j,-1+k,3)+U(1+i,j,k,3))/hx
      endif
        dwdxt=trandert(i,j,1+k,3,1)
        dudye=trandere(1+i,j,k,1,2)
      if((masks(-1+i,-1+j,k).gt.0).or.(masks(i,-1+j,k).gt.0)) then 
        dudyw=5.d-1*(5.d-1*(-3.d0*U(-1+i,j,k,1)+4.d0*U(-1+
     &  i,1+j,k,1)-U(-1+i,2+j,k,1))/hy+5.d-1*(-3.d0*U(i,j,
     &  k,1)+4.d0*U(i,1+j,k,1)-U(i,2+j,k,1))/hy)
      else
        dudyw=2.5d-1*(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(
     &  i,-1+j,k,1)+U(i,1+j,k,1))/hy
      endif
      if((masks(i,-1+j,-1+k).gt.0).or.(masks(i,-1+j,k).gt.0)) then 
        dwdyb=5.d-1*(5.d-1*(-3.d0*U(i,j,-1+k,3)+4.d0*U(i,1
     &  +j,-1+k,3)-U(i,2+j,-1+k,3))/hy+5.d-1*(-3.d0*U(i,j,
     &  k,3)+4.d0*U(i,1+j,k,3)-U(i,2+j,k,3))/hy)
      else
        dwdyb=2.5d-1*(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,
     &  1+j,-1+k,3)+U(i,1+j,k,3))/hy
      endif
        dwdyt=trandert(i,j,1+k,3,2)
      if((maskt(-1+i,j,1+k).gt.0).or.(maskt(i,j,1+k).gt.0)) then 
        dudzw=5.d-1*(5.d-1*(U(-1+i,j,-2+k,1)-4.d0*U(-1+i,j
     &  ,-1+k,1)+3.d0*U(-1+i,j,k,1))/hz+5.d-1*(U(i,j,-2+k,
     &  1)-4.d0*U(i,j,-1+k,1)+3.d0*U(i,j,k,1))/hz)
      else
        dudzw=2.5d-1*(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(
     &  i,j,-1+k,1)+U(i,j,1+k,1))/hz
      endif
        dudze=trandere(1+i,j,k,1,3)
        dvdzs=tranders(i,-1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              xflux(i+1,j,k,1)=2.d0*dudxe*muX(1+i,j,k)
              xflux(i+1,j,k,2)=(dudye+dvdxe)*muX(1+i,j,k)
              xflux(i+1,j,k,3)=(dudze+dwdxe)*muX(1+i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

              zflux(i,j,k+1,1)=(dudzt+dwdxt)*muZ(i,j,1+k)
              zflux(i,j,k+1,2)=(dvdzt+dwdyt)*muZ(i,j,1+k)
              zflux(i,j,k+1,3)=2.d0*dwdzt*muZ(i,j,1+k)


c ::: case 22, top-south-west corner
c ::: Null
c ::: tbsw
           i=lo(1)
           j=lo(2)
           k=hi(3)
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzt = (-U(i,j,k,3)+U(i,j,1+k,3))/hz
               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzt = (-U(i,j,k,1)+U(i,j,1+k,1))/hz
               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz
               dvdzt = (-U(i,j,k,2)+U(i,j,1+k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxs=tranders(i,-1+j,k,2,1)
      if((maskw(-1+i,j,-1+k).gt.0).or.(maskw(-1+i,j,k).gt.0)) then 
        dwdxb=5.d-1*(5.d-1*(-3.d0*U(i,j,-1+k,3)+4.d0*U(1+i
     &  ,j,-1+k,3)-U(2+i,j,-1+k,3))/hx+5.d-1*(-3.d0*U(i,j,
     &  k,3)+4.d0*U(1+i,j,k,3)-U(2+i,j,k,3))/hx)
      else
        dwdxb=2.5d-1*(-U(-1+i,j,-1+k,3)-U(-1+i,j,k,3)+U(1+
     &  i,j,-1+k,3)+U(1+i,j,k,3))/hx
      endif
        dwdxt=trandert(i,j,1+k,3,1)
        dudyw=tranderw(-1+i,j,k,1,2)
      if((masks(i,-1+j,-1+k).gt.0).or.(masks(i,-1+j,k).gt.0)) then 
        dwdyb=5.d-1*(5.d-1*(-3.d0*U(i,j,-1+k,3)+4.d0*U(i,1
     &  +j,-1+k,3)-U(i,2+j,-1+k,3))/hy+5.d-1*(-3.d0*U(i,j,
     &  k,3)+4.d0*U(i,1+j,k,3)-U(i,2+j,k,3))/hy)
      else
        dwdyb=2.5d-1*(-U(i,-1+j,-1+k,3)-U(i,-1+j,k,3)+U(i,
     &  1+j,-1+k,3)+U(i,1+j,k,3))/hy
      endif
        dwdyt=trandert(i,j,1+k,3,2)
        dudzw=tranderw(-1+i,j,k,1,3)
        dvdzs=tranders(i,-1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

              zflux(i,j,k+1,1)=(dudzt+dwdxt)*muZ(i,j,1+k)
              zflux(i,j,k+1,2)=(dvdzt+dwdyt)*muZ(i,j,1+k)
              zflux(i,j,k+1,3)=2.d0*dwdzt*muZ(i,j,1+k)

c ::: case 23, bottom-north-east corner
c ::: Null
c ::: bnsew
           i=hi(1)
           j=hi(2)
           k=lo(3)
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if((maske(1+i,-1+j,k).gt.0).or.(maske(1+i,j,k).gt.0)) then 
        dvdxs=5.d-1*(5.d-1*(U(-2+i,-1+j,k,2)-4.d0*U(-1+i,-
     &  1+j,k,2)+3.d0*U(i,-1+j,k,2))/hx+5.d-1*(U(-2+i,j,k,
     &  2)-4.d0*U(-1+i,j,k,2)+3.d0*U(i,j,k,2))/hx)
      else
        dvdxs=2.5d-1*(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+
     &  i,-1+j,k,2)+U(1+i,j,k,2))/hx
      endif
        dwdxb=tranderb(i,j,-1+k,3,1)
        dudye=trandere(1+i,j,k,1,2)
      if((maskn(-1+i,1+j,k).gt.0).or.(maskn(i,1+j,k).gt.0)) then 
        dudyw=5.d-1*(5.d-1*(U(-1+i,-2+j,k,1)-4.d0*U(-1+i,-
     &  1+j,k,1)+3.d0*U(-1+i,j,k,1))/hy+5.d-1*(U(i,-2+j,k,
     &  1)-4.d0*U(i,-1+j,k,1)+3.d0*U(i,j,k,1))/hy)
      else
        dudyw=2.5d-1*(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(
     &  i,-1+j,k,1)+U(i,1+j,k,1))/hy
      endif
        dwdyb=tranderb(i,j,-1+k,3,2)
      if((maskb(-1+i,j,-1+k).gt.0).or.(maskb(i,j,-1+k).gt.0)) then 
        dudzw=5.d-1*(5.d-1*(-3.d0*U(-1+i,j,k,1)+4.d0*U(-1+
     &  i,j,1+k,1)-U(-1+i,j,2+k,1))/hz+5.d-1*(-3.d0*U(i,j,
     &  k,1)+4.d0*U(i,j,1+k,1)-U(i,j,2+k,1))/hz)
      else
        dudzw=2.5d-1*(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(
     &  i,j,-1+k,1)+U(i,j,1+k,1))/hz
      endif
        dudze=trandere(1+i,j,k,1,3)
      if((maskb(i,-1+j,-1+k).gt.0).or.(maskb(i,j,-1+k).gt.0)) then 
        dvdzs=5.d-1*(5.d-1*(-3.d0*U(i,-1+j,k,2)+4.d0*U(i,-
     &  1+j,1+k,2)-U(i,-1+j,2+k,2))/hz+5.d-1*(-3.d0*U(i,j,
     &  k,2)+4.d0*U(i,j,1+k,2)-U(i,j,2+k,2))/hz)
      else
        dvdzs=2.5d-1*(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(
     &  i,j,-1+k,2)+U(i,j,1+k,2))/hz
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              xflux(i+1,j,k,1)=2.d0*dudxe*muX(1+i,j,k)
              xflux(i+1,j,k,2)=(dudye+dvdxe)*muX(1+i,j,k)
              xflux(i+1,j,k,3)=(dudze+dwdxe)*muX(1+i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              yflux(i,j+1,k,1)=(dudyn+dvdxn)*muY(i,1+j,k)
              yflux(i,j+1,k,2)=2.d0*dvdyn*muY(i,1+j,k)
              yflux(i,j+1,k,3)=(dvdzn+dwdyn)*muY(i,1+j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

c ::: case 24, bottom-north-west corner
c ::: Null
c ::: bnsw
           i=lo(1)
           j=hi(2)
           k=lo(3)
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdyn = (-U(i,j,k,2)+U(i,1+j,k,2))/hy
               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudyn = (-U(i,j,k,1)+U(i,1+j,k,1))/hy
               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy
               dwdyn = (-U(i,j,k,3)+U(i,1+j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxn=trandern(i,1+j,k,2,1)
      if((maskw(-1+i,-1+j,k).gt.0).or.(maskw(-1+i,j,k).gt.0)) then 
        dvdxs=5.d-1*(5.d-1*(-3.d0*U(i,-1+j,k,2)+4.d0*U(1+i
     &  ,-1+j,k,2)-U(2+i,-1+j,k,2))/hx+5.d-1*(-3.d0*U(i,j,
     &  k,2)+4.d0*U(1+i,j,k,2)-U(2+i,j,k,2))/hx)
      else
        dvdxs=2.5d-1*(-U(-1+i,-1+j,k,2)-U(-1+i,j,k,2)+U(1+
     &  i,-1+j,k,2)+U(1+i,j,k,2))/hx
      endif
        dwdxb=tranderb(i,j,-1+k,3,1)
        dudyw=tranderw(-1+i,j,k,1,2)
        dwdyb=tranderb(i,j,-1+k,3,2)
        dudzw=tranderw(-1+i,j,k,1,3)
      if((maskb(i,-1+j,-1+k).gt.0).or.(maskb(i,j,-1+k).gt.0)) then 
        dvdzs=5.d-1*(5.d-1*(-3.d0*U(i,-1+j,k,2)+4.d0*U(i,-
     &  1+j,1+k,2)-U(i,-1+j,2+k,2))/hz+5.d-1*(-3.d0*U(i,j,
     &  k,2)+4.d0*U(i,j,1+k,2)-U(i,j,2+k,2))/hz)
      else
        dvdzs=2.5d-1*(-U(i,-1+j,-1+k,2)+U(i,-1+j,1+k,2)-U(
     &  i,j,-1+k,2)+U(i,j,1+k,2))/hz
      endif
        dvdzn=trandern(i,1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              yflux(i,j+1,k,1)=(dudyn+dvdxn)*muY(i,1+j,k)
              yflux(i,j+1,k,2)=2.d0*dvdyn*muY(i,1+j,k)
              yflux(i,j+1,k,3)=(dvdzn+dwdyn)*muY(i,1+j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)


c ::: case 25, bottom-south-east corner
c ::: Null
c ::: bsew
           i=hi(1)
           j=lo(2)
           k=lo(3)
c     :::  get expressions for normal derivatives
               dudxe = (-U(i,j,k,1)+U(1+i,j,k,1))/hx
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxe = (-U(i,j,k,2)+U(1+i,j,k,2))/hx
               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx
               dwdxe = (-U(i,j,k,3)+U(1+i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxb=tranderb(i,j,-1+k,3,1)
        dudye=trandere(1+i,j,k,1,2)
      if((masks(-1+i,-1+j,k).gt.0).or.(masks(i,-1+j,k).gt.0)) then 
        dudyw=5.d-1*(5.d-1*(-3.d0*U(-1+i,j,k,1)+4.d0*U(-1+
     &  i,1+j,k,1)-U(-1+i,2+j,k,1))/hy+5.d-1*(-3.d0*U(i,j,
     &  k,1)+4.d0*U(i,1+j,k,1)-U(i,2+j,k,1))/hy)
      else
        dudyw=2.5d-1*(-U(-1+i,-1+j,k,1)+U(-1+i,1+j,k,1)-U(
     &  i,-1+j,k,1)+U(i,1+j,k,1))/hy
      endif
        dwdyb=tranderb(i,j,-1+k,3,2)
      if((maskb(-1+i,j,-1+k).gt.0).or.(maskb(i,j,-1+k).gt.0)) then 
        dudzw=5.d-1*(5.d-1*(-3.d0*U(-1+i,j,k,1)+4.d0*U(-1+
     &  i,j,1+k,1)-U(-1+i,j,2+k,1))/hz+5.d-1*(-3.d0*U(i,j,
     &  k,1)+4.d0*U(i,j,1+k,1)-U(i,j,2+k,1))/hz)
      else
        dudzw=2.5d-1*(-U(-1+i,j,-1+k,1)+U(-1+i,j,1+k,1)-U(
     &  i,j,-1+k,1)+U(i,j,1+k,1))/hz
      endif
        dudze=trandere(1+i,j,k,1,3)
        dvdzs=tranders(i,-1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              xflux(i+1,j,k,1)=2.d0*dudxe*muX(1+i,j,k)
              xflux(i+1,j,k,2)=(dudye+dvdxe)*muX(1+i,j,k)
              xflux(i+1,j,k,3)=(dudze+dwdxe)*muX(1+i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)

c ::: case 26, bottom-south-west corner
c ::: Null
           i=lo(1)
           j=lo(2)
           k=lo(3)
c     :::  get expressions for normal derivatives
               dudxw = (-U(-1+i,j,k,1)+U(i,j,k,1))/hx

               dvdys = (-U(i,-1+j,k,2)+U(i,j,k,2))/hy

               dwdzb = (-U(i,j,-1+k,3)+U(i,j,k,3))/hz

               dudys = (-U(i,-1+j,k,1)+U(i,j,k,1))/hy

               dvdxw = (-U(-1+i,j,k,2)+U(i,j,k,2))/hx

               dudzb = (-U(i,j,-1+k,1)+U(i,j,k,1))/hz

               dvdzb = (-U(i,j,-1+k,2)+U(i,j,k,2))/hz

               dwdxw = (-U(-1+i,j,k,3)+U(i,j,k,3))/hx

               dwdys = (-U(i,-1+j,k,3)+U(i,j,k,3))/hy

c     ::: get expressions for tangential derivatives
        dvdxs=tranders(i,-1+j,k,2,1)
        dwdxb=tranderb(i,j,-1+k,3,1)
        dudyw=tranderw(-1+i,j,k,1,2)
        dwdyb=tranderb(i,j,-1+k,3,2)
        dudzw=tranderw(-1+i,j,k,1,3)
        dvdzs=tranders(i,-1+j,k,2,3)

c     :::  evaluate expression
              xflux(i,j,k,1)=2.d0*dudxw*muX(i,j,k)
              xflux(i,j,k,2)=(dudyw+dvdxw)*muX(i,j,k)
              xflux(i,j,k,3)=(dudzw+dwdxw)*muX(i,j,k)
              
              yflux(i,j,k,1)=(dudys+dvdxs)*muY(i,j,k)
              yflux(i,j,k,2)=2.d0*dvdys*muY(i,j,k)
              yflux(i,j,k,3)=(dvdzs+dwdys)*muY(i,j,k)

              zflux(i,j,k,1)=(dudzb+dwdxb)*muZ(i,j,k)
              zflux(i,j,k,2)=(dvdzb+dwdyb)*muZ(i,j,k)
              zflux(i,j,k,3)=2.d0*dwdzb*muZ(i,j,k)


#endif


      end
