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
c     tau is a two by two tensor
c     tau = | t_xx     t_xy |
c           | t_xy     t_yy |
c
c     t_xx = 2*mu * u_x
c     t_xy = mu*(u_y + v_x)
c     t_yy = 2*mu * v_y

c     and alpha is a diagonal 2x2 matrix, says rick

c ::: define standard replacements used by Mathematica
c ::: see file visc2d.ma
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
     $     out, DIMS(out),
     $     maskn,DIMS(maskn),
     $     maske,DIMS(maske),
     $     maskw,DIMS(maskw),
     $     masks,DIMS(masks),
     $     trandern,DIMS(trandern),
     $     trandere,DIMS(trandere),
     $     tranderw,DIMS(tranderw),
     $     tranders,DIMS(tranders),
     $     lo,hi,h
     $     )
      REAL_T alpha, beta
      integer DIMDEC(u)
      integer DIMDEC(a)
      integer DIMDEC(muX)
      integer DIMDEC(muY)
      integer DIMDEC(out)
      integer DIMDEC(maskn)
      integer DIMDEC(maske)
      integer DIMDEC(maskw)
      integer DIMDEC(masks)
      integer DIMDEC(trandern)
      integer DIMDEC(trandere)
      integer DIMDEC(tranderw)
      integer DIMDEC(tranders)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      REAL_T h(BL_SPACEDIM)
      REAL_T u(DIMV(u),2)
      REAL_T a(DIMV(a),2)
      REAL_T muX(DIMV(muX))
      REAL_T muY(DIMV(muY))
      REAL_T out(DIMV(out),2)
      integer  maskn(DIMV(maskn))
      integer  maske(DIMV(maske))
      integer  maskw(DIMV(maskw))
      integer  masks(DIMV(masks))

      REAL_T trandern(DIMV(trandern),2)
      REAL_T trandere(DIMV(trandere),2)
      REAL_T tranderw(DIMV(tranderw),2)
      REAL_T tranders(DIMV(tranders),2)

      integer i,j
      REAL_T tauxxw,tauxxe,tauyyn,tauyys
      REAL_T tauxye,tauxyw,tauxyn,tauxys
      REAL_T hx,hy
      REAL_T dudxe,dudxw
      REAL_T dvdyn,dvdys
      REAL_T dudye,dudyw,dudyn,dudys
      REAL_T dvdxe,dvdxw,dvdxn,dvdxs
      logical False, True
      parameter( False=.false.)
      parameter( True=.true. )

      hx = h(1)
      hy = h(2)

c ::: make some mathematica definitions
c ::: allDerivAllUV makes replacements for the deriv expressions and
c ::: maps U and V to the correct values
c ::: Null

c ::: cases:
c     ::: 1) interior.  All normal derivative expressions
c     ::: 2) edge but not on corner.  check mask value to see if use outside
c     :::    value or use interpolations
c     ::: 3) corners.  Must check mask values to see if use outside values
c     :::    values or use interpolations ON TWO SIDES.  But must also
c     :::    check to see if need to do one-sided derivs for tangential
c     :::    derivatives.

c ::: case 1
      do j=lo(2)+1,hi(2)-1
         do i=lo(1)+1,hi(1)-1
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudye =         (-U(i,-1 + j,1) + U(i,1 + j,1) - U(1 + i,-1 + j,1) + 
     -    U(1 + i,1 + j,1))/(4.*hy)
            dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)

            dvdxn =         (-U(-1 + i,j,2) - U(-1 + i,1 + j,2) + U(1 + i,j,2) + 
     -    U(1 + i,1 + j,2))/(4.*hx)
            dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            out(i,j,1) =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            out(i,j,2) =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)
         enddo
      enddo

c ::: ::: case 2a: north side
c ::: Null
c ::: Null
c :::
      j = hi(2)
         do i=lo(1)+1,hi(1)-1
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            if( maskn(i,1 + j).gt.0.or.maskn(1 + i,1 + j).gt.0 ) then
               dudye =         ((U(i,-2 + j,1) - 4*U(i,-1 + j,1) + 3*U(i,j,1))/(2.*hy) + 
     -    (U(1 + i,-2 + j,1) - 4*U(1 + i,-1 + j,1) + 3*U(1 + i,j,1))/(2.*hy))/
     -  2.
            else
               dudye =         (-U(i,-1 + j,1) + U(i,1 + j,1) - U(1 + i,-1 + j,1) + 
     -    U(1 + i,1 + j,1))/(4.*hy)
            endif
            if( maskn(-1 + i,1 + j).gt.0.or.maskn(i,1 + j).gt.0 ) then
               dudyw =         ((U(-1 + i,-2 + j,1) - 4*U(-1 + i,-1 + j,1) + 3*U(-1 + i,j,1))/
     -     (2.*hy) + (U(i,-2 + j,1) - 4*U(i,-1 + j,1) + 3*U(i,j,1))/(2.*hy))/
     -  2.
            else
               dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)
            endif

            dvdxn = trandern(i,j+1,2)
            dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            out(i,j,1) =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            out(i,j,2) =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)
         enddo

c ::: case 2b south side
c ::: Null
c ::: Null
      j=lo(2)
         do i=lo(1)+1,hi(1)-1
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            if( masks(i,-1 + j).gt.0.or.masks(1 + i,-1 + j).gt.0 ) then
               dudye =         ((-3*U(i,j,1) + 4*U(i,1 + j,1) - U(i,2 + j,1))/(2.*hy) + 
     -    (-3*U(1 + i,j,1) + 4*U(1 + i,1 + j,1) - U(1 + i,2 + j,1))/(2.*hy))/
     -  2.
            else
               dudye =         (-U(i,-1 + j,1) + U(i,1 + j,1) - U(1 + i,-1 + j,1) + 
     -    U(1 + i,1 + j,1))/(4.*hy)
            endif
            if( masks(-1 + i,-1 + j).gt.0.or.masks(i,-1 + j).gt.0 ) then
               dudyw =         ((-3*U(-1 + i,j,1) + 4*U(-1 + i,1 + j,1) - U(-1 + i,2 + j,1))/
     -     (2.*hy) + (-3*U(i,j,1) + 4*U(i,1 + j,1) - U(i,2 + j,1))/(2.*hy))/2.
            else
               dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)
            endif

            dvdxn =         (-U(-1 + i,j,2) - U(-1 + i,1 + j,2) + U(1 + i,j,2) + 
     -    U(1 + i,1 + j,2))/(4.*hx)
            dvdxs = tranders(i,j-1,2)
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            out(i,j,1) =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            out(i,j,2) =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)
         enddo


c ::: case 2c: west side story
c ::: Null
c ::: Null
      i=lo(1)
      do j=lo(2)+1,hi(2)-1
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudye =         (-U(i,-1 + j,1) + U(i,1 + j,1) - U(1 + i,-1 + j,1) + 
     -    U(1 + i,1 + j,1))/(4.*hy)
            dudyw = tranderw(i-1,j,1)

            if( maskw(-1 + i,j).gt.0.or.maskw(-1 + i,1 + j).gt.0 ) then
               dvdxn =         ((-3*U(i,j,2) + 4*U(1 + i,j,2) - U(2 + i,j,2))/(2.*hx) + 
     -    (-3*U(i,1 + j,2) + 4*U(1 + i,1 + j,2) - U(2 + i,1 + j,2))/(2.*hx))/
     -  2.
            else
               dvdxn =         (-U(-1 + i,j,2) - U(-1 + i,1 + j,2) + U(1 + i,j,2) + 
     -    U(1 + i,1 + j,2))/(4.*hx)
            endif
            if( maskw(-1 + i,-1 + j).gt.0.or.maskw(-1 + i,j).gt.0 ) then
               dvdxs =         ((-3*U(i,-1 + j,2) + 4*U(1 + i,-1 + j,2) - U(2 + i,-1 + j,2))/
     -     (2.*hx) + (-3*U(i,j,2) + 4*U(1 + i,j,2) - U(2 + i,j,2))/(2.*hx))/2.
            else
               dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            endif
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            out(i,j,1) =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            out(i,j,2) =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)
      enddo

c ::: case 2d: east side story
c ::: Null
c ::: Null
      i=hi(1)
      do j=lo(2)+1,hi(2)-1
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudye = trandere(i+1,j,1)
            dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)

            if( maske(1 + i,j).gt.0.or.maske(1 + i,1 + j).gt.0 ) then
               dvdxn =         ((U(-2 + i,j,2) - 4*U(-1 + i,j,2) + 3*U(i,j,2))/(2.*hx) + 
     -    (U(-2 + i,1 + j,2) - 4*U(-1 + i,1 + j,2) + 3*U(i,1 + j,2))/(2.*hx))/
     -  2.
            else
               dvdxn =         (-U(-1 + i,j,2) - U(-1 + i,1 + j,2) + U(1 + i,j,2) + 
     -    U(1 + i,1 + j,2))/(4.*hx)
            endif
            if( maske(1 + i,-1 + j).gt.0.or.maske(1 + i,j).gt.0 ) then
               dvdxs =         ((U(-2 + i,-1 + j,2) - 4*U(-1 + i,-1 + j,2) + 3*U(i,-1 + j,2))/
     -     (2.*hx) + (U(-2 + i,j,2) - 4*U(-1 + i,j,2) + 3*U(i,j,2))/(2.*hx))/
     -  2.
            else
               dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            endif
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            out(i,j,1) =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            out(i,j,2) =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)
         enddo

c ::: case 3a: north-east corner
c ::: Null
c ::: Null
      j=hi(2)
      i=hi(1)
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudye = trandere(i+1,j,1)
            if( maskn(-1 + i,1 + j).gt.0.or.maskn(i,1 + j).gt.0 ) then
               dudyw =         ((U(-1 + i,-2 + j,1) - 4*U(-1 + i,-1 + j,1) + 3*U(-1 + i,j,1))/
     -     (2.*hy) + (U(i,-2 + j,1) - 4*U(i,-1 + j,1) + 3*U(i,j,1))/(2.*hy))/
     -  2.
            else
               dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)
            endif

            dvdxn = trandern(i,j+1,2)
            if( maske(1 + i,-1 + j).gt.0.or.maske(1 + i,j).gt.0 ) then
               dvdxs =         ((U(-2 + i,-1 + j,2) - 4*U(-1 + i,-1 + j,2) + 3*U(i,-1 + j,2))/
     -     (2.*hx) + (U(-2 + i,j,2) - 4*U(-1 + i,j,2) + 3*U(i,j,2))/(2.*hx))/
     -  2.
            else
               dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            endif
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            out(i,j,1) =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            out(i,j,2) =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)

c ::: case 3b: south-east corner, plus i, minus j, plus i, plus j
c ::: Null
      j=lo(2)
      i=hi(1)
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudye = trandere(i+1,j,1)
            if( masks(-1 + i,-1 + j).gt.0.or.masks(i,-1 + j).gt.0 ) then 
               dudyw =         ((-3*U(-1 + i,j,1) + 4*U(-1 + i,1 + j,1) - U(-1 + i,2 + j,1))/
     -     (2.*hy) + (-3*U(i,j,1) + 4*U(i,1 + j,1) - U(i,2 + j,1))/(2.*hy))/2.
            else
               dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)
            endif

            if( maske(1 + i,j).gt.0.or.maske(1 + i,1 + j).gt.0 ) then
               dvdxn =         ((U(-2 + i,j,2) - 4*U(-1 + i,j,2) + 3*U(i,j,2))/(2.*hx) + 
     -    (U(-2 + i,1 + j,2) - 4*U(-1 + i,1 + j,2) + 3*U(i,1 + j,2))/(2.*hx))/
     -  2. 
            else
               dvdxn =         (-U(-1 + i,j,2) - U(-1 + i,1 + j,2) + U(1 + i,j,2) + 
     -    U(1 + i,1 + j,2))/(4.*hx)
            endif
            dvdxs = tranders(i,j-1,2)
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            out(i,j,1) =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            out(i,j,2) =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)

c ::: case 3c: south-west corner, minus i, minus j
c ::: Null
      j=lo(2)
      i=lo(1)
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            if( masks(i,-1 + j).gt.0.or.masks(1 + i,-1 + j).gt.0 ) then
               dudye =         ((-3*U(i,j,1) + 4*U(i,1 + j,1) - U(i,2 + j,1))/(2.*hy) + 
     -    (-3*U(1 + i,j,1) + 4*U(1 + i,1 + j,1) - U(1 + i,2 + j,1))/(2.*hy))/
     -  2.
            else
               dudye =         (-U(i,-1 + j,1) + U(i,1 + j,1) - U(1 + i,-1 + j,1) + 
     -    U(1 + i,1 + j,1))/(4.*hy)
            endif
            dudyw = tranderw(i-1,j,1)

            if( maskw(-1 + i,j).gt.0.or.maskw(-1 + i,1 + j).gt.0 ) then
               dvdxn =         ((-3*U(i,j,2) + 4*U(1 + i,j,2) - U(2 + i,j,2))/(2.*hx) + 
     -    (-3*U(i,1 + j,2) + 4*U(1 + i,1 + j,2) - U(2 + i,1 + j,2))/(2.*hx))/
     -  2.
            else
               dvdxn =         (-U(-1 + i,j,2) - U(-1 + i,1 + j,2) + U(1 + i,j,2) + 
     -    U(1 + i,1 + j,2))/(4.*hx)
            endif
            dvdxs = tranders(i,j-1,2)
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            out(i,j,1) =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            out(i,j,2) =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)

c ::: case 3d: north-west corner, minus i, plus j
c ::: Null
      j=hi(2)
      i=lo(1)
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            if( maskn(i,1 + j).gt.0.or.maskn(1 + i,1 + j).gt.0 ) then
               dudye =         ((U(i,-2 + j,1) - 4*U(i,-1 + j,1) + 3*U(i,j,1))/(2.*hy) + 
     -    (U(1 + i,-2 + j,1) - 4*U(1 + i,-1 + j,1) + 3*U(1 + i,j,1))/(2.*hy))/
     -  2.
            else
               dudye =         (-U(i,-1 + j,1) + U(i,1 + j,1) - U(1 + i,-1 + j,1) + 
     -    U(1 + i,1 + j,1))/(4.*hy)
            endif
            dudyw = tranderw(i-1,j,1)

            dvdxn = trandern(i,j+1,2)
            if( maskw(-1 + i,-1 + j).gt.0.or.maskw(-1 + i,j).gt.0 ) then
               dvdxs =         ((-3*U(i,-1 + j,2) + 4*U(1 + i,-1 + j,2) - U(2 + i,-1 + j,2))/
     -     (2.*hx) + (-3*U(i,j,2) + 4*U(1 + i,j,2) - U(2 + i,j,2))/(2.*hx))/2.
            else
               dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            endif
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            out(i,j,1) =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            out(i,j,2) =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)



      return
      end


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
     $     xflux,DIMS(xflux),
     $     yflux,DIMS(yflux),
     $     maskn,DIMS(maskn),
     $     maske,DIMS(maske),
     $     maskw,DIMS(maskw),
     $     masks,DIMS(masks),
     $     trandern,DIMS(trandern),
     $     trandere,DIMS(trandere),
     $     tranderw,DIMS(tranderw),
     $     tranders,DIMS(tranders),
     $     lo,hi,h
     $     )
      REAL_T alpha, beta
      integer DIMDEC(u)
      integer DIMDEC(a)
      integer DIMDEC(muX)
      integer DIMDEC(muY)
      integer DIMDEC(xflux)
      integer DIMDEC(yflux)
      integer DIMDEC(maskn)
      integer DIMDEC(maske)
      integer DIMDEC(maskw)
      integer DIMDEC(masks)
      integer DIMDEC(trandern)
      integer DIMDEC(trandere)
      integer DIMDEC(tranderw)
      integer DIMDEC(tranders)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      REAL_T h(BL_SPACEDIM)
      REAL_T u(DIMV(u),2)
      REAL_T a(DIMV(a),2)
      REAL_T muX(DIMV(muX))
      REAL_T muY(DIMV(muY))
      REAL_T xflux(DIMV(xflux),2)
      REAL_T yflux(DIMV(yflux),2)
      integer  maskn(DIMV(maskn))
      integer  maske(DIMV(maske))
      integer  maskw(DIMV(maskw))
      integer  masks(DIMV(masks))

      REAL_T trandern(DIMV(trandern),2)
      REAL_T trandere(DIMV(trandere),2)
      REAL_T tranderw(DIMV(tranderw),2)
      REAL_T tranders(DIMV(tranders),2)

      integer i,j
      REAL_T tauxxw,tauxxe,tauyyn,tauyys
      REAL_T tauxye,tauxyw,tauxyn,tauxys
      REAL_T hx,hy
      REAL_T dudxe,dudxw
      REAL_T dvdyn,dvdys
      REAL_T dudye,dudyw,dudyn,dudys
      REAL_T dvdxe,dvdxw,dvdxn,dvdxs
      logical False, True
      parameter( False=.false.)
      parameter( True=.true. )

      hx = h(1)
      hy = h(2)

c ::: make some mathematica definitions
c ::: allDerivAllUV makes replacements for the deriv expressions and
c ::: maps U and V to the correct values
c ::: Null

c ::: cases:
c     ::: 1) interior.  All normal derivative expressions
c     ::: 2) edge but not on corner.  check mask value to see if use outside
c     :::    value or use interpolations
c     ::: 3) corners.  Must check mask values to see if use outside values
c     :::    values or use interpolations ON TWO SIDES.  But must also
c     :::    check to see if need to do one-sided derivs for tangential
c     :::    derivatives.

c ::: case 1
      do j=lo(2)+1,hi(2)
         do i=lo(1)+1,hi(1)
c     :::  get expressions for appropriate derivatives
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)

            dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            xflux(i,j,1) = -(2*dudxw*muX(i,j))
            xflux(i,j,2) = -((dudyw + dvdxw)*muX(i,j))
            yflux(i,j,1) = -((dudys + dvdxs)*muY(i,j))
            yflux(i,j,2) = -(2*dvdys*muY(i,j))
         enddo
      enddo

c ::: ::: case 2a: north side
c ::: Null
c ::: Null
c :::
      j = hi(2)
         do i=lo(1)+1,hi(1)-1
c     :::  get expressions for appropriate derivatives
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            if( maskn(-1 + i,1 + j).gt.0.or.maskn(i,1 + j).gt.0 ) then
               dudyw =         ((U(-1 + i,-2 + j,1) - 4*U(-1 + i,-1 + j,1) + 3*U(-1 + i,j,1))/
     -     (2.*hy) + (U(i,-2 + j,1) - 4*U(i,-1 + j,1) + 3*U(i,j,1))/(2.*hy))/
     -  2.
            else
               dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)
            endif

            dvdxn = trandern(i,j+1,2)
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            xflux(i,j,1) = -(2*dudxw*muX(i,j))
            xflux(i,j,2) = -((dudyw + dvdxw)*muX(i,j))
            yflux(i,j+1,1) = -((dudyn + dvdxn)*muY(i,1 + j))
            yflux(i,j+1,2) = -(2*dvdyn*muY(i,1 + j))
         enddo

c ::: case 2b south side
c ::: Null
c ::: Null
      j=lo(2)
         do i=lo(1)+1,hi(1)-1
c     :::  get expressions for appropriate derivatives
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            if( masks(-1 + i,-1 + j).gt.0.or.masks(i,-1 + j).gt.0 ) then
               dudyw =         ((-3*U(-1 + i,j,1) + 4*U(-1 + i,1 + j,1) - U(-1 + i,2 + j,1))/
     -     (2.*hy) + (-3*U(i,j,1) + 4*U(i,1 + j,1) - U(i,2 + j,1))/(2.*hy))/2.
            else
               dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)
            endif

            dvdxs = tranders(i,j-1,2)
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            xflux(i,j,1) = -(2*dudxw*muX(i,j))
            xflux(i,j,2) = -((dudyw + dvdxw)*muX(i,j))
            yflux(i,j,1) = -((dudys + dvdxs)*muY(i,j))
            yflux(i,j,2) = -(2*dvdys*muY(i,j))
         enddo


c ::: case 2c: west side story
c ::: Null
c ::: Null
      i=lo(1)
      do j=lo(2)+1,hi(2)-1
c     :::  get expressions for appropriate derivatives
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudyw = tranderw(i-1,j,1)

            if( maskw(-1 + i,-1 + j).gt.0.or.maskw(-1 + i,j).gt.0 ) then
               dvdxs =         ((-3*U(i,-1 + j,2) + 4*U(1 + i,-1 + j,2) - U(2 + i,-1 + j,2))/
     -     (2.*hx) + (-3*U(i,j,2) + 4*U(1 + i,j,2) - U(2 + i,j,2))/(2.*hx))/2.
            else
               dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            endif
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            xflux(i,j,1) = -(2*dudxw*muX(i,j))
            xflux(i,j,2) = -((dudyw + dvdxw)*muX(i,j))
            yflux(i,j,1) = -((dudys + dvdxs)*muY(i,j))
            yflux(i,j,2) = -(2*dvdys*muY(i,j))
      enddo

c ::: case 2d: east side story
c ::: Null
c ::: Null
      i=hi(1)
      do j=lo(2)+1,hi(2)-1
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx

            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudye = trandere(i+1,j,1)

            if( maske(1 + i,-1 + j).gt.0.or.maske(1 + i,j).gt.0 ) then
               dvdxs =         ((U(-2 + i,-1 + j,2) - 4*U(-1 + i,-1 + j,2) + 3*U(i,-1 + j,2))/
     -     (2.*hx) + (U(-2 + i,j,2) - 4*U(-1 + i,j,2) + 3*U(i,j,2))/(2.*hx))/
     -  2.
            else
               dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            endif
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
c     :::  evaluate expression
            xflux(i+1,j,1) = -(2*dudxe*muX(1 + i,j))
            xflux(i+1,j,2) = -((dudye + dvdxe)*muX(1 + i,j))
            yflux(i,j,1) = -((dudys + dvdxs)*muY(i,j))
            yflux(i,j,2) = -(2*dvdys*muY(i,j))
         enddo

c ::: case 3a: north-east corner
c ::: Null
c ::: Null
      j=hi(2)
      i=hi(1)
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudye = trandere(i+1,j,1)
            if( maskn(-1 + i,1 + j).gt.0.or.maskn(i,1 + j).gt.0 ) then
               dudyw =         ((U(-1 + i,-2 + j,1) - 4*U(-1 + i,-1 + j,1) + 3*U(-1 + i,j,1))/
     -     (2.*hy) + (U(i,-2 + j,1) - 4*U(i,-1 + j,1) + 3*U(i,j,1))/(2.*hy))/
     -  2.
            else
               dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)
            endif

            dvdxn = trandern(i,j+1,2)
            if( maske(1 + i,-1 + j).gt.0.or.maske(1 + i,j).gt.0 ) then
               dvdxs =         ((U(-2 + i,-1 + j,2) - 4*U(-1 + i,-1 + j,2) + 3*U(i,-1 + j,2))/
     -     (2.*hx) + (U(-2 + i,j,2) - 4*U(-1 + i,j,2) + 3*U(i,j,2))/(2.*hx))/
     -  2.
            else
               dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            endif
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            xflux(i,j,1) = -(2*dudxw*muX(i,j))
            xflux(i,j,2) = -((dudyw + dvdxw)*muX(i,j))
            xflux(i+1,j,1) = -(2*dudxe*muX(1 + i,j))
            xflux(i+1,j,2) = -((dudye + dvdxe)*muX(1 + i,j))
            yflux(i,j,1) = -((dudys + dvdxs)*muY(i,j))
            yflux(i,j,2) = -(2*dvdys*muY(i,j))
            yflux(i,j+1,1) = -((dudyn + dvdxn)*muY(i,1 + j))
            yflux(i,j+1,2) = -(2*dvdyn*muY(i,1 + j))

c ::: case 3b: south-east corner, plus i, minus j, plus i, plus j
c ::: Null
      j=lo(2)
      i=hi(1)
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudye = trandere(i+1,j,1)
            if( masks(-1 + i,-1 + j).gt.0.or.masks(i,-1 + j).gt.0 ) then 
               dudyw =         ((-3*U(-1 + i,j,1) + 4*U(-1 + i,1 + j,1) - U(-1 + i,2 + j,1))/
     -     (2.*hy) + (-3*U(i,j,1) + 4*U(i,1 + j,1) - U(i,2 + j,1))/(2.*hy))/2.
            else
               dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)
            endif

            if( maske(1 + i,j).gt.0.or.maske(1 + i,1 + j).gt.0 ) then
               dvdxn =         ((U(-2 + i,j,2) - 4*U(-1 + i,j,2) + 3*U(i,j,2))/(2.*hx) + 
     -    (U(-2 + i,1 + j,2) - 4*U(-1 + i,1 + j,2) + 3*U(i,1 + j,2))/(2.*hx))/
     -  2. 
            else
               dvdxn =         (-U(-1 + i,j,2) - U(-1 + i,1 + j,2) + U(1 + i,j,2) + 
     -    U(1 + i,1 + j,2))/(4.*hx)
            endif
            dvdxs = tranders(i,j-1,2)
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            xflux(i,j,1) = -(2*dudxw*muX(i,j))
            xflux(i,j,2) = -((dudyw + dvdxw)*muX(i,j))
            xflux(i+1,j,1) = -(2*dudxe*muX(1 + i,j))
            xflux(i+1,j,2) = -((dudye + dvdxe)*muX(1 + i,j))
            yflux(i,j,1) = -((dudys + dvdxs)*muY(i,j))
            yflux(i,j,2) = -(2*dvdys*muY(i,j))
            yflux(i,j+1,1) = -((dudyn + dvdxn)*muY(i,1 + j))
            yflux(i,j+1,2) = -(2*dvdyn*muY(i,1 + j))

c ::: case 3c: south-west corner, minus i, minus j
c ::: Null
      j=lo(2)
      i=lo(1)
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            if( masks(i,-1 + j).gt.0.or.masks(1 + i,-1 + j).gt.0 ) then
               dudye =         ((-3*U(i,j,1) + 4*U(i,1 + j,1) - U(i,2 + j,1))/(2.*hy) + 
     -    (-3*U(1 + i,j,1) + 4*U(1 + i,1 + j,1) - U(1 + i,2 + j,1))/(2.*hy))/
     -  2.
            else
               dudye =         (-U(i,-1 + j,1) + U(i,1 + j,1) - U(1 + i,-1 + j,1) + 
     -    U(1 + i,1 + j,1))/(4.*hy)
            endif
            dudyw = tranderw(i-1,j,1)

            if( maskw(-1 + i,j).gt.0.or.maskw(-1 + i,1 + j).gt.0 ) then
               dvdxn =         ((-3*U(i,j,2) + 4*U(1 + i,j,2) - U(2 + i,j,2))/(2.*hx) + 
     -    (-3*U(i,1 + j,2) + 4*U(1 + i,1 + j,2) - U(2 + i,1 + j,2))/(2.*hx))/
     -  2.
            else
               dvdxn =         (-U(-1 + i,j,2) - U(-1 + i,1 + j,2) + U(1 + i,j,2) + 
     -    U(1 + i,1 + j,2))/(4.*hx)
            endif
            dvdxs = tranders(i,j-1,2)
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            xflux(i,j,1) = -(2*dudxw*muX(i,j))
            xflux(i,j,2) = -((dudyw + dvdxw)*muX(i,j))
            xflux(i+1,j,1) = -(2*dudxe*muX(1 + i,j))
            xflux(i+1,j,2) = -((dudye + dvdxe)*muX(1 + i,j))
            yflux(i,j,1) = -((dudys + dvdxs)*muY(i,j))
            yflux(i,j,2) = -(2*dvdys*muY(i,j))
            yflux(i,j+1,1) = -((dudyn + dvdxn)*muY(i,1 + j))
            yflux(i,j+1,2) = -(2*dvdyn*muY(i,1 + j))

c ::: case 3d: north-west corner, minus i, plus j
c ::: Null
      j=hi(2)
      i=lo(1)
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            if( maskn(i,1 + j).gt.0.or.maskn(1 + i,1 + j).gt.0 ) then
               dudye =         ((U(i,-2 + j,1) - 4*U(i,-1 + j,1) + 3*U(i,j,1))/(2.*hy) + 
     -    (U(1 + i,-2 + j,1) - 4*U(1 + i,-1 + j,1) + 3*U(1 + i,j,1))/(2.*hy))/
     -  2.
            else
               dudye =         (-U(i,-1 + j,1) + U(i,1 + j,1) - U(1 + i,-1 + j,1) + 
     -    U(1 + i,1 + j,1))/(4.*hy)
            endif
            dudyw = tranderw(i-1,j,1)

            dvdxn = trandern(i,j+1,2)
            if( maskw(-1 + i,-1 + j).gt.0.or.maskw(-1 + i,j).gt.0 ) then
               dvdxs =         ((-3*U(i,-1 + j,2) + 4*U(1 + i,-1 + j,2) - U(2 + i,-1 + j,2))/
     -     (2.*hx) + (-3*U(i,j,2) + 4*U(1 + i,j,2) - U(2 + i,j,2))/(2.*hx))/2.
            else
               dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            endif
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  evaluate expression
            xflux(i,j,1) = -(2*dudxw*muX(i,j))
            xflux(i,j,2) = -((dudyw + dvdxw)*muX(i,j))
            xflux(i+1,j,1) = -(2*dudxe*muX(1 + i,j))
            xflux(i+1,j,2) = -((dudye + dvdxe)*muX(1 + i,j))
            yflux(i,j,1) = -((dudys + dvdxs)*muY(i,j))
            yflux(i,j,2) = -(2*dvdys*muY(i,j))
            yflux(i,j+1,1) = -((dudyn + dvdxn)*muY(i,1 + j))
            yflux(i,j+1,2) = -(2*dvdyn*muY(i,1 + j))



      return
      end



c-----------------------------------------------------------------------
      subroutine FORT_APPLYBC (
     $     flagden, flagbc, maxorder,
     $     u,   DIMS(u),
     $     cdir, bct, bcl,
     $     bcval, DIMS(bcval),
     $     mask,  DIMS(mask),
     $     maskp, DIMS(maskp),
     $     maskm, DIMS(maskm),
     $     den,   DIMS(den),
     $     exttd, DIMS(exttd),
     $     tander,DIMS(tander),
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
c     The bc type = LO_REFLECT_ODD is a special type of boundary condition.

c ::: the old FORT_APPLYBC put values out in ghost cells.  The new
c ::: FORT_APPLYBC makes a distinction between normal and tangential derivs.
c ::: the normal derivatives are still evaluated with ghost cells ( and
c ::: the den array for adjusting the diagonal element).  Tangential
c ::: derivatives are NOT computed from the ghost cells (note: we are
c ::: discussing here the tangential derivates which are centered on the
c ::: external cells walls of the rectangular domain.  Internal tangential
c ::: derivatives whose stencil extends outside the rectangular domain have
c ::: still ANOTHER calculational trick, one sided derivatives, which is 
c ::: implemented in the apply operator).  For these tangential derivatives,
c ::: if some of the cells in the stencil are uncovered, the derivative is
c ::: computed half from cells which are internal, plus an externally supplied
c ::: tangential derivative centered at location bcl.
c ??? exttd is the externally supplied tangential derivative which is
c ??? centered at bcl from wall.  Since the location of exttd is
c ??? indeterminate (determined by bcl), we arbitrarily index it in the
c ??? cell just outside the wall.
c ??? tander is the calculated tangential derivative which is centered at
c ??? at the wall.  It is indexed as if it were a cell centered quantity
c ??? just outside the wall.  This is anomalous since it is really edge
c ??? centered on the outer wall.
c      
c ::: other notes since previous developers didn't bother to document
c ::: cdir is mnemonic for coordinate direction, i.e. which side this is
c ::: cdir==0 -> left wall
c ::: cdir==2 -> right
c ::: cdir==1 -> bottom
c ::: cdir==3 -> top

c ::: notes about mask
c ::: since we now require other masks, there are 3.
c ::: 1) mask -- the main one covering the cells off the side of valid domain.
c ::: 2) maskp -- the plus end of transverse direction
c ::: 3) maskn -- the minus end of transverse direction

      integer maxorder
      integer nc, cdir, flagden, flagbc
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      integer DIMDEC(u)
      REAL_T u(DIMV(u),nc)
      integer DIMDEC(den)
      REAL_T den(DIMV(den),nc)
      integer DIMDEC(exttd)
      REAL_T exttd(DIMV(exttd),nc)
      integer DIMDEC(tander)
      REAL_T tander(DIMV(tander),nc)
      integer DIMDEC(bcval)
      REAL_T bcval(DIMV(bcval),nc)
      integer DIMDEC(mask)
      integer mask(DIMV(mask))
      integer DIMDEC(maskp)
      integer maskp(DIMV(maskp))
      integer DIMDEC(maskm)
      integer maskm(DIMV(maskm))
      integer bct(nc)
      REAL_T bcl
      REAL_T h(BL_SPACEDIM)
c
      REAL_T hx,hy,lambda
      REAL_T innder,outder
      REAL_T innloc,outloc,edgloc
      integer i
      integer j
      integer n
      logical is_dirichlet
      logical is_neumann
      logical is_odd
      REAL_T xb
c
      integer lenx
      integer leny
      integer m
c
      integer Lmaxorder
      integer maxmaxorder
      parameter(maxmaxorder=4)
      REAL_T x(-1:maxmaxorder-2)
      REAL_T coef(-1:maxmaxorder-2)
      REAL_T xInt
c
      logical False, True
      parameter( False=.false.)
      parameter( True=.true. )
c
      is_dirichlet(i) = ( i .eq. LO_DIRICHLET   )
      is_neumann(i)   = ( i .eq. LO_NEUMANN )
      hx = h(1)
      hy = h(2)
c
      if ( maxorder .eq. -1 ) then
         Lmaxorder = maxmaxorder
      else
         Lmaxorder = MIN(maxorder,maxmaxorder)
      endif
      lenx = MIN(hi(1)-lo(1), Lmaxorder-2)
      leny = MIN(hi(2)-lo(2), Lmaxorder-2)
c
c     TODO:
c     In order for this to work with growing multigrid, must
c     sort xa[] because it is possible for the xb value to lay
c     within this range.
c     
c     The Left face of the grid
c
      if(cdir .eq. 0) then
         do n = 1, nc
            if (is_neumann(bct(n))) then
               do j = lo(2), hi(2)
                  u(lo(1)-1,j,n) = cvmgt(
     $                 u(lo(1),j,n),
     $                 u(lo(1)-1,j,n),
     $                 mask(lo(1)-1,j) .gt. 0)
               enddo
               if ( flagden .eq. 1) then
                  do j = lo(2), hi(2)
                     den(lo(1),j,n) = 1.0
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: Null
c     :::::: interior part of left side
               i = lo(1)
               do j=lo(2)+1,hi(2)-1
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                  tander(i-1,j,n) = innder
               enddo
c     :::::: now end points
c     :::::: lower leftside
               j = lo(2)
               if( maskm(i,-1 + j).gt.0 ) then
                    innder = (-3*U(i,j,n) + 4*U(i,1 + j,n) - U(i,2 + j,n))/(2.*hy)
               else
                   innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
               endif
               tander(i-1,j,n) = innder
c     :::::: upper left side
                j = hi(2)
                if( maskp(i,1 + j).gt.0 ) then
                   innder = (U(i,-2 + j,n) - 4*U(i,-1 + j,n) + 3*U(i,j,n))/(2.*hy)
                else
                   innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                endif
                tander(i-1,j,n) = innder
            else if (is_dirichlet(bct(n))) then
               do m=0,lenx
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(1)
               xInt = - 0.5
               call polyInterpCoeff(xInt, x, lenx+2, coef)
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     u(lo(1)-1, j, n) = cvmgt(
     $                    bcval(lo(1)-1,j,n)*coef(-1),
     $                    u(lo(1)-1, j, n),
     $                    mask(lo(1)-1,j) .gt. 0)
                  enddo
               else
                  do j = lo(2), hi(2)
                     u(lo(1)-1, j, n) = cvmgt(
     $                    0.0,
     $                    u(lo(1)-1, j, n),
     $                    mask(lo(1)-1,j) .gt. 0)
                  enddo
               endif
               do m = 0, lenx
                  do j = lo(2), hi(2)
                     u(lo(1)-1,j,n) = cvmgt(
     $                    u(lo(1)-1,j,n)
     $                    + u(lo(1)+m, j, n)*coef(m),
     $                    u(lo(1)-1,j,n),
     $                    mask(lo(1)-1,j) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(lo(1),j,n) = cvmgt(coef(0), 0.0,
     $                    mask(lo(1)-1,j) .gt. 0)
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: Null
c     :::::: interior part of left side
               edgloc = 0.
               i = lo(1)
               do j=lo(2)+1,hi(2)-1
                  if( mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-U(-1 + i,-1 + j,n) + U(-1 + i,1 + j,n))/(2.*hy)
                  elseif( mask(-1 + i,j).eq.0.and.mask(-1 + i,1 + j).eq.0.and.mask(-1 + i,2 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(-1 + i,j,n) + 4*U(-1 + i,1 + j,n) - U(-1 + i,2 + j,n))/(2.*hy)
                  elseif( mask(-1 + i,-2 + j).eq.0.and.mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(-1 + i,-2 + j,n) - 4*U(-1 + i,-1 + j,n) + 3*U(-1 + i,j,n))/(2.*hy)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i-1,j,n)
                     else
                        outder = 0.
                     endif
                  endif
                  innloc = 0.5
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
            
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  tander(i-1,j,n) = lambda*innder+(1-lambda)*outder
               enddo
c     :::::: now end points
c     :::::: lower leftside
               j = lo(2)
               if( mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0 ) then
                   outloc = -0.5
                   outder = (-U(-1 + i,-1 + j,n) + U(-1 + i,1 + j,n))/(2.*hy)
               elseif( mask(-1 + i,j).eq.0.and.mask(-1 + i,1 + j).eq.0.and.mask(-1 + i,2 + j).eq.0 ) then
                   outloc = -0.5
                   outder = (-3*U(-1 + i,j,n) + 4*U(-1 + i,1 + j,n) - U(-1 + i,2 + j,n))/(2.*hy)
               else
                   outloc = x(-1)
                   if( flagbc .eq. 1) then
                      outder = exttd(i-1,j,n)
                   else
                      outder = 0.
                   endif
               endif
               if( maskm(i,-1 + j).gt.0 ) then
                    innder = (-3*U(i,j,n) + 4*U(i,1 + j,n) - U(i,2 + j,n))/(2.*hy)
               else
                   innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i-1,j,n) = lambda*innder+(1-lambda)*outder
c     :::::: upper left side
                j = hi(2)
                if( mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0 ) then
                   outloc = -0.5
                   outder = (-U(-1 + i,-1 + j,n) + U(-1 + i,1 + j,n))/(2.*hy)
                elseif( mask(-1 + i,-2 + j).eq.0.and.mask(-1 + i,-1 + j).eq.0.and.mask(-1 + i,j).eq.0 ) then
                   outloc = -0.5
                   outder = (U(-1 + i,-2 + j,n) - 4*U(-1 + i,-1 + j,n) + 3*U(-1 + i,j,n))/(2.*hy)
                else
                   outloc = x(-1)
                   if( flagbc .eq. 1) then
                      outder = exttd(i-1,j,n)
                   else
                      outder = 0.
                   endif
                endif
                if( maskp(i,1 + j).gt.0 ) then
                   innder = (U(i,-2 + j,n) - 4*U(i,-1 + j,n) + 3*U(i,j,n))/(2.*hy)
                else
                   innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                endif
                innloc = 0.5
                lambda = (edgloc-outloc)/(innloc-outloc)
                tander(i-1,j,n) = lambda*innder+(1-lambda)*outder
               
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               
               do j = lo(2), hi(2)
                  u(lo(1)-1, j, n) = cvmgt(
     $                 -u(lo(1),j,n),
     $                 u(lo(1)-1, j, n),
     $                 mask(lo(1)-1,j) .gt. 0)
                  tander(lo(1)-1,j,n) = 0.0
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(lo(1),j,n) = cvmgt(-1.0, 0.0,
     $                    mask(lo(1)-1,j) .gt. 0)
                  enddo
               endif
               
            else 
               print *,'UNKNOWN BC ON LEFT FACE IN APPLYBC'
               stop
            endif
         enddo
      endif
c     
c     The Right face of the grid
c
      if(cdir .eq. 2) then
         do n = 1, nc
            if(is_neumann(bct(n))) then
               do j = lo(2), hi(2)
                  u(hi(1)+1,j,n) = cvmgt(
     $                 u(hi(1), j, n),
     $                 u(hi(1)+1, j, n),
     $                 mask(hi(1)+1,j) .gt. 0)
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(hi(1),j,n) = 1.0
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: Null
c     :::::: interior part of right side
               i = hi(1)
               do j=lo(2)+1,hi(2)-1
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
                  tander(i+1,j,n) = innder
               enddo
c :::::: now end points
c :::::: lower right side
               j = lo(2)
               if( maskm(i,-1 + j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(i,1 + j,n) - U(i,2 + j,n))/(2.*hy)
               else
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
               endif
               tander(i+1,j,n) = innder
c     :::::: upper right side
               j = hi(2)
               if( maskp(i,1 + j).gt.0 ) then
                  innder = (U(i,-2 + j,n) - 4*U(i,-1 + j,n) + 3*U(i,j,n))/(2.*hy)
               else
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
               endif
               tander(i+1,j,n) = innder
            else if (is_dirichlet(bct(n))) then
               do m=0,lenx
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(1)
               xInt = - 0.5
               call polyInterpCoeff(xInt, x, lenx+2, coef)
               if ( flagbc .eq. 1 ) then
                  do j = lo(2), hi(2)
                     u(hi(1)+1,j,n) = cvmgt(
     $                    bcval(hi(1)+1,j,n)*coef(-1),
     $                    u(hi(1)+1,j,n),
     $                    mask(hi(1)+1,j) .gt. 0)
                  enddo
               else
                  do j = lo(2), hi(2)
                     u(hi(1)+1,j,n) = cvmgt(
     $                    0.0,
     $                    u(hi(1)+1,j,n),
     $                    mask(hi(1)+1,j) .gt. 0)
                  enddo
               endif
               do m = 0, lenx
                  do j = lo(2), hi(2)
                     u(hi(1)+1,j,n) = cvmgt(
     $                    u(hi(1)+1,j,n)
     $                    + u(hi(1)-m,j,n)*coef(m),
     $                    u(hi(1)+1,j,n),
     $                    mask(hi(1)+1,j) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(hi(1),j,n)   = cvmgt(coef(0), 0.0,
     $                    mask(hi(1)+1,j) .gt. 0)
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: Null
c     :::::: interior part of right side
               edgloc = 0.
               i = hi(1)
               do j=lo(2)+1,hi(2)-1
                  if( mask(1 + i,-1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-U(1 + i,-1 + j,n) + U(1 + i,1 + j,n))/(2.*hy)
                  elseif( mask(1 + i,j).eq.0.and.mask(1 + i,1 + j).eq.0.and.mask(1 + i,2 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(1 + i,j,n) + 4*U(1 + i,1 + j,n) - U(1 + i,2 + j,n))/(2.*hy)
                  elseif( mask(1 + i,-2 + j).eq.0.and.mask(1 + i,-1 + j).eq.0.and.mask(1 + i,j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(1 + i,-2 + j,n) - 4*U(1 + i,-1 + j,n) + 3*U(1 + i,j,n))/(2.*hy)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i+1,j,n)
                     else
                        outder = 0.
                     endif
                  endif
                  innloc = 0.5
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)

                  lambda = (edgloc-outloc)/(innloc-outloc)
                  tander(i+1,j,n) = lambda*innder+(1-lambda)*outder
               enddo
c :::::: now end points
c :::::: lower right side
               j = lo(2)
               if( mask(1 + i,-1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(1 + i,-1 + j,n) + U(1 + i,1 + j,n))/(2.*hy)
               elseif( mask(1 + i,j).eq.0.and.mask(1 + i,1 + j).eq.0.and.mask(1 + i,2 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-3*U(1 + i,j,n) + 4*U(1 + i,1 + j,n) - U(1 + i,2 + j,n))/(2.*hy)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i+1,j,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskm(i,-1 + j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(i,1 + j,n) - U(i,2 + j,n))/(2.*hy)
               else
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i+1,j,n) = lambda*innder+(1-lambda)*outder
c     :::::: upper right side
               j = hi(2)
               if( mask(1 + i,-1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(1 + i,-1 + j,n) + U(1 + i,1 + j,n))/(2.*hy)
               elseif( mask(1 + i,-2 + j).eq.0.and.mask(1 + i,-1 + j).eq.0.and.mask(1 + i,j).eq.0 ) then
                  outloc = -0.5
                  outder = (U(1 + i,-2 + j,n) - 4*U(1 + i,-1 + j,n) + 3*U(1 + i,j,n))/(2.*hy)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i+1,j,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskp(i,1 + j).gt.0 ) then
                  innder = (U(i,-2 + j,n) - 4*U(i,-1 + j,n) + 3*U(i,j,n))/(2.*hy)
               else
                  innder = (-U(i,-1 + j,n) + U(i,1 + j,n))/(2.*hy)
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i+1,j,n) = lambda*innder+(1-lambda)*outder
               
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               
               do j = lo(2), hi(2)
                  u(hi(1)+1, j, n) = cvmgt(
     $                 -u(hi(1),j,n),
     $                 u(hi(1)+1, j, n),
     $                 mask(hi(1)+1,j) .gt. 0)
c                 ::: tangential der is trivial
                  tander(hi(1)+1,j,n) = 0.0
               enddo
               if ( flagden .eq. 1 ) then
                  do j = lo(2), hi(2)
                     den(hi(1),j,n) = cvmgt(-1.0, 0.0,
     $                    mask(hi(1)+1,j) .gt. 0)
                  enddo
               endif
               
            else
               print *,'UNKNOWN BC ON RIGHT FACE IN APPLYBC'
               stop
            endif
         enddo
      endif
c
c     The Bottom of the Grid
c
      if(cdir .eq. 1) then
         do n = 1, nc
            if(is_neumann(bct(n))) then
               do i = lo(1),hi(1)
                  u(i,lo(2)-1,n) = cvmgt(
     $                 u(i,lo(2),n),
     $                 u(i,lo(2)-1,n),
     $                 mask(i,lo(2)-1) .gt. 0)
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1),hi(1)
                     den(i,lo(2),n)   = 1.0
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: interior part of bottom side
               j = lo(2)
               do i=lo(1)+1,hi(1)-1
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  tander(i,j-1,n) = innder
               enddo
c :::::: now end points
c :::::: left bottom side
               i = lo(1)
               if( maskm(-1 + i,j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(1 + i,j,n) - U(2 + i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               tander(i,j-1,n) = innder
c :::::: right bottom side
               i = hi(1)
               if( maskp(1 + i,j).gt.0 ) then
                  innder = (U(-2 + i,j,n) - 4*U(-1 + i,j,n) + 3*U(i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               tander(i,j-1,n) = innder
            else if (is_dirichlet(bct(n))) then
               do m=0,leny
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(2)
               xInt = - 0.5
               call polyInterpCoeff(xInt, x, leny+2, coef)
               if ( flagbc .eq. 1 ) then
                  do i = lo(1), hi(1)
                     u(i,lo(2)-1,n) = cvmgt(
     $                    bcval(i,lo(2)-1,n)*coef(-1),
     $                    u(i,lo(2)-1,n),
     $                    mask(i,lo(2)-1) .gt. 0)
                  enddo
               else
                  do i = lo(1), hi(1)
                     u(i,lo(2)-1,n) = cvmgt(
     $                    0.0,
     $                    u(i,lo(2)-1,n),
     $                    mask(i,lo(2)-1) .gt. 0)
                  enddo
               endif
               do m = 0, leny
                  do i = lo(1), hi(1)
                     u(i, lo(2)-1, n) = cvmgt(
     $                    u(i, lo(2)-1,n)
     $                    + u(i, lo(2)+m,n)*coef(m),
     $                    u(i, lo(2)-1, n),
     $                    mask(i, lo(2)-1) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i, lo(2),n)   = cvmgt(coef(0), 0.0,
     $                    mask(i, lo(2)-1) .gt. 0)
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: interior part of bottom side
               edgloc = 0.
               j = lo(2)
               do i=lo(1)+1,hi(1)-1
                  if(  mask(-1 + i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-U(-1 + i,-1 + j,n) + U(1 + i,-1 + j,n))/(2.*hx)
                  else if(  mask(i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0.and.mask(2 + i,-1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(i,-1 + j,n) + 4*U(1 + i,-1 + j,n) - U(2 + i,-1 + j,n))/(2.*hx)
                  else if(  mask(-2 + i,-1 + j).eq.0.and.mask(-1 + i,-1 + j).eq.0.and.mask(i,-1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(-2 + i,-1 + j,n) - 4*U(-1 + i,-1 + j,n) + 3*U(i,-1 + j,n))/(2.*hx)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i,j-1,n)
                     else
                        outder = 0.
                     endif
                  endif
                  innloc = 0.5
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)

                  lambda = (edgloc-outloc)/(innloc-outloc)
                  tander(i,j-1,n) = lambda*innder+(1-lambda)*outder
               enddo
c :::::: now end points
c :::::: left bottom side
               i = lo(1)
               if(  mask(-1 + i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(-1 + i,-1 + j,n) + U(1 + i,-1 + j,n))/(2.*hx)
               else if( mask(i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0.and.mask(2 + i,-1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(i,-1 + j,n) + 4*U(1 + i,-1 + j,n) - U(2 + i,-1 + j,n))/(2.*hx)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i,j-1,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskm(-1 + i,j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(1 + i,j,n) - U(2 + i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i,j-1,n) = lambda*innder+(1-lambda)*outder
c :::::: right bottom side
               i = hi(1)
               if( mask(-1 + i,-1 + j).eq.0.and.mask(1 + i,-1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(-1 + i,-1 + j,n) + U(1 + i,-1 + j,n))/(2.*hx)
               else if( mask(-2 + i,-1 + j).eq.0.and.mask(-1 + i,-1 + j).eq.0.and.mask(i,-1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(-2 + i,-1 + j,n) - 4*U(-1 + i,-1 + j,n) + 3*U(i,-1 + j,n))/(2.*hx)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i,j-1,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskp(1 + i,j).gt.0 ) then
                  innder = (U(-2 + i,j,n) - 4*U(-1 + i,j,n) + 3*U(i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i,j-1,n) = lambda*innder+(1-lambda)*outder
               
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               
               do i = lo(1), hi(1)
                  u(i,lo(2)-1,n) = cvmgt(
     $                 -u(i,lo(2),n),
     $                 u(i,lo(2)-1,n),
     $                 mask(i,lo(2)-1) .gt. 0)
                  tander(i,lo(2)-1,n) = 0.0
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i,lo(2),n) = cvmgt(-1.0, 0.0,
     $                    mask(i,lo(2)-1) .gt. 0)
                  enddo
               endif
               
            else
               print *,'UNKNOWN BC ON BOTTOM FACE IN APPLYBC'
               stop
            endif
        enddo
      endif
c     
c     The top of the grid
c
      if (cdir .eq. 3) then
         do n = 1, nc
            if(is_neumann(bct(n))) then
               do i = lo(1), hi(1)
                  u(i,hi(2)+1,n) = cvmgt(
     $                 u(i,hi(2),n),
     $                 u(i,hi(2)+1,n),
     $                 mask(i,hi(2)+1) .gt. 0)
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i,hi(2),n)   = 1.0
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: interior part of top side
               j = hi(2)
               do i=lo(1)+1,hi(1)-1
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  tander(i,j+1,n) = innder
               enddo
c     :::::: now end points
c     :::::: left top side
               i = lo(1)
               if( maskm(-1 + i,j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(1 + i,j,n) - U(2 + i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               tander(i,j+1,n) = innder
c     :::::: right top side
               i = hi(1)
               if( maskp(1 + i,j).gt.0 ) then
                  innder = (U(-2 + i,j,n) - 4*U(-1 + i,j,n) + 3*U(i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               tander(i,j+1,n) = innder
            else if (is_dirichlet(bct(n))) then
               if ( bct(n) .eq. LO_REFLECT_ODD ) leny = 0
               do m=0,leny
                  x(m) = m + 0.5
               enddo
               x(-1) = - bcl/h(2)
               xInt = - 0.5
               call polyInterpCoeff(xInt, x, leny+2, coef)
               if ( flagbc .eq. 1 ) then
                  do i = lo(1), hi(1)
                     u(i,hi(2)+1,n) = cvmgt(
     $                    bcval(i,hi(2)+1,n)*coef(-1),
     $                    u(i,hi(2)+1,n),
     $                    mask(i,hi(2)+1) .gt. 0)
                  enddo
               else
                  do i = lo(1), hi(1)
                     u(i,hi(2)+1,n) = cvmgt(
     $                    0.0,
     $                    u(i,hi(2)+1,n),
     $                    mask(i,hi(2)+1) .gt. 0)
                  enddo
               endif
               do m = 0, leny
                  do i = lo(1), hi(1)
                     u(i, hi(2)+1,n) = cvmgt(
     $                    u(i,hi(2)+1,n)
     $                    + u(i, hi(2)-m,n)*coef(m),
     $                    u(i,hi(2)+1,n),
     $                    mask(i,hi(2)+1) .gt. 0)
                  enddo
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i,hi(2),n)   = cvmgt(coef(0), 0.0,
     $                    mask(i,hi(2)+1) .gt. 0)
                  enddo
               endif
c     :::::: now do the tangential derivative part
c     :::::: interior part of top side
               edgloc = 0.
               j = hi(2)
               do i=lo(1)+1,hi(1)-1
                  if(  mask(-1 + i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-U(-1 + i,1 + j,n) + U(1 + i,1 + j,n))/(2.*hx)
                  else if(  mask(i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0.and.mask(2 + i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(i,1 + j,n) + 4*U(1 + i,1 + j,n) - U(2 + i,1 + j,n))/(2.*hx)
                  else if( mask(-2 + i,1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0.and.mask(i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(-2 + i,1 + j,n) - 4*U(-1 + i,1 + j,n) + 3*U(i,1 + j,n))/(2.*hx)
                  else
                     outloc = x(-1)
                     if( flagbc .eq. 1) then
                        outder = exttd(i,j+1,n)
                     else
                        outder = 0.
                     endif
                  endif
                  innloc = 0.5
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
                  
                  lambda = (edgloc-outloc)/(innloc-outloc)
                  tander(i,j+1,n) = lambda*innder+(1-lambda)*outder
               enddo
c     :::::: now end points
c     :::::: left top side
               i = lo(1)
               if( mask(-1 + i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(-1 + i,1 + j,n) + U(1 + i,1 + j,n))/(2.*hx)
               else if( mask(i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0.and.mask(2 + i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (-3*U(i,1 + j,n) + 4*U(1 + i,1 + j,n) - U(2 + i,1 + j,n))/(2.*hx)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i,j+1,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskm(-1 + i,j).gt.0 ) then
                  innder = (-3*U(i,j,n) + 4*U(1 + i,j,n) - U(2 + i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i,j+1,n) = lambda*innder+(1-lambda)*outder
c     :::::: right top side
               i = hi(1)
               if( mask(-1 + i,1 + j).eq.0.and.mask(1 + i,1 + j).eq.0 ) then
                  outloc = -0.5
                  outder = (-U(-1 + i,1 + j,n) + U(1 + i,1 + j,n))/(2.*hx)
               else if(  mask(-2 + i,1 + j).eq.0.and.mask(-1 + i,1 + j).eq.0.and.mask(i,1 + j).eq.0 ) then
                     outloc = -0.5
                     outder = (U(-2 + i,1 + j,n) - 4*U(-1 + i,1 + j,n) + 3*U(i,1 + j,n))/(2.*hx)
               else
                  outloc = x(-1)
                  if( flagbc .eq. 1) then
                     outder = exttd(i,j+1,n)
                  else
                     outder = 0.
                  endif
               endif
               if( maskp(1 + i,j).gt.0 ) then
                  innder = (U(-2 + i,j,n) - 4*U(-1 + i,j,n) + 3*U(i,j,n))/(2.*hx)
               else
                  innder = (-U(-1 + i,j,n) + U(1 + i,j,n))/(2.*hx)
               endif
               innloc = 0.5
               lambda = (edgloc-outloc)/(innloc-outloc)
               tander(i,j+1,n) = lambda*innder+(1-lambda)*outder
               
            else if ( bct(n) .eq. LO_REFLECT_ODD ) then
               
               do i = lo(1), hi(1)
                  u(i,hi(2)+1,n) = cvmgt(
     $                 -u(i,hi(2),n),
     $                 u(i,hi(2)+1,n),
     $                 mask(i,hi(2)+1) .gt. 0)
                  tander(i,hi(2)+1,n) = 0.0
               enddo
               if ( flagden .eq. 1 ) then
                  do i = lo(1), hi(1)
                     den(i,hi(2),n) = cvmgt(-1.0, 0.0,
     $                    mask(i,hi(2)+1) .gt. 0)
                  enddo
               endif
               
            else
               print *,'UNKNOWN BC ON TOP FACE IN APPLYBC'
               stop
            endif
         enddo
      endif
c
      end


      subroutine FORT_GSRB (
     $     u, DIMS(u),
     $     rhs, DIMS(rhs),
     $     alpha, beta,
     $     a, DIMS(a),
     $     muX, DIMS(muX),
     $     muY, DIMS(muY),
     $     maskn,DIMS(maskn),
     $     fn, DIMS(fn),
     $     maske,DIMS(maske),
     $     fe, DIMS(fe),
     $     maskw,DIMS(maskw),
     $     fw, DIMS(fw),
     $     masks,DIMS(masks),
     $     fs, DIMS(fs),
     $     trandern,DIMS(trandern),
     $     trandere,DIMS(trandere),
     $     tranderw,DIMS(tranderw),
     $     tranders,DIMS(tranders),
     $     lo,hi,h,nc,phaseflag
     $     )

      REAL_T alpha, beta
      integer DIMDEC(u)
      integer DIMDEC(rhs)
      integer DIMDEC(a)
      integer DIMDEC(muX)
      integer DIMDEC(muY)
      integer DIMDEC(maskn)
      integer DIMDEC(maske)
      integer DIMDEC(maskw)
      integer DIMDEC(masks)
      integer DIMDEC(fn)
      integer DIMDEC(fe)
      integer DIMDEC(fw)
      integer DIMDEC(fs)
      integer DIMDEC(trandern)
      integer DIMDEC(trandere)
      integer DIMDEC(tranderw)
      integer DIMDEC(tranders)
      integer lo(BL_SPACEDIM), hi(BL_SPACEDIM)
      REAL_T h(BL_SPACEDIM)
      integer nc,phaseflag

      REAL_T u(DIMV(u),2)
      REAL_T rhs(DIMV(rhs),2)
      REAL_T a(DIMV(a),2)
      REAL_T muX(DIMV(muX))
      REAL_T muY(DIMV(muY))
      integer maskn(DIMV(maskn))
      integer maske(DIMV(maske))
      integer maskw(DIMV(maskw))
      integer masks(DIMV(masks))
      REAL_T fn(DIMV(fn),nc)
      REAL_T fe(DIMV(fe),nc)
      REAL_T fw(DIMV(fw),nc)
      REAL_T fs(DIMV(fs),nc)
      REAL_T trandern(DIMV(trandern),2)
      REAL_T trandere(DIMV(trandere),2)
      REAL_T tranderw(DIMV(tranderw),2)
      REAL_T tranders(DIMV(tranders),2)

      integer i,j,modx,mody
      integer istart,jstart
      REAL_T tauxxw,tauxxe,tauyyn,tauyys
      REAL_T tauxye,tauxyw,tauxyn,tauxys
      REAL_T hx,hy
      REAL_T dudxe,dudxw
      REAL_T dvdyn,dvdys
      REAL_T dudye,dudyw,dudyn,dudys
      REAL_T dvdxe,dvdxw,dvdxn,dvdxs
      REAL_T diagu, diagv
      REAL_T operu, operv
      logical False, True
      parameter( False=.false.)
      parameter( True=.true. )

      hx = h(1)
      hy = h(2)
      if( nc .ne. 2 ) then
         write(6,*)'FORT_GSRB: bad nc',nc
         stop
      endif

c ::: convert phase flag into modx and mody
      if(    phaseflag.eq.0 ) then
         modx = 0
         mody = 0
      elseif(phaseflag.eq.1) then
         modx = 1
         mody = 0
      elseif(phaseflag.eq.2) then
         modx = 0
         mody = 1
      elseif(phaseflag.eq.3) then
         modx = 1
         mody = 1
      else
         write(6,*)'FORT_GSRB: bad phaseflag', phaseflag
         stop
      endif


c ::: cases:
c     ::: 1) interior.  All normal derivative expressions
c     ::: 2) edge but not on corner.  check mask value to see if use outside
c     :::    value or use interpolations
c     ::: 3) corners.  Must check mask values to see if use outside values
c     :::    values or use interpolations ON TWO SIDES.  But must also
c     :::    check to see if need to do one-sided derivs for tangential
c     :::    derivatives.

c ::: case 1
      istart = lo(1)+1
      if( mod(istart,2) .ne. modx ) istart = istart+1
      jstart = lo(2)+1
      if( mod(jstart,2) .ne. mody ) jstart = jstart+1
      do j = jstart,hi(2)-1,2
         do i = istart,hi(1)-1,2
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudye =         (-U(i,-1 + j,1) + U(i,1 + j,1) - U(1 + i,-1 + j,1) + 
     -    U(1 + i,1 + j,1))/(4.*hy)
            dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)

            dvdxn =         (-U(-1 + i,j,2) - U(-1 + i,1 + j,2) + U(1 + i,j,2) + 
     -    U(1 + i,1 + j,2))/(4.*hx)
            dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  get diagonal element
            diagu = alpha*a(i,j,1) + 
     $           2*beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $             beta*(muY(i,j)+muY(i,j+1))/(hy**2)
            diagv = alpha*a(i,j,2) +
     $             beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $           2*beta*(muY(i,j)+muY(i,j+1))/(hy**2)
c     :::  evaluate operator
            operu =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            operv =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)
c     :::  do relaxation
            u(i,j,1) = (rhs(i,j,1)-operu)/diagu + u(i,j,1)
            u(i,j,2) = (rhs(i,j,2)-operv)/diagv + u(i,j,2)
         enddo
      enddo

c ::: case 2a) north side
c ::: Null
c ::: Null
      j = hi(2)
      if( mod(j,2) .eq. mody ) then
         do i = istart,hi(1)-1,2
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            if( maskn(i,1 + j).gt.0.or.maskn(1 + i,1 + j).gt.0 ) then
               dudye =         ((U(i,-2 + j,1) - 4*U(i,-1 + j,1) + 3*U(i,j,1))/(2.*hy) + 
     -    (U(1 + i,-2 + j,1) - 4*U(1 + i,-1 + j,1) + 3*U(1 + i,j,1))/(2.*hy))/
     -  2.
            else
               dudye =         (-U(i,-1 + j,1) + U(i,1 + j,1) - U(1 + i,-1 + j,1) + 
     -    U(1 + i,1 + j,1))/(4.*hy)
            endif
            if( maskn(-1 + i,1 + j).gt.0.or.maskn(i,1 + j).gt.0 ) then
               dudyw =         ((U(-1 + i,-2 + j,1) - 4*U(-1 + i,-1 + j,1) + 3*U(-1 + i,j,1))/
     -     (2.*hy) + (U(i,-2 + j,1) - 4*U(i,-1 + j,1) + 3*U(i,j,1))/(2.*hy))/
     -  2.
            else
               dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)
            endif

            dvdxn = trandern(i,j+1,2)
            dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  get diagonal element
            diagu = alpha*a(i,j,1) + 
     $           2*beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $             beta*(muY(i,j)+muY(i,j+1))/(hy**2)
            diagv = alpha*a(i,j,2) +
     $             beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $           2*beta*(muY(i,j)+muY(i,j+1))/(hy**2)
	    if( maskn(i,j+1).gt. 0) then
	       diagu = diagu - fn(i,j,1)*beta*muY(i,j+1)/hy**2
	       diagv = diagv - fn(i,j,2)*2*beta*muY(i,j+1)/hy**2
	    endif
c     :::  evaluate operator
            operu =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            operv =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)
c     :::  do relaxation
            u(i,j,1) = (rhs(i,j,1)-operu)/diagu + u(i,j,1)
            u(i,j,2) = (rhs(i,j,2)-operv)/diagv + u(i,j,2)
         enddo
      endif



c ::: case 2b south side
c ::: Null
c ::: Null

      j=lo(2)
      if( mod(j,2) .eq. mody ) then
         do i=istart,hi(1)-1,2
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            if( masks(i,-1 + j).gt.0.or.masks(1 + i,-1 + j).gt.0 ) then
               dudye =         ((-3*U(i,j,1) + 4*U(i,1 + j,1) - U(i,2 + j,1))/(2.*hy) + 
     -    (-3*U(1 + i,j,1) + 4*U(1 + i,1 + j,1) - U(1 + i,2 + j,1))/(2.*hy))/
     -  2.
            else
               dudye =         (-U(i,-1 + j,1) + U(i,1 + j,1) - U(1 + i,-1 + j,1) + 
     -    U(1 + i,1 + j,1))/(4.*hy)
            endif
            if( masks(-1 + i,-1 + j).gt.0.or.masks(i,-1 + j).gt.0 ) then
               dudyw =         ((-3*U(-1 + i,j,1) + 4*U(-1 + i,1 + j,1) - U(-1 + i,2 + j,1))/
     -     (2.*hy) + (-3*U(i,j,1) + 4*U(i,1 + j,1) - U(i,2 + j,1))/(2.*hy))/2.
            else
               dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)
            endif

            dvdxn =         (-U(-1 + i,j,2) - U(-1 + i,1 + j,2) + U(1 + i,j,2) + 
     -    U(1 + i,1 + j,2))/(4.*hx)
            dvdxs = tranders(i,j-1,2)
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  get diagonal element
            diagu = alpha*a(i,j,1) + 
     $           2*beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $             beta*(muY(i,j)+muY(i,j+1))/(hy**2)
            diagv = alpha*a(i,j,2) +
     $             beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $           2*beta*(muY(i,j)+muY(i,j+1))/(hy**2)
	    if( masks(i,j-1).gt. 0) then
	       diagu = diagu - fs(i,j,1)*beta*muY(i,j)/hy**2
	       diagv = diagv - fs(i,j,2)*2*beta*muY(i,j)/hy**2
	    endif
c     :::  evaluate operator
            operu =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            operv =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)
c     :::  do relaxation
            u(i,j,1) = (rhs(i,j,1)-operu)/diagu + u(i,j,1)
            u(i,j,2) = (rhs(i,j,2)-operv)/diagv + u(i,j,2)
         enddo
      endif


c ::: case 2c: west side story
c ::: Null
c ::: Null
      i=lo(1)
      if( mod(i,2) .eq. modx ) then
         do j=jstart,hi(2)-1,2
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx
            
            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy
            
            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudye =         (-U(i,-1 + j,1) + U(i,1 + j,1) - U(1 + i,-1 + j,1) + 
     -    U(1 + i,1 + j,1))/(4.*hy)
            dudyw = tranderw(i-1,j,1)
            
            if( maskw(-1 + i,j).gt.0.or.maskw(-1 + i,1 + j).gt.0 ) then
               dvdxn =         ((-3*U(i,j,2) + 4*U(1 + i,j,2) - U(2 + i,j,2))/(2.*hx) + 
     -    (-3*U(i,1 + j,2) + 4*U(1 + i,1 + j,2) - U(2 + i,1 + j,2))/(2.*hx))/
     -  2.
            else
               dvdxn =         (-U(-1 + i,j,2) - U(-1 + i,1 + j,2) + U(1 + i,j,2) + 
     -    U(1 + i,1 + j,2))/(4.*hx)
            endif
            if( maskw(-1 + i,-1 + j).gt.0.or.maskw(-1 + i,j).gt.0 ) then
               dvdxs =         ((-3*U(i,-1 + j,2) + 4*U(1 + i,-1 + j,2) - U(2 + i,-1 + j,2))/
     -     (2.*hx) + (-3*U(i,j,2) + 4*U(1 + i,j,2) - U(2 + i,j,2))/(2.*hx))/2.
            else
               dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            endif
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  get diagonal element
            diagu = alpha*a(i,j,1) + 
     $           2*beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $             beta*(muY(i,j)+muY(i,j+1))/(hy**2)
            diagv = alpha*a(i,j,2) +
     $             beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $           2*beta*(muY(i,j)+muY(i,j+1))/(hy**2)
	    if( maskw(i-1,j).gt. 0) then
	       diagu = diagu - fw(i,j,1)*2*beta*muX(i,j)/hx**2
	       diagv = diagv - fw(i,j,2)*  beta*muX(i,j)/hy**2
	    endif
c     :::  evaluate operator
            operu =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            operv =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)
c     :::  do relaxation
            u(i,j,1) = (rhs(i,j,1)-operu)/diagu + u(i,j,1)
            u(i,j,2) = (rhs(i,j,2)-operv)/diagv + u(i,j,2)
         enddo
      endif


c ::: case 2d: east side story
c ::: Null
c ::: Null
      i=hi(1)
      if( mod(i,2) .eq. modx ) then
         do j=jstart,hi(2)-1,2
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudye = trandere(i+1,j,1)
            dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)

            if( maske(1 + i,j).gt.0.or.maske(1 + i,1 + j).gt.0 ) then
               dvdxn =         ((U(-2 + i,j,2) - 4*U(-1 + i,j,2) + 3*U(i,j,2))/(2.*hx) + 
     -    (U(-2 + i,1 + j,2) - 4*U(-1 + i,1 + j,2) + 3*U(i,1 + j,2))/(2.*hx))/
     -  2.
            else
               dvdxn =         (-U(-1 + i,j,2) - U(-1 + i,1 + j,2) + U(1 + i,j,2) + 
     -    U(1 + i,1 + j,2))/(4.*hx)
            endif
            if( maske(1 + i,-1 + j).gt.0.or.maske(1 + i,j).gt.0 ) then
               dvdxs =         ((U(-2 + i,-1 + j,2) - 4*U(-1 + i,-1 + j,2) + 3*U(i,-1 + j,2))/
     -     (2.*hx) + (U(-2 + i,j,2) - 4*U(-1 + i,j,2) + 3*U(i,j,2))/(2.*hx))/
     -  2.
            else
               dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            endif
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  get diagonal element
            diagu = alpha*a(i,j,1) + 
     $           2*beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $             beta*(muY(i,j)+muY(i,j+1))/(hy**2)
            diagv = alpha*a(i,j,2) +
     $             beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $           2*beta*(muY(i,j)+muY(i,j+1))/(hy**2)
	    if( maske(i+1,j).gt. 0) then
	       diagu = diagu - fe(i,j,1)*2*beta*muX(i+1,j)/hx**2
	       diagv = diagv - fe(i,j,2)*  beta*muX(i+1,j)/hy**2
	    endif
c     :::  evaluate operator
            operu =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            operv =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)
c     :::  do relaxation
            u(i,j,1) = (rhs(i,j,1)-operu)/diagu + u(i,j,1)
            u(i,j,2) = (rhs(i,j,2)-operv)/diagv + u(i,j,2)
         enddo
      endif

c ::: case 3a: north-east corner, plus i, plus j
c ::: Null
c ::: Null
      j=hi(2)
      i=hi(1)
      if( mod(i,2).eq.modx .and. mod(j,2).eq.mody) then
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudye = trandere(i+1,j,1)
            if( maskn(-1 + i,1 + j).gt.0.or.maskn(i,1 + j).gt.0 ) then
               dudyw =         ((U(-1 + i,-2 + j,1) - 4*U(-1 + i,-1 + j,1) + 3*U(-1 + i,j,1))/
     -     (2.*hy) + (U(i,-2 + j,1) - 4*U(i,-1 + j,1) + 3*U(i,j,1))/(2.*hy))/
     -  2.
            else
               dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)
            endif

            dvdxn = trandern(i,j+1,2)
            if( maske(1 + i,-1 + j).gt.0.or.maske(1 + i,j).gt.0 ) then
               dvdxs =         ((U(-2 + i,-1 + j,2) - 4*U(-1 + i,-1 + j,2) + 3*U(i,-1 + j,2))/
     -     (2.*hx) + (U(-2 + i,j,2) - 4*U(-1 + i,j,2) + 3*U(i,j,2))/(2.*hx))/
     -  2.
            else
               dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            endif
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  get diagonal element
            diagu = alpha*a(i,j,1) + 
     $           2*beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $             beta*(muY(i,j)+muY(i,j+1))/(hy**2)
            diagv = alpha*a(i,j,2) +
     $             beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $           2*beta*(muY(i,j)+muY(i,j+1))/(hy**2)
	    if( maske(i+1,j).gt. 0) then
	       diagu = diagu - fe(i,j,1)*2*beta*muX(i+1,j)/hx**2
	       diagv = diagv - fe(i,j,2)*  beta*muX(i+1,j)/hy**2
	    endif
	    if( maskn(i,j+1).gt. 0) then
	       diagu = diagu - fn(i,j,1)*beta*muY(i,j+1)/hy**2
	       diagv = diagv - fn(i,j,2)*2*beta*muY(i,j+1)/hy**2
	    endif
c     :::  evaluate operator
            operu =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            operv =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)
c     :::  do relaxation
            u(i,j,1) = (rhs(i,j,1)-operu)/diagu + u(i,j,1)
            u(i,j,2) = (rhs(i,j,2)-operv)/diagv + u(i,j,2)
      endif


c ::: case 3b: south-east corner, plus i, minus j, plus i, plus j
c ::: Null
      j=lo(2)
      i=hi(1)
      
      if( mod(i,2).eq.modx .and. mod(j,2).eq.mody) then
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            dudye = trandere(i+1,j,1)
            if( masks(-1 + i,-1 + j).gt.0.or.masks(i,-1 + j).gt.0 ) then 
               dudyw =         ((-3*U(-1 + i,j,1) + 4*U(-1 + i,1 + j,1) - U(-1 + i,2 + j,1))/
     -     (2.*hy) + (-3*U(i,j,1) + 4*U(i,1 + j,1) - U(i,2 + j,1))/(2.*hy))/2.
            else
               dudyw =         (-U(-1 + i,-1 + j,1) + U(-1 + i,1 + j,1) - U(i,-1 + j,1) + 
     -    U(i,1 + j,1))/(4.*hy)
            endif

            if( maske(1 + i,j).gt.0.or.maske(1 + i,1 + j).gt.0 ) then
               dvdxn =         ((U(-2 + i,j,2) - 4*U(-1 + i,j,2) + 3*U(i,j,2))/(2.*hx) + 
     -    (U(-2 + i,1 + j,2) - 4*U(-1 + i,1 + j,2) + 3*U(i,1 + j,2))/(2.*hx))/
     -  2. 
            else
               dvdxn =         (-U(-1 + i,j,2) - U(-1 + i,1 + j,2) + U(1 + i,j,2) + 
     -    U(1 + i,1 + j,2))/(4.*hx)
            endif
            dvdxs = tranders(i,j-1,2)
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  get diagonal element
            diagu = alpha*a(i,j,1) + 
     $           2*beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $             beta*(muY(i,j)+muY(i,j+1))/(hy**2)
            diagv = alpha*a(i,j,2) +
     $             beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $           2*beta*(muY(i,j)+muY(i,j+1))/(hy**2)
	    if( maske(i+1,j).gt. 0) then
	       diagu = diagu - fe(i,j,1)*2*beta*muX(i+1,j)/hx**2
	       diagv = diagv - fe(i,j,2)*  beta*muX(i+1,j)/hy**2
	    endif
	    if( masks(i,j-1).gt. 0) then
	       diagu = diagu - fs(i,j,1)*beta*muY(i,j)/hy**2
	       diagv = diagv - fs(i,j,2)*2*beta*muY(i,j)/hy**2
	    endif
c     :::  evaluate operator
            operu =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            operv =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)
c     :::  do relaxation
            u(i,j,1) = (rhs(i,j,1)-operu)/diagu + u(i,j,1)
            u(i,j,2) = (rhs(i,j,2)-operv)/diagv + u(i,j,2)
      endif

c ::: case 3c: south-west corner, minus i, minus j
c ::: Null
      j=lo(2)
      i=lo(1)

      if( mod(i,2).eq.modx .and. mod(j,2).eq.mody) then
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            if( masks(i,-1 + j).gt.0.or.masks(1 + i,-1 + j).gt.0 ) then
               dudye =         ((-3*U(i,j,1) + 4*U(i,1 + j,1) - U(i,2 + j,1))/(2.*hy) + 
     -    (-3*U(1 + i,j,1) + 4*U(1 + i,1 + j,1) - U(1 + i,2 + j,1))/(2.*hy))/
     -  2.
            else
               dudye =         (-U(i,-1 + j,1) + U(i,1 + j,1) - U(1 + i,-1 + j,1) + 
     -    U(1 + i,1 + j,1))/(4.*hy)
            endif
            dudyw = tranderw(i-1,j,1)

            if( maskw(-1 + i,j).gt.0.or.maskw(-1 + i,1 + j).gt.0 ) then
               dvdxn =         ((-3*U(i,j,2) + 4*U(1 + i,j,2) - U(2 + i,j,2))/(2.*hx) + 
     -    (-3*U(i,1 + j,2) + 4*U(1 + i,1 + j,2) - U(2 + i,1 + j,2))/(2.*hx))/
     -  2.
            else
               dvdxn =         (-U(-1 + i,j,2) - U(-1 + i,1 + j,2) + U(1 + i,j,2) + 
     -    U(1 + i,1 + j,2))/(4.*hx)
            endif
            dvdxs = tranders(i,j-1,2)
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  get diagonal element
            diagu = alpha*a(i,j,1) + 
     $           2*beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $             beta*(muY(i,j)+muY(i,j+1))/(hy**2)
            diagv = alpha*a(i,j,2) +
     $             beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $           2*beta*(muY(i,j)+muY(i,j+1))/(hy**2)
	    if( maskw(i-1,j).gt. 0) then
	       diagu = diagu - fw(i,j,1)*2*beta*muX(i,j)/hx**2
	       diagv = diagv - fw(i,j,2)*  beta*muX(i,j)/hy**2
	    endif
	    if( masks(i,j-1).gt. 0) then
	       diagu = diagu - fs(i,j,1)*beta*muY(i,j)/hy**2
	       diagv = diagv - fs(i,j,2)*2*beta*muY(i,j)/hy**2
	    endif
c     :::  evaluate operator
            operu =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            operv =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)
c     :::  do relaxation
            u(i,j,1) = (rhs(i,j,1)-operu)/diagu + u(i,j,1)
            u(i,j,2) = (rhs(i,j,2)-operv)/diagv + u(i,j,2)
      endif



c ::: case 3d: north-west corner, minus i, plus j
c ::: Null
      j=hi(2)
      i=lo(1)
      if( mod(i,2).eq.modx .and. mod(j,2).eq.mody) then
c     :::  get expressions for appropriate derivatives
            dudxe = (-U(i,j,1) + U(1 + i,j,1))/hx
            dudxw = (-U(-1 + i,j,1) + U(i,j,1))/hx

            dvdyn = (-U(i,j,2) + U(i,1 + j,2))/hy
            dvdys = (-U(i,-1 + j,2) + U(i,j,2))/hy

            dudyn = (-U(i,j,1) + U(i,1 + j,1))/hy
            dudys = (-U(i,-1 + j,1) + U(i,j,1))/hy
            if( maskn(i,1 + j).gt.0.or.maskn(1 + i,1 + j).gt.0 ) then
               dudye =         ((U(i,-2 + j,1) - 4*U(i,-1 + j,1) + 3*U(i,j,1))/(2.*hy) + 
     -    (U(1 + i,-2 + j,1) - 4*U(1 + i,-1 + j,1) + 3*U(1 + i,j,1))/(2.*hy))/
     -  2.
            else
               dudye =         (-U(i,-1 + j,1) + U(i,1 + j,1) - U(1 + i,-1 + j,1) + 
     -    U(1 + i,1 + j,1))/(4.*hy)
            endif
            dudyw = tranderw(i-1,j,1)

            dvdxn = trandern(i,j+1,2)
            if( maskw(-1 + i,-1 + j).gt.0.or.maskw(-1 + i,j).gt.0 ) then
               dvdxs =         ((-3*U(i,-1 + j,2) + 4*U(1 + i,-1 + j,2) - U(2 + i,-1 + j,2))/
     -     (2.*hx) + (-3*U(i,j,2) + 4*U(1 + i,j,2) - U(2 + i,j,2))/(2.*hx))/2.
            else
               dvdxs =         (-U(-1 + i,-1 + j,2) - U(-1 + i,j,2) + U(1 + i,-1 + j,2) + 
     -    U(1 + i,j,2))/(4.*hx)
            endif
            dvdxe = (-U(i,j,2) + U(1 + i,j,2))/hx
            dvdxw = (-U(-1 + i,j,2) + U(i,j,2))/hx
c     :::  get diagonal element
            diagu = alpha*a(i,j,1) + 
     $           2*beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $             beta*(muY(i,j)+muY(i,j+1))/(hy**2)
            diagv = alpha*a(i,j,2) +
     $             beta*(muX(i,j)+muX(i+1,j))/(hx**2)+
     $           2*beta*(muY(i,j)+muY(i,j+1))/(hy**2)
	    if( maskw(i-1,j).gt. 0) then
	       diagu = diagu - fw(i,j,1)*2*beta*muX(i,j)/hx**2
	       diagv = diagv - fw(i,j,2)*  beta*muX(i,j)/hy**2
	    endif
	    if( maskn(i,j+1).gt. 0) then
	       diagu = diagu - fn(i,j,1)*beta*muY(i,j+1)/hy**2
	       diagv = diagv - fn(i,j,2)*2*beta*muY(i,j+1)/hy**2
	    endif
c     :::  evaluate operator
            operu =         -((beta*(hy*(-2*dudxw*muX(i,j) + 2*dudxe*muX(1 + i,j)) + 
     -         hx*(-((dudys + dvdxs)*muY(i,j)) + (dudyn + dvdxn)*muY(i,1 + j))
     -         ))/(hx*hy)) + alpha*a(i,j,1)*u(i,j,1)
            operv =         -((beta*(hy*(-((dudyw + dvdxw)*muX(i,j)) + 
     -            (dudye + dvdxe)*muX(1 + i,j)) + 
     -         hx*(-2*dvdys*muY(i,j) + 2*dvdyn*muY(i,1 + j))))/(hx*hy)) + 
     -  alpha*a(i,j,2)*u(i,j,2)
c     :::  do relaxation
            u(i,j,1) = (rhs(i,j,1)-operu)/diagu + u(i,j,1)
            u(i,j,2) = (rhs(i,j,2)-operv)/diagv + u(i,j,2)
      endif


      return
      end
