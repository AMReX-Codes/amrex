
#ifndef _LAPLACIAN_H_
#define _LAPLACIAN_H_

//
// $Id: Laplacian.H,v 1.6 2000-10-02 20:51:16 lijewski Exp $
//

#include <LinOp.H>

//@Man:
/*@Memo:
  A Laplacian tailors the description of a linear operator to apply
  the second-order central difference approximation to the laplacian
  differential operator on a cell-centered MultiFab of data
  */        
/*@Doc:
  A Laplacian tailors a LinOp (a virtual base class for general linear
  operators) to compute the second-order central difference
  approximation to the laplacian differential operator,
  L(phi) = div(grad(phi)), on a cell-centered MultiFab, phi.  This class
  provides the necessary implementations for applying the operator and
  for smoothing a guessed solution for the linear
  equation system, L(phi) = rhs (where rhs is another cell-centered
  MultiFab).  These customizations are designed for 2D and 3D, with
  uniform, constant mesh spacing in all directions.

  Application of the operator and smoothing are
  "level" operations, and therefore must access "hidden" level data
  (such as boundary conditions, etc) as required, in concert with the
  semantics of the underlying LinOp defintions.  Accordingly, these
  implementations are "protected" and called only by the publically
  accessible members of the virtual base class, LinOp.  Note that these
  implementations may therefore assume that boundary (ghost) nodes and
  other pertinent information are filled properly by the base class
  functions prior to call.
        
  This class does NOT provide a copy constructor or assignment operator.
  */

class Laplacian
    :
    public LinOp
{
public:
    //
    //@ManDoc: constructor (merely falls through to LinOp constructor)
    //
    Laplacian (const BndryData& bd,
               Real             _h)
        :
        LinOp(bd,_h) {}
    //
    //@ManDoc: destructor
    //
    ~Laplacian() {}
    //
    //@ManDoc: Compute flux associated with the op
    //
    virtual void compFlux (D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
			   MultiFab& in, const BC_Mode& bc_mode);
    
    virtual Real norm (int nm = 0, int level = 0);
protected:
    //
    //@ManDoc: compute out=L(in) at level=level
    //
    virtual void Fapply (MultiFab&       out,
                         const MultiFab& in,
                         int             level);
    //
    //@ManDoc: apply smoother to improve residual to L(solnL)=rhsL
    //
    virtual void Fsmooth (MultiFab&       solnL,
                          const MultiFab& rhsL,
                          int             level,
                          int             rgbflag);
};

#endif /*_LAPLACIAN_H_*/
