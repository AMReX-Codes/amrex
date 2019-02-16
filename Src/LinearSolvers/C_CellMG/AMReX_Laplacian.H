
#ifndef AMREX_LAPLACIAN_H_
#define AMREX_LAPLACIAN_H_

#include <AMReX_LinOp.H>

namespace amrex {

/**
  \brief A Laplacian tailors the description of a linear operator to apply
  the second-order central difference approximation to the laplacian
  differential operator on a cell-centered MultiFab of data

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

    /**
    * \brief constructor;
    *
    * \param bd
    * \param _h
    */
    Laplacian (const BndryData& bd,
               Real             _h);

    /**
    * \brief destructor
    */
    virtual ~Laplacian() override;

    /**
    * \brief Compute flux associated with the op
    *
    * \paramxflux
    * \paramyflux
    * \paramzflux)
    * \param in
    * \param bc_mode
    * \param sComp
    * \param dComp
    * \param nComp
    * \param bndComp
    */
    virtual void compFlux (AMREX_D_DECL(MultiFab &xflux, MultiFab &yflux, MultiFab &zflux),
			   MultiFab& in, const BC_Mode& bc_mode=LinOp::Inhomogeneous_BC,
			   int sComp=0, int dComp=0, int nComp=1, int bndComp=0) override;

    virtual Real norm (int nm = 0, int level = 0, const bool local = false) override;

protected:

    /**
    * \brief compute out=L(in) at level=level
    *
    * \param out
    * \param in
    * \param level
    */
    virtual void Fapply (MultiFab&       out,
                         const MultiFab& in,
                         int             level) override;
    virtual void Fapply (MultiFab&       out,
			 int             dst_comp,
                         const MultiFab& in,
			 int             src_comp,
			 int             num_comp,
                         int             level) override;

    /**
    * \brief apply GSRB smoother to improve residual to ``L(solnL)=rhsL``
    *
    * \param solnL
    * \param rhsL
    * \param level
    * \param rgbflag
    */
    virtual void Fsmooth (MultiFab&       solnL,
                          const MultiFab& rhsL,
                          int             level,
                          int             rgbflag) override;

    /**
    * \brief apply Jacobi smoother to improve residual to ``L(solnL)=rhsL``
    *
    * \param solnL
    * \param rhsL
    * \param level
    */
    virtual void Fsmooth_jacobi (MultiFab&       solnL,
                                 const MultiFab& rhsL,
                                 int             level) override;
};

}

#endif /*_LAPLACIAN_H_*/
