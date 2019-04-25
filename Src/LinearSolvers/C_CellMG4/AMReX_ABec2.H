
#ifndef AMREX_ABec2_H_
#define AMREX_ABec2_H_

#include <AMReX_Array.H>
#include <AMReX_LinOp.H>
#include <AMReX_ABecLaplacian.H>

namespace amrex {

class ABec2
    : public ABecLaplacian
{
public:

  //! constructor for box array, boundary data and vector mesh spacing
  ABec2 (const BndryData& bd,
         const Real*      h)
    : ABecLaplacian(bd,h) {}

  virtual ~ABec2 () {}

  void altSmooth (MultiFab&       solnL,
                  const MultiFab& resL,
                  int             level,
                  int             redBlackFlag);

  virtual void smooth (MultiFab&       solnL,
		       const MultiFab& rhsL,
		       int             level,
		       LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);

  virtual void jacobi_smooth (MultiFab&       solnL,
			      const MultiFab& rhsL,
			      int             level,
			      LinOp::BC_Mode  bc_mode = LinOp::Inhomogeneous_BC);

protected:

  void altApplyBC (int  level,
                   bool local);

  virtual void Fsmooth_jacobi (MultiFab&       solnL,
			       const MultiFab& resL,
			       int             level);

  virtual void Fsmooth (MultiFab&       solnL,
			const MultiFab& resL,
			int             level,
			int             redBlackFlag);

  //! Disallow copy constructors
  ABec2 (const ABec2&);
  ABec2& operator= (const ABec2&);
};

}

#endif /*_ABec2_H_*/
