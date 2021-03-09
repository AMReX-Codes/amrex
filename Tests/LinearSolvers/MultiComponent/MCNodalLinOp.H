//
// Tutorial:    MultiComponent Linear Solve
//
// File:        MCNodalLinOp.H
//
// Author:      Brandon Runnels
//              University of Colorado Colorado Springs
//              brunnels@uccs.edu
//              solids.uccs.edu
//
// Date:        September 3, 2019
// 
// Description: This file defines the multi-component linear operator
//
// See also:    ./MCNodalLinOp.cpp (for implementation of class)
//              ./main.cpp (for execution of the tutorial)
// 

#ifndef MCNodalLinOp_H
#define MCNodalLinOp_H

#include <AMReX_MLNodeLinOp.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BaseFab.H>

using namespace amrex;

class MCNodalLinOp : public MLNodeLinOp
{

public :
	MCNodalLinOp () {}
	virtual ~MCNodalLinOp ();
	void define (const Vector<Geometry>& a_geom, const Vector<BoxArray>& a_grids,
		     const Vector<DistributionMapping>& a_dmap,
		     const LPInfo& a_info = LPInfo(),
		     const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

	// Set the number of components that the operator will work with.
	// Default = 1
	void setNComp(int a_ncomp) {ncomp = a_ncomp;}
	// Set the coefficients for the operator. The input (a_coeff) must have
	// ncomp x ncomp ( x ncomp) components. It is the matrix of coefficients
	// used in the definition of the operator.
	void setCoeff(const Vector<Real> a_coeff) {coeff = a_coeff;}

protected:
	// Apply the operator, which is defined to be (in index notation)
	//    D(phi)_i = coeff_{ij} * Laplacian(phi_j)
	// where coeff_{ij} is stored as a vector, passed in using 
	// the setCoeff function
	virtual void Fapply (int amrlev, int mglev,MultiFab& out,const MultiFab& in) const override ;
	// Compute the diagonal of the operator
	void Diag (int amrlev, int mglev, MultiFab& diag);
protected:
	void Diagonal (bool recompute=false);
	// Run a Gauss-Siedel operation on the solution. 
	// (Note: this implementation uses just Diag and Fapply; no operator-specific code needed)
	virtual void Fsmooth (int amrlev, int mglev, MultiFab& x,const MultiFab& b) const override;
	// Divide out by the diagonal.
	// (Note: this implementation uses just Diag, no operator-specific code needed)
	virtual void normalize (int amrlev, int mglev, MultiFab& mf) const override;
	// Resolve the residual at the coarse/fine boundary.
	// (Note: unlike other nodal solvers, this is NOT operator-specific. Instead is just performs a
	//        restriction from the fine level to the coarse. MUST be used with CFStrategy::ghostnodes)
	virtual void reflux (int crse_amrlev, MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
			     MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const override;
protected:
	// Usual MG restriction
	virtual void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final;
	// Usual MG interpolation
	virtual void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const override final;
	// Usual averageDownSolutionRHS
	virtual void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs, const MultiFab& fine_sol, const MultiFab& fine_rhs) final;
	// Usual prepareForSolve
	virtual void prepareForSolve () final;
	// isSingular always returns 0 (singular operators not supported)
	virtual bool isSingular (int amrlev) const final {return (amrlev == 0) ? m_is_bottom_singular : false; }
	// isBottomSingular always returns 0 (singular operators not supported)
	virtual bool isBottomSingular () const final { return m_is_bottom_singular; }
	// Usual applyBC
	virtual void applyBC (int amrlev, int mglev, MultiFab& phi, BCMode bc_mode, amrex::MLLinOp::StateMode /**/, bool skip_fillboundary=false) const final;
	// Usual fixUpResidualMask
	virtual void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final;
public:
	// Return the number of ghost cells
	// Note: MUST return 2!
	virtual int getNGrow() const override final {return 2;}
	// Return the number of components
	virtual int getNComp() const override final {return ncomp;}
	// Overrides solutionResidual to add call to realFillBoundary
	virtual void solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
				       const MultiFab* crse_bcdata=nullptr) override final;
	// Overrides correctionResidual to add call to realFillBoundary
	virtual void correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
					 BCMode bc_mode, const MultiFab* crse_bcdata=nullptr) override final;
private:
	// Usual buildMasks
	void buildMasks ();

	bool m_is_bottom_singular = false;
	bool m_masks_built = false;
	
	int ncomp = 1;
	Vector<Real> coeff = {1.0};
	
protected:
	bool m_diagonal_computed = false;
	amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab> > > m_diag;
};

#endif
