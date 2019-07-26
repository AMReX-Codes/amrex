#ifndef MCNodalLinOp_H
#define MCNodalLinOp_H

#include <AMReX_MLNodeLinOp.H>
#include <AMReX_MLNodeLap_F.H>
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

	void setNComp(int a_ncomp) {ncomp = a_ncomp;}
	void setCoeff(const Vector<Real> a_coeff) {coeff = a_coeff;}
	

protected:
	virtual void Fapply (int amrlev, int mglev,MultiFab& out,const MultiFab& in) const override ;
	void Diag (int amrlev, int mglev, MultiFab& diag);
	void averageDownCoeffs ();
	//
	// Virtual: you SHOULD override these functions
	//
protected:
	void Diagonal (bool recompute=false);
	//
	// Virtual: you CAN override these functions (but probably don't need to)
	//
	virtual void Fsmooth (int amrlev, int mglev, MultiFab& x,const MultiFab& b) const override;
	virtual void normalize (int amrlev, int mglev, MultiFab& mf) const override;
	virtual void reflux (int crse_amrlev, MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
			     MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const override;
	//
	// Virtual: you CAN'T override these functions (they take care of AMReX business) and you SHOULDN'T
	//
protected:
	virtual void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final;
	virtual void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const override final;
	virtual void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs, const MultiFab& fine_sol, const MultiFab& fine_rhs) final;
	virtual void prepareForSolve () final;
	virtual bool isSingular (int amrlev) const final {return (amrlev == 0) ? m_is_bottom_singular : false; }
	virtual bool isBottomSingular () const final { return m_is_bottom_singular; }
	virtual void applyBC (int amrlev, int mglev, MultiFab& phi, BCMode bc_mode, amrex::MLLinOp::StateMode /**/, bool skip_fillboundary=false) const final;
	virtual void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final;
public:
	virtual int getNGrow() const override final {return 2;}
	virtual int getNComp() const override final {return ncomp;}
	virtual void solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
				       const MultiFab* crse_bcdata=nullptr) override final;
	virtual void correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
					 BCMode bc_mode, const MultiFab* crse_bcdata=nullptr) override final;
	//
	// Private: internal machinery that you should never need to mess with
	//
private:
	void buildMasks ();
	void realFillBoundary(MultiFab &phi, const Geometry &geom) const;

	bool m_is_bottom_singular = false;
	bool m_masks_built = false;

	
	int ncomp = 1;
	Vector<Real> coeff = {1.0};
	
protected:
	bool m_diagonal_computed = false;
	amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab> > > m_diag;
};

#endif
