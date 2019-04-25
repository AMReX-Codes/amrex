#ifndef AMREX_MLEBABECLAP_H_
#define AMREX_MLEBABECLAP_H_

#include <AMReX_EBFabFactory.H>
#include <AMReX_MLCellABecLap.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

// (alpha * a - beta * (del dot b grad)) phi

class MLEBABecLap
    : public MLCellABecLap
{
public:

    MLEBABecLap () {}
    MLEBABecLap (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory);

    virtual ~MLEBABecLap ();

    MLEBABecLap (const MLEBABecLap&) = delete;
    MLEBABecLap (MLEBABecLap&&) = delete;
    MLEBABecLap& operator= (const MLEBABecLap&) = delete;
    MLEBABecLap& operator= (MLEBABecLap&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory);

    void setScalars (Real a, Real b);
    void setACoeffs (int amrlev, const MultiFab& alpha);
    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta);

    // Tells the solver that EB boundaries have Dirichlet bc's specified by "phi"
    void setEBDirichlet      (int amrlev, const MultiFab& phi, const MultiFab& beta);

    // Tells the solver that EB boundaries have homogeneous Dirichlet bc's
    void setEBHomogDirichlet (int amrlev,                      const MultiFab& beta);

    // Tells the solver to use the Higher Order Dirichlet stencil at EB boundaries 
    void setEBHODirichlet    ();

    virtual bool needsUpdate () const final override {
        return (m_needs_update || MLCellABecLap::needsUpdate());
    }
    virtual void update () final override;

    virtual std::unique_ptr<FabFactory<FArrayBox> > makeFactory (int amrlev, int mglev) const final override;

    virtual bool isCrossStencil () const override { return false; }

    virtual void applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode, StateMode s_mode,
                          const MLMGBndry* bndry=nullptr, bool skip_fillboundary=false) const final override;
    virtual void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                        StateMode s_mode, const MLMGBndry* bndry=nullptr) const final override;
    virtual void compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                           MultiFab& sol, Location loc) const final override;

    virtual void prepareForSolve () final override;
    virtual bool isSingular (int amrlev) const final override { return m_is_singular[amrlev]; }
    virtual bool isBottomSingular () const final override { return m_is_singular[0]; }
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final override;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const final override;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, Location loc,
                        const int face_only=0) const final override;

    virtual void normalize (int amrlev, int mglev, MultiFab& mf) const final override;

    virtual Real getAScalar () const final override { return m_a_scalar; }
    virtual Real getBScalar () const final override { return m_b_scalar; }
    virtual MultiFab const* getACoeffs (int amrlev, int mglev) const final override
        { return &(m_a_coeffs[amrlev][mglev]); }
    virtual Array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int amrlev, int mglev) const final override
        { return amrex::GetArrOfConstPtrs(m_b_coeffs[amrlev][mglev]); }

    virtual std::unique_ptr<MLLinOp> makeNLinOp (int grid_size) const final override {
        amrex::Abort("MLABecLaplacian::makeNLinOp: Not implmented");
        return std::unique_ptr<MLLinOp>{};
    }

    virtual void restriction (int, int, MultiFab& crse, MultiFab& fine) const final override;

    virtual void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final override;

    virtual void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs) final override;

#ifdef AMREX_USE_HYPRE
    virtual std::unique_ptr<Hypre> makeHypre (Hypre::Interface hypre_interface) const override;
#endif

#ifdef AMREX_USE_PETSC
    virtual std::unique_ptr<PETScABecLap> makePETSc () const override;
#endif

protected:

    bool m_needs_update = true;

private:

    Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
    Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
    Vector<Vector<MultiFab> > m_a_coeffs;
    Vector<Vector<Array<MultiFab,AMREX_SPACEDIM> > > m_b_coeffs;
    Vector<Vector<iMultiFab> > m_cc_mask;

    Vector<std::unique_ptr<MultiFab> > m_eb_phi;
    Vector<Vector<std::unique_ptr<MultiFab> > > m_eb_b_coeffs;

    Vector<int> m_is_singular;

    mutable int m_is_eb_inhomog;
    mutable int m_is_eb_ho_dir = false;

    //
    // functions
    //
    bool isEBDirichlet   () const noexcept { return m_eb_phi[0] != nullptr; }
    bool isHOEBDirichlet () const noexcept { return m_is_eb_ho_dir;}

    void averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a,
                                        Vector<Array<MultiFab,AMREX_SPACEDIM> >& b,
                                        const Vector<MultiFab*>& b_eb);
    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);
};

}

#endif
