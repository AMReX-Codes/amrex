#ifndef AMREX_MLEBABECLAP_H_
#define AMREX_MLEBABECLAP_H_
#include <AMReX_Config.H>

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

    MLEBABecLap () = default;
    MLEBABecLap (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory,
                 int a_ncomp = 1);

    ~MLEBABecLap () override;

    MLEBABecLap (const MLEBABecLap&) = delete;
    MLEBABecLap (MLEBABecLap&&) = delete;
    MLEBABecLap& operator= (const MLEBABecLap&) = delete;
    MLEBABecLap& operator= (MLEBABecLap&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory,
                 int a_ncomp = 1);

    void setPhiOnCentroid ();

    void setScalars (Real a, Real b);
    void setACoeffs (int amrlev, const MultiFab& alpha);
    void setACoeffs (int amrlev, Real alpha);

    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta,
                     Location a_beta_loc);
    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta)
        {setBCoeffs (amrlev, beta, Location::FaceCenter);}

    void setBCoeffs (int amrlev, Real beta);
    void setBCoeffs (int amrlev, Vector<Real> const& beta);

    // Tells the solver that EB boundaries have Dirichlet bc's specified by "phi"
    void setEBDirichlet      (int amrlev, const MultiFab& phi, const MultiFab& beta);
    void setEBDirichlet      (int amrlev, const MultiFab& phi, Real beta);
    void setEBDirichlet      (int amrlev, const MultiFab& phi, Vector<Real> const& beta);

    // Tells the solver that EB boundaries have homogeneous Dirichlet bc's
    void setEBHomogDirichlet (int amrlev,                      const MultiFab& beta);
    void setEBHomogDirichlet (int amrlev,                      Real beta);
    void setEBHomogDirichlet (int amrlev,                      Vector<Real> const& beta);

    int getNComp () const override { return m_ncomp; }

    bool needsUpdate () const override {
        return (m_needs_update || MLCellABecLap::needsUpdate());
    }
    void update () override;

    std::unique_ptr<FabFactory<FArrayBox> > makeFactory (int amrlev, int mglev) const final;

    bool isCrossStencil () const override { return false; }

    void applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode, StateMode s_mode,
                          const MLMGBndry* bndry=nullptr, bool skip_fillboundary=false) const final;
    void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                        StateMode s_mode, const MLMGBndry* bndry=nullptr) const override;
    void compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                           MultiFab& sol, Location loc) const final;

    void prepareForSolve () override;
    bool isSingular (int amrlev) const override { return m_is_singular[amrlev]; }
    bool isBottomSingular () const override { return m_is_singular[0]; }
    void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const final;
    void FFlux (int amrlev, const MFIter& mfi,
                        const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, Location loc,
                        int face_only=0) const final;

    void normalize (int amrlev, int mglev, MultiFab& mf) const final;

    Real getAScalar () const final { return m_a_scalar; }
    Real getBScalar () const final { return m_b_scalar; }
    MultiFab const* getACoeffs (int amrlev, int mglev) const final
        { return &(m_a_coeffs[amrlev][mglev]); }
    Array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int amrlev, int mglev) const final
        { return amrex::GetArrOfConstPtrs(m_b_coeffs[amrlev][mglev]); }

    std::unique_ptr<MLLinOp> makeNLinOp (int /*grid_size*/) const final {
        amrex::Abort("MLABecLaplacian::makeNLinOp: Not implemented");
        return std::unique_ptr<MLLinOp>{};
    }

    void restriction (int, int, MultiFab& crse, MultiFab& fine) const final;

    void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final;

    void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs) final;

    void getEBFluxes (const Vector<MultiFab*>& a_flux,
                              const Vector<MultiFab*>& a_sol) const override;

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    [[nodiscard]] std::unique_ptr<Hypre> makeHypre (Hypre::Interface hypre_interface) const override;
#endif

#ifdef AMREX_USE_PETSC
    [[nodiscard]] std::unique_ptr<PETScABecLap> makePETSc () const override;
#endif

protected:

    int m_ncomp = 1;

    bool m_needs_update = true;

    Location m_beta_loc; // Location of coefficients: face centers or face centroids
    Location m_phi_loc;  // Location of solution variable: cell centers or cell centroids

    Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
    Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
    Vector<Vector<MultiFab> > m_a_coeffs;
    Vector<Vector<Array<MultiFab,AMREX_SPACEDIM> > > m_b_coeffs;
    Vector<Vector<iMultiFab> > m_cc_mask;

    Vector<std::unique_ptr<MultiFab> > m_eb_phi;
    Vector<Vector<std::unique_ptr<MultiFab> > > m_eb_b_coeffs;

    Vector<int> m_is_singular;

    mutable int m_is_eb_inhomog;

    //
    // functions
    //
    bool isEBDirichlet   () const noexcept { return m_eb_phi[0] != nullptr; }

    void averageDownCoeffsSameAmrLevel (int amrlev, Vector<MultiFab>& a,
                                        Vector<Array<MultiFab,AMREX_SPACEDIM> >& b,
                                        const Vector<MultiFab*>& b_eb);
    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);
};

}

#endif
