#ifndef AMREX_MLNODEABECLAPLACIAN_H_
#define AMREX_MLNODEABECLAPLACIAN_H_
#include <AMReX_Config.H>

#include <AMReX_MLNodeLinOp.H>

namespace amrex {

// (alpha * a - beta * (del dot b grad)) phi = rhs
// a, phi and rhs are nodal. b is cell-centered.

class MLNodeABecLaplacian
    : public MLNodeLinOp
{
public:

    MLNodeABecLaplacian () = default;
    MLNodeABecLaplacian (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info = LPInfo(),
                         const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
    ~MLNodeABecLaplacian () override = default;

    MLNodeABecLaplacian (const MLNodeABecLaplacian&) = delete;
    MLNodeABecLaplacian (MLNodeABecLaplacian&&) = delete;
    MLNodeABecLaplacian& operator= (const MLNodeABecLaplacian&) = delete;
    MLNodeABecLaplacian& operator= (MLNodeABecLaplacian&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    std::string name () const override { return std::string("MLNodeABecLaplacian"); }

    void setScalars (Real a, Real b) {
        m_a_scalar = a;
        m_b_scalar = b;
    }

    void setACoeffs (int amrlev, Real a_acoef);
    void setACoeffs (int amrlev, const MultiFab& a_acoef);

    void setBCoeffs (int amrlev, Real a_bcoef);
    void setBCoeffs (int amrlev, const MultiFab& a_bcoef);

    void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const final;

    void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final;

    bool isSingular (int /*amrlev*/) const final { return false; }
    bool isBottomSingular () const final { return false; }

    void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final;
    void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final;
    void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                 const MultiFab& fine_sol, const MultiFab& fine_rhs) final;

    void reflux (int crse_amrlev,
                 MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                 MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const final;

    void prepareForSolve () final;

    [[nodiscard]] bool needsUpdate () const final { return m_needs_update; }

    void update () final;

    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);
    void averageDownCoeffsSameAmrLevel (int amrlev);

private:

    bool m_needs_update = true;

    Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
    Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
    Vector<Vector<MultiFab> > m_a_coeffs;
    Vector<Vector<MultiFab> > m_b_coeffs;
};

}

#endif
