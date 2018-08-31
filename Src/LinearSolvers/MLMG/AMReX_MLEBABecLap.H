#ifndef AMREX_MLEBABECLAP_H_
#define AMREX_MLEBABECLAP_H_

#include <AMReX_EBFabFactory.H>
#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

// (alpha * a - beta * (del dot b grad)) phi

class MLEBABecLap
    : public MLCellLinOp
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

    void setEBDirichlet (int amrlev, const MultiFab& phi, const MultiFab& beta);

    virtual bool needsUpdate () const final {
        return (m_needs_update || MLCellLinOp::needsUpdate());
    }
    virtual void update () final;

    bool m_needs_update = true;

protected:

    virtual std::unique_ptr<FabFactory<FArrayBox> > makeFactory (int amrlev, int mglev) const final;

    virtual bool isCrossStencil () const { return false; }

    virtual void applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode,
                          const MLMGBndry* bndry=nullptr, bool skip_fillboundary=false) const final;
    virtual void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                        const MLMGBndry* bndry=nullptr) const final;

    virtual void prepareForSolve () final;
    virtual bool isSingular (int amrlev) const final { return m_is_singular[amrlev]; }
    virtual bool isBottomSingular () const final { return m_is_singular[0]; }
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const final;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only=0) const final;

    virtual void normalize (int amrlev, int mglev, MultiFab& mf) const final;

    virtual Real getAScalar () const final { return m_a_scalar; }
    virtual Real getBScalar () const final { return m_b_scalar; }
    virtual MultiFab const* getACoeffs (int amrlev, int mglev) const final
        { return &(m_a_coeffs[amrlev][mglev]); }
    virtual Array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int amrlev, int mglev) const final
        { return amrex::GetArrOfConstPtrs(m_b_coeffs[amrlev][mglev]); }

    virtual std::unique_ptr<MLLinOp> makeNLinOp (int grid_size) const final {
        amrex::Abort("MLABecLaplacian::makeNLinOp: Not implmented");
        return std::unique_ptr<MLLinOp>{};
    }

    virtual void restriction (int, int, MultiFab& crse, MultiFab& fine) const final;

    virtual void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final;

    virtual void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs) final;

private:

    Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
    Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
    Vector<Vector<MultiFab> > m_a_coeffs;
    Vector<Vector<Array<MultiFab,AMREX_SPACEDIM> > > m_b_coeffs;
    Vector<Vector<iMultiFab> > m_cc_mask;

    Vector<std::unique_ptr<MultiFab> > m_eb_phi;
    Vector<Vector<std::unique_ptr<MultiFab> > > m_eb_b_coeffs;

    Vector<int> m_is_singular;

    mutable int m_is_inhomog;

    //
    // functions
    //
    bool isEBDirichlet () const { return m_eb_phi[0] != nullptr; }

    void averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a,
                                        Vector<Array<MultiFab,AMREX_SPACEDIM> >& b,
                                        const Vector<MultiFab*>& b_eb);
    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);
};

}

#endif
