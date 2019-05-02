#ifndef AMREX_ML_ABECLAPLACIAN_H_
#define AMREX_ML_ABECLAPLACIAN_H_

#include <AMReX_MLCellABecLap.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

// (alpha * a - beta * (del dot b grad)) phi

class MLABecLaplacian
    : public MLCellABecLap
{
public:

    MLABecLaplacian () {}
    MLABecLaplacian (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info = LPInfo(),
                     const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
    virtual ~MLABecLaplacian ();

    MLABecLaplacian (const MLABecLaplacian&) = delete;
    MLABecLaplacian (MLABecLaplacian&&) = delete;
    MLABecLaplacian& operator= (const MLABecLaplacian&) = delete;
    MLABecLaplacian& operator= (MLABecLaplacian&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    void setScalars (Real a, Real b) noexcept;
    void setACoeffs (int amrlev, const MultiFab& alpha);
    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta);

    virtual bool needsUpdate () const override {
        return (m_needs_update || MLCellABecLap::needsUpdate());
    }
    virtual void update () override;

    virtual void prepareForSolve () override;
    virtual bool isSingular (int amrlev) const override { return m_is_singular[amrlev]; }
    virtual bool isBottomSingular () const override { return m_is_singular[0]; }
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final override;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const final override;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, Location /* loc */,
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

    void averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a,
                                        Vector<Array<MultiFab,AMREX_SPACEDIM> >& b);
    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);

    void applyMetricTermsCoeffs ();

    static void FFlux (Box const& box, Real const* dxinv, Real bscalar,
                       Array<FArrayBox const*, AMREX_SPACEDIM> const& bcoef,
                       Array<FArrayBox*,AMREX_SPACEDIM> const& flux,
                       FArrayBox const& sol, int face_only, int ncomp);

protected:

    bool m_needs_update = true;

    Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
    Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
    Vector<Vector<MultiFab> > m_a_coeffs;
    Vector<Vector<Array<MultiFab,AMREX_SPACEDIM> > > m_b_coeffs;

    Vector<int> m_is_singular;
};

}

#endif
