#ifndef AMREX_ML_ABECLAPLACIAN_H_
#define AMREX_ML_ABECLAPLACIAN_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

// (alpha * a - beta * (del dot b grad)) phi

class MLABecLaplacian
    : public MLCellLinOp
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

    void setScalars (Real a, Real b);
    void setACoeffs (int amrlev, const MultiFab& alpha);
    void setBCoeffs (int amrlev, const std::array<MultiFab const*,AMREX_SPACEDIM>& beta);

protected:

    virtual void prepareForSolve () final;
    virtual bool isSingular (int amrlev) const final { return m_is_singular[amrlev]; }
    virtual bool isBottomSingular () const final { return m_is_singular[0]; }
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const final;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only=0) const final;

    virtual void normalize (int amrlev, int mglev, MultiFab& mf) const final;

    virtual Real getAScalar () const final { return m_a_scalar; }
    virtual Real getBScalar () const final { return m_b_scalar; }
    virtual MultiFab const* getACoeffs (int amrlev, int mglev) const final
        { return &(m_a_coeffs[amrlev][mglev]); }
    virtual std::array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int amrlev, int mglev) const final
        { return amrex::GetArrOfConstPtrs(m_b_coeffs[amrlev][mglev]); }

    virtual std::unique_ptr<MLLinOp> makeNLinOp (int grid_size) const final {
        amrex::Abort("MLABecLaplacian::makeNLinOp: Not implmented");
        return std::unique_ptr<MLLinOp>{};
    }

private:

    Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
    Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
    Vector<Vector<MultiFab> > m_a_coeffs;
    Vector<Vector<std::array<MultiFab,AMREX_SPACEDIM> > > m_b_coeffs;

    Vector<int> m_is_singular;

    //
    // functions
    //

    void averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a,
                                        Vector<std::array<MultiFab,AMREX_SPACEDIM> >& b);
    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);

    void applyMetricTermsCoeffs ();
};

}

#endif
