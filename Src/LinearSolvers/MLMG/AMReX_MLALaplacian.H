#ifndef AMREX_MLALAPLACIAN_H_
#define AMREX_MLALAPLACIAN_H_

#include <AMReX_MLCellABecLap.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

class MLALaplacian
    : public MLCellABecLap
{
public:

    MLALaplacian () {}
    MLALaplacian (const Vector<Geometry>& a_geom,
                  const Vector<BoxArray>& a_grids,
                  const Vector<DistributionMapping>& a_dmap,
                  const LPInfo& a_info = LPInfo(),
                  const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
    virtual ~MLALaplacian ();

    MLALaplacian (const MLALaplacian&) = delete;
    MLALaplacian (MLALaplacian&&) = delete;
    MLALaplacian& operator= (const MLALaplacian&) = delete;
    MLALaplacian& operator= (MLALaplacian&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    void setScalars (Real a, Real b) noexcept;
    void setACoeffs (int amrlev, const MultiFab& alpha);

    virtual void prepareForSolve () final override;
    virtual bool isSingular (int amrlev) const final override { return m_is_singular[amrlev]; }
    virtual bool isBottomSingular () const final override { return m_is_singular[0]; }
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final override;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const final override;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, Location /* loc */,
                        const int face_only=0) const final override;

    virtual void normalize (int marlve, int mglev, MultiFab& mf) const final override;

    virtual Real getAScalar () const final override { return m_a_scalar; }
    virtual Real getBScalar () const final override { return m_b_scalar; }
    virtual MultiFab const* getACoeffs (int amrlev, int mglev) const final override
        { return &(m_a_coeffs[amrlev][mglev]); }
    virtual Array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int amrlev, int mglev) const final override
        { return { AMREX_D_DECL(nullptr,nullptr,nullptr)}; }

    virtual std::unique_ptr<MLLinOp> makeNLinOp (int grid_size) const final override {
        amrex::Abort("MLALaplacian::makeNLinOp: Not implmented");
        return std::unique_ptr<MLLinOp>{};
    }

    void averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a);
    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);

private:

    Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
    Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
    Vector<Vector<MultiFab> > m_a_coeffs;

    Vector<int> m_is_singular;
};

}

#endif
