#ifndef AMREX_ML_ABECLAPLACIAN_H_
#define AMREX_ML_ABECLAPLACIAN_H_

#include <AMReX_MLLinOp.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

class MLABecLaplacian
    : public MLLinOp
{
public:

    MLABecLaplacian () {}
    MLABecLaplacian (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap);
    ~MLABecLaplacian ();

    MLABecLaplacian (const MLABecLaplacian&) = delete;
    MLABecLaplacian (MLABecLaplacian&&) = delete;
    MLABecLaplacian& operator= (const MLABecLaplacian&) = delete;
    MLABecLaplacian& operator= (MLABecLaplacian&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap);

    void setScalars (Real a, Real b) { m_a_scalar = a; m_b_scalar = b; }
    void setACoeffs (int amrlev, const MultiFab& alpha);
    void setBCoeffs (int amrlev,
                     const std::array<MultiFab const*,AMREX_SPACEDIM>& beta);

    virtual void prepareForSolve () final;
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const final;
    virtual void FFlux (int amrlev, const MFIter& mfi, std::array<FArrayBox,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only=0) const final;

private:

    Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
    Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
    Vector<Vector<MultiFab> > m_a_coeffs;
    Vector<Vector<std::array<MultiFab,AMREX_SPACEDIM> > > m_b_coeffs;

    void averageDownCoeffs ();

    static void averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a,
                                               Vector<std::array<MultiFab,AMREX_SPACEDIM> >& b);
    void averageDownCoeffsToCoarseAmrLevel (int flev);
};

}

#endif
