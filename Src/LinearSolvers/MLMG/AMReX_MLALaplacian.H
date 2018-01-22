#ifndef AMREX_MLALAPLACIAN_H_
#define AMREX_MLALAPLACIAN_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

// del dot grad phi

class MLALaplacian
    : public MLCellLinOp
{
public:

    MLALaplacian () {}
    MLALaplacian (const Vector<Geometry>& a_geom,
                  const Vector<BoxArray>& a_grids,
                  const Vector<DistributionMapping>& a_dmap,
                  const LPInfo& a_info = LPInfo());
    virtual ~MLALaplacian ();

    MLALaplacian (const MLALaplacian&) = delete;
    MLALaplacian (MLALaplacian&&) = delete;
    MLALaplacian& operator= (const MLALaplacian&) = delete;
    MLALaplacian& operator= (MLALaplacian&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo());

    void setScalars (Real a, Real b);
    void setACoeffs (int amrlev, const MultiFab& alpha);

protected:

    virtual void prepareForSolve () final;
    virtual bool isSingular (int amrlev) const final { return m_is_singular[amrlev]; }
    virtual bool isBottomSingular () const final { return m_is_singular[0]; }
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const final;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only=0) const final;
    virtual Real getBScalar () const final { return m_b_scalar; }

    virtual std::unique_ptr<MLLinOp> makeNLinOp (int grid_size) const final {
        amrex::Abort("MLALaplacian::makeNLinOp: Not implmented");
        return std::unique_ptr<MLLinOp>{};
    }

private:

    Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
    Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
    Vector<Vector<MultiFab> > m_a_coeffs;

    Vector<int> m_is_singular;

    //
    // functions
    //

    void averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a);
    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);
};

}

#endif
