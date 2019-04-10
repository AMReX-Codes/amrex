#ifndef AMREX_MLPOISSON_H_
#define AMREX_MLPOISSON_H_

#include <AMReX_MLCellABecLap.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

// del dot grad phi

class MLPoisson
    : public MLCellABecLap
{
public:

    MLPoisson () noexcept {}
    MLPoisson (const Vector<Geometry>& a_geom,
               const Vector<BoxArray>& a_grids,
               const Vector<DistributionMapping>& a_dmap,
               const LPInfo& a_info = LPInfo(),
               const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
    virtual ~MLPoisson ();

    MLPoisson (const MLPoisson&) = delete;
    MLPoisson (MLPoisson&&) = delete;
    MLPoisson& operator= (const MLPoisson&) = delete;
    MLPoisson& operator= (MLPoisson&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    virtual void prepareForSolve () final override;
    virtual bool isSingular (int amrlev) const final override { return m_is_singular[amrlev]; }
    virtual bool isBottomSingular () const final override { return m_is_singular[0]; }
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final override;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const final override;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, Location loc, const int face_only=0) const final override;

    virtual void normalize (int amrlev, int mglev, MultiFab& mf) const final override;

    virtual Real getAScalar () const final override { return  0.0; }
    virtual Real getBScalar () const final override { return -1.0; }
    virtual MultiFab const* getACoeffs (int amrlev, int mglev) const final override { return nullptr; }
    virtual Array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int amrlev, int mglev) const final override
        { return { AMREX_D_DECL(nullptr,nullptr,nullptr)}; }

    virtual std::unique_ptr<MLLinOp> makeNLinOp (int grid_size) const final override;

private:

    Vector<int> m_is_singular;
};

}

#endif
