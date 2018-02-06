#ifndef AMREX_MLPOISSON_H_
#define AMREX_MLPOISSON_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

// del dot grad phi

class MLPoisson
    : public MLCellLinOp
{
public:

    MLPoisson () {}
    MLPoisson (const Vector<Geometry>& a_geom,
               const Vector<BoxArray>& a_grids,
               const Vector<DistributionMapping>& a_dmap,
               const LPInfo& a_info = LPInfo(),
               const Vector<FactoryType const*>& a_factory = {});
    virtual ~MLPoisson ();

    MLPoisson (const MLPoisson&) = delete;
    MLPoisson (MLPoisson&&) = delete;
    MLPoisson& operator= (const MLPoisson&) = delete;
    MLPoisson& operator= (MLPoisson&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FactoryType const*>& a_factory = {});

protected:

    virtual void prepareForSolve () final;
    virtual bool isSingular (int amrlev) const final { return m_is_singular[amrlev]; }
    virtual bool isBottomSingular () const final { return m_is_singular[0]; }
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const final;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only=0) const final;

    virtual Real getAScalar () const final { return  0.0; }
    virtual Real getBScalar () const final { return -1.0; }
    virtual MultiFab const* getACoeffs (int amrlev, int mglev) const final { return nullptr; }
    virtual std::array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int amrlev, int mglev) const final
        { return { AMREX_D_DECL(nullptr,nullptr,nullptr)}; }

    virtual std::unique_ptr<MLLinOp> makeNLinOp (int grid_size) const final;

private:

    Vector<int> m_is_singular;
};

}

#endif
