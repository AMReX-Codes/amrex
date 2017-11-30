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
               const LPInfo& a_info = LPInfo());
    virtual ~MLPoisson ();

    MLPoisson (const MLPoisson&) = delete;
    MLPoisson (MLPoisson&&) = delete;
    MLPoisson& operator= (const MLPoisson&) = delete;
    MLPoisson& operator= (MLPoisson&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo());

protected:

    virtual void prepareForSolve () final;
    virtual bool isSingular (int amrlev) const final { return m_is_singular[amrlev]; }
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const final;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only=0) const final;
    virtual Real Anorm (int amrlev, int mglev) const final;
    virtual Real getBScalar () const final { return 1.0; }

private:

    Vector<int> m_is_singular;
};

}

#endif
