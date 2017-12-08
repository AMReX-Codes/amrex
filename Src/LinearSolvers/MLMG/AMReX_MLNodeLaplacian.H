#ifndef AMREX_ML_NODE_LAPLACIAN_H_
#define AMREX_ML_NODE_LAPLACIAN_H_

#include <AMReX_MLNodeLinOp.H>

namespace amrex {

// del dot (sigma grah phi) = rhs
// where phi and rhs are nodal, and sigma is cell-centered.

class MLNodeLaplacian
    : public MLNodeLinOp
{
public :

    MLNodeLaplacian () {}
    MLNodeLaplacian (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info = LPInfo());
    virtual ~MLNodeLaplacian ();

    MLNodeLaplacian (const MLNodeLaplacian&) = delete;
    MLNodeLaplacian (MLNodeLaplacian&&) = delete;
    MLNodeLaplacian& operator= (const MLNodeLaplacian&) = delete;
    MLNodeLaplacian& operator= (MLNodeLaplacian&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo());

    void setSigma (int amrlev, const MultiFab& a_sigma);

    void compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel,
                  const Vector<const MultiFab*>& rhnd,
                  const Vector<const MultiFab*>& rhcc) const;

    void setGaussSeidel (bool flag) { m_use_gauss_seidel = flag; }

protected:

    virtual void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final;
    virtual void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final;

    virtual void prepareForSolve () final;
    virtual bool isSingular (int amrlev) const final { return false; }
    virtual void applyBC (int amrlev, int mglev, MultiFab& phi) const final;
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const final;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only=0) const final;
    virtual Real Anorm (int amrlev, int mglev) const final;
    virtual Real getBScalar () const final { return 1.0; }

private:

    Vector<Vector<std::array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> > > m_sigma;
    bool m_use_gauss_seidel = false;

    //
    // functions
    //
    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);
    void averageDownCoeffsSameAmrLevel (int amrlev);

    static void FillBoundaryCoeff (MultiFab& sigma, const Geometry& geom);
};

}

#endif
