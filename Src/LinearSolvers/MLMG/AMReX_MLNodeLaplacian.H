#ifndef AMREX_ML_NODE_LAPLACIAN_H_
#define AMREX_ML_NODE_LAPLACIAN_H_

#include <AMReX_MLNodeLinOp.H>
#ifdef AMREX_USE_HYPRE
#include <AMReX_HypreNodeLap.H>
#endif

namespace amrex {

// del dot (sigma grah phi) = rhs
// where phi and rhs are nodal, and sigma is cell-centered.

class MLNodeLaplacian
    : public MLNodeLinOp
{
public :

    enum struct CoarseningStrategy : int { Sigma, RAP };

    MLNodeLaplacian () noexcept {}
    MLNodeLaplacian (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info = LPInfo(),
                     const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
#ifdef AMREX_USE_EB
    MLNodeLaplacian (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info,
                     const Vector<EBFArrayBoxFactory const*>& a_factory);
#endif
    virtual ~MLNodeLaplacian ();

    MLNodeLaplacian (const MLNodeLaplacian&) = delete;
    MLNodeLaplacian (MLNodeLaplacian&&) = delete;
    MLNodeLaplacian& operator= (const MLNodeLaplacian&) = delete;
    MLNodeLaplacian& operator= (MLNodeLaplacian&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

#ifdef AMREX_USE_EB
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory);
#endif

    virtual std::string name () const override { return std::string("MLNodeLaplacian"); }

    void setRZCorrection (bool rz) noexcept { m_is_rz = rz; }

    void setSigma (int amrlev, const MultiFab& a_sigma);

    void compDivergence (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel);

    void compRHS (const Vector<MultiFab*>& rhs, const Vector<MultiFab*>& vel,
                  const Vector<const MultiFab*>& rhnd,
                  const Vector<MultiFab*>& rhcc);

    void updateVelocity (const Vector<MultiFab*>& vel, const Vector<MultiFab const*>& sol) const;

    void compSyncResidualCoarse (MultiFab& sync_resid, const MultiFab& phi,
                                 const MultiFab& vold, const MultiFab* rhcc,
                                 const BoxArray& fine_grids, const IntVect& ref_ratio);

    void compSyncResidualFine (MultiFab& sync_resid, const MultiFab& phi, const MultiFab& vold,
                               const MultiFab* rhcc);

    void setGaussSeidel (bool flag) noexcept { m_use_gauss_seidel = flag; }
    void setHarmonicAverage (bool flag) noexcept { m_use_harmonic_average = flag; }

    void setCoarseningStrategy (CoarseningStrategy cs) noexcept { m_coarsening_strategy = cs; }

    virtual BottomSolver getDefaultBottomSolver () const final override {
        return (m_coarsening_strategy == CoarseningStrategy::RAP) ?
            BottomSolver::bicgcg : BottomSolver::bicgstab;
    }

    virtual void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final override;
    virtual void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final override;
    virtual void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs) final override;

    virtual void reflux (int crse_amrlev,
                         MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                         MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const final override;

    virtual void prepareForSolve () final override;
    virtual bool isSingular (int amrlev) const final override
        { return (amrlev == 0) ? m_is_bottom_singular : false; }
    virtual bool isBottomSingular () const final override { return m_is_bottom_singular; }
    virtual void applyBC (int amrlev, int mglev, MultiFab& phi, BCMode bc_mode, StateMode s_mode,
                          bool skip_fillboundary=false) const final override;
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final override;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const final override;
    virtual void normalize (int amrlev, int mglev, MultiFab& mf) const final override;

    virtual void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final override;

    virtual void getFluxes (const Vector<MultiFab*>& a_flux,
                            const Vector<MultiFab*>& a_sol) const final override;

    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);
    void averageDownCoeffsSameAmrLevel (int amrlev);

    void restrictInteriorNodes (int camrlev, MultiFab& crhs, MultiFab& frhs) const;

    void FillBoundaryCoeff (MultiFab& sigma, const Geometry& geom);

    void buildMasks ();

    void buildStencil ();

#ifdef AMREX_USE_EB
    void buildIntegral ();
#endif

#ifdef AMREX_USE_HYPRE
    virtual std::unique_ptr<HypreNodeLap> makeHypreNodeLap (int bottom_verbose) const override;

    void fillIJMatrix (MFIter const& mfi, Array4<HypreNodeLap::Int const> const& nid,
                       Array4<int const> const& owner,
                       Vector<HypreNodeLap::Int>& ncols, Vector<HypreNodeLap::Int>& rows,
                       Vector<HypreNodeLap::Int>& cols, Vector<Real>& mat) const;
#endif

private:

    int m_is_rz = 0;

#ifdef AMREX_USE_EB
    CoarseningStrategy m_coarsening_strategy = CoarseningStrategy::RAP;
#else
    CoarseningStrategy m_coarsening_strategy = CoarseningStrategy::Sigma;
#endif

    Vector<Vector<Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> > > m_sigma;
    Vector<Vector<std::unique_ptr<MultiFab> > > m_stencil;

#ifdef AMREX_USE_EB
    // they could be MultiCutFab
    Vector<std::unique_ptr<MultiFab> > m_integral;
    bool m_integral_built = false;
#endif

    bool m_use_gauss_seidel = true;
    bool m_use_harmonic_average = false;

    bool m_is_bottom_singular = false;
    bool m_masks_built = false;

    virtual void checkPoint (std::string const& file_name) const final;
};

}

#endif
