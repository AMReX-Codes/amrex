#ifndef AMREX_ML_CURL_CURL_H_
#define AMREX_ML_CURL_CURL_H_
#include <AMReX_Config.H>

#include <AMReX_MLLinOp.H>
#include <AMReX_MLCurlCurl_K.H>

namespace amrex {

/**
 * \brief curl (alpha curl E) + beta E = rhs
 *
 * Here E is an Array of 3 MultiFabs on staggered grid, alpha is a positive
 * scalar, and beta is a non-negative scalar.
 *
 * It's the caller's responsibility to make sure rhs has consistent nodal
 * data. If needed, one could call prepareRHS for this.
 *
 * The smoother is based on the 4-color Gauss-Seidel smoother of Li
 * et. al. 2020.  "An Efficient Preconditioner for 3-D Finite Difference
 * Modeling of the Electromagnetic Diffusion Process in the Frequency
 * Domain", IEEE Transactions on Geoscience and Remote Sensing, 58, 500-509.
 */
class MLCurlCurl
    : public MLLinOpT<Array<MultiFab,3> >
{
public:
    using MF        = Array<MultiFab,3>;
    using RT        = typename MLLinOpT<MF>::RT;
    using BCType    = typename MLLinOpT<MF>::BCType;
    using BCMode    = typename MLLinOpT<MF>::BCMode;
    using StateMode = typename MLLinOpT<MF>::StateMode;
    using Location  = typename MLLinOpT<MF>::Location;

    MLCurlCurl () = default;
    MLCurlCurl (const Vector<Geometry>& a_geom,
                const Vector<BoxArray>& a_grids,
                const Vector<DistributionMapping>& a_dmap,
                const LPInfo& a_info = LPInfo());

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo());

    void setScalars (RT a_alpha, RT a_beta) noexcept;

    //! Synchronize RHS on nodal points. If the user can guarantee it, this
    //! function does not need to be called.
    void prepareRHS (Vector<MF*> const& rhs) const;

    void setDirichletNodesToZero (int amrlev, int mglev, MF& a_mf) const override;

    [[nodiscard]] std::string name () const override {
        return std::string("curl of curl");
    }

    void setLevelBC (int amrlev, const MF* levelbcdata,
                     const MF* robinbc_a = nullptr,
                     const MF* robinbc_b = nullptr,
                     const MF* robinbc_f = nullptr) override;

    void restriction (int amrlev, int cmglev, MF& crse, MF& fine) const override;

    void interpolation (int amrlev, int fmglev, MF& fine, const MF& crse) const override;

    void apply (int amrlev, int mglev, MF& out, MF& in, BCMode bc_mode,
                StateMode s_mode, const MLMGBndryT<MF>* bndry=nullptr) const override;

    void smooth (int amrlev, int mglev, MF& sol, const MF& rhs,
                 bool skip_fillboundary=false) const override;

    void solutionResidual (int amrlev, MF& resid, MF& x, const MF& b,
                           const MF* crse_bcdata=nullptr) override;

    void correctionResidual (int amrlev, int mglev, MF& resid, MF& x,
                             const MF& b, BCMode bc_mode,
                             const MF* crse_bcdata=nullptr) override;

    void prepareForSolve () override;

    [[nodiscard]] bool isSingular (int /*amrlev*/) const override { return false; }
    [[nodiscard]] bool isBottomSingular () const override { return false; }

    RT xdoty (int amrlev, int mglev, const MF& x, const MF& y, bool local) const override;

    [[nodiscard]] RT normInf (int amrlev, MF const& mf, bool local) const override;

    void averageDownAndSync (Vector<MF>& sol) const override;

    [[nodiscard]] IntVect getNGrowVectRestriction () const override {
        return IntVect(1);
    }

    void make (Vector<Vector<MF> >& mf, IntVect const& ng) const override;

    [[nodiscard]] MF make (int amrlev, int mglev, IntVect const& ng) const override;

    [[nodiscard]] MF makeAlias (MF const& mf) const override;

    [[nodiscard]] MF makeCoarseMG (int amrlev, int mglev, IntVect const& ng) const override;

    [[nodiscard]] MF makeCoarseAmr (int famrlev, IntVect const& ng) const override;

// public for cuda

    void smooth4 (int amrlev, int mglev, MF& sol, MF const& rhs, int color) const;

    void compresid (int amrlev, int mglev, MF& resid, MF const& b) const;

    void applyPhysBC (int amrlev, int mglev, MultiFab& mf, CurlCurlStateType type) const;

private:

    void applyBC (int amrlev, int mglev, MF& in, CurlCurlStateType type) const;

    [[nodiscard]] iMultiFab const& getDotMask (int amrlev, int mglev, int idim) const;

    [[nodiscard]] CurlCurlDirichletInfo getDirichletInfo (int amrlev, int mglev) const;
    [[nodiscard]] CurlCurlSymmetryInfo getSymmetryInfo (int amrlev, int mglev) const;

    RT m_alpha = std::numeric_limits<RT>::lowest();
    RT m_beta  = std::numeric_limits<RT>::lowest();

    Array<IntVect,3> m_etype
#if (AMREX_SPACEDIM == 3)
         {IntVect(0,1,1), IntVect(1,0,1), IntVect(1,1,0)};
#else
         {IntVect(0,1), IntVect(1,0), IntVect(1,1)};
#endif

    mutable Vector<Vector<Array<std::unique_ptr<iMultiFab>,3>>> m_dotmask;
    static constexpr int m_ncomp = 1;
    Vector<Vector<std::unique_ptr<Gpu::DeviceScalar
                                  <LUSolver<AMREX_SPACEDIM*2,RT>>>>> m_lusolver;
};

}

#endif
