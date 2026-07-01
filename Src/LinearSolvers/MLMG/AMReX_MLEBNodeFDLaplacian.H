#ifndef AMREX_MLEBNODEFDLAPLACIAN_H_
#define AMREX_MLEBNODEFDLAPLACIAN_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif
#include <AMReX_MLNodeLinOp.H>

#include <limits>

namespace amrex {

/**
 * \file AMReX_MLEBNodeFDLaplacian.H
 *
 * Nodal Laplacian with optional embedded boundaries, radial metrics,
 * and tensor-diagonal conductivity inputs; derived from `MLNodeLinOp`.
 */

// Although the class has EB in the name, it works for non-EB build too.
//
// del dot (sigma grad phi) = rhs, for non-RZ
// where phi and rhs are nodal multifab, and sigma is a tensor constant
// with only diagonal components.  The EB is assumed to be Dirichlet.
//
// del dot (sigma grad phi) - alpha/r^2 phi = rhs, for RZ where alpha is a
// scalar constant that is zero by default.
//
// New feature: sigma can also be a single-component cell-centered multifab.

/**
 * \brief Nodal finite-difference Laplacian with optional embedded boundaries.
 *
 * Handles both Cartesian and RZ geometries, supports tensor-diagonal conductivity
 * (`sigma`) provided either as scalars or MultiFabs, and enforces Dirichlet EB
 * values when factories expose geometry information.
 *
 * \ingroup amrex_eb
 */
class MLEBNodeFDLaplacian
    : public MLNodeLinOp
{
public:

    MLEBNodeFDLaplacian () = default;

#ifdef AMREX_USE_EB
    //! Construct directly from EB-aware factories.
    MLEBNodeFDLaplacian (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info,
                         const Vector<EBFArrayBoxFactory const*>& a_factory);
#endif

    //! Construct without EB data (works in non-EB builds as well).
    MLEBNodeFDLaplacian (const Vector<Geometry>& a_geom,
                         const Vector<BoxArray>& a_grids,
                         const Vector<DistributionMapping>& a_dmap,
                         const LPInfo& a_info);

    ~MLEBNodeFDLaplacian () override = default;

    MLEBNodeFDLaplacian (const MLEBNodeFDLaplacian&) = delete;
    MLEBNodeFDLaplacian (MLEBNodeFDLaplacian&&) = delete;
    MLEBNodeFDLaplacian& operator= (const MLEBNodeFDLaplacian&) = delete;
    MLEBNodeFDLaplacian& operator= (MLEBNodeFDLaplacian&&) = delete;

    //! Assign constant diagonal conductivity tensor `sigma`.
    void setSigma (Array<Real,AMREX_SPACEDIM> const& a_sigma) noexcept;

    //! Provide conductivity as a cell-centered MultiFab on AMR level \p amrlev.
    void setSigma (int amrlev, MultiFab const& a_sigma);

    //! Enable/disable RZ corrections (radial metrics).
    void setRZ (bool flag);

    //! Set the radial `alpha/r^2` term used when RZ is enabled.
    void setAlpha (Real a_alpha);

#ifdef AMREX_USE_EB

    //! Override phi on embedded boundaries (constant value).
    void setEBDirichlet (Real a_phi_eb);
    /**
     * \brief Override EB Dirichlet values using a coordinate-dependent callable.
     *
     * The callable must provide `Real operator()(Real x, Real y, Real z)` (with unused
     * coordinates ignored in lower dimensions) and returns the desired Dirichlet value.
     */
    template <typename F>
    requires (IsCallableR<Real,F,AMREX_D_DECL(Real,Real,Real)>::value)
    void setEBDirichlet (F const& f);

    /**
     * \brief Define the hierarchy using EB factories (captures cut-cell layout).
     *
     * \param a_geom    Per-level geometries.
     * \param a_grids   Per-level grids.
     * \param a_dmap    Distribution mappings.
     * \param a_info    LPInfo overrides.
     * \param a_factory EB factories per level.
     */
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory);

    //! EB-aware factory allocator for (\p amrlev,\p mglev).
    [[nodiscard]] std::unique_ptr<FabFactory<FArrayBox> > makeFactory (int amrlev, int mglev) const final;

    [[nodiscard]] bool scaleRHS (int amrlev, MultiFab* rhs) const final;

#endif

    /**
     * \brief Define the hierarchy when embedded boundaries are absent.
     *
     * \param a_geom  Per-level geometries.
     * \param a_grids Per-level grids.
     * \param a_dmap  Distribution mappings.
     * \param a_info  Optional LPInfo overrides.
     */
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info);

    [[nodiscard]] std::string name () const override { return std::string("MLEBNodeFDLaplacian"); }

    /**
     * \brief Restrict nodal data from fine to coarse MG levels.
     *
     * \param amrlev AMR level index.
     * \param cmglev Coarse MG level receiving the restricted data.
     * \param crse   Destination MultiFab on (\p amrlev,\p cmglev).
     * \param fine   Source MultiFab on the next finer MG level.
     */
    void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const final;
    /**
     * \brief Add the prolongation of coarse data onto the fine grid
     * (\p fine += prolong(\p crse)).
     *
     * \param amrlev AMR level index.
     * \param fmglev Fine MG level index.
     * \param fine   Destination fine MultiFab (updated in place).
     * \param crse   Source coarse MultiFab.
     */
    void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const final;

    [[nodiscard]] bool needsUpdate () const override {
        return (m_needs_update || MLNodeLinOp::needsUpdate());
    }
    void update () override;

    //! Finalize sigma/alpha/EB data prior to invoking MLMG.
    void prepareForSolve () final;
    //! Apply the nodal operator to \p in and write to \p out.
    void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    /**
     * \brief Perform the nodal smoother on (\p amrlev,\p mglev).
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param sol    Correction updated in place.
     * \param rhs    Right-hand side.
     */
    void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs) const final;
    //! Divide \p mf by the diagonal of the operator (used by CG-family bottom solvers).
    void normalize (int amrlev, int mglev, MultiFab& mf) const final;

    //! Adjust the residual mask \p resmsk to honor EB Dirichlet nodes (not yet implemented).
    void fixUpResidualMask (int amrlev, iMultiFab& resmsk) final;

    [[nodiscard]] bool isSingular (int) const final { return false; }
    [[nodiscard]] bool isBottomSingular () const final { return false; }

    /**
     * \brief Compute gradients of \p sol into \p grad.
     *
     * \param amrlev AMR level index.
     * \param grad   Destination gradients per direction.
     * \param sol    Solution supplying nodal values.
     */
    void compGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grad,
                   MultiFab& sol, Location /*loc*/) const override;

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    void fillIJMatrix (MFIter const& mfi,
                       Array4<HypreNodeLap::AtomicInt const> const& gid,
                       Array4<int const> const& lid,
                       HypreNodeLap::Int* ncols,
                       HypreNodeLap::Int* cols,
                       Real* mat) const override;

    void fillRHS (MFIter const& mfi,
                  Array4<int const> const& lid,
                  Real* rhs,
                  Array4<Real const> const& bfab) const override;
#endif

    //! Post-process the solution hierarchy.
    void postSolve (Vector<MultiFab*> const& sol) const override;

private:
    GpuArray<Real,AMREX_SPACEDIM> m_sigma{{AMREX_D_DECL(1_rt,1_rt,1_rt)}};
    Vector<Vector<std::unique_ptr<MultiFab>>> m_sigma_mf;
    bool m_has_sigma_mf = false;
    bool m_needs_update = true;
    Real m_s_phi_eb = std::numeric_limits<Real>::lowest();
    Vector<MultiFab> m_phi_eb;
    int m_rz = false;
    Real m_rz_alpha = 0._rt;

    void update_sigma ();
};

#ifdef AMREX_USE_EB

template <typename F>
requires (IsCallableR<Real,F,AMREX_D_DECL(Real,Real,Real)>::value)
void MLEBNodeFDLaplacian::setEBDirichlet (F const& f)
{
    m_phi_eb.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        auto const* factory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory[amrlev][0].get());
        if (factory) {
            Geometry const& geom = m_geom[amrlev][0];
            auto const problo = geom.ProbLoArray();
            auto const cellsize = geom.CellSizeArray();
            if (m_phi_eb[amrlev].empty()) {
                m_phi_eb[amrlev].define(amrex::convert(m_grids[amrlev][0],IntVect(1)),
                                        m_dmap[amrlev][0], 1, 1);
                m_phi_eb[amrlev].setVal(0.0);
            }
            auto const& flags = factory->getMultiEBCellFlagFab();
            auto const& levset = factory->getLevelSet();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(m_phi_eb[amrlev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& ndbx = mfi.growntilebox();
                const auto& flag = flags[mfi];
                if (flag.getType() != FabType::regular) {
                    Array4<Real const> const lstarr = levset.const_array(mfi);
                    Array4<Real> const& phi = m_phi_eb[amrlev].array(mfi);
                    AMREX_HOST_DEVICE_FOR_3D(ndbx, i, j, k,
                    {
                        if (lstarr(i,j,k) >= Real(0.0)) {
                            phi(i,j,k) = f(AMREX_D_DECL(problo[0]+i*cellsize[0],
                                                        problo[1]+j*cellsize[1],
                                                        problo[2]+k*cellsize[2]));
                        }
                    });
                }
            }
        }
    }
}

#endif

}

#endif
