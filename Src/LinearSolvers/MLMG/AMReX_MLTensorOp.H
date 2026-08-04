#ifndef AMREX_ML_TENSOR_OP_H_
#define AMREX_ML_TENSOR_OP_H_
#include <AMReX_Config.H>

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_Array.H>

namespace amrex {

/**
 * \file AMReX_MLTensorOp.H
 *
 * Cell-centered tensor operator ``alpha a v - beta div(tau)``.
 */

// Tensor solver for high Reynolds flows with small gradient in viscosity.
// The system it solves is
//
//   alpha a v - beta div dot tau = rhs
//
// where tau = eta [grad v + (grad v)^T] +  (kappa-(2/3)eta) (div v) I.
// Here eta and kappa are shear and bulk viscosity, and I is identity tensor.
//
// The user needs to provide `a` by `setACoeffs`, eta by `setShearViscosity`,
// and kappa by `setBulkViscosity`.  If `setBulkViscosity` is not called,
// kappa is set to zero.
//
// The scalars alpha and beta can be set with `setScalars(Real, Real)`.  If
// they are not set, their default value is 1.

/// \brief Tensor solver with shear/bulk viscosity inputs for high-Reynolds flows.
class MLTensorOp
    : public MLABecLaplacian
{
public:

    //! Construct an empty tensor operator; call define() before solving.
    MLTensorOp ();
    /**
     * \brief Convenience constructor without overset masks (for cell-only domains).
     *
     * Internally calls define() with the provided metadata.
     */
    MLTensorOp (const Vector<Geometry>& a_geom,
                const Vector<BoxArray>& a_grids,
                const Vector<DistributionMapping>& a_dmap,
                const LPInfo& a_info = LPInfo(),
                const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
    /**
     * \brief Convenience constructor that includes overset-mask support (1 = unknown, 0 = known).
     */
    MLTensorOp (const Vector<Geometry>& a_geom,
                const Vector<BoxArray>& a_grids,
                const Vector<DistributionMapping>& a_dmap,
                const Vector<iMultiFab const*>& a_overset_mask, // 1: unknown, 0: known
                const LPInfo& a_info = LPInfo(),
                const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
    ~MLTensorOp () override = default;

    MLTensorOp (const MLTensorOp&) = delete;
    MLTensorOp (MLTensorOp&&) = delete;
    MLTensorOp& operator= (const MLTensorOp&) = delete;
    MLTensorOp& operator= (MLTensorOp&&) = delete;

    /**
     * \brief Define the hierarchy for standard tensor solves.
     *
     * \param a_geom    Per-level geometries (fine to coarse).
     * \param a_grids   Per-level grids.
     * \param a_dmap    Distribution mappings.
     * \param a_info    Optional LPInfo overrides.
     * \param a_factory Optional FAB factories.
     */
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    /**
     * \brief Define the hierarchy with overset masks (1 = unknown, 0 = known).
     *
     * \param a_geom         Per-level geometries.
     * \param a_grids        Per-level grids.
     * \param a_dmap         Distribution mappings.
     * \param a_overset_mask Overset masks per AMR level.
     * \param a_info         Optional LPInfo overrides.
     * \param a_factory      Optional FAB factories.
     */
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const Vector<iMultiFab const*>& a_overset_mask,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    /**
     * \brief Provide shear viscosity `eta` as face-centered MultiFabs.
     *
     * \param amrlev AMR level receiving the data.
     * \param eta    Array of MultiFabs, one per spatial direction.
     */
    void setShearViscosity (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& eta);
    /**
     * \brief Populate `eta` with a constant scalar.
     *
     * \param amrlev AMR level receiving the value.
     * \param eta    Scalar shear viscosity.
     */
    void setShearViscosity (int amrlev, Real eta);
    /**
     * \brief Provide bulk viscosity `kappa` as face-centered MultiFabs.
     *
     * \param amrlev AMR level receiving the data.
     * \param kappa  Array of MultiFabs, one per spatial direction.
     */
    void setBulkViscosity (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& kappa);
    /**
     * \brief Populate `kappa` with a constant scalar.
     *
     * \param amrlev AMR level receiving the value.
     * \param kappa  Scalar bulk viscosity.
     */
    void setBulkViscosity (int amrlev, Real kappa);

    [[nodiscard]] int getNComp () const final { return AMREX_SPACEDIM; }

    [[nodiscard]] bool isCrossStencil () const final { return false; }
    [[nodiscard]] bool isTensorOp () const final { return true; }

    /**
     * \brief True if coefficients have changed since the last prepareForSolve().
     *
     * The operator may be reused for any number of solves as long as its
     * coefficients are left untouched.
     */
    [[nodiscard]] bool needsUpdate () const final {
        return (m_needs_update || MLABecLaplacian::needsUpdate());
    }
    /**
     * \brief Not implemented for tensor operators.
     *
     * Changing viscosity or `a` coefficients after prepareForSolve() is not
     * supported; build a new operator instead.
     */
    void update () final {
        amrex::Abort("MLTensorOp: update TODO");
    }

    //! Assemble viscosity data/BC caches before solving.
    void prepareForSolve () final;
    [[nodiscard]] bool isSingular (int /*armlev*/) const final { return false; }
    [[nodiscard]] bool isBottomSingular () const final { return false; }

    /**
     * \brief Apply the tensor operator on AMR level \p amrlev / MG level \p mglev.
     *
     * \param amrlev AMR level targeted by the smoothing/apply stage.
     * \param mglev  Multigrid level within that AMR level.
     * \param out    Destination MultiFab receiving ``L(in)``.
     * \param in     Source MultiFab holding the field being operated on.
     * \param bc_mode Boundary-condition mode (homogeneous vs inhomogeneous).
     * \param s_mode  State mode (solution vs correction).
     * \param bndry  Optional boundary metadata supplied by MLMG.
     */
    void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                StateMode s_mode, const MLMGBndry* bndry=nullptr) const final;

    /**
     * \brief Assemble viscous fluxes on AMR level \p amrlev.
     *
     * \param amrlev Level whose fluxes are built.
     * \param fluxes Destination arrays (one per direction).
     * \param sol    Solution values defining the flux.
     * \param loc    Sampling location for the stencil (faces, centroids, etc.).
     */
    void compFlux (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes,
                   MultiFab& sol, Location loc) const override;

    /**
     * \brief Compute velocity gradients (symmetrized) at the requested location.
     *
     * \param amrlev AMR level index.
     * \param fluxes Destination gradient MultiFabs per direction.
     * \param sol    Solution whose gradients are computed.
     * \param loc    Sampling location (cell/face center, centroid, etc.).
     */
    void compVelGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes,
                      MultiFab& sol, Location loc) const;

    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta) = delete;

protected:

    bool m_has_kappa = false;
    Vector<Vector<Array<MultiFab,AMREX_SPACEDIM> > > m_kappa;

private:

    bool m_needs_update = true;

public: // for cuda

    /**
     * \brief Apply tensor-specific boundary conditions on (\p amrlev,\p mglev).
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param vel    MultiFab updated in place with boundary values.
     * \param bc_mode Boundary-condition mode requested by MLMG.
     * \param s_mode  State mode (solution/correction).
     * \param bndry   Boundary metadata supplied by MLMG (may be null).
     */
    void applyBCTensor (int amrlev, int mglev, MultiFab& vel,
                        BCMode bc_mode, StateMode s_mode, const MLMGBndry* bndry
                        ) const;
};

}

#endif
