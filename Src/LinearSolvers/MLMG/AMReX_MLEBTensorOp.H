#ifndef AMREX_ML_EB_TENSOROP_H_
#define AMREX_ML_EB_TENSOROP_H_
#include <AMReX_Config.H>

#include <AMReX_MLEBABecLap.H>

namespace amrex {

/**
 * \file AMReX_MLEBTensorOp.H
 *
 * Embedded-boundary tensor operator for `alpha a v - beta div(tau)` where
 * viscosities may vary spatially on both cells and EB faces.
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
// kappa is set to zero.  The user must also call `setEBShearViscosity` to set
// viscosity on EB, and `setEBBulkViscosity` to set bulk viscosity on EB.
// `setBulkViscosity` and `setEBBulkViscosity` must either both be called or
// both be omitted.  Passing a value of zero still counts as calling them.
//
// The scalars alpha and beta can be set with `setScalar(Real, Real)`.  If
// they are not set, their default value is 1.

/// \brief Embedded-boundary tensor solver mirroring MLTensorOp inside the domain.
class MLEBTensorOp
    : public MLEBABecLap
{
public:

    //! Construct an empty EB tensor operator; call define() before using.
    MLEBTensorOp ();
    /**
     * \brief Convenience constructor that forwards to define().
     *
     * Requires EB factories for each AMR level.
     */
    MLEBTensorOp (const Vector<Geometry>& a_geom,
                  const Vector<BoxArray>& a_grids,
                  const Vector<DistributionMapping>& a_dmap,
                  const LPInfo& a_info,
                  const Vector<EBFArrayBoxFactory const*>& a_factory);

    ~MLEBTensorOp () override;

    MLEBTensorOp (const MLEBTensorOp&) = delete;
    MLEBTensorOp (MLEBTensorOp&&) = delete;
    MLEBTensorOp& operator= (const MLEBTensorOp&) = delete;
    MLEBTensorOp& operator= (MLEBTensorOp&&) = delete;

    /**
     * \brief Define the AMR hierarchy using EB-aware factories.
     *
     * \param a_geom    Per-level geometries (fine to coarse).
     * \param a_grids   Per-level grids.
     * \param a_dmap    Distribution mappings.
     * \param a_info    Optional LPInfo overrides.
     * \param a_factory EB factories per AMR level.
     */
    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<EBFArrayBoxFactory const*>& a_factory);

    /**
     * \brief Provide shear viscosity `eta` as face-centered MultiFabs.
     *
     * \param amrlev AMR level receiving the data.
     * \param eta    Array of MultiFabs (per spatial direction).
     * \param a_beta_loc Geometric location of the supplied coefficients.
     */
    void setShearViscosity (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& eta,
                            Location a_beta_loc);
    /**
     * \brief Populate `eta` with a scalar value on AMR level \p amrlev.
     *
     * \param amrlev AMR level receiving the constant viscosity.
     * \param eta    Scalar viscosity applied to every face.
     */
    void setShearViscosity (int amrlev, Real eta);

    /**
     * \brief Provide bulk viscosity `kappa` as face-centered MultiFabs.
     *
     * \param amrlev AMR level receiving the data.
     * \param kappa  Array of MultiFabs (per spatial direction).
     */
    void setBulkViscosity (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& kappa);
    /**
     * \brief Populate `kappa` with a scalar value.
     *
     * \param amrlev AMR level receiving the constant.
     * \param kappa  Scalar bulk viscosity.
     */
    void setBulkViscosity (int amrlev, Real kappa);

    /**
     * \brief Provide EB-face shear viscosity.
     *
     * \param amrlev AMR level receiving the data.
     * \param eta    EB-face MultiFab describing the viscosity.
     */
    void setEBShearViscosity (int amrlev, MultiFab const& eta);
    /**
     * \brief Set EB shear viscosity to a constant scalar.
     *
     * \param amrlev AMR level receiving the value.
     * \param eta    Scalar EB viscosity.
     */
    void setEBShearViscosity (int amrlev, Real eta);
    /**
     * \brief Provide EB shear viscosity plus inflow velocity corrections.
     *
     * \param amrlev AMR level receiving the data.
     * \param eta    EB-face viscosity MultiFab.
     * \param eb_vel EB-face inflow correction velocities.
     */
    void setEBShearViscosityWithInflow (int amrlev, MultiFab const& eta, MultiFab const& eb_vel);

    /**
     * \brief Provide EB bulk viscosity from a MultiFab.
     *
     * \param amrlev AMR level receiving the data.
     * \param kappa  EB-face bulk viscosity.
     */
    void setEBBulkViscosity (int amrlev, MultiFab const& kappa);
    /**
     * \brief Set EB bulk viscosity to a constant value.
     *
     * Passing zero still counts as providing the EB bulk viscosity, so
     * `setBulkViscosity` must be called as well.
     *
     * \param amrlev AMR level receiving the value.
     * \param kappa  Scalar EB bulk viscosity.
     */
    void setEBBulkViscosity (int amrlev, Real kappa);

    [[nodiscard]] int getNComp () const final { return AMREX_SPACEDIM; }

    [[nodiscard]] bool isCrossStencil () const final { return false; }
    [[nodiscard]] bool isTensorOp () const final { return true; }

    //! True if coefficients/viscosities must be refreshed before apply().
    [[nodiscard]] bool needsUpdate () const final {
        return (m_needs_update || MLEBABecLap::needsUpdate());
    }
    //! Refresh coefficient averages (not yet implemented).
    void update () final {
        amrex::Abort("MLEBTensorOp: update TODO");
    }

    //! Prepare shear/bulk coefficient caches and BC data before solving.
    void prepareForSolve () final;
    [[nodiscard]] bool isSingular (int /*armlev*/) const final { return false; }
    [[nodiscard]] bool isBottomSingular () const final { return false; }

    /**
     * \brief Apply the EB tensor operator at (\p amrlev,\p mglev).
     *
     * \param amrlev AMR level targeted by the apply.
     * \param mglev  Multigrid level within that AMR level.
     * \param out    Destination MultiFab receiving ``L(in)``.
     * \param in     Source MultiFab storing the vector field being operated on.
     * \param bc_mode Boundary-condition handling (homogeneous vs inhomogeneous).
     * \param s_mode  State interpretation (solution vs correction).
     * \param bndry   Optional boundary metadata supplied by MLMG.
     */
    void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                StateMode s_mode, const MLMGBndry* bndry=nullptr) const final;
    /**
     * \brief Assemble viscous fluxes on AMR level \p amrlev.
     *
     * \param amrlev Level whose fluxes are being computed.
     * \param fluxes Destination arrays (one per spatial direction).
     * \param sol    Solution whose fluxes are being computed.
     * \param loc    Geometric location for the stencil (face center, centroid, etc.).
     */
    void compFlux (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& fluxes,
                   MultiFab& sol, Location loc) const override;

    /**
     * \brief Compute velocity gradients at the requested location.
     *
     * \param amrlev AMR level index.
     * \param grads  Destination gradient MultiFabs per direction.
     * \param sol    Solution whose gradients are computed.
     * \param loc    Location (cell/face center, centroid, etc.).
     */
    void compVelGrad (int amrlev, const Array<MultiFab*,AMREX_SPACEDIM>& grads,
                      MultiFab& sol, Location loc) const;

    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta,
                     Location a_beta_loc) = delete;
    void setEBDirichlet (int amrlev, const MultiFab& phi, const MultiFab& beta) = delete;
    void setEBHomogDirichlet (int amrlev, const MultiFab& beta) = delete;

protected:

    bool m_has_kappa = false;
    bool m_has_eb_kappa = false;
    Vector<Vector<Array<MultiFab,AMREX_SPACEDIM> > > m_kappa;
    Vector<Vector<MultiFab> > m_eb_kappa;
    mutable Vector<Vector<Array<MultiFab,AMREX_SPACEDIM> > > m_tauflux;

private:

    bool m_needs_update = true;

public: // for cuda

    /**
     * \brief Apply tensor boundary conditions prior to an operator call.
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param vel    Velocity field updated in place.
     * \param bc_mode Boundary-condition mode supplied by MLMG.
     * \param s_mode  State interpretation (solution/correction).
     * \param bndry   Boundary metadata (may be null to rebuild on the fly).
     */
    void applyBCTensor (int amrlev, int mglev, MultiFab& vel,
                        BCMode bc_mode, StateMode s_mode, const MLMGBndry* bndry) const;
    /**
     * \brief Compute cross-derivative stress terms for cached BC data.
     *
     * \param amrlev AMR level index.
     * \param mglev  Multigrid level index.
     * \param mf     MultiFab providing the field whose derivatives are needed.
     * \param bndry  Boundary metadata describing where to store the terms.
     */
    void compCrossTerms(int amrlev, int mglev, MultiFab const& mf,
                        const MLMGBndry* bndry) const;
};

}

#endif

// add ghost cells to kappa
