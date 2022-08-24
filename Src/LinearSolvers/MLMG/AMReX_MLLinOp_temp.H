#ifndef AMREX_MLLINOP_TEMP_H_
#define AMREX_MLLINOP_TEMP_H_

//! This is a template for writing your own linear operator class for Ax=b.

#include <AMReX_MLLinOp.H>

namespace amrex_temp
{

class MLLinOpTemp
    : public amrex::MLLinOp
{
public:

    //! In this example, there are 3 edge based MultiFabs.
    using Container = amrex::Array<amrex::MultiFab,3>;

    MLLinOpTemp () {}

    virtual ~MLLinOpTemp () {}

    MLLinOpTemp (const MLLinOpTemp&) = delete;
    MLLinOpTemp (MLLinOpTemp&&) = delete;
    MLLinOpTemp& operator= (const MLLinOpTemp&) = delete;
    MLLinOpTemp& operator= (MLLinOpTemp&&) = delete;

    MLLinOpTemp (const amrex::Vector<amrex::Geometry>& a_geom,
                 const amrex::Vector<amrex::BoxArray>& a_grids,
                 const amrex::Vector<amrex::DistributionMapping>& a_dmap,
                 const amrex::LPInfo& a_info = amrex::LPInfo(),
                 const amrex::Vector<amrex::FabFactory<amrex::FArrayBox> const*>& a_factory = {})
    {
        define(a_geom, a_grids, a_dmap, a_info, a_factory);
    }

    void define (const amrex::Vector<amrex::Geometry>& a_geom,
                 const amrex::Vector<amrex::BoxArray>& a_grids,
                 const amrex::Vector<amrex::DistributionMapping>& a_dmap,
                 const amrex::LPInfo& a_info = amrex::LPInfo(),
                 const amrex::Vector<amrex::FabFactory<amrex::FArrayBox> const*>& a_factory = {})
    {
        amrex::MLLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);
    }

    /**
     * \brief Return the default solver at the bottom of MG cycles.  By
     * default, MLLinOp uses a BiCGStab solver implemented in
     * AMReX::MLCGSolver.  However, it only supports a single MultiFab.
     * Since our data type is different, we use a smoother instead. In the
     * future we can try to generalize MLCGSolver.
     */
    virtual amrex::BottomSolver getDefaultBottomSolver () const override {
        return amrex::BottomSolver::smoother;
    }

    /**
     * \brief Make data container (e.g., MultiFabs stored in Any) for given level.
     *
     * \param amrlev AMR level.  Note that the lowest level is always 0.
     * \param mglev MG level. Note that mglev+1 is one level coarser than mglev.
     * \param ng number of ghost cells.
     */
    virtual amrex::Any AnyMake (int amrlev, int mglev, amrex::IntVect const& ng) const override
    {
        auto const& ba = m_grids[amrlev][mglev];
        auto const& dm = m_dmap [amrlev][mglev];
        auto const& fc = *m_factory[amrlev][mglev];
        return amrex::Any(Container{amrex::MultiFab(amrex::convert(ba,amrex::IntVect(0,1,1)),
                                                    dm, 1, ng, amrex::MFInfo(), fc),
                                    amrex::MultiFab(amrex::convert(ba,amrex::IntVect(1,0,1)),
                                                    dm, 1, ng, amrex::MFInfo(), fc),
                                    amrex::MultiFab(amrex::convert(ba,amrex::IntVect(1,1,0)),
                                                    dm, 1, ng, amrex::MFInfo(), fc)});
    }

    /**
     * \brief Make data container with coarsened BoxArray and
     * DistributionMapping of the give MG level.
     *
     * \param amrlev AMR level.  Note that the lowest level is always 0.
     * \param mglev MG level.  The coarser level is mglev+1.
     * \param ng number of ghost cells.
     */
    virtual amrex::Any AnyMakeCoarseMG (int amrlev, int mglev, amrex::IntVect const& ng) const override
    {
        auto ratio = (amrlev > 0) ? amrex::IntVect(2) : this->mg_coarsen_ratio_vec[mglev];
        auto const& ba = amrex::coarsen(m_grids[amrlev][mglev], ratio);
        auto const& dm = m_dmap[amrlev][mglev];
        return amrex::Any(Container{amrex::MultiFab(amrex::convert(ba,amrex::IntVect(0,1,1)),
                                                    dm, 1, ng),
                                    amrex::MultiFab(amrex::convert(ba,amrex::IntVect(1,0,1)),
                                                    dm, 1, ng),
                                    amrex::MultiFab(amrex::convert(ba,amrex::IntVect(1,1,0)),
                                                    dm, 1, ng)});
    }

    /**
     * \brief Make data container with coarsened BoxArray and
     * DistributionMapping of the given AMR level.
     *
     * \param famrlev AMR level.  The coarser AMR level is famrlev-1.
     * \param ng number of ghost cells.
     */
    virtual amrex::Any AnyMakeCoarseAmr (int famrlev, amrex::IntVect const& ng) const override
    {
        amrex::IntVect ratio(this->AMRRefRatio(famrlev-1));
        auto const& ba = amrex::coarsen(m_grids[famrlev][0], ratio);
        auto const& dm = m_dmap[famrlev][0];
        return amrex::Any(Container{amrex::MultiFab(amrex::convert(ba,amrex::IntVect(0,1,1)),
                                                    dm, 1, ng),
                                    amrex::MultiFab(amrex::convert(ba,amrex::IntVect(1,0,1)),
                                                    dm, 1, ng),
                                    amrex::MultiFab(amrex::convert(ba,amrex::IntVect(1,1,0)),
                                                    dm, 1, ng)});
    }

    /**
     * \brief Make an alias of the given Any without deepcopying.
     *
     * \param a an Any object.
     */
    virtual amrex::Any AnyMakeAlias (amrex::Any const& a) const override
    {
        auto const& rhs = a.get<Container>();
        return amrex::Any(Container{amrex::MultiFab(rhs[0], amrex::make_alias, 0, 1),
                                    amrex::MultiFab(rhs[1], amrex::make_alias, 0, 1),
                                    amrex::MultiFab(rhs[2], amrex::make_alias, 0, 1)});
    }

    /**
     * \brief Retuen the number of ghost cells in the given Any.
     *
     * \param a an Any object.
     */
    virtual amrex::IntVect AnyGrowVect (amrex::Any const& a) const override
    {
        auto const& mfs = a.get<Container>();
        return mfs[0].nGrowVect();
    }

    /**
     * \brief Copy data from source Any to destination Any.
     *
     * \param dst destination Any.
     * \param src source Any.
     * \param ng number of ghost cells included in the operation.
     */
    virtual void AnyCopy (amrex::Any& dst, amrex::Any const& src, amrex::IntVect const& ng) const override
    {
        auto& dmf = dst.get<Container>();
        auto const& smf = src.get<Container>();
        for (int idim=0; idim < 3; ++idim) {
            amrex::MultiFab::Copy(dmf[idim], smf[idim], 0, 0, 1, ng);
        }
    }

    /**
     * \brief Add data from source Any to destination Any.
     *
     * \param dst destination Any.
     * \param src source Any.
     * \param ng number of ghost cells included in the operation.
     */
    virtual void AnyAdd (amrex::Any& dst, amrex::Any const& src, amrex::IntVect const& ng) const override
    {
        auto& dmf = dst.get<Container>();
        auto const& smf = src.get<Container>();
        for (int idim=0; idim < 3; ++idim) {
            amrex::MultiFab::Add(dmf[idim], smf[idim], 0, 0, 1, ng);
        }
    }

    /**
     * \brief Set the given Any to zero.
     *
     * \param a an Any object.
     */
    virtual void AnySetToZero (amrex::Any& a) const override
    {
        auto& mfs = a.get<Container>();
        for (int idim=0; idim < 3; ++idim) {
            mfs[idim].setVal(amrex::Real(0.0));
        }
    }

    /**
     * \brief Set boundary (i.e., ghost cells) the given Any to zero.
     *
     * \param a an Any object.
     */
    virtual void AnySetBndryToZero (amrex::Any& a) const override
    {
        auto& mfs = a.get<Container>();
        for (int idim=0; idim < 3; ++idim) {
            mfs[idim].setBndry(amrex::Real(0.0), 0, 1);
        }
    }

#ifdef AMREX_USE_EB
    /**
     * \brief Set covered region of the given Any to zero.
     *
     * \param a an Any object.
     */
    virtual void AnySetCoveredToZero (amrex::Any& a) const override
    {
        auto& mfs = a.get<Container>();
        for (int idim=0; idim < 3; ++idim) {
            amrex::EB_set_covered(mfs[idim], 0, 1, 0, amrex::Real(0.0));
        }
    }
#endif

    /**
     * \brief ParallelCopy from source Any ot destination Any.
     *
     * \param dst destination Any.
     * \param src source Any.
     * \param src_nghost number of ghost cells in the source included in the operation.
     * \param dst_nghost number of ghost cells in the destination included in the operation.
     * \param period Periodicity.
     */
    virtual void AnyParallelCopy (amrex::Any& dst, amrex::Any const& src,
                                  amrex::IntVect const& src_nghost, amrex::IntVect const& dst_nghost,
                                  amrex::Periodicity const& period = amrex::Periodicity::NonPeriodic()) const override
    {
        auto& dmf = dst.get<Container>();
        auto const& smf = src.get<Container>();
        for (int idim=0; idim < 3; ++idim) {
            dmf[idim].ParallelCopy_nowait(smf[idim], 0, 0, 1, src_nghost, dst_nghost, period);
        }
        for (int idim=0; idim < 3; ++idim) {
            dmf[idim].ParallelCopy_finish();
        }
    }

    /**
     * \brief Return the infinity norm of the given Any.
     *
     * \param a an Any object.
     */
    virtual amrex::Real AnyNormInf (amrex::Any& a) const override
    {
        auto& mfs = a.get<Container>();
        amrex::Real r = amrex::Real(0.0);
        for (int idim=0; idim < 3; ++idim) {
            auto tmp = mfs[idim].norminf(0, 0, true);
            r = std::max(r, tmp);
        }
        amrex::ParallelAllReduce::Max(r, amrex::ParallelContext::CommunicatorSub());
        return r;
    }

    /**
     * \brief Return the infinity norm of the masked region of the given Any.
     *
     * For a composite solve with multiple AMR levels, the region covered by
     * finer AMR levels are not included in the operation.
     *
     * \parame amrlev AMR level.
     * \param a an Any object.
     * \parame local determines if the reduction is local (i.e., no MPI communication) or not.
     */
    virtual amrex::Real AnyNormInfMask (int amrlev, amrex::Any const& a, bool local) const override
    {
        amrex::ignore_unused(amrlev, a, local);
        amrex::Abort("TODO: AnyNormInfMask");
        // This is only needed for multi-level composite solve
        return amrex::Real(0.0);
    }

    /**
     * \brief Compute residual of the original form, r = b - Ax.
     *
     * \param amrlev AMR level
     * \param resid residual
     * \param x the solution x
     * \param b the RHS b
     * \param crse_bcdata provides Dirichlet BC at AMR coarse/fine interface.
     *                    It's a nullptr for single level solve.
     */
    virtual void AnySolutionResidual (int amrlev, amrex::Any& resid, amrex::Any& x, amrex::Any const& b,
                                      amrex::Any const* crse_bcdata = nullptr) override
    {
        amrex::ignore_unused(amrlev, resid, x, b, crse_bcdata);
        amrex::Abort("TODO: AnySolutionResidual");
    }

    /**
     * \brief Compute residual of the residual correction form, r = b - Ax.
     *
     * \param amrlev AMR level.
     * \param resid residual of the residual correction form.
     * \param x the correction.
     * \param b the RHS for the residual correction form (i.e., the residual of the original form.
     * \param bc_mode is either Homogeneous or Inhomogeneous.
     * \param crse_bcdata provides inhomogenous Dirichlet BC at AMR coarse/fine interface.
     *                    It's ignored for homogeneous Dirichlet BC.
     */
    virtual void AnyCorrectionResidual (int amrlev, int mglev, amrex::Any& resid, amrex::Any& x,
                                        const amrex::Any& b, MLLinOp::BCMode bc_mode,
                                        const amrex::Any* crse_bcdata=nullptr) override
    {
        amrex::ignore_unused(amrlev, mglev, resid, x, b, bc_mode, crse_bcdata);
        amrex::Abort("TODO: AnyCorrectionResidual");
    }

    /**
     * \brief Reflux
     *
     * This modifies the coarse level residual at the coarse/fine interface.
     *
     * \param crse_amrlev coarse AMR level.
     * \param res coarse level residual.
     * \param crse_sol coarse level x.
     * \param crse_rhs coarse level b.
     * \param fine_res fine level residual.  This may not be needed depending on the coarse/fine stencil.
     * \param fine_sol fine level x.
     * \param fine_rhs fine level b.
     */
    virtual void AnyReflux (int crse_amrlev,
                            amrex::Any& res, const amrex::Any& crse_sol, const amrex::Any& crse_rhs,
                            amrex::Any& fine_res, amrex::Any& fine_sol, const amrex::Any& fine_rhs) override
    {
        amrex::ignore_unused(crse_amrlev, res, crse_sol, crse_rhs, fine_res, fine_sol, fine_rhs);
        amrex::Abort("TODO: AnyReflux");
        // This is only needed for multi-level composite solve
    }

    /**
     * \brief Average down residual from fine to coarse AMR level.
     *
     * \param clev coarse ARR level.
     * \param cres coarse level residual.
     * \param fres fine level residual.
     */
    virtual void AnyAvgDownResAmr (int clev, amrex::Any& cres, amrex::Any const& fres) const override
    {
        amrex::ignore_unused(clev, cres, fres);
        amrex::Abort("TODO: AnyAvgDownResAmr");
        // This is only needed for mulit-level composite solve.
        // And maybe there is nothing neeed to be done here, like in the nodal projection solver.
    }

    /**
     * \brief Average down residual from fine to coarse MG level.
     *
     * This is only needed for MG F-cycle, and we don't need to implement this for V-cycle.
     *
     * \param clev coarse MG level.
     * \param cres coarse level residual.
     * \param fres fine level residual.
     */
    virtual void AnyAvgDownResMG (int clev, amrex::Any& cres, amrex::Any const& fres) const override
    {
        amrex::ignore_unused(clev, cres, fres);
        amrex::Abort("TODO: AnyAvgDownResMG"); // Not needed for V-cycle.
    }

    /**
     * \brief Smooth the given level.
     *
     * \param amrlev AMR level.  Note that the lowest level is always 0.
     * \param mglev MG level.  Note that mglev+1 is one level coarser than mglev.
     * \param sol x
     * \param rhs b
     * \param skip_fillboundary a flag for if we need to fill ghost cells in this function.
     */
    virtual void AnySmooth (int amrlev, int mglev, amrex::Any& sol, const amrex::Any& rhs,
                            bool skip_fillboundary=false) const override
    {
        amrex::ignore_unused(amrlev, mglev, sol, rhs, skip_fillboundary);
        amrex::Abort("TODO: AnySmooth");
    }

    /**
     * \brief Restriction from fine to coarse MG level.
     *
     * \param amrlev AMR level.
     * \param cmglev coarse MG level.  The fine MG level is cmglev-1.
     * \param crse coarse data.
     * \param fine fine data.  This is not const& because we may need to fill its ghost cells.
     */
    virtual void AnyRestriction (int amrlev, int cmglev, amrex::Any& crse, amrex::Any& fine) const override
    {
        amrex::ignore_unused(amrlev, cmglev, crse, fine);
        amrex::Abort("TODO: AnyRestriction");
    }

    /**
     * \brief Add interpolated coarse data onto the fine MG level.
     *
     * Note that it's an ADD operation.
     *
     * \param amrlev AMR level.
     * \param fmglev fine MG level.  The coarse MG level is fmglev+1.
     * \param fine fine MG level data.
     * \param crse coarse MG level data.
     */
    virtual void AnyInterpolationMG (int amrlev, int fmglev, amrex::Any& fine, const amrex::Any& crse) const override
    {
        amrex::ignore_unused(amrlev, fmglev, fine, crse);
        amrex::Abort("TODO: AnyInterpolationMG");
    }

    /**
     * \brief Assign (i.e., copy) interpolated coarse data onto the fine MG level.
     *
     * Note that it's an ASSIGN operation.  This is used in MG F-cycle, and
     * does not need to be implemented for V-cycle.
     *
     * \param amrlev AMR level.
     * \param fmglev fine MG level.  The coarse MG level is fmglev+1.
     * \param fine fine MG level data.
     * \param crse coarse MG level data.
     */
    virtual void AnyInterpAssignMG (int amrlev, int fmglev, amrex::Any& fine, amrex::Any& crse) const override
    {
        amrex::ignore_unused(amrlev, fmglev, fine, crse);
        amrex::Abort("TODO: AnyInterpAssignMG"); // not needed for V-cycle.
    }

    /**
     * \brief Interpolate data from coarse to fine AMR level.
     *
     * \param famrlev fine AMR level. The coarse AMR level is famrlev-1.
     * \param fine data on fine AMR level.
     * \param crse data on coarse AMR level.
     */
    virtual void AnyInterpolationAmr (int famrlev, amrex::Any& fine, const amrex::Any& crse,
                                      amrex::IntVect const& /*nghost*/) const override
    {
        amrex::ignore_unused(famrlev, fine, crse);
        // This is only needed for multi-level composite solve
        amrex::Abort("TODO: AnyInterpolationAmr");
    }

    /**
     * \brief Average down x and b from fine to coarse AMR level.
     *
     * This is called before V-cycle to make data on AMR levels consistent.
     *
     * \param camrlev coarse AMR level.  The fine level is camrlev+1.
     * \param crse_sol x on coarse level.
     * \param crse_rhs b on coarse level.
     * \param fine_sol x on fine level.
     * \param fine_rhs b on fine level.
     */
    virtual void AnyAverageDownSolutionRHS (int camrlev, amrex::Any& crse_sol, amrex::Any& crse_rhs,
                                            const amrex::Any& fine_sol, const amrex::Any& fine_rhs) override
    {
        amrex::ignore_unused(camrlev, crse_sol, crse_rhs, fine_sol, fine_rhs);
        // This is only needed for multi-level composite solve
        amrex::Abort("AnyAverageDownSolutionRHS");
    }

    /**
     * \brief Average down and synchronize AMR data.
     *
     * Synchronize the data on each level.  That is the nodal data in the
     * same MultiFab needs to be synchronized.  This function also needs to
     * average down the data from fine to coarse AMR levels.
     *
     * \param sol data on all AMR levels.
     */
    virtual void AnyAverageDownAndSync (amrex::Vector<amrex::Any>& sol) const override
    {
        amrex::ignore_unused(sol);
        // Even for single level, we shoudl synchronize the data on level 0.
        amrex::Abort("TODO: AnyAverageDownAndSync");
    }

    /**
     * \brief Prepare the solver for MG cycle.
     */
    virtual void prepareForSolve () override
    {
        amrex::Abort("TODO: prepareForSolve");
    }
};

}


#endif
