#ifndef AMREX_ML_LINOP_H_
#define AMREX_ML_LINOP_H_
#include <AMReX_Config.H>

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
#include <AMReX_Hypre.H>
#include <AMReX_HypreNodeLap.H>
#endif

#if defined(AMREX_USE_PETSC) && (AMREX_SPACEDIM > 1)
#include <AMReX_PETSc.H>
#endif

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiCutFab.H>
#endif

#include <AMReX_Any.H>
#include <AMReX_BndryRegister.H>
#include <AMReX_FabDataType.H>
#include <AMReX_MLMGBndry.H>
#include <AMReX_MultiFabUtil.H>

#include <algorithm>
#include <string>

namespace amrex {

enum class BottomSolver : int {
    Default, smoother, bicgstab, cg, bicgcg, cgbicg, hypre, petsc
};

struct LPInfo
{
    bool do_agglomeration = true;
    bool do_consolidation = true;
    bool do_semicoarsening = false;
    int agg_grid_size = -1;
    int con_grid_size = -1;
    int con_ratio = 2;
    int con_strategy = 3;
    bool has_metric_term = true;
    int max_coarsening_level = 30;
    int max_semicoarsening_level = 0;
    int semicoarsening_direction = -1;
    int hidden_direction = -1;

    LPInfo& setAgglomeration (bool x) noexcept { do_agglomeration = x; return *this; }
    LPInfo& setConsolidation (bool x) noexcept { do_consolidation = x; return *this; }
    LPInfo& setSemicoarsening (bool x) noexcept { do_semicoarsening = x; return *this; }
    LPInfo& setAgglomerationGridSize (int x) noexcept { agg_grid_size = x; return *this; }
    LPInfo& setConsolidationGridSize (int x) noexcept { con_grid_size = x; return *this; }
    LPInfo& setConsolidationRatio (int x) noexcept { con_ratio = x; return *this; }
    LPInfo& setConsolidationStrategy (int x) noexcept { con_strategy = x; return *this; }
    LPInfo& setMetricTerm (bool x) noexcept { has_metric_term = x; return *this; }
    LPInfo& setMaxCoarseningLevel (int n) noexcept { max_coarsening_level = n; return *this; }
    LPInfo& setMaxSemicoarseningLevel (int n) noexcept { max_semicoarsening_level = n; return *this; }
    LPInfo& setSemicoarseningDirection (int n) noexcept { semicoarsening_direction = n; return *this; }
    LPInfo& setHiddenDirection (int n) noexcept { hidden_direction = n; return *this; }

    [[nodiscard]] bool hasHiddenDimension () const noexcept {
        return hidden_direction >=0 && hidden_direction < AMREX_SPACEDIM;
    }

    static constexpr int getDefaultAgglomerationGridSize () {
#ifdef AMREX_USE_GPU
        return 32;
#else
        return AMREX_D_PICK(32, 16, 8);
#endif
    }

    static constexpr int getDefaultConsolidationGridSize () {
#ifdef AMREX_USE_GPU
        return 32;
#else
        return AMREX_D_PICK(32, 16, 8);
#endif
    }
};

struct LinOpEnumType
{
    enum struct BCMode { Homogeneous, Inhomogeneous };
    enum struct StateMode { Solution, Correction };
    enum struct Location { FaceCenter, FaceCentroid, CellCenter, CellCentroid };
};

template <typename T> class MLMGT;
template <typename T> class MLCGSolverT;
template <typename T> class MLPoissonT;
template <typename T> class MLABecLaplacianT;
template <typename T> class GMRESMLMGT;

template <typename MF>
class MLLinOpT
{
public:

    template <typename T> friend class MLMGT;
    template <typename T> friend class MLCGSolverT;
    template <typename T> friend class MLPoissonT;
    template <typename T> friend class MLABecLaplacianT;
    template <typename T> friend class GMRESMLMGT;

    using MFType = MF;
    using FAB = typename FabDataType<MF>::fab_type;
    using RT  = typename FabDataType<MF>::value_type;

    using BCType = LinOpBCType;
    using BCMode    = LinOpEnumType::BCMode;
    using StateMode = LinOpEnumType::StateMode;
    using Location  = LinOpEnumType::Location;

    MLLinOpT () = default;
    virtual ~MLLinOpT () = default;

    MLLinOpT (const MLLinOpT<MF>&) = delete;
    MLLinOpT (MLLinOpT<MF>&&) = delete;
    MLLinOpT<MF>& operator= (const MLLinOpT<MF>&) = delete;
    MLLinOpT<MF>& operator= (MLLinOpT<MF>&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<FabFactory<FAB> const*>& a_factory,
                 bool eb_limit_coarsening = true);

    [[nodiscard]] virtual std::string name () const { return std::string("Unspecified"); }

    /**
     * \brief Boundary of the whole domain.
     *
     * This functions must be called, and must be called before other bc
     * functions. This version is for single-component solve or when all the
     * components have the same BC types.
     *
     * \param lobc lower boundaries
     * \param hibc upper boundaries
     */
    void setDomainBC (const Array<BCType,AMREX_SPACEDIM>& lobc,
                      const Array<BCType,AMREX_SPACEDIM>& hibc) noexcept;

    /**
     * \brief Boundary of the whole domain.
     *
     * This functions must be called, and must be called before other bc
     * functions. This version is for multi-component solve.
     *
     * \param lobc lower boundaries
     * \param hibc upper boundaries
     */
    void setDomainBC (const Vector<Array<BCType,AMREX_SPACEDIM> >& lobc,
                      const Vector<Array<BCType,AMREX_SPACEDIM> >& hibc) noexcept;

    /**
     * \brief Set location of domain boundaries.
     *
     * By default, domain BC is on the domain face.  If that's the
     * case, this function doesn't need to be called.  However, one
     * could use this function to set non-zero domain BC locations.
     * Note all values should be >= 0.  If this function is called,
     * it MUST be called before setLevelBC.
     */
    void setDomainBCLoc (const Array<Real,AMREX_SPACEDIM>& lo_bcloc,
                         const Array<Real,AMREX_SPACEDIM>& hi_bcloc) noexcept;

    /**
    * \brief Needs coarse data for bc?
    *
    * If the lowest level grids does not cover the entire domain, coarse
    * level data are needed for supplying Dirichlet bc at coarse/fine
    * boundary, even when the domain bc is not Dirichlet.
    */
    [[nodiscard]] bool needsCoarseDataForBC () const noexcept { return m_needs_coarse_data_for_bc; }

    /**
    * \brief Set coarse/fine boundary conditions. For cell-centered solves
    * only.
    *
    * If we want to do a linear solve where the boundary conditions on the
    * coarsest AMR level of the solve come from a coarser level (e.g. the
    * base AMR level of the solve is > 0 and does not cover the entire
    * domain), we must explicitly provide the coarser data.  Boundary
    * conditions from a coarser level are Dirichlet by default.  The MultiFab
    * crse does not need to have ghost cells and is at a coarser resolution
    * than the coarsest AMR level of the solve; it is used to supply
    * (interpolated) boundary conditions for the solve.  NOTE: If this is
    * called, it must be called before `setLevelBC`.  If crse is nullptr,
    * then the bc values are assumed to be zero. The coarse/fine BC type
    * can be changed to homogeneous Neumann by the bc_type argument. In that
    * case, use nullptr for the crse argument.
    *
    * \param crse       the coarse AMR level data
    * \param crse_ratio the coarsening ratio between fine and coarse AMR levels.
    * \param bc_type    optional. It's Dirichlet by default, and can be Neumann.
    */
    void setCoarseFineBC (const MF* crse, int crse_ratio,
                          LinOpBCType bc_type = LinOpBCType::Dirichlet) noexcept;

    void setCoarseFineBC (const MF* crse, IntVect const& crse_ratio,
                          LinOpBCType bc_type = LinOpBCType::Dirichlet) noexcept;

    template <typename AMF, std::enable_if_t<!std::is_same_v<MF,AMF>,int> = 0>
    void setCoarseFineBC (const AMF* crse, int crse_ratio,
                          LinOpBCType bc_type = LinOpBCType::Dirichlet) noexcept;

    template <typename AMF, std::enable_if_t<!std::is_same_v<MF,AMF>,int> = 0>
    void setCoarseFineBC (const AMF* crse, IntVect const& crse_ratio,
                          LinOpBCType bc_type = LinOpBCType::Dirichlet) noexcept;

    /**
    * \brief Set boundary conditions for given level. For cell-centered
    * solves only.
    *
    * This must be called for each level.  Argument `levelbcdata` is used to
    * supply Dirichlet or Neumann bc at the physical domain; if those data
    * are homogeneous we can pass nullptr instead of levelbcdata.
    * Regardless, this function must be called.  If used, the MultiFab
    * levelbcdata must have one ghost cell.  Only the data outside the
    * physical domain will be used.  It is assumed that the data in those
    * ghost cells outside the domain live exactly on the face of the
    * physical domain.  Argument `amrlev` is relative level such that the
    * lowest to the solver is always 0.  The optional arguments
    * robinbc_[a|b|f] provide Robin boundary condition `a*phi + b*dphi/dn =
    * f`.  Note that `d./dn` is `d./dx` at the upper boundary and `-d./dx`
    * at the lower boundary, for Robin BC. However, for inhomogeneous
    * Neumann BC, the value in leveldata is assumed to be `d./dx`.
    */
    virtual void setLevelBC (int /*amrlev*/, const MF* /*levelbcdata*/,
                             const MF* /*robinbc_a*/ = nullptr,
                             const MF* /*robinbc_b*/ = nullptr,
                             const MF* /*robinbc_f*/ = nullptr) = 0;

    template <typename AMF, std::enable_if_t<!std::is_same_v<MF,AMF>,int> = 0>
    void setLevelBC (int amrlev, const AMF* levelbcdata,
                     const AMF* robinbc_a = nullptr,
                     const AMF* robinbc_b = nullptr,
                     const AMF* robinbc_f = nullptr);

    //! Set verbosity
    void setVerbose (int v) noexcept { verbose = v; }

    //! Set order of interpolation at coarse/fine boundary
    void setMaxOrder (int o) noexcept { maxorder = o; }
    //! Get order of interpolation at coarse/fine boundary
    [[nodiscard]] int getMaxOrder () const noexcept { return maxorder; }

    //! Set the flag for whether the solver should try to make singular
    //! problem solvable, which is on by default.
    void setEnforceSingularSolvable (bool o) noexcept { enforceSingularSolvable = o; }
    //! Get the flag for whether the solver should try to make singular
    //! problem solvable.
    [[nodiscard]] bool getEnforceSingularSolvable () const noexcept { return enforceSingularSolvable; }

    [[nodiscard]] virtual BottomSolver getDefaultBottomSolver () const { return BottomSolver::bicgstab; }

    //! Return number of components
    [[nodiscard]] virtual int getNComp () const { return 1; }

    [[nodiscard]] virtual int getNGrow (int /*a_lev*/ = 0, int /*mg_lev*/ = 0) const { return 0; }

    //! Does it need update if it's reused?
    [[nodiscard]] virtual bool needsUpdate () const { return false; }
    //! Update for reuse.
    virtual void update () {}

    /**
     * \brief Restriction onto coarse MG level
     *
     * \param amrlev AMR level
     * \param cmglev coarse MG level
     * \param crse   coarse data. This is the output.
     * \param fine   fine data. This is the input. Some operators might need to fill ghost cells.
     *               This is why it's not a const reference.
     */
    virtual void restriction (int amrlev, int cmglev, MF& crse, MF& fine) const = 0;

    /**
     * \brief Add interpolated coarse MG level data to fine MG level data
     *
     * \param amrlev AMR level
     * \param fmglev fine MG level
     * \param crse   coarse data.
     * \param fine   fine data.
     */
    virtual void interpolation (int amrlev, int fmglev, MF& fine, const MF& crse) const = 0;

    /**
     * \brief Overwrite fine MG level data with interpolated coarse data.
     *
     * \param amrlev AMR level
     * \param fmglev fine MG level
     * \param fine   fine MG level data
     * \param crse   coarse MG level data
     */
    virtual void interpAssign (int amrlev, int fmglev, MF& fine, MF& crse) const
    {
        amrex::ignore_unused(amrlev, fmglev, fine, crse);
        amrex::Abort("MLLinOpT::interpAssign: Must be implemented for FMG cycle");
    }

    /**
     * \brief Interpolation between AMR levels
     *
     * \param famrlev fine AMR level
     * \param fine    fine level data
     * \param crse    coarse level data
     * \param nghost number of ghost cells
     */
    virtual void interpolationAmr (int famrlev, MF& fine, const MF& crse,
                                   IntVect const& nghost) const
    {
        amrex::ignore_unused(famrlev, fine, crse, nghost);
        amrex::Abort("MLLinOpT::interpolationAmr: Must be implemented for composite solves across multiple AMR levels");
    }

    /**
     * \brief Average-down data from fine AMR level to coarse AMR level.
     *
     * \param camrlev  coarse AMR level
     * \param crse_sol solutoin on coarse AMR level
     * \param crse_rhs RHS on coarse AMR level
     * \param fine_sol solution on fine AMR level
     * \param fine_rhs RHS on fine AMR level
     */
    virtual void averageDownSolutionRHS (int camrlev, MF& crse_sol, MF& crse_rhs,
                                         const MF& fine_sol, const MF& fine_rhs)
    {
        amrex::ignore_unused(camrlev, crse_sol, crse_rhs, fine_sol, fine_rhs);
        amrex::Abort("MLLinOpT::averageDownSolutionRHS: Must be implemented for composite solves across multiple AMR levels");
    }

    /**
     * \brief Apply the linear operator, out = L(in)
     *
     * \param amrlev  AMR level
     * \param mglev   MG level
     * \param out     output
     * \param in      input
     * \param bc_mode Is the BC homogeneous or inhomogeneous?
     * \param s_mode  Are data data solution or correction?
     * \param bndry   object for handling coarse/fine and physical boundaries
     */
    virtual void apply (int amrlev, int mglev, MF& out, MF& in, BCMode bc_mode,
                        StateMode s_mode, const MLMGBndryT<MF>* bndry=nullptr) const = 0;

    /**
     * \brief Smooth
     *
     * \param amrlev            AMR level
     * \param mglev             MG level
     * \param sol               unknowns
     * \param rhs               RHS
     * \param skip_fillboundary flag controlling whether ghost cell filling can be skipped.
     */
    virtual void smooth (int amrlev, int mglev, MF& sol, const MF& rhs,
                         bool skip_fillboundary=false) const = 0;

    //! Divide mf by the diagonal component of the operator. Used by bicgstab.
    virtual void normalize (int amrlev, int mglev, MF& mf) const {
        amrex::ignore_unused(amrlev, mglev, mf);
    }

    /**
     * \brief Compute residual for solution
     *
     * \param amrlev       AMR level
     * \param resid        residual
     * \param x            solution
     * \param b            RHS
     * \param crse_bc_data optional argument providing BC at coarse/fine boundary.
     */
    virtual void solutionResidual (int amrlev, MF& resid, MF& x, const MF& b,
                                   const MF* crse_bcdata=nullptr) = 0;

    virtual void prepareForFluxes (int /*amrlev*/, const MF* /*crse_bcdata*/ = nullptr) {}

    /**
     * \brief Compute residual for the residual-correction form, resid = b - L(x)
     *
     * \param amrlev       AMR level
     * \param mglev        MG level
     * \param resid        residual
     * \param x            unknown in the residual-correction form
     * \param b            RHS in the residual-correction form
     * \param bc_mode      Is the BC homogeneous or inhomogeneous?
     * \param crse_bc_data optional argument providing BC at coarse/fine boundary.
     */
    virtual void correctionResidual (int amrlev, int mglev, MF& resid, MF& x, const MF& b,
                                     BCMode bc_mode, const MF* crse_bcdata=nullptr) = 0;

    /**
     * \brief Reflux at AMR coarse/fine boundary
     *
     * \param crse_amrlev coarse AMR level
     * \param res         coarse level residual
     * \param crse_sol    coarse level solution
     * \param crse_rhs    coarse level RHS
     * \param fine_res    fine level residual
     * \param fine_sol    fine level solution
     * \param fine_rhs    fine level RHS
     */
    virtual void reflux (int crse_amrlev,
                         MF& res, const MF& crse_sol, const MF& crse_rhs,
                         MF& fine_res, MF& fine_sol, const MF& fine_rhs) const
    {
        amrex::ignore_unused(crse_amrlev, res, crse_sol, crse_rhs, fine_res,
                             fine_sol, fine_rhs);
        amrex::Abort("MLLinOpT::reflux: Must be implemented for composite solves across multiple AMR levels");
    }

    /**
     * \brief Compute fluxes
     *
     * \param amrlev AMR level
     * \param fluxes fluxes
     * \param sol    solution
     * \param loc    location of the fluxes
     */
    virtual void compFlux (int /*amrlev*/, const Array<MF*,AMREX_SPACEDIM>& /*fluxes*/,
                           MF& /*sol*/, Location /*loc*/) const
    {
        amrex::Abort("AMReX_MLLinOp::compFlux::How did we get here?");
    }

    /**
     * \brief Compute gradients of the solution
     *
     * \param amrlev AMR level
     * \param grad   grad(sol)
     * \param sol    solution
     * \param loc    location of the gradients
     */
    virtual void compGrad (int /*amrlev*/, const Array<MF*,AMREX_SPACEDIM>& /*grad*/,
                           MF& /*sol*/, Location /*loc*/) const
    {
        amrex::Abort("AMReX_MLLinOp::compGrad::How did we get here?");
    }

    //! apply metric terms if there are any
    virtual void applyMetricTerm (int /*amrlev*/, int /*mglev*/, MF& /*rhs*/) const {}
    //! unapply metric terms if there are any
    virtual void unapplyMetricTerm (int /*amrlev*/, int /*mglev*/, MF& /*rhs*/) const {}

    //! This is needed for our nodal projection solver
    virtual void unimposeNeumannBC (int /*amrlev*/, MF& /*rhs*/) const {}

    //! Extra terms introduced when we treat inhomogeneous Nuemann BC as homogeneous.
    virtual void applyInhomogNeumannTerm (int /*amrlev*/, MF& /*rhs*/) const {}

    //! for overset solver only
    virtual void applyOverset (int /*amlev*/, MF& /*rhs*/) const {}

    //! scale RHS to fix solvability
    [[nodiscard]] virtual bool scaleRHS (int /*amrlev*/, MF* /*rhs*/) const {
        return false;
    }

    //! get offset for fixing solvability
    virtual Vector<RT> getSolvabilityOffset (int /*amrlev*/, int /*mglev*/,
                                               MF const& /*rhs*/) const { return {}; }

    //! fix solvability by subtracting offset from RHS
    virtual void fixSolvabilityByOffset (int /*amrlev*/, int /*mglev*/, MF& /*rhs*/,
                                         Vector<RT> const& /*offset*/) const {}

    virtual void prepareForSolve () = 0;

    virtual void preparePrecond () {}

    //! For GMRES to work, this might need to be implemented to mask out
    //! Dirichlet nodes or cells (e.g., EB covered cells or overset cells)
    virtual void setDirichletNodesToZero (int /*amrlev*/, int /*mglev*/,
                                          MF& /*mf*/) const
    {
        amrex::Warning("This function might need to be implemented for GMRES to work with this LinOp.");
    }

    //! Is it singular on given AMR level?
    [[nodiscard]] virtual bool isSingular (int amrlev) const = 0;
    //! Is the bottom of MG singular?
    [[nodiscard]] virtual bool isBottomSingular () const = 0;

    //! x dot y, used by the bottom solver
    virtual RT xdoty (int amrlev, int mglev, const MF& x, const MF& y, bool local) const = 0;

    virtual RT dotProductPrecond (Vector<MF const*> const& x, Vector<MF const*> const& y) const;

    virtual RT norm2Precond (Vector<MF const*> const& x) const;

    virtual std::unique_ptr<MLLinOpT<MF>> makeNLinOp (int /*grid_size*/) const
    {
        amrex::Abort("MLLinOp::makeNLinOp: N-Solve not supported");
        return nullptr;
    }

    virtual void getFluxes (const Vector<Array<MF*,AMREX_SPACEDIM> >& /*a_flux*/,
                            const Vector<MF*>& /*a_sol*/,
                            Location /*a_loc*/) const {
        amrex::Abort("MLLinOp::getFluxes: How did we get here?");
    }
    virtual void getFluxes (const Vector<MF*>& /*a_flux*/,
                            const Vector<MF*>& /*a_sol*/) const {
        amrex::Abort("MLLinOp::getFluxes: How did we get here?");
    }

#ifdef AMREX_USE_EB
    virtual void getEBFluxes (const Vector<MF*>& /*a_flux*/,
                              const Vector<MF*>& /*a_sol*/) const {
        amrex::Abort("MLLinOp::getEBFluxes: How did we get here?");
    }
#endif

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    [[nodiscard]] virtual std::unique_ptr<Hypre> makeHypre (Hypre::Interface /*hypre_interface*/) const {
        amrex::Abort("MLLinOp::makeHypre: How did we get here?");
        return {nullptr};
    }
    [[nodiscard]] virtual std::unique_ptr<HypreNodeLap> makeHypreNodeLap(
        int /*bottom_verbose*/,
        const std::string& /* options_namespace */) const
    {
        amrex::Abort("MLLinOp::makeHypreNodeLap: How did we get here?");
        return {nullptr};
    }
#endif

#if defined(AMREX_USE_PETSC) && (AMREX_SPACEDIM > 1)
    [[nodiscard]] virtual std::unique_ptr<PETScABecLap> makePETSc () const {
        amrex::Abort("MLLinOp::makePETSc: How did we get here?");
        return {nullptr};
    }
#endif

    [[nodiscard]] virtual bool supportNSolve () const { return false; }

    virtual void copyNSolveSolution (MF&, MF const&) const {}

    virtual void postSolve (Vector<MF*> const& /*sol*/) const {}

    [[nodiscard]] virtual RT normInf (int amrlev, MF const& mf, bool local) const = 0;

    virtual void averageDownAndSync (Vector<MF>& sol) const = 0;

    virtual void avgDownResAmr (int clev, MF& cres, MF const& fres) const
    {
        amrex::ignore_unused(clev, cres, fres);
        amrex::Abort("MLLinOpT::avgDownResAmr: Must be implemented for composite solves across multiple AMR levels");
    }

    // This function is needed for FMG cycle, but not V-cycle.
    virtual void avgDownResMG (int clev, MF& cres, MF const& fres) const;

    virtual void beginPrecondBC () { m_precond_mode = true; }
    virtual void endPrecondBC () { m_precond_mode = false; }

    [[nodiscard]] bool isMFIterSafe (int amrlev, int mglev1, int mglev2) const;

    //! Return the number of AMR levels
    [[nodiscard]] int NAMRLevels () const noexcept { return m_num_amr_levels; }

    //! Return the number of MG levels at given AMR level
    [[nodiscard]] int NMGLevels (int amrlev) const noexcept { return m_num_mg_levels[amrlev]; }

    [[nodiscard]] const Geometry& Geom (int amr_lev, int mglev=0) const noexcept { return m_geom[amr_lev][mglev]; }

    // BC
    Vector<Array<BCType, AMREX_SPACEDIM> > m_lobc;
    Vector<Array<BCType, AMREX_SPACEDIM> > m_hibc;
    // Need to save the original copy because we change the BC type to
    // Neumann for inhomogeneous Neumann and Robin.
    Vector<Array<BCType, AMREX_SPACEDIM> > m_lobc_orig;
    Vector<Array<BCType, AMREX_SPACEDIM> > m_hibc_orig;

protected:

    static constexpr int mg_coarsen_ratio = 2;
    static constexpr int mg_box_min_width = 2;
    int mg_domain_min_width = 2;

    LPInfo info;

    int verbose = 0;

    int maxorder = 3;

    bool enforceSingularSolvable = true;

    int m_num_amr_levels = 0;
    Vector<int> m_amr_ref_ratio;

    Vector<int> m_num_mg_levels;
    const MLLinOpT<MF>* m_parent = nullptr;

    IntVect m_ixtype;

    bool m_do_agglomeration = false;
    bool m_do_consolidation = false;

    bool m_do_semicoarsening = false;
    Vector<IntVect> mg_coarsen_ratio_vec;

    //! first Vector is for amr level and second is mg level
    Vector<Vector<Geometry> >            m_geom;
    Vector<Vector<BoxArray> >            m_grids;
    Vector<Vector<DistributionMapping> > m_dmap;
    Vector<Vector<std::unique_ptr<FabFactory<FAB> > > > m_factory;
    Vector<int>                          m_domain_covered;

    MPI_Comm m_default_comm = MPI_COMM_NULL;
    MPI_Comm m_bottom_comm = MPI_COMM_NULL;
    struct CommContainer {
        MPI_Comm comm;
        CommContainer (MPI_Comm m) noexcept : comm(m) {}
        CommContainer (const CommContainer&) = delete;
        CommContainer (CommContainer&&) = delete;
        void operator= (const CommContainer&) = delete;
        void operator= (CommContainer&&) = delete;
        ~CommContainer () { // NOLINT(modernize-use-equals-default)
#ifdef BL_USE_MPI
            if (comm != MPI_COMM_NULL) { MPI_Comm_free(&comm); }
#endif
        }
    };
    std::unique_ptr<CommContainer> m_raii_comm;

    Array<Real, AMREX_SPACEDIM> m_domain_bloc_lo {{AMREX_D_DECL(0._rt,0._rt,0._rt)}};
    Array<Real, AMREX_SPACEDIM> m_domain_bloc_hi {{AMREX_D_DECL(0._rt,0._rt,0._rt)}};

    bool m_needs_coarse_data_for_bc = false;
    LinOpBCType m_coarse_fine_bc_type = LinOpBCType::Dirichlet;
    IntVect m_coarse_data_crse_ratio = IntVect(-1);
    RealVect m_coarse_bc_loc;
    const MF* m_coarse_data_for_bc = nullptr;
    MF m_coarse_data_for_bc_raii;

    bool m_precond_mode = false;

    //! Return AMR refinement ratios
    [[nodiscard]] const Vector<int>& AMRRefRatio () const noexcept { return m_amr_ref_ratio; }

    //! Return AMR refinement ratio at given AMR level
    [[nodiscard]] int AMRRefRatio (int amr_lev) const noexcept { return m_amr_ref_ratio[amr_lev]; }

    [[nodiscard]] FabFactory<FAB> const* Factory (int amr_lev, int mglev=0) const noexcept {
        return m_factory[amr_lev][mglev].get();
    }

    [[nodiscard]] GpuArray<BCType,AMREX_SPACEDIM> LoBC (int icomp = 0) const noexcept {
        return GpuArray<BCType,AMREX_SPACEDIM>{{AMREX_D_DECL(m_lobc[icomp][0],
                                                             m_lobc[icomp][1],
                                                             m_lobc[icomp][2])}};
    }
    [[nodiscard]] GpuArray<BCType,AMREX_SPACEDIM> HiBC (int icomp = 0) const noexcept {
        return GpuArray<BCType,AMREX_SPACEDIM>{{AMREX_D_DECL(m_hibc[icomp][0],
                                                             m_hibc[icomp][1],
                                                             m_hibc[icomp][2])}};
    }

    [[nodiscard]] bool hasBC (BCType bct) const noexcept;
    [[nodiscard]] bool hasInhomogNeumannBC () const noexcept;
    [[nodiscard]] bool hasRobinBC () const noexcept;

    [[nodiscard]] virtual bool supportRobinBC () const noexcept { return false; }
    [[nodiscard]] virtual bool supportInhomogNeumannBC () const noexcept { return false; }

#ifdef BL_USE_MPI
    [[nodiscard]] bool isBottomActive () const noexcept { return m_bottom_comm != MPI_COMM_NULL; }
#else
    [[nodiscard]] bool isBottomActive () const noexcept { return true; }
#endif
    [[nodiscard]] MPI_Comm BottomCommunicator () const noexcept { return m_bottom_comm; }
    [[nodiscard]] MPI_Comm Communicator () const noexcept { return m_default_comm; }

    void setCoarseFineBCLocation (const RealVect& cloc) noexcept { m_coarse_bc_loc = cloc; }

    [[nodiscard]] bool doAgglomeration () const noexcept { return m_do_agglomeration; }
    [[nodiscard]] bool doConsolidation () const noexcept { return m_do_consolidation; }
    [[nodiscard]] bool doSemicoarsening () const noexcept { return m_do_semicoarsening; }

    [[nodiscard]] bool isCellCentered () const noexcept { return m_ixtype == 0; }

    [[nodiscard]] virtual IntVect getNGrowVectRestriction () const {
        return isCellCentered() ? IntVect(0) : IntVect(1);
    }

    virtual void make (Vector<Vector<MF> >& mf, IntVect const& ng) const;

    [[nodiscard]] virtual MF make (int amrlev, int mglev, IntVect const& ng) const;

    [[nodiscard]] virtual MF makeAlias (MF const& mf) const;

    [[nodiscard]] virtual MF makeCoarseMG (int amrlev, int mglev, IntVect const& ng) const;

    [[nodiscard]] virtual MF makeCoarseAmr (int famrlev, IntVect const& ng) const;

    [[nodiscard]] virtual std::unique_ptr<FabFactory<FAB> > makeFactory (int /*amrlev*/, int /*mglev*/) const {
        return std::make_unique<DefaultFabFactory<FAB>>();
    }

    virtual void resizeMultiGrid (int new_size);

    [[nodiscard]] bool hasHiddenDimension () const noexcept { return info.hasHiddenDimension(); }
    [[nodiscard]] int hiddenDirection () const noexcept { return info.hidden_direction; }
    [[nodiscard]] Box compactify (Box const& b) const noexcept;

    template <typename T>
    [[nodiscard]] Array4<T> compactify (Array4<T> const& a) const noexcept
    {
        if (info.hidden_direction == 0) {
            return Array4<T>(a.dataPtr(), {a.begin.y,a.begin.z,0}, {a.end.y,a.end.z,1}, a.nComp());
        } else if (info.hidden_direction == 1) {
            return Array4<T>(a.dataPtr(), {a.begin.x,a.begin.z,0}, {a.end.x,a.end.z,1}, a.nComp());
        } else if (info.hidden_direction == 2) {
            return Array4<T>(a.dataPtr(), {a.begin.x,a.begin.y,0}, {a.end.x,a.end.y,1}, a.nComp());
        } else {
            return a;
        }
    }

    template <typename T>
    [[nodiscard]] T get_d0 (T const& dx, T const& dy, T const&) const noexcept
    {
        if (info.hidden_direction == 0) {
            return dy;
        } else {
            return dx;
        }
    }

    template <typename T>
    [[nodiscard]] T get_d1 (T const&, T const& dy, T const& dz) const noexcept
    {
        if (info.hidden_direction == 0 || info.hidden_direction == 1) {
            return dz;
        } else {
            return dy;
        }
    }

private:

    void defineGrids (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap,
                      const Vector<FabFactory<FAB> const*>& a_factory);
    void defineBC ();
    static void makeAgglomeratedDMap (const Vector<BoxArray>& ba, Vector<DistributionMapping>& dm);
    static void makeConsolidatedDMap (const Vector<BoxArray>& ba, Vector<DistributionMapping>& dm,
                                      int ratio, int strategy);
    [[nodiscard]] MPI_Comm makeSubCommunicator (const DistributionMapping& dm);

    virtual void checkPoint (std::string const& /*file_name*/) const {
        amrex::Abort("MLLinOp:checkPoint: not implemented");
    }

    Vector<std::unique_ptr<MF>> levelbc_raii;
    Vector<std::unique_ptr<MF>> robin_a_raii;
    Vector<std::unique_ptr<MF>> robin_b_raii;
    Vector<std::unique_ptr<MF>> robin_f_raii;
};

template <typename MF>
void
MLLinOpT<MF>::define (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap,
                      const LPInfo& a_info,
                      const Vector<FabFactory<FAB> const*>& a_factory,
                      [[maybe_unused]] bool eb_limit_coarsening)
{
    BL_PROFILE("MLLinOp::define()");

    info = a_info;
#ifdef AMREX_USE_GPU
    if (Gpu::notInLaunchRegion())
    {
        if (info.agg_grid_size <= 0) { info.agg_grid_size = AMREX_D_PICK(32, 16, 8); }
        if (info.con_grid_size <= 0) { info.con_grid_size = AMREX_D_PICK(32, 16, 8); }
    }
    else
#endif
    {
        if (info.agg_grid_size <= 0) { info.agg_grid_size = LPInfo::getDefaultAgglomerationGridSize(); }
        if (info.con_grid_size <= 0) { info.con_grid_size = LPInfo::getDefaultConsolidationGridSize(); }
    }

#ifdef AMREX_USE_EB
    if (!a_factory.empty() && eb_limit_coarsening) {
        const auto *f = dynamic_cast<EBFArrayBoxFactory const*>(a_factory[0]);
        if (f) {
            info.max_coarsening_level = std::min(info.max_coarsening_level,
                                                 f->maxCoarseningLevel());
        }
    }
#endif
    defineGrids(a_geom, a_grids, a_dmap, a_factory);
    defineBC();
}

template <typename MF>
void
MLLinOpT<MF>::defineGrids (const Vector<Geometry>& a_geom,
                           const Vector<BoxArray>& a_grids,
                           const Vector<DistributionMapping>& a_dmap,
                           const Vector<FabFactory<FAB> const*>& a_factory)
{
    BL_PROFILE("MLLinOp::defineGrids()");

#ifdef AMREX_USE_EB
    if ( ! a_factory.empty() ) {
        auto const* ebf = dynamic_cast<EBFArrayBoxFactory const*>(a_factory[0]);
        if (ebf && !(ebf->isAllRegular())) { // Has non-trivial EB
            mg_domain_min_width = 4;
        }
    }
#endif

    m_num_amr_levels = 0;
    for (int amrlev = 0; amrlev < a_geom.size(); amrlev++) {
        if (!a_grids[amrlev].empty()) {
            m_num_amr_levels++;
        }
    }

    m_amr_ref_ratio.resize(m_num_amr_levels);
    m_num_mg_levels.resize(m_num_amr_levels);

    m_geom.resize(m_num_amr_levels);
    m_grids.resize(m_num_amr_levels);
    m_dmap.resize(m_num_amr_levels);
    m_factory.resize(m_num_amr_levels);

    m_default_comm = ParallelContext::CommunicatorSub();

    const RealBox& rb = a_geom[0].ProbDomain();
    const int coord = a_geom[0].Coord();
    const Array<int,AMREX_SPACEDIM>& is_per = a_geom[0].isPeriodic();

    IntVect mg_coarsen_ratio_v(mg_coarsen_ratio);
    IntVect mg_box_min_width_v(mg_box_min_width);
    IntVect mg_domain_min_width_v(mg_domain_min_width);
    if (hasHiddenDimension()) {
        AMREX_ASSERT_WITH_MESSAGE(AMREX_SPACEDIM == 3 && m_num_amr_levels == 1,
                                  "Hidden direction only supported for 3d level solve");
        mg_coarsen_ratio_v[info.hidden_direction] = 1;
        mg_box_min_width_v[info.hidden_direction] = 0;
        mg_domain_min_width_v[info.hidden_direction] = 0;
    }

    // fine amr levels
    for (int amrlev = m_num_amr_levels-1; amrlev > 0; --amrlev)
    {
        m_num_mg_levels[amrlev] = 1;
        m_geom[amrlev].push_back(a_geom[amrlev]);
        m_grids[amrlev].push_back(a_grids[amrlev]);
        m_dmap[amrlev].push_back(a_dmap[amrlev]);
        if (amrlev < a_factory.size()) {
            m_factory[amrlev].emplace_back(a_factory[amrlev]->clone());
        } else {
            m_factory[amrlev].push_back(std::make_unique<DefaultFabFactory<FAB>>());
        }

        IntVect rr = mg_coarsen_ratio_v;
        const Box& dom = a_geom[amrlev].Domain();
        for (int i = 0; i < 2; ++i)
        {
            if (!dom.coarsenable(rr)) { amrex::Abort("MLLinOp: Uncoarsenable domain"); }

            const Box& cdom = amrex::coarsen(dom,rr);
            if (cdom == a_geom[amrlev-1].Domain()) { break; }

            ++(m_num_mg_levels[amrlev]);

            m_geom[amrlev].emplace_back(cdom, rb, coord, is_per);

            m_grids[amrlev].push_back(a_grids[amrlev]);
            AMREX_ASSERT(m_grids[amrlev].back().coarsenable(rr));
            m_grids[amrlev].back().coarsen(rr);

            m_dmap[amrlev].push_back(a_dmap[amrlev]);

            rr *= mg_coarsen_ratio_v;
        }

#if (AMREX_SPACEDIM > 1)
        if (hasHiddenDimension()) {
            m_amr_ref_ratio[amrlev-1] = rr[AMREX_SPACEDIM-info.hidden_direction];
        } else
#endif
        {
            m_amr_ref_ratio[amrlev-1] = rr[0];
        }
    }

    // coarsest amr level
    m_num_mg_levels[0] = 1;
    m_geom[0].push_back(a_geom[0]);
    m_grids[0].push_back(a_grids[0]);
    m_dmap[0].push_back(a_dmap[0]);
    if (!a_factory.empty()) {
        m_factory[0].emplace_back(a_factory[0]->clone());
    } else {
        m_factory[0].push_back(std::make_unique<DefaultFabFactory<FAB>>());
    }

    m_domain_covered.resize(m_num_amr_levels, false);
    auto npts0 = m_grids[0][0].numPts();
    m_domain_covered[0] = (npts0 == compactify(m_geom[0][0].Domain()).numPts());
    for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
    {
        if (!m_domain_covered[amrlev-1]) { break; }
        m_domain_covered[amrlev] = (m_grids[amrlev][0].numPts() ==
                                    compactify(m_geom[amrlev][0].Domain()).numPts());
    }

    Box aggbox;
    bool aggable = false;

    if (m_grids[0][0].size() > 1 && info.do_agglomeration)
    {
        if (m_domain_covered[0])
        {
            aggbox = m_geom[0][0].Domain();
            if (hasHiddenDimension()) {
                aggbox.makeSlab(hiddenDirection(), m_grids[0][0][0].smallEnd(hiddenDirection()));
            }
            aggable = true;
        }
        else
        {
            aggbox = m_grids[0][0].minimalBox();
            aggable = (aggbox.numPts() == npts0);
        }
    }

    bool agged = false;
    bool coned = false;
    int agg_lev = 0, con_lev = 0;

    AMREX_ALWAYS_ASSERT( ! (info.do_semicoarsening && info.hasHiddenDimension())
                         && info.semicoarsening_direction >= -1
                         && info.semicoarsening_direction < AMREX_SPACEDIM );

    if (info.do_agglomeration && aggable)
    {
        Box dbx = m_geom[0][0].Domain();
        Box bbx = aggbox;
        Real const nbxs = static_cast<Real>(m_grids[0][0].size());
        Real const threshold_npts = static_cast<Real>(AMREX_D_TERM(info.agg_grid_size,
                                                                  *info.agg_grid_size,
                                                                  *info.agg_grid_size));
        Vector<Box> domainboxes{dbx};
        Vector<Box> boundboxes{bbx};
        Vector<int> agg_flag{false};
        Vector<IntVect> accum_coarsen_ratio{IntVect(1)};
        int numsclevs = 0;

        for (int lev = 0; lev < info.max_coarsening_level; ++lev)
        {
            IntVect rr_level = mg_coarsen_ratio_v;
            bool const do_semicoarsening_level = info.do_semicoarsening
                && numsclevs < info.max_semicoarsening_level;
            if (do_semicoarsening_level
                && info.semicoarsening_direction != -1)
            {
                rr_level[info.semicoarsening_direction] = 1;
            }
            IntVect is_coarsenable;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                IntVect rr_dir(1);
                rr_dir[idim] = rr_level[idim];
                is_coarsenable[idim] = dbx.coarsenable(rr_dir, mg_domain_min_width_v)
                    && bbx.coarsenable(rr_dir, mg_box_min_width_v);
                if (!is_coarsenable[idim] && do_semicoarsening_level
                    && info.semicoarsening_direction == -1)
                {
                    is_coarsenable[idim] = true;
                    rr_level[idim] = 1;
                }
            }
            if (is_coarsenable != IntVect(1) || rr_level == IntVect(1)) {
                break;
            }
            if (do_semicoarsening_level && info.semicoarsening_direction == -1) {
                // make sure there is at most one direction that is not coarsened
                int n_ones = AMREX_D_TERM(  static_cast<int>(rr_level[0] == 1),
                                          + static_cast<int>(rr_level[1] == 1),
                                          + static_cast<int>(rr_level[2] == 1));
                if (n_ones > 1) { break; }
            }
            if (rr_level != mg_coarsen_ratio_v) {
                ++numsclevs;
            }

            accum_coarsen_ratio.push_back(accum_coarsen_ratio.back()*rr_level);
            domainboxes.push_back(dbx.coarsen(rr_level));
            boundboxes.push_back(bbx.coarsen(rr_level));
            bool to_agg = (bbx.d_numPts() / nbxs) < 0.999*threshold_npts;
            agg_flag.push_back(to_agg);
        }

        for (int lev = 1, nlevs = static_cast<int>(domainboxes.size()); lev < nlevs; ++lev) {
            if (!agged && !agg_flag[lev] &&
                a_grids[0].coarsenable(accum_coarsen_ratio[lev], mg_box_min_width_v))
            {
                m_grids[0].push_back(amrex::coarsen(a_grids[0], accum_coarsen_ratio[lev]));
                m_dmap[0].push_back(a_dmap[0]);
            } else {
                IntVect cr = domainboxes[lev-1].length() / domainboxes[lev].length();
                if (!m_grids[0].back().coarsenable(cr)) {
                    break; // average_down would fail if fine boxarray is not coarsenable.
                }
                m_grids[0].emplace_back(boundboxes[lev]);
                IntVect max_grid_size(info.agg_grid_size);
                if (info.do_semicoarsening && info.max_semicoarsening_level >= lev
                    && info.semicoarsening_direction != -1)
                {
                    IntVect blen = amrex::enclosedCells(boundboxes[lev]).size();
                    AMREX_D_TERM(int mgs_0 = (max_grid_size[0]+blen[0]-1) / blen[0];,
                                 int mgs_1 = (max_grid_size[1]+blen[1]-1) / blen[1];,
                                 int mgs_2 = (max_grid_size[2]+blen[2]-1) / blen[2]);
                    max_grid_size[info.semicoarsening_direction]
                        *= AMREX_D_TERM(mgs_0, *mgs_1, *mgs_2);
                }
                m_grids[0].back().maxSize(max_grid_size);
                m_dmap[0].push_back(DistributionMapping());
                if (!agged) {
                    agged = true;
                    agg_lev = lev;
                }
            }
            m_geom[0].emplace_back(domainboxes[lev],rb,coord,is_per);
        }
    }
    else
    {
        Long consolidation_threshold = 0;
        Real avg_npts = 0.0;
        if (info.do_consolidation) {
            avg_npts = static_cast<Real>(a_grids[0].d_numPts()) / static_cast<Real>(ParallelContext::NProcsSub());
            consolidation_threshold = AMREX_D_TERM(Long(info.con_grid_size),
                                                       *info.con_grid_size,
                                                       *info.con_grid_size);
        }

        Box const& dom0 = a_geom[0].Domain();
        IntVect rr_vec(1);
        int numsclevs = 0;
        for (int lev = 0; lev < info.max_coarsening_level; ++lev)
        {
            IntVect rr_level = mg_coarsen_ratio_v;
            bool do_semicoarsening_level = info.do_semicoarsening
                && numsclevs < info.max_semicoarsening_level;
            if (do_semicoarsening_level
                && info.semicoarsening_direction != -1)
            {
                rr_level[info.semicoarsening_direction] = 1;
            }
            IntVect is_coarsenable;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                IntVect rr_dir(1);
                rr_dir[idim] = rr_vec[idim] * rr_level[idim];
                is_coarsenable[idim] = dom0.coarsenable(rr_dir, mg_domain_min_width_v)
                    && a_grids[0].coarsenable(rr_dir, mg_box_min_width_v);
                if (!is_coarsenable[idim] && do_semicoarsening_level
                    && info.semicoarsening_direction == -1)
                {
                    is_coarsenable[idim] = true;
                    rr_level[idim] = 1;
                }
            }
            if (is_coarsenable != IntVect(1) || rr_level == IntVect(1)) {
                break;
            }
            if (do_semicoarsening_level && info.semicoarsening_direction == -1) {
                // make sure there is at most one direction that is not coarsened
                int n_ones = AMREX_D_TERM(  static_cast<int>(rr_level[0] == 1),
                                          + static_cast<int>(rr_level[1] == 1),
                                          + static_cast<int>(rr_level[2] == 1));
                if (n_ones > 1) { break; }
            }
            if (rr_level != mg_coarsen_ratio_v) {
                ++numsclevs;
            }
            rr_vec *= rr_level;

            m_geom[0].emplace_back(amrex::coarsen(dom0, rr_vec), rb, coord, is_per);
            m_grids[0].push_back(amrex::coarsen(a_grids[0], rr_vec));

            if (info.do_consolidation)
            {
                if (avg_npts/static_cast<Real>(AMREX_D_TERM(rr_vec[0], *rr_vec[1], *rr_vec[2]))
                    < Real(0.999)*static_cast<Real>(consolidation_threshold))
                {
                    coned = true;
                    con_lev = m_dmap[0].size();
                    m_dmap[0].push_back(DistributionMapping());
                }
                else
                {
                    m_dmap[0].push_back(m_dmap[0].back());
                }
            }
            else
            {
                m_dmap[0].push_back(a_dmap[0]);
            }
        }
    }

    m_num_mg_levels[0] = m_grids[0].size();

    for (int mglev = 0; mglev < m_num_mg_levels[0] - 1; mglev++){
        const Box& fine_domain = m_geom[0][mglev].Domain();
        const Box& crse_domain = m_geom[0][mglev+1].Domain();
        mg_coarsen_ratio_vec.push_back(fine_domain.length()/crse_domain.length());
    }

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev) {
        if (AMRRefRatio(amrlev) == 4 && mg_coarsen_ratio_vec.empty()) {
            mg_coarsen_ratio_vec.push_back(IntVect(2));
        }
    }

    if (agged)
    {
        makeAgglomeratedDMap(m_grids[0], m_dmap[0]);
    }
    else if (coned)
    {
        makeConsolidatedDMap(m_grids[0], m_dmap[0], info.con_ratio, info.con_strategy);
    }

    if (agged || coned)
    {
        m_bottom_comm = makeSubCommunicator(m_dmap[0].back());
    }
    else
    {
        m_bottom_comm = m_default_comm;
    }

    m_do_agglomeration = agged;
    m_do_consolidation = coned;

    if (verbose > 1) {
        if (agged) {
            Print() << "MLLinOp::defineGrids(): agglomerated AMR level 0 starting at MG level "
                    << agg_lev << " of " << m_num_mg_levels[0] << "\n";
        } else if (coned) {
            Print() << "MLLinOp::defineGrids(): consolidated AMR level 0 starting at MG level "
                    << con_lev << " of " << m_num_mg_levels[0]
                    << " (ratio = " << info.con_ratio << ")" << "\n";
        } else {
            Print() << "MLLinOp::defineGrids(): no agglomeration or consolidation of AMR level 0\n";
        }
    }

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_factory[amrlev].emplace_back(makeFactory(amrlev,mglev));
        }
    }

    for (int amrlev = 1; amrlev < m_num_amr_levels; ++amrlev)
    {
        AMREX_ASSERT_WITH_MESSAGE(m_grids[amrlev][0].coarsenable(m_amr_ref_ratio[amrlev-1]),
                                  "MLLinOp: grids not coarsenable between AMR levels");
    }
}

template <typename MF>
void
MLLinOpT<MF>::defineBC ()
{
    m_needs_coarse_data_for_bc = !m_domain_covered[0];

    levelbc_raii.resize(m_num_amr_levels);
    robin_a_raii.resize(m_num_amr_levels);
    robin_b_raii.resize(m_num_amr_levels);
    robin_f_raii.resize(m_num_amr_levels);
}

template <typename MF>
void
MLLinOpT<MF>::setDomainBC (const Array<BCType,AMREX_SPACEDIM>& a_lobc,
                           const Array<BCType,AMREX_SPACEDIM>& a_hibc) noexcept
{
    const int ncomp = getNComp();
    setDomainBC(Vector<Array<BCType,AMREX_SPACEDIM> >(ncomp,a_lobc),
                Vector<Array<BCType,AMREX_SPACEDIM> >(ncomp,a_hibc));
}

template <typename MF>
void
MLLinOpT<MF>::setDomainBC (const Vector<Array<BCType,AMREX_SPACEDIM> >& a_lobc,
                           const Vector<Array<BCType,AMREX_SPACEDIM> >& a_hibc) noexcept
{
    const int ncomp = getNComp();
    AMREX_ASSERT_WITH_MESSAGE(ncomp == a_lobc.size() && ncomp == a_hibc.size(),
                              "MLLinOp::setDomainBC: wrong size");
    m_lobc = a_lobc;
    m_hibc = a_hibc;
    m_lobc_orig = m_lobc;
    m_hibc_orig = m_hibc;
    for (int icomp = 0; icomp < ncomp; ++icomp) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (m_geom[0][0].isPeriodic(idim)) {
                AMREX_ALWAYS_ASSERT(m_lobc[icomp][idim] == BCType::Periodic &&
                                    m_hibc[icomp][idim] == BCType::Periodic);
            } else {
                AMREX_ALWAYS_ASSERT(m_lobc[icomp][idim] != BCType::Periodic &&
                                    m_hibc[icomp][idim] != BCType::Periodic);
            }

            if (m_lobc[icomp][idim] == LinOpBCType::inhomogNeumann ||
                m_lobc[icomp][idim] == LinOpBCType::Robin)
            {
                m_lobc[icomp][idim] = LinOpBCType::Neumann;
            }

            if (m_hibc[icomp][idim] == LinOpBCType::inhomogNeumann ||
                m_hibc[icomp][idim] == LinOpBCType::Robin)
            {
                m_hibc[icomp][idim] = LinOpBCType::Neumann;
            }
        }
    }

    if (hasHiddenDimension()) {
        const int hd = hiddenDirection();
        for (int n = 0; n < ncomp; ++n) {
            m_lobc[n][hd] = LinOpBCType::Neumann;
            m_hibc[n][hd] = LinOpBCType::Neumann;
        }
    }

    if (hasInhomogNeumannBC() && !supportInhomogNeumannBC()) {
        amrex::Abort("Inhomogeneous Neumann BC not supported");
    }
    if (hasRobinBC() && !supportRobinBC()) {
        amrex::Abort("Robin BC not supported");
    }
}

template <typename MF>
bool
MLLinOpT<MF>::hasBC (BCType bct) const noexcept
{
    int ncomp = m_lobc_orig.size();
    for (int n = 0; n < ncomp; ++n) {
        for (int idim = 0; idim <AMREX_SPACEDIM; ++idim) {
            if (m_lobc_orig[n][idim] == bct || m_hibc_orig[n][idim] == bct) {
                return true;
            }
        }
    }
    return false;
}

template <typename MF>
bool
MLLinOpT<MF>::hasInhomogNeumannBC () const noexcept
{
    return hasBC(BCType::inhomogNeumann);
}

template <typename MF>
bool
MLLinOpT<MF>::hasRobinBC () const noexcept
{
    return hasBC(BCType::Robin);
}

template <typename MF>
Box
MLLinOpT<MF>::compactify (Box const& b) const noexcept
{
#if (AMREX_SPACEDIM == 3)
    if (info.hasHiddenDimension()) {
        const auto& lo = b.smallEnd();
        const auto& hi = b.bigEnd();
        if (info.hidden_direction == 0) {
            return Box(IntVect(lo[1],lo[2],0), IntVect(hi[1],hi[2],0), b.ixType());
        } else if (info.hidden_direction == 1) {
            return Box(IntVect(lo[0],lo[2],0), IntVect(hi[0],hi[2],0), b.ixType());
        } else {
            return Box(IntVect(lo[0],lo[1],0), IntVect(hi[0],hi[1],0), b.ixType());
        }
    } else
#endif
    {
        return b;
    }
}

template <typename MF>
void
MLLinOpT<MF>::makeAgglomeratedDMap (const Vector<BoxArray>& ba,
                                    Vector<DistributionMapping>& dm)
{
    BL_PROFILE("MLLinOp::makeAgglomeratedDMap");

    BL_ASSERT(!dm[0].empty());
    for (int i = 1, N=static_cast<int>(ba.size()); i < N; ++i)
    {
        if (dm[i].empty())
        {
            const std::vector< std::vector<int> >& sfc = DistributionMapping::makeSFC(ba[i]);

            const int nprocs = ParallelContext::NProcsSub();
            AMREX_ASSERT(static_cast<int>(sfc.size()) == nprocs);

            Vector<int> pmap(ba[i].size());
            for (int iproc = 0; iproc < nprocs; ++iproc) {
                int grank = ParallelContext::local_to_global_rank(iproc);
                for (int ibox : sfc[iproc]) {
                    pmap[ibox] = grank;
                }
            }
            dm[i].define(std::move(pmap));
        }
    }
}

template <typename MF>
void
MLLinOpT<MF>::makeConsolidatedDMap (const Vector<BoxArray>& ba,
                                    Vector<DistributionMapping>& dm,
                                    int ratio, int strategy)
{
    BL_PROFILE("MLLinOp::makeConsolidatedDMap()");

    int factor = 1;
    BL_ASSERT(!dm[0].empty());
    for (int i = 1, N=static_cast<int>(ba.size()); i < N; ++i)
    {
        if (dm[i].empty())
        {
            factor *= ratio;

            const int nprocs = ParallelContext::NProcsSub();
            const auto& pmap_fine = dm[i-1].ProcessorMap();
            Vector<int> pmap(pmap_fine.size());
            ParallelContext::global_to_local_rank(pmap.data(), pmap_fine.data(), static_cast<int>(pmap.size()));
            if (strategy == 1) {
                for (auto& x: pmap) {
                    x /= ratio;
                }
            } else if (strategy == 2) {
                int nprocs_con = static_cast<int>(std::ceil(static_cast<Real>(nprocs)
                                                            / static_cast<Real>(factor)));
                for (auto& x: pmap) {
                    auto d = std::div(x,nprocs_con);
                    x = d.rem;
                }
            } else if (strategy == 3) {
                if (factor == ratio) {
                    const std::vector< std::vector<int> >& sfc = DistributionMapping::makeSFC(ba[i]);
                    for (int iproc = 0; iproc < nprocs; ++iproc) {
                        for (int ibox : sfc[iproc]) {
                            pmap[ibox] = iproc;
                        }
                    }
                }
                for (auto& x: pmap) {
                    x /= ratio;
                }
            }

            if (ParallelContext::CommunicatorSub() == ParallelDescriptor::Communicator()) {
                dm[i].define(std::move(pmap));
            } else {
                Vector<int> pmap_g(pmap.size());
                ParallelContext::local_to_global_rank(pmap_g.data(), pmap.data(), static_cast<int>(pmap.size()));
                dm[i].define(std::move(pmap_g));
            }
        }
    }
}

template <typename MF>
MPI_Comm
MLLinOpT<MF>::makeSubCommunicator (const DistributionMapping& dm)
{
    BL_PROFILE("MLLinOp::makeSubCommunicator()");

#ifdef BL_USE_MPI

    Vector<int> newgrp_ranks = dm.ProcessorMap();
    std::sort(newgrp_ranks.begin(), newgrp_ranks.end());
    auto last = std::unique(newgrp_ranks.begin(), newgrp_ranks.end());
    newgrp_ranks.erase(last, newgrp_ranks.end());

    MPI_Comm newcomm;
    MPI_Group defgrp, newgrp;
    MPI_Comm_group(m_default_comm, &defgrp);
    if (ParallelContext::CommunicatorSub() == ParallelDescriptor::Communicator()) {
        MPI_Group_incl(defgrp, static_cast<int>(newgrp_ranks.size()), newgrp_ranks.data(), &newgrp);
    } else {
        Vector<int> local_newgrp_ranks(newgrp_ranks.size());
        ParallelContext::global_to_local_rank(local_newgrp_ranks.data(),
                                              newgrp_ranks.data(), static_cast<int>(newgrp_ranks.size()));
        MPI_Group_incl(defgrp, static_cast<int>(local_newgrp_ranks.size()), local_newgrp_ranks.data(), &newgrp);
    }

    MPI_Comm_create(m_default_comm, newgrp, &newcomm);

    m_raii_comm = std::make_unique<CommContainer>(newcomm);

    MPI_Group_free(&defgrp);
    MPI_Group_free(&newgrp);

    return newcomm;
#else
    amrex::ignore_unused(dm);
    return m_default_comm;
#endif
}

template <typename MF>
void
MLLinOpT<MF>::setDomainBCLoc (const Array<Real,AMREX_SPACEDIM>& lo_bcloc,
                              const Array<Real,AMREX_SPACEDIM>& hi_bcloc) noexcept
{
    m_domain_bloc_lo = lo_bcloc;
    m_domain_bloc_hi = hi_bcloc;
}

template <typename MF>
void
MLLinOpT<MF>::setCoarseFineBC (const MF* crse, int crse_ratio,
                               LinOpBCType bc_type) noexcept
{
    setCoarseFineBC(crse, IntVect(crse_ratio), bc_type);
}

template <typename MF>
void
MLLinOpT<MF>::setCoarseFineBC (const MF* crse, IntVect const& crse_ratio,
                               LinOpBCType bc_type) noexcept
{
    m_coarse_data_for_bc = crse;
    m_coarse_data_crse_ratio = crse_ratio;
    m_coarse_fine_bc_type = bc_type;
}

template <typename MF>
template <typename AMF, std::enable_if_t<!std::is_same_v<MF,AMF>,int>>
void
MLLinOpT<MF>::setCoarseFineBC (const AMF* crse, int crse_ratio,
                               LinOpBCType bc_type) noexcept
{
    setCoarseFineBC(crse, IntVect(crse_ratio), bc_type);
}

template <typename MF>
template <typename AMF, std::enable_if_t<!std::is_same_v<MF,AMF>,int>>
void
MLLinOpT<MF>::setCoarseFineBC (const AMF* crse, IntVect const& crse_ratio,
                               LinOpBCType bc_type) noexcept
{
    m_coarse_data_for_bc_raii = MF(crse->boxArray(), crse->DistributionMap(),
                                   crse->nComp(), crse->nGrowVect());
    m_coarse_data_for_bc_raii.LocalCopy(*crse, 0, 0, crse->nComp(),
                                        crse->nGrowVect());
    m_coarse_data_for_bc = &m_coarse_data_for_bc_raii;
    m_coarse_data_crse_ratio = crse_ratio;
    m_coarse_fine_bc_type = bc_type;
}

template <typename MF>
void
MLLinOpT<MF>::make (Vector<Vector<MF> >& mf, IntVect const& ng) const
{
    mf.clear();
    mf.resize(m_num_amr_levels);
    for (int alev = 0; alev < m_num_amr_levels; ++alev) {
        mf[alev].resize(m_num_mg_levels[alev]);
        for (int mlev = 0; mlev < m_num_mg_levels[alev]; ++mlev) {
            mf[alev][mlev] = make(alev, mlev, ng);
        }
    }
}

template <typename MF>
MF
MLLinOpT<MF>::make (int amrlev, int mglev, IntVect const& ng) const
{
    if constexpr (IsMultiFabLike_v<MF>) {
        return MF(amrex::convert(m_grids[amrlev][mglev], m_ixtype),
                  m_dmap[amrlev][mglev], getNComp(), ng, MFInfo(),
                  *m_factory[amrlev][mglev]);
    } else {
        amrex::ignore_unused(amrlev, mglev, ng);
        amrex::Abort("MLLinOpT::make: how did we get here?");
        return {};
    }
}

template <typename MF>
MF
MLLinOpT<MF>::makeAlias (MF const& mf) const
{
    if constexpr (IsMultiFabLike_v<MF>) {
        return MF(mf, amrex::make_alias, 0, mf.nComp());
    } else {
        amrex::ignore_unused(mf);
        amrex::Abort("MLLinOpT::makeAlias: how did we get here?");
        return {};
    }
}

template <typename MF>
MF
MLLinOpT<MF>::makeCoarseMG (int amrlev, int mglev, IntVect const& ng) const
{
    if constexpr (IsMultiFabLike_v<MF>) {
        BoxArray cba = m_grids[amrlev][mglev];
        IntVect ratio = (amrlev > 0) ? IntVect(2) : mg_coarsen_ratio_vec[mglev];
        cba.coarsen(ratio);
        cba.convert(m_ixtype);
        return MF(cba, m_dmap[amrlev][mglev], getNComp(), ng);
    } else {
        amrex::ignore_unused(amrlev, mglev, ng);
        amrex::Abort("MLLinOpT::makeCoarseMG: how did we get here?");
        return {};
    }
}

template <typename MF>
MF
MLLinOpT<MF>::makeCoarseAmr (int famrlev, IntVect const& ng) const
{
    if constexpr (IsMultiFabLike_v<MF>) {
        BoxArray cba = m_grids[famrlev][0];
        IntVect ratio(AMRRefRatio(famrlev-1));
        cba.coarsen(ratio);
        cba.convert(m_ixtype);
        return MF(cba, m_dmap[famrlev][0], getNComp(), ng);
    } else {
        amrex::ignore_unused(famrlev, ng);
        amrex::Abort("MLLinOpT::makeCoarseAmr: how did we get here?");
        return {};
    }
}

template <typename MF>
void
MLLinOpT<MF>::resizeMultiGrid (int new_size)
{
    if (new_size <= 0 || new_size >= m_num_mg_levels[0]) { return; }

    m_num_mg_levels[0] = new_size;

    m_geom[0].resize(new_size);
    m_grids[0].resize(new_size);
    m_dmap[0].resize(new_size);
    m_factory[0].resize(new_size);

    if (m_bottom_comm != m_default_comm) {
        m_bottom_comm = makeSubCommunicator(m_dmap[0].back());
    }
}

template <typename MF>
void
MLLinOpT<MF>::avgDownResMG (int clev, MF& cres, MF const& fres) const
{
    amrex::ignore_unused(clev, cres, fres);
    if constexpr (amrex::IsFabArray<MF>::value) {
        const int ncomp = this->getNComp();
#ifdef AMREX_USE_EB
        if (!fres.isAllRegular()) {
            if constexpr (std::is_same<MF,MultiFab>()) {
                amrex::EB_average_down(fres, cres, 0, ncomp,
                                       mg_coarsen_ratio_vec[clev-1]);
            } else {
                amrex::Abort("EB_average_down only works with MultiFab");
            }
        } else
#endif
        {
            amrex::average_down(fres, cres, 0, ncomp, mg_coarsen_ratio_vec[clev-1]);
        }
    } else {
        amrex::Abort("For non-FabArray, MLLinOpT<MF>::avgDownResMG should be overridden.");
    }
}

template <typename MF>
bool
MLLinOpT<MF>::isMFIterSafe (int amrlev, int mglev1, int mglev2) const
{
    return m_dmap[amrlev][mglev1] == m_dmap[amrlev][mglev2]
        && BoxArray::SameRefs(m_grids[amrlev][mglev1], m_grids[amrlev][mglev2]);
}

template <typename MF>
template <typename AMF, std::enable_if_t<!std::is_same_v<MF,AMF>,int>>
void
MLLinOpT<MF>::setLevelBC (int amrlev, const AMF* levelbcdata,
                          const AMF* robinbc_a, const AMF* robinbc_b,
                          const AMF* robinbc_f)
{
    const int ncomp = this->getNComp();
    if (levelbcdata) {
        levelbc_raii[amrlev] = std::make_unique<MF>(levelbcdata->boxArray(),
                                                    levelbcdata->DistributionMap(),
                                                    ncomp, levelbcdata->nGrowVect());
        levelbc_raii[amrlev]->LocalCopy(*levelbcdata, 0, 0, ncomp,
                                        levelbcdata->nGrowVect());
    } else {
        levelbc_raii[amrlev].reset();
    }

    if (robinbc_a) {
        robin_a_raii[amrlev] = std::make_unique<MF>(robinbc_a->boxArray(),
                                                    robinbc_a->DistributionMap(),
                                                    ncomp, robinbc_a->nGrowVect());
        robin_a_raii[amrlev]->LocalCopy(*robinbc_a, 0, 0, ncomp,
                                        robinbc_a->nGrowVect());
    } else {
        robin_a_raii[amrlev].reset();
    }

    if (robinbc_b) {
        robin_b_raii[amrlev] = std::make_unique<MF>(robinbc_b->boxArray(),
                                                    robinbc_b->DistributionMap(),
                                                    ncomp, robinbc_b->nGrowVect());
        robin_b_raii[amrlev]->LocalCopy(*robinbc_b, 0, 0, ncomp,
                                        robinbc_b->nGrowVect());
    } else {
        robin_b_raii[amrlev].reset();
    }

    if (robinbc_f) {
        robin_f_raii[amrlev] = std::make_unique<MF>(robinbc_f->boxArray(),
                                                    robinbc_f->DistributionMap(),
                                                    ncomp, robinbc_f->nGrowVect());
        robin_f_raii[amrlev]->LocalCopy(*robinbc_f, 0, 0, ncomp,
                                        robinbc_f->nGrowVect());
    } else {
        robin_f_raii[amrlev].reset();
    }

    this->setLevelBC(amrlev, levelbc_raii[amrlev].get(), robin_a_raii[amrlev].get(),
                     robin_b_raii[amrlev].get(), robin_f_raii[amrlev].get());
}

template <typename MF>
auto
MLLinOpT<MF>::dotProductPrecond (Vector<MF const*> const& x, Vector<MF const*> const& y) const -> RT
{
    AMREX_ALWAYS_ASSERT(NAMRLevels() == 1);
    return xdoty(0,0,*x[0],*y[0],false);
}

template <typename MF>
auto
MLLinOpT<MF>::norm2Precond (Vector<MF const*> const& x) const -> RT
{
    AMREX_ALWAYS_ASSERT(NAMRLevels() == 1);
    auto r = xdoty(0,0,*x[0],*x[0],false);
    return std::sqrt(r);
}

extern template class MLLinOpT<MultiFab>;

using MLLinOp = MLLinOpT<MultiFab>;

}

#endif
