#ifndef AMREX_ML_LINOP_H_
#define AMREX_ML_LINOP_H_
#include <AMReX_Config.H>

#include <AMReX_Any.H>
#include <AMReX_SPACE.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BndryRegister.H>
#include <AMReX_YAFluxRegister.H>
#include <AMReX_MLMGBndry.H>
#include <AMReX_VisMF.H>

#ifdef AMREX_USE_EB
#include <AMReX_MultiCutFab.H>
#endif

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
#include <AMReX_Hypre.H>
#include <AMReX_HypreNodeLap.H>
#endif

namespace amrex {

enum class BottomSolver : int {
    Default, smoother, bicgstab, cg, bicgcg, cgbicg, hypre, petsc
};

#ifdef AMREX_USE_PETSC
class PETScABecLap;
#endif

class MLMG;

struct LPInfo
{
    bool do_agglomeration = true;
    bool do_consolidation = true;
    bool do_semicoarsening = false;
    int agg_grid_size = -1;
    int con_grid_size = -1;
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
    LPInfo& setMetricTerm (bool x) noexcept { has_metric_term = x; return *this; }
    LPInfo& setMaxCoarseningLevel (int n) noexcept { max_coarsening_level = n; return *this; }
    LPInfo& setMaxSemicoarseningLevel (int n) noexcept { max_semicoarsening_level = n; return *this; }
    LPInfo& setSemicoarseningDirection (int n) noexcept { semicoarsening_direction = n; return *this; }
    LPInfo& setHiddenDirection (int n) noexcept { hidden_direction = n; return *this; }

    bool hasHiddenDimension () const noexcept {
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

class MLLinOp
{
public:

    friend class MLMG;
    friend class MLCGSolver;
    friend class MLPoisson;
    friend class MLABecLaplacian;

    enum struct BCMode { Homogeneous, Inhomogeneous };
    using BCType = LinOpBCType;

    enum struct StateMode { Solution, Correction };

    enum struct Location { FaceCenter, FaceCentroid, CellCenter, CellCentroid };

    static void Initialize ();
    static void Finalize ();

    MLLinOp ();
    virtual ~MLLinOp ();

    MLLinOp (const MLLinOp&) = delete;
    MLLinOp (MLLinOp&&) = delete;
    MLLinOp& operator= (const MLLinOp&) = delete;
    MLLinOp& operator= (MLLinOp&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info,
                 const Vector<FabFactory<FArrayBox> const*>& a_factory,
                 bool eb_limit_coarsening = true);

    virtual std::string name () const { return std::string("Unspecified"); }

    /**
    * \brief Boundary of the whole domain. This functions must be called,
    * and must be called before other bc functions.
    *
    * \param lobc
    * \param hibc
    */
    // For single component
    void setDomainBC (const Array<BCType,AMREX_SPACEDIM>& lobc,
                      const Array<BCType,AMREX_SPACEDIM>& hibc) noexcept;

    // For multi-component
    void setDomainBC (const Vector<Array<BCType,AMREX_SPACEDIM> >& lobc,
                      const Vector<Array<BCType,AMREX_SPACEDIM> >& hibc) noexcept;

    // By default, domain BC is on the domain face.  If that's the
    // case, this function doesn't need to be called.  However, one
    // could use this function to set non-zero domain BC locations.
    // Note all values should be >= 0.  If this function is called,
    // it MUST be called before setLevelBC.
    void setDomainBCLoc (const Array<Real,AMREX_SPACEDIM>& lo_bcloc,
                         const Array<Real,AMREX_SPACEDIM>& hi_bcloc) noexcept;

    /**
    * \brief Needs coarse data for bc?  If the lowest level grids does not
    * cover the entire domain, coarse level data are needed for
    * supplying Dirichlet bc at coarse/fine boundary, even when the
    * domain bc is not Dirichlet.
    */
    bool needsCoarseDataForBC () const noexcept { return m_needs_coarse_data_for_bc; }


    /**
    * \brief For cell-centered solves only: if we want to do a linear solve
    * where the boundary conditions on the coarsest AMR level of the solve
    * come from a coarser level (e.g. the base AMR level of the solve is > 0
    * and does not cover the entire domain), we must explicitly provide the coarser
    * data.  Boundary conditions from a coarser level are always Dirichlet.
    * The MultiFab crse does not need to have ghost cells and is at a coarser
    * resolution than the coarsest AMR level of the solve; it is used to supply
    * (interpolated) boundary conditions for the solve.
    * NOTE: If this is called, it must be called before `setLevelBC`.
    * If crse is nullptr, then the bc values are assumed to be zero.
    *
    * \param crse
    * \param crse_ratio
    */
    void setCoarseFineBC (const MultiFab* crse, int crse_ratio) noexcept;

    /**
    * \brief For cell-centered solves only: this must be called for each
    * level.  Argument `levelbcdata` is used to supply Dirichlet or Neumann
    * bc at the physical domain; if those data are homogeneous we can pass
    * nullptr instead of levelbcdata.  Regardless, this function must be
    * called.  If used, the MultiFab levelbcdata must have one ghost cell.
    * Only the data outside the physical domain will be used.  It is assumed
    * that the data in those ghost cells outside the domain live exactly on
    * the face of the physical domain.  Argument `amrlev` is relative level
    * such that the lowest to the solver is always 0.  The optional
    * arguments robinbc_[a|b|f] provide Robin boundary condition `a*phi +
    * b*dphi/dn = f`.  Note that `d./dn` is `d./dx` at the upper boundary
    * and `-d./dx` at the lower boundary, for Robin BC. However, for
    * inhomogeneous Neumann BC, the value in leveldata is assumed to be
    * `d./dx`.
    */
    virtual void setLevelBC (int /*amrlev*/, const MultiFab* /*levelbcdata*/,
                             const MultiFab* /*robinbc_a*/ = nullptr,
                             const MultiFab* /*robinbc_b*/ = nullptr,
                             const MultiFab* /*robinbc_f*/ = nullptr) {}

    void setVerbose (int v) noexcept { verbose = v; }

    void setMaxOrder (int o) noexcept { maxorder = o; }
    int getMaxOrder () const noexcept { return maxorder; }

    void setEnforceSingularSolvable (bool o) noexcept { enforceSingularSolvable = o; }
    bool getEnforceSingularSolvable () const noexcept { return enforceSingularSolvable; }

    virtual BottomSolver getDefaultBottomSolver () const { return BottomSolver::bicgstab; }
    virtual int getNComp () const { return 1; }
    virtual int getNGrow (int /*a_lev*/ = 0, int /*mg_lev*/ = 0) const { return 0; }

    virtual bool needsUpdate () const { return false; }
    virtual void update () {}

    virtual void restriction (int /*amrlev*/, int /*cmglev*/, MultiFab& /*crse*/, MultiFab& /*fine*/) const {}
    virtual void interpolation (int /*amrlev*/, int /*fmglev*/, MultiFab& /*fine*/, const MultiFab& /*crse*/) const {}
    virtual void interpAssign (int /*amrlev*/, int /*fmglev*/, MultiFab& /*fine*/, MultiFab& /*crse*/) const {}
    virtual void averageDownSolutionRHS (int /*camrlev*/, MultiFab& /*crse_sol*/, MultiFab& /*crse_rhs*/,
                                         const MultiFab& /*fine_sol*/, const MultiFab& /*fine_rhs*/) {}

    virtual void apply (int /*amrlev*/, int /*mglev*/, MultiFab& /*out*/, MultiFab& /*in*/, BCMode /*bc_mode*/,
                        StateMode /*s_mode*/, const MLMGBndry* /*bndry*/=nullptr) const {}
    virtual void smooth (int /*amrlev*/, int /*mglev*/, MultiFab& /*sol*/, const MultiFab& /*rhs*/,
                         bool /*skip_fillboundary*/=false) const {}

    // Divide mf by the diagonal component of the operator. Used by bicgstab.
    virtual void normalize (int /*amrlev*/, int /*mglev*/, MultiFab& /*mf*/) const {}

    virtual void solutionResidual (int /*amrlev*/, MultiFab& /*resid*/, MultiFab& /*x*/, const MultiFab& /*b*/,
                                   const MultiFab* /*crse_bcdata*/=nullptr) {}
    virtual void correctionResidual (int /*amrlev*/, int /*mglev*/, MultiFab& /*resid*/, MultiFab& /*x*/, const MultiFab& /*b*/,
                                     BCMode /*bc_mode*/, const MultiFab* /*crse_bcdata*/=nullptr) {}

    virtual void reflux (int /*crse_amrlev*/,
                         MultiFab& /*res*/, const MultiFab& /*crse_sol*/, const MultiFab& /*crse_rhs*/,
                         MultiFab& /*fine_res*/, MultiFab& /*fine_sol*/, const MultiFab& /*fine_rhs*/) const {}
    virtual void compFlux (int /*amrlev*/, const Array<MultiFab*,AMREX_SPACEDIM>& /*fluxes*/,
                           MultiFab& /*sol*/, Location /*loc*/) const {}
    virtual void compGrad (int /*amrlev*/, const Array<MultiFab*,AMREX_SPACEDIM>& /*grad*/,
                           MultiFab& /*sol*/, Location /*loc*/) const {}

    virtual void applyMetricTerm (int /*amrlev*/, int /*mglev*/, Any& /*rhs*/) const {}
    virtual void unapplyMetricTerm (int /*amrlev*/, int /*mglev*/, MultiFab& /*rhs*/) const {}

    virtual void unimposeNeumannBC (int /*amrlev*/, Any& /*rhs*/) const {} // only nodal solver might need it
    virtual void applyInhomogNeumannTerm (int /*amrlev*/, Any& /*rhs*/) const {}
    virtual void applyOverset (int /*amlev*/, Any& /*rhs*/) const {}
    virtual void scaleRHS (int /*amrlev*/, Any& /*rhs*/) const {}
    virtual Vector<Real> getSolvabilityOffset (int /*amrlev*/, int /*mglev*/,
                                               Any const& /*rhs*/) const { return {}; }
    virtual void fixSolvabilityByOffset (int /*amrlev*/, int /*mglev*/, Any& /*rhs*/,
                                         Vector<Real> const& /*offset*/) const {}

    virtual void prepareForSolve () = 0;
    virtual bool isSingular (int /*amrlev*/) const { return false; }
    virtual bool isBottomSingular () const { return false; }
    virtual Real xdoty (int /*amrlev*/, int /*mglev*/, const MultiFab& /*x*/, const MultiFab& /*y*/, bool /*local*/) const { return 0._rt; }

    virtual std::unique_ptr<MLLinOp> makeNLinOp (int /*grid_size*/) const { return {nullptr}; }

    virtual void getFluxes (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& /*a_flux*/,
                            const Vector<MultiFab*>& /*a_sol*/,
                            Location /*a_loc*/) const {
        amrex::Abort("MLLinOp::getFluxes: How did we get here?");
    }
    virtual void getFluxes (const Vector<MultiFab*>& /*a_flux*/,
                            const Vector<MultiFab*>& /*a_sol*/) const {
        amrex::Abort("MLLinOp::getFluxes: How did we get here?");
    }

#ifdef AMREX_USE_EB
    virtual void getEBFluxes (const Vector<MultiFab*>& /*a_flux*/,
                              const Vector<MultiFab*>& /*a_sol*/) const {
        amrex::Abort("MLLinOp::getEBFluxes: How did we get here?");
    }
#endif

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)
    virtual std::unique_ptr<Hypre> makeHypre (Hypre::Interface /*hypre_interface*/) const {
        amrex::Abort("MLLinOp::makeHypre: How did we get here?");
        return {nullptr};
    }
    virtual std::unique_ptr<HypreNodeLap> makeHypreNodeLap(
        int /*bottom_verbose*/,
        const std::string& /* options_namespace */) const
    {
        amrex::Abort("MLLinOp::makeHypreNodeLap: How did we get here?");
        return {nullptr};
    }
#endif

#ifdef AMREX_USE_PETSC
    virtual std::unique_ptr<PETScABecLap> makePETSc () const;
#endif

    virtual bool supportNSolve () const { return false; }

    virtual void copyNSolveSolution (MultiFab&, MultiFab const&) const {}

    virtual Any AnyMake (int amrlev, int mglev, IntVect const& ng) const;
    virtual Any AnyMakeCoarseMG (int amrlev, int mglev, IntVect const& ng) const;
    virtual Any AnyMakeCoarseAmr (int famrlev, IntVect const& ng) const;
    virtual Any AnyMakeAlias (Any const& a) const;
    virtual IntVect AnyGrowVect (Any const& a) const;
    virtual void AnyCopy (Any& dst, Any const& src, IntVect const& ng) const;
    virtual void AnyAdd (Any& dst, Any const& src, IntVect const& ng) const;
    virtual void AnySetToZero (Any& a) const;
    virtual void AnySetBndryToZero (Any& a) const;
#ifdef AMREX_USE_EB
    virtual void AnySetCoveredToZero (Any& a) const;
#endif
    virtual void AnyParallelCopy (Any& dst, Any const& src,
                                  IntVect const& src_nghost, IntVect const& dst_nghost,
                                  Periodicity const& period = Periodicity::NonPeriodic()) const;

    virtual Real AnyNormInf (Any& a) const;

    virtual Real AnyNormInfMask (int amrlev, Any const& a, bool local) const = 0;

    virtual void AnySolutionResidual (int amrlev, Any& resid, Any& x, Any const& b,
                                      Any const* crse_bcdata = nullptr);
    virtual void AnyCorrectionResidual (int amrlev, int mglev, Any& resid, Any& x,
                                        const Any& b, BCMode bc_mode,
                                        const Any* crse_bcdata=nullptr);
    virtual void AnyReflux (int crse_amrlev,
                            Any& res, const Any& crse_sol, const Any& crse_rhs,
                            Any& fine_res, Any& fine_sol, const Any& fine_rhs);

    virtual void AnyAvgDownResAmr (int clev, Any& cres, Any const& fres) const = 0;
    virtual void AnyAvgDownResMG (int clev, Any& cres, Any const& fres) const;

    virtual void AnySmooth (int amrlev, int mglev, Any& sol, const Any& rhs,
                            bool skip_fillboundary=false) const;

    virtual void AnyRestriction (int amrlev, int cmglev, Any& crse, Any& fine) const;

    virtual void AnyInterpolationMG (int amrlev, int fmglev, Any& fine, const Any& crse) const;
    virtual void AnyInterpAssignMG (int amrlev, int fmglev, Any& fine, Any& crse) const;
    virtual void AnyInterpolationAmr (int famrlev, Any& fine, const Any& crse,
                                      IntVect const& /*nghost*/) const = 0;

    virtual void AnyAverageDownSolutionRHS (int camrlev, Any& crse_sol, Any& crse_rhs,
                                            const Any& fine_sol, const Any& fine_rhs);

    virtual void AnyAverageDownAndSync (Vector<Any>& sol) const = 0;

    Real MFNormInf (MultiFab const& mf, iMultiFab const* fine_mask, bool local) const;

    bool isMFIterSafe (int amrlev, int mglev1, int mglev2) const;

protected:

    static constexpr int mg_coarsen_ratio = 2;
    static constexpr int mg_box_min_width = 2;
#ifdef AMREX_USE_EB
    static constexpr int mg_domain_min_width = 4;
#else
    static constexpr int mg_domain_min_width = 2;
#endif

    LPInfo info;

    int verbose = 0;

    int maxorder = 3;

    bool enforceSingularSolvable = true;

    int m_num_amr_levels;
    Vector<int> m_amr_ref_ratio;

    Vector<int> m_num_mg_levels;
    const MLLinOp* m_parent = nullptr;

    IntVect m_ixtype;

    bool m_do_agglomeration = false;
    bool m_do_consolidation = false;

    bool m_do_semicoarsening = false;
    Vector<IntVect> mg_coarsen_ratio_vec;

    //! first Vector is for amr level and second is mg level
    Vector<Vector<Geometry> >            m_geom;
    Vector<Vector<BoxArray> >            m_grids;
    Vector<Vector<DistributionMapping> > m_dmap;
    Vector<Vector<std::unique_ptr<FabFactory<FArrayBox> > > > m_factory;
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
        ~CommContainer () {
#ifdef BL_USE_MPI
            if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
#endif
        }
    };
    std::unique_ptr<CommContainer> m_raii_comm;

    // BC
    Vector<Array<BCType, AMREX_SPACEDIM> > m_lobc;
    Vector<Array<BCType, AMREX_SPACEDIM> > m_hibc;
    // Need to save the original copy because we change the BC type to
    // Neumann for inhomogeneous Neumann and Robin.
    Vector<Array<BCType, AMREX_SPACEDIM> > m_lobc_orig;
    Vector<Array<BCType, AMREX_SPACEDIM> > m_hibc_orig;

    Array<Real, AMREX_SPACEDIM> m_domain_bloc_lo {{AMREX_D_DECL(0.,0.,0.)}};
    Array<Real, AMREX_SPACEDIM> m_domain_bloc_hi {{AMREX_D_DECL(0.,0.,0.)}};

    bool m_needs_coarse_data_for_bc;
    int m_coarse_data_crse_ratio = -1;
    RealVect m_coarse_bc_loc;
    const MultiFab* m_coarse_data_for_bc = nullptr;

    /**
    * \brief functions
    */
    int NAMRLevels () const noexcept { return m_num_amr_levels; }
    int NMGLevels (int amrlev) const noexcept { return m_num_mg_levels[amrlev]; }
    const Vector<int>& AMRRefRatio () const noexcept { return m_amr_ref_ratio; }
    int AMRRefRatio (int amr_lev) const noexcept { return m_amr_ref_ratio[amr_lev]; }

    const Geometry& Geom (int amr_lev, int mglev=0) const noexcept { return m_geom[amr_lev][mglev]; }
    FabFactory<FArrayBox> const* Factory (int amr_lev, int mglev=0) const noexcept {
        return m_factory[amr_lev][mglev].get();
    }

    GpuArray<BCType,AMREX_SPACEDIM> LoBC (int icomp = 0) const noexcept {
        return GpuArray<BCType,AMREX_SPACEDIM>{{AMREX_D_DECL(m_lobc[icomp][0],
                                                             m_lobc[icomp][1],
                                                             m_lobc[icomp][2])}};
    }
    GpuArray<BCType,AMREX_SPACEDIM> HiBC (int icomp = 0) const noexcept {
        return GpuArray<BCType,AMREX_SPACEDIM>{{AMREX_D_DECL(m_hibc[icomp][0],
                                                             m_hibc[icomp][1],
                                                             m_hibc[icomp][2])}};
    }

    bool hasInhomogNeumannBC () const noexcept;
    bool hasRobinBC () const noexcept;

    virtual bool supportRobinBC () const noexcept { return false; }
    virtual bool supportInhomogNeumannBC () const noexcept { return false; }

#ifdef BL_USE_MPI
    bool isBottomActive () const noexcept { return m_bottom_comm != MPI_COMM_NULL; }
#else
    bool isBottomActive () const noexcept { return true; }
#endif
    MPI_Comm BottomCommunicator () const noexcept { return m_bottom_comm; }
    MPI_Comm Communicator () const noexcept { return m_default_comm; }

    void setCoarseFineBCLocation (const RealVect& cloc) noexcept { m_coarse_bc_loc = cloc; }

    bool doAgglomeration () const noexcept { return m_do_agglomeration; }
    bool doConsolidation () const noexcept { return m_do_consolidation; }
    bool doSemicoarsening () const noexcept { return m_do_semicoarsening; }

    bool isCellCentered () const noexcept { return m_ixtype == 0; }

    void make (Vector<Vector<Any> >& mf, IntVect const& ng) const;

    virtual std::unique_ptr<FabFactory<FArrayBox> > makeFactory (int /*amrlev*/, int /*mglev*/) const {
        return std::make_unique<FArrayBoxFactory>();
    }

    virtual void resizeMultiGrid (int new_size);

    bool hasHiddenDimension () const noexcept { return info.hasHiddenDimension(); }
    int hiddenDirection () const noexcept { return info.hidden_direction; }
    Box compactify (Box const& b) const noexcept;

    template <typename T>
    Array4<T> compactify (Array4<T> const& a) const noexcept
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
    T get_d0 (T const& dx, T const& dy, T const&) const noexcept
    {
        if (info.hidden_direction == 0) {
            return dy;
        } else {
            return dx;
        }
    }

    template <typename T>
    T get_d1 (T const&, T const& dy, T const& dz) const noexcept
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
                      const Vector<FabFactory<FArrayBox> const*>& a_factory);
    void defineAuxData ();
    void defineBC ();
    static void makeAgglomeratedDMap (const Vector<BoxArray>& ba, Vector<DistributionMapping>& dm);
    static void makeConsolidatedDMap (const Vector<BoxArray>& ba, Vector<DistributionMapping>& dm,
                                      int ratio, int strategy);
    MPI_Comm makeSubCommunicator (const DistributionMapping& dm);
    void remapNeighborhoods (Vector<DistributionMapping> & dms);

    virtual void checkPoint (std::string const& /*file_name*/) const {
        amrex::Abort("MLLinOp:checkPoint: not implemented");
    }
};

}

#endif
