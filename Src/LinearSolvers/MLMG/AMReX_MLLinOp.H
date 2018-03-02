#ifndef AMREX_ML_LINOP_H_
#define AMREX_ML_LINOP_H_

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

namespace amrex {

struct LPInfo
{
    bool do_agglomeration = true;
    bool do_consolidation = true;
    int agg_grid_size = AMREX_D_PICK(32, 16, 8);
    int con_grid_size = AMREX_D_PICK(32, 16, 8);
    bool has_metric_term = true;

    LPInfo& setAgglomeration (bool x) { do_agglomeration = x; return *this; }
    LPInfo& setConsolidation (bool x) { do_consolidation = x; return *this; }
    LPInfo& setAgglomerationGridSize (int x) { agg_grid_size = x; return *this; }
    LPInfo& setConsolidationGridSize (int x) { con_grid_size = x; return *this; }
    LPInfo& setMetricTerm (bool x) { has_metric_term = x; return *this; }
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
                 const Vector<FabFactory<FArrayBox> const*>& a_factory);

    // Boundary of the whole domain. This functions must be called,
    // and must be called before other bc functions.
    void setDomainBC (const std::array<BCType,AMREX_SPACEDIM>& lobc,
                      const std::array<BCType,AMREX_SPACEDIM>& hibc);

    // Needs coarse data for bc?  If the lowest level grids does not
    // cover the entire domain, coarse level data are needed for
    // supplying Dirichlet bc at coarse/fine boundary, even when the
    // domain bc is not Dirichlet.
    bool needsCoarseDataForBC () const { return m_needs_coarse_data_for_bc; }

    // If needs coarse data for bc, this should be called.  MultiFab
    // crse does not need to have ghost cells.  The data are at coarse
    // resolution.  If this is called, it *MUST* be called before
    // `setLevelBC`.  If crse is nullptr, then bc value is assumed to
    // be zero.
    void setCoarseFineBC (const MultiFab* crse, int crse_ratio);

    // Must be called for each level.  Argument `levelbcdata` is used
    // to supply Dirichlet bc at the physical domain.  However, even
    // if there is no Dirichlet bc, this funcion must still be called.
    // In that case, argument `levelbcdata` will be ignored and thus
    // could be nullptr.  The supplied MultiFab must have one ghost
    // cell.  Only the data outside the physical domain at Dirichlet
    // boundary will be used.  It is assumed that the data in those
    // ghost cells outside the domain live exactly on the face of the
    // phsyical domain.  Argument `amrlev` is relative level such that
    // the lowest to the solver is always 0.
    virtual void setLevelBC (int amrlev, const MultiFab* levelbcdata) = 0;

    void setVerbose (int v) { verbose = v; }

    void setMaxOrder (int o) { maxorder = o; }

protected:

    static constexpr int mg_coarsen_ratio = 2;
    static constexpr int mg_box_min_width = 128;

    LPInfo info;

    int verbose = 0;

    int maxorder = 3;

    int m_num_amr_levels;
    Vector<int> m_amr_ref_ratio;

    Vector<int> m_num_mg_levels;
    const MLLinOp* m_parent = nullptr;

    IntVect m_ixtype;

    bool m_do_agglomeration = false;
    bool m_do_consolidation = false;

    // first Vector is for amr level and second is mg level
    Vector<Vector<Geometry> >            m_geom;
    Vector<Vector<BoxArray> >            m_grids;
    Vector<Vector<DistributionMapping> > m_dmap;
    Vector<Vector<std::unique_ptr<FabFactory<FArrayBox> > > > m_factory;
    Vector<int>                          m_domain_covered;

    MPI_Comm m_default_comm = MPI_COMM_NULL;
    MPI_Comm m_bottom_comm = MPI_COMM_NULL;
    struct CommContainer {
        MPI_Comm comm;
        CommContainer (MPI_Comm m) : comm(m) {}
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
    std::array<BCType, AMREX_SPACEDIM> m_lobc {{AMREX_D_DECL(BCType::bogus,BCType::bogus,BCType::bogus)}};
    std::array<BCType, AMREX_SPACEDIM> m_hibc {{AMREX_D_DECL(BCType::bogus,BCType::bogus,BCType::bogus)}};

    bool m_needs_coarse_data_for_bc;
    int m_coarse_data_crse_ratio = -1;
    RealVect m_coarse_bc_loc;
    const MultiFab* m_coarse_data_for_bc = nullptr;

    //
    // functions
    //

    int NAMRLevels () const { return m_num_amr_levels; }
    int NMGLevels (int amrlev) const { return m_num_mg_levels[amrlev]; }
    const Vector<int>& AMRRefRatio () const { return m_amr_ref_ratio; }
    int AMRRefRatio (int amr_lev) const { return m_amr_ref_ratio[amr_lev]; }

    const Geometry& Geom (int amr_lev, int mglev=0) const { return m_geom[amr_lev][mglev]; }

#ifdef BL_USE_MPI
    bool isBottomActive () const { return m_bottom_comm != MPI_COMM_NULL; }
#else
    bool isBottomActive () const { return true; }
#endif
    MPI_Comm BottomCommunicator () const { return m_bottom_comm; }
    MPI_Comm Communicator (int amrlev, int mglev) const {
        if (amrlev == 0 && mglev == NMGLevels(0)-1) {
            return m_bottom_comm;
        } else {
            return m_default_comm;
        }
    }

    void setCoarseFineBCLocation (const RealVect& cloc) { m_coarse_bc_loc = cloc; }

    bool doAgglomeration () const { return m_do_agglomeration; }
    bool doConsolidation () const { return m_do_consolidation; }

    bool isCellCentered () const { return m_ixtype == 0; }

    void make (Vector<Vector<MultiFab> >& mf, int nc, int ng) const;

    virtual void restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const = 0;
    virtual void interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const = 0;
    virtual void averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab& crse_rhs,
                                         const MultiFab& fine_sol, const MultiFab& fine_rhs) = 0;

    virtual void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                        const MLMGBndry* bndry=nullptr) const = 0;
    virtual void smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                         bool skip_fillboundary=false) const = 0;

    // Divide mf by the diagonal component of the operator. Used by bicgstab.
    virtual void normalize (int amrlev, int mglev, MultiFab& mf) const {}

    virtual void solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                   const MultiFab* crse_bcdata=nullptr) = 0;
    virtual void correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                                     BCMode bc_mode, const MultiFab* crse_bcdata=nullptr) = 0;

    virtual void reflux (int crse_amrlev,
                         MultiFab& res, const MultiFab& crse_sol, const MultiFab& crse_rhs,
                         MultiFab& fine_res, MultiFab& fine_sol, const MultiFab& fine_rhs) const = 0;
    virtual void compFlux (int amrlev, const std::array<MultiFab*,AMREX_SPACEDIM>& fluxes, MultiFab& sol) const = 0;
    virtual void compGrad (int amrlev, const std::array<MultiFab*,AMREX_SPACEDIM>& grad, MultiFab& sol) const = 0;
    
    virtual void applyMetricTerm (int amrlev, int mglev, MultiFab& rhs) const = 0;
    virtual void unapplyMetricTerm (int amrlev, int mglev, MultiFab& rhs) const = 0;
    virtual void fillSolutionBC (int amrlev, MultiFab& sol, const MultiFab* crse_bcdata=nullptr) = 0;

    virtual void prepareForSolve () = 0;
    virtual bool isSingular (int amrlev) const = 0;
    virtual bool isBottomSingular () const = 0;
    virtual Real xdoty (int amrlev, int mglev, const MultiFab& x, const MultiFab& y, bool local) const = 0;

    virtual Real getAScalar () const = 0;
    virtual Real getBScalar () const = 0;
    virtual MultiFab const* getACoeffs (int amrlev, int mglev) const = 0;
    virtual std::array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int amrlev, int mglev) const = 0;

    virtual void fixUpResidualMask (int amrlev, iMultiFab& resmsk) { }
    virtual void nodalSync (int amrlev, int mglev, MultiFab& mf) const {}

    virtual std::unique_ptr<MLLinOp> makeNLinOp (int grid_size) const = 0;

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
};

}

#endif
