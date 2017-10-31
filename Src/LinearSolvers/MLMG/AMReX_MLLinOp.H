#ifndef AMREX_ML_LINOP_H_
#define AMREX_ML_LINOP_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_MacBndry.H>
#include <AMReX_YAFluxRegister.H>

namespace amrex {

class MLLinOp
{
public:

    friend class MLMG;
    friend class MLCGSolver;

    enum struct BCMode { Homogeneous, Inhomogeneous };
    enum struct BCType { Dirichlet, Neumann, Periodic, Bogus };

    MLLinOp () {}
    MLLinOp (const Vector<Geometry>& a_geom,
             const Vector<BoxArray>& a_grids,
             const Vector<DistributionMapping>& a_dmap);
    virtual ~MLLinOp ();

    MLLinOp (const MLLinOp&) = delete;
    MLLinOp (MLLinOp&&) = delete;
    MLLinOp& operator= (const MLLinOp&) = delete;
    MLLinOp& operator= (MLLinOp&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap);

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
    void setBCWithCoarseData (const MultiFab* crse, int crse_ratio);

    // Must be called for each level.  Argument `levelbcdata` is used
    // to supply Dirichlet bc at the physical domain.  However, even
    // if there is no Dirichlet bc, this funcion must still be called.
    // In that case, argument `levelbcdata` will be ignored and thus
    // could be nullptr.  The supplied MultiFab must have one ghost
    // cell.  Only the data outside the physical domain at Dirichlet
    // boundary will be used.  It is assumed that the data in those
    // ghost cells outside the domain live exactly on the face of the
    // phsyical domain.  Argumen `amrlev` is relative level such that
    // the lowest to the solver is always 0.
    void setLevelBC (int amrlev, const MultiFab* levelbcdata);

    void setMaxOrder (int o) { maxorder = o; }
    void setAgglomeration (int f) { do_agglomeration = f; }
    void setConsolidation (int f) { do_consolidation = f; }
    void setAgglomerationGridSize (int n) { agg_grid_size = n; }

protected:

    static constexpr int mg_coarsen_ratio = 2;
    static constexpr int mg_box_min_width = 2;

    int maxorder = 3;

    int do_agglomeration = 1;
    int do_consolidation = 1;
#if (AMREX_SPACEDIM == 1)
    int agg_grid_size = 32;
#elif (AMREX_SPACEDIM == 2)
    int agg_grid_size = 16;
#elif (AMREX_SPACEDIM == 3)
    int agg_grid_size = 8;
#endif

    int m_num_amr_levels;
    Vector<int> m_amr_ref_ratio;

    Vector<int> m_num_mg_levels;

    // first Vector is for amr level and second is mg level
    Vector<Vector<Geometry> >            m_geom;
    Vector<Vector<BoxArray> >            m_grids;
    Vector<Vector<DistributionMapping> > m_dmap;
    Vector<int>                          m_domain_covered;

    // BC
    std::array<BCType, AMREX_SPACEDIM> m_lobc {{AMREX_D_DECL(BCType::Bogus,BCType::Bogus,BCType::Bogus)}};
    std::array<BCType, AMREX_SPACEDIM> m_hibc {{AMREX_D_DECL(BCType::Bogus,BCType::Bogus,BCType::Bogus)}};

    bool m_needs_coarse_data_for_bc;
    int m_coarse_data_crse_ratio;
    const MultiFab* m_coarse_data_for_bc = nullptr;

    Vector<std::unique_ptr<MacBndry> >   m_bndry_sol;
    Vector<std::unique_ptr<BndryRegister> > m_crse_sol_br;

    Vector<std::unique_ptr<MacBndry> > m_bndry_cor;
    Vector<std::unique_ptr<BndryRegister> > m_crse_cor_br;

    // In case of agglomeration, coarse MG grids on amr level 0 are
    // not simply coarsened from fine MG grids.  So we need to build
    // bcond and bcloc for each MG level.
    using RealTuple = std::array<Real,2*BL_SPACEDIM>;
    using BCTuple   = std::array<BoundCond,2*BL_SPACEDIM>;
    class BndryCondLoc
    {
    public:
        BndryCondLoc (const BoxArray& ba, const DistributionMapping& dm);
        void setBndryConds (const Geometry& geom, const BCRec& phys_bc, int ratio);
        const BCTuple& bndryConds (const MFIter& mfi) const {
            return bcond[mfi];
        }
        const RealTuple& bndryLocs (const MFIter& mfi) const {
            return bcloc[mfi];
        }
    private:
        LayoutData<BCTuple>   bcond;
        LayoutData<RealTuple> bcloc;
    };
    Vector<Vector<std::unique_ptr<BndryCondLoc> > > m_bcondloc;

    // used to save interpolation coefficients of the first interior cells
    mutable Vector<Vector<BndryRegister> > m_undrrelxr;

    // boundary cell flags for covered, not_covered, outside_domain
    Vector<Vector<std::array<MultiMask,2*AMREX_SPACEDIM> > > m_maskvals;

    mutable Vector<YAFluxRegister> m_fluxreg;

    //
    // functions
    //

    int NAMRLevels () const { return m_num_amr_levels; }
    int NMGLevels (int amrlev) const { return m_num_mg_levels[amrlev]; }
    const Vector<int>& AMRRefRatio () const { return m_amr_ref_ratio; }
    int AMRRefRatio (int amr_lev) const { return m_amr_ref_ratio[amr_lev]; }

    const Geometry& Geom (int amr_lev, int mglev=0) const { return m_geom[amr_lev][mglev]; }

    void make (Vector<Vector<MultiFab> >& mf, int nc, int ng) const;

    void updateSolBC (int amrlev, const MultiFab& crse_bcdata) const;
    void updateCorBC (int amrlev, const MultiFab& crse_bcdata) const;

    void solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                           const MultiFab* crse_bcdata=nullptr);
    void correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
                             BCMode bc_mode, const MultiFab* crse_bcdata=nullptr);
    void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                const MacBndry* bndry=nullptr) const;
    void applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode,
                  const MacBndry* bndry=nullptr, bool skip_fillboundary=false) const;
    void smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs,
                 bool skip_fillboundary=false) const;
    // The assumption is crse_sol's boundary has been filled, but not fine_sol.
    void reflux (int crse_amrlev, MultiFab& res, const MultiFab& crse_sol, MultiFab& fine_sol) const;
    void compFlux (int amrlev, const std::array<MultiFab*,AMREX_SPACEDIM>& fluxes, MultiFab& sol) const;

    virtual void prepareForSolve ();
    virtual bool isSingular (int amrlev) const = 0;
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const = 0;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const = 0;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only=0) const = 0;
    virtual Real Anorm (int amrlev, int mglev) const = 0;

private:

    void defineGrids (const Vector<Geometry>& a_geom,
                      const Vector<BoxArray>& a_grids,
                      const Vector<DistributionMapping>& a_dmap);
    void defineAuxData ();
    void defineBC ();

};

}

#endif
