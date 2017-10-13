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

    enum struct BCMode { Homogeneous, Inhomogeneous };

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

    // bcdata is expected to have one ghost cell filled.
    // crse_bcdata, if !nullptr, doesn't need to have ghost cells.
    void setDirichletBC (int amrlev, const MultiFab& bcdata, const MultiFab* crse_bcdata=nullptr);

protected:

    static constexpr int mg_coarsen_ratio = 2;
    static constexpr int mg_box_min_width = 2;
    static constexpr int maxorder = 2;

    int m_num_amr_levels;
    Vector<int> m_amr_ref_ratio;

    Vector<int> m_num_mg_levels;

    // first Vector is for amr level and second is mg level
    Vector<Vector<Geometry> >            m_geom;
    Vector<Vector<BoxArray> >            m_grids;
    Vector<Vector<DistributionMapping> > m_dmap;

    // used to save interpolation coefficients of the first interior cells
    mutable Vector<Vector<BndryRegister> > m_undrrelxr;

    // boundary cell flags for covered, not_covered, outside_domain
    Vector<Vector<std::array<MultiMask,2*AMREX_SPACEDIM> > > m_maskvals;

    mutable Vector<YAFluxRegister> m_fluxreg;

    Vector<std::unique_ptr<MacBndry> >   m_bndry_sol;
    Vector<std::unique_ptr<BndryRegister> > m_crse_sol_br;

    Vector<std::unique_ptr<MacBndry> > m_bndry_cor;
    Vector<std::unique_ptr<BndryRegister> > m_crse_cor_br;

    //
    // functions
    //

    int NAMRLevels () const { return m_num_amr_levels; }
    int NMGLevels (int amrlev) const { return m_num_mg_levels[amrlev]; }
    const Vector<int>& AMRRefRatio () const { return m_amr_ref_ratio; }
    int AMRRefRatio (int amr_lev) const { return m_amr_ref_ratio[amr_lev]; }

    const Geometry& Geom (int amr_lev) const { return m_geom[amr_lev][0]; }

    void make (Vector<Vector<MultiFab> >& mf, int nc, int ng) const;

    // TODO: need a flag to specify Dirichlet/Neumann, ...

    void updateSolBC (int amrlev, const MultiFab& crse_bcdata);
    void updateCorBC (int amrlev, const MultiFab& crse_bcdata);

    void residual (int amrlev, int mglev,
                   MultiFab& resid, MultiFab& sol, const MultiFab& rhs,
                   BCMode bc_mode) const;
    void correctionResidual (int amrlev, MultiFab& resid, MultiFab& sol, const MultiFab& rhs) const;
    void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode,
                const MacBndry* bndry=nullptr) const;
    void applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode,
                  const MacBndry* bndry=nullptr) const;
    void smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, BCMode bc_mode) const;
    void reflux (int crse_amrlev, MultiFab& res, const MultiFab& crse_sol, const MultiFab& fine_sol) const;

    virtual void prepareForSolve () = 0;
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const = 0;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const = 0;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        std::array<FArrayBox,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, const int face_only=0) const = 0;

};

}

#endif
