#ifndef AMREX_ML_LINOP_H_
#define AMREX_ML_LINOP_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_MacBndry.H>

namespace amrex {

class MLLinOp
{
public:

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

    int NAMRLevels () const { return m_num_amr_levels; }
    int NMGLevels (int amrlev) const { return m_num_mg_levels[amrlev]; }
    const Vector<int>& AMRRefRatio () const { return m_amr_ref_ratio; }

    void make (Vector<Vector<MultiFab> >& mf, int nc, int ng) const;

    // TODO: need a flag to specify Dirichlet/Neumann, ...

    // The provided bcdata are expected to have ghost cells already filled.
    // Do we need crse_data's ghost cells filled?
    void setBC (int amrlev, const MultiFab* bcdata=nullptr, const MultiFab* crse_bcdata=nullptr);

    void updateBC (int amrlev, const MultiFab& crse_bcdata);
    
    virtual void prepareForSolve () = 0;
    virtual void residual (int amrlev, int mglev,
                           MultiFab& resid, MultiFab& sol, const MultiFab& rhs,
                           BCMode bc_mode); // const?
    virtual void apply (int amrlev, int mglev, MultiFab& out, MultiFab& in, BCMode bc_mode);
    virtual void applyBC (int amrlev, int mglev, MultiFab& in, BCMode bc_mode);
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) = 0;

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
    Vector<Vector<BndryRegister> >       m_undrrelxr;

    // boundary cell flags for covered, not_covered, outside_domain
    Vector<Vector<std::array<MultiMask,2*AMREX_SPACEDIM> > > m_maskvals;

    Vector<std::unique_ptr<MacBndry> >   m_macbndry;

};

}

#endif
