#ifndef AMREX_ML_LINOP_H_
#define AMREX_ML_LINOP_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

namespace amrex {

class MLLinOp
{
public:

    MLLinOp () {}
    MLLinOp (const Vector<Geometry>& a_geom,
             const Vector<BoxArray>& a_grids,
             const Vector<DistributionMapping>& a_dmap);
    ~MLLinOp ();

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

    virtual void prepareForSolve () = 0;

protected:

    static constexpr int mg_coarsen_ratio = 2;
    static constexpr int mg_box_min_width = 2;

    int m_num_amr_levels;
    Vector<int> m_amr_ref_ratio;

    Vector<int> m_num_mg_levels;

    // first Vector is for amr level and second is mg level
    Vector<Vector<Geometry> >            m_geom;
    Vector<Vector<BoxArray> >            m_grids;
    Vector<Vector<DistributionMapping> > m_dmap;

};

}

#endif
