#ifndef AMREX_EB2_MULTIGFAB_H_
#define AMREX_EB2_MULTIGFAB_H_

#include <AMReX_EB2_Graph.H>
#include <AMReX_Array.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MultiFab.H>

namespace amrex { namespace EB2 {

class GFab
{
public:

    GFab () = default;

    static constexpr int ngrow_levelset = 2;
    static constexpr int ngrow_celltype = 2;
    static constexpr int ngrow_facetype = 1;
    static constexpr int ngrow_edgetype = 1;
    
    void define(const Box& vbx) {
        m_validbox = vbx;
        const Box& bxg1 = amrex::grow(vbx,1);
        const Box& bxg2 = amrex::grow(vbx,2);
        m_levelset.resize(amrex::surroundingNodes(bxg2));
        m_celltype.resize(bxg2,1);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            m_facetype[idim].resize(amrex::convert(bxg1,IntVect::TheDimensionVector(idim)));
            m_edgetype[idim].resize(amrex::convert(bxg1,IntVect::TheNodeVector()
                                                       -IntVect::TheDimensionVector(idim)));
        }
    }

    Box validbox () const { return m_validbox; }

    const BaseFab<Real>& getLevelSet () const { return m_levelset; }
    BaseFab<Real>& getLevelSet () { return m_levelset; }

    const BaseFab<Type_t>& getCellType () const { return m_celltype; }
    BaseFab<Type_t>& getCellType () { return m_celltype; }

    const Array<BaseFab<Type_t>,AMREX_SPACEDIM>& getFaceType () const { return m_facetype; }
    Array<BaseFab<Type_t>,AMREX_SPACEDIM>& getFaceType () { return m_facetype; }

    const Array<BaseFab<Type_t>,AMREX_SPACEDIM>& getEdgeType () const { return m_edgetype; }
    Array<BaseFab<Type_t>,AMREX_SPACEDIM>& getEdgeType () { return m_edgetype; }

    const Graph& getGraph () const { return m_graph; }
    Graph& getGraph () { return m_graph; }

    void buildTypes ();

private:

    Box m_validbox;
    BaseFab<Real> m_levelset;
    BaseFab<Type_t> m_celltype;
    Array<BaseFab<Type_t>, AMREX_SPACEDIM> m_facetype;
    Array<BaseFab<Type_t>, AMREX_SPACEDIM> m_edgetype;
    Graph m_graph;
};

class MultiGFab
    : public LayoutData<GFab>
{
public:

    MultiGFab () = default;

    MultiGFab (const BoxArray& ba, const DistributionMapping& dm) {
        define(ba,dm);
    }

    void define (const BoxArray& ba, const DistributionMapping& dm) {
        LayoutData<GFab>::define(ba, dm);
        for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto& gfab = (*this)[mfi];
            gfab.define(bx);
        }
    }

    MultiFab getLevelSet ();
};

}}

#endif
