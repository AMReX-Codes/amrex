#ifndef AMREX_EB2_MULTIGFAB_H_
#define AMREX_EB2_MULTIGFAB_H_

#include <AMReX_EB2_Graph.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_Array.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MultiFab.H>

namespace amrex { namespace EB2 {

class GFab
{
public:

    GFab () = default;

    static constexpr int ng = 2;

    void define(const Box& vbx) {
        m_validbox = vbx;
        const Box& bxg2 = amrex::grow(vbx,ng);
        m_levelset.resize(amrex::surroundingNodes(bxg2));
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            m_facetype[idim].resize(amrex::convert(bxg2,IntVect::TheDimensionVector(idim)));
#if (AMREX_SPACEDIM == 3)
            m_edgetype[idim].resize(amrex::convert(bxg2,IntVect::TheNodeVector()
                                                       -IntVect::TheDimensionVector(idim)));
#endif
        }
    }

    Box validbox () const { return m_validbox; }

    const BaseFab<Real>& getLevelSet () const { return m_levelset; }
    BaseFab<Real>& getLevelSet () { return m_levelset; }

    const Array<BaseFab<Type_t>,AMREX_SPACEDIM>& getFaceType () const { return m_facetype; }
    Array<BaseFab<Type_t>,AMREX_SPACEDIM>& getFaceType () { return m_facetype; }

#if (AMREX_SPACEDIM == 3)
    const Array<BaseFab<Type_t>,AMREX_SPACEDIM>& getEdgeType () const { return m_edgetype; }
    Array<BaseFab<Type_t>,AMREX_SPACEDIM>& getEdgeType () { return m_edgetype; }
#endif

    const Graph& getGraph () const { return m_graph; }
    Graph& getGraph () { return m_graph; }

    void buildTypes (EBCellFlagFab& celltype);

private:

    Box m_validbox;
    BaseFab<Real> m_levelset;
    Array<BaseFab<Type_t>, AMREX_SPACEDIM> m_facetype;
#if (AMREX_SPACEDIM == 3)
    Array<BaseFab<Type_t>, AMREX_SPACEDIM> m_edgetype;
#endif
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
