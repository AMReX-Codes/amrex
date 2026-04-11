#ifndef AMREX_EB2_MULTIGFAB_H_
#define AMREX_EB2_MULTIGFAB_H_
#include <AMReX_Config.H>

#include <AMReX_EB2_Graph.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_Array.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MultiFab.H>

/**
 * \file AMReX_EB2_MultiGFab.H
 *
 * Buffers backing EB graph construction (level sets, face/edge flags, etc.).
 */

namespace amrex::EB2 {

//! Object that stores level-set, face/edge types, and EB graph per box.
class GFab
{
public:

    GFab () = default;

    static constexpr int ng = 2;

    /**
     * \brief Allocate level-set and type arrays for valid box \p vbx.
     */
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

    //! Valid box covered by this GFab.
    Box validbox () const { return m_validbox; }

    //! Access the nodal level-set field (const and mutable).
    const BaseFab<Real>& getLevelSet () const { return m_levelset; }
    BaseFab<Real>& getLevelSet () { return m_levelset; }

    //! Access face-type classification arrays.
    const Array<BaseFab<Type_t>,AMREX_SPACEDIM>& getFaceType () const { return m_facetype; }
    Array<BaseFab<Type_t>,AMREX_SPACEDIM>& getFaceType () { return m_facetype; }

#if (AMREX_SPACEDIM == 3)
    //! Access edge-type arrays (3-D only).
    const Array<BaseFab<Type_t>,AMREX_SPACEDIM>& getEdgeType () const { return m_edgetype; }
    Array<BaseFab<Type_t>,AMREX_SPACEDIM>& getEdgeType () { return m_edgetype; }
#endif

    //! Accessor for the EB graph structure stored within this tile.
    const Graph& getGraph () const { return m_graph; }
    Graph& getGraph () { return m_graph; }

    /**
     * \brief Derive face/edge types from the supplied EB cell flags.
     *
     * \param celltype EBCellFlag data covering the grown region.
     */
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

//! LayoutData wrapper that allocates a GFab per FAB/Box.
class MultiGFab
    : public LayoutData<GFab>
{
public:

    MultiGFab () = default;

    //! Construct and define the layout immediately.
    MultiGFab (const BoxArray& ba, const DistributionMapping& dm) {
        define(ba,dm);
    }

    /**
     * \brief Define the layout and allocate per-box GFab storage.
     */
    void define (const BoxArray& ba, const DistributionMapping& dm) {
        LayoutData<GFab>::define(ba, dm);
        for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto& gfab = (*this)[mfi];
            gfab.define(bx);
        }
    }

    //! Collect the level-set data from each GFab into a MultiFab.
    MultiFab getLevelSet ();
};

}

#endif
