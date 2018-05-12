#ifndef AMREX_EB2_ISLEVEL_H_
#define AMREX_EB2_ISLEVEL_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_EB2_Graph.H>

namespace amrex { namespace EB2 {

template <typename G>
class ISLevel
{
public:

    ISLevel (G const& gshop, const Geometry& geom, int max_grid_size);

private:

    Geometry m_geom;    
    BoxArray m_grids;
    DistributionMapping m_dmap;
    FabArray<Graph> m_graph;
};

template <typename G>
ISLevel<G>::ISLevel (G const& gshop, const Geometry& geom, int max_grid_size)
    : m_geom(geom)
{
    m_grids.define(geom.Domain());
    m_grids.maxSize(max_grid_size);
    m_dmap.define(m_grids);

    LayoutData<int> box_type(m_grids, m_dmap);
    for (MFIter mfi(box_type); mfi.isValid(); ++mfi)
    {
        const Box& bx = amrex::surroundingNodes(mfi.validbox());
        box_type[mfi] = gshop.getBoxType(bx, geom);
    }


    // we should simplify boxarray

//    FabArray<BaseFab<uint32_t> > cell_type(m_grids, m_dmap, 1, 2);
    m_graph.define(m_grids, m_dmap, 1, 2, MFInfo(), DefaultFabFactory<Graph>());
}


}}

#endif
