#ifndef AMREX_EB2_ISLEVEL_H_
#define AMREX_EB2_ISLEVEL_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_EB2_MultiGFab.H>

namespace amrex { namespace EB2 {

template <typename G>
class ISLevel
{
public:

    ISLevel (G const& gshop, const Geometry& geom, int max_grid_size);

private:

    bool m_allregular = false;
    Geometry m_geom;    
    BoxArray m_grids;
    DistributionMapping m_dmap;
    MultiGFab m_mgf;
};

template <typename G>
ISLevel<G>::ISLevel (G const& gshop, const Geometry& geom, int max_grid_size)
    : m_geom(geom)
{
    m_grids.define(geom.Domain());
    m_grids.maxSize(max_grid_size);
    m_dmap.define(m_grids);

    Box pdom = geom.Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            pdom.grow(idim,2);
        }
    }

    Vector<Box> cut_boxes;

    for (MFIter mfi(m_grids, m_dmap); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const Box& gbx = amrex::surroundingNodes(amrex::grow(vbx,1) & pdom);
        int box_type = gshop.getBoxType(gbx, geom);
        if (box_type > 0) { // contains irregular cells
            cut_boxes.push_back(vbx);
        }
    }

    amrex::AllGatherBoxes(cut_boxes);

    if (cut_boxes.empty()) {
        m_allregular = true;
        m_grids = BoxArray();
        m_dmap = DistributionMapping();
        return;
    }

    m_grids = BoxArray(BoxList(std::move(cut_boxes)));
    m_dmap = DistributionMapping(m_grids);

    m_mgf.define(m_grids, m_dmap);
    
    BaseFab<int> inout;
    for (MFIter mfi(m_mgf); mfi.isValid(); ++mfi)
    {
        auto& gfab = m_mgf[mfi];

        auto& levelset = gfab.getLevelSet();
        gshop.fillFab(levelset, geom);

        gfab.buildTypes();

        // set cell_type on +1 region
        // set face_type
        // need to fix regular right next to covered cells. on +1 region
        
        // compute volume stuff in irregular cells
        // compute 
        
        // add cutcells in +1 region to the graph and note that in ghost region
        // information is incomplete.
        
        // 
    }

}


}}

#endif
