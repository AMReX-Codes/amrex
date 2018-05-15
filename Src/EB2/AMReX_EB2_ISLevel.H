#ifndef AMREX_EB2_ISLEVEL_H_
#define AMREX_EB2_ISLEVEL_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_VisMF.H>
#include <AMReX_EB2_MultiGFab.H>

namespace amrex { namespace EB2 {

template <typename G>
class ISLevel
{
public:

    ISLevel (G const& gshop, const Geometry& geom, int max_grid_size);

private:

    Geometry m_geom;    
    BoxArray m_grids;
    BoxArray m_covered_grids;
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
    Vector<Box> covered_boxes;

    for (MFIter mfi(m_grids, m_dmap); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const Box& gbx = amrex::surroundingNodes(amrex::grow(vbx,1) & pdom);
        int box_type = gshop.getBoxType(gbx, geom);
        if (box_type == gshop.allcovered) {
            covered_boxes.push_back(vbx);
        } else if (box_type == gshop.mixedcells) {
            cut_boxes.push_back(vbx);
        }
    }

    amrex::AllGatherBoxes(cut_boxes);
    amrex::AllGatherBoxes(covered_boxes);

    if (!covered_boxes.empty()) {
        m_covered_grids = BoxArray(BoxList(std::move(covered_boxes)));
    }

    if (cut_boxes.empty()) {
        m_grids = BoxArray();
        m_dmap = DistributionMapping();
        return;
    }

    m_grids = BoxArray(BoxList(std::move(cut_boxes)));
    m_dmap = DistributionMapping(m_grids);

    m_mgf.define(m_grids, m_dmap);

    using faces_t = Array<Vector<std::pair<IntVect,
                                           Array<Real,AMREX_SPACEDIM-1> > >,
                          AMREX_SPACEDIM>;
    faces_t irreg_faces;
    for (MFIter mfi(m_mgf); mfi.isValid(); ++mfi)
    {
        auto& gfab = m_mgf[mfi];

        auto& levelset = gfab.getLevelSet();
        gshop.fillFab(levelset, geom);

        gfab.buildTypes();

        auto& facetype = gfab.getFaceType();
        
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            auto& irreg_face_vec = irreg_faces[idim];
            irreg_face_vec.clear();
            const auto& facetype_fab = facetype[idim];
            const Box& b = facetype_fab.box();
            facetype_fab.ForEachIV(b, 0, 1,
                [&irreg_face_vec, &gshop, &geom, &idim] (const Type_t& t, const IntVect& iv) {
                       if (t == Type::irregular) {
                           Array<Real,AMREX_SPACEDIM-1> intercept = gshop.getFaceIntercept(idim, iv, geom);
                           irreg_face_vec.emplace_back(std::make_pair(iv,intercept));
                       }
                });
        }


        // compute volume stuff in irregular cells
        // compute 
        
        // add cutcells in +1 region to the graph and note that in ghost region
        // information is incomplete.
        
        // 
    }

    VisMF::Write(m_mgf.getLevelSet(), "levelset");
}


}}

#endif
