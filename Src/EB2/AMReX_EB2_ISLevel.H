#ifndef AMREX_EB2_ISLEVEL_H_
#define AMREX_EB2_ISLEVEL_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_VisMF.H>
#include <AMReX_EB2_MultiGFab.H>

#include <unordered_map>

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
    MultiFab m_volfrac;
    MultiFab m_centroid;
    MultiFab m_bndrycent;
    MultiFab m_bndrynorm;
    Array<MultiFab,AMREX_SPACEDIM> m_areafrac;
    Array<MultiFab,AMREX_SPACEDIM> m_facecent;
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
    m_volfrac.define(m_grids, m_dmap, 1, 0);
    m_centroid.define(m_grids, m_dmap, AMREX_SPACEDIM, 0);
    m_bndrycent.define(m_grids, m_dmap, AMREX_SPACEDIM, 0);
    m_bndrynorm.define(m_grids, m_dmap, AMREX_SPACEDIM, 0);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_areafrac[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)),
                                m_dmap, 1, 0);
        m_facecent[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)),
                                m_dmap, AMREX_SPACEDIM, 0);
    }

    for (MFIter mfi(m_mgf); mfi.isValid(); ++mfi)
    {
        auto& gfab = m_mgf[mfi];
        const Box& vbx = gfab.validbox();

        auto& levelset = gfab.getLevelSet();
        gshop.fillFab(levelset, geom);

        gfab.buildTypes();

//        auto& facetype = gfab.getFaceType();
        
#if 0
        Array<std::unordered_map<IntVect,Array<Real,AMREX_SPACEDIM-1>,IntVect::hasher>,
              AMREX_SPACEDIM> irreg_faces;

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            auto& irreg_face_map = irreg_faces[idim];
            const auto& facetype_fab = facetype[idim];
            const Box& b = facetype_fab.box();
            facetype_fab.ForEachIV(b, 0, 1,
                [&irreg_face_map, &gshop, &geom, &idim] (const Type_t& t, const IntVect& iv) {
                    if (t == Type::irregular) {
                        const auto& intercept = gshop.getFaceIntercept(idim, iv, geom);
//                        irreg_face_map[iv] = 
                    }
                });
        }

        auto& celltype_fab = gfab.getCellType();
        auto& graph_fab = gfab.getGraph();
        celltype_fab.ForEachIV(vbx, 0, 1,
            [] (const Type_t& t, const IntVect& iv) {
                if (t == Type::irregular) {
                    amrex::Print() << "irreg cell " << iv << "\n";
                    
                }
            });

        // compute volume stuff in irregular cells
        // compute 
        
        // add cutcells in +1 region to the graph and note that in ghost region
        // information is incomplete.
        
        // 
#endif
    }

    VisMF::Write(m_mgf.getLevelSet(), "levelset");
}


}}

#endif
