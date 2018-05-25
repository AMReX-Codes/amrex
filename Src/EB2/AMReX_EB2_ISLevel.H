#ifndef AMREX_EB2_ISLEVEL_H_
#define AMREX_EB2_ISLEVEL_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_VisMF.H>
#include <AMReX_EB2_MultiGFab.H>
#include <AMReX_EB2_F.H>

#include <unordered_map>
#include <limits>

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
    m_volfrac.define(m_grids, m_dmap, 1, 1);
    m_centroid.define(m_grids, m_dmap, AMREX_SPACEDIM, 1);
    m_bndrycent.define(m_grids, m_dmap, AMREX_SPACEDIM, 1);
    m_bndrynorm.define(m_grids, m_dmap, AMREX_SPACEDIM, 1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_areafrac[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)),
                                m_dmap, 1, 1);
        m_facecent[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)),
                                m_dmap, AMREX_SPACEDIM, 1);
    }

    Array<BaseFab<Real>, AMREX_SPACEDIM> intercept;
    const Real* dx = geom.CellSize();
    const Real* dxinv = geom.InvCellSize();
    const Real* problo = geom.ProbLo();

    for (MFIter mfi(m_mgf); mfi.isValid(); ++mfi)
    {
        auto& gfab = m_mgf[mfi];
        const Box& vbx = gfab.validbox();

        auto& levelset = gfab.getLevelSet();
        gshop.fillFab(levelset, geom);

        gfab.buildTypes();

#if (AMREX_SPACEDIM == 3)
        auto& celltype = gfab.getCellType();
        auto& facetype = gfab.getFaceType();
        auto& edgetype = gfab.getEdgeType();
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const auto& edgetype_fab = edgetype[idim];
            const Box& b = edgetype_fab.box();
            auto& inter = intercept[idim];
            inter.resize(b);
            inter.setVal(std::numeric_limits<Real>::quiet_NaN());
            edgetype_fab.ForEachIV(b, 0, 1,
                [&inter, &gshop, &geom, &idim] (const Type_t& t, const IntVect& iv) {
                    if (t == Type::irregular) {
                        inter(iv) = gshop.getIntercept(idim, iv, geom);
                    }
                });
        }

        amrex_eb2_build_faces(BL_TO_FORTRAN_BOX(vbx),
                              BL_TO_FORTRAN_ANYD(celltype),
                              BL_TO_FORTRAN_ANYD(facetype[0]),
                              BL_TO_FORTRAN_ANYD(facetype[1]),
                              BL_TO_FORTRAN_ANYD(facetype[2]),
                              BL_TO_FORTRAN_ANYD(edgetype[0]),
                              BL_TO_FORTRAN_ANYD(edgetype[1]),
                              BL_TO_FORTRAN_ANYD(edgetype[2]),
                              BL_TO_FORTRAN_ANYD(levelset),
                              BL_TO_FORTRAN_ANYD(intercept[0]),
                              BL_TO_FORTRAN_ANYD(intercept[1]),
                              BL_TO_FORTRAN_ANYD(intercept[2]),
                              BL_TO_FORTRAN_ANYD(m_areafrac[0][mfi]),
                              BL_TO_FORTRAN_ANYD(m_areafrac[1][mfi]),
                              BL_TO_FORTRAN_ANYD(m_areafrac[2][mfi]),
                              BL_TO_FORTRAN_ANYD(m_facecent[0][mfi]),
                              BL_TO_FORTRAN_ANYD(m_facecent[1][mfi]),
                              BL_TO_FORTRAN_ANYD(m_facecent[2][mfi]),
                              dx, dxinv, problo);
                              
#endif

    }

    VisMF::Write(m_areafrac[0], "apx");
    VisMF::Write(m_facecent[0], "fcx");
    VisMF::Write(m_mgf.getLevelSet(), "levelset");
}


}}

#endif
