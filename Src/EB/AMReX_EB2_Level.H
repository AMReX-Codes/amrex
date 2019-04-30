#ifndef AMREX_EB2_LEVEL_H_
#define AMREX_EB2_LEVEL_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_VisMF.H>
#include <AMReX_Array.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EB2_MultiGFab.H>
#include <AMReX_EB2_F.H>
#include <AMReX_EB2_IF_AllRegular.H>

#include <unordered_map>
#include <limits>
#include <cmath>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex { namespace EB2 {

class IndexSpace;

class Level
{
public:

    bool isAllRegular () const noexcept { return m_allregular; }
    bool isOK () const noexcept { return m_ok; }
    void fillEBCellFlag (FabArray<EBCellFlagFab>& cellflag, const Geometry& geom) const;
    void fillVolFrac (MultiFab& vfrac, const Geometry& geom) const;
    void fillCentroid (MultiCutFab& centroid, const Geometry& geom) const;
    void fillCentroid (   MultiFab& centroid, const Geometry& geom) const;
    void fillBndryArea (MultiCutFab& bndryarea, const Geometry& geom) const;
    void fillBndryArea (   MultiFab& bndryarea, const Geometry& geom) const;
    void fillBndryCent (MultiCutFab& bndrycent, const Geometry& geom) const;
    void fillBndryCent (   MultiFab& bndrycent, const Geometry& geom) const;
    void fillBndryNorm (MultiCutFab& bndrynorm, const Geometry& geom) const;
    void fillBndryNorm (   MultiFab& bndrynorm, const Geometry& geom) const;
    void fillAreaFrac (Array<MultiCutFab*,AMREX_SPACEDIM> const& areafrac, const Geometry& geom) const;
    void fillAreaFrac (Array<   MultiFab*,AMREX_SPACEDIM> const& areafrac, const Geometry& geom) const;
    void fillFaceCent (Array<MultiCutFab*,AMREX_SPACEDIM> const& facefrac, const Geometry& geom) const;
    void fillFaceCent (Array<   MultiFab*,AMREX_SPACEDIM> const& facefrac, const Geometry& geom) const;
    void fillLevelSet (MultiFab& levelset, const Geometry& geom) const;

    const BoxArray& boxArray () const noexcept { return m_grids; }
    const DistributionMapping& DistributionMap () const noexcept { return m_dmap; }

    Level (IndexSpace const* is, const Geometry& geom) : m_geom(geom), m_parent(is) {}
    void prepareForCoarsening (const Level& rhs, int max_grid_size, IntVect ngrow);

    const Geometry& Geom () const noexcept { return m_geom; }
    IndexSpace const* getEBIndexSpace () const noexcept { return m_parent; }

protected:

    Level (Level && rhs) = default;

    Level (Level const& rhs) = delete;
    Level& operator=(Level const& rhs) = delete;
    Level& operator=(Level && rhs) = delete;

    int coarsenFromFine (Level& fineLevel, bool fill_boundary);
    void buildCellFlag ();

    Geometry m_geom;
    IntVect  m_ngrow;
    BoxArray m_grids;
    BoxArray m_covered_grids;
    DistributionMapping m_dmap;
    MultiGFab m_mgf;
    MultiFab m_levelset;
    FabArray<EBCellFlagFab> m_cellflag;
    MultiFab m_volfrac;
    MultiFab m_centroid;
    MultiFab m_bndryarea;
    MultiFab m_bndrycent;
    MultiFab m_bndrynorm;
    Array<MultiFab,AMREX_SPACEDIM> m_areafrac;
    Array<MultiFab,AMREX_SPACEDIM> m_facecent;
    bool m_allregular = false;
    bool m_ok = false;
    IndexSpace const* m_parent;
};

template <typename G>
class GShopLevel
    : public Level
{
public:
    GShopLevel (IndexSpace const* is, G const& gshop, const Geometry& geom, int max_grid_size, int ngrow);
    GShopLevel (IndexSpace const* is, int ilev, int max_grid_size, int ngrow,
                const Geometry& geom, GShopLevel<G>& fineLevel);
};

template <typename G>
GShopLevel<G>::GShopLevel (IndexSpace const* is, G const& gshop, const Geometry& geom,
                           int max_grid_size, int ngrow)
    : Level(is, geom)
{
    if (std::is_same<typename G::FunctionType, AllRegularIF>::value) {
        m_allregular = true;
        m_ok = true;
        return;
    }

    BL_PROFILE("EB2::GShopLevel()-fine");

    // make sure ngrow is multiple of 16
    m_ngrow = IntVect{static_cast<int>(std::ceil(ngrow/16.)) * 16};

    Box domain_grown = geom.Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (Geometry::isPeriodic(idim)) {
            m_ngrow[idim] = 0;
        } else {
            m_ngrow[idim] = std::min(m_ngrow[idim], domain_grown.length(idim));
        }
    }
    domain_grown.grow(m_ngrow);

    m_grids.define(domain_grown);
    m_grids.maxSize(max_grid_size);
    m_dmap.define(m_grids);

    Vector<Box> cut_boxes;
    Vector<Box> covered_boxes;

    for (MFIter mfi(m_grids, m_dmap); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        const Box& gbx = amrex::surroundingNodes(amrex::grow(vbx,1));
        int box_type = gshop.getBoxType(gbx, geom);
        if (box_type == gshop.allcovered) {
            covered_boxes.push_back(vbx);
        } else if (box_type == gshop.mixedcells) {
            cut_boxes.push_back(vbx);
        }
    }

    amrex::AllGatherBoxes(cut_boxes);
    amrex::AllGatherBoxes(covered_boxes);

    if ( cut_boxes.empty() && 
        !covered_boxes.empty()) 
    {
        amrex::Abort("AMReX_EB2_Level.H: Domain is completely covered");
    }

    if (!covered_boxes.empty()) {
        m_covered_grids = BoxArray(BoxList(std::move(covered_boxes)));
    }

    if (cut_boxes.empty()) {
        m_grids = BoxArray();
        m_dmap = DistributionMapping();
        m_allregular = true;
        m_ok = true;
        return;
    }

    m_grids = BoxArray(BoxList(std::move(cut_boxes)));
    m_dmap = DistributionMapping(m_grids);

    m_mgf.define(m_grids, m_dmap);
    const int ng = 2;
    m_cellflag.define(m_grids, m_dmap, 1, ng);
    m_volfrac.define(m_grids, m_dmap, 1, ng);
    m_centroid.define(m_grids, m_dmap, AMREX_SPACEDIM, ng);
    m_bndryarea.define(m_grids, m_dmap, 1, ng);
    m_bndrycent.define(m_grids, m_dmap, AMREX_SPACEDIM, ng);
    m_bndrynorm.define(m_grids, m_dmap, AMREX_SPACEDIM, ng);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_areafrac[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)),
                                m_dmap, 1, ng);
        m_facecent[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)),
                                m_dmap, AMREX_SPACEDIM-1, ng);
    }

    const Real* dx = geom.CellSize();
    const Real* dxinv = geom.InvCellSize();
    const Real* problo = geom.ProbLo();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<BaseFab<Real>, AMREX_SPACEDIM> intercept;
        Array<BaseFab<Real>, AMREX_SPACEDIM> M2;

        for (MFIter mfi(m_mgf); mfi.isValid(); ++mfi)
        {
            auto& gfab = m_mgf[mfi];
            const Box& vbx = gfab.validbox();

            auto& levelset = gfab.getLevelSet();
            gshop.fillFab(levelset, geom);

            auto& cellflag = m_cellflag[mfi];

            gfab.buildTypes(cellflag);

#if (AMREX_SPACEDIM == 3)
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

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const Box& b = facetype[idim].box();
                M2[idim].resize(b,3);
            }

            // regular by default
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_areafrac[idim][mfi].setVal(1.0);
                m_facecent[idim][mfi].setVal(0.0);
            }

            amrex_eb2_build_faces(BL_TO_FORTRAN_BOX(vbx),
                                  BL_TO_FORTRAN_ANYD(cellflag),
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
                                  BL_TO_FORTRAN_ANYD(M2[0]),
                                  BL_TO_FORTRAN_ANYD(M2[1]),
                                  BL_TO_FORTRAN_ANYD(M2[2]),
                                  dx, dxinv, problo);

            // regular by default
            m_volfrac[mfi].setVal(1.0);
            m_centroid[mfi].setVal(0.0);
            m_bndryarea[mfi].setVal(0.0);
            m_bndrycent[mfi].setVal(-1.0);
            m_bndrynorm[mfi].setVal(0.0);

            amrex_eb2_build_cells(BL_TO_FORTRAN_BOX(vbx),
                                  BL_TO_FORTRAN_ANYD(cellflag),
                                  BL_TO_FORTRAN_ANYD(facetype[0]),
                                  BL_TO_FORTRAN_ANYD(facetype[1]),
                                  BL_TO_FORTRAN_ANYD(facetype[2]),
                                  BL_TO_FORTRAN_ANYD(m_areafrac[0][mfi]),
                                  BL_TO_FORTRAN_ANYD(m_areafrac[1][mfi]),
                                  BL_TO_FORTRAN_ANYD(m_areafrac[2][mfi]),
                                  BL_TO_FORTRAN_ANYD(m_facecent[0][mfi]),
                                  BL_TO_FORTRAN_ANYD(m_facecent[1][mfi]),
                                  BL_TO_FORTRAN_ANYD(m_facecent[2][mfi]),
                                  BL_TO_FORTRAN_ANYD(M2[0]),
                                  BL_TO_FORTRAN_ANYD(M2[1]),
                                  BL_TO_FORTRAN_ANYD(M2[2]),
                                  BL_TO_FORTRAN_ANYD(m_volfrac[mfi]),
                                  BL_TO_FORTRAN_ANYD(m_centroid[mfi]),
                                  BL_TO_FORTRAN_ANYD(m_bndryarea[mfi]),
                                  BL_TO_FORTRAN_ANYD(m_bndrycent[mfi]),
                                  BL_TO_FORTRAN_ANYD(m_bndrynorm[mfi]));

#elif (AMREX_SPACEDIM == 2)
            auto& facetype = gfab.getFaceType();
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const auto& facetype_fab = facetype[idim];
                const Box& b = facetype_fab.box();
                auto& inter = intercept[idim];
                inter.resize(b);
                inter.setVal(std::numeric_limits<Real>::quiet_NaN());
                facetype_fab.ForEachIV(b, 0, 1,
                    [&inter, &gshop, &geom, &idim] (const Type_t& t, const IntVect& iv) {
                        if (t == Type::irregular) {
                            inter(iv) = gshop.getIntercept(1-idim, iv, geom);
                        }
                    });
            }

            // regular by default
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                m_areafrac[idim][mfi].setVal(1.0);
                m_facecent[idim][mfi].setVal(0.0);
            }

            amrex_eb2_build_faces(BL_TO_FORTRAN_BOX(vbx),
                                  BL_TO_FORTRAN_ANYD(cellflag),
                                  BL_TO_FORTRAN_ANYD(facetype[0]),
                                  BL_TO_FORTRAN_ANYD(facetype[1]),
                                  BL_TO_FORTRAN_ANYD(levelset),
                                  BL_TO_FORTRAN_ANYD(intercept[0]),
                                  BL_TO_FORTRAN_ANYD(intercept[1]),
                                  BL_TO_FORTRAN_ANYD(m_areafrac[0][mfi]),
                                  BL_TO_FORTRAN_ANYD(m_areafrac[1][mfi]),
                                  BL_TO_FORTRAN_ANYD(m_facecent[0][mfi]),
                                  BL_TO_FORTRAN_ANYD(m_facecent[1][mfi]),
                                  dx, dxinv, problo);

            // regular by default
            m_volfrac[mfi].setVal(1.0);
            m_centroid[mfi].setVal(0.0);
            m_bndryarea[mfi].setVal(0.0);
            m_bndrycent[mfi].setVal(-1.0);
            m_bndrynorm[mfi].setVal(0.0);

            amrex_eb2_build_cells(BL_TO_FORTRAN_BOX(vbx),
                                  BL_TO_FORTRAN_ANYD(cellflag),
                                  BL_TO_FORTRAN_ANYD(facetype[0]),
                                  BL_TO_FORTRAN_ANYD(facetype[1]),
                                  BL_TO_FORTRAN_ANYD(m_areafrac[0][mfi]),
                                  BL_TO_FORTRAN_ANYD(m_areafrac[1][mfi]),
                                  BL_TO_FORTRAN_ANYD(m_volfrac[mfi]),
                                  BL_TO_FORTRAN_ANYD(m_centroid[mfi]),
                                  BL_TO_FORTRAN_ANYD(m_bndryarea[mfi]),
                                  BL_TO_FORTRAN_ANYD(m_bndrycent[mfi]),
                                  BL_TO_FORTRAN_ANYD(m_bndrynorm[mfi]));
#endif
        }
    }

    m_levelset = m_mgf.getLevelSet();

    m_ok = true;
}


template <typename G>
GShopLevel<G>::GShopLevel (IndexSpace const* is, int ilev, int max_grid_size, int /*ngrow*/,
                           const Geometry& geom, GShopLevel<G>& fineLevel)
    : Level(is, geom)
{
    if (fineLevel.isAllRegular()) {
        m_allregular = true;
        m_ok = true;
        return;
    }

    BL_PROFILE("EB2::GShopLevel()-coarse");

    const BoxArray& fine_grids = fineLevel.m_grids;
    const BoxArray& fine_covered_grids = fineLevel.m_covered_grids;

    const int coarse_ratio = 2;
    const int min_width = 8;
    bool coarsenable = fine_grids.coarsenable(coarse_ratio, min_width)
        && (fine_covered_grids.empty() || fine_covered_grids.coarsenable(coarse_ratio));

    m_ngrow = amrex::coarsen(fineLevel.m_ngrow,2);
    if (amrex::scale(m_ngrow,2) != fineLevel.m_ngrow) {
        m_ngrow = IntVect::TheZeroVector();
    }

    if (coarsenable)
    {
        int ierr = coarsenFromFine(fineLevel, true);
        m_ok = (ierr == 0);
    }
    else
    {
        Level fine_level_2(is, fineLevel.m_geom);
        fine_level_2.prepareForCoarsening(fineLevel, max_grid_size, amrex::scale(m_ngrow,2));
        int ierr = coarsenFromFine(fine_level_2, false);
        m_ok = (ierr == 0);
    }
}

}}

#endif
