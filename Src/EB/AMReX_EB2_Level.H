#ifndef AMREX_EB2_LEVEL_H_
#define AMREX_EB2_LEVEL_H_

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_VisMF.H>
#include <AMReX_Array.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EB2_MultiGFab.H>
#include <AMReX_EB2_C.H>
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

public: // for cuda
    int coarsenFromFine (Level& fineLevel, bool fill_boundary);
    void buildCellFlag ();
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

    Real small_volfrac = 1.e-14;
    {
        ParmParse pp("eb2");
        pp.query("small_volfrac", small_volfrac);
    }

    // make sure ngrow is multiple of 16
    m_ngrow = IntVect{static_cast<int>(std::ceil(ngrow/16.)) * 16};

    Box domain_grown = geom.Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
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
        int box_type = gshop.getBoxType(gbx, geom, RunOn::Gpu);
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
    MFInfo mf_info;
    mf_info.SetTag("EB2::Level");
    m_cellflag.define(m_grids, m_dmap, 1, ng, mf_info);
    m_volfrac.define(m_grids, m_dmap, 1, ng, mf_info);
    m_centroid.define(m_grids, m_dmap, AMREX_SPACEDIM, ng, mf_info);
    m_bndryarea.define(m_grids, m_dmap, 1, ng, mf_info);
    m_bndrycent.define(m_grids, m_dmap, AMREX_SPACEDIM, ng, mf_info);
    m_bndrynorm.define(m_grids, m_dmap, AMREX_SPACEDIM, ng, mf_info);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_areafrac[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)),
                                m_dmap, 1, ng, mf_info);
        m_facecent[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)),
                                m_dmap, AMREX_SPACEDIM-1, ng, mf_info);
    }

    const auto dx = geom.CellSizeArray();
    const auto dxinv = geom.InvCellSizeArray();
    const auto problo = geom.ProbLoArray();

    RunOn gshop_run_on = (Gpu::inLaunchRegion() && gshop.isGPUable())
        ? RunOn::Gpu : RunOn::Cpu;

    bool hybrid = Gpu::inLaunchRegion() and (gshop_run_on == RunOn::Cpu);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        Array<BaseFab<Real>, AMREX_SPACEDIM> intercept;
        Array<BaseFab<Real>, AMREX_SPACEDIM> M2;
        EBCellFlagFab cellflagtmp;

        for (MFIter mfi(m_mgf); mfi.isValid(); ++mfi)
        {
            auto& gfab = m_mgf[mfi];
            const Box& vbx = gfab.validbox();

            auto& levelset = gfab.getLevelSet();
            gshop.fillFab(levelset, geom, gshop_run_on);

            if (hybrid) levelset.prefetchToDevice();

            auto& cellflag = m_cellflag[mfi];

            gfab.buildTypes(cellflag);

            Array4<Real const> const& lst = levelset.const_array();
            Array4<EBCellFlag> const& cfg = m_cellflag.array(mfi);
            Array4<Real> const& vfr = m_volfrac.array(mfi);
            Array4<Real> const& ctr = m_centroid.array(mfi);
            Array4<Real> const& bar = m_bndryarea.array(mfi);
            Array4<Real> const& bct = m_bndrycent.array(mfi);
            Array4<Real> const& bnm = m_bndrynorm.array(mfi);
            AMREX_D_TERM(Array4<Real> const& apx = m_areafrac[0].array(mfi);,
                         Array4<Real> const& apy = m_areafrac[1].array(mfi);,
                         Array4<Real> const& apz = m_areafrac[2].array(mfi););
            AMREX_D_TERM(Array4<Real> const& fcx = m_facecent[0].array(mfi);,
                         Array4<Real> const& fcy = m_facecent[1].array(mfi);,
                         Array4<Real> const& fcz = m_facecent[2].array(mfi););

            auto& facetype = gfab.getFaceType();
            AMREX_D_TERM(Array4<Type_t> const& ftx = facetype[0].array();,
                         Array4<Type_t> const& fty = facetype[1].array();,
                         Array4<Type_t> const& ftz = facetype[2].array(););

#if (AMREX_SPACEDIM == 3)
            auto& edgetype = gfab.getEdgeType();
            Array4<Type_t> const& xdg = edgetype[0].array();
            Array4<Type_t> const& ydg = edgetype[1].array();
            Array4<Type_t> const& zdg = edgetype[2].array();

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                intercept[idim].resize(edgetype[idim].box());
            }
            Elixir xinter_eli = intercept[0].elixir();
            Elixir yinter_eli = intercept[1].elixir();
            Elixir zinter_eli = intercept[2].elixir();
            Array4<Real> const& xip = intercept[0].array();
            Array4<Real> const& yip = intercept[1].array();
            Array4<Real> const& zip = intercept[2].array();

            if (hybrid) {
                Gpu::streamSynchronize();
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    edgetype[idim].prefetchToHost();
                    intercept[idim].prefetchToHost();
                }
            }

            gshop.getIntercept(intercept, edgetype, geom, gshop_run_on);

            if (hybrid) {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    edgetype[idim].prefetchToDevice();
                    intercept[idim].prefetchToDevice();
                }
            }

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const Box& b = facetype[idim].box();
                M2[idim].resize(b,3);
            }
            Elixir xm2_eli = M2[0].elixir();
            Elixir ym2_eli = M2[1].elixir();
            Elixir zm2_eli = M2[2].elixir();
            Array4<Real> const& xm2 = M2[0].array();
            Array4<Real> const& ym2 = M2[1].array();
            Array4<Real> const& zm2 = M2[2].array();

            build_faces(vbx, cfg, ftx, fty, ftz, xdg, ydg, zdg, lst,
                        xip, yip, zip, apx, apy, apz, fcx, fcy, fcz,
                        xm2, ym2, zm2, dx, problo);

            cellflagtmp.resize(m_cellflag[mfi].box());
            Elixir cellflagtmp_eli = cellflagtmp.elixir();
            Array4<EBCellFlag> const& cfgtmp = cellflagtmp.array();

            build_cells(vbx, cfg, ftx, fty, ftz, apx, apy, apz,
                        fcx, fcy, fcz, xm2, ym2, zm2, vfr, ctr,
                        bar, bct, bnm, cfgtmp, small_volfrac);

#elif (AMREX_SPACEDIM == 2)
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                intercept[idim].resize(facetype[idim].box());
            }
            Elixir xinter_eli = intercept[0].elixir();
            Elixir yinter_eli = intercept[1].elixir();
            Array4<Real> const& xip = intercept[0].array();
            Array4<Real> const& yip = intercept[1].array();

            if (hybrid) {
                Gpu::streamSynchronize();
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    facetype[idim].prefetchToHost();
                    intercept[idim].prefetchToHost();
                }
            }

            gshop.getIntercept(intercept, facetype, geom, gshop_run_on);

            if (hybrid) {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    facetype[idim].prefetchToDevice();
                    intercept[idim].prefetchToDevice();
                }
            }

            build_faces(vbx, cfg, ftx, fty, lst, xip, yip, apx, apy, fcx, fcy, dx, problo);

            build_cells(vbx, cfg, ftx, fty, apx, apy, vfr, ctr, bar, bct, bnm, small_volfrac);
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
