#ifndef AMREX_EB2_LEVEL_H_
#define AMREX_EB2_LEVEL_H_
#include <AMReX_Config.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_LayoutData.H>
#include <AMReX_VisMF.H>
#include <AMReX_Array.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_MultiCutFab.H>
#if (AMREX_SPACEDIM == 3)
#  include <AMReX_MarchingCubes.H>
#endif
#include <AMReX_EB2_MultiGFab.H>
#include <AMReX_EB2_C.H>
#include <AMReX_EB2_IF_AllRegular.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <cmath>
#include <limits>
#include <map>
#include <type_traits>
#include <unordered_map>

namespace amrex
{
    class STLtools;
}

namespace amrex::EB2 {

class IndexSpace;
template <typename G> class GShopLevel;

class Level
{
public:

    //! \return True if every cell in this level is regular (no cut cells).
    bool isAllRegular () const noexcept { return m_allregular; }
    //! \return True if the level successfully built its EB description.
    bool isOK () const noexcept { return m_ok; }
    /**
     * \brief Fill \p cellflag with EBCellFlag data using geometry \p geom.
     *
     * \param cellflag Destination flag array.
     * \param geom     Geometry describing the sampling box.
     */
    void fillEBCellFlag (FabArray<EBCellFlagFab>& cellflag, const Geometry& geom) const;
    /**
     * \brief Populate volume fractions into \p vfrac using \p geom.
     */
    void fillVolFrac (MultiFab& vfrac, const Geometry& geom) const;
    /**
     * \brief Populate cut-cell centroids into a MultiCutFab.
     */
    void fillCentroid (MultiCutFab& centroid, const Geometry& geom) const;
    //! \overload
    void fillCentroid (   MultiFab& centroid, const Geometry& geom) const;
    /**
     * \brief Populate boundary-face areas into \p bndryarea.
     */
    void fillBndryArea (MultiCutFab& bndryarea, const Geometry& geom) const;
    //! \overload
    void fillBndryArea (   MultiFab& bndryarea, const Geometry& geom) const;
    /**
     * \brief Populate boundary centroids into \p bndrycent.
     */
    void fillBndryCent (MultiCutFab& bndrycent, const Geometry& geom) const;
    //! \overload
    void fillBndryCent (   MultiFab& bndrycent, const Geometry& geom) const;
    /**
     * \brief Populate boundary normals into \p bndrynorm.
     */
    void fillBndryNorm (MultiCutFab& bndrynorm, const Geometry& geom) const;
    //! \overload
    void fillBndryNorm (   MultiFab& bndrynorm, const Geometry& geom) const;
    /**
     * \brief Populate face area fractions for each direction.
     */
    void fillAreaFrac (Array<MultiCutFab*,AMREX_SPACEDIM> const& areafrac, const Geometry& geom) const;
    //! \overload
    void fillAreaFrac (Array<   MultiFab*,AMREX_SPACEDIM> const& areafrac, const Geometry& geom) const;
    /**
     * \brief Populate face centroids for each direction.
     */
    void fillFaceCent (Array<MultiCutFab*,AMREX_SPACEDIM> const& facecent, const Geometry& geom) const;
    //! \overload
    void fillFaceCent (Array<   MultiFab*,AMREX_SPACEDIM> const& facecent, const Geometry& geom) const;
    /**
     * \brief Populate edge centroids for each direction.
     */
    void fillEdgeCent (Array<MultiCutFab*,AMREX_SPACEDIM> const& edgecent, const Geometry& geom) const;
    //! \overload
    void fillEdgeCent (Array<   MultiFab*,AMREX_SPACEDIM> const& edgecent, const Geometry& geom) const;
    /**
     * \brief Write the implicit function values into \p levelset.
     */
    void fillLevelSet (MultiFab& levelset, const Geometry& geom) const;

    //! \return Grid layout associated with this level.
    const BoxArray& boxArray () const noexcept { return m_grids; }
    //! \return Distribution mapping describing how grids are assigned to ranks.
    const DistributionMapping& DistributionMap () const noexcept { return m_dmap; }

    //! Construct a Level bound to EB index space \p is and Geometry \p geom.
    Level (IndexSpace const* is, const Geometry& geom) : m_geom(geom), m_parent(is) {}
    /**
     * \brief Prepare for coarsening by copying from \p rhs, respecting max grid size and grow vector.
     */
    void prepareForCoarsening (const Level& rhs, int max_grid_size, IntVect const& ngrow);

    //! \return Reference to the level Geometry.
    const Geometry& Geom () const noexcept { return m_geom; }
    //! \return Pointer to the EB IndexSpace.
    IndexSpace const* getEBIndexSpace () const noexcept { return m_parent; }

    //! \return EB grow-vector.
    IntVect const& nGrowVect () const noexcept { return m_ngrow; }

    /**
     * \brief Write this level’s EB data to a checkpoint file \p fname.
     *
     * \param fname              Target checkpoint file name.
     * \param extend_domain_face Whether domain faces were extended.
     * \param max_grid_size      Grid blocking factor used during build.
     */
    void write_to_chkpt_file (const std::string& fname, bool extend_domain_face, int max_grid_size) const;

    //! \return True if EB data (flags, area fractions, etc.) are available.
    bool hasEBInfo () const noexcept { return m_has_eb_info; }
    //! Populate the cut-cell mask into \p cutcellmask.
    void fillCutCellMask (iMultiFab& cutcellmask, const Geometry& geom) const;

    /**
     * \brief Shift this level by \p ncells cells along direction \p direction (for periodic tiling).
     */
    void setShift (int direction, int ncells);

// public: // for cuda
    /**
     * \brief Coarsen EB data from \p fineLevel into this level.
     *
     * \param fineLevel Source level.
     * \param fill_boundary Whether to fill boundary data after coarsening.
     * \return 0 on success.
     */
    int coarsenFromFine (Level& fineLevel, bool fill_boundary);
    //! Build EBCellFlag data from the current implicit function field.
    void buildCellFlag ();
    //! Build the cut-cell mask by coarsening from \p fine_level.
    void buildCutCellMask (Level const& fine_level);

protected:

    Geometry m_geom;
    IntVect  m_ngrow;
    BoxArray m_grids;
    BoxArray m_covered_grids;
    DistributionMapping m_dmap;
    MultiFab m_sdf;
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
    Array<MultiFab,AMREX_SPACEDIM> m_edgecent;
    iMultiFab m_cutcellmask;
    bool m_allregular = false;
    bool m_ok = false;
    bool m_has_eb_info = true;
    IntVect m_shift{0};
    IndexSpace const* m_parent;

private:
    template <typename G> friend class GShopLevel;

    // Need this function to work around a gcc bug.
    void setRegularLevel () {
        m_allregular = true;
        m_ok = true;
        m_has_eb_info = false;
    }
};

template <typename G>
class GShopLevel
    : public Level
{
public:
    /**
     * \brief Build a level directly from a GeometryShop object.
     */
    GShopLevel (IndexSpace const* is, G const& gshop, const Geometry& geom, int max_grid_size,
                int ngrow, bool extend_domain_face, int num_crse_opt);
    /**
     * \brief Construct a coarse level by coarsening \p fineLevel.
     */
    GShopLevel (IndexSpace const* is, int ilev, int max_grid_size, int ngrow,
                const Geometry& geom, GShopLevel<G>& fineLevel);
    //! Construct an empty level (usually regular) bound to \p geom.
    GShopLevel (IndexSpace const* is, const Geometry& geom);
    //! Define data for the finest level using the supplied GeometryShop.
    void define_fine (G const& gshop, const Geometry& geom,
                      int max_grid_size, int ngrow, bool extend_domain_face, int num_crse_opt);

#if (AMREX_SPACEDIM == 3)
    /**
     * \brief Define a fine level for multi-valued multi-cut STL geometries (3-D only)
     */
    template <typename GS = G>
    requires (std::same_as<GS,STLtools>)
    void define_fine_marching_cubes (GS const& gshop, const Geometry& geom,
                                     int max_grid_size, int ngrow,
                                     bool extend_domain_face, int num_crse_opt);
#endif

    /**
     * \brief Build a regular (no EB) level.
     */
    static GShopLevel<G>
    makeAllRegular(IndexSpace const* is, const Geometry& geom)
    {
        GShopLevel<G> r(is, geom);
        r.setRegularLevel();
        return r;
    }

private:

    Box m_bounding_box;

    void prepare_grids (G const& gshop, Geometry const& geom, int max_grid_size, int ngrow,
                        bool extend_domain_face, int num_crse_opt);

};

template <typename G>
GShopLevel<G>::GShopLevel (IndexSpace const* is, const Geometry& geom)
    : Level(is, geom)
{}

template <typename G>
GShopLevel<G>::GShopLevel (IndexSpace const* is, G const& gshop, const Geometry& geom,
                           int max_grid_size, int ngrow, bool extend_domain_face, int num_crse_opt)
    : Level(is, geom)
{
    if (std::is_same_v<typename G::FunctionType, AllRegularIF>) {
        m_allregular = true;
        m_ok = true;
        return;
    }

    define_fine(gshop, geom, max_grid_size, ngrow, extend_domain_face, num_crse_opt);
}

template <typename G>
void
GShopLevel<G>::prepare_grids (G const& gshop, const Geometry& geom,
                              int max_grid_size, int ngrow, bool extend_domain_face, int num_crse_opt)
{
    if (amrex::Verbose() > 0 && extend_domain_face == false) {
        amrex::Print() << "AMReX WARNING: extend_domain_face=false is not recommended!\n";
    }

    BL_PROFILE("EB2::GShopLevel()-prepare_grids");

    // make sure ngrow is multiple of 16
    m_ngrow = IntVect{static_cast<int>(std::ceil(ngrow/16.)) * 16};

    Box const& domain = geom.Domain();
    Box domain_grown = domain;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            m_ngrow[idim] = 0;
        } else {
            m_ngrow[idim] = std::min(m_ngrow[idim], domain_grown.length(idim));
        }
    }
    domain_grown.grow(m_ngrow);
    m_bounding_box = (extend_domain_face) ? domain : domain_grown;
    m_bounding_box.surroundingNodes();

    BoxList cut_boxes;
    BoxList covered_boxes;

    const int nprocs = ParallelDescriptor::NProcs();
    const int iproc = ParallelDescriptor::MyProc();

    num_crse_opt = std::max(0,std::min(8,num_crse_opt));
    for (int clev = num_crse_opt; clev >= 0; --clev) {
        IntVect crse_ratio(1 << clev);
        if (domain.coarsenable(crse_ratio)) {
            Box const& crse_bounding_box = amrex::coarsen(m_bounding_box, crse_ratio);
            Geometry const& crse_geom = amrex::coarsen(geom, crse_ratio);
            BoxList test_boxes;
            if (cut_boxes.isEmpty()) {
                covered_boxes.clear();
                test_boxes = BoxList(crse_geom.Domain());
                test_boxes.maxSize(max_grid_size);
            } else {
                test_boxes.swap(cut_boxes);
                test_boxes.coarsen(crse_ratio);
                test_boxes.maxSize(max_grid_size);
            }

            const Long nboxes = test_boxes.size();
            const auto& boxes = test_boxes.data();
            for (Long i = iproc; i < nboxes; i += nprocs) {
                const Box& vbx = boxes[i];
                const Box& gbx = amrex::surroundingNodes(amrex::grow(vbx,1));
                auto box_type = gshop.getBoxType(gbx&crse_bounding_box,crse_geom,RunOn::Gpu);
                if (box_type == gshop.allcovered) {
                    covered_boxes.push_back(amrex::refine(vbx, crse_ratio));
                } else if (box_type == gshop.mixedcells) {
                    cut_boxes.push_back(amrex::refine(vbx, crse_ratio));
                }
            }

            amrex::AllGatherBoxes(cut_boxes.data());
        }
    }

    amrex::AllGatherBoxes(covered_boxes.data());

    if (m_ngrow != 0) {
        auto grow_at_domain_boundary = [&] (BoxList& bl)
        {
            const IntVect& domlo = domain.smallEnd();
            const IntVect& domhi = domain.bigEnd();
            for (auto& b : bl) {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    if (m_ngrow[idim] != 0) {
                        if (b.smallEnd(idim) == domlo[idim]) {
                            b.growLo(idim,m_ngrow[idim]);
                        }
                        if (b.bigEnd(idim) == domhi[idim]) {
                            b.growHi(idim,m_ngrow[idim]);
                        }
                    }
                }
            }
        };
        grow_at_domain_boundary(covered_boxes);
        grow_at_domain_boundary(cut_boxes);
    }

    if ( cut_boxes.isEmpty() &&
        !covered_boxes.isEmpty())
    {
        amrex::Abort("AMReX_EB2_Level.H: Domain is completely covered");
    }

    if (!covered_boxes.isEmpty()) {
        if (num_crse_opt > 2) { // don't want the box too big
            covered_boxes.maxSize(max_grid_size*4);
        }
        m_covered_grids = BoxArray(std::move(covered_boxes));
    }

    if (cut_boxes.isEmpty()) {
        m_grids = BoxArray();
        m_dmap = DistributionMapping();
        m_allregular = true;
        m_ok = true;
    } else {
        m_grids = BoxArray(std::move(cut_boxes));
        m_dmap = DistributionMapping(m_grids);
    }
}

template <typename G>
void
GShopLevel<G>::define_fine (G const& gshop, const Geometry& geom,
                            int max_grid_size, int ngrow, bool extend_domain_face, int num_crse_opt)
{
    BL_PROFILE("EB2::GShopLevel()-fine");

#ifdef AMREX_USE_FLOAT
    Real small_volfrac = 1.e-5_rt;
#else
    Real small_volfrac = 1.e-14;
#endif
    bool cover_multiple_cuts = false;
    int maxiter = 32;
    {
        ParmParse pp("eb2");
        pp.queryAdd("small_volfrac", small_volfrac);
        pp.queryAdd("cover_multiple_cuts", cover_multiple_cuts);
        pp.queryAdd("maxiter", maxiter);
    }
    maxiter = std::min(100000, maxiter);

    prepare_grids(gshop, geom, max_grid_size, ngrow, extend_domain_face, num_crse_opt);

    if (m_allregular) { return; }

    m_mgf.define(m_grids, m_dmap);
    const int ng = GFab::ng;
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
        IntVect edge_type{1}; edge_type[idim] = 0;
        m_edgecent[idim].define(amrex::convert(m_grids, edge_type), m_dmap, 1, ng, mf_info);
    }

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    auto bounding_box = m_bounding_box;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (!extend_domain_face || geom.isPeriodic(idim)) {
            bounding_box.grow(idim,GFab::ng);
        }
    }

    RunOn gshop_run_on = (Gpu::inLaunchRegion() && gshop.isGPUable())
        ? RunOn::Gpu : RunOn::Cpu;

    bool hybrid = Gpu::inLaunchRegion() && (gshop_run_on == RunOn::Cpu);
    amrex::ignore_unused(hybrid);

    int iter = 0;
    for (; iter < maxiter; ++iter)
    {
        int nsmallcells = 0;
        int nmulticuts = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(+:nsmallcells,nmulticuts)
#endif
        {
#if (AMREX_SPACEDIM == 3)
            Array<BaseFab<Real>, AMREX_SPACEDIM> M2;
            EBCellFlagFab cellflagtmp;
#endif
            for (MFIter mfi(m_mgf); mfi.isValid(); ++mfi)
            {
                auto& gfab = m_mgf[mfi];
                const Box& vbx = gfab.validbox();

                auto& levelset = gfab.getLevelSet();
                if (iter == 0) {
                    gshop.fillFab(levelset, geom, gshop_run_on, bounding_box);
#ifdef AMREX_USE_GPU
                    if (hybrid) {
                        levelset.prefetchToDevice();
                    }
#endif
                }

                auto& cellflag = m_cellflag[mfi];

                gfab.buildTypes(cellflag);

                Array4<Real const> const& clst = levelset.const_array();
                Array4<Real      > const&  lst = levelset.array();
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

                int nmc = 0;
                int nsm = 0;

#if (AMREX_SPACEDIM == 3)
                auto& edgetype = gfab.getEdgeType();
                Array4<Type_t const> const& xdg = edgetype[0].const_array();
                Array4<Type_t const> const& ydg = edgetype[1].const_array();
                Array4<Type_t const> const& zdg = edgetype[2].const_array();

                Array4<Real> const& xip = m_edgecent[0].array(mfi);
                Array4<Real> const& yip = m_edgecent[1].array(mfi);
                Array4<Real> const& zip = m_edgecent[2].array(mfi);

                if (iter == 0) {
#ifdef AMREX_USE_GPU
                    if (hybrid) {
                        Gpu::streamSynchronize();
                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                            edgetype[idim].prefetchToHost();
                            m_edgecent[idim][mfi].prefetchToHost();
                        }
                    }
#endif

                    gshop.getIntercept({xip,yip,zip}, {xdg,ydg,zdg}, clst,
                                       geom, gshop_run_on, bounding_box);

#ifdef AMREX_USE_GPU
                    if (hybrid) {
                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                            edgetype[idim].prefetchToDevice();
                            m_edgecent[idim][mfi].prefetchToDevice();
                        }
                    }
#endif
                } else {
                    gshop.updateIntercept({xip,yip,zip}, {xdg,ydg,zdg}, clst, geom);
                }

                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    const Box& b = facetype[idim].box();
                    M2[idim].resize(b,3);
                }
                Array4<Real> const& xm2 = M2[0].array();
                Array4<Real> const& ym2 = M2[1].array();
                Array4<Real> const& zm2 = M2[2].array();

                nmc = build_faces(vbx, cfg, ftx, fty, ftz, xdg, ydg, zdg, lst,
                                  xip, yip, zip, apx, apy, apz, fcx, fcy, fcz,
                                  xm2, ym2, zm2, dx, problo, cover_multiple_cuts);

                cellflagtmp.resize(m_cellflag[mfi].box());
                Elixir cellflagtmp_eli = cellflagtmp.elixir();
                Array4<EBCellFlag> const& cfgtmp = cellflagtmp.array();

                build_cells(vbx, cfg, ftx, fty, ftz, apx, apy, apz,
                            fcx, fcy, fcz, xm2, ym2, zm2, dx, vfr, ctr,
                            bar, bct, bnm, cfgtmp, lst,
                            small_volfrac, geom, extend_domain_face, cover_multiple_cuts,
                            nsm, nmc);

                // Because it is used in a synchronous reduction kernel in
                // build_cells, we do not need to worry about M2's lifetime.
                // But we still need to use Elixir to extend the life of
                // cellflagtmp.

#elif (AMREX_SPACEDIM == 2)
                Array4<Real> const& xip = m_edgecent[0].array(mfi);
                Array4<Real> const& yip = m_edgecent[1].array(mfi);

                if (iter == 0) {
#ifdef AMREX_USE_GPU
                    if (hybrid) {
                        Gpu::streamSynchronize();
                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                            facetype[idim].prefetchToHost();
                            m_edgecent[idim][mfi].prefetchToHost();
                        }
                    }
#endif

                    //                         yes, factype[1] and then [0]
                    gshop.getIntercept({xip,yip},
                                       {facetype[1].const_array(), facetype[0].const_array()},
                                       clst, geom, gshop_run_on, bounding_box);

#ifdef AMREX_USE_GPU
                    if (hybrid) {
                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                            facetype[idim].prefetchToDevice();
                            m_edgecent[idim][mfi].prefetchToDevice();
                        }
                    }
#endif
                } else {
                    gshop.updateIntercept({xip,yip},
                                          {facetype[1].const_array(), facetype[0].const_array()},
                                          clst, geom);
                }

                nmc = build_faces(vbx, cfg, ftx, fty, lst, xip, yip, apx, apy, fcx, fcy,
                                  dx, problo, cover_multiple_cuts, nsm);

                build_cells(vbx, cfg, ftx, fty, apx, apy, dx, vfr, ctr,
                            bar, bct, bnm, lst, small_volfrac, geom, extend_domain_face,
                            nsm, nmc);
#endif
                nsmallcells += nsm;
                nmulticuts  += nmc;
            }
        }

        ParallelAllReduce::Sum<int>({nsmallcells,nmulticuts}, ParallelContext::CommunicatorSub());
        if (nsmallcells == 0 && nmulticuts == 0) {
            break;
        } else {
            auto ls = m_mgf.getLevelSet();
            // This is an alias MultiFab, therefore FillBoundary on it is fine.
            ls.FillBoundary(geom.periodicity());
            if (amrex::Verbose() > 0) {
                if (nsmallcells) {
                    amrex::Print() << "AMReX EB: Iter. " << iter+1 << " fixed " << nsmallcells
                                   << " small cells" << '\n';
                }
                if (nmulticuts) {
                    amrex::Print() << "AMReX EB: Iter. " << iter+1 << " fixed " << nmulticuts
                                   << " multicuts" << '\n';
                }
            }
        }
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < maxiter, "EB: failed to fix small cells");

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(m_mgf); mfi.isValid(); ++mfi)
    {
        auto& gfab = m_mgf[mfi];
        auto const& levelset = gfab.getLevelSet();
        Array4<Real const> const& clst = levelset.const_array();

        AMREX_D_TERM(Array4<Real> const& xip = m_edgecent[0].array(mfi);,
                     Array4<Real> const& yip = m_edgecent[1].array(mfi);,
                     Array4<Real> const& zip = m_edgecent[2].array(mfi);)
#if (AMREX_SPACEDIM == 3)
        auto const& edgetype = gfab.getEdgeType();
        Array4<Type_t const> const& xdg = edgetype[0].const_array();
        Array4<Type_t const> const& ydg = edgetype[1].const_array();
        Array4<Type_t const> const& zdg = edgetype[2].const_array();

        intercept_to_edge_centroid(xip, yip, zip, xdg, ydg, zdg, clst, dx, problo);

#elif (AMREX_SPACEDIM == 2)
        auto& facetype = gfab.getFaceType();
        Array4<Type_t const> const& ftx = facetype[0].const_array();
        Array4<Type_t const> const& fty = facetype[1].const_array();
        //                                  fty then ftx
        intercept_to_edge_centroid(xip, yip, fty, ftx, clst, dx, problo);
#endif
    }

    m_levelset = m_mgf.getLevelSet();

    m_ok = true;
}

#if (AMREX_SPACEDIM == 3)
template <typename G>
template <typename GS>
    requires(std::same_as<GS, STLtools>)
void GShopLevel<G>::define_fine_marching_cubes (GS const& gshop, const Geometry& geom,
                                                int max_grid_size, int ngrow,
                                                bool extend_domain_face, int num_crse_opt)
{
    BL_PROFILE("EB2::GShopLevel()-fine-marching-cubes");

    auto const cell_size = geom.CellSizeArray();
    Real const max_cell_size = amrex::max(cell_size[0], amrex::max(cell_size[1], cell_size[2]));
    Real const cubic_tolerance = 16.0_rt * std::numeric_limits<Real>::epsilon() * max_cell_size;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        geom.Coord() == CoordSys::cartesian &&
            std::abs(cell_size[0] - cell_size[1]) <= cubic_tolerance &&
            std::abs(cell_size[0] - cell_size[2]) <= cubic_tolerance,
        "Marching-cubes EB construction requires a 3D Cartesian grid with dx == "
        "dy == dz");

#ifdef AMREX_USE_FLOAT
    Real small_volfrac = 1.e-5_rt;
#else
    Real small_volfrac = 1.e-14_rt;
#endif
    int maxiter = 32;
    bool cover_multiple_cuts = false;
    {
        ParmParse pp("eb2");
        pp.queryAdd("small_volfrac", small_volfrac);
        pp.queryAdd("maxiter", maxiter);
        pp.queryAdd("cover_multiple_cuts", cover_multiple_cuts);
    }
    maxiter = std::min(100000, maxiter);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(maxiter > 0, "eb2.maxiter must be positive");

    // Keep only the cut-box discovery from the legacy path. Every EB field
    // below is allocated and populated from the MC working SDF.
    prepare_grids(gshop, geom, max_grid_size, ngrow, extend_domain_face, num_crse_opt);
    if (m_allregular) {
        return;
    }

    MC::Initialize();

    int const ng = GFab::ng;
    MFInfo mf_info;
    mf_info.SetTag("EB2::Level-MC");
    m_cellflag.define(m_grids, m_dmap, 1, ng, mf_info);
    m_volfrac.define(m_grids, m_dmap, 1, ng, mf_info);
    m_centroid.define(m_grids, m_dmap, AMREX_SPACEDIM, ng, mf_info);
    m_bndryarea.define(m_grids, m_dmap, 1, ng, mf_info);
    m_bndrycent.define(m_grids, m_dmap, AMREX_SPACEDIM, ng, mf_info);
    m_bndrynorm.define(m_grids, m_dmap, AMREX_SPACEDIM, ng, mf_info);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_areafrac[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)), m_dmap,
                                1, ng, mf_info);
        m_facecent[idim].define(amrex::convert(m_grids, IntVect::TheDimensionVector(idim)), m_dmap,
                                AMREX_SPACEDIM - 1, ng, mf_info);
        IntVect edge_type{1};
        edge_type[idim] = 0;
        m_edgecent[idim].define(amrex::convert(m_grids, edge_type), m_dmap, 1, ng, mf_info);
    }

    BoxArray const nodal_grids = amrex::convert(m_grids, IntVect::TheNodeVector());
    m_sdf.define(nodal_grids, m_dmap, 1, IntVect(2), mf_info);
    gshop.fillMarchingCubesLevelSet(m_sdf, m_sdf.nGrowVect(), geom);

    std::map<int, std::unique_ptr<MC::MCFab>> marching_cubes;
    for (MFIter mfi(m_sdf, MFItInfo().DisableDeviceSync()); mfi.isValid(); ++mfi) {
        marching_cubes[mfi.index()] = std::make_unique<MC::MCFab>();
        auto& mc_fab = *marching_cubes[mfi.index()];
        mc_fab.defineEdgeIntersections(m_sdf[mfi].box());
        gshop.fillMarchingCubesEdgeIntersections(m_sdf[mfi], geom, mc_fab);
    }
    auto const& mc_fabs = marching_cubes;

    iMultiFab rejected_cells(m_grids, m_dmap, 1, 2, mf_info);
    Array<iMultiFab, AMREX_SPACEDIM> rejected_faces;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        rejected_faces[idim].define(m_areafrac[idim].boxArray(), m_dmap, 1, 2, mf_info);
    }
    Box const domain = geom.Domain();
    GpuArray<int, 3> const is_periodic{geom.isPeriodic(0), geom.isPeriodic(1), geom.isPeriodic(2)};

    bool converged = false;
    for (int iter = 0; iter < maxiter; ++iter) {
        rejected_cells.setVal(0);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            rejected_faces[idim].setVal(0);
        }

        MFItInfo info{};
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
        info.SetDynamic(true);
#pragma omp parallel
#endif
        for (MFIter mfi(m_sdf, info); mfi.isValid(); ++mfi) {
            MC::marching_cubes(geom, m_sdf[mfi], *mc_fabs.at(mfi.index()));
        }

        m_volfrac.setVal(1.0_rt, 0, 1, m_volfrac.nGrowVect());
        m_centroid.setVal(0.0_rt, 0, AMREX_SPACEDIM, m_centroid.nGrowVect());
        m_bndryarea.setVal(0.0_rt, 0, 1, m_bndryarea.nGrowVect());
        m_bndrycent.setVal(-1.0_rt, 0, AMREX_SPACEDIM, m_bndrycent.nGrowVect());
        m_bndrynorm.setVal(0.0_rt, 0, AMREX_SPACEDIM, m_bndrynorm.nGrowVect());

        int face_geometry_errors = 0;
        int cell_geometry_errors = 0;
        int face_rejections = 0;
        int topology_rejections = 0;
        int small_cell_rejections = 0;
        int extension_rejections = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())                                                 \
    reduction(+ : face_geometry_errors, cell_geometry_errors, face_rejections,                     \
                  topology_rejections, small_cell_rejections, extension_rejections)
#endif
        for (MFIter mfi(m_sdf, info); mfi.isValid(); ++mfi) {
            Box const vbx = amrex::enclosedCells(mfi.validbox());
            Box const gbx = amrex::grow(vbx, 1);
            auto const sdf = m_sdf.const_array(mfi);
            auto const vfrac = m_volfrac.array(mfi);
            auto const vcent = m_centroid.array(mfi);
            auto const barea = m_bndryarea.array(mfi);
            auto const bcent = m_bndrycent.array(mfi);
            auto const bnorm = m_bndrynorm.array(mfi);

            ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                int nfluid = 0;
                nfluid += sdf(i, j, k) > 0.0_rt;
                nfluid += sdf(i + 1, j, k) > 0.0_rt;
                nfluid += sdf(i, j + 1, k) > 0.0_rt;
                nfluid += sdf(i + 1, j + 1, k) > 0.0_rt;
                nfluid += sdf(i, j, k + 1) > 0.0_rt;
                nfluid += sdf(i + 1, j, k + 1) > 0.0_rt;
                nfluid += sdf(i, j + 1, k + 1) > 0.0_rt;
                nfluid += sdf(i + 1, j + 1, k + 1) > 0.0_rt;
                vfrac(i, j, k) = nfluid == 0 ? 0.0_rt : (nfluid == 8 ? 1.0_rt : -1.0_rt);
                barea(i, j, k) = 0.0_rt;
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    vcent(i, j, k, n) = 0.0_rt;
                    bcent(i, j, k, n) = -1.0_rt;
                    bnorm(i, j, k, n) = 0.0_rt;
                }
            });

            face_geometry_errors += MC::build_face_fractions(
                gbx, *mc_fabs.at(mfi.index()), m_sdf[mfi], m_areafrac[0][mfi], m_areafrac[1][mfi],
                m_areafrac[2][mfi], m_facecent[0][mfi], m_facecent[1][mfi], m_facecent[2][mfi]);
            MC::build_edge_centroids(gbx, *mc_fabs.at(mfi.index()), m_sdf[mfi], m_edgecent[0][mfi],
                                     m_edgecent[1][mfi], m_edgecent[2][mfi]);
            cell_geometry_errors += MC::build_cell_fractions(
                vbx, geom, *mc_fabs.at(mfi.index()), m_sdf[mfi], m_areafrac[0][mfi],
                m_areafrac[1][mfi], m_areafrac[2][mfi], m_volfrac[mfi], m_centroid[mfi],
                m_bndryarea[mfi], m_bndrycent[mfi], m_bndrynorm[mfi]);

            face_rejections += MC::mark_faces_for_cleanup(
                vbx, *mc_fabs.at(mfi.index()), m_sdf[mfi], rejected_faces[0][mfi],
                rejected_faces[1][mfi], rejected_faces[2][mfi]);
            GpuArray<int, 2> const counts =
                MC::mark_cells_for_cleanup(vbx, *mc_fabs.at(mfi.index()), m_sdf[mfi],
                                           m_volfrac[mfi], small_volfrac, rejected_cells[mfi]);
            topology_rejections += counts[0];
            small_cell_rejections += counts[1];
            if (extend_domain_face) {
                extension_rejections += MC::mark_cells_for_domain_extension(
                    vbx, domain, is_periodic, m_volfrac[mfi], rejected_cells[mfi]);
            }
        }

        ParallelAllReduce::Sum<int>({face_geometry_errors, cell_geometry_errors, face_rejections,
                                     topology_rejections, small_cell_rejections,
                                     extension_rejections},
                                    ParallelContext::CommunicatorSub());
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(face_geometry_errors == 0,
                                         "Marching-cubes EB encountered " +
                                             std::to_string(face_geometry_errors) +
                                             " inconsistent Cartesian-face geometry records");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            cell_geometry_errors == 0 || topology_rejections > 0,
            "Marching-cubes EB could not map an invalid cell-moment record to "
            "the nodal repair set");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            cover_multiple_cuts || face_rejections + topology_rejections == 0,
            "Marching-cubes EB found " + std::to_string(face_rejections + topology_rejections) +
                " faces or cells whose topology cannot be represented by one VoF; "
                "set eb2.cover_multiple_cuts=1 to enable legacy-style nodal "
                "repair");

        if (face_rejections + topology_rejections + small_cell_rejections + extension_rejections ==
            0) {
            converged = true;
            break;
        }

        rejected_cells.FillBoundary(geom.periodicity());
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            rejected_faces[idim].OverrideSync(geom.periodicity());
            rejected_faces[idim].FillBoundary(geom.periodicity());
        }

        Long changed_nodes = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(+ : changed_nodes)
#endif
        for (MFIter mfi(m_sdf, info); mfi.isValid(); ++mfi) {
            changed_nodes += MC::zero_nodes_for_cleanup(
                mfi.validbox(), rejected_cells[mfi], rejected_faces[0][mfi], rejected_faces[1][mfi],
                rejected_faces[2][mfi], m_sdf[mfi]);
        }
        ParallelAllReduce::Sum(changed_nodes, ParallelContext::CommunicatorSub());
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(changed_nodes > 0,
                                         "Marching-cubes EB repair found rejected "
                                         "geometry but changed no fluid nodes");

        m_sdf.OverrideSync(geom.periodicity());
        m_sdf.FillBoundary(geom.periodicity());
        if (amrex::Verbose() > 0) {
            if (small_cell_rejections > 0) {
                amrex::Print() << "AMReX MC EB: Iter. " << iter + 1 << " fixed "
                               << small_cell_rejections << " small cells\n";
            }
            if (face_rejections + topology_rejections > 0) {
                amrex::Print() << "AMReX MC EB: Iter. " << iter + 1 << " fixed "
                               << face_rejections + topology_rejections
                               << " unsupported faces/cells\n";
            }
            if (extension_rejections > 0) {
                amrex::Print() << "AMReX MC EB: Iter. " << iter + 1 << " extended "
                               << extension_rejections
                               << " covered cells beyond physical domain faces\n";
            }
            if (cell_geometry_errors > 0) {
                amrex::Print() << "AMReX MC EB: Iter. " << iter + 1 << " included "
                               << cell_geometry_errors
                               << " invalid cell-moment records in the repair set\n";
            }
        }
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(converged,
                                     "Marching-cubes EB failed to fix small or unsupported cells");

    m_volfrac.FillBoundary(geom.periodicity());
    m_centroid.FillBoundary(geom.periodicity());
    m_bndryarea.FillBoundary(geom.periodicity());
    m_bndrycent.FillBoundary(geom.periodicity());
    m_bndrynorm.FillBoundary(geom.periodicity());
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_areafrac[idim].FillBoundary(geom.periodicity());
        m_facecent[idim].FillBoundary(geom.periodicity());
        m_edgecent[idim].FillBoundary(geom.periodicity());
    }

    MFItInfo final_info{};
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
    final_info.SetDynamic(true);
#endif
    int geometry_errors = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(+ : geometry_errors)
#endif
    for (MFIter mfi(m_sdf, final_info); mfi.isValid(); ++mfi) {
        Box const vbx = amrex::enclosedCells(mfi.validbox());
        m_cellflag[mfi].setVal(EBCellFlag::TheDefaultCell());
        geometry_errors += MC::build_cell_topology(
            vbx, *mc_fabs.at(mfi.index()), m_sdf[mfi], m_cellflag[mfi], m_volfrac[mfi],
            m_areafrac[0][mfi], m_areafrac[1][mfi], m_areafrac[2][mfi]);
        m_cellflag[mfi].resetType(GFab::ng);
    }
    ParallelAllReduce::Sum(geometry_errors, ParallelContext::CommunicatorSub());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geometry_errors == 0,
                                     "Final repaired marching-cubes topology is inconsistent");

    // The converged triangulation and the public EB data now describe the same
    // repaired domain.
    std::string stl_output;
    ParmParse("eb2").query("mc_stl_file", stl_output);
    if (!stl_output.empty()) {
        MC::write_stl(stl_output, marching_cubes, m_grids);
    }
    marching_cubes.clear();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(m_sdf, final_info); mfi.isValid(); ++mfi) {
        auto const phi = m_sdf.array(mfi);
        ParallelFor(m_sdf[mfi].box(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // EB2's public convention is negative in fluid.
            phi(i, j, k) = -phi(i, j, k);
        });
    }
    m_levelset = std::move(m_sdf);

    m_ok = true;
}
#endif

template <typename G>
GShopLevel<G>::GShopLevel (IndexSpace const* is, int /*ilev*/, int max_grid_size, int /*ngrow*/,
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

}

#endif
