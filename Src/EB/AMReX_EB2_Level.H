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

extern AMREX_EXPORT int max_grid_size;

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

/**
     * \brief Face-centered EB data storage structure.
     *
     * Stores EB data for a specific face direction (x, y, or z).
     * Data is derived from refined (2x finer) cell-centered grids.
     */
    struct FCData {
        FabArray<EBCellFlagFab> m_cellflag_fc;
        MultiFab m_volfrac_fc;
        MultiFab m_centroid_fc;
        MultiFab m_bndryarea_fc;
        MultiFab m_bndrycent_fc;
        MultiFab m_bndrynorm_fc;
        Array<MultiFab,AMREX_SPACEDIM> m_areafrac_fc;
        Array<MultiFab,AMREX_SPACEDIM> m_facecent_fc;
        Array<MultiFab,AMREX_SPACEDIM> m_edgecent_fc;
        bool m_built = false;  //!< Lazy building flag
    };

    //! Check if face-centered data is available for face_dir (0=x, 1=y, 2=z).
    [[nodiscard]] bool hasFCData (int face_dir) const noexcept {
        return face_dir >= 0 && face_dir < AMREX_SPACEDIM &&
               m_fc_data[face_dir] && m_fc_data[face_dir]->m_built;
    }

    //! Get reference to face-centered data for face_dir (0=x, 1=y, 2=z).
    [[nodiscard]] const FCData& getFCData (int face_dir) const {
        AMREX_ASSERT(hasFCData(face_dir));
        return *m_fc_data[face_dir];
    }

    //! Fill face-centered volume fractions for face_dir.
    void fillVolFracFC (MultiFab& vfrac, int face_dir, const Geometry& geom) const;
    //! Fill face-centered area fractions for face_dir.
    void fillAreaFracFC (Array<   MultiFab*,AMREX_SPACEDIM> const& areafrac, int face_dir, const Geometry& geom) const;
    void fillAreaFracFC (Array<MultiCutFab*,AMREX_SPACEDIM> const& areafrac, int face_dir, const Geometry& geom) const;
    //! Fill face-centered centroids for face_dir.
    void fillCentroidFC (MultiFab& centroid, int face_dir, const Geometry& geom) const;
    void fillCentroidFC (MultiCutFab& centroid, int face_dir, const Geometry& geom) const;
    //! Fill face-centered boundary areas for face_dir.
    void fillBndryAreaFC (MultiFab& bndryarea, int face_dir, const Geometry& geom) const;
    void fillBndryAreaFC (MultiCutFab& bndryarea, int face_dir, const Geometry& geom) const;
    //! Fill face-centered boundary centroids for face_dir.
    void fillBndryCentFC (MultiFab& bndrycent, int face_dir, const Geometry& geom) const;
    void fillBndryCentFC (MultiCutFab& bndrycent, int face_dir, const Geometry& geom) const;
    //! Fill face-centered boundary normals for face_dir.
    void fillBndryNormFC (MultiFab& bndrynorm, int face_dir, const Geometry& geom) const;
    void fillBndryNormFC (MultiCutFab& bndrynorm, int face_dir, const Geometry& geom) const;
    //! Fill face-centered EBCellFlags for face_dir.
    void fillEBCellFlagFC (FabArray<EBCellFlagFab>& cellflag, int face_dir, const Geometry& geom) const;
    //! Fill face-centered face centroids for face_dir.
    void fillFaceCentFC (Array<   MultiFab*,AMREX_SPACEDIM> const& facecent, int face_dir, const Geometry& geom) const;
    void fillFaceCentFC (Array<MultiCutFab*,AMREX_SPACEDIM> const& facecent, int face_dir, const Geometry& geom) const;
    //! Fill face-centered edge centroids for face_dir.
    void fillEdgeCentFC (Array<   MultiFab*,AMREX_SPACEDIM> const& edgecent, int face_dir, const Geometry& geom) const;
    void fillEdgeCentFC (Array<MultiCutFab*,AMREX_SPACEDIM> const& edgecent, int face_dir, const Geometry& geom) const;

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
#if (AMREX_SPACEDIM == 3)
    std::map<int,std::unique_ptr<MC::MCFab>> m_marching_cubes;
#endif
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

    //! Face-centered data storage: [0]=x-face, [1]=y-face, [2]=z-face
    Array<std::unique_ptr<FCData>, AMREX_SPACEDIM> m_fc_data;

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
    void define_fine_mvmc (GS const& gshop, const Geometry& geom,
                           int max_grid_size, int ngrow, bool extend_domain_face, int num_crse_opt);
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

    //! Build face-centered EB data for face_dir. If fine_level provided, coarsen from it; otherwise build refined CC transiently.
    void buildFCData (G const& gshop, int face_dir, Level const* fine_level = nullptr);

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
requires (std::same_as<GS,STLtools>)
void
GShopLevel<G>::define_fine_mvmc (GS const& gshop, const Geometry& geom,
                                 int max_grid_size, int ngrow, bool extend_domain_face, int num_crse_opt)
{
    BL_PROFILE("EB2::GShopLevel()-fine-mvmc");

    prepare_grids(gshop, geom, max_grid_size, ngrow, extend_domain_face, num_crse_opt);

    if (m_allregular) { return; }

    MC::Initialize();

    m_sdf.define(amrex::convert(m_grids,IntVect(1)), m_dmap, 1, IntVect(1));
    gshop.fillSignedDistance(m_sdf, m_sdf.nGrowVect(), geom);

    for (MFIter mfi(m_sdf,MFItInfo().DisableDeviceSync()); mfi.isValid(); ++mfi) {
        m_marching_cubes[mfi.index()] = std::make_unique<MC::MCFab>();
    }

    MFItInfo info{};

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
    info.SetDynamic(true);
#pragma omp parallel
#endif
    for (MFIter mfi(m_sdf,info); mfi.isValid(); ++mfi) {
        amrex::MC::marching_cubes(geom, m_sdf[mfi], *m_marching_cubes[mfi.index()]);
    }

    MC::write_stl("test.stl", m_marching_cubes);

    m_ok = false; // xxxxx TODO
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

template <typename G>
void
GShopLevel<G>::buildFCData (G const& gshop, int face_dir, Level const* fine_level)
{
    BL_PROFILE("EB2::GShopLevel::buildFCData");

    AMREX_ASSERT(face_dir >= 0 && face_dir < AMREX_SPACEDIM);

    if (m_allregular) { return; }

    // Allocate FC data structure if needed
    if (!m_fc_data[face_dir]) {
        m_fc_data[face_dir] = std::make_unique<FCData>();
    }

    if (m_fc_data[face_dir]->m_built) { return; }

    // **********************************************************************
    // STEP1: Get or build the one-level-finer CC data
    // **********************************************************************

    Level const* refined_cc_level = nullptr;
    std::unique_ptr<GShopLevel<G>> transient_level;

    if (fine_level != nullptr) {
        // Reuse existing finer CC level from hierarchy
        refined_cc_level = fine_level;
    } else {
        // Finest FC level: build transient refined CC level (2x finer)
        Geometry refined_geom = amrex::refine(m_geom, 2);
        transient_level = std::make_unique<GShopLevel<G>>(
            m_parent, gshop, refined_geom,
            EB2::max_grid_size, m_ngrow[0] > 0 ? m_ngrow[0] : 4, true, 0);

        if (!transient_level->isOK()) {
            amrex::Abort("GShopLevel::buildFCData: failed to build refined level");
        }
        refined_cc_level = transient_level.get();

        // Need FillBoundary here?
    }

    // **********************************************************************
    // STEP2: Check coarsenable
    // **********************************************************************

    const BoxArray& fine_grids = refined_cc_level->m_grids;
    const BoxArray& fine_covered_grids = refined_cc_level->m_covered_grids;
    const int coarse_ratio = 2;
    const int min_width = 8;
    bool coarsenable = fine_grids.coarsenable(coarse_ratio, min_width)
        && (fine_covered_grids.empty() || fine_covered_grids.coarsenable(coarse_ratio));

    // prepareForCoarsening can be called if not coarsenable.
    // if (!coarsenable) {
    //     Level refined_cc_level_2(is, refined_cc_level.m_geom);
    //     refined_cc_level_2.prepareForCoarsening(refined_cc_level, max_grid_size, amrex::scale(m_ngrow,2));
    // }

    // **********************************************************************
    // STEP3: Check multi-cut error
    // **********************************************************************

    auto const& f_levelset = refined_cc_level->m_levelset;
    m_levelset.define(amrex::convert(m_grids,IntVect::TheNodeVector()), m_dmap, 1, 0);

    int mvmc_error = 0;

    if (Gpu::notInLaunchRegion())
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:mvmc_error)
#endif
        for (MFIter mfi(m_levelset,true); mfi.isValid(); ++mfi)
        {
            const Box& ccbx = mfi.tilebox(IntVect::TheCellVector());
            const Box& ndbx = mfi.tilebox();
            auto const& crse = m_levelset.array(mfi);
            auto const& fine = f_levelset.const_array(mfi);

            amrex::LoopConcurrentOnCpu(ndbx,
            [=] (int i, int j, int k) noexcept
            {
                crse(i,j,k) = fine(2*i,2*j,2*k);
            });

            int tile_error = 0;
            amrex::LoopOnCpu(ccbx,
            [=,&tile_error] (int i, int j, int k) noexcept
            {
                int ierror = check_mvmc(i,j,k,fine);
                tile_error = std::max(tile_error,ierror);
            });

            mvmc_error = std::max(mvmc_error, tile_error);
        }
    }
    else
    {
        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(m_levelset); mfi.isValid(); ++mfi)
        {
            const Box& ndbx = mfi.validbox();
            const Box& ccbx = amrex::enclosedCells(ndbx);
            auto const& crse = m_levelset.array(mfi);
            auto const& fine = f_levelset.const_array(mfi);
            reduce_op.eval(ndbx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                crse(i,j,k) = fine(2*i,2*j,2*k);
                int ierror;
                if (ccbx.contains(IntVect{AMREX_D_DECL(i,j,k)})) {
                    ierror = check_mvmc(i,j,k,fine);
                } else {
                    ierror = 0;
                }
                return {ierror};
            });
        }
        ReduceTuple rv = reduce_data.value(reduce_op);
        mvmc_error = amrex::max(0, amrex::get<0>(rv));
    }

    {
        bool b = mvmc_error;
        ParallelDescriptor::ReduceBoolOr(b);
        mvmc_error = b;
    }
    if (mvmc_error) { 
        amrex::Abort("GShopLevel::buildFCData: MVMC error");
    }

    // **********************************************************************
    // STEP3: Define FC data structures
    // **********************************************************************

    // FC MultiFabs
    auto& fc_cellflag = m_fc_data[face_dir]->m_cellflag_fc;
    MultiFab& fc_volfrac = m_fc_data[face_dir]->m_volfrac_fc;
    MultiFab& fc_centroid = m_fc_data[face_dir]->m_centroid_fc;
    MultiFab& fc_bndryarea = m_fc_data[face_dir]->m_bndryarea_fc;
    MultiFab& fc_bndrycent = m_fc_data[face_dir]->m_bndrycent_fc;
    MultiFab& fc_bndrynorm = m_fc_data[face_dir]->m_bndrynorm_fc;
    auto& fc_areafrac = m_fc_data[face_dir]->m_areafrac_fc;
    auto& fc_facecent = m_fc_data[face_dir]->m_facecent_fc;
    auto& fc_edgecent = m_fc_data[face_dir]->m_edgecent_fc;

    // Define FC MultiFabs
    const int ng = 0;  // No ghost cells for FC data
    fc_cellflag.define(m_grids, m_dmap, 1, ng);
    fc_volfrac.define(m_grids, m_dmap, 1, ng);
    fc_centroid.define(m_grids, m_dmap, AMREX_SPACEDIM, ng);
    fc_bndryarea.define(m_grids, m_dmap, 1, ng);
    fc_bndrycent.define(m_grids, m_dmap, AMREX_SPACEDIM, ng);
    fc_bndrynorm.define(m_grids, m_dmap, AMREX_SPACEDIM, ng);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        BoxArray faceba = amrex::convert(m_grids, IntVect::TheDimensionVector(idim));
        fc_areafrac[idim].define(faceba, m_dmap, 1, ng);
        fc_facecent[idim].define(faceba, m_dmap, AMREX_SPACEDIM-1, ng);
        IntVect edge_type{1}; edge_type[idim] = 0;
        fc_edgecent[idim].define(amrex::convert(m_grids, edge_type), m_dmap, 1, ng);
    }

    // FillBoundary ?

    // **********************************************************************
    // STEP4: Coarsen from refined CC to FC with staggered averaging
    // **********************************************************************

    // Fine CC MultiFabs
    auto const& cc_cellflag = refined_cc_level->m_cellflag;
    MultiFab const& cc_volfrac = refined_cc_level->m_volfrac;
    MultiFab const& cc_centroid = refined_cc_level->m_centroid;
    MultiFab const& cc_bndryarea = refined_cc_level->m_bndryarea;
    MultiFab const& cc_bndrycent = refined_cc_level->m_bndrycent;
    MultiFab const& cc_bndrynorm = refined_cc_level->m_bndrynorm;
    auto const& cc_areafrac = refined_cc_level->m_areafrac;
    auto const& cc_facecent = refined_cc_level->m_facecent;
    auto const& cc_edgecent = refined_cc_level->m_edgecent;

    int error = 0;

    if (Gpu::notInLaunchRegion())
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:error)
#endif
        for (MFIter mfi(fc_volfrac,true); mfi.isValid(); ++mfi)
        {
            // Coarse FC arrays
            auto const& cvol = fc_volfrac.array(mfi);
            auto const& ccent = fc_centroid.array(mfi);
            auto const& cba = fc_bndryarea.array(mfi);
            auto const& cbc = fc_bndrycent.array(mfi);
            auto const& cbn = fc_bndrynorm.array(mfi);
            AMREX_D_TERM(auto const& capx = fc_areafrac[0].array(mfi);,
                         auto const& capy = fc_areafrac[1].array(mfi);,
                         auto const& capz = fc_areafrac[2].array(mfi););
            AMREX_D_TERM(auto const& cfcx = fc_facecent[0].array(mfi);,
                         auto const& cfcy = fc_facecent[1].array(mfi);,
                         auto const& cfcz = fc_facecent[2].array(mfi););
            AMREX_D_TERM(auto const& cecx = fc_edgecent[0].array(mfi);,
                         auto const& cecy = fc_edgecent[1].array(mfi);,
                         auto const& cecz = fc_edgecent[2].array(mfi););
            auto const& cflag = fc_cellflag.array(mfi);

            // Fine CC arrays
            auto const& fvol = cc_volfrac.const_array(mfi);
            auto const& fcent = cc_centroid.const_array(mfi);
            auto const& fba = cc_bndryarea.const_array(mfi);
            auto const& fbc = cc_bndrycent.const_array(mfi);
            auto const& fbn = cc_bndrynorm.const_array(mfi);
            AMREX_D_TERM(auto const& fapx = cc_areafrac[0].const_array(mfi);,
                         auto const& fapy = cc_areafrac[1].const_array(mfi);,
                         auto const& fapz = cc_areafrac[2].const_array(mfi););
            AMREX_D_TERM(auto const& ffcx = cc_facecent[0].const_array(mfi);,
                         auto const& ffcy = cc_facecent[1].const_array(mfi);,
                         auto const& ffcz = cc_facecent[2].const_array(mfi););
            AMREX_D_TERM(auto const& fecx = cc_edgecent[0].const_array(mfi);,
                         auto const& fecy = cc_edgecent[1].const_array(mfi);,
                         auto const& fecz = cc_edgecent[2].const_array(mfi););
            auto const& fflag = cc_cellflag.const_array(mfi);

            Box const& bx = mfi.validbox();
            Box const& ndgbx = mfi.grownnodaltilebox(-1,2);

            int tile_error = 0;
            amrex::LoopOnCpu(ndgbx,
            [=,&tile_error] (int i, int j, int k) noexcept
            {
                amrex::ignore_unused(j,k);
                int ierr = coarsen_from_fine_staggered(AMREX_D_DECL(i,j,k), bx, 2, face_dir,
                                             cvol,ccent,cba,cbc,cbn,
                                             AMREX_D_DECL(capx,capy,capz),
                                             AMREX_D_DECL(cfcx,cfcy,cfcz),
                                             AMREX_D_DECL(cecx,cecy,cecz),
                                             cflag,fvol,fcent,fba,fbc,fbn,
                                             AMREX_D_DECL(fapx,fapy,fapz),
                                             AMREX_D_DECL(ffcx,ffcy,ffcz),
                                             AMREX_D_DECL(fecx,fecy,fecz),
                                             fflag);
                tile_error = std::max(tile_error,ierr);
            });

            error = std::max(error,tile_error);
        }
    }
    else
    {
        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(fc_volfrac); mfi.isValid(); ++mfi)
        {
            auto const& cvol = fc_volfrac.array(mfi);
            auto const& ccent = fc_centroid.array(mfi);
            auto const& cba = fc_bndryarea.array(mfi);
            auto const& cbc = fc_bndrycent.array(mfi);
            auto const& cbn = fc_bndrynorm.array(mfi);
            AMREX_D_TERM(auto const& capx = fc_areafrac[0].array(mfi);,
                         auto const& capy = fc_areafrac[1].array(mfi);,
                         auto const& capz = fc_areafrac[2].array(mfi););
            AMREX_D_TERM(auto const& cfcx = fc_facecent[0].array(mfi);,
                         auto const& cfcy = fc_facecent[1].array(mfi);,
                         auto const& cfcz = fc_facecent[2].array(mfi););
            // We use MultiArray4 instead of Array4 here
            // to be below the 2 kB kernel parameter limit for SYCL
            AMREX_D_TERM(auto const& cecx = fc_edgecent[0].arrays();,
                         auto const& cecy = fc_edgecent[1].arrays();,
                         auto const& cecz = fc_edgecent[2].arrays(););
            auto const& cflag = fc_cellflag.array(mfi);

            auto const& fvol = cc_volfrac.const_array(mfi);
            auto const& fcent = cc_centroid.const_array(mfi);
            auto const& fba = cc_bndryarea.const_array(mfi);
            auto const& fbc = cc_bndrycent.const_array(mfi);
            auto const& fbn = cc_bndrynorm.const_array(mfi);
            AMREX_D_TERM(auto const& fapx = cc_areafrac[0].const_array(mfi);,
                         auto const& fapy = cc_areafrac[1].const_array(mfi);,
                         auto const& fapz = cc_areafrac[2].const_array(mfi););
            AMREX_D_TERM(auto const& ffcx = cc_facecent[0].const_array(mfi);,
                         auto const& ffcy = cc_facecent[1].const_array(mfi);,
                         auto const& ffcz = cc_facecent[2].const_array(mfi););
            AMREX_D_TERM(auto const& fecx = cc_edgecent[0].const_arrays();,
                         auto const& fecy = cc_edgecent[1].const_arrays();,
                         auto const& fecz = cc_edgecent[2].const_arrays(););
            auto const& fflag = cc_cellflag.const_array(mfi);

            const int lidx = mfi.LocalIndex();
            Box const& bx = mfi.validbox();
            Box const& gbx = amrex::grow(bx,2);
            Box const& ndgbx = amrex::surroundingNodes(gbx);

            reduce_op.eval(ndgbx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                amrex::ignore_unused(j,k);
                int ierr = coarsen_from_fine_staggered(AMREX_D_DECL(i,j,k), bx, 2, face_dir,
                                             cvol,ccent,cba,cbc,cbn,
                                             AMREX_D_DECL(capx,capy,capz),
                                             AMREX_D_DECL(cfcx,cfcy,cfcz),
                                             AMREX_D_DECL(cecx[lidx],cecy[lidx],cecz[lidx]),
                                             cflag,fvol,fcent,fba,fbc,fbn,
                                             AMREX_D_DECL(fapx,fapy,fapz),
                                             AMREX_D_DECL(ffcx,ffcy,ffcz),
                                             AMREX_D_DECL(fecx[lidx],fecy[lidx],fecz[lidx]),
                                             fflag);
                return {ierr};
            });
        }

        ReduceTuple rv = reduce_data.value(reduce_op);
        error = amrex::max(0, amrex::get<0>(rv));
    }

    {
        bool b = error;
        ParallelDescriptor::ReduceBoolOr(b);
        error = b;
    }

    if (!error) {
        buildCellFlag();
    }

    // return error;

    m_fc_data[face_dir]->m_built = true;
    // If transient_level exists, it's automatically destroyed here
}

}

#endif
