#ifndef BL_AMRMESH_H_
#define BL_AMRMESH_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_RealBox.H>
#include <AMReX_IntVect.H>
#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <AMReX_TagBox.H>

/**
 * \file AMReX_AmrMesh.H
 * \brief Core mesh hierarchy infrastructure shared by AmrCore and AmrLevel.
 */

#ifdef AMREX_USE_BITTREE
#include <Bittree_BittreeAmr.h>
#endif

namespace amrex {

struct AmrInfo {
    int verbose = 0;
    //! Maximum allowed level.
    int max_level = 0;
    //! Refinement ratios
    Vector<IntVect> ref_ratio {{IntVect(2)}};
    //! Blocking factor in grid generation (by level).
    Vector<IntVect> blocking_factor {{IntVect(8)}};
    //! Maximum allowable grid size (by level).
#if defined(AMREX_USE_GPU)
    Vector<IntVect> max_grid_size {{IntVect(AMREX_D_PICK(128,128,64))}};
#else
    Vector<IntVect> max_grid_size {{IntVect(AMREX_D_PICK(128,128,32))}};
#endif
    //! Buffer cells around each tagged cell.
    Vector<IntVect> n_error_buf {{IntVect(1)}};
    //! Grid efficiency.
    Real grid_eff = static_cast<Real>(0.7);
    //! Cells required for proper nesting.
    int n_proper = 1;
    int use_fixed_upto_level = 0;
    bool use_fixed_coarse_grids = false;

    /**
     * Split the number grids until the number of grids is no less than the
     * number of procs. Overridden by refine_grid_layout_dims.
     */
    bool refine_grid_layout = true;

    /**
     * For each specified dimension, chop in that dimension so that
     * the number of grids is no less than the number of procs.
     */
    IntVect refine_grid_layout_dims = IntVect(1);

    bool check_input = true;
    bool use_new_chop = false;
    bool iterate_on_new_grids = true;
    int max_grid_iterations = 4;
};

/// \ingroup amrex_amrcore
class AmrMesh
    : protected AmrInfo
{
public:

    /**
     * \brief Stream a summary of the mesh hierarchy (level extents, cell counts, etc.).
     *
     * \param os        Destination stream.
     * \param amr_mesh  Mesh to describe.
     * \return The destination stream for chaining.
     */
    friend std::ostream& operator<< (std::ostream& os, AmrMesh const& amr_mesh);

    //! Construct an empty mesh.
    AmrMesh ();

    /**
     * \brief Construct using legacy inputs (same signature as AmrCore).
     *
     * \param rb        Optional domain definition; nullptr uses defaults.
     * \param max_level_in Maximum level count.
     * \param n_cell_in Level-0 cell counts.
     * \param coord     Coordinate system index (defaults to Cartesian if not
     *                  specified in ParmParse parameter geometry.coord_sys).
     * \param refrat    Optional refinement ratios per level.
     * \param is_per    Optional periodicity flags (defaults to non-periodic if
     *                  not specified in ParmParse parameter geometry.is_periodic).
     */
    AmrMesh (const RealBox* rb, int max_level_in,
             const Vector<int>& n_cell_in, int coord=-1,
             Vector<IntVect> refrat = Vector<IntVect>(),
             const int* is_per = nullptr);

    /**
     * \brief Construct with explicit domain and periodicity array.
     *
     * \param rb        Level-0 physical domain.
     * \param max_level_in Maximum level count.
     * \param n_cell_in Level-0 cell counts.
     * \param coord     Coordinate system index.
     * \param a_refrat  Per-level refinement ratios.
     * \param is_per    Array of periodicity flags (0 = nonperiodic, 1 = periodic).
     */
    AmrMesh (const RealBox& rb, int max_level_in,
             const Vector<int>& n_cell_in, int coord,
             Vector<IntVect> const& a_refrat,
             Array<int,AMREX_SPACEDIM> const& is_per);

    /**
     * \brief Construct from a pre-built Geometry and AmrInfo bundle.
     *
     * \param level_0_geom Geometry for level 0.
     * \param amr_info     AMR configuration.
     */
    AmrMesh (Geometry const& level_0_geom, AmrInfo const& amr_info);

    AmrMesh (const AmrMesh& rhs) = delete;
    AmrMesh& operator= (const AmrMesh& rhs) = delete;

    AmrMesh (AmrMesh&& rhs) = default;
    AmrMesh& operator= (AmrMesh&& rhs) = default;

    virtual ~AmrMesh () = default;

    [[nodiscard]] int Verbose () const noexcept { return verbose; }

    //! Return the max level
    [[nodiscard]] int maxLevel () const noexcept { return max_level; }

    //! Return the finest level
    [[nodiscard]] int finestLevel () const noexcept { return finest_level; }

    //! Return the refinement ratio for level lev
    [[nodiscard]] IntVect refRatio (int lev) const noexcept { return ref_ratio[lev]; }

    //! Return the maximum refinement ratio in any direction.
    [[nodiscard]] int MaxRefRatio (int lev) const noexcept;

    //! Return refinement ratios between all levels.
    [[nodiscard]] const Vector<IntVect>& refRatio () const noexcept { return ref_ratio; }

    [[nodiscard]] const Vector<Geometry>& Geom () const noexcept { return geom; }
    [[nodiscard]] const Vector<DistributionMapping>& DistributionMap () const noexcept { return dmap; }
    [[nodiscard]] const Vector<BoxArray>& boxArray () const noexcept { return grids; }

    [[nodiscard]] const Geometry& Geom (int lev) const noexcept { return geom[lev]; }
    [[nodiscard]] const DistributionMapping& DistributionMap (int lev) const noexcept { return dmap[lev]; }
    [[nodiscard]] const BoxArray& boxArray (int lev) const noexcept { return grids[lev]; }

    [[nodiscard]] Vector<Geometry> Geom (int a_coarsest_lev, int a_finest_lev) const noexcept {
        Vector<Geometry> r;
        r.reserve(a_finest_lev-a_coarsest_lev+1);
        for (int lev = a_coarsest_lev; lev <= a_finest_lev; ++lev) {
            r.push_back(geom[lev]);
        }
        return r;
    }
    [[nodiscard]] Vector<BoxArray> boxArray (int a_coarsest_lev, int a_finest_lev) const noexcept {
        Vector<BoxArray> r;
        r.reserve(a_finest_lev-a_coarsest_lev+1);
        for (int lev = a_coarsest_lev; lev <= a_finest_lev; ++lev) {
            r.push_back(grids[lev]);
        }
        return r;
    }
    [[nodiscard]] Vector<DistributionMapping> DistributionMap (int a_coarsest_lev, int a_finest_lev) const noexcept {
        Vector<DistributionMapping> r;
        r.reserve(a_finest_lev-a_coarsest_lev+1);
        for (int lev = a_coarsest_lev; lev <= a_finest_lev; ++lev) {
            r.push_back(dmap[lev]);
        }
        return r;
    }

    Vector<Geometry>& Geom () noexcept { return geom; }
    Geometry& Geom (int lev) noexcept { return geom[lev]; }

    /**
     * \brief Override the max grid size constraint for all levels (scalar version).
     *
     * \param new_mgs Maximum length in each direction.
     */
    void SetMaxGridSize (int new_mgs) noexcept {
        max_grid_size.assign(max_level+1, IntVect{AMREX_D_DECL(new_mgs,new_mgs,new_mgs)});
    }
    //! Override the max grid size constraint using explicit IntVect \p new_mgs.
    void SetMaxGridSize (const IntVect& new_mgs) noexcept {
        max_grid_size.assign(max_level+1, new_mgs);
    }
    //! Override per-level max grid sizes using scalar-per-level vector \p new_mgs.
    void SetMaxGridSize (const Vector<int>& new_mgs) noexcept {
        max_grid_size.resize(max_level+1);
        for (int i = 0; i <= max_level; ++i) {
            max_grid_size[i] = IntVect{AMREX_D_DECL(new_mgs[i],new_mgs[i],new_mgs[i])};
        }
    }
    //! Override per-level max grid sizes using explicit IntVect-per-level vector \p new_mgs.
    void SetMaxGridSize (const Vector<IntVect>& new_mgs) noexcept {
        max_grid_size.assign(new_mgs.cbegin(), new_mgs.cbegin()+max_level+1);
    }

    /**
     * \brief Override the blocking factor constraint for all levels (scalar version).
     *
     * \param new_bf Blocking factor to enforce in every direction.
     */
    void SetBlockingFactor (int new_bf) noexcept {
        blocking_factor.assign(max_level+1, IntVect{AMREX_D_DECL(new_bf,new_bf,new_bf)});
    }
    //! Override blocking factor using explicit IntVect \p new_bf.
    void SetBlockingFactor (const IntVect& new_bf) noexcept {
        blocking_factor.assign(max_level+1, new_bf);
    }
    //! Override per-level blocking factors using scalar-per-level vector \p new_bf.
    void SetBlockingFactor (const Vector<int>& new_bf) noexcept {
        blocking_factor.resize(max_level+1);
        for (int i = 0; i <= max_level; ++i) {
            blocking_factor[i] = IntVect{AMREX_D_DECL(new_bf[i],new_bf[i],new_bf[i])};
        }
    }
    //! Override per-level blocking factors using explicit IntVect vector \p new_bf.
    void SetBlockingFactor (const Vector<IntVect>& new_bf) noexcept {
        AMREX_ALWAYS_ASSERT(new_bf.size() >= (max_level+1));
        blocking_factor.assign(new_bf.cbegin(), new_bf.cbegin()+max_level+1);
    }

    //! Set the minimum acceptable grid efficiency to \p eff.
    void SetGridEff (Real eff) noexcept { grid_eff = eff; }

    //! Set the number of cell blocks required for proper nesting to \p n.
    void SetNProper (int n) noexcept { n_proper = n; }

    // There is no set for ref_ratio because it would require rebuilding Geometry objects.

    //! Update the recorded finest level index to \p new_finest_level (without reallocating).
    void SetFinestLevel (int new_finest_level) noexcept { finest_level = new_finest_level; }
    /**
     * \brief Install a DistributionMapping for level \p lev.
     *
     * \param lev     Level index being updated.
     * \param dmap_in Distribution mapping to store.
     */
    void SetDistributionMap (int lev, const DistributionMapping& dmap_in) noexcept;

    /**
     * \brief Install a BoxArray for level \p lev.
     *
     * \param lev   Level index.
     * \param ba_in Grid layout to store.
     */
    void SetBoxArray (int lev, const BoxArray& ba_in) noexcept;

    /**
     * \brief Install a Geometry object for level \p lev.
     *
     * \param lev      Level index.
     * \param geom_in  Geometry to store (copied).
     */
    void SetGeometry (int lev, const Geometry& geom_in) noexcept;

    //! Given domain box, return AMR level.  Return -1 if there is no match.
    [[nodiscard]] int GetLevel (Box const& domain) const noexcept;

    //! Clear the stored DistributionMapping for level \p lev.
    void ClearDistributionMap (int lev) noexcept;
    //! Clear the stored BoxArray for level \p lev.
    void ClearBoxArray (int lev) noexcept;

    //! Return the number of buffer cells (as a single integer) in error estimator.
    [[nodiscard]] int nErrorBuf (int lev, int direction = 0) const noexcept { return n_error_buf[lev][direction]; }

    //! Return the number of buffer cells (as an IntVect) in error estimator.
    [[nodiscard]] const IntVect& nErrorBufVect (int lev) const noexcept { return n_error_buf[lev]; }

    //! Return the minimum allowable grid efficiency.
    [[nodiscard]] Real gridEff () const noexcept { return grid_eff; }

    //! Return the number of cell blocks to define proper nesting
    [[nodiscard]] int nProper () const noexcept { return n_proper; }

    //! Return the blocking factor at level lev
    [[nodiscard]] const IntVect& blockingFactor (int lev) const noexcept { return blocking_factor[lev]; }

    //! Return the largest allowable grid.
    [[nodiscard]] const IntVect& maxGridSize (int lev) const noexcept { return max_grid_size[lev]; }

    //! Return true if valid grids/dmaps exist for level \p lev.
    [[nodiscard]] bool LevelDefined (int lev) const noexcept;

    //! Should we keep the coarser grids fixed (and not regrid those levels) at all?
    [[nodiscard]] bool useFixedCoarseGrids () const noexcept { return use_fixed_coarse_grids; }

    //! Up to what level should we keep the coarser grids fixed (and not regrid those levels)?
    [[nodiscard]] int useFixedUpToLevel () const noexcept { return use_fixed_upto_level; }

    /**
     * \brief Try to chop grids on level \p lev until BoxArray \p ba has at least \p target_size boxes.
     *
     * \param lev         Level whose blocking and efficiency rules apply.
     * \param ba          BoxArray to refine in place.
     * \param target_size Desired minimum number of boxes.
     */
    void ChopGrids (int lev, BoxArray& ba, int target_size) const;

    //! Make a level-0 BoxArray covering the whole domain (without installing it).
    [[nodiscard]] BoxArray MakeBaseGrids () const;

    /**
    * \brief Make new grids based on error estimates.  This function
    * expects that valid BoxArrays exist in this->grids from level
    * lbase to level this->finest_level (the current finest level).
    * new_grids.  On return, the new finest level is stored in
    * new_finest, and the new grids are stored in new_grids from Array
    * element lbase+1 to new_finest_level (unless fixed grids are
    * used).  Note that this function adds at most one more level to
    * the existing levels, and it may remove all levels above the base
    * level.  This function does not change the value of
    * this->finest_level, nor does it modifies any BoxArrays stored in
    * this->grids.  It also does not modify new_grids's elements
    * outside the range [lbase+1,new_finest_level].
    */
    /**
     * \brief Build new grids above level \p lbase using current error tags.
     *
     * \param lbase      Coarsest level to keep fixed.
     * \param time       Simulation time used when tagging.
     * \param new_finest Output: resulting finest level index.
     * \param new_grids  Output: new BoxArrays for levels > lbase (vector sized finest+1).
     */
    void MakeNewGrids (int lbase, Real time, int& new_finest, Vector<BoxArray>& new_grids);

    //! Convenience overload that rebuilds every level (including level 0) at time \p time.
    void MakeNewGrids (Real time = 0.0);

    /**
     * \brief Hook invoked by MakeNewGrids() to build a fresh level.
     *
     * \param lev Level index being created.
     * \param time Simulation time used for initialization.
     * \param ba  BoxArray describing the level's grids.
     * \param dm  Distribution mapping describing ownership.
     */
    virtual void MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) {
        amrex::ignore_unused(lev, time, ba, dm);
    }

    /**
     * \brief Hook for tagging cells on level \p lev.
     *
     * \param lev   Level index.
     * \param tags  Tag boxes (built on level \p lev).
     * \param time  Simulation time.
     * \param ngrow Number of ghost layers filled in \p tags.
     */
    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) {
        amrex::ignore_unused(lev, tags, time, ngrow);
    }

    /**
     * \brief Manual tagging hook. Tags are defined on level \p lev grids coarsened by \p bf_lev[lev].
     *
     * \param lev    Level index.
     * \param tags   Tag boxes to update.
     * \param bf_lev Buffer sizes per level.
     */
    virtual void ManualTagsPlacement (int lev, TagBoxArray& tags, const Vector<IntVect>& bf_lev) {
        amrex::ignore_unused(lev, tags, bf_lev);
    }

    /**
     * \brief Allow derived classes to edit the base-level grids before installing them.
     *
     * For example, use this function if you want to remove covered grids on the coarsest refinement level.
     *
     * This is called by MakeNewGrids() after computing the coarsest-level BoxArray and
     * before invoking MakeNewLevelFromScratch().
     *
     * \param box_array Base-level boxes (modifiable in place).
     */
    virtual void PostProcessBaseGrids(BoxArray& box_array) const { amrex::ignore_unused(box_array); }

    /**
     * \brief Return a BoxArray that should never be tagged on level \p lev (e.g., buffer zones).
     */
    [[nodiscard]] virtual BoxArray GetAreaNotToTag (int lev) { (void)lev; return BoxArray(); }

    //! Count the total number of cells at level \p lev.
    [[nodiscard]] Long CountCells (int lev) const noexcept;

    //! Build a distribution map for level \p lev grid layout \p ba (override to customize load balancing).
    [[nodiscard]] virtual DistributionMapping MakeDistributionMap (int lev, BoxArray const& ba);

protected:

    int finest_level;    //!< Current finest level.
    Vector<Geometry>            geom;
    Vector<DistributionMapping> dmap;
    Vector<BoxArray>            grids;

#ifdef AMREX_USE_BITTREE
    bool use_bittree = false;
    std::unique_ptr<bittree::BittreeAmr> btmesh;
#endif

    unsigned int num_setdm = 0;
    unsigned int num_setba = 0;

    //! Validate AmrInfo inputs (throws if inconsistent).
    void checkInput();

    //! Disable the iterate-on-new-grids loop for the next MakeNewGrids call.
    void SetIterateToFalse () noexcept { iterate_on_new_grids = false; }

    //! Force MakeNewGrids to use the newer chopping algorithm.
    void SetUseNewChop () noexcept { use_new_chop = true; }

private:
    void InitAmrMesh (int max_level_in, const Vector<int>& n_cell_in,
                      Vector<IntVect> refrat = Vector<IntVect>(),
                      const RealBox* rb = nullptr, int coord = -1,
                      const int* is_per = nullptr);

    static void ProjPeriodic (BoxList& blout, const Box& domain,
                              Array<int,AMREX_SPACEDIM> const& is_per);
};

//! Stream helper; forwards to the friend declared inside AmrMesh.
std::ostream& operator<< (std::ostream& os, AmrMesh const& amr_mesh);

}

#endif
