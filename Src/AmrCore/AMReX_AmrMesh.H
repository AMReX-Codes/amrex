#ifndef BL_AMRMESH_H_
#define BL_AMRMESH_H_

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_RealBox.H>
#include <AMReX_IntVect.H>
#include <AMReX_Geometry.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <AMReX_TagBox.H>

namespace amrex {

struct AmrInfo {
    int verbose = 0;
    // Maximum allowed level.
    int max_level = 0;
    // Refinement ratios
    Vector<IntVect> ref_ratio {{IntVect(2)}};
    // Blocking factor in grid generation (by level).
    Vector<IntVect> blocking_factor {{IntVect(8)}};
    // Maximum allowable grid size (by level).
    Vector<IntVect> max_grid_size {{IntVect(AMREX_D_PICK(128,128,32))}};
    // Buffer cells around each tagged cell.
    Vector<IntVect> n_error_buf {{IntVect(1)}};
    // Grid efficiency.
    Real grid_eff = 0.7;
    // Cells required for proper nesting.
    int n_proper = 1;
    int use_fixed_upto_level = 0;
    bool use_fixed_coarse_grids = false;
    //chop up grids to have the number of grids no less the number of procs
    bool refine_grid_layout = true;
    bool check_input = true;
    bool use_new_chop = false;
    bool iterate_on_new_grids = true;
};

class AmrMesh
    : protected AmrInfo
{
public:

    friend std::ostream& operator<< (std::ostream& os, AmrMesh const& amr_mesh);

    AmrMesh ();

    AmrMesh (const RealBox* rb, int max_level_in,
             const Vector<int>& n_cell_in, int coord=-1,
             Vector<IntVect> refrat = Vector<IntVect>(),
             const int* is_per = nullptr);

    AmrMesh (const RealBox& rb, int max_level_in,
             const Vector<int>& n_cell_in, int coord,
             Vector<IntVect> const& ref_ratios,
             Array<int,AMREX_SPACEDIM> const& is_per);

    AmrMesh (Geometry const& level_0_geom, AmrInfo const& amr_info);

    AmrMesh (const AmrMesh& rhs) = delete;
    AmrMesh& operator= (const AmrMesh& rhs) = delete;

    AmrMesh (AmrMesh&& rhs) = default;
    AmrMesh& operator= (AmrMesh&& rhs) = default;

    virtual ~AmrMesh () = default;

    int Verbose () const noexcept { return verbose; }

    //! Return the max level
    int maxLevel () const noexcept { return max_level; }

    //! Return the finest level
    int finestLevel () const noexcept { return finest_level; }

    //! Return the refinement ratio for level lev
    IntVect refRatio (int lev) const noexcept { return ref_ratio[lev]; }

    //! Return the maximum refinement ratio in any direction.
    int MaxRefRatio (int lev) const noexcept;

    //! Return refinement ratios between all levels.
    const Vector<IntVect>& refRatio () const noexcept { return ref_ratio; }

    const Vector<Geometry>& Geom () const noexcept { return geom; }
    const Vector<DistributionMapping>& DistributionMap () const noexcept { return dmap; }
    const Vector<BoxArray>& boxArray () const noexcept { return grids; }

    const Geometry& Geom (int lev) const noexcept { return geom[lev]; }
    const DistributionMapping& DistributionMap (int lev) const noexcept { return dmap[lev]; }
    const BoxArray& boxArray (int lev) const noexcept { return grids[lev]; }

    Vector<Geometry> Geom (int a_coarsest_lev, int a_finest_lev) const noexcept {
        Vector<Geometry> r;
        r.reserve(a_finest_lev-a_coarsest_lev+1);
        for (int lev = a_coarsest_lev; lev <= a_finest_lev; ++lev) {
            r.push_back(geom[lev]);
        }
        return r;
    }
    Vector<BoxArray> boxArray (int a_coarsest_lev, int a_finest_lev) const noexcept {
        Vector<BoxArray> r;
        r.reserve(a_finest_lev-a_coarsest_lev+1);
        for (int lev = a_coarsest_lev; lev <= a_finest_lev; ++lev) {
            r.push_back(grids[lev]);
        }
        return r;
    }
    Vector<DistributionMapping> DistributionMap (int a_coarsest_lev, int a_finest_lev) const noexcept {
        Vector<DistributionMapping> r;
        r.reserve(a_finest_lev-a_coarsest_lev+1);
        for (int lev = a_coarsest_lev; lev <= a_finest_lev; ++lev) {
            r.push_back(dmap[lev]);
        }
        return r;
    }

    Vector<Geometry>& Geom () noexcept { return geom; }
    Geometry& Geom (int lev) noexcept { return geom[lev]; }

    void SetMaxGridSize (int new_mgs) noexcept {
        max_grid_size.assign(max_level+1, IntVect{AMREX_D_DECL(new_mgs,new_mgs,new_mgs)});
    }
    void SetMaxGridSize (const IntVect& new_mgs) noexcept {
        max_grid_size.assign(max_level+1, new_mgs);
    }
    void SetMaxGridSize (const Vector<int>& new_mgs) noexcept {
        max_grid_size.resize(max_level+1);
        for (int i = 0; i <= max_level; ++i) {
            max_grid_size[i] = IntVect{AMREX_D_DECL(new_mgs[i],new_mgs[i],new_mgs[i])};
        }
    }
    void SetMaxGridSize (const Vector<IntVect>& new_mgs) noexcept {
        max_grid_size.assign(new_mgs.cbegin(), new_mgs.cbegin()+max_level+1);
    }

    void SetBlockingFactor (int new_bf) noexcept {
        blocking_factor.assign(max_level+1, IntVect{AMREX_D_DECL(new_bf,new_bf,new_bf)});
    }
    void SetBlockingFactor (const IntVect& new_bf) noexcept {
        blocking_factor.assign(max_level+1, new_bf);
    }
    void SetBlockingFactor (const Vector<int>& new_bf) noexcept {
        blocking_factor.resize(max_level+1);
        for (int i = 0; i <= max_level; ++i) {
            blocking_factor[i] = IntVect{AMREX_D_DECL(new_bf[i],new_bf[i],new_bf[i])};
        }
    }
    void SetBlockingFactor (const Vector<IntVect>& new_bf) noexcept {
        blocking_factor.assign(new_bf.cbegin(), new_bf.cend()+max_level+1);
    }

    void SetGridEff (Real eff) noexcept { grid_eff = eff; }
    void SetNProper (int n) noexcept { n_proper = n; }

    //! Set ref_ratio would require rebuiling Geometry objects.

    void SetFinestLevel (int new_finest_level) noexcept { finest_level = new_finest_level; }
    void SetDistributionMap (int lev, const DistributionMapping& dmap_in) noexcept;
    void SetBoxArray (int lev, const BoxArray& ba_in) noexcept;
    void SetGeometry (int lev, const Geometry& geom_in) noexcept;

    //! Given domain box, return AMR level.  Return -1 if there is no match.
    int GetLevel (Box const& domain) noexcept;

    void ClearDistributionMap (int lev) noexcept;
    void ClearBoxArray (int lev) noexcept;

    //! Return the number of buffer cells (as a single integer) in error estimator.
    int nErrorBuf (int lev, int direction = 0) const noexcept { return n_error_buf[lev][direction]; }

    //! Return the number of buffer cells (as an IntVect) in error estimator.
    const IntVect& nErrorBufVect (int lev) const noexcept { return n_error_buf[lev]; }

    //! Return the minimum allowable grid efficiency.
    Real gridEff () const noexcept { return grid_eff; }

    //! Return the number of cells to define proper nesting
    int nProper () const noexcept { return n_proper; }

    //! Return the blocking factor at level lev
    const IntVect& blockingFactor (int lev) const noexcept { return blocking_factor[lev]; }

    //! Return the largest allowable grid.
    const IntVect& maxGridSize (int lev) const noexcept { return max_grid_size[lev]; }

    bool LevelDefined (int lev) noexcept;

    //! Should we keep the coarser grids fixed (and not regrid those levels) at all?
    bool useFixedCoarseGrids () const noexcept { return use_fixed_coarse_grids; }

    //! Up to what level should we keep the coarser grids fixed (and not regrid those levels)?
    int useFixedUpToLevel () const noexcept { return use_fixed_upto_level; }

    //! "Try" to chop up grids so that the number of boxes in the BoxArray is greater than the target_size.
    void ChopGrids (int lev, BoxArray& ba, int target_size) const;

    //! Make a level 0 grids covering the whole domain.  It does NOT install the new grids.
    BoxArray MakeBaseGrids () const;

    /**
    * \brief Make new grids based on error estimates.  This functin
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
    void MakeNewGrids (int lbase, Real time, int& new_finest, Vector<BoxArray>& new_grids);

    //! This function makes new grid for all levels (including level 0).
    void MakeNewGrids (Real time = 0.0);

    //! This function is called by the second version of MakeNewGrids.
    //! Make a new level from scratch using provided BoxArray and DistributionMapping.
    //! Only used during initialization.
    virtual void MakeNewLevelFromScratch (int /*lev*/, Real /*time*/, const BoxArray& /*ba*/, const DistributionMapping& /*dm*/) {}

    //! Tag cells for refinement.  TagBoxArray tags is built on level lev grids.
    virtual void ErrorEst (int /*lev*/, TagBoxArray& /*tags*/, Real /*time*/, int /*ngrow*/) {}

    //! Manually tag.  Note that tags is built on level lev grids coarsened by bf_lev[lev].
    virtual void ManualTagsPlacement (int /*lev*/, TagBoxArray& /*tags*/, const Vector<IntVect>& /*bf_lev*/) {}

    //! Apply some user-defined changes the to base grids.
    //!
    //! This function is only called by MakeNewGrids after computing a box array for the coarsest level
    //! and before calling MakeNewLevelFromScratch.
    //!
    //! For example, use this function if you want to remove covered grids on the coarsest refinement level.
    virtual void PostProcessBaseGrids(BoxArray& /*box_array*/) const {}

    virtual BoxArray GetAreaNotToTag (int /*lev*/) { return BoxArray(); }

    long CountCells (int lev) noexcept;

protected:

    int finest_level;    //!< Current finest level.
    Vector<Geometry>            geom;
    Vector<DistributionMapping> dmap;
    Vector<BoxArray>            grids;

    unsigned int num_setdm = 0;

    void checkInput();

    void SetIterateToFalse () noexcept { iterate_on_new_grids = false; }
    void SetUseNewChop () noexcept { use_new_chop = true; }

private:
    void InitAmrMesh (int max_level_in, const Vector<int>& n_cell_in,
                      Vector<IntVect> refrat = Vector<IntVect>(),
                      const RealBox* rb = nullptr, int coord = -1,
                      const int* is_per = nullptr);

    static void ProjPeriodic (BoxList& bd, const Box& domain,
                              Array<int,AMREX_SPACEDIM> const& is_per);
};

std::ostream& operator<< (std::ostream& os, AmrMesh const& amr_mesh);

}

#endif
