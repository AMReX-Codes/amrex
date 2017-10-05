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

class AmrMesh
{
public:
    AmrMesh ();
    AmrMesh (const RealBox* rb, int max_level_in, const Vector<int>& n_cell_in, int coord=-1,
             std::vector<int> refrat = std::vector<int>());

    AmrMesh (const AmrMesh& rhs) = delete;
    AmrMesh& operator= (const AmrMesh& rhs) = delete;

    virtual ~AmrMesh ();

    int Verbose () const { return verbose; }

    //! Return the max level
    int maxLevel () const { return max_level; }

    //! Return the finest level
    int finestLevel () const { return finest_level; }

    //! Return the refinement ratio for level lev
    IntVect refRatio (int lev) const { return ref_ratio[lev]; }

    //! Return the maximum refinement ratio in any direction.
    int MaxRefRatio (int lev) const;

    //! Return refinement ratios between all levels.
    const Vector<IntVect>& refRatio () const { return ref_ratio; }

    const Vector<Geometry>& Geom () const { return geom; }
    const Vector<DistributionMapping>& DistributionMap () const { return dmap; }
    const Vector<BoxArray>& boxArray () const { return grids; }

    const Geometry& Geom (int lev) const { return geom[lev]; }
    const DistributionMapping& DistributionMap (int lev) const { return dmap[lev]; }
    const BoxArray& boxArray (int lev) const { return grids[lev]; }

    Vector<Geometry>& Geom () { return geom; }
    Geometry& Geom (int lev) { return geom[lev]; }

    void SetMaxGridSize (int new_mgs) {
        max_grid_size.assign(max_level+1, IntVect{D_DECL(new_mgs,new_mgs,new_mgs)});
    }
    void SetMaxGridSize (const IntVect& new_mgs) {
        max_grid_size.assign(max_level+1, new_mgs);
    }
    void SetMaxGridSize (const Vector<int>& new_mgs) { 
        max_grid_size.resize(max_level+1);
        for (int i = 0; i <= max_level; ++i) {
            max_grid_size[i] = IntVect{D_DECL(new_mgs[i],new_mgs[i],new_mgs[i])};
        }
    }
    void SetMaxGridSize (const Vector<IntVect>& new_mgs) { 
        max_grid_size.assign(new_mgs.cbegin(), new_mgs.cbegin()+max_level+1);
    }

    void SetBlockingFactor (int new_bf) {
        blocking_factor.assign(max_level+1, IntVect{D_DECL(new_bf,new_bf,new_bf)});
    }
    void SetBlockingFactor (const IntVect& new_bf) {
        blocking_factor.assign(max_level+1, new_bf);
    }
    void SetBlockingFactor (const Vector<int>& new_bf) {
        blocking_factor.resize(max_level+1);
        for (int i = 0; i <= max_level; ++i) {
            blocking_factor[i] = IntVect{D_DECL(new_bf[i],new_bf[i],new_bf[i])};
        }
    }
    void SetBlockingFactor (const Vector<IntVect>& new_bf) {
        blocking_factor.assign(new_bf.cbegin(), new_bf.cend()+max_level+1);
    }

    void SetGridEff (Real eff) { grid_eff = eff; }
    void SetNProper (int n) { n_proper = n; }

    // Set ref_ratio would require rebuiling Geometry objects.

    void SetFinestLevel (int new_finest_level) { finest_level = new_finest_level; }
    void SetDistributionMap (int lev, const DistributionMapping& dmap_in);
    void SetBoxArray (int lev, const BoxArray& ba_in);

    void ClearDistributionMap (int lev);
    void ClearBoxArray (int lev);

    //! Return the number of buffer cells in error estimator.
    int nErrorBuf (int lev) const { return n_error_buf[lev]; }

    //! Return the minimum allowable grid efficiency.
    Real gridEff () const { return grid_eff; }

    //! Return the number of cells to define proper nesting 
    int nProper () const { return n_proper; }

    //! Return the blocking factor at level lev
    const IntVect& blockingFactor (int lev) const { return blocking_factor[lev]; }

    //! Return the largest allowable grid.
    const IntVect& maxGridSize (int lev) const { return max_grid_size[lev]; }

    bool LevelDefined (int lev);

    //! Should we keep the coarser grids fixed (and not regrid those levels) at all?
    bool useFixedCoarseGrids () const { return use_fixed_coarse_grids; }

    //! Up to what level should we keep the coarser grids fixed (and not regrid those levels)?
    int useFixedUpToLevel () const { return use_fixed_upto_level; }

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
    virtual void MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) {}

    //! Tag cells for refinement.  TagBoxArray tags is built on level lev grids.
    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) {}

    //! Manually tag.  Note that tags is built on level lev grids coarsened by bf_lev[lev].
    virtual void ManualTagsPlacement (int lev, TagBoxArray& tags, const Vector<IntVect>& bf_lev) {}

    virtual BoxArray GetAreaNotToTag (int lev) { return BoxArray(); }

    static void Initialize ();
    static void Finalize ();

protected:
    
    int              verbose;
    int              max_level;    // Maximum allowed level.
    Vector<IntVect>   ref_ratio;    // Refinement ratios [0:finest_level-1]

    int              finest_level; // Current finest level.

    Vector<int>       n_error_buf;     // Buffer cells around each tagged cell.
    Vector<IntVect>   blocking_factor; // Blocking factor in grid generation (by level).
    Vector<IntVect>   max_grid_size;   // Maximum allowable grid size (by level).
    Real             grid_eff;        // Grid efficiency.
    int              n_proper;        // # cells required for proper nesting.

    bool use_fixed_coarse_grids;
    int  use_fixed_upto_level;
    bool refine_grid_layout; // chop up grids to have the number of grids no less the number of procs
    bool check_input;

    Vector<Geometry>            geom;
    Vector<DistributionMapping> dmap;
    Vector<BoxArray>            grids;    

    void checkInput();

private:
  void InitAmrMesh (int max_level_in, const Vector<int>& n_cell_in,
                    std::vector<int> refrat = std::vector<int>());

    static void ProjPeriodic (BoxList& bd, const Geometry& geom);
};

}

#endif

