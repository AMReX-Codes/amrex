#ifndef AMRCORE_H_
#define AMRCORE_H_

#include <memory>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_TagBox.H>

namespace amrex {

#ifdef USE_PARTICLES
class AmrParGDB;
#endif

/**
 * \brief Provide basic functionalities to set up an AMR hierarchy
 *
 * Because it does not own any floating point data, AmrCore used
 * virtual functions to allocate, initialize and delete data.  It also
 * requires the derived class to tag cells for refinement.
 */

class AmrCore
{
public:

    AmrCore ();
    AmrCore (const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord);

    virtual ~AmrCore ();

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
    const Array<IntVect>& refRatio () const { return ref_ratio; }

    //! Should we keep the coarser grids fixed (and not regrid those levels) at all?
    bool useFixedCoarseGrids () const { return use_fixed_coarse_grids; }

    //! Up to what level should we keep the coarser grids fixed (and not regrid those levels)?
    int useFixedUpToLevel () const { return use_fixed_upto_level; }

    const Array<Geometry>& Geom () const { return geom; }
    const Array<DistributionMapping>& DistributionMap () const { return dmap; }
    const Array<BoxArray>& boxArray () const { return grids; }

    const Geometry& Geom (int lev) const { return geom[lev]; }
    const DistributionMapping& DistributionMap (int lev) const { return dmap[lev]; }
    const BoxArray& boxArray (int lev) const { return grids[lev]; }

    Array<Geometry>& Geom () { return geom; }
    Geometry& Geom (int lev) { return geom[lev]; }

    void SetDistributionMap (int lev, const DistributionMapping& dmap_in);
    void SetBoxArray (int lev, const BoxArray& ba_in);

    void ClearDistributionMap (int lev);
    void ClearBoxArray (int lev);

#ifdef USE_PARTICLES
    AmrParGDB* GetParGDB () const { return m_gdb.get(); }
#endif

    //! Return the number of buffer cells in error estimator.
    int nErrorBuf (int lev) const { return n_error_buf[lev]; }

    //! Return the minimum allowable grid efficiency.
    Real gridEff () const { return grid_eff; }

    //! Return the number of cells to define proper nesting 
    int nProper () const { return n_proper; }

    //! Return the blocking factor at level lev
    int blockingFactor (int lev) const { return blocking_factor[lev]; }

    //! Return the largest allowable grid.
    int maxGridSize (int lev) const { return max_grid_size[lev]; }

    bool LevelDefined (int lev);

    //! "Try" to chop up grids so that the number of boxes in the BoxArray is greater than the target_size.
    void ChopGrids (int lev, BoxArray& ba, int target_size) const;

    //! Make a level 0 grids covering the whole domain.  It does NOT install the new grids.
    BoxArray MakeBaseGrids () const;

    /**
    * \brief  Make new grids based on error estimates.  Note that new_finest returns the new finest level.
    * Also note that this function adds at most one more level to the existing levels.
    * The current finest level is stored in this->finest_level.  This function does not
    * change its value, nor does it modifies any BoxArrays stored in this->grids.
    */
    void MakeNewGrids (int lbase, Real time, int& new_finest, Array<BoxArray>& new_grids);

    /**
     * \brief Initialize BoxArray, DistributionMapping and data from scratch.
     * Calling this function requires the derive class implement its own MakeNewLevelFromScratch
     * to allocate and intialize data. 
     * Also note usually one needs to average the fine data down to coarse level after this.
     */
    void InitFromScratch (Real time);

    //! Rebuild levels finer than lbase
    virtual void regrid (int lbase, Real time, bool initial=false);

    static void Initialize ();
    static void Finalize ();

protected:

    //! Make a new level from scratch using provided BoxArray and DistributionMapping.
    //! Only used during initialization.
    virtual void MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) = 0;

    //! Make a new level using provided BoxArray and DistributionMapping and fill with interpolated coarse level data.
    virtual void MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) = 0;

    //! Remake an existing level using provided BoxArray and DistributionMapping and fill with existing fine and coarse data.
    virtual void RemakeLevel (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm) = 0;

    //! Delete level data
    virtual void ClearLevel (int lev) = 0;

    //! Tagging cells for refinement
    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) = 0;

    virtual BoxArray GetAreaNotToTag (int lev) { return BoxArray(); }
    virtual void ManualTagsPlacement (int lev, TagBoxArray& tags, Array<IntVect>& bf_lev) {}

    static void ProjPeriodic (BoxList& bd, const Geometry& geom);
    
    int              verbose;
    int              max_level;    // Maximum allowed level.
    Array<IntVect>   ref_ratio;    // Refinement ratios [0:finest_level-1]

    int              finest_level; // Current finest level.

    Array<int>       n_error_buf;     // Buffer cells around each tagged cell.
    Array<int>       blocking_factor; // Blocking factor in grid generation (by level).
    Array<int>       max_grid_size;   // Maximum allowable grid size (by level).
    Real             grid_eff;        // Grid efficiency.
    int              n_proper;        // # cells required for proper nesting.

    bool use_fixed_coarse_grids;
    int use_fixed_upto_level;

    bool refine_grid_layout;

#ifdef USE_PARTICLES
    std::unique_ptr<AmrParGDB> m_gdb;
#endif

    Array<Geometry>            geom;
    Array<DistributionMapping> dmap;
    Array<BoxArray>            grids;    

private:
    //
    // Disallowed
    //
    AmrCore (const AmrCore& rhs);
    AmrCore& operator= (const AmrCore& rhs);

    void InitAmrCore (int max_level_in, const Array<int>& n_cell_in);
};

}

#endif
