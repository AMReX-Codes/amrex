#ifndef AMRCORE_H_
#define AMRCORE_H_

#include <memory>

#include <BoxArray.H>
#include <DistributionMapping.H>
#include <Geometry.H>
#include <TagBox.H>

#ifdef USE_PARTICLES
class AmrParGDB;
#endif

class AmrCore
{
public:

    AmrCore ();
    AmrCore (const RealBox* rb, int max_level_in, const Array<int>& n_cell_in, int coord);

    virtual ~AmrCore ();

    int Verbose () const { return verbose; }

    int maxLevel () const { return max_level; }
    int finestLevel () const { return finest_level; }
    // Refinement ratio for level lev
    IntVect refRatio (int lev) const { return ref_ratio[lev]; }
    // Maximum refinement ratio in any direction.
    int MaxRefRatio (int lev) const;
    // Refinement ratios between all levels.
    const Array<IntVect>& refRatio () const { return ref_ratio; }
    //
    bool useFixedCoarseGrids () const { return use_fixed_coarse_grids; }
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

    // Number of buffer cells in error estimator.
    int nErrorBuf (int lev) const { return n_error_buf[lev]; }
    // Minimum allowable grid efficiency.
    Real gridEff () const { return grid_eff; }
    // Number of cells to define proper nesting 
    int nProper () const { return n_proper; }
    int blockingFactor (int lev) const { return blocking_factor[lev]; }
    // Largest allowable grid.
    int maxGridSize (int lev) const { return max_grid_size[lev]; }

    bool LevelDefined (int lev);

    // "Try" to chop up grids so that the number of boxes in the BoxArray is greater than
    // the target_size.
    void ChopGrids (int lev, BoxArray& ba, int target_size) const;

    // Make a level 0 grids covering the whole domain.  It does NOT install the new grids.
    BoxArray MakeBaseGrids () const;

    // Make new grids based on error estimates.  Note that new_finest returns the new finest level.
    // Also note that this function adds at most one more level to the existing levels.
    // The current finest level is stored in this->finest_level.  This function does not
    // change its value, nor does it modifies any BoxArrays stored in this->grids.
    void MakeNewGrids (int lbase, Real time, int& new_finest, Array<BoxArray>& new_grids);

    static void Initialize ();
    static void Finalize ();

protected:

    static void ProjPeriodic (BoxList& bd, const Geometry& geom);
    
    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) = 0;
    virtual BoxArray GetAreaNotToTag (int lev) { return BoxArray(); }
    virtual void ManualTagsPlacement (int lev, TagBoxArray& tags, Array<IntVect>& bf_lev) {}

    int              verbose;
    int              max_level;    // Maximum allowed level.
    Array<IntVect>   ref_ratio;    // Refinement ratios [0:finest_level-1]

    int              finest_level; // Current finest level.

    Array<int>       n_error_buf;     // Buffer cells around each tagged cell.
    Array<int>       blocking_factor; // Blocking factor in grid generation (by level).
    Array<int>       max_grid_size;   // Maximum allowable grid size (by level).
    Real             grid_eff;        // Grid efficiency.
    int              n_proper;        // # cells required for proper nesting.

    // Should we keep the coarser grids fixed (and not regrid those levels) at all?
    bool use_fixed_coarse_grids;
    // Up to what level should we keep the coarser grids fixed (and not regrid those levels)?
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

#endif
