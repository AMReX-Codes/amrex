#ifndef AMRCORE_H_
#define AMRCORE_H_

#include <memory>

#include <BoxArray.H>
#include <DistributionMapping.H>
#include <Geometry.H>

#ifdef USE_PARTICLES
class AmrParGDB;
#endif

class AmrCore
{
public:

    AmrCore ();
    AmrCore (const RealBox* rb, int max_level_in, Array<int> n_cell_in, int coord);

    ~AmrCore ();

    int maxLevel () const { return max_level; }
    int finestLevel () const { return finest_level; }
    // Refinement ratio for level lev
    IntVect refRatio (int lev) const { return ref_ratio[lev]; }
    // Maximum refinement ratio in any direction.
    int MaxRefRatio (int lev) const;
    // Refinement ratios between all levels.
    const Array<IntVect>& refRatio () const { return ref_ratio; }

    const Geometry& Geom (int lev) const { return geom[lev]; }

    const DistributionMapping& DistributionMap (int lev) const { return dmap[lev]; }
    const BoxArray& boxArray (int lev) const { return grids[lev]; }

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

    static void Initialize ();
    static void Finalize ();

protected:

    int              verbose;
    int              max_level;    // Maximum allowed level.
    Array<IntVect>   ref_ratio;    // Refinement ratios [0:finest_level-1]

    int              finest_level; // Current finest level.

    Array<Geometry>            geom;
    Array<DistributionMapping> dmap;
    Array<BoxArray>            grids;

    Array<int>       n_error_buf;     // Buffer cells around each tagged cell.
    Array<int>       blocking_factor; // Blocking factor in grid generation (by level).
    Array<int>       max_grid_size;   // Maximum allowable grid size (by level).
    Real             grid_eff;        // Grid efficiency.
    int              n_proper;        // # cells required for proper nesting.

#ifdef USE_PARTICLES
    std::unique_ptr<AmrParGDB> m_gdb;
#endif

private:
    //
    // Disallowed
    //
    AmrCore (const AmrCore& rhs);
    AmrCore& operator= (const AmrCore& rhs);

    void InitAmrCore (int max_level_in, Array<int> n_cell_in);
};

#endif
