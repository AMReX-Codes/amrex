#ifndef BL_BITTREE_H_
#define BL_BITTREE_H_

#include <AMReX_IntVect.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>

#include <Bittree_BittreeAmr.h>

using namespace bittree;

namespace amrex {

/* 
Include in Make.local:
BITTREE_PATH = /path/to/bittree/installation
INCLUDE_LOCATIONS += $(BITTREE_PATH)/include
LIBRARY_LOCATIONS += $(BITTREE_PATH)/lib
LIBRARIES += -lbittree

Include in inputs:
amr.use_bittree = true
amr.bt_derefine = true
amr.infer_bt_grids = false

 */

class btUnit {
  // Functions used in AmrMesh
  public:
    static int  btRefine(BittreeAmr* const mesh,
                         std::vector<int>& btTags,
                         MPI_Comm comm);
    static void btCalculateGrids(BittreeAmr* const mesh,
                                 int lbase, Real time, int& new_finest,
                                 Vector<BoxArray>& new_grids,
                                 Vector<IntVect>& max_grid_size);
    static void btCalculateLevel(BittreeAmr* const mesh,
                                 int lev, Real time,
                                 BoxArray& ba,
                                 IntVect& max_grid_size);
  // Utils
  public:
    static int getBitid(BittreeAmr* const mesh, bool updated,
                        int lev, int idx_on_lev);
    static int getIndex(BittreeAmr* const mesh, bool updated,
                        int lev, int bitid);

  // Functions to implement strict octree logic
  private:
    static void btCheckRefine(BittreeAmr* const mesh,
                              std::vector<int>& btTags,
                              MPI_Comm comm);
    static void btCheckDerefine(BittreeAmr* const mesh,
                                std::vector<int>& btTags,
                                MPI_Comm comm);

    // Utility Functions
    static bool checkNeighborsRefine(BittreeAmr* const mesh,
                                     MortonTree::Block b);
    static std::vector<int> neighIntCoords(BittreeAmr* const mesh,
                                           unsigned lev, unsigned* lcoord,
                                           int* gCell);

  public:
    // Represents whether domain has periodic BC in each direction
    // true = Periodic, false = Non-Periodic
    static bool bcPeriodic[AMREX_SPACEDIM];
};


}
#endif
