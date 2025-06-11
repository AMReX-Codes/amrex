#ifndef BL_BITTREE_H_
#define BL_BITTREE_H_

#include <AMReX_IntVect.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <Bittree_BittreeAmr.h>

namespace amrex {

/*
Include in Make.local:
BITTREE_PATH = /path/to/bittree/installation
INCLUDE_LOCATIONS += $(BITTREE_PATH)/include
LIBRARY_LOCATIONS += $(BITTREE_PATH)/lib
LIBRARIES += -lbittree

Include in inputs:
amr.use_bittree = true
*/
class btUnit {
  // Functions used in AmrMesh
  public:
    static int  btRefine (bittree::BittreeAmr* mesh,
                          std::vector<int>& btTags,
                          int max_crse, int lbase,
                          Vector<BoxArray>& grids, Vector<DistributionMapping>& dmap,
                          MPI_Comm comm);
    static void btCalculateGrids (bittree::BittreeAmr* mesh,
                                  int lbase,
                                  int& new_finest,
                                  Vector<BoxArray>& new_grids,
                                  Vector<IntVect> const& max_grid_size);
    static void btCalculateLevel (bittree::BittreeAmr* mesh,
                                  int lev,
                                  BoxArray& ba,
                                  IntVect const& max_grid_size);
  // Utils
    static int getBitid (bittree::BittreeAmr* mesh, bool updated,
                         int lev, int idx_on_lev);
    static int getIndex (bittree::BittreeAmr* mesh, bool updated,
                         int lev, int bitid);

  // Functions to implement strict octree logic
  private:
    static void btCheckRefine (bittree::BittreeAmr* mesh,
                               std::vector<int>& btTags,
                               int max_crse, int lbase,
                               Vector<BoxArray>& grids, Vector<DistributionMapping>& dmap,
                               MPI_Comm comm);

    static void btCheckDerefine (bittree::BittreeAmr* mesh,
                                 std::vector<int>& btTags,
                                 int max_crse, int lbase,
                                 Vector<BoxArray>& grids, Vector<DistributionMapping>& dmap,
                                 MPI_Comm comm);

    // Utility Functions
    static bool checkNeighborsRefine (bittree::BittreeAmr* mesh,
                                      bittree::MortonTree::Block b);
    static std::vector<int> neighIntCoords (bittree::BittreeAmr* mesh,
                                            unsigned lev, unsigned const* lcoord,
                                            int const* gCell);

  public:
    // Represents whether domain has periodic BC in each direction
    // true = Periodic, false = Non-Periodic
    static bool bcPeriodic[AMREX_SPACEDIM];
};


}
#endif
