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

 */

class btUnit {
  // Functions corresponding to the three steps of the Bittree version of MakeNewGrids
  public:
    static int  btRefine( std::shared_ptr<BittreeAmr> mesh,
                          std::vector<int>& btTags,
                          MPI_Comm comm);
    static void btCalculateGrids( std::shared_ptr<BittreeAmr> mesh,
                                  int lbase, Real time, int& new_finest,
                                  Vector<BoxArray>& new_grids,
                                  Vector<DistributionMapping>& new_dm,
                                  Vector<IntVect>& max_grid_size);
    static void btCalculateLevel( std::shared_ptr<BittreeAmr> mesh,
                                  int lev, Real time,
                                  BoxArray& ba,
                                  DistributionMapping& dm,
                                  IntVect& max_grid_size);
  // Utils
  public:
    static int getBitid(std::shared_ptr<BittreeAmr> mesh, bool updated,
                        int lev, int idx_on_lev);
    static int getIndex(std::shared_ptr<BittreeAmr> mesh, bool updated,
                        int lev, int bitid);

  // Functions for octree logic, adapted from PARAMESH
  private:
    static void btCheckRefine( std::shared_ptr<BittreeAmr> mesh,
                               std::vector<int>& btTags,
                               MPI_Comm comm);
    static void btCheckDerefine( std::shared_ptr<BittreeAmr> mesh,
                                 std::vector<int>& btTags,
                                 MPI_Comm comm);

    // Utility Functions
    static bool checkNeighborsRefine( std::shared_ptr<BittreeAmr> mesh,
                                MortonTree::Block b);
    static std::vector<int> neighIntCoords(std::shared_ptr<BittreeAmr> mesh,
                                           unsigned lev, unsigned* lcoord,
                                           int* gCell);
};


}
#endif
