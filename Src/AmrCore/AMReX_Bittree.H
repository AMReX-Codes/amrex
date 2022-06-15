#ifndef BL_BITTREE_H_
#define BL_BITTREE_H_


#include <AMReX_IntVect.H>
#include <AMReX_BoxList.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>


#include <AMReX_BoxArray.H>
#include <Bittree_Amr.h>

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
    static int btRefine( std::shared_ptr<BittreeAmr> mesh, std::vector<int>& btTags, MPI_Comm comm);
    static void btCalculateGrids( std::shared_ptr<BittreeAmr> mesh, int lbase,
                                Real time, int& new_finest,
                                Vector<BoxArray>& new_grids,
                                Vector<DistributionMapping>& new_dm,
                                Vector<IntVect>& max_grid_size);

  // Functions imported from Flash-X
  private:
    static void btCheckRefine( std::shared_ptr<BittreeAmr> mesh, std::vector<int>& btTags, MPI_Comm comm);
    static void btCheckDerefine( std::shared_ptr<BittreeAmr> mesh );
    static std::vector<int> calcNeighIntCoords(unsigned lev, unsigned* lcoord, int* gCell, std::shared_ptr<BittreeAmr> mesh);

  // cache of mesh data
  public:
    // 1 for refine, -1 for derefine, 0 for neither
    static std::vector<int> btTags;
};


}
#endif
