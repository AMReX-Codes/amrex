#ifndef BL_BITTREE_H_
#define BL_BITTREE_H_

#include <AMReX_BoxArray.H>
#include <Bittree_Amr.h>

using namespace bittree;

class btUnit {
  // Functions
  public:
    static void btErrorEst( std::shared_ptr<BittreeAmr> mesh );
    static void btRefine( std::shared_ptr<BittreeAmr> mesh );
    static void btMakeNewGrids( std::shared_ptr<BittreeAmr> mesh, int lbase,
                                Real time, int& new_finest,
                                Vector<BoxArray>& new_grids,
                                Vector<DistributionMapping>& new_dm);

  private:
    static void btCheckRefine( std::shared_ptr<BittreeAmr> mesh );
    static void btCheckDerefine( std::shared_ptr<BittreeAmr> mesh );
    static std::vector<int> calcNeighIntCoords(unsigned lev, unsigned* lcoord, int* gCell, std::shared_ptr<BittreeAmr> mesh);

  // cache of mesh data
  private:
    static std::vector<bool> refine;
    static std::vector<bool> derefine;
    static std::vector<std::vector<unsigned>> lcoord;
    static std::vector<unsigned> lev;
    static std::vector<bool> is_par;
    static std::vector<unsigned> bitid;
    static std::vector<std::vector<double>> error;
    static std::vector<std::vector<double>> error_par;

};


#endif
