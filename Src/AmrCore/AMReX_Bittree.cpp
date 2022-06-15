#include <AMReX_Bittree.H>

#include <iostream>
#include <functional>

#define NVARS 2

namespace amrex {

/*
NOTE: Bittree object is created in AmrMesh::MakeNewGrids (Real time)
   with
      `mesh = std::make_shared<BittreeAmr>(top,includes);`

The functions here are called in the BT version of MakeNewGrids which has three steps:
    1. Error Estimation and tagging - btTagging

    2. Bitree's actual bitmap generated/updated - btRefine

    3. AMReX updates grids based on bitree - btCalculateGrids
*/


/** New Bittree mesh is generated.
  *
  * This makes use of BT library functions and as well as routines adapted
  * from Flash-X that enforce Octree logic.
  */
int btUnit::btRefine( std::shared_ptr<BittreeAmr> mesh, std::vector<int>& btTags, MPI_Comm comm) {
    // Tree before refinement. With only one rank, lnblocks = nblocks.
    auto tree0 = mesh->getTree();
    unsigned id0 = tree0->level_id0(0);
    unsigned id1 = tree0->id_upper_bound();

//--Mark leaves to be refined
    for( unsigned id = id0; id < id1; ++id) {
      if (btTags[id]==1) {
        mesh->refine_mark(id, true);
      }
    }
    mesh->refine_reduce(comm);

    btCheckRefine(mesh,btTags,comm);

//--Mark leaves for derefinement (on BT, parent is marked for nodetype change)
    for( unsigned id = id0; id < id1; ++id) {
      if (btTags[id]==-1) {
        auto b = tree0->locate(id);
        unsigned parCoord[AMREX_SPACEDIM];
        for(unsigned d=0; d<AMREX_SPACEDIM; ++d) parCoord[d] = b.coord[d]/2;
        auto p = tree0->identify(b.level-1,parCoord);

        mesh->refine_mark(p.id, true);
      }
    }
    mesh->refine_reduce(comm);

    //btCheckDerefine(mesh);

//--Generate updated Bittree
    mesh->refine_update();

    //TODO sort (distribute over procs)

    // return delta count
    return static_cast<int>( mesh->delta_count() );
}

/** Creates new box arrays to match the new Bittree mesh.
  * TODO: also calculate an appropriate DistributionMapping.
  */
void btUnit::btCalculateGrids(std::shared_ptr<BittreeAmr> mesh, int lbase,
                            Real time,int& new_finest,
                            Vector<BoxArray>& new_grids,
                            Vector<DistributionMapping>& new_dm,
                            Vector<IntVect>& max_grid_size) {
    auto tree1 = mesh->getTree(true);
    int nlevs = tree1->levels();
    new_finest = nlevs - 1;

//--Calculate the new grid layout and distribution map based on Bittree
    for(int lev=lbase; lev<nlevs; ++lev) {
      //Bittree has its own indices for blocks which I call bitid; get
      //the range of bitids for the level being made. Bitid range is
      //contiguous for each level.
      int id0 = tree1->level_id0(lev);
      int id1 = tree1->level_id1(lev);
      int nblocks = tree1->level_blocks(lev);

      BoxList bl;
      Vector<int> pmap;

      for(int i=id0; i<id1; ++i) {
        //Get coordinates and morton index.
        auto b = tree1->locate(i);

        if(b.level != lev) {
            std::string msg = "Error identifying block in btCalculateGrids";
            //TODO throw error
        }

        IntVect coordVec{AMREX_D_DECL(static_cast<int>(b.coord[0]),
                                      static_cast<int>(b.coord[1]),
                                      static_cast<int>(b.coord[2]))};
        IntVect lo = max_grid_size[lev]*coordVec;
        IntVect hi = max_grid_size[lev]*(coordVec+1) - 1;
        bl.push_back( Box{lo,hi} );

        //TODO Calculate the processor based on position in the global Morton curve.
        int proc = 0;
        pmap.push_back(proc);

      }

      new_grids[lev].define(bl);
      new_dm[lev].define(pmap);

    }
}



//---------------------------------------------------------------------
// Local Routines
//---------------------------------------------------------------------

/** Implements the logic which ensures the generated Bittree adheres
  * to a strict octree structure with no more than one level difference
  * between surrounding leaf blocks.
  */
void btUnit::btCheckRefine( std::shared_ptr<BittreeAmr> mesh, std::vector<int>& btTags, MPI_Comm comm ) {
    int cIdx[3];
    int nIdx[3];

    // Tree before refinement. With only one rank, lnblocks = nblocks.
    // TODO only check local blocks?
    auto tree0 = mesh->getTree();
    unsigned lnblocks = tree0->blocks();
    unsigned id0 = tree0->level_id0(0);
    unsigned id1 = tree0->id_upper_bound();
    std::vector<int> ref_test( tree0->id_upper_bound());


    bool repeat = true;

//--Repeat is left true if another round is needed
    while(repeat) {
      std::fill(ref_test.begin(),ref_test.end(),0);

//----Check adjacent children of neighbors of leaf blocks to see if they
//----are marked for refinement.

      for( int id = id0; id < id1; ++id) {
        auto b = tree0->locate(id);
        if( !b.is_parent && btTags[id]!=1 ) {

//--------Loop over neighbors
          for(nIdx[2]= -K3D; nIdx[2]<= K3D; ++nIdx[2]) {
          for(nIdx[1]= -K2D; nIdx[1]<= K2D; ++nIdx[1]) {
          for(nIdx[0]= -K1D; nIdx[0]<= K1D; ++nIdx[0]) {
            std::vector<int> neighCoord = calcNeighIntCoords(b.level, b.coord, nIdx, mesh);

            if(AMREX_D_TERM(neighCoord[0]<0, || neighCoord[1]<0, || neighCoord[2]<0 )) {
                continue;
            }

            unsigned neighCoord_u[AMREX_SPACEDIM];
            for(unsigned d=0; d<AMREX_SPACEDIM; ++d) neighCoord_u[d] = static_cast<unsigned>(neighCoord[d]);
            auto n = tree0->identify(b.level, neighCoord_u);

            if (b.level == n.level) { //neighbor exists
              if(n.is_parent) {
                // Check its children
                unsigned childCoord_u[AMREX_SPACEDIM];

                for(cIdx[2]= 0; cIdx[2]<= K3D; ++cIdx[2]) {
                for(cIdx[1]= 0; cIdx[1]<= K2D; ++cIdx[1]) {
                for(cIdx[0]= 0; cIdx[0]<= K1D; ++cIdx[0]) {
                  for(unsigned d=0; d<AMREX_SPACEDIM; ++d) {
                    childCoord_u[d] = neighCoord_u[d]*2 + cIdx[d];
                  }
                  auto c = tree0->identify(n.level+1, childCoord_u);
                  if(c.is_parent or btTags[c.id]==1) {
                    ref_test[id] = 1;
                    //goto endloops;
                  }
                }}}

              }
            }
          }}} // neighbor loop
//endloops:
          1+1;
        } // if leaf block not already marked
      } // local block loop

      repeat = false;

      for( unsigned id = id0; id < id1; ++id) {
          if( ref_test[id] && btTags[id]!=1 ) {
              repeat = true;
              btTags[id] = 1;
              mesh->refine_mark(id,true);
          }
      }

      // TODO Check all processors to see if a repeat is necessary
      // if(repeat) mesh->refine_reduce(comm);

    } // while repeat
}


/** Implements the logic which ensures the generated Bittree adheres
  * to a strict octree structure with no more than one level difference
  * between surrounding leaf blocks.
  */
//void btUnit::btCheckDerefine( std::shared_ptr<BittreeAmr> mesh ) {
//    // Tree before refinement. With only one rank, lnblocks = nblocks.
//    auto tree = mesh->getTree();
//    unsigned lnblocks = tree->blocks();
//
//
//    bool repeat = true;
//
////--Repeat is left true if another round is needed
//    while(repeat) {
//      //Turn OFF deref_test if block can't be derefined
//      std::vector<bool> deref_test = derefine;
//
////----Check neighbors - if any neighbor is either a parent 
////----or marked for refinement, do not derefine. 
//      for( unsigned lb = 0; lb<lnblocks; ++lb) {
//        if( derefine[lb] ) {
//
////--------Loop over neighbors
//          for(int lk= -K3D; lk<= K3D; ++lk) {
//          for(int lj= -K2D; lj<= K2D; ++lj) {
//          for(int li= -K1D; li<= K1D; ++li) {
//            std::vector<int> gCell{li,lj,lk};
//
//            auto neighCoord = calcNeighIntCoords(lev[lb], lcoord[lb].data(), gCell.data(), mesh);
//            // TODO continue if any(neighCoords<0)
//            unsigned neighCoord_u[AMREX_SPACEDIM];
//            for(unsigned d=0; d<AMREX_SPACEDIM; ++d) neighCoord_u[d] = static_cast<unsigned>(neighCoord[d]);
//            auto b = tree->identify(lev[lb], neighCoord_u);
//
//            if (b.level == lev[lb]) { //neighbor exists
//              bool neighRefine = mesh->check_refine_bit(b.id);
//              if(neighRefine) {
//                deref_test[lb] = false;
//                // TODO break neighbor loop
//              }
//            }
//
//            /* TODO
//            if( b.is_parent) {
//              if( any(lcoord[lb]/2 != neighCoord/2) ) {
//                bool deref_mark = mesh->check_refine_bit(b.id);
//                if(!deref_mark) {
//                  deref_test[lb] = false;
//                  // break
//                }
//              } else {
//                  deref_test[lb] = false;
//                  // break
//              }
//            } // neighbor is parent */
//          }}} // neighbor loop
//        } // if leaf block not already marked
//      } // local block loop
//
//      repeat = false;
//
//      for(unsigned lb=0; lb<lnblocks; ++lb) {
//          if( !deref_test[lb] && derefine[lb] ) {
//              repeat = true;
//              derefine[lb] = false;
//
//              // Unmark for derefinement
//          }
//      }
//
//      //mesh->refine_reduce_and(meshComm);
//
////----If any blocks are still marked for derefinement, check to make
////----sure their parents are still marked on bittree. This ensures blocks
////----only derefine if ALL siblings are marked for derefine.
//      /*
//      for(unsigned lb=0; lb<lnblocks; ++lb) {
//        if( derefine[lb] ) {
//          auto p = tree->(lev[lb]-1,lcoord/2);
//          bool deref_mark = mesh->check_refine_bit(p.id);
//          if(!deref_mark) {
//            repeat = true;
//            derefine[lb] = false;
//          }
//        }
//      }*/
//
//
//      // TODO Check all processors to see if a repeat is necessary
//
//    } // while repeat
//}

/** Calculate integer coordinates of neighbors, taking into acount BCs.
  * Currently assuming Periodic in all directions.
  */
std::vector<int> btUnit::calcNeighIntCoords(unsigned lev, unsigned* lcoord, int* gCell, std::shared_ptr<BittreeAmr> mesh) {
    auto tree = mesh->getTree();

    std::vector<int> neighCoord(AMREX_SPACEDIM);

//--Calculate integer coordinates of neighbor in direction
    for(unsigned d=0;d<AMREX_SPACEDIM;++d)
      neighCoord[d] = static_cast<int>(lcoord[d]) + gCell[d];

//--Make sure not out-of-bounds. If periodic BCs, apply modulo
    std::vector<int> maxcoord(AMREX_SPACEDIM);
    for(unsigned d=0;d<AMREX_SPACEDIM;++d)
      maxcoord[d] = static_cast<int>(tree->top_size(d)) << lev;

    constexpr unsigned PERIODIC = 0;
    constexpr unsigned REFLECTING = 1;
    constexpr unsigned OUTFLOW = 2;
    std::vector<unsigned> bcLo(AMREX_SPACEDIM, PERIODIC);
    std::vector<unsigned> bcHi(AMREX_SPACEDIM, PERIODIC);

    for(unsigned d=0;d<AMREX_SPACEDIM;++d) {
      if (neighCoord[d] < 0 ) {
        if ( bcLo[d] == PERIODIC )
          neighCoord[d] = neighCoord[d] + maxcoord[d];
        else if ( bcLo[d] == REFLECTING )
          neighCoord[d] = static_cast<int>(lcoord[d]) - gCell[d];
        else
          neighCoord[d] = -1;
      }

      if (neighCoord[d] >= maxcoord[d]) {
        if ( bcHi[d] == PERIODIC )
          neighCoord[d] = neighCoord[d] + maxcoord[d];
        else if ( bcHi[d] == REFLECTING )
          neighCoord[d] = static_cast<int>(lcoord[d]) - gCell[d];
        else
          neighCoord[d] = -1;
      }

    }

    return neighCoord;
}

}
