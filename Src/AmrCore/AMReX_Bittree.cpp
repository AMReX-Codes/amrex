#include <AMReX_Bittree.H>

#include <iostream>
#include <functional>

#define NVARS 2

namespace amrex {

std::vector<bool> btUnit::refine;
std::vector<bool> btUnit::derefine;
std::vector<std::vector<unsigned>> btUnit::lcoord;
std::vector<unsigned> btUnit::lev;
std::vector<unsigned> btUnit::bitid;
std::vector<bool> btUnit::is_par;
std::vector<std::vector<double>> btUnit::error;
std::vector<std::vector<double>> btUnit::error_par;


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

/*
NOTE: Bittree object is created in AmrMesh::MakeNewGrids (Real time)
   with
      `mesh = std::make_shared<BittreeAmr>(top,includes);`

The functions here are called in the BT version of MakeNewGrids which has three steps:
    1. Error Estimation and tagging - btTagging

    2. Bitree's actual bitmap generated/updated - btRefine

    3. AMReX updates grids based on bitree - btCalculateGrids
*/

//---------------------------------------------------------------------
//  Routines that need revising
//---------------------------------------------------------------------
/** Error estimation and tagging for refinement/derefinement.
  * This is implemented below in a paramesh-like algorithm
  * but this could be all changed to fit AMReX.
  *
  * The requirement is that two lists are created: `refine` and `derefine`, with
  * flags for all local blocks marking them for refinement/derefinement.
  */
void btUnit::btErrorEst( std::shared_ptr<BittreeAmr> mesh ) {
    // Tree before refinement. With only one rank, lnblocks = nblocks.
    auto tree0 = mesh->getTree();
    unsigned lnblocks = tree0->blocks();

    // Initialize the mesh data caches to appropriate values
    refine    = std::vector<bool>(lnblocks, false);
    derefine  = std::vector<bool>(lnblocks, false);
    lcoord    = std::vector<std::vector<unsigned>>(
                      lnblocks, std::vector<unsigned>(AMREX_SPACEDIM,0u) );
    lev       = std::vector<unsigned>(lnblocks, 0 );
    is_par    = std::vector<bool>(lnblocks, false);
    bitid     = std::vector<unsigned>(lnblocks, 0 );
    error     = std::vector<std::vector<double>>(
                      lnblocks, std::vector<double>(NVARS,0.0) );
    error_par = std::vector<std::vector<double>>(
                      lnblocks, std::vector<double>(NVARS,0.0) );


    // Loop over local blocks to cache metadata in order of morton curve.
    unsigned id0 = tree0->level_id0(0);
    unsigned id1 = tree0->id_upper_bound();
    for( unsigned lb = id0; lb < id1; ++lb) {
        // In this case I used bittree functionality to actually get the
        // coordinates and level of the list of blocks, but
        // in AMReX this data can be gotten elsewhere obviously.
        MortonTree::Block actual = tree0->locate(lb);

        // Bittree calcualtes blkID - such that the pair (rank, blkID)
        // uniquely identifies a block (grid). For one rank, this is just equal
        // to the index of the block in the morton curve.
        MortonTree::Block b = tree0->identify( actual.level, actual.coord);
        unsigned blkID = b.mort;

        // Should probably check to make sure actual.lev == b.lev

        // Cache lev, lcoord, and is_par for later reference
        for(unsigned d=0; d<AMREX_SPACEDIM; ++d)
            lcoord[blkID][d] = actual.coord[d];
        lev[blkID] = actual.level;
        is_par[blkID] = tree0->block_is_parent(b.id);
        bitid[blkID] = b.id;

        // Estimate error
        double error_calc_result;
        for(unsigned v=0; v<NVARS; ++v) {
            // TODO actually do error estimation and tagging
            error[blkID][v] = 0.5;
        }
    }


    // Exchange error between parents and children. On one rank, this is
    // just having each parent give its error to all of its children.
    for( unsigned v=0; v<NVARS; ++v) {
      for( unsigned lb = 0; lb<lnblocks; ++lb) { //lb now is local blk num
        if(is_par[lb]) {
          unsigned ch[3];
          unsigned lev_ch;
          unsigned coord_ch[AMREX_SPACEDIM];
          for(ch[2]=0; ch[2]<=K3D; ++ch[2]) {
          for(ch[1]=0; ch[1]<=K2D; ++ch[1]) {
          for(ch[0]=0; ch[0]<=K1D; ++ch[0]) {
              lev_ch = lev[lb] + 1;
              for(unsigned d=0; d<AMREX_SPACEDIM; ++d) {
                  coord_ch[d] = lcoord[lb][d]*2 + ch[d];
              }
              auto b = tree0->identify(lev_ch, coord_ch);

              error_par[b.mort][v] = error[lb][v];
          }}}
        }
      }
    }

    // Actual marking for refinement/derefinement
    double refineCutoff, derefineCutoff;
    bool allVarsDeref;
    unsigned maxLevel = 2;
    unsigned minLevel = 0;

    for( unsigned lb = 0; lb<lnblocks; ++lb) {
      allVarsDeref = true;
      for(unsigned v=0; v<NVARS; ++v) { 

        // These are application parameters, could be different per var
        // Since error is hard coded to 0.5 above, nothing is marked for
        // refinement or derefinement.
        refineCutoff = 0.99;
        derefineCutoff = 0.0; 

        // If block's error is too large for any single refinement variable,
        // then the block should be refined. The block's error is too small
        // for ALL of the variables, the block should be derefined.
        if( lev[lb]<maxLevel &&
            !is_par[lb] &&
            (error[lb][v]>refineCutoff || lev[lb]<minLevel)) {
          refine[lb] = true;
          derefine[lb] = false; 
        }

        if(refine[lb]) break; // no need to do other vars

        if( lev[lb] > minLevel  && !is_par[lb] ) {
          if((error[lb][v] <= derefineCutoff  &&
              error_par[lb][v] <= derefineCutoff && allVarsDeref)
             || lev[lb] > maxLevel ) {
            derefine[lb] = true;
          }
        }

        if( error[lb][v] > derefineCutoff ||
            error_par[lb][v] > derefineCutoff) {
           allVarsDeref = false;
           derefine[lb] = false;
        }
      }
    }
}

/** New Bittree mesh is generated.
  *
  * This makes use of BT library functions and as well as routines adapted
  * from Flash-X that enforce Octree logic.
  */
void btUnit::btRefine( std::shared_ptr<BittreeAmr> mesh ) {
    // Tree before refinement. With only one rank, lnblocks = nblocks.
    auto tree0 = mesh->getTree();
    unsigned lnblocks = tree0->blocks();

//--Mark leaves to be refined
    for( unsigned lb = 0; lb<lnblocks; ++lb) {
      if (refine[lb]) {
        mesh->refine_mark(bitid[lb], true);
      }
    }
    //TODO mesh->refine_reduce(meshComm);

    btCheckRefine(mesh);

//--Mark leaves for derefinement (on BT, parent is marked for nodetype change)
    for( unsigned lb = 0; lb<lnblocks; ++lb) {
      if (derefine[lb]) {
        unsigned parCoord[AMREX_SPACEDIM];
        for(unsigned d=0; d<AMREX_SPACEDIM; ++d) parCoord[d] = lcoord[lb][d]/2;
        
        auto p = tree0->identify(lev[lb]-1,parCoord);
        mesh->refine_mark(p.id, true);
      }
    }
    //TODO mesh->refine_reduce(meshComm);

    btCheckDerefine(mesh);

//--Generate updated Bittree
    mesh->refine_update();

    //TODO sort (distribute over procs)

}


//---------------------------------------------------------------------
// Local Routines
//---------------------------------------------------------------------

/** Implements the logic which ensures the generated Bittree adheres
  * to a strict octree structure with no more than one level difference
  * between surrounding leaf blocks.
  */
void btUnit::btCheckRefine( std::shared_ptr<BittreeAmr> mesh ) {
    // Tree before refinement. With only one rank, lnblocks = nblocks.
    auto tree0 = mesh->getTree();
    unsigned lnblocks = tree0->blocks();


    bool repeat = true;

//--Repeat is left true if another round is needed
    while(repeat) {
      std::vector<bool> ref_test( refine.size(), false );

//----Check adjacent children of neighbors of local leaf blocks to see if they
//----are marked for refinement.
      for( unsigned lb = 0; lb<lnblocks; ++lb) {
        if( !is_par[lb] && !refine[lb] ) {

//--------Loop over neighbors
          for(int lk= -K3D; lk<= K3D; ++lk) {
          for(int lj= -K2D; lj<= K2D; ++lj) {
          for(int li= -K1D; li<= K1D; ++li) {
            std::vector<int> gCell{li,lj,lk};

            auto neighCoord = calcNeighIntCoords(lev[lb], lcoord[lb].data(), gCell.data(), mesh);
            // TODO skip if any(neighCoords<0)
            unsigned neighCoord_u[AMREX_SPACEDIM];
            for(unsigned d=0; d<AMREX_SPACEDIM; ++d) neighCoord_u[d] = static_cast<unsigned>(neighCoord[d]);
            auto b = tree0->identify(lev[lb], neighCoord_u);
            if (b.level == lev[lb]) { //neighbor exists
              if(b.is_parent) {
                // Check its children

                // TODO loop over children
                //  { auto c = tree->identify(lev[lb]+1, childCoord);
                //    if(c.is_par) ref_test[lb] = true; }
              }
            }
          }}} // neighbor loop
        } // if leaf block not already marked
      } // local block loop

      repeat = false;

      for(unsigned lb=0; lb<lnblocks; ++lb) {
          if( ref_test[lb] && !refine[lb] ) {
              repeat = true;
              refine[lb] = true;
              derefine[lb] = false;

              mesh->refine_mark(bitid[lb],true);
          }
      }

      // TODO Check all processors to see if a repeat is necessary
      // TODO if(repeat) mesh->refine_reduce(meshComm);

    } // while repeat
}


/** Implements the logic which ensures the generated Bittree adheres
  * to a strict octree structure with no more than one level difference
  * between surrounding leaf blocks.
  */
void btUnit::btCheckDerefine( std::shared_ptr<BittreeAmr> mesh ) {
    // Tree before refinement. With only one rank, lnblocks = nblocks.
    auto tree = mesh->getTree();
    unsigned lnblocks = tree->blocks();


    bool repeat = true;

//--Repeat is left true if another round is needed
    while(repeat) {
      //Turn OFF deref_test if block can't be derefined
      std::vector<bool> deref_test = derefine;

//----Check neighbors - if any neighbor is either a parent 
//----or marked for refinement, do not derefine. 
      for( unsigned lb = 0; lb<lnblocks; ++lb) {
        if( derefine[lb] ) {

//--------Loop over neighbors
          for(int lk= -K3D; lk<= K3D; ++lk) {
          for(int lj= -K2D; lj<= K2D; ++lj) {
          for(int li= -K1D; li<= K1D; ++li) {
            std::vector<int> gCell{li,lj,lk};

            auto neighCoord = calcNeighIntCoords(lev[lb], lcoord[lb].data(), gCell.data(), mesh);
            // TODO continue if any(neighCoords<0)
            unsigned neighCoord_u[AMREX_SPACEDIM];
            for(unsigned d=0; d<AMREX_SPACEDIM; ++d) neighCoord_u[d] = static_cast<unsigned>(neighCoord[d]);
            auto b = tree->identify(lev[lb], neighCoord_u);

            if (b.level == lev[lb]) { //neighbor exists
              bool neighRefine = mesh->check_refine_bit(b.id);
              if(neighRefine) {
                deref_test[lb] = false;
                // TODO break neighbor loop
              }
            }

            /* TODO
            if( b.is_parent) {
              if( any(lcoord[lb]/2 != neighCoord/2) ) {
                bool deref_mark = mesh->check_refine_bit(b.id);
                if(!deref_mark) {
                  deref_test[lb] = false;
                  // break
                }
              } else {
                  deref_test[lb] = false;
                  // break
              }
            } // neighbor is parent */
          }}} // neighbor loop
        } // if leaf block not already marked
      } // local block loop

      repeat = false;

      for(unsigned lb=0; lb<lnblocks; ++lb) {
          if( !deref_test[lb] && derefine[lb] ) {
              repeat = true;
              derefine[lb] = false;

              // Unmark for derefinement
          }
      }

      //mesh->refine_reduce_and(meshComm);

//----If any blocks are still marked for derefinement, check to make
//----sure their parents are still marked on bittree. This ensures blocks
//----only derefine if ALL siblings are marked for derefine.
      /*
      for(unsigned lb=0; lb<lnblocks; ++lb) {
        if( derefine[lb] ) {
          auto p = tree->(lev[lb]-1,lcoord/2);
          bool deref_mark = mesh->check_refine_bit(p.id);
          if(!deref_mark) {
            repeat = true;
            derefine[lb] = false;
          }
        }
      }*/


      // TODO Check all processors to see if a repeat is necessary

    } // while repeat
}

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
    std::vector<unsigned> bcLo(AMREX_SPACEDIM, PERIODIC);
    std::vector<unsigned> bcHi(AMREX_SPACEDIM, PERIODIC);

    for(unsigned d=0;d<AMREX_SPACEDIM;++d) {
      if (neighCoord[d] < 0 ) {
        if ( bcLo[d] == PERIODIC )
          neighCoord[d] = neighCoord[d] % maxcoord[d];
        else if ( bcLo[d] == REFLECTING )
          neighCoord[d] = static_cast<int>(lcoord[d]) - gCell[d];
        else
          neighCoord[d] = -1;
      }

      if (neighCoord[d] >= maxcoord[d]) {
        if ( bcHi[d] == PERIODIC )
          neighCoord[d] = neighCoord[d] % maxcoord[d];
        else if ( bcHi[d] == REFLECTING )
          neighCoord[d] = static_cast<int>(lcoord[d]) - gCell[d];
        else
          neighCoord[d] = -1;
      }

    }

    return neighCoord;
}

}
