#include <AMReX_Bittree.H>

#include <iostream>
#include <functional>

// static variable initialization
std::vector<bool> btUnit::refine;
std::vector<bool> btUnit::derefine;
std::vector<std::vector<unsigned>> btUnit::lcoord;
std::vector<unsigned> btUnit::lev;
std::vector<unsigned> btUnit::bitid;
std::vector<bool> btUnit::is_par;
std::vector<std::vector<double>> btUnit::error;
std::vector<std::vector<double>> btUnit::error_par;

//NOTE: Create bittree in AmrMesh::MakeNewGrids (Real time)
//   with
//      `mesh = std::make_shared<BittreeAmr>(top,includes);`

/* Structure of AmrMesh::MakeNewGrids (int lbase, Real time, int& new_finest, Vector<BoxArray>& new_grids)
The BT version of MakeNewGrids should have three steps:
    - Error Estimation and tagging - this is TBD, can use either the
      below routine or a more AMReX-style tagging

    - Bitree mesh generated - btRefineInitalize

    - AMReX updates grids based on bitree - adapt gr_btMakeNewGridsCallback.F90
*/

/** Error estimation and tagging for refinement/derefinement.
  * This is application specific, needs a callback to some error
  * estimation routine.
  *
  * Requirements: two lists are created: `refine` and `derefine`, with
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
                      lnblocks, std::vector<unsigned>(NDIM,0u) );
    lev       = std::vector<unsigned>(lnblocks, 0 );
    is_par    = std::vector<bool>(lnblocks, false);
    bitid     = std::vector<unsigned>(lnblocks, 0 );
    error     = std::vector<std::vector<double>>(
                      lnblocks, std::vector<double>(NVARS,0.0) );
    error_par = std::vector<std::vector<double>>(
                      lnblocks, std::vector<double>(NVARS,0.0) );


    // Loop over local blocks to cache metadata in order of morton curve.
    // In this case I used bittree functionality to actually get the
    // coordinates and level of the list of blocks, but
    // in AMReX this data can be gotten elsewhere obviously.
    unsigned id0 = tree0->level_id0(0);
    unsigned id1 = tree0->id_upper_bound();
    for( unsigned lb = id0; lb < id1; ++lb) {
        MortonTree::Block actual = tree0->locate(lb);

        // Bittree calcualtes blkID - the local ID the block would have
        // in paramesh
        MortonTree::Block b = tree0->identify( actual.level, actual.coord);
        unsigned blkID = b.mort;

        // Should probably check to make sure actual.lev == b.lev

        // Cache lev, lcoord, and is_par for later reference
        for(unsigned d=0; d<NDIM; ++d)
            lcoord[blkID][d] = actual.coord[d];
        lev[blkID] = actual.level;
        is_par[blkID] = tree0->block_is_parent(b.id);
        bitid[blkID] = b.id;

        // Estimate error
        double error_calc_result;
        for(unsigned v=0; v<NVARS; ++v) {
            // TODO replace with actual calculation
            error_calc_result = 0; 
            for(unsigned d=0; d<NDIM; ++d)
              if(lcoord[blkID][d]==0) error_calc_result += 1.0/NDIM;

            error[blkID][v] = error_calc_result;
        }
    }


    // Exchange error between parents and children. On one rank, this is
    // just having each parent give its error to all of its children.
    for( unsigned v=0; v<NVARS; ++v) {
      for( unsigned lb = 0; lb<lnblocks; ++lb) { //lb now is local blk num
        if(is_par[lb]) {
          unsigned ch[3];
          unsigned lev_ch;
          unsigned coord_ch[NDIM];
          for(ch[2]=0; ch[2]<=K3D; ++ch[2]) {
          for(ch[1]=0; ch[1]<=K2D; ++ch[1]) {
          for(ch[0]=0; ch[0]<=K1D; ++ch[0]) {
              lev_ch = lev[lb] + 1;
              for(unsigned d=0; d<NDIM; ++d) {
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
    unsigned maxLevel = 4;
    unsigned minLevel = 0;

    for( unsigned lb = 0; lb<lnblocks; ++lb) {
      allVarsDeref = true;
      for(unsigned v=0; v<NVARS; ++v) { 

        // These are application parameters, could be different per var
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
  * This routine should be more-or-less copy pasted into AMReX.
  */
void btUnit::btRefine( std::shared_ptr<BittreeAmr> mesh ) {
    // Tree before refinement. With only one rank, lnblocks = nblocks.
    auto tree0 = mesh->getTree();
    unsigned lnblocks = tree0->blocks();

//--Initialize bittree refinement and mark leaves to be refined
    //mesh->refine_init();

    for( unsigned lb = 0; lb<lnblocks; ++lb) {
      if (refine[lb]) {
        mesh->refine_mark(bitid[lb], true);
      }
    }
    //TODO mesh->refine_reduce(meshComm);

    btCheckRefine(mesh);

    // mark for derefinement
    for( unsigned lb = 0; lb<lnblocks; ++lb) {
      if (derefine[lb]) {
        unsigned parCoord[NDIM];
        for(unsigned d=0; d<NDIM; ++d) parCoord[d] = lcoord[lb][d]/2;
        
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

/** Creates new box arrays to match the new Bittree mesh.
  */
void btUnit::btMakeNewGrids(std::shared_ptr<BittreeAmr> mesh, int lbase,
                            Real time,int& new_finest,
                            Vector<BoxArray>& new_grids,
                            Vector<DistributionMapping>& new_dm) {
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

      //int lev_idx = 0
      for(int i=id0; i<id1; ++i) {
        //Get coordinates and morton index.
        auto b = tree1->locate(i);

        if(b.level != lev) {
            std::string msg = "Error identifying block in btMakeNewGrids";
        }

        //Calculate the processor based on position in the global Morton curve.
        //(This also calculates a local index that Paramesh uses, ignored here.)
        //call calc_proc_locblk(nprocs,localMortUB,mort,proc,locblk)
        int proc = 0;

        lo = {NXB*b.coord[0], NYB*b.coord[1], NZB*b.coord[2]};
        hi = {NXB*(b.coord[0]+1)-1, NYB*(b.coord[1]+1)-1, NZB*(b.coord[2]+1)-1};
        bl.push_back( Box{lo,hi} );
        pmap.push_back(proc);

        //lev_idx = lev_idx + 1;
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
            // TODO continue if any(neighCoords<0)
            unsigned neighCoord_u[NDIM];
            for(unsigned d=0; d<NDIM; ++d) neighCoord_u[d] = static_cast<unsigned>(neighCoord[d]);
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
            unsigned neighCoord_u[NDIM];
            for(unsigned d=0; d<NDIM; ++d) neighCoord_u[d] = static_cast<unsigned>(neighCoord[d]);
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
  */
std::vector<int> btUnit::calcNeighIntCoords(unsigned lev, unsigned* lcoord, int* gCell, std::shared_ptr<BittreeAmr> mesh) {
    auto tree = mesh->getTree();

    std::vector<int> neighCoord(NDIM);

//--Calculate integer coordinates of neighbor in direction
    for(unsigned d=0;d<NDIM;++d)
      neighCoord[d] = static_cast<int>(lcoord[d]) + gCell[d];

//--Make sure not out-of-bounds. If periodic BCs, apply modulo
    std::vector<int> maxcoord(NDIM);
    for(unsigned d=0;d<NDIM;++d)
      maxcoord[d] = static_cast<int>(tree->top_size(d)) << lev;

    constexpr unsigned PERIODIC = 0;
    constexpr unsigned REFLECTING = 1;
    std::vector<unsigned> bcLo(NDIM, PERIODIC);
    std::vector<unsigned> bcHi(NDIM, PERIODIC);

    for(unsigned d=0;d<NDIM;++d) {
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
