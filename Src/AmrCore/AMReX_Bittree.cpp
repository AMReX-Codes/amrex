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
    mesh->refine_update();

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
    mesh->refine_update();

    //btCheckDerefine(mesh);


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
    // Tree before refinement.
    auto tree0 = mesh->getTree();

    // Ref test is marked 1 if block needs a tag (and doesn't have one).
    unsigned id0 = tree0->level_id0(0);
    unsigned id1 = tree0->id_upper_bound();
    std::vector<int> ref_test(tree0->id_upper_bound());

    // Repeat is made true if another round is needed
    bool repeat = false;

    do {
        // Clear out ref_test
        std::fill(ref_test.begin(),ref_test.end(),0);

        // Check adjacent children of neighbors of leaf blocks to see if they
        // are marked for refinement.
        // TODO only check local blocks?
        for( unsigned id = id0; id < id1; ++id) {
            auto b = tree0->locate(id);
            if( !b.is_parent && btTags[id]!=1 ) {
                bool needsTag = checkNeighborsRefine( mesh, b);
                if(needsTag) {
                    ref_test[id] = 1;
                }
            }
        }

        repeat = false;
        for( unsigned id = id0; id < id1; ++id) {
            if( ref_test[id]==1 && btTags[id]!=1 ) {
                repeat = true;
                btTags[id] = 1;
                mesh->refine_mark(id,true);
            }
        }

        // TODO Check all processors to see if a repeat is necessary
        //if(repeat) mesh->refine_reduce(comm);
        mesh->refine_update();

    } while(repeat);
}


/** Implements the logic which ensures the generated Bittree adheres
  * to a strict octree structure with no more than one level difference
  * between surrounding leaf blocks.
  */
void btUnit::btCheckDerefine( std::shared_ptr<BittreeAmr> mesh, std::vector<int>& btTags, MPI_Comm comm ) {
    // Tree before refinement. With only one rank, lnblocks = nblocks.
    auto tree0 = mesh->getTree();

    unsigned id0 = tree0->level_id0(0);
    unsigned id1 = tree0->id_upper_bound();
    std::vector<int> deref_test(tree0->id_upper_bound());

    // Repeat is made true if another round is needed
    bool repeat = false;

//--Repeat is left true if another round is needed
    do {
        //Turn deref_test to default 0 if block can't be derefined
        deref_test = btTags;

//------Check neighbors - if any neighbor is either a parent
//------or marked for refinement, do not derefine.
        for( unsigned id = id0; id < id1; ++id) {
            auto b = tree0->locate(id);
            if( !b.is_parent && btTags[id]==-1 ) {
                bool canDeref = checkNeighborsDerefine( mesh, b);
                if(!canDeref) {
                    deref_test[id] = 0;
                }
            }
        }

        repeat = false;
        for( unsigned id = id0; id < id1; ++id) {
            if( deref_test[id]==0 && btTags[id]==-1 ) {
                repeat = true;
                btTags[id] = 0;
                // Unmark for derefinement
                // TODO: for optimization have BT library find id of parent
                auto b = tree0->locate(id);
                unsigned parCoord[AMREX_SPACEDIM];
                for(unsigned d=0; d<AMREX_SPACEDIM; ++d)
                    parCoord[d] = b.coord[d]/2;
                auto p = tree0->identify(b.level-1,parCoord);

                mesh->refine_mark(p.id, false);
            }
        }

//------If any blocks are still marked for derefinement, check to make
//------sure their parents are still marked on bittree. This ensures blocks
//------only derefine if ALL siblings are marked for derefine.
        for( unsigned id = id0; id < id1; ++id) {
            if( btTags[id]==-1 ) {
                // TODO: for optimization have BT library find id of parent
                auto b = tree0->locate(id);
                unsigned parCoord[AMREX_SPACEDIM];
                for(unsigned d=0; d<AMREX_SPACEDIM; ++d)
                    parCoord[d] = b.coord[d]/2;
                auto p = tree0->identify(b.level-1,parCoord);

                bool deref_mark = mesh->check_refine_bit(p.id);
                if(!deref_mark) {
                  repeat = true;
                  btTags[id] = 0;
                }
            }
        }


        // TODO Check all processors to see if a repeat is necessary
        //if(repeat) mesh->refine_reduce_and(comm);
        mesh->refine_update();

    } while(repeat);
}


// Check all neighbors to see if their children are marked for refinement.
// If so, the block in question also needs to be refined.
bool btUnit::checkNeighborsRefine( std::shared_ptr<BittreeAmr> mesh, MortonTree::Block b) {
    auto tree0 = mesh->getTree();
    auto tree1 = mesh->getTree(true);
    int nIdx[3];
    unsigned cIdx[3];
    unsigned childCoord_u[AMREX_SPACEDIM];

    // Loop over neighbors
    for(nIdx[2]= -K3D; nIdx[2]<= K3D; ++nIdx[2]) {
    for(nIdx[1]= -K2D; nIdx[1]<= K2D; ++nIdx[1]) {
    for(nIdx[0]= -K1D; nIdx[0]<= K1D; ++nIdx[0]) {
        std::vector<int> nCoord = neighIntCoords(mesh, b.level, b.coord, nIdx);

        // If neighbor is outside domain or otherwise invalid, continue.
        if(AMREX_D_TERM(nCoord[0]<0, || nCoord[1]<0, || nCoord[2]<0 )) {
            continue;
        }

        // Identify neighbor from Bittree.
        unsigned neighCoord_u[AMREX_SPACEDIM];
        for(unsigned d=0; d<AMREX_SPACEDIM; ++d) {
            neighCoord_u[d] = static_cast<unsigned>(nCoord[d]);
        }
        auto n = tree0->identify(b.level, neighCoord_u);

        if(b.level==n.level && n.is_parent) {
            // Loop over children of neighbor.
            for(cIdx[2]= 0; cIdx[2]<= K3D; ++cIdx[2]) {
            for(cIdx[1]= 0; cIdx[1]<= K2D; ++cIdx[1]) {
            for(cIdx[0]= 0; cIdx[0]<= K1D; ++cIdx[0]) {

                // Only check adjacent children
                if (( ((1-nIdx[0])/2)==cIdx[0] || nIdx[0] == 0 ) &&
                    ( ((1-nIdx[1])/2)==cIdx[1] || nIdx[1] == 0 ) &&
                    ( ((1-nIdx[2])/2)==cIdx[2] || nIdx[2] == 0 )) {

                    // Identify child on updated tree
                    for(unsigned d=0; d<AMREX_SPACEDIM; ++d) {
                      childCoord_u[d] = neighCoord_u[d]*2 + cIdx[d];
                    }
                    auto c = tree1->identify(n.level+1, childCoord_u);

                    // If child WILL be parent, need to tag block
                    if( c.level==(b.level+1) && c.is_parent) {
                        return true;
                    }
                }
            }}}
        }
    }}}

    // Don't need to tag block
    return false;
}

// Check all neighbors to see if it is a parent or marked for refinement.
// If so, the block in question cannot be derefined.
bool btUnit::checkNeighborsDerefine( std::shared_ptr<BittreeAmr> mesh, MortonTree::Block b) {
    auto tree0 = mesh->getTree();
    auto tree1 = mesh->getTree(true);
    int nIdx[3];
    unsigned cIdx[3];
    unsigned childCoord_u[AMREX_SPACEDIM];

    // Loop over neighbors
    for(nIdx[2]= -K3D; nIdx[2]<= K3D; ++nIdx[2]) {
    for(nIdx[1]= -K2D; nIdx[1]<= K2D; ++nIdx[1]) {
    for(nIdx[0]= -K1D; nIdx[0]<= K1D; ++nIdx[0]) {
        std::vector<int> nCoord = neighIntCoords(mesh, b.level, b.coord, nIdx);

        // If neighbor is outside domain or otherwise invalid, continue.
        if(AMREX_D_TERM(nCoord[0]<0, || nCoord[1]<0, || nCoord[2]<0 )) {
            continue;
        }

        // Identify neighbor on updated tree
        unsigned neighCoord_u[AMREX_SPACEDIM];
        for(unsigned d=0; d<AMREX_SPACEDIM; ++d) {
            neighCoord_u[d] = static_cast<unsigned>(nCoord[d]);
        }
        auto n = tree1->identify(b.level, neighCoord_u);

        // If neighbor WILL be parent, cannot derefine
        if(b.level==n.level && n.is_parent) {
            return false;
        }
    }}}

    // Allowed to derefine
    return true;
}


/** Calculate integer coordinates of neighbors, taking into acount BCs.
  * Currently assuming Periodic in all directions.
  */
std::vector<int> btUnit::neighIntCoords(std::shared_ptr<BittreeAmr> mesh,
                                   unsigned lev, unsigned* lcoord, int* gCell) {
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
