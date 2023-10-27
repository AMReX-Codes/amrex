#include <AMReX_Bittree.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MFIter.H>
#include <functional>

using namespace bittree;

namespace amrex {
static constexpr auto K1D = int(AMREX_SPACEDIM>=1);
static constexpr auto K2D = int(AMREX_SPACEDIM>=2);
static constexpr auto K3D = int(AMREX_SPACEDIM>=3);


bool btUnit::bcPeriodic[AMREX_SPACEDIM];

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
  * from Flash-X that enforce Octree nesting.
  */
int btUnit::btRefine (BittreeAmr* const mesh, std::vector<int>& btTags,
                      int max_crse, int lbase,
                      Vector<BoxArray>& grids, Vector<DistributionMapping>& dmap, MPI_Comm comm)
{
    BL_PROFILE("Bittree-btRefine");

    // Tree before refinement. With only one rank, lnblocks = nblocks.
    auto tree0 = mesh->getTree();

    // Mark leaves to be refined
    for (int lev=max_crse; lev>=lbase; --lev) {
        for (MFIter mfi(grids[lev], dmap[lev]); mfi.isValid(); ++mfi) {
            int id = getBitid(mesh,false,lev,mfi.index());
            if (btTags[id]==1) {
                if(!tree0->block_is_parent(id)) {
                    mesh->refine_mark(id, true);
                }
            }
        }
    }

    mesh->refine_reduce(comm);
    mesh->refine_update();

    btCheckRefine(mesh, btTags, max_crse, lbase, grids, dmap, comm);

    // Mark derefinement (parents who will nodetype change to leaf)
    for (int lev=max_crse; lev>=lbase; --lev) {
        for (MFIter mfi(grids[lev], dmap[lev]); mfi.isValid(); ++mfi) {
            int id = getBitid(mesh,false,lev,mfi.index());
            if (btTags[id]==-1) {
                if(tree0->block_is_parent(id)) {
                    mesh->refine_mark(id, true);
                }
            }
        }
    }

    mesh->refine_reduce(comm);
    mesh->refine_update();

    btCheckDerefine(mesh, btTags, max_crse, lbase, grids, dmap, comm);

    // return delta count
    return static_cast<int>( mesh->delta_count() );
}

/** Creates new box arrays to match the new Bittree mesh.
  */
void btUnit::btCalculateGrids (BittreeAmr* const mesh, int lbase,
                               int& new_finest,
                               Vector<BoxArray>& new_grids,
                               Vector<IntVect> const& max_grid_size)
{
    BL_PROFILE("Bittree-btCalculateGrids");

    auto tree1 = mesh->getTree(true);
    auto nlevs = tree1->levels();
    new_finest = int(nlevs - 1);

//--Calculate the new grid layout and distribution map based on Bittree
    for(int lev=lbase; lev<=new_finest; ++lev) {
        btCalculateLevel(mesh, lev, new_grids[lev],
                         max_grid_size[lev]);
    }
}

/** Creates a box array based on Bittree.
  */
void btUnit::btCalculateLevel (BittreeAmr* const mesh, int lev,
                               BoxArray& ba,
                               IntVect const& max_grid_size)
{
    auto tree1 = mesh->getTree(true);

    //Bittree has its own indices for blocks which I call bitid; get
    //the range of bitids for the level being made. Bitid range is
    //contiguous for each level.
    auto id0 = tree1->level_id0(lev);
    auto id1 = tree1->level_id1(lev);
    // int nblocks = tree1->level_blocks(lev);

    BoxList bl;

    for(auto i=id0; i<id1; ++i) {
      //Get coordinates and morton index.
      auto b = tree1->locate(i);

      if(b.level != lev) {
          std::string msg = "Error identifying block in btCalculateGrids";
          //throw error?
      }

      IntVect coordVec{AMREX_D_DECL(static_cast<int>(b.coord[0]),
                                    static_cast<int>(b.coord[1]),
                                    static_cast<int>(b.coord[2]))};
      IntVect lo = max_grid_size*coordVec;
      IntVect hi = max_grid_size*(coordVec+1) - 1;
      bl.push_back( Box{lo,hi} );
    }

    ba = BoxArray(bl);
}

int btUnit::getBitid (BittreeAmr* const mesh, bool updated,
                      int lev, int idx_on_lev)
{
    return idx_on_lev + int(mesh->getTree(updated)->level_id0(lev));
}

int btUnit::getIndex (BittreeAmr* const mesh, bool updated,
                      int lev, int bitid)
{
    return bitid - int(mesh->getTree(updated)->level_id0(lev));
}



//---------------------------------------------------------------------
// Local Routines
//---------------------------------------------------------------------

/** Implements the logic which ensures the generated Bittree adheres
  * to a strict octree structure with no more than one level difference
  * between surrounding leaf blocks.
  */
void btUnit::btCheckRefine (BittreeAmr* const mesh, std::vector<int>& btTags,
                            int max_crse, int lbase,
                            Vector<BoxArray>& grids,
                            Vector<DistributionMapping>& dmap, MPI_Comm comm)
{
    BL_PROFILE("Bittree-btCheckRefine");

    // Tree before refinement.
    auto tree0 = mesh->getTree();

    // Ref test is marked 1 if block needs a tag (and doesn't have one).
    std::vector<int> ref_test(tree0->id_upper_bound());

    // Repeat is made true if another round is needed
    bool repeat = false;

    do {
        // Clear out ref_test
        std::fill(ref_test.begin(),ref_test.end(),0);

        // Check neighbors - if any adjacent child of a neighbor is either a parent
        // or marked for refinement, this block needs to be refined.
        for (int lev=max_crse; lev>=lbase; --lev) {
            for (MFIter mfi(grids[lev], dmap[lev]); mfi.isValid(); ++mfi) {
                int id = getBitid(mesh,false,lev,mfi.index());
                auto b = tree0->locate(id);
                if( !b.is_parent && btTags[id]!=1 ) {
                    bool needsTag = checkNeighborsRefine( mesh, b);
                    //amrex::Print() << "needsTag for " << id << " : " << needsTag <<std::endl;
                    if(needsTag) {
                        ref_test[id] = 1;
                    }
                }
            }
        }

        // Mark blocks who need to be refined (as per above check).
        repeat = false;
        for (int lev=max_crse; lev>=lbase; --lev) {
            for (MFIter mfi(grids[lev], dmap[lev]); mfi.isValid(); ++mfi) {
                int id = getBitid(mesh,false,lev,mfi.index());
                if( ref_test[id]==1 && btTags[id]!=1 ) {
                    repeat = true;
                    btTags[id] = 1;
                    mesh->refine_mark(id,true);
                }
            }
        }

        // If only processing local blocks, check all processors to see if
        // a repeat is necessary, then reduce bittree to update on all ranks.
        ParallelDescriptor::ReduceBoolOr(repeat);

        if(repeat) {
            mesh->refine_reduce(comm);
            mesh->refine_update();
        }

    } while(repeat);
}


/** Implements the logic which ensures the generated Bittree adheres
  * to a strict octree structure with no more than one level difference
  * between surrounding leaf blocks.
  */
void btUnit::btCheckDerefine (BittreeAmr* const mesh, std::vector<int>& btTags,
                              int max_crse, int lbase,
                              Vector<BoxArray>& grids,
                              Vector<DistributionMapping>& dmap, MPI_Comm comm)
{
    BL_PROFILE("Bittree-btCheckDerefine");

    // Tree before refinement. With only one rank, lnblocks = nblocks.
    auto tree0 = mesh->getTree();

    std::vector<int> deref_test(tree0->id_upper_bound());

    // Repeat is made true if another round is needed
    bool repeat = false;

    // Repeat is left true if another round is needed
    do {
        // Turn deref_test to default 0 if block can't be derefined
        deref_test = btTags;

        // Check neighbors - if any adjacent child of neighbor is either a parent
        // or marked for refinement, do not derefine.
        for (int lev=max_crse; lev>=lbase; --lev) {
            for (MFIter mfi(grids[lev], dmap[lev]); mfi.isValid(); ++mfi) {
                int id = getBitid(mesh,false,lev,mfi.index());
                auto b = tree0->locate(id);
                if( btTags[id]==-1 ) {
                    bool cantDeref = checkNeighborsRefine( mesh, b);
                    if(cantDeref) {
                        deref_test[id] = 0;
                    }
                }
            }
        }

        // Unmark any blocks who cannot derefine (as per above check).
        repeat = false;
        for (int lev=max_crse; lev>=lbase; --lev) {
            for (MFIter mfi(grids[lev], dmap[lev]); mfi.isValid(); ++mfi) {
                int id = getBitid(mesh,false,lev,mfi.index());
                if( deref_test[id]==0 && btTags[id]==-1 ) {
                    repeat = true;
                    btTags[id] = 0;

                    // Unmark for derefinement
                    mesh->refine_mark(id, false);
                }
            }
        }

        // If only processing local blocks, check all processors to see if
        // a repeat is necessary, then reduce bittree to update on all ranks.
        ParallelDescriptor::ReduceBoolOr(repeat);

        if(repeat) {
            mesh->refine_reduce_and(comm);
            mesh->refine_update();
        }

    } while(repeat);
}


// Check all neighbors to see if their adjacent children are parents or marked for refinement.
bool btUnit::checkNeighborsRefine (BittreeAmr* const mesh, MortonTree::Block b)
{
    BL_PROFILE("Bittree-checkNeighborsRefine");

    auto tree0 = mesh->getTree();
    auto tree1 = mesh->getTree(true);
    int nIdx[3], cIdx[3];
    unsigned childCoord_u[AMREX_SPACEDIM];

    // Loop over neighbors
    for(nIdx[2]= -1*K3D; nIdx[2]<= K3D; ++nIdx[2]) {
    for(nIdx[1]= -1*K2D; nIdx[1]<= K2D; ++nIdx[1]) {
    for(nIdx[0]= -1*K1D; nIdx[0]<= K1D; ++nIdx[0]) {
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
                      childCoord_u[d] = neighCoord_u[d]*2 + static_cast<unsigned>(cIdx[d]);
                    }
                    auto c = tree1->identify(n.level+1, childCoord_u);

                    // If child WILL be parent, return true
                    if( c.level==(b.level+1) && c.is_parent) {
                        return true;
                    }
                }
            }}}
        }
    }}}

    // Return false otherwise
    return false;
}

/** Calculate integer coordinates of neighbors, taking into account BCs.
  * Currently assuming Periodic in all directions.
  */
std::vector<int> btUnit::neighIntCoords (BittreeAmr* const mesh,
                                         unsigned lev, unsigned const* lcoord,
                                         int const* gCell)
{
    auto tree = mesh->getTree();

    std::vector<int> neighCoord(AMREX_SPACEDIM);

//--Calculate integer coordinates of neighbor in direction
    for(unsigned d=0;d<AMREX_SPACEDIM;++d) {
        neighCoord[d] = static_cast<int>(lcoord[d]) + gCell[d];
    }

//--Make sure not out-of-bounds. If periodic BCs, apply modulo
    std::vector<int> maxcoord(AMREX_SPACEDIM);
    for(unsigned d=0;d<AMREX_SPACEDIM;++d) {
        maxcoord[d] = static_cast<int>(tree->top_size(d)) << lev;
    }

    for(unsigned d=0;d<AMREX_SPACEDIM;++d) {
        if (neighCoord[d] < 0 ) {
            if ( bcPeriodic[d] == true ) {
                neighCoord[d] = neighCoord[d] + maxcoord[d];
            } else {
                neighCoord[d] = -1;
            }
        }

        if (neighCoord[d] >= maxcoord[d]) {
            if ( bcPeriodic[d] == true ) {
                neighCoord[d] = neighCoord[d] - maxcoord[d];
            } else {
                neighCoord[d] = -1;
            }
        }
    }

    return neighCoord;
}

}
