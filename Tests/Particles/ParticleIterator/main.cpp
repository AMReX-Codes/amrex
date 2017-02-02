#include <iostream>
#include <map>
#include <vector>

#include "AMReX_Array.H"
#include "AMReX_Particles.H"

using namespace amrex;

class Iterator
{

public:
  Iterator (const BoxArray& _box_array, bool _do_tiling, const IntVect& _tile_size)
    : box_array(_box_array), do_tiling(_do_tiling), tile_size(_tile_size)
  {
    Initialize();
  }
  
  bool isValid () {return current_grid_index < num_boxes;}

  void operator++ () {
    current_tile_index++;
    if (current_tile_index >= num_tiles[current_grid_index]) {
      current_grid_index++;
      current_tile_index = 0;
    }
  }

  int gridIndex () const {return current_grid_index;}

  int tileIndex () const {return current_tile_index;}

  Box gridBox () {return box_array[current_grid_index];}

  Box tileBox () {
    if (do_tiling) {
      const Box& box = box_array[current_grid_index];
      const IntVect& small = box.smallEnd();
      const IntVect& size  = box.size();
      IntVect ntiles(D_DECL((size[0]+tile_size[0]-1)/tile_size[0],
			    (size[1]+tile_size[1]-1)/tile_size[1],
			    (size[2]+tile_size[2]-1)/tile_size[2]));
      
      IntVect index = { D_DECL(0, 0, 0) };
      int linear_index = current_tile_index;
      for (int i = 0; i < BL_SPACEDIM; i++) {
	index[i] = linear_index % ntiles[i];
	linear_index /= ntiles[i];
      }
      
      IntVect tile_lo_corner = small + index*tile_size; 
      IntVect tile_hi_corner = tile_lo_corner + tile_size - IntVect::TheUnitVector(); 
      Box tile_box(tile_lo_corner, tile_hi_corner);
      tile_box &= box;
      return tile_box;
    }
    
    return box_array[current_grid_index];
  }

private:

  void Initialize()
  {
    current_grid_index = 0;
    current_tile_index = 0;
    num_boxes = box_array.size();
    num_tiles.resize(num_boxes, 0);
    if (do_tiling) {
      for (int i = 0; i < num_boxes; i++) {
	const Box& box = box_array[i];
	const IntVect& small = box.smallEnd();
	const IntVect& size  = box.size();
	IntVect ntiles(D_DECL((size[0]+tile_size[0]-1)/tile_size[0],
			      (size[1]+tile_size[1]-1)/tile_size[1],
			      (size[2]+tile_size[2]-1)/tile_size[2]));
	int ntottiles = D_TERM(ntiles[0],*ntiles[1],*ntiles[2]);
	num_tiles[i] = ntottiles;
      }
    }
  }

  bool do_tiling;
  IntVect tile_size;

  const BoxArray& box_array;

  int current_grid_index;
  int current_tile_index;
  int num_boxes;
  Array<int> num_tiles;
};

int main(int argc, char* argv[])
{

  amrex::Initialize(argc,argv);

  int ncell = 48;
  int max_grid_size = 32;
  int nlevs = 1;
  int coord = 0;

  bool do_tiling = true;
  IntVect tile_size = { D_DECL(102400, 8, 8) };

  RealBox real_box;
  for (int n = 0; n < BL_SPACEDIM; n++)
    {
      real_box.setLo(n,0.0);
      real_box.setHi(n,1.0);
    }

  IntVect domain_lo(0 , 0, 0); 
  IntVect domain_hi(ncell-1, ncell-1, ncell-1); 
 
  const Box domain(domain_lo, domain_hi);

  Array<int> rr(nlevs-1);
  for (int lev = 1; lev < nlevs; lev++)
    rr[lev-1] = 2;
 
  int is_per[BL_SPACEDIM];
  for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 1;

  Array<Geometry> geom(nlevs);
  geom[0].define(domain, &real_box, coord, is_per);

  Array<BoxArray> ba(nlevs);
  ba[0].define(domain);  

  for (int lev = 0; lev < nlevs; lev++)
    ba[lev].maxSize(max_grid_size);

  for (Iterator it(ba[0], do_tiling, tile_size); it.isValid(); ++it) {
    if (ParallelDescriptor::IOProcessor()) {
      Box tile_box = it.tileBox();
      Box grid_box = it.gridBox();
      std::cout << "Grid number: "  << it.gridIndex() << " " << grid_box << " " << grid_box.size();
      std::cout << " Tile number: " << it.tileIndex() << " " << tile_box << " " << tile_box.size();
      std::cout << std::endl;
    }
  }

  Array<DistributionMapping> dmap(nlevs);
  for (int lev = 0; lev < nlevs; lev++)
    dmap[lev].define(ba[lev]);
    
  amrex::Finalize();
}
