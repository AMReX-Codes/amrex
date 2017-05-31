#include <iostream>
#include <random>

#include <AMReX.H>
#include "ShortRangeParticleContainer.H"

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    
    ParmParse pp;
    
    int size, max_step, max_grid_size;
    bool write_particles;
    Real dt;

    pp.get("size", size);
    pp.get("max_step", max_step);
    pp.get("max_grid_size", max_grid_size);
    pp.get("write_particles", write_particles);
    pp.get("dt", dt);

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, size);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(size - 1, size - 1, size - 1));
    const Box domain(domain_lo, domain_hi);
    
    int coord = 0;
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) 
        is_per[i] = 0; 
    Geometry geom(domain, &real_box, coord, is_per);
    
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    
    DistributionMapping dmap(ba);
   
    int num_neighbor_cells = 1;
    ShortRangeParticleContainer myPC(geom, dmap, ba, num_neighbor_cells);

    myPC.InitParticles();

    for (int i = 0; i < max_step; i++) {
        if (write_particles) myPC.writeParticles(i);
        
        myPC.fillNeighbors();
        myPC.computeForces();
        myPC.clearNeighbors();

        myPC.moveParticles(dt);

        myPC.Redistribute();
    }

    if (write_particles) myPC.writeParticles(max_step);
    
    amrex::Finalize();
}
