#include <iostream>

#include <AMReX.H>
#include "NeighborListParticleContainer.H"

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
    
    ParmParse pp;
    
    int size, max_step, max_grid_size, nlevs;
    bool write_particles, do_nl;
    Real dt;

    pp.get("size", size);
    pp.get("max_step", max_step);
    pp.get("max_grid_size", max_grid_size);
    pp.get("write_particles", write_particles);
    pp.get("dt", dt);
    pp.get("do_nl", do_nl);
    pp.get("nlevs", nlevs);

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, size);
    }

    RealBox fine_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
       fine_box.setLo(n,0.25*size);
       fine_box.setHi(n,0.75*size);
    }
    
    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(size - 1, size - 1, size - 1));
    const Box domain(domain_lo, domain_hi);

    Vector<int> rr(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev-1] = 2;
    
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) 
        is_per[i] = 1;

    // This defines a Geometry object which is useful for writing the plotfiles  
    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
	geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
			 &real_box, CoordSys::cartesian, is_per);
    }

    Vector<BoxArray> ba(nlevs);
    ba[0].define(domain);
    
    // Now we make the refined level be the center eighth of the domain
    if (nlevs > 1) {
        int n_fine = size*rr[0];
        IntVect refined_lo(D_DECL(n_fine/4,n_fine/4,n_fine/4)); 
        IntVect refined_hi(D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba[1].define(refined_patch);
    }
    
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(max_grid_size);
    }

    Vector<DistributionMapping> dmap(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap[lev] = DistributionMapping{ba[lev]};
    }
    
    int num_neighbor_cells = 1;

    NeighborListParticleContainer myPC(geom, dmap, ba, rr, num_neighbor_cells);

    myPC.InitParticles();

    for (int i = 0; i < max_step; i++) {
        if (write_particles) myPC.writeParticles(i);
        
        myPC.fillNeighbors();

        if (do_nl) { myPC.computeForcesNL(); } 
        else {       myPC.computeForces();   }

        myPC.clearNeighbors();

        myPC.moveParticles(dt);

        myPC.Redistribute();
    }

    if (write_particles) myPC.writeParticles(max_step);

    }
    
    amrex::Finalize();
}
